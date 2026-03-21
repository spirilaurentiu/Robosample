#include "World.hpp"
#include "Compound.h"
#include "Constraint.h"
#include "OpenMM.hpp"
#include "Sampler.hpp"
#include "TopologyElements.hpp"
#include "bgeneral.hpp"
#include "common.h"
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <ostream>

void World::setAtomTargetLocationsToState(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets)
{
	// Update cache
	// TODO so far we only need to update the cache for testing purposes only
	atomTargetLocaltionsCache = atomTargets;

	// Match Compound and DuMM coordinates
	for(std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++)
	{
		Topology& topology = topologies[topoIx];

		// Use Molmodel's Compound match functions to set the new conf
		topology.matchAtomTargetLocations(atomTargets[topoIx]);

		// Get the Ground to Top Transform
		const SimTK::Transform G_X_T = topology.getTopLevelTransform();
		SimTK::Vec3 locationInMobod = SimTK::Vec3(0);

		// Set atoms' stations on body
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topology.getNumAtoms(); ++aIx) {
			// Get previous location in mobod
			const SimTK::Vec3& locInMobod = topology.getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, *forceField);
			
			if(locInMobod == 0){
				// Atom is at body's origin
				locationInMobod = SimTK::Vec3(0);
			} else {
				// Atom is not at body's origin
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);

				const std::pair<int, SimTK::Compound::AtomIndex>& topoRootAtomPair = getMobodRootAtomIndex(mbx);
				SimTK::Compound::AtomIndex mobodRootAIx = topoRootAtomPair.second;

				const SimTK::Transform& T_X_root = topology.getTopTransform(mobodRootAIx);
				SimTK::Transform G_X_root = G_X_T * T_X_root;

				const SimTK::Vec3& G_vchild = atomTargets[topoIx].at(aIx);
				SimTK::Transform root_X_child = alignFlipAndTranslateFrameAlongXAxis(G_X_root, G_vchild);

				topology.bsetFrameInMobilizedBodyFrame(aIx, root_X_child);
				locationInMobod = root_X_child.p();
			}

			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
			SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(aIx);

			// Set station_B
			forceField->bsetAtomStationOnBody(dAIx, locationInMobod);
			forceField->bsetAllAtomStationOnBody(dAIx, locationInMobod);

			// Set included atom
			forceField->updIncludedAtomStation(dAIx) = locationInMobod;
			forceField->updAllAtomStation(dAIx) = locationInMobod;

			// Atom placements in clusters
			forceField->bsetAtomPlacementStation(dAIx, mbx, locationInMobod);
		}
	}

	// Set default child mobod inboard (X_PF) and outboard (X_BM) frames
	// This method only uses internal coordinates per molecule bonds info
	updateFramesFromTopologies();

	// Set every mobod's mass properties
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx)
	{
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		SimTK::DuMM::ClusterIndex clusterIx = forceField->bgetMobodClusterIndex(mbx);
		SimTK::MassProperties massProperties = forceField->calcClusterMassProperties(clusterIx);
		mobod.setDefaultMassProperties(massProperties);
	}

	// TODO add in the documentation that we advance to Position stage in this function
	const SimTK::State& state = compoundSystem->realizeTopology();
	compoundSystem->realize(state, SimTK::Stage::Position);
	
	integrator->updAdvancedState() = state;
	
	if (testing) {
		coordinateTransferErrors.push_back(checkCoordinateTransfer(atomTargets));
	}
}

CoordinateTransferError World::checkCoordinateTransfer(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets)
{
	auto wrapAngle = [](SimTK::Real x) {
        return std::atan2(std::sin(x), std::cos(x));
    };

	CoordinateTransferError error;
	int numAtoms = 0,
		numBonds = 0,
		numAngles = 0,
		numPropers = 0,
		numImpropers = 0;

	// Check topology atom target location matching residuals
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		SimTK::Real matchError = topologies[topoIx].getMatchError(atomTargets[topoIx]);
		error.matchResiduals.push_back(matchError);
	}

	// Collect new atom coordinates
	std::vector<std::vector<SimTK::Vec3>> newAtomtargets(topologies.size());

	const SimTK::State& state = integrator->updAdvancedState();
	compoundSystem->realize(state, SimTK::Stage::Position);

	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (SimTK::Compound::AtomIndex cAIx(0); cAIx < topologies[topoIx].getNumAtoms(); ++cAIx) {
			const auto& computedLoc = topologies[topoIx].calcAtomLocationInGroundFrameThroughSimbody(cAIx, *forceField, *matter, state);
			newAtomtargets[topoIx].push_back(computedLoc);
		}
	}

	// Run through all molecules
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {

		// Check per-atom cartesian displacements
		for (SimTK::Compound::AtomIndex cAIx(0); cAIx < topologies[topoIx].getNumAtoms(); ++cAIx) {
			const auto& computedLoc = newAtomtargets[topoIx][cAIx];
			const auto& targetLoc = atomTargets[topoIx].at(cAIx);

			const SimTK::Real diffNorm = (targetLoc - computedLoc).norm();
			error.cartesian += diffNorm;
			error.cartesianMax = std::max(error.cartesianMax, diffNorm);

			numAtoms++;
		}

		// BAT check - bond
		for (const auto& bond : topologies[topoIx].getBonds()) {
			const SimTK::Compound::AtomIndex parentAIx = bond.compoundAtomIndices[0];
			const SimTK::Compound::AtomIndex childAIx = bond.compoundAtomIndices[1];

			const auto& computedParent = newAtomtargets[topoIx][parentAIx];
			const auto& computedChild = newAtomtargets[topoIx][childAIx];

			const auto& targetParent = atomTargets[topoIx].at(parentAIx);
			const auto& targetChild = atomTargets[topoIx].at(childAIx);

			const SimTK::Real computedBondLength = (computedChild - computedParent).norm();
			const SimTK::Real targetBondLength = (targetChild - targetParent).norm();

			const SimTK::Real diffBondLength = std::abs(targetBondLength - computedBondLength);
			error.bonds += diffBondLength;
			error.bondsMax = std::max(error.bondsMax, diffBondLength);

			numBonds++;
		}

		// BAT check - angle
		for (const auto& angle : topologies[topoIx].getAngles()) {
			const SimTK::Compound::AtomIndex aIx1 = angle.compoundAtomIndices[0];
			const SimTK::Compound::AtomIndex aIx2 = angle.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex aIx3 = angle.compoundAtomIndices[2];

			const auto& computedPos1 = newAtomtargets[topoIx][aIx1];
			const auto& computedPos2 = newAtomtargets[topoIx][aIx2];
			const auto& computedPos3 = newAtomtargets[topoIx][aIx3];

			const auto& targetPos1 = atomTargets[topoIx].at(aIx1);
			const auto& targetPos2 = atomTargets[topoIx].at(aIx2);
			const auto& targetPos3 = atomTargets[topoIx].at(aIx3);

			SimTK::Real computedAngle = bAngle(computedPos1, computedPos2, computedPos3);
			SimTK::Real targetAngle = bAngle(targetPos1, targetPos2, targetPos3);

			const SimTK::Real diffAngle = std::abs(wrapAngle(targetAngle - computedAngle));
			error.angles += diffAngle;
			error.anglesMax = std::max(error.anglesMax, diffAngle);

			numAngles++;
		}

		// BAT chekck -	dihedrals
		for (const auto& torsion : topologies[topoIx].getPeriodicTorsions()) {
			const auto& ids = torsion.compoundAtomIndices;
		
			// Calculate computed dihedral
			SimTK::Real computedDihedral = bDihedral(
				newAtomtargets[topoIx][ids[0]], 
				newAtomtargets[topoIx][ids[1]], 
				newAtomtargets[topoIx][ids[2]], 
				newAtomtargets[topoIx][ids[3]]
			);

			// Calculate target dihedral
			SimTK::Real targetDihedral = bDihedral(
				atomTargets[topoIx].at(ids[0]), 
				atomTargets[topoIx].at(ids[1]), 
				atomTargets[topoIx].at(ids[2]), 
				atomTargets[topoIx].at(ids[3])
			);

			const SimTK::Real diffDihedral = std::abs(wrapAngle(targetDihedral - computedDihedral));
			if (torsion.improper) {
				error.improperDihedrals += diffDihedral;
				error.improperDihedralsMax = std::max(error.improperDihedralsMax, diffDihedral);
				numImpropers++;
			} else {
				error.properDihedrals += diffDihedral;
				error.properDihedralsMax = std::max(error.properDihedralsMax, diffDihedral);
				numPropers++;
			}
		}

		for (const auto& torsion : topologies[topoIx].getImproperHarmonicTorsions()) {
			const auto& ids = torsion.compoundAtomIndices;
		
			// Calculate computed dihedral
			SimTK::Real computedDihedral = bDihedral(
				newAtomtargets[topoIx][ids[0]], 
				newAtomtargets[topoIx][ids[1]], 
				newAtomtargets[topoIx][ids[2]], 
				newAtomtargets[topoIx][ids[3]]
			);

			// Calculate target dihedral
			SimTK::Real targetDihedral = bDihedral(
				atomTargets[topoIx].at(ids[0]), 
				atomTargets[topoIx].at(ids[1]), 
				atomTargets[topoIx].at(ids[2]), 
				atomTargets[topoIx].at(ids[3])
			);

			const SimTK::Real diffDihedral = std::abs(wrapAngle(targetDihedral - computedDihedral));
			error.improperDihedrals += diffDihedral;
			error.improperDihedralsMax = std::max(error.improperDihedralsMax, diffDihedral);
			numImpropers++;
		}
	}

	if (numAtoms) error.cartesian /= numAtoms;
	if (numBonds) error.bonds /= numBonds;
	if (numAngles) error.angles /= numAngles;
	if (numPropers) error.properDihedrals /= numPropers;
	if (numImpropers) error.improperDihedrals /= numImpropers;

	return error;
}

void World::updateFramesFromTopologies()
{
	// Iterate molecules
	for (const auto& bond : rigidBodyAtomBonds) {

		const auto& topology = topologies[bond.topologyIndex];
		SimTK::MobilizedBody& childAtomMobod = matter->updMobilizedBody(bond.childMBIx);
		SimTK::MobilizedBody& parentAtomMobod = matter->updMobilizedBody(bond.parentMBIx);

		// Bound to Ground
		if(parentAtomMobod.isGround()){
			const SimTK::Transform& G_X_T = topology.getTopLevelTransform();
			const SimTK::Transform& T_X_base = topology.getTopTransform(SimTK::Compound::AtomIndex(0));

			// Create transforms for child default inboard frame (XPF) and default outboard frame (XBM)
			const SimTK::Transform XPF = G_X_T * T_X_base;
			const SimTK::Transform XBM = SimTK::Transform();

			childAtomMobod.setDefaultInboardFrame(XPF);
			childAtomMobod.setDefaultOutboardFrame(XBM);

			continue;
		}

		// Get parent-child BondCenters relationship
		const SimTK::Transform X_parentBC_childBC = topology.getDefaultBondCenterFrameInOtherBondCenterFrame(bond.childCAIx, bond.parentCAIx);
		const SimTK::Transform X_childBC_parentBC = ~X_parentBC_childBC;

		// Get Top frame
		const SimTK::Transform& T_X_root = topology.getTopTransform(bond.childCAIx);

		// Origin of the parent mobod
		const SimTK::Transform& T_X_Proot = topology.getTopTransform(bond.parentMobodRootCAIx);
		const SimTK::Transform Proot_X_T = ~T_X_Proot;
		const SimTK::Transform Proot_X_root = Proot_X_T * T_X_root;

		// Create transforms for child default inboard frame (XPF) and default outboard frame (XBM)
		SimTK::Transform XPF, XBM;
		switch (bond.mobility) {
			case SimTK::BondMobility::Mobility::AnglePin:
			case SimTK::BondMobility::Mobility::Slider:
			case SimTK::BondMobility::Mobility::BendStretch: {
				const SimTK::Transform& B_X_M_anglePin = X_parentBC_childBC;
				const SimTK::Transform P_X_F_anglePin = Proot_X_root * B_X_M_anglePin;

				XPF = P_X_F_anglePin;
				XBM = B_X_M_anglePin;
				break;
			}

			case SimTK::BondMobility::Mobility::Torsion:
			case SimTK::BondMobility::Mobility::Cylinder: {
				const SimTK::Transform B_X_M_pin = X_parentBC_childBC * X_to_Z;
				const SimTK::Transform P_X_F_pin = Proot_X_root * B_X_M_pin;

				XPF = P_X_F_pin;
				XBM = B_X_M_pin;
				break;
			}

			case SimTK::BondMobility::Mobility::BallM:
			case SimTK::BondMobility::Mobility::Rigid:
			case SimTK::BondMobility::Mobility::Translation: {
				const SimTK::Transform& B_X_M = X_to_Z; // Samuel Flores' terminology aka M_X_pin = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis)
				const SimTK::Transform P_X_F = Proot_X_root * B_X_M;
				
				XPF = P_X_F;
				XBM = B_X_M;
				break;
			}

			case SimTK::BondMobility::Mobility::Spherical: {
				const SimTK::Transform P_X_F_spheric = SimTK::Transform() * Proot_X_root  * X_parentBC_childBC * X_to_Y * Y_to_Z; 
				const SimTK::Transform M_X_B_spheric = SimTK::Transform() * Z_to_Y * Y_to_X * X_childBC_parentBC;
				const SimTK::Transform B_X_M_spheric = ~M_X_B_spheric;

				XPF = P_X_F_spheric;
				XBM = B_X_M_spheric;
				break;
			}

			case SimTK::BondMobility::Mobility::OrthoSpherical: {
				// Get parent-child BC transform
				const SimTK::Transform X_parentAtom_BCpar = topology.calcDefaultBondCenterFrameInParentAtomFrame(bond.parentCAIx, bond.childCAIx);
				const SimTK::Transform X_childAtom_BCchi = topology.calcDefaultBondCenterFrameInChildAtomFrame(bond.parentCAIx, bond.childCAIx);
				const SimTK::Transform X_BCchi_childAtom = ~X_childAtom_BCchi;
			
				const SimTK::Transform P_X_F_orthospheric = X_parentAtom_BCpar; // BAT from Compound
				const SimTK::Transform M_X_B_orthospheric = X_parentBC_childBC * X_BCchi_childAtom; // BAT from Compound
				const SimTK::Transform B_X_M_orthospheric = ~M_X_B_orthospheric; // X_childAtom_BC * X_childBC_parentBC;
			
				XPF = P_X_F_orthospheric;
				XBM = B_X_M_orthospheric;
				break;
			}

			default:
				std::string errorMsg = "Mobility type '" + std::string(getBondMobilityName(bond.mobility)) + "' not supported in World::setFramesFromTopologies().";
				throw std::runtime_error(errorMsg);
		}

		childAtomMobod.setDefaultInboardFrame(XPF);
		childAtomMobod.setDefaultOutboardFrame(XBM);
	}

	// Handle bonds involving root atoms separately
	for (const auto& bond : rootAtomBonds) {
		const auto& topology = topologies[bond.topologyIndex];
		// SimTK_ASSERT_ALWAYS(bond.topologyIndex == 0, "Root atom bonds only supported for a single molecule. See bond.parentCAIx referencing below fo understand what i mean.");

		const SimTK::Transform& G_X_T = topology.getTopLevelTransform();
		
		if(topology.getAtoms()[bond.parentCAIx].connectivity.root){
			const SimTK::Transform& T_X_base = topology.getTopTransform(bond.parentCAIx);
			const SimTK::Transform G_X_base = G_X_T * T_X_base;

			SimTK::MobilizedBody& parentAtomMobod = matter->updMobilizedBody(bond.parentMBIx);
			parentAtomMobod.setDefaultInboardFrame(G_X_base);
			parentAtomMobod.setDefaultOutboardFrame(SimTK::Transform());
		}
		
		if(topology.getAtoms()[bond.childCAIx].connectivity.root){
			const SimTK::Transform& T_X_base = topology.getTopTransform(bond.childCAIx);
			const SimTK::Transform G_X_base = G_X_T * T_X_base;

			SimTK::MobilizedBody& childAtomMobod = matter->updMobilizedBody(bond.childMBIx);
			childAtomMobod.setDefaultInboardFrame(G_X_base);
			childAtomMobod.setDefaultOutboardFrame(SimTK::Transform());				
		}
	}
}

// TODO write a pdb writer for all the Compounds
// TODO move this in Topology since they work only for one Compound
void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  SimTK::Real mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  std::string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  SimTK::Real mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  std::string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, SimTK::Real aTime)
{
  SimTK::Real mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  std::string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, SimTK::Real aTime)
{
  SimTK::Real mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  std::string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::PdbStructure pdb, const char *FN)
{
  std::filebuf fb;
  fb.open(FN, std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); //automatically multiplies by ten (nm to A)
  fb.close();
}

void World::generateDummParams(
	const std::vector<RoboAtom>& atoms,
	const std::vector<RoboBond>& bonds,
	const std::vector<RoboAngle>& angles,
	const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
	const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions)
{
	// Make a counter that checks if the atom class index already exists since Molmodel does not check for and does not allow re-definitions
	std::vector<bool> atomClassDefined(atoms.size(), false);
	std::vector<bool> chargedAtomTypeDefined(atoms.size(), false);

	// Define atom classes and charged atom types
	for (auto& atom : atoms) {
		if (!atomClassDefined[atom.identity.atomClassIndex]) {
			atomClassDefined[atom.identity.atomClassIndex] = true;

			forceField->defineAtomClass(
				atom.identity.atomClassIndex,
				atom.identity.atomClassName.c_str(),
				atom.elementInfo.atomicNumber,
				atom.connectivity.numBondsInvolved,
				atom.physics.vdwRadiusInNm,
				atom.physics.vdwWellDepthInKJ
			);
		}

		if (!chargedAtomTypeDefined[atom.identity.chargedAtomTypeIndex]) {
			chargedAtomTypeDefined[atom.identity.chargedAtomTypeIndex] = true;

			forceField->defineChargedAtomType(
				atom.identity.chargedAtomTypeIndex,
				atom.identity.chargedAtomTypeName.c_str(),
				atom.identity.atomClassIndex,
				atom.physics.chargeInE
			);
			forceField->setBiotypeChargedAtomType(atom.identity.chargedAtomTypeIndex, atom.biotypeIndex);
		}
	}

	// Define bonds parameters using atom classes
	// DuMM canonicalizes the order of atom class indices internally (smallest first)
	// DuMM also checks for re-definitions and ignores them if found
	// However, it will fail if we try to define a bond with the same atom classes but different parameters
	std::map<BondStretchKey, BondStretchValue> definedBondStretches;

	for (auto& bond : bonds) {
		const auto aCIx1 = atoms[bond.globalIndices[0]].identity.atomClassIndex;
		const auto aCIx2 = atoms[bond.globalIndices[1]].identity.atomClassIndex;

		const auto key = BondStretchKey(aCIx1, aCIx2);
		const BondStretchValue newParams{
			{ atoms[bond.globalIndices[0]].identity.globalIndex, atoms[bond.globalIndices[1]].identity.globalIndex },
			bond.stiffnessInKJPerNmSq,
			bond.nominalLengthInNm
		};
		
		const auto it = definedBondStretches.find(key);
		if (it != definedBondStretches.end()) {
			const BondStretchValue& existingParams = it->second;
			if (existingParams != newParams) {
				std::string message = "Error in World::generateDummParams(): Conflicting bond stretch parameters.\n";
				message += "Existing paramters defined for atoms:\n";
				for (const auto i : existingParams.globalAtomIndices)
					message += "  " + getAtomDescription(atoms[i]) + "\n";
				message += "  Stiffness: " + std::to_string(existingParams.stiffness) + " kJ/(nm^2 mol), Length: " + std::to_string(existingParams.length) + " nm\n";
				message += "New parameters defined for atoms:\n";
				for (const auto i : newParams.globalAtomIndices)
					message += "  " + getAtomDescription(atoms[i]) + "\n";
				message += "  Stiffness: " + std::to_string(newParams.stiffness) + " kJ/(nm^2 mol), Length: " + std::to_string(newParams.length) + " nm\n";
				SimTK_ASSERT_ALWAYS(existingParams == newParams, message.c_str());
			}
		} else {
			definedBondStretches.insert(it, std::make_pair(key, newParams));
		}

		forceField->defineBondStretch(aCIx1, aCIx2, bond.stiffnessInKJPerNmSq, bond.nominalLengthInNm);
	}

	// Define angles
	// DuMM canonicalizes the order of atom class indices internally (smallest first)
	// DuMM also checks for re-definitions and ignores them if found
	// However, it will fail if we try to define an angle with the same atom classes but different parameters
	std::map<BondBendKey, BondBendValue> definedBondBends;

	for (const auto& angle : angles) {
		const auto aCIx1 = atoms[angle.globalIndices[0]].identity.atomClassIndex;
		const auto aCIx2 = atoms[angle.globalIndices[1]].identity.atomClassIndex;
		const auto aCIx3 = atoms[angle.globalIndices[2]].identity.atomClassIndex;

		const BondBendKey key(aCIx1, aCIx2, aCIx3);
		const BondBendValue newParams{
			{ atoms[angle.globalIndices[0]].identity.globalIndex, atoms[angle.globalIndices[1]].identity.globalIndex, atoms[angle.globalIndices[2]].identity.globalIndex },
			angle.stiffnessInKJPerRadSq,
			angle.nominalAngleInDeg
		};

		const auto it = definedBondBends.find(key);
		if (it != definedBondBends.end()) {
			const BondBendValue& existingParams = it->second;
			if (existingParams != newParams) {
				std::string message = "Error in World::generateDummParams(): Conflicting bond bend parameters.\n";
				message += "Existing paramters defined for atoms:\n";
				for (const auto i : existingParams.globalAtomIndices)
					message += "  " + getAtomDescription(atoms[i]) + "\n";
				message += "  Stiffness: " + std::to_string(existingParams.stiffness) + " kJ/(rad^2 mol), Angle: " + std::to_string(existingParams.angleDeg) + " deg\n";
				message += "New parameters defined for atoms:\n";
				for (const auto i : newParams.globalAtomIndices)
					message += "  " + getAtomDescription(atoms[i]) + "\n";
				message += "  Stiffness: " + std::to_string(newParams.stiffness) + " kJ/(rad^2 mol), Angle: " + std::to_string(newParams.angleDeg) + " deg\n";
				SimTK_ASSERT_ALWAYS(existingParams == newParams, message.c_str());
			}
		} else {
			definedBondBends.insert(it, std::make_pair(key, newParams));
		}

		forceField->defineBondBend(
			aCIx1, aCIx2, aCIx3,
			angle.stiffnessInKJPerRadSq,
			angle.nominalAngleInDeg
		);
	}

	// Define AMBER-style proper periodic torsions
	std::map<PeriodicTorsionKey, PeriodicTorsionValue> definedPeriodicProperTorsions;
	std::map<PeriodicTorsionKey, PeriodicTorsionValue> definedPeriodicImproperTorsions;

	for (const auto& torsion : properPeriodicTorsions) {
		const auto aCIx1 = atoms[torsion.globalIndices[0]].identity.atomClassIndex;
		const auto aCIx2 = atoms[torsion.globalIndices[1]].identity.atomClassIndex;
		const auto aCIx3 = atoms[torsion.globalIndices[2]].identity.atomClassIndex;
		const auto aCIx4 = atoms[torsion.globalIndices[3]].identity.atomClassIndex;

		// Dont canonicalize improper torsions
		const PeriodicTorsionKey key (aCIx1, aCIx2, aCIx3, aCIx4, !torsion.improper);
		const PeriodicTorsionValue newParams {
			{ atoms[torsion.globalIndices[0]].identity.globalIndex, atoms[torsion.globalIndices[1]].identity.globalIndex, atoms[torsion.globalIndices[2]].identity.globalIndex, atoms[torsion.globalIndices[3]].identity.globalIndex },
			{ torsion.terms[0].periodicity, torsion.terms[1].periodicity, torsion.terms[2].periodicity, torsion.terms[3].periodicity, torsion.terms[4].periodicity },
			{ torsion.terms[0].amplitudeKJ, torsion.terms[1].amplitudeKJ, torsion.terms[2].amplitudeKJ, torsion.terms[3].amplitudeKJ, torsion.terms[4].amplitudeKJ },
			{ torsion.terms[0].phaseDeg, torsion.terms[1].phaseDeg, torsion.terms[2].phaseDeg, torsion.terms[3].phaseDeg, torsion.terms[4].phaseDeg },
			torsion.numTerms
		};

		bool found = false;
		std::map<PeriodicTorsionKey, PeriodicTorsionValue>::iterator it;
		if (torsion.improper) {
			it = definedPeriodicImproperTorsions.find(key);
			found = definedPeriodicImproperTorsions.end() != it;
		} else {
			it = definedPeriodicProperTorsions.find(key);
			found = definedPeriodicProperTorsions.end() != it;
		}

		if (found) {
			const PeriodicTorsionValue& existingParams = it->second;
			if (existingParams != newParams) {
				std::string message = "Error in World::generateDummParams(): Conflicting " + std::string(torsion.improper ? "improper" : "proper") + " periodic torsion parameters for atom classes ";
				message += "Existing parameters defined for atoms:\n";
				for (const auto i : existingParams.globalAtomIndices) {
					message += "\t\t" + getAtomDescription(atoms[i]) + "\n";
				}
				message += "Existing parameters:\n";
				for (size_t i = 0; i < existingParams.numTerms; i++) {
					message += "\t\tTerm " + std::to_string(i) + ": periodicity=" + std::to_string(existingParams.periodicity[i]) +
							   ", amplitude=" + std::to_string(existingParams.amplitude[i]) +
							   ", phase=" + std::to_string(existingParams.phase[i]) + "\n";
				}

				message += "New parameters defined for atoms:\n";
				for (const auto i : newParams.globalAtomIndices) {
					message += "\t\t" + getAtomDescription(atoms[i]) + "\n";
				}
				message += "New parameters:\n";
				for (size_t i = 0; i < newParams.numTerms; i++) {
					message += "\t\tTerm " + std::to_string(i) + ": periodicity=" + std::to_string(newParams.periodicity[i]) +
							   ", amplitude=" + std::to_string(newParams.amplitude[i]) +
							   ", phase=" + std::to_string(newParams.phase[i]) + "\n";
				}
				SimTK_ASSERT_ALWAYS(existingParams == newParams, message.c_str());
			}
		} else {
			if (torsion.improper) {
				definedPeriodicImproperTorsions.insert(it, std::make_pair(key, newParams));
			} else {
				definedPeriodicProperTorsions.insert(it, std::make_pair(key, newParams));
			}
		}

		if (torsion.improper) {
			switch (torsion.numTerms)
			{
			case 1:
				forceField->defineAmberImproperTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg);
				break;
			case 2:
				forceField->defineAmberImproperTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg);
				break;
			case 3:
				forceField->defineAmberImproperTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg,
					torsion.terms[2].periodicity, torsion.terms[2].amplitudeKJ, torsion.terms[2].phaseDeg);
				break;

			default:
				SimTK_ASSERT_ALWAYS(false, "Error in World::generateDummParams: Improper torsions can only have 1, 2, or 3 terms.");
				break;
			}
		} else {
			switch (torsion.numTerms)
			{
			case 1:
				forceField->defineBondTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg);
				break;
			case 2:
				forceField->defineBondTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg);
				break;
			case 3:
				forceField->defineBondTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg,
					torsion.terms[2].periodicity, torsion.terms[2].amplitudeKJ, torsion.terms[2].phaseDeg);
				break;
			case 4:
				forceField->defineBondTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg,
					torsion.terms[2].periodicity, torsion.terms[2].amplitudeKJ, torsion.terms[2].phaseDeg,
					torsion.terms[3].periodicity, torsion.terms[3].amplitudeKJ, torsion.terms[3].phaseDeg);
				break;
			case 5:
				forceField->defineBondTorsion(
					aCIx1, aCIx2, aCIx3, aCIx4,
					torsion.terms[0].periodicity, torsion.terms[0].amplitudeKJ, torsion.terms[0].phaseDeg,
					torsion.terms[1].periodicity, torsion.terms[1].amplitudeKJ, torsion.terms[1].phaseDeg,
					torsion.terms[2].periodicity, torsion.terms[2].amplitudeKJ, torsion.terms[2].phaseDeg,
					torsion.terms[3].periodicity, torsion.terms[3].amplitudeKJ, torsion.terms[3].phaseDeg,
					torsion.terms[4].periodicity, torsion.terms[4].amplitudeKJ, torsion.terms[4].phaseDeg);
				break;
			
			default:
				SimTK_ASSERT_ALWAYS(false, "Error in World::generateDummParams: Proper torsions can only have 1 to 5 terms.");
				break;
			}
			
		}
	}

	// Define CHARMM-style improper harmonic torsions
	for (const auto& torsion : harmonicImproperTorsions) {
		const auto aCIx1 = atoms[torsion.globalIndices[1]].identity.atomClassIndex;
		const auto aCIx2 = atoms[torsion.globalIndices[2]].identity.atomClassIndex;
		const auto aCIx3 = atoms[torsion.globalIndices[3]].identity.atomClassIndex;
		const auto aCIx4 = atoms[torsion.globalIndices[4]].identity.atomClassIndex;

		forceField->defineCustomBondTorsion(aCIx1, aCIx2, aCIx3, aCIx4, new HarmonicImproperTorsionForce(torsion.stiffnessInKJPerRadSq, torsion.nominalAngleInRad));
	}
}

bool World::isOverconstrained() const {

	std::cout << "[INFO] Checking for overconstraints in world " << ownWorldIndex << ":" << std::endl;

	// Get the advanced state from the integrator to check for overconstraints
	const auto& state = integrator->getAdvancedState();
	
	// Realize will calculate the forces and potential energy using OpenMM
	// OPENMM::get().setActiveForceGroup(world.getOwnIndex());
	compoundSystem->realize(state, SimTK::Stage::Acceleration);

	// Structure to hold mobilized bodies info
	struct MobilizedBodyInfo {
		std::unordered_set<std::string> atoms;
		std::unordered_set<std::string> atomsInRingClosingBonds;
		std::vector<SimTK::Real> q;
		std::vector<SimTK::Real> u;
		std::vector<SimTK::Real> qdot;
		std::vector<SimTK::Real> udot;
		std::vector<SimTK::Real> qdotdot;
		bool invalid { false };
		bool hasRingClosingBond { false };
	};

	const int NU = compoundSystem->getMatterSubsystem().getNU(state);
	const int numBodies = compoundSystem->getMatterSubsystem().getNumBodies();
	std::vector<MobilizedBodyInfo> mobilizedBodiesInfo(numBodies);
	std::map<std::pair<SimTK::MobilizedBodyIndex, SimTK::MobilizedBodyIndex>, std::string> dihedralTypes;

	// Get atom names of each mobilized body
	for (const auto& topology : topologies) {
		for (const auto& bond : topology.getBonds()) {
			const auto mbx1 = topology.getAtomMobilizedBodyIndexThroughDumm(bond.compoundAtomIndices[0], *forceField);
			const auto mbx2 = topology.getAtomMobilizedBodyIndexThroughDumm(bond.compoundAtomIndices[1], *forceField);

			const auto& atom1_name = topology.getAtoms()[bond.compoundAtomIndices[0]].identity.uniqueAtomName;
			const auto& atom2_name = topology.getAtoms()[bond.compoundAtomIndices[1]].identity.uniqueAtomName;

			// Save dihedral type
			if (mbx1 != mbx2) {
				dihedralTypes[std::make_pair(mbx1, mbx2)] = bond.dihedralType;
				dihedralTypes[std::make_pair(mbx2, mbx1)] = bond.dihedralType;
			}

			// Save atoms involved in ring-closing bonds
			if (bond.ringClosing) {
				mobilizedBodiesInfo[mbx1].atomsInRingClosingBonds.insert(atom1_name);
				mobilizedBodiesInfo[mbx2].atomsInRingClosingBonds.insert(atom2_name);

				mobilizedBodiesInfo[mbx1].hasRingClosingBond = true;
				mobilizedBodiesInfo[mbx2].hasRingClosingBond = true;
			}
		}

		for (const auto& atom : topology.getAtoms()) {
			const auto mbIx = topology.getAtomMobilizedBodyIndexThroughDumm(atom.identity.compoundAtomIndex, *forceField);
        	const auto& name = atom.identity.uniqueAtomName;

			if (mobilizedBodiesInfo[mbIx].atomsInRingClosingBonds.find(name) == mobilizedBodiesInfo[mbIx].atomsInRingClosingBonds.end()) {
				mobilizedBodiesInfo[mbIx].atoms.insert(name);
			}
		}
	}

	auto isOutlier = [](SimTK::Real val, bool strict) {
		if (std::isnan(val) || std::isinf(val)) return true;
		return std::abs(val) > 1e3;
	};

	bool anyInvalid = false;
	for (SimTK::MobilizedBodyIndex mbIx(1); mbIx < numBodies; ++mbIx) {
		const auto& mobod = compoundSystem->getMatterSubsystem().getMobilizedBody(mbIx);
		const bool strict = (!mobilizedBodiesInfo[mbIx].hasRingClosingBond);
		bool invalid = false;

		for (const SimTK::Real q: mobod.getQAsVector(state)) {
			mobilizedBodiesInfo[mbIx].q.push_back(q);
			if (isOutlier(q, strict)) invalid = true;
		};
		for (const SimTK::Real u: mobod.getUAsVector(state)) {
			mobilizedBodiesInfo[mbIx].u.push_back(u);
			if (isOutlier(u, strict)) invalid = true;
		};
		for (const SimTK::Real qdot: mobod.getQDotAsVector(state)) {
			mobilizedBodiesInfo[mbIx].qdot.push_back(qdot);
			if (isOutlier(qdot, strict)) invalid = true;
		};
		for (const SimTK::Real udot: mobod.getUDotAsVector(state)) {
			mobilizedBodiesInfo[mbIx].udot.push_back(udot);
			if (isOutlier(udot, strict)) invalid = true;
		};
		for (const SimTK::Real qdotdot: mobod.getQDotDotAsVector(state)) {
			mobilizedBodiesInfo[mbIx].qdotdot.push_back(qdotdot);
			if (isOutlier(qdotdot, strict)) invalid = true;
		};

		mobilizedBodiesInfo[mbIx].invalid = invalid;
		anyInvalid = anyInvalid || invalid;
	}

	// Stop here if all mobilized bodies are valid
	// We don't want to print a huge table if everything is fine
	if (!anyInvalid) {
		std::cout << "\tNo overconstraints detected." << std::endl;
		return false;
	}

	// Print mobilized bodies info
	struct Row {
		std::string mbx, parent_mbx, dihedral, invalid, q, u, qdot, udot, qddot, atoms, ring_atoms;
	};

	// Initialize widths
	std::vector<Row> rows;
	size_t w_mbIx = 4, w_parent = 10, w_dihedral = 13, w_invalid = 7, 
		w_q = 1, w_u = 1, w_qdot = 1, w_udot = 1, w_qddot = 1, 
		w_atoms = 5, w_ring = 10;

	for (SimTK::MobilizedBodyIndex mbIx(1); mbIx < numBodies; ++mbIx) {
		const auto& mobod = matter->getMobilizedBody(mbIx);
		const auto& info = mobilizedBodiesInfo[mbIx];
		const auto parentIx = mobod.getParentMobilizedBody().getMobilizedBodyIndex();

		Row r;
		r.mbx        = std::to_string(mbIx);
		r.parent_mbx = std::to_string(parentIx);
		
		// Look up dihedral type between this body and its parent
		auto it = dihedralTypes.find(std::make_pair(mbIx, parentIx));
		r.dihedral   = (it != dihedralTypes.end()) ? it->second : "ground";

		r.invalid    = info.invalid ? "invalid" : "";
		r.q          = vecToString(info.q);
		r.u          = vecToString(info.u);
		r.qdot       = vecToString(info.qdot);
		r.udot       = vecToString(info.udot);
		r.qddot      = vecToString(info.qdotdot);
		r.atoms      = atomsToString(info.atoms);
		// New column for ring closing atoms
		r.ring_atoms = atomsToString(info.atomsInRingClosingBonds);

		// Update widths
		w_mbIx     = std::max(w_mbIx,     r.mbx.size());
		w_parent   = std::max(w_parent,   r.parent_mbx.size());
		w_dihedral = std::max(w_dihedral, r.dihedral.size());
		w_q        = std::max(w_q,        r.q.size());
		w_u        = std::max(w_u,        r.u.size());
		w_qdot     = std::max(w_qdot,     r.qdot.size());
		w_udot     = std::max(w_udot,     r.udot.size());
		w_qddot    = std::max(w_qddot,    r.qddot.size());
		w_atoms    = std::max(w_atoms,    r.atoms.size());
		w_ring     = std::max(w_ring,     r.ring_atoms.size());

		rows.push_back(std::move(r));
	}

	// Print Header
	std::cerr << std::left
			<< std::setw(w_mbIx + 2)     << "mbx"
			<< std::setw(w_parent + 2)   << "parent"
			<< std::setw(w_dihedral + 2) << "dihedral"
			<< std::setw(w_invalid + 2)  << "invalid"
			<< std::setw(w_q + 2)        << "q"
			<< std::setw(w_u + 2)        << "u"
			<< std::setw(w_qdot + 2)     << "qdot"
			<< std::setw(w_udot + 2)     << "udot"
			<< std::setw(w_qddot + 2)    << "qddot"
			<< std::setw(w_atoms + 2)    << "atoms"
			<< "ring_atoms"
			<< std::endl;

	// Recalculate total width for the separator line (10 columns with padding)
	size_t total_width = w_mbIx + w_parent + w_dihedral + w_invalid + w_q + w_u + 
						w_qdot + w_udot + w_qddot + w_atoms + w_ring + (2 * 10);

	std::cerr << std::string(total_width, '-') << std::endl;

	// Print Rows
	for (const auto& r : rows) {
		std::cerr << std::left
				<< std::setw(w_mbIx + 2)     << r.mbx
				<< std::setw(w_parent + 2)   << r.parent_mbx
				<< std::setw(w_dihedral + 2) << r.dihedral
				<< std::setw(w_invalid + 2)  << r.invalid
				<< std::setw(w_q + 2)        << r.q
				<< std::setw(w_u + 2)        << r.u
				<< std::setw(w_qdot + 2)     << r.qdot
				<< std::setw(w_udot + 2)     << r.udot
				<< std::setw(w_qddot + 2)    << r.qddot
				<< std::setw(w_atoms + 2)    << r.atoms
				<< r.ring_atoms
				<< std::endl;
	}

	return true;
}

bool World::hasRigidBodyViolations(SimTK::Real timeStep, int numSteps) {
	if (samplers.empty()) {
		throw std::runtime_error("World::hasRigidBodyViolations() requires at least one sampler to be defined.");
	}

	auto wrapAngle = [](SimTK::Real x) {
		return std::atan2(std::sin(x), std::cos(x));
	};

	std::stringstream nullStream;
	bool hasViolations = false;

	// The currently implemented tests check if rigid bonds, angles and torsions are correctly ignored
	// This is only supported for rigid or torsional bonds
	std::vector<RigidBond> rigidBonds;
	std::vector<RigidAngle> rigidAngles;
	std::vector<RigidTorsion> rigidProperTorsions, rigidImproperTorsions;

	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		const auto& topology = topologies[topoIx];

		// Add bonds (stretches)
		for (const auto& bond : topology.getBonds()) {
			const SimTK::Compound::AtomIndex childCAIx = bond.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex parentCAIx = bond.compoundAtomIndices[0];

			const SimTK::MobilizedBodyIndex mbxChild = topology.getAtomMobilizedBodyIndex(childCAIx);
			const SimTK::MobilizedBodyIndex mbxParent = topology.getAtomMobilizedBodyIndex(parentCAIx);

			// // If the two atoms are not in the same mobilized body, then this bond cannot be rigid
			// if (mbxChild != mbxParent) {
			// 	continue;
			// }

			RigidBond rigidBond;
			rigidBond.topologyIndex = topoIx;
			rigidBond.childCAIx = childCAIx;
			rigidBond.parentCAIx = parentCAIx;
			rigidBond.ringClosing = bond.ringClosing;
			rigidBonds.push_back(rigidBond);
		}

		// Add angles (bends)
		for (const auto& angle :topology.getAngles()) {
			const SimTK::Compound::AtomIndex aIx1 = angle.compoundAtomIndices[0];
			const SimTK::Compound::AtomIndex aIx2 = angle.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex aIx3 = angle.compoundAtomIndices[2];
			const SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(aIx1);
			const SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(aIx2);
			const SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(aIx3);

			// // Ignore angles that are not between two atoms in the same mobilized body
			// if (mbx1 != mbx2 || mbx2 != mbx3)
			// 	continue;

			RigidAngle rigidAngle;
			rigidAngle.topologyIndex = topoIx;
			rigidAngle.cAIx1 = aIx1;
			rigidAngle.cAIx2 = aIx2;
			rigidAngle.cAIx3 = aIx3;
			rigidAngle.ringClosing =
				topology.getBondByCompoundAtomIndex(aIx1, aIx2).ringClosing ||
				topology.getBondByCompoundAtomIndex(aIx2, aIx3).ringClosing;
			rigidAngles.push_back(rigidAngle);
		}

		// Process periodic dihedrals - proper and improper
		for (const auto& torsion :topology.getPeriodicTorsions()) {
			const SimTK::Compound::AtomIndex aIx1 = torsion.compoundAtomIndices[0];
			const SimTK::Compound::AtomIndex aIx2 = torsion.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex aIx3 = torsion.compoundAtomIndices[2];
			const SimTK::Compound::AtomIndex aIx4 = torsion.compoundAtomIndices[3];
			const SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(aIx1);
			const SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(aIx2);
			const SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(aIx3);
			const SimTK::MobilizedBodyIndex mbx4 = topology.getAtomMobilizedBodyIndex(aIx4);

			RigidTorsion rigidTorsion;
			rigidTorsion.topologyIndex = topoIx;
			rigidTorsion.cAIx1 = aIx1;
			rigidTorsion.cAIx2 = aIx2;
			rigidTorsion.cAIx3 = aIx3;
			rigidTorsion.cAIx4 = aIx4;

			if (torsion.improper) {
				if (mbx1 == mbx2 && mbx2 == mbx3 && mbx3 == mbx4) {
					rigidTorsion.ringClosing = 
						topology.getBondByCompoundAtomIndex(aIx1, aIx3).ringClosing ||
						topology.getBondByCompoundAtomIndex(aIx2, aIx3).ringClosing ||
						topology.getBondByCompoundAtomIndex(aIx3, aIx4).ringClosing;
					rigidImproperTorsions.push_back(rigidTorsion);
				}
			} else {
				// For proper torsions, the mobile axis is the 2-3 bond.
				if (mbx2 == mbx3) {
					rigidTorsion.ringClosing = 
						topology.getBondByCompoundAtomIndex(aIx1, aIx2).ringClosing ||
						topology.getBondByCompoundAtomIndex(aIx2, aIx3).ringClosing ||
						topology.getBondByCompoundAtomIndex(aIx3, aIx4).ringClosing;
					rigidProperTorsions.push_back(rigidTorsion);
				}
			}
		}

		// Process improper harmonic torsions
		for (const auto& torsion :topology.getImproperHarmonicTorsions()) {
			const SimTK::Compound::AtomIndex aIx1 = torsion.compoundAtomIndices[0];
			const SimTK::Compound::AtomIndex aIx2 = torsion.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex aIx3 = torsion.compoundAtomIndices[2];
			const SimTK::Compound::AtomIndex aIx4 = torsion.compoundAtomIndices[3];
			const SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndex(aIx1);
			const SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndex(aIx2);
			const SimTK::MobilizedBodyIndex mbx3 = topology.getAtomMobilizedBodyIndex(aIx3);
			const SimTK::MobilizedBodyIndex mbx4 = topology.getAtomMobilizedBodyIndex(aIx4);

			RigidTorsion rigidTorsion;
			rigidTorsion.topologyIndex = topoIx;
			rigidTorsion.cAIx1 = aIx1;
			rigidTorsion.cAIx2 = aIx2;
			rigidTorsion.cAIx3 = aIx3;
			rigidTorsion.cAIx4 = aIx4;
			rigidTorsion.ringClosing = 
				topology.getBondByCompoundAtomIndex(aIx1, aIx3).ringClosing ||
				topology.getBondByCompoundAtomIndex(aIx2, aIx3).ringClosing ||
				topology.getBondByCompoundAtomIndex(aIx3, aIx4).ringClosing;

			if (mbx1 == mbx2 && mbx2 == mbx3 && mbx3 == mbx4) {
				rigidImproperTorsions.push_back(rigidTorsion);
			}
		}
	}

	// Deep copy the initial state to avoid modifying it during sampling
	const auto initialState = integrator->getAdvancedState();
	SimTK::State& state = integrator->updAdvancedState();

	// Set up the sampler for testing
	setTemperature(300);
	setBoostTemperature(300);
	samplers[0]->setMDStepsPerSample(numSteps);
	samplers[0]->setBoostMDSteps(numSteps);
	samplers[0]->setTimestep(timeStep, false);
	samplers[0]->setAcceptRejectMode(AcceptRejectMode::AlwaysAccept);
	samplers[0]->reinitialize(state, nullStream, false);

	// Save initial atom target locations before sampling
	atomTargetLocaltionsCacheOld.resize(topologies.size());
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (SimTK::Compound::AtomIndex cAIx = SimTK::Compound::AtomIndex(0); cAIx < topologies[topoIx].getAtoms().size(); cAIx++) {
			const SimTK::Vec3 location = topologies[topoIx].calcAtomLocationInGroundFrameThroughSimbody(cAIx, *forceField, *matter, state);
			atomTargetLocaltionsCacheOld[topoIx][cAIx] = location;
		}
	}

	// Get one sample
	samplers[0]->sample_iteration(state, nullStream, false);

	// Update atom target locations cache after sampling
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (SimTK::Compound::AtomIndex cAIx = SimTK::Compound::AtomIndex(0); cAIx < topologies[topoIx].getAtoms().size(); cAIx++) {
			const SimTK::Vec3 location = topologies[topoIx].calcAtomLocationInGroundFrameThroughSimbody(cAIx, *forceField, *matter, state);
			atomTargetLocaltionsCache[topoIx][cAIx] = location;
		}
	}
	
	// Begin
	std::cout << "[INFO] Checking for rigid body violations in world " << getOwnIndex() << ":" << std::endl;

	// Not necessarily a rigid body check, but if the RMSD is very low while always accepting, then something is wrong
	SimTK::Real sumSquaredDist = 0.0;
	int numAtoms = 0;
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (SimTK::Compound::AtomIndex cAIx = SimTK::Compound::AtomIndex(0); cAIx < topologies[topoIx].getAtoms().size(); cAIx++) {
			const SimTK::Vec3& newLocation = atomTargetLocaltionsCache[topoIx].at(cAIx);
			const SimTK::Vec3& oldLocation = atomTargetLocaltionsCacheOld[topoIx].at(cAIx);
			sumSquaredDist += (newLocation - oldLocation).normSqr();
			numAtoms++;
		}
	}
	const SimTK::Real rmsd = numAtoms > 0 ? std::sqrt(sumSquaredDist / numAtoms) : 0.0;
	if (rmsd < 1e-6) {
		std::cerr << "\t[ERROR] RMSD after 1 step is " << rmsd << ", which is below the threshold." << std::endl;
		hasViolations = true;
	}

	// Bonds inside rigid bodies should not change lengths
	for (const auto& rigidBond : rigidBonds) {
		const auto topoIx = rigidBond.topologyIndex;

		const SimTK::Vec3& childOld = atomTargetLocaltionsCacheOld[topoIx].at(rigidBond.childCAIx);
		const SimTK::Vec3& parentOld = atomTargetLocaltionsCacheOld[topoIx].at(rigidBond.parentCAIx);
		const SimTK::Real oldDistance = (childOld - parentOld).norm();

		const SimTK::Vec3& childNew = atomTargetLocaltionsCache[topoIx].at(rigidBond.childCAIx);
		const SimTK::Vec3& parentNew = atomTargetLocaltionsCache[topoIx].at(rigidBond.parentCAIx);
		const SimTK::Real newDistance = (childNew - parentNew).norm();

		const SimTK::Real diff = std::abs(newDistance - oldDistance);

		const auto atom1Name = topologies[topoIx].getAtoms()[rigidBond.childCAIx].identity.uniqueAtomName;
		const auto atom2Name = topologies[topoIx].getAtoms()[rigidBond.parentCAIx].identity.uniqueAtomName;

		if (rigidBond.ringClosing) {
			// For ring-closing bonds, the constraint residual satisfies max_deviation <= C * sqrt(constraint_tolerance) where C is a geometry-dependent amplification factor.
			// In multibody systems, C reflects how generalized coordinate errors map to Cartesian bond errors.
			// For loop closures spanning multiple bodies, this amplification can be O(10 - 100).
			// Empirically, our systems exhibit C ≈ 50 (proline, disulfide bridges), consistent with the depth and conditioning of the kinematic chain.
			// This value reflects the expected conditioning of the constraint projection in generalized coordinates.
			if (diff > integrator->getConstraintToleranceInUse() * 50) {
				std::cerr << "\t[ERROR] Bond between " << atom1Name << " - " << atom2Name
					<< " changed length from " << oldDistance * 10 << " A to " << newDistance * 10 << " A (ring closing: yes)." << std::endl;
				hasViolations = true;
			}
		} else {
			if (diff > 1e-6) {
				std::cerr << "\t[ERROR] Bond between " << atom1Name << " - " << atom2Name
					<< " changed length from " << oldDistance * 10 << " A to " << newDistance * 10 << " A (ring closing: no)." << std::endl;
				hasViolations = true;
			}
		}
	}

	// Angles inside rigid bodies should not change values
	for (const auto& rigidAngle : rigidAngles) {
		const auto topoIx = rigidAngle.topologyIndex;

		const SimTK::Vec3& atom1Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidAngle.cAIx1);
		const SimTK::Vec3& atom2Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidAngle.cAIx2);
		const SimTK::Vec3& atom3Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidAngle.cAIx3);
		const SimTK::Real oldAngle = bAngle(atom1Old, atom2Old, atom3Old);

		const SimTK::Vec3& atom1New = atomTargetLocaltionsCache[topoIx].at(rigidAngle.cAIx1);
		const SimTK::Vec3& atom2New = atomTargetLocaltionsCache[topoIx].at(rigidAngle.cAIx2);
		const SimTK::Vec3& atom3New = atomTargetLocaltionsCache[topoIx].at(rigidAngle.cAIx3);
		const SimTK::Real newAngle = bAngle(atom1New, atom2New, atom3New);

		const SimTK::Real diff = std::abs(wrapAngle(newAngle - oldAngle)) * SimTK_RADIAN_TO_DEGREE;

		const auto atom1Name = topologies[topoIx].getAtoms()[rigidAngle.cAIx1].identity.uniqueAtomName;
		const auto atom2Name = topologies[topoIx].getAtoms()[rigidAngle.cAIx2].identity.uniqueAtomName;
		const auto atom3Name = topologies[topoIx].getAtoms()[rigidAngle.cAIx3].identity.uniqueAtomName;

		if (rigidAngle.ringClosing) {
			if (diff > 1) {
				std::cerr << "\t[ERROR] Angle between " << atom1Name << " - " << atom2Name << " - " << atom3Name
					<< " changed from " << oldAngle * SimTK_RADIAN_TO_DEGREE << " deg to " << newAngle * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: yes)." << std::endl;
				hasViolations = true;
			}
		} else {
			if (diff > 1e-6) {
				std::cerr << "\t[ERROR] Angle between " << atom1Name << " - " << atom2Name << " - " << atom3Name
					<< " changed from " << oldAngle * SimTK_RADIAN_TO_DEGREE << " deg to " << newAngle * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: no)." << std::endl;
				hasViolations = true;
			}
		}
	}

	// Proper and improper dihedrals inside rigid bodies should not change values
	for (const auto& rigidProperTorsion : rigidProperTorsions) {
		const auto topoIx = rigidProperTorsion.topologyIndex;

		const SimTK::Vec3& atom1Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidProperTorsion.cAIx1);
		const SimTK::Vec3& atom2Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidProperTorsion.cAIx2);
		const SimTK::Vec3& atom3Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidProperTorsion.cAIx3);
		const SimTK::Vec3& atom4Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidProperTorsion.cAIx4);
		const SimTK::Real oldDihedral = bDihedral(atom1Old, atom2Old, atom3Old, atom4Old);

		const SimTK::Vec3& atom1New = atomTargetLocaltionsCache[topoIx].at(rigidProperTorsion.cAIx1);
		const SimTK::Vec3& atom2New = atomTargetLocaltionsCache[topoIx].at(rigidProperTorsion.cAIx2);
		const SimTK::Vec3& atom3New = atomTargetLocaltionsCache[topoIx].at(rigidProperTorsion.cAIx3);
		const SimTK::Vec3& atom4New = atomTargetLocaltionsCache[topoIx].at(rigidProperTorsion.cAIx4);
		const SimTK::Real newDihedral = bDihedral(atom1New, atom2New, atom3New, atom4New);

		const SimTK::Real diff = std::abs(wrapAngle(newDihedral - oldDihedral)) * SimTK_RADIAN_TO_DEGREE;

		const auto atom1Name = topologies[topoIx].getAtoms()[rigidProperTorsion.cAIx1].identity.uniqueAtomName;
		const auto atom2Name = topologies[topoIx].getAtoms()[rigidProperTorsion.cAIx2].identity.uniqueAtomName;
		const auto atom3Name = topologies[topoIx].getAtoms()[rigidProperTorsion.cAIx3].identity.uniqueAtomName;
		const auto atom4Name = topologies[topoIx].getAtoms()[rigidProperTorsion.cAIx4].identity.uniqueAtomName;

		if (rigidProperTorsion.ringClosing) {
			if (diff > 1) {
				std::cerr << "\t[ERROR] Proper torsion between " << atom1Name << " - " << atom2Name << " - " << atom3Name << " - " << atom4Name
					<< " changed from " << oldDihedral * SimTK_RADIAN_TO_DEGREE << " deg to " << newDihedral * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: yes)." << std::endl;
				hasViolations = true;
			}
		} else {
			if (diff > 1e-6) {
				std::cerr << "\t[ERROR] Proper torsion between " << atom1Name << " - " << atom2Name << " - " << atom3Name << " - " << atom4Name
					<< " changed from " << oldDihedral * SimTK_RADIAN_TO_DEGREE << " deg to " << newDihedral * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: no)." << std::endl;
				hasViolations = true;
			}
		}
	}

	for (const auto& rigidImproperTorsion : rigidImproperTorsions) {
		const auto topoIx = rigidImproperTorsion.topologyIndex;

		const SimTK::Vec3& atom1Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidImproperTorsion.cAIx1);
		const SimTK::Vec3& atom2Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidImproperTorsion.cAIx2);
		const SimTK::Vec3& atom3Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidImproperTorsion.cAIx3);
		const SimTK::Vec3& atom4Old = atomTargetLocaltionsCacheOld[topoIx].at(rigidImproperTorsion.cAIx4);
		const SimTK::Real oldDihedral = bDihedral(atom1Old, atom2Old, atom3Old, atom4Old);

		const SimTK::Vec3& atom1New = atomTargetLocaltionsCache[topoIx].at(rigidImproperTorsion.cAIx1);
		const SimTK::Vec3& atom2New = atomTargetLocaltionsCache[topoIx].at(rigidImproperTorsion.cAIx2);
		const SimTK::Vec3& atom3New = atomTargetLocaltionsCache[topoIx].at(rigidImproperTorsion.cAIx3);
		const SimTK::Vec3& atom4New = atomTargetLocaltionsCache[topoIx].at(rigidImproperTorsion.cAIx4);
		const SimTK::Real newDihedral = bDihedral(atom1New, atom2New, atom3New, atom4New);

		const SimTK::Real diff = std::abs(wrapAngle(newDihedral - oldDihedral)) * SimTK_RADIAN_TO_DEGREE;
		
		const auto atom1Name = topologies[topoIx].getAtoms()[rigidImproperTorsion.cAIx1].identity.uniqueAtomName;
		const auto atom2Name = topologies[topoIx].getAtoms()[rigidImproperTorsion.cAIx2].identity.uniqueAtomName;
		const auto atom3Name = topologies[topoIx].getAtoms()[rigidImproperTorsion.cAIx3].identity.uniqueAtomName;
		const auto atom4Name = topologies[topoIx].getAtoms()[rigidImproperTorsion.cAIx4].identity.uniqueAtomName;

		if (rigidImproperTorsion.ringClosing) {
			if (diff > 1) {
				std::cerr << "\t[ERROR] Improper torsion between " << atom1Name << " - " << atom2Name << " - " << atom3Name << " - " << atom4Name
					<< " changed from " << oldDihedral * SimTK_RADIAN_TO_DEGREE << " deg to " << newDihedral * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: yes)." << std::endl;
				hasViolations = true;
			}
		} else {
			if (diff > 1e-6) {
				std::cerr << "\t[ERROR] Improper torsion between " << atom1Name << " - " << atom2Name << " - " << atom3Name << " - " << atom4Name
					<< " changed from " << oldDihedral * SimTK_RADIAN_TO_DEGREE << " deg to " << newDihedral * SimTK_RADIAN_TO_DEGREE << " deg (ring closing: no)." << std::endl;
				hasViolations = true;
			}
		}
	}

	// Reset state
	integrator->updAdvancedState() = initialState;

	if (!hasViolations) {
		std::cout << "\tNo rigid body violations detected." << std::endl;
	}

	return hasViolations;
}


/** Print a Compound Cartesian coordinates as given by
	 * Compound::calcAtomLocationInGroundFrame **/
void World::printPoss(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vertex;
	std::cout<<"Positions:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
			//c.calcAtomLocationInGroundFrameThroughSimbody(
			//	aIx, forceField, matter, advanced);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
	}
}

void World::printVels(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vel;
	std::cout<<"Velocities:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
}

void World::printPossVels(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vertex, vel;
	std::cout<<"Positions:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Velocities:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
}


//==============================================================================
//                   CLASS TaskSpace
//==============================================================================
/**
 *  Contains a Symbody task space and additional data
 **/
StationTaskLaurentiu::StationTaskLaurentiu(void)
{

}

//==============================================================================
//                             1. STRUCTURAL FUNCTIONS
//==============================================================================

/** Constructor. Initializes the following objects:
 *  - CompoundSystem,
 *	  - SimbodyMatterSubsystem, GeneralForceSubsystem, DecorationSubsystem,
 *		Visualizer, Visualizer::Reporter, DuMMForceFieldSubsystem,
 *  - Integrator with a TimeStepper on top **/
World::World(int worldIndex, Span<Topology> topo, bool testing, const ZMatrix& _zMatrix, bool isVisual, SimTK::Real visualizerFrequency) :
	ownWorldIndex(worldIndex),
	topologies(topo),
	testing(testing),
	zMatrix(_zMatrix)
{
	compoundSystem = std::make_unique<SimTK::CompoundSystem>();
	matter = std::make_unique<SimTK::SimbodyMatterSubsystem>(*compoundSystem);
	forces = std::make_unique<SimTK::GeneralForceSubsystem>(*compoundSystem);
	forceField = std::make_unique<SimTK::DuMMForceFieldSubsystem>(*compoundSystem);
	integrator = std::make_unique<SimTK::VerletIntegrator>(*compoundSystem);
	timeStepper = std::make_unique<SimTK::TimeStepper>(*compoundSystem, *integrator);

	// Thermodynamics
	this->temperature = -1; // this leads to unusal behaviour hopefully

	// Allocate atom target locations cache
	for (const auto& t : topologies) {
		atomTargetLocaltionsCache.emplace_back();
		numMolecules++;
		for (const auto& a : t.getAtoms()) {
			atomTargetLocaltionsCache.back().emplace(std::make_pair(a.identity.compoundAtomIndex, a.position));
			numAtoms++;
		}
	}

	// Contact system
	/*
	tracker = std::make_unique<ContactTrackerSubsystem>(compoundSystem);
	contactForces = std::make_unique<CompliantContactSubsystem>(compoundSystem, *tracker);
	contactForces->setTrackDissipatedEnergy(true);
	contactForces->setTransitionVelocity(1e-2);
    	clique1 = ContactSurface::createNewContactClique();
	*/


	// const SimTK::MultibodySystem& mbs = forces->getMul tibodySystem();
	
	// // Set the visual flag and if true initialize a Decorations Subsystem,
	// // a Visualizer and a Simbody EventReporter which interacts with the
	// // Visualizer
	// this->visual = isVisual;
	// if(visual){

	// 	decorations = std::make_unique<SimTK::DecorationSubsystem>(compoundSystem);
	// 	visualizer = std::make_unique<SimTK::Visualizer>(compoundSystem);
	// 	visualizerReporter = std::make_unique<SimTK::Visualizer::Reporter>(
	// 		*visualizer, std::abs(visualizerFrequency));

	// 	compoundSystem->addEventReporter(visualizerReporter.get());

	// 	if(contactForces){
	// 		std::cout << "[WARNING] Victor check Teodor's contacts." << std::endl;
	// 		visualizer->addDecorationGenerator(
	// 			new ForceArrowGenerator(mbs, *contactForces));
	// 	}else{
	// 		std::cout << "[WARNING] Teodor's contacts." << std::endl;
	// 	}
		

	// 	// Initialize a DecorationGenerator
	// 	paraMolecularDecorator = std::make_unique<ParaMolecularDecorator>(
	// 		compoundSystem->get(),
	// 		matter->get(),
	// 		forceField->get(),
	// 		forces->get()
	// 	);

	// 	visualizer->addDecorationGenerator(paraMolecularDecorator.get());
	// }

	// currStage = this->getSimbodyMatterSubsystem()->getStage(
	// 	this->getSimbodyMatterSubsystem()->getSystem().getDefaultState()
	// );

}

void World::modelTopologies() {
	for (std::size_t topoIx = 0; topoIx < this->topologies.size(); topoIx++) {
		Topology& topology = topologies[topoIx];

		compoundSystem->adoptCompound(topology);
		compoundSystem->modelOneCompound(SimTK::CompoundSystem::CompoundIndex(topoIx), topology.updAtomFrameCache(), "Rigid");

		// Iterate through atoms and get their MobilizedBodyIndeces
		for (const auto& atom : topology.getAtoms()) {
			SimTK::Compound::AtomIndex aIx = atom.identity.compoundAtomIndex;
			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
			
			// Map mbx2aIx contains only atoms at the origin of mobods
			if (topology.getAtomLocationInMobilizedBodyFrame(aIx) == 0) {
				std::pair<int, SimTK::Compound::AtomIndex> topoAtomPair(topoIx, aIx);
				mbx2aIx.insert(std::make_pair(mbx, topoAtomPair));
			}

			// TODO desk_mass_related - i think this is already done in CompoundSystem::modelOneCompound()
			SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(aIx);
			forceField->setDuMMAtomMass(dAIx, atom.physics.massInDaltons);
		}

		// Cache two types of special bonds: bonds invovling the root atoms and bonds between two rigid bodies
		for (const auto& bond : topology.getBonds()) {
			const SimTK::Compound::AtomIndex childCAIx = bond.compoundAtomIndices[1];
			const SimTK::Compound::AtomIndex parentCAIx = bond.compoundAtomIndices[0];
			const SimTK::MobilizedBodyIndex childMBIx = topology.getAtomMobilizedBodyIndex(childCAIx);
			const SimTK::MobilizedBodyIndex parentMBIx = topology.getAtomMobilizedBodyIndex(parentCAIx);

			// Check if one of the atoms is a root atom
			// This kind of bond can have any mobility or be ring closing
			if (topology.getAtoms()[bond.compoundAtomIndices[1]].connectivity.root || topology.getAtoms()[bond.compoundAtomIndices[0]].connectivity.root) {
				RootAtomBond rootAtomBond;
				rootAtomBond.topologyIndex = topoIx;
				rootAtomBond.childCAIx = childCAIx;
				rootAtomBond.parentCAIx = parentCAIx;
				rootAtomBond.childMBIx = childMBIx;
				rootAtomBond.parentMBIx = parentMBIx;

				rootAtomBonds.push_back(rootAtomBond);
			}

			// Now we start looking for bonds between two rigid bodies
			// TODO can ring closing bonds be between to rigid bodies? i think not
			if (bond.ringClosing) continue;

			// We only want flexible bonds
			if (bond.getBondMobility(ownWorldIndex) == SimTK::BondMobility::Mobility::Rigid) continue;

			// We have a bond between two rigid bodies
			RigidBodyAtomBond rigidBodyAtomBond;
			rigidBodyAtomBond.topologyIndex = topoIx;
			rigidBodyAtomBond.mobility = bond.getBondMobility(ownWorldIndex);

			// Save global atom indices
			rigidBodyAtomBond.childAtomGlobalIndex = bond.globalIndices[1];
			rigidBodyAtomBond.parentAtomGlobalIndex = bond.globalIndices[0];
	
			// Save compound atom indices. Note that they are local to each topology
			rigidBodyAtomBond.childCAIx = childCAIx;
			rigidBodyAtomBond.parentCAIx = parentCAIx;

			// Save mobilized body indices
			rigidBodyAtomBond.childMBIx = childMBIx;
			rigidBodyAtomBond.parentMBIx = parentMBIx;

			// Save compound atom index of the root atom in the parent rigid body
			rigidBodyAtomBond.parentMobodRootCAIx = mbx2aIx.at(rigidBodyAtomBond.parentMBIx).second;

			// Check if we have grand parents. Three atoms are needed to define angles
			// If parent atom is 0, then no grand parent
			if (rigidBodyAtomBond.parentCAIx > 0) {
				rigidBodyAtomBond.grandParentCAIx = topology.getInboardAtomIndex(rigidBodyAtomBond.parentCAIx);
				rigidBodyAtomBond.grandParentMBIx = topology.getAtomMobilizedBodyIndex(rigidBodyAtomBond.grandParentCAIx);
				// rigidBodyAtomBond.grandParentAtomGlobalIndex = topology.getDuMMAtomIndex(rigidBodyAtomBond.grandParentCAIx);
			}

			rigidBodyAtomBonds.push_back(rigidBodyAtomBond);
		}
	}

	compoundSystem->realizeTopology();
}

// Print recommended timesteps. We need and advanced State here
SimTK::Real World::getRecommendedTimesteps(void)
{
	SimTK::State& someState = integrator->updAdvancedState();
	int nu = matter->getNU(someState);

	SimTK::Real minTimeStep;
    SimTK::Real prevMinTimeStep = SimTK::Infinity;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        minTimeStep = 0.0007 * std::sqrt(mobod.getBodyMass(someState));
        if(minTimeStep < prevMinTimeStep){
            prevMinTimeStep = minTimeStep;
        }
    }

	return prevMinTimeStep;
}

//==============================================================================
//                   TaskSpace Functions
//==============================================================================

/**
 * Allocate memory for a task space consisting of a set of body indeces
 * station on the bodies expresed in both guest (target) and host 
 * and the difference between them
*/
void World::addTaskSpaceLS(void)
{
	// //StationTaskLaurentiu stationTask;

	// int guestTopology = 1;
	// std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	// int topi = -1;
	// for(auto& topology : topologies){
	// 	topi++;

	// 	if(topi == guestTopology){

	// 		// Guest atoms iteration
	// 		for (int bAtomIx : bAtomIxs_guest) {
	// 			SimTK::Compound::AtomIndex aIx = (topology.subAtomList[bAtomIx]).identity.compoundAtomIndex;
	// 			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, forceField);

	// 			onBodyB.emplace_back(mbx);
	// 			taskStationPInGuest.emplace_back(SimTK::Vec3());
	// 			taskStationPInHost.emplace_back(SimTK::Vec3());
	// 			taskDeltaStationP.emplace_back(SimTK::Vec3());

	// 			if(this->visual == true){
	// 				paraMolecularDecorator->loadArrow(SimTK::Vec3(0), SimTK::Vec3(0));
	// 			}

	// 		}
	// 	}
	// }
}

/**
 * Update target task space
*/
void World::updateTaskSpace(const SimTK::State& someState)
{

	// int hostTopology = 0;
	// int guestTopology = 1;

	// std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	// std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	// // Get stations in host
	// int topi = -1;
	// for(auto& topology : topologies){
	// 	topi++;
	// 	if(topi == hostTopology){

	// 		// Atoms
	// 		int tz = -1;
	// 		for (int bAtomIx : bAtomIxs_host) {
	// 			tz++;
	// 			SimTK::Compound::AtomIndex aIx = (topology.subAtomList[bAtomIx]).identity.compoundAtomIndex;
	// 			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, forceField);
	// 			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

	// 			SimTK::Transform X_GB = mobod.getBodyTransform(someState);
	// 			SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
	// 			taskStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);
	// 		}
	// 	}
	// }

	// // Get stations in guest
	// topi = -1;
	// for(auto& topology : topologies){
	// 	topi++;

	// 	if(topi == guestTopology){

	// 		// Atoms
	// 		int tz = -1;
	// 		for (int bAtomIx : bAtomIxs_guest) {
	// 			tz++;
	// 			SimTK::Compound::AtomIndex aIx = (topology.subAtomList[bAtomIx]).identity.compoundAtomIndex;
	// 			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, forceField);
	// 			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

	// 			SimTK::Transform X_GB = mobod.getBodyTransform(someState);
	// 			SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
	// 			taskStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);

	// 			taskDeltaStationP[tz] = taskStationPInHost[tz] - taskStationPInGuest[tz];

	// 			if(this->visual == true){
	// 				paraMolecularDecorator->updateArrow(tz, taskStationPInGuest[tz], taskStationPInGuest[tz] + taskDeltaStationP[tz]);
	// 			}
	// 		}
	// 	}
	// }	
}

/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceStationPInGuest(void)
{
	return taskStationPInGuest;
}
/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceStationPInHost(void)
{
	return taskStationPInHost;
}

/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceDeltaStationP(void)
{
	return taskDeltaStationP;
}

/**
 * Calc Station Jacobian JS
*/
void World::calcStationJacobian(
	const SimTK::State&                           someState,
	SimTK::Matrix_<SimTK::Vec3>&                      JS) const
{
		matter->calcStationJacobian(someState, onBodyB, taskStationPInGuest, JS);

		std::cout << "Task Bodies ";
		std::cout << onBodyB << std::endl;
		std::cout << "Task Stations ";
		std::cout << taskStationPInGuest << std::endl;
		std::cout << "Station Jacobian ";
		std::cout << JS << std::endl;

		//matter->calcBiasForStationJacobian(someState, onBodyB, stationPInB, JSDotu);
}


//=============================================================================
//                   CONSTRAINTS
//=============================================================================

/**
 * Add contact constraints to specific bodies.
 **/
void World::addRodConstraint(SimTK::State& someState)
{

	int hostTopology = 0;
	int guestTopology = 1;

	std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology
	rodBodies.emplace_back(std::make_pair(SimTK::MobilizedBody(), SimTK::MobilizedBody()));
	conStationPInGuest.emplace_back(SimTK::Vec3());
	conStationPInHost.emplace_back(SimTK::Vec3());
	conDeltaStationP.emplace_back(SimTK::Vec3());

	// Get stations in host
	int topi = -1;
	for(auto& topology : topologies){
		topi++;
		if(topi == hostTopology){

			/* // Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_host) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, forceField);

				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				rodBodies[0].first = mbx;

				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				conStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);
			} */
		}
	}

	// Get stations in guest
	topi = -1;
	for(auto& topology : topologies){
		topi++;

		if(topi == guestTopology){

			/* // Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_guest) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, forceField);
				
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				rodBodies[0].second = mbx;
				
				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				conStationPInGuest[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);

			} */
		}
	}

/* 	std::cout << "Rod bodies: " << rodBodies[0].first
		<< " " << rodBodies[0].second << std::endl << std::flush;

	rodConstraints.emplace_back( SimTK::Constraint::Rod(
		matter->updMobilizedBody(rodBodies[0].first),  SimTK::Vec3(),
		matter->updMobilizedBody(rodBodies[0].second), SimTK::Vec3(), 0.1) ); */

	//rodConstraints.back().enable(someState);

}

/** Add speed constraints to specific bodies.
TODO:use number of mobilities. TODO: Solve if **/
const SimTK::State& World::addSpeedConstraint(int prmtopIndex)
{

	int hostTopology = 0;
	int guestTopology = 1;

	std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	if(prmtopIndex >= 0){
		std::cout << "Adding constraint to atom with prmtop index "
			<< prmtopIndex << "\n" ;
		SimTK::MobilizedBodyIndex mbx =
			topologies[0].getAtomMobilizedBodyIndexThroughDumm(
				SimTK::Compound::AtomIndex(prmtopIndex), *forceField
		);
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		SimTK::Constraint::ConstantSpeed B3291ConstraintU1(
			mobod, SimTK::MobilizerUIndex(0), 0);
		if(matter->getNumBodies() > 5000){
			SimTK::Constraint::ConstantSpeed B3291ConstraintU2(
				mobod, SimTK::MobilizerUIndex(1), 0);
			SimTK::Constraint::ConstantSpeed B3291ConstraintU3(
				mobod, SimTK::MobilizerUIndex(2), 0);
		}
	}

	const SimTK::State& returnState = compoundSystem->realizeTopology();
	return returnState;
}

//=============================================================================
//                   CONTACTS Functions
//=============================================================================

/** Add contact surfaces to bodies **/
void World::addContacts(const std::vector<int>& prmtopIndex, const int topologyIx, 
						const SimTK::ContactCliqueId cliqueId)
{
	// // Iterate through the prmtopIndex and add the appropriate spheres

	// for(int prmtopIx : (prmtopIndex)){
	// 	if (prmtopIx == -1) {
	// 		std::cout << "World::addContacts: prmtop index -1 found. Skipping topologyIx " << topologyIx << std::endl;
	// 		break;
	// 	}

	// 	std::cout << "Adding contacts with membrane to atom with prmtop index " <<
	// 		prmtopIx << " to topology " << topologyIx << " in clique " << cliqueId << std::endl;



	// 	const Real stiffness = 10000; // stiffness in pascals
	// 	const Real dissipation = 0.0; 
	// 	SimTK::Real staticFriction = 0.0;
	// 	SimTK::Real dynamicFriction = 0.0;
	// 	SimTK::Real viscousFriction = 0.0;

	// 	// Map from AmberAtomIndex to CompoundAtomIndex 
	// 	SimTK::Compound::AtomIndex cAIx = (topologies[topologyIx]).
	// 									subAtomList[prmtopIx].identity.compoundAtomIndex;

	// 	SimTK::MobilizedBodyIndex
	// 	mbx = (topologies[topologyIx]).getAtomMobilizedBodyIndexThroughDumm(
	// 		cAIx, forceField);
	// 	SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		
	// 	// Need to get vdw radius for said atom
	// 	const SimTK::DuMM::AtomIndex dAIx = (topologies[topologyIx]).getDuMMAtomIndex(cAIx);
	// 	Real vdwRadius = forceField->getAtomRadius(dAIx);

	// 	std::cout << "Atom element: " << forceField->getAtomElement(dAIx)
	// 			  << " Atom radius (nm): " << forceField->getAtomRadius(dAIx) << std::endl;

	// 	PolygonalMesh sphereMesh;
    // 	sphereMesh = PolygonalMesh::createSphereMesh(vdwRadius,1);
			
	// 	// We need to account for the position of the actual atom, so we
	// 	// must also apply a translation from the Origin of the MoBod to the
	// 	// location of the station coressponding to said atom.

	// 	SimTK::Vec3 atomPos = (topologies[topologyIx]).getAtomLocationInMobilizedBodyFrameThroughDumm
	// 							(cAIx,getForceField());


	// 	ContactGeometry::TriangleMesh sphere(sphereMesh);
	// 	mobod.updBody().addContactSurface(Transform(atomPos), 
	// 		ContactSurface(sphere,
	// 			ContactMaterial(stiffness, dissipation,
	// 		staticFriction, dynamicFriction, viscousFriction),
	// 			0.1).joinClique(cliqueId));

	// 	if (visual == true) {
	// 		DecorativeMesh showSphere(sphere.createPolygonalMesh());
	// 		mobod.updBody().addDecoration(Transform(atomPos), 
	// 			showSphere.setColor(Cyan).setOpacity(0.2));
	// 		//TODO: Remove this in production
	// 		mobod.updBody().addDecoration(Transform(atomPos), 
	// 			showSphere.setColor(Gray).setRepresentation(DecorativeGeometry::DrawWireframe));
	// 	}
	// }
}

// CONTACT DEBUG
/*
int numForces = updWorld(currentWorldIx)->contactForces->getNumContactForces(
	currentAdvancedState);
SimTK::Real dissEnergy = updWorld(currentWorldIx)->contactForces->
	getDissipatedEnergy(currentAdvancedState);
bool hasDefaultForceGenerator =
	updWorld(currentWorldIx)->contactForces->hasDefaultForceGenerator();

const MultibodySystem & mbs =
	updWorld(currentWorldIx)->contactForces->getMultibodySystem();
int nofMobods = mbs.getMatterSubsystem().getNumBodies();

const ContactTrackerSubsystem & cts =
	updWorld(currentWorldIx)->contactForces->getContactTrackerSubsystem();
int ctsNofSurfaces = cts.getNumSurfaces();

std::cout << "CONTACT INFO:"
	<< " #forces= " << numForces	const SimTK::State& addContacts(int prmtopIndex);

	const SimTK::State& realizeTopology();

	<< " dissEnergy= " << dissEnergy
	<< " hasDefaultForceGenerator= " << hasDefaultForceGenerator
	<< " #mobods= " << nofMobods
	<< " ctsNofSurfaces= " << ctsNofSurfaces
<< std::endl;
*/
// CONTACT DEBUG enD
/** Assign a scale factor for generalized velocities to every mobilized
body **/

//=============================================================================
//                   MEMBRANE Functions
//=============================================================================

/** Add a membrane represented by a contact surface **/
void World::addMembrane(const SimTK::Real halfThickness)
{
	// SimTK::Real stiffness = 10000;
	// SimTK::Real dissipation = 0;
	// SimTK::Real staticFriction = 0.0;
	// SimTK::Real dynamicFriction = 0.0;
	// SimTK::Real viscousFriction = 0.0;

 	// matter->Ground().updBody().addContactSurface(
	// 	Transform(Rotation(-0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, halfThickness)),
	// 	ContactSurface(
	// 	ContactGeometry::HalfSpace(),
	// 	ContactMaterial(stiffness, dissipation,
	// 		staticFriction, dynamicFriction, viscousFriction))
	// 		.joinClique(SimTK::ContactCliqueId(1))
	// 		.joinClique(SimTK::ContactCliqueId(2))
	// 		.joinClique(SimTK::ContactCliqueId(3))
	// 		);

	// matter->Ground().updBody().addContactSurface(
	// 	Transform(Rotation(-0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, -halfThickness)),
	// 	ContactSurface(
	// 	ContactGeometry::HalfSpace(),
	// 	ContactMaterial(stiffness, dissipation,
	// 		staticFriction, dynamicFriction, viscousFriction))
	// 		.joinClique(SimTK::ContactCliqueId(0))
	// 		.joinClique(SimTK::ContactCliqueId(2))
	// 		.joinClique(SimTK::ContactCliqueId(3)));
	
	// matter->Ground().updBody().addContactSurface(
	// 	Transform(Rotation(0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, halfThickness)),
	// 	ContactSurface(
	// 	ContactGeometry::HalfSpace(),
	// 	ContactMaterial(stiffness, dissipation,
	// 		staticFriction, dynamicFriction, viscousFriction))
	// 		.joinClique(SimTK::ContactCliqueId(0))
	// 		.joinClique(SimTK::ContactCliqueId(1))
	// 		.joinClique(SimTK::ContactCliqueId(3))
	// 		);

	// matter->Ground().updBody().addContactSurface(
	// 	Transform(Rotation(0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, -halfThickness)),
	// 	ContactSurface(
	// 	ContactGeometry::HalfSpace(),
	// 	ContactMaterial(stiffness, dissipation,
	// 		staticFriction, dynamicFriction, viscousFriction))
	// 		.joinClique(SimTK::ContactCliqueId(0))
	// 		.joinClique(SimTK::ContactCliqueId(1))
	// 		.joinClique(SimTK::ContactCliqueId(2))
	// 		);


	// // if (visual == true) {
	// // 	DecorativeFrame contactGeometryDecoFrame;
	// // 	matter->Ground().updBody().addDecoration(
	// // 	Transform(),
    // //     DecorativeBrick(Vec3(10,10,halfThickness)).setColor(Orange).setOpacity(0.25));
	// // }
}

/*! 
 * <!-- Assign a scale factor for generalized velocities to every mobilized
 * body -->
 */
void World::setUScaleFactorsToMobods(void)
{
	// //for(auto& topology : topologies){ // SAFE
	// for(auto& topology : topologies){ // DANGER
	// 	// Iterate bonds

	// 	//for(const auto& AtomList : topology.bAtomList){
	// 	for(const auto& bond : topology.getBonds()){
	// 		bond.getBondGlobalIndex();
	// 		topology.getAtomINdex
	// 		SimTK::Compound::AtomIndex aIx1 = topology.subAtomList[Bond.i].identity.compoundAtomIndex;
	// 		SimTK::Compound::AtomIndex aIx2 = topology.subAtomList[Bond.j].identity.compoundAtomIndex;

	// 		SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndexThroughDumm(aIx1, forceField);
	// 		SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndexThroughDumm(aIx2, forceField);

	// 		std::cout
	// 		<< "DEBUG World::setUScaleFactorsToMobods aIx1 aIx2 mbx1 mbx2 "
	// 		<< aIx1 << " " << aIx2 << " "
	// 		<< mbx1 << " " << mbx2 << std::endl;

	// 		const SimTK::MobilizedBody& mobod1 = matter->getMobilizedBody(mbx1);
	// 		const SimTK::MobilizedBody& mobod2 = matter->getMobilizedBody(mbx2);

	// 		int level1 = mobod1.getLevelInMultibodyTree();
	// 		int level2 = mobod2.getLevelInMultibodyTree();

	// 		if(level1 > level2){
	// 			mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx1, Bond.getUScaleFactor(ownWorldIndex)));
	// 		}else if(level2 > level1){
	// 			mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx2, Bond.getUScaleFactor(ownWorldIndex)));
	// 		}else{
	// 			if(Bond.getUScaleFactor(ownWorldIndex) != 0){
	// 				std::cout << "World::setUScaleFactorsToMobods Warning: Trying to scale a bond inside a rigid body\n";
	// 			}
	// 		}
	// 	}
	// }
}

// Get the (potential) energy transfer
// If any of the Q, U or tau is actively modifyied by the sampler
// the Jacobian of that transformation will be included too
SimTK::Real World::getWorkOrHeat(void)
{
	// Accumulate in this variable
	SimTK::Real retValue = 0.0;

	// Get the energy transfer from all the samplers
	for(auto& sampler : this->samplers){

		// Get the potential energy difference
		retValue += getSampler(0)->getCurrentEnergy().potential - getSampler(0)->getPreviousEnergy().potential;

		// Get Fixman potential difference
		retValue += getSampler(0)->getCurrentEnergy().fixman - getSampler(0)->getPreviousEnergy().fixman;

		/* // Get the Q modifying samplers Jacobians
		if(sampler->getDistortOpt() < 0){
			retValue -= 
				sampler->getDistortJacobianDetLog();
		} */
	}
	
	return retValue;
}

// Get the (potential) energy transfer in the form of work
// If any of the Q, U or tau is actively modifyied by the sampler
// the Jacobian of that transformation will be included too
SimTK::Real World::getWork(void) const
{
	// Accumulate in this variable
	SimTK::Real retValue = 0.0;

	// Get the energy transfer from all the samplers
	for(auto& sampler : this->samplers){

		if(sampler->getDistortOpt() < 0){

			// Get the potential energy difference
			retValue += getSampler(0)->getCurrentEnergy().potential - getSampler(0)->getPreviousEnergy().potential;

			// Get Fixman potential difference
			retValue += getSampler(0)->getCurrentEnergy().fixman - getSampler(0)->getPreviousEnergy().fixman;
			/* // Get the Jacobians
			retValue -= 
				sampler->getDistortJacobianDetLog(); */
		}
		
	}
	
	return retValue;
}

/*
 * Shift all the generalized coordinates
 */
void World::getTransformsStatistics(SimTK::State& someState)
{
	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		// Get BAT coordinate "angle"
		/*
		/ cos t
		| 
		|
		\
		*/
		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		/* std::cout << "World " << ownWorldIndex << " " 
			//<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			//<< "angleMean " << int(mbx) - 1 << " "
			//<< angleMean * (180 / SimTK::Pi) << " "
			//<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */

	}

}

// Print bond lengths and angle bends
void World::traceBendStretch(SimTK::State& someState){
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){
		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		
		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
		//PrintTransform(X_PF, 4, "X_PF " + std::to_string(int(mbx)));

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
		//PrintTransform(X_BM, 4, "X_BM " + std::to_string(int(mbx)));

		// Get BAT coordinates B and A
		SimTK::Vec3 bondVector = X_BM.p();
		trace( bondVector.norm(),  std::acos(X_PF.R()(0)(0)));
		//trace("X_FM");
		//PrintTransform(X_FM, 10);
	
	}
}

// Print X_PF means
void World::PrintAcosX_PFs(void)
{
	int i = -1;
	for(const auto &xpf : acosX_PF00 ){
		i += 1;
		//std::cout << "X_PF " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
		std::cout << "acosX_PF " << i << " " << xpf << std::endl;

	}
}

// Print X_PF means
void World::PrintNormX_BMs(void)
{
	int i = -1;
	for(const auto &xbm : normX_BMp ){
		i += 1;
		std::cout << "normX_BM " << i << " " << xbm << std::endl;
	}
}

// Print X_PF means
void World::PrintAcosX_PFMeans(void)
{
	int i = -1;
	for(const auto &xpf : acosX_PF00_means ){
		i += 1;
		//std::cout << "X_PFMean " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
		std::cout << "acosX_PFMean " << i << " " << xpf << std::endl;
	}
}

// Print X_PF means
void World::PrintNormX_BMMeans(void)
{
	int i = -1;
	for(const auto &xbm : normX_BMp_means ){
		i += 1;
		std::cout << "normX_BMMean " << i << " " << xbm << std::endl;
	}
}

// REORIENT

SimTK::Transform& World::getReorientTransformInAnotherBody(
	const SimTK::State &someState,
	const SimTK::MobilizedBody &inBodyA, const SimTK::MobilizedBody &ofBodyB,
	const SimTK::Transform &reorientAB,
	SimTK::Transform& X_FMprim)
{
	SimTK::Transform X_MB = ~(ofBodyB.getOutboardFrame(someState));
	SimTK::Transform X_FM = ofBodyB.getMobilizerTransform(someState);
	SimTK::Transform X_AB = 
		ofBodyB.findBodyTransformInAnotherBody(someState, inBodyA);

	SimTK::Transform X_BBprim = (~X_AB) * reorientAB;
	X_FMprim = X_FM * X_MB * X_BBprim * X_MB;

	return X_FMprim;
}

//...............

/**
 * Set X_PF, X_BM means
*/
void World::setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
const std::vector<SimTK::Real>& givenX_BM)
{
	// Update acosX_PF00 means
	int i = -1;
	for(auto &xpf : acosX_PF00_means ){
		i += 1;
		xpf = givenX_PF[i]; 
	}
	
	// Update normX_BMp means
	i = -1;
	for(auto &xbm : normX_BMp_means ){
		i += 1;
		xbm = givenX_BM[i];
	}
}

/**
 * Set X_PF, X_BM means to initial values
*/
void World::setTransformsMeansToIni(void)
{
	const SimTK::State& defaultState = matter->getSystem().getDefaultState();

	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(defaultState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(defaultState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(defaultState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		/* SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		std::cout 
			<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			<< "angleMean " << int(mbx) - 1 << " "
			<< angleMean * (180 / SimTK::Pi) << " "
			<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */
	}
}

/*
 * Shift all the generalized coordinates
 */
void World::setTransformsMeansToCurrent(SimTK::State& someState)
{
	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00_means[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp_means[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		/* SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		std::cout 
			<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			<< "angleMean " << int(mbx) - 1 << " "
			<< angleMean * (180 / SimTK::Pi) << " "
			<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */

	}

}

/**
 * Update X_PF, X_BM means
*/
void World::updateTransformsMeans(SimTK::State& someState)
{
	int nofSamples = getNofSamples() + 1;
	//std::cout << "Nof samples " << nofSamples << std::endl;

	// Useful vars
	SimTK::Real N_1overN = 9999, NInv = 9999;

	if(nofSamples == 1){
		for(unsigned int k = 0; k < acosX_PF00_means.size(); k++){
			acosX_PF00_means[k] = acosX_PF00[k];
		}
		for(unsigned int k = 0; k < normX_BMp_means.size(); k++){
			normX_BMp_means[k] = normX_BMp[k];
			//std::cout << "World " << ownWorldIndex << " bondUpdMean " << k << " " << normX_BMp_means[k] << std::endl;
		}
	}else{
		if(nofSamples == 2){
			N_1overN = NInv = 0.5;
		}else{
			// Update useful vars
			SimTK::Real N_1 = nofSamples - 1.0;
			N_1overN = N_1 / nofSamples;
			NInv = 1.0 / nofSamples;
		}
		//std::cout << "updateX_PFMeans check " << " "
		//	<<  N_1overN << " " <<  NInv  << " " << std::flush;

		// Update acosX_PF00 means
		int i = -1;
		for(auto &xpf : acosX_PF00_means ){
			i += 1;
			xpf = (N_1overN * xpf) + (NInv * acosX_PF00.at(i)); 
		}
		
		// Update normX_BMp means
		i = -1;
		for(auto &xbm : normX_BMp_means ){
			i += 1;
			xbm = (N_1overN * xbm) + (NInv * normX_BMp.at(i));
			//std::cout << "World " << ownWorldIndex << " bondUpdMean " << i << " " << xbm << std::endl;
		}
	}

}

// Get X_PF means
std::vector<SimTK::Real>& World::getX_PFMeans(void)
{
	return acosX_PF00_means;
}

// Get X_BM means
std::vector<SimTK::Real>& World::getX_BMMeans(void)
{
	return normX_BMp_means;
}

/**
 * Calculate bond length and angle deviations from their means
*/ 
void World::calcBendStretchDeviations(SimTK::State& someState,
	std::vector<SimTK::Real>& X_PFdiffs,
	std::vector<SimTK::Real>& X_BMdiffs
)
{
	// Make sure it has 
	X_PFdiffs.resize(this->acosX_PF00_means.size(), 0.0);
	X_BMdiffs.resize(this->normX_BMp_means.size(), 0.0);

	// 
	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		X_PFdiffs[k] = this->acosX_PF00[k] - this->acosX_PF00_means[k];
	}
	for(unsigned int k = 0; k < X_BMdiffs.size(); k++){
		X_BMdiffs[k] = this->normX_BMp[k] - this->normX_BMp_means[k];
		/* std::cout << "World " << ownWorldIndex << " bondDiff " << k << " "
		//<< this->normX_BMp[k] << " " << this->normX_BMp_means[k] << " "
		<< X_BMdiffs[k] << std::endl; */
	}

}

/** Get MobilizedBody to AtomIndex map **/
std::map< SimTK::MobilizedBodyIndex, std::pair<int, SimTK::Compound::AtomIndex>>&
World::getMbx2aIx(){
	return mbx2aIx;
}


std::size_t World::getNofMobilizedBodies() const{
	return mbx2aIx.size();
}

/** Get U scale factor for the mobilized body **/
SimTK::Real World::getMobodUScaleFactor(SimTK::MobilizedBodyIndex& mbx) const
{
	if(!mbx2uScale.empty()){
		if(mbx2uScale.find(mbx) != mbx2uScale.end()){
			return mbx2uScale.at(mbx);
		}else{
			//std::cout << "Warning: U scale factor for mobod " << int(mbx) << " not found.\n";
			return 1;
		}
	}else{
		return 1;
	}
}

/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index maps **/
void World::printMaps(void)
{
	for (auto& topology : topologies){
		topology.printMaps();
	}
}

SimTK::Vec3 World::calcAtomLocationInGroundFrameThroughOMM(const SimTK::DuMM::AtomIndex& )
{
	assert("NOT IMPLEMENTED!");
	return SimTK::Vec3(-1, -1, -1);
}

//==============================================================================
//                             2. Inter-world functions.
//==============================================================================
// Pass configurations between Worlds

/*! <!--  --> */
void World::PrintBATFromSimbody() const
{
	const SimTK::State& advState = integrator->getAdvancedState();

	bool parFlag = false;
	bool parParFlag = false;
	int childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;
	bool printTransforms = true;

	std::vector<SimTK::Real> BONDLengths(matter->getNumBodies() - 1, -99999);
	std::vector<SimTK::Real> ANGLEBends(matter->getNumBodies() - 1, -99999);
	std::vector<SimTK::Real> TORSIONAngles(matter->getNumBodies() - 1, -99999);
	std::vector<std::vector<int>> ZMatrix(matter->getNumBodies() - 1, std::vector<int>(4, -99999));

	for (SimTK::MobilizedBodyIndex childMbx(1); childMbx < matter->getNumBodies(); ++childMbx){
		childIx = int(childMbx);

		const SimTK::MobilizedBody& childMobod = matter->getMobilizedBody(childMbx);
		const SimTK::MobilizedBody& parentMobod = childMobod.getParentMobilizedBody();
		const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
		parentIx = int(parentMbx);

		if(int(childMbx) > 1){
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();
			grandIx = int(grandMbx);
		}

		if(int(childMbx) > 2){
			const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandGrandMbx = grandGrandMobod.getMobilizedBodyIndex();
			grandGrandIx = int(grandGrandMbx);
		}

		// Print out the indices
		std::cout << "World " << ownWorldIndex << " " << "child " << childIx << " " << "parent " << parentIx << " " << "parPar " << grandIx << " " << "grandGrandIx " << grandGrandIx << " " << std::endl;

		// BOND ==============
		const SimTK::Transform& B_X_Fb = childMobod.getInboardFrame(advState);
		const SimTK::Transform& C_X_Mb = childMobod.getOutboardFrame(advState);
		const SimTK::Transform& Fb_X_Mb = childMobod.getMobilizerTransform(advState);

		const SimTK::Transform& G_X_C = childMobod.getBodyTransform(advState);
		const SimTK::Transform& G_X_B = parentMobod.getBodyTransform(advState);
		SimTK::Transform G_X_Fb = G_X_B * B_X_Fb;
		SimTK::Transform G_X_Mb = G_X_C * C_X_Mb;

		SimTK::Transform B_X_C = B_X_Fb * Fb_X_Mb * (~C_X_Mb);
		//BONDLengths[int(childMbx) - 1] = B_X_C.p().norm(); // correct
		BONDLengths[int(childMbx) - 1] = C_X_Mb.p().norm(); // correct correct
		
		ZMatrix[int(childMbx) - 1][0] = int(childMbx);
		ZMatrix[int(childMbx) - 1][1] = int(parentMbx);

		if(printTransforms){
			//SimTK::Test::PrintTransform(G_X_C, 6, "G_X_C", "G_X_C:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			// SimTK::Test::PrintTransform(G_X_Mb, 6, "G_X_Mb", "G_X_Mb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			// SimTK::Test::PrintTransform(G_X_Fb, 6, "G_X_Fb", "G_X_Fb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(G_X_B, 6, "G_X_B", "G_X_B:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			SimTK::Test::PrintTransform(B_X_C, 6, "B_X_C", "B_X_C:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(B_X_Fb, 6, "B_X_Fb", "B_X_Fb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			SimTK::Test::PrintTransform(Fb_X_Mb, 6, "Fb_X_Mb", "Fb_X_Mb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			SimTK::Test::PrintTransform(C_X_Mb, 6, "C_X_Mb", "C_X_Mb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
		}

		if(int(childMbx) > 1){ // ANGLE ========================
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

			const SimTK::Transform& A_X_Fa = parentMobod.getInboardFrame(advState); // A_X_Fa
			const SimTK::Transform& B_X_Ma = parentMobod.getOutboardFrame(advState); // B_X_Ma
			const SimTK::Transform& Fa_X_Ma = parentMobod.getMobilizerTransform(advState); // Fa_X_Ma
			SimTK::Transform A_X_B = A_X_Fa * Fa_X_Ma * (~B_X_Ma);
			SimTK::Transform A_X_C = A_X_B * B_X_C;

			const SimTK::Transform& G_X_A = grandMobod.getBodyTransform(advState); // G_X_A
			SimTK::Transform G_X_Fa = G_X_A * A_X_Fa;
			SimTK::Transform G_X_Ma = G_X_B * B_X_Ma;

			// Checks
			SimTK::Vec3 pAXB_A = A_X_B.p();
			SimTK::Vec3 pAXB_B = ~(A_X_B.R()) * pAXB_A;
			SimTK::Vec3 pBXC_B = B_X_C.p();


			SimTK::Vec3 pAXB_G = ~(G_X_A.R()) * pAXB_A;
			SimTK::Vec3 pBXC_G = ~(G_X_B.R()) * pBXC_B;

			// WORK ==========================================
			SimTK::Vec3 xAXB_A = A_X_B.R()(0);
			SimTK::Vec3 xAXB_G = ~(G_X_A.R()) * xAXB_A;

			SimTK::Vec3 xBXC_B = B_X_C.R()(0);
			SimTK::Vec3 xBXC_G = ~(G_X_B.R()) * xBXC_B;

			// std::cout << "check: " << int(childMbx)
			// 	<<" "<< std::acos(SimTK::dot(pAXB_A.normalize(), pBXC_B.normalize()))
			// 	<<" "<< std::acos(SimTK::dot(xAXB_G, xBXC_G))
			// 	<<std::endl;

			SimTK::Vec3 xBXFb_B = (B_X_Fb.R())(0);
			SimTK::Vec3 xBXMa_B = (B_X_Ma.R())(0);
			SimTK::Vec3 _xBXMa_B = -1 * xBXMa_B;

			//ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot(pAXB_A.normalize(), pBXC_B.normalize()));
			//ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot(pAXB_G.normalize(), pBXC_G.normalize()));
			ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot((-1 * pAXB_B).normalize(), pBXC_B.normalize())); // correct
			//ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot( xBXFb_B, _xBXMa_B )); // correct correct
			//ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot(xAXB_G, xBXC_G));
			// ================================================

			ZMatrix[int(childMbx) - 1][2] = int(grandMbx);

			if(printTransforms){
				// SimTK::Test::PrintTransform(G_X_Ma, 6, "G_X_Ma", "G_X_Ma:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintTransform(G_X_Fa, 6, "G_X_Fa", "G_X_Fa:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(G_X_A, 6, "G_X_A", "G_X_A:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandIx)));
				SimTK::Test::PrintTransform(A_X_B, 6, "A_X_B", "A_X_B:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(A_X_Fa, 6, "A_X_Fa", "A_X_Fa:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandMbx)));
				SimTK::Test::PrintTransform(Fa_X_Ma, 6, "Fa_X_Ma", "Fa_X_Ma:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandMbx)));
				SimTK::Test::PrintTransform(B_X_Ma, 6, "B_X_Ma", "B_X_Ma:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			}

			if(int(childMbx) > 2){ // TORSION =======================
				const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();

				const SimTK::Transform& T_X_Ft = grandMobod.getInboardFrame(advState); // T_X_Ft
				const SimTK::Transform& A_X_Mt = grandMobod.getOutboardFrame(advState); // A_X_Mt
				const SimTK::Transform& G_X_T = grandGrandMobod.getBodyTransform(advState); // G_X_T
				const SimTK::Transform& Ft_X_Mt = grandMobod.getMobilizerTransform(advState); // Ft_X_Mt
				SimTK::Transform T_X_A = T_X_Ft * Ft_X_Mt * (~A_X_Mt);

				SimTK::Transform G_X_Ft = G_X_T * T_X_Ft;
				SimTK::Transform G_X_Mt = G_X_A * A_X_Mt;

				// Checks

				// WORK ==========================================

				SimTK::Vec3 xAXFa_A = (A_X_Fa.R())(0);
				SimTK::Vec3 xAXMt_A = (A_X_Mt.R())(0);
				
				SimTK::Vec3 v2_B = SimTK::cross(xBXFb_B, xBXMa_B);
				SimTK::Vec3 v1_A = SimTK::cross(xAXFa_A, xAXMt_A);

				SimTK::Vec3 v2_B_hat = v2_B.normalize();
				SimTK::Vec3 v1_A_hat = v1_A.normalize();

				SimTK::Vec3 v1_B_hat = (~(A_X_B.R())) * v1_A_hat;

				SimTK::Vec3 v3_B = SimTK::cross(v2_B_hat, v1_B_hat);
				SimTK::Vec3 v3_B_hat = v3_B.normalize();

				SimTK::Real tors_cos = SimTK::dot(v1_B_hat, v2_B_hat);
				SimTK::Real tors_sin = SimTK::dot(v3_B, _xBXMa_B);

				// SimTK::Test::PrintVec3(v1_A, 6, "v1_A", "v1_A:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintVec3(v2_B, 6, "v2_B", "v2_B:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintVec3(v3_B, 6, "v3_B", "v3_B:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));

				// SimTK::Test::PrintVec3(v1_B_hat, 6, "v1_B_hat", "v1_B_hat:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintVec3(v2_B_hat, 6, "v2_B_hat", "v2_B_hat:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintVec3(v3_B_hat, 6, "v3_B_hat", "v3_B_hat:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));

				//std::cout << " tors cos sin " << tors_cos <<" "<< tors_sin << std::endl;

				TORSIONAngles[int(childMbx) - 1] = std::atan2(tors_sin, tors_cos);
				// ==============================================

				ZMatrix[int(childMbx) - 1][3] = int(grandGrandIx);

				if(printTransforms){
					//SimTK::Test::PrintTransform(G_X_Mt, 6, "G_X_Mt", "G_X_Mt:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					//SimTK::Test::PrintTransform(G_X_Ft, 6, "G_X_Ft", "G_X_Ft:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(G_X_T, 6, "G_X_T", "G_X_T:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(T_X_A, 6, "T_X_A", "T_X_A:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(T_X_Ft, 6, "T_X_Ft", "T_X_Ft:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(Ft_X_Mt, 6, "Ft_X_Mt", "Ft_X_Mt:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(A_X_Mt, 6, "A_X_Mt", "A_X_Mt:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
				}


			}
		}

		childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;

		/*
		const SimTK::Transform& X_PF = mobod.getInboardFrame(advState);
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(advState);
		SimTK::Test::PrintTransform(X_PF, 6, "X_PF", "X_PF:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
		SimTK::Test::PrintTransform(X_BM, 6, "X_BM", "X_BM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));

		const SimTK::MobilizedBody& parentMobod = mobod.getParentMobilizedBody();
		const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

		if(int(parentMbx) != 0){
			parFlag = true;

			const SimTK::Transform& parX_PF = parentMobod.getInboardFrame(advState);
			const SimTK::Transform& parX_BM = parentMobod.getOutboardFrame(advState);
			SimTK::Test::PrintTransform(parX_PF, 6, "parX_PF", "parX_PF:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			SimTK::Test::PrintTransform(parX_BM, 6, "parX_BM", "parX_BM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));

			const SimTK::MobilizedBody& parParMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex parParMbx = parParMobod.getMobilizedBodyIndex();

			if(int(parParMbx) != 0){
				parParFlag = true;

				const SimTK::Transform& parParX_PF = parParMobod.getInboardFrame(advState);
				const SimTK::Transform& parParX_BM = parParMobod.getOutboardFrame(advState);
				SimTK::Test::PrintTransform(parParX_PF, 6, "parParX_PF", "parParX_PF:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parParMbx)));
				SimTK::Test::PrintTransform(parParX_BM, 6, "parParX_BM", "parParX_BM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parParMbx)));

				const SimTK::MobilizedBody& parParParMobod = parParMobod.getParentMobilizedBody();
				const SimTK::MobilizedBodyIndex parParParMbx = parParParMobod.getMobilizedBodyIndex();

			}else{ // _end_ if parent parent
				;
			} // _end_ if parent parent else

		}else{ // _end_ if parent
			;
		} // _end_ if parent else

		parFlag = false;
		parParFlag = false;
		child = -1, parent = -1, parPar = -1, grandGrandIx = -2; */

	} // _end_ for mbx

	// Print
	for (int BOIx = 0; BOIx < BONDLengths.size(); BOIx++){
		std::cout << "ZMatrixBATSimbody:"
			<<" " << ZMatrix[BOIx][0] << " " << ZMatrix[BOIx][1] << " " << ZMatrix[BOIx][2] << " " << ZMatrix[BOIx][3]
			<<" "<< BONDLengths[BOIx] << " " << ANGLEBends[BOIx] << " " << TORSIONAngles[BOIx] << std::endl;
	}


}


/*! <!--  --> */
void World::calcSimbodyBAT(std::vector<std::vector<int>>& ZMatrix, std::vector<SimTK::Real>& BONDLengths, std::vector<SimTK::Real>& ANGLEBends, std::vector<SimTK::Real>& TORSIONAngles)
{
	SimTK::State& advState = integrator->updAdvancedState();

	bool parFlag = false;
	bool parParFlag = false;
	int childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;

	if(BONDLengths.size() != matter->getNumBodies() -1){
		BONDLengths.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(ANGLEBends.size() != matter->getNumBodies() -1){
		ANGLEBends.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(TORSIONAngles.size() != matter->getNumBodies() -1){
		TORSIONAngles.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(ZMatrix.size() != matter->getNumBodies() -1){
		ZMatrix.resize(matter->getNumBodies() - 1, std::vector<int>(4, -1));
	}

	bool printTransforms = false;

	for (SimTK::MobilizedBodyIndex childMbx(1); childMbx < matter->getNumBodies(); ++childMbx){
		childIx = int(childMbx);

		const SimTK::MobilizedBody& childMobod = matter->getMobilizedBody(childMbx);
		const SimTK::MobilizedBody& parentMobod = childMobod.getParentMobilizedBody();
		const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
		parentIx = int(parentMbx);

		if(int(childMbx) > 1){
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();
			grandIx = int(grandMbx);
		}

		if(int(childMbx) > 2){
			const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandGrandMbx = grandGrandMobod.getMobilizedBodyIndex();
			grandGrandIx = int(grandGrandMbx);
		}

		// Print out the indices
		if(printTransforms){
			std::cout << "World " << ownWorldIndex << " " << "child " << childIx << " " << "parent " << parentIx << " " << "parPar " << grandIx << " " << "grandGrandIx " << grandGrandIx << " " << std::endl << std::flush;
		}

		// BOND ==============
		const SimTK::Transform& B_X_Fb = childMobod.getInboardFrame(advState);
		const SimTK::Transform& C_X_Mb = childMobod.getOutboardFrame(advState);
		const SimTK::Transform& Fb_X_Mb = childMobod.getMobilizerTransform(advState);

		const SimTK::Transform& G_X_C = childMobod.getBodyTransform(advState);
		const SimTK::Transform& G_X_B = parentMobod.getBodyTransform(advState);
		SimTK::Transform G_X_Fb = G_X_B * B_X_Fb;
		SimTK::Transform G_X_Mb = G_X_C * C_X_Mb;

		SimTK::Transform B_X_C = B_X_Fb * Fb_X_Mb * (~C_X_Mb);
		BONDLengths[int(childMbx) - 1] = B_X_C.p().norm(); // correct
		
		ZMatrix[int(childMbx) - 1][0] = int(childMbx);
		ZMatrix[int(childMbx) - 1][1] = int(parentMbx);

		if(printTransforms){
			//SimTK::Test::PrintTransform(G_X_C, 6, "G_X_C", "G_X_C:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			//SimTK::Test::PrintTransform(G_X_Mb, 6, "G_X_Mb", "G_X_Mb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			//SimTK::Test::PrintTransform(G_X_Fb, 6, "G_X_Fb", "G_X_Fb:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(G_X_B, 6, "G_X_B", "G_X_B:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
			SimTK::Test::PrintTransform(B_X_C, 6, "B_X_C", "B_X_C:" + std::to_string(int(parentMbx)) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(B_X_Fb, 6, "B_X_Fb", "B_X_Fb:" + std::to_string(int(parentMbx)) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(Fb_X_Mb, 6, "Fb_X_Mb", "Fb_X_Mb:" + std::to_string(int(parentMbx)) + ":" + std::to_string(int(childMbx)));
			SimTK::Test::PrintTransform(C_X_Mb, 6, "C_X_Mb", "C_X_Mb:" + std::to_string(int(parentMbx)) + ":" + std::to_string(int(childMbx)));
		}

		if(int(childMbx) > 1){ // ANGLE ========================
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

			const SimTK::Transform& A_X_Fa = parentMobod.getInboardFrame(advState); // A_X_Fa
			const SimTK::Transform& B_X_Ma = parentMobod.getOutboardFrame(advState); // B_X_Ma
			const SimTK::Transform& Fa_X_Ma = parentMobod.getMobilizerTransform(advState); // Fa_X_Ma
			SimTK::Transform A_X_B = A_X_Fa * Fa_X_Ma * (~B_X_Ma);
			SimTK::Transform B_X_A = ~A_X_B;
			SimTK::Transform A_X_C = A_X_B * B_X_C;

			const SimTK::Transform& G_X_A = grandMobod.getBodyTransform(advState); // G_X_A
			SimTK::Transform G_X_Fa = G_X_A * A_X_Fa;
			SimTK::Transform G_X_Ma = G_X_B * B_X_Ma;

			// Checks
			SimTK::Vec3 pAXB_A = A_X_B.p();
			//SimTK::Vec3 pAXB_B = ~(A_X_B.R()) * pAXB_A;
			SimTK::Vec3 pBXA_B = B_X_A.p();
			SimTK::Vec3 pBXC_B = B_X_C.p();

			ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot((-1 * pBXA_B).normalize(), pBXC_B.normalize())); // correct
			// ================================================

			ZMatrix[int(childMbx) - 1][2] = int(grandMbx);

			if(printTransforms){
				// SimTK::Test::PrintTransform(G_X_Ma, 6, "G_X_Ma", "G_X_Ma:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				// SimTK::Test::PrintTransform(G_X_Fa, 6, "G_X_Fa", "G_X_Fa:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(G_X_A, 6, "G_X_A", "G_X_A:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandIx)));
				SimTK::Test::PrintTransform(A_X_B, 6, "A_X_B", "A_X_B:" + std::to_string(int(grandMbx)) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(A_X_Fa, 6, "A_X_Fa", "A_X_Fa:" + std::to_string(int(grandMbx)) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(Fa_X_Ma, 6, "Fa_X_Ma", "Fa_X_Ma:" + std::to_string(int(grandMbx)) + ":" + std::to_string(int(parentMbx)));
				SimTK::Test::PrintTransform(B_X_Ma, 6, "B_X_Ma", "B_X_Ma:" + std::to_string(int(grandMbx)) + ":" + std::to_string(int(parentMbx)));
			}

			if(int(childMbx) > 2){ // TORSION =======================
				const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();

				const SimTK::Transform& T_X_Ft = grandMobod.getInboardFrame(advState); // T_X_Ft
				const SimTK::Transform& A_X_Mt = grandMobod.getOutboardFrame(advState); // A_X_Mt
				const SimTK::Transform& G_X_T = grandGrandMobod.getBodyTransform(advState); // G_X_T
				const SimTK::Transform& Ft_X_Mt = grandMobod.getMobilizerTransform(advState); // Ft_X_Mt
				SimTK::Transform T_X_A = T_X_Ft * Ft_X_Mt * (~A_X_Mt);
				SimTK::Transform T_X_B = T_X_A * A_X_B;

				//Transform G_X_Ft = G_X_T * T_X_Ft;
				//Transform G_X_Mt = G_X_A * A_X_Mt;

				// Checks

				// WORK ==========================================
				SimTK::Vec3 b1_B = (T_X_B.R()) * T_X_A.p();
				SimTK::Vec3 b2_B = -1.0 * B_X_A.p();
				SimTK::Vec3 b3_B = B_X_C.p();

				SimTK::Vec3 b1_B_hat = b1_B.normalize();
				SimTK::Vec3 b2_B_hat = b2_B.normalize();
				SimTK::Vec3 b3_B_hat = b3_B.normalize();

				SimTK::Vec3 n1_B_hat = (SimTK::cross(b1_B, b2_B)).normalize();
				SimTK::Vec3 n2_B_hat = (SimTK::cross(b2_B, b3_B)).normalize();

				// tors_cos
				SimTK::Real tors_cos = SimTK::dot(n1_B_hat, n2_B_hat);

				//tors_sin
				SimTK::Vec3 m1_B = SimTK::cross(n1_B_hat, b2_B_hat);
				SimTK::Real tors_sin = SimTK::dot(m1_B, n2_B_hat);

				//std::cout << " tors cos sin " << tors_cos <<" "<< tors_sin << std::endl;

				TORSIONAngles[int(childMbx) - 1] = std::atan2(tors_sin, tors_cos);
				// ==============================================

				ZMatrix[int(childMbx) - 1][3] = int(grandGrandIx);

				if(printTransforms){
					//SimTK::Test::PrintTransform(G_X_Mt, 6, "G_X_Mt", "G_X_Mt:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					//SimTK::Test::PrintTransform(G_X_Ft, 6, "G_X_Ft", "G_X_Ft:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(G_X_T, 6, "G_X_T", "G_X_T:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
					SimTK::Test::PrintTransform(T_X_A, 6, "T_X_A", "T_X_A:" + std::to_string(int(grandGrandIx)) + ":" + std::to_string(int(grandMbx)));
					SimTK::Test::PrintTransform(T_X_Ft, 6, "T_X_Ft", "T_X_Ft:" + std::to_string(int(grandGrandIx)) + ":" + std::to_string(int(grandMbx)));
					SimTK::Test::PrintTransform(Ft_X_Mt, 6, "Ft_X_Mt", "Ft_X_Mt:" + std::to_string(int(grandGrandIx)) + ":" + std::to_string(int(grandMbx)));
					SimTK::Test::PrintTransform(A_X_Mt, 6, "A_X_Mt", "A_X_Mt:" + std::to_string(int(grandGrandIx)) + ":" + std::to_string(int(grandMbx)));
				}

			}
		}

		childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;


	} // _end_ for mbx

	// Print
	bool printZmatBAT = false;
	if(printZmatBAT){
		for (int BOIx = 0; BOIx < BONDLengths.size(); BOIx++){
			std::cout << "ZMatrixBATSimbody:"
				<<" " << ZMatrix[BOIx][0] << " " << ZMatrix[BOIx][1] << " " << ZMatrix[BOIx][2] << " " << ZMatrix[BOIx][3]
				<<" "<< BONDLengths[BOIx] << " " << ANGLEBends[BOIx] << " " << TORSIONAngles[BOIx] << std::endl;
		}
	}


}


/*! <!--  --> */
void World::calcSimbodyBAT_TODEL(
	std::vector<std::vector<int>>& ZMatrix,
	std::vector<SimTK::Real>& BONDLengths,
	std::vector<SimTK::Real>& ANGLEBends,
	std::vector<SimTK::Real>& TORSIONAngles)
{

	/*
	SimTK::State& advState = integrator->updAdvancedState();

	bool parFlag = false;
	bool parParFlag = false;
	int childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;

	if(BONDLengths.size() != matter->getNumBodies() -1){
		BONDLengths.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(ANGLEBends.size() != matter->getNumBodies() -1){
		ANGLEBends.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(TORSIONAngles.size() != matter->getNumBodies() -1){
		TORSIONAngles.resize(matter->getNumBodies() - 1, SimTK::NaN);
	}
	if(ZMatrix.size() != matter->getNumBodies() -1){
		ZMatrix.resize(matter->getNumBodies() - 1, std::vector<int>(4, -1));
	}

	for (SimTK::MobilizedBodyIndex childMbx(1); childMbx < matter->getNumBodies(); ++childMbx){
		childIx = int(childMbx);

		const SimTK::MobilizedBody& childMobod = matter->getMobilizedBody(childMbx);
		const SimTK::MobilizedBody& parentMobod = childMobod.getParentMobilizedBody();
		const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
		parentIx = int(parentMbx);

		if(int(childMbx) > 1){
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();
			grandIx = int(grandMbx);
		}

		if(int(childMbx) > 2){
			const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandGrandMbx = grandGrandMobod.getMobilizedBodyIndex();
			grandGrandIx = int(grandGrandMbx);
		}


		// BOND ==============
		const SimTK::Transform& B_X_Fb = childMobod.getInboardFrame(advState);
		const SimTK::Transform& Fb_X_Mb = childMobod.getMobilizerTransform(advState);
		const SimTK::Transform& C_X_Mb = childMobod.getOutboardFrame(advState);

		BONDLengths[int(childMbx) - 1] = C_X_Mb.p().norm(); // correct correct
		
		ZMatrix[int(childMbx) - 1][0] = int(childMbx);
		ZMatrix[int(childMbx) - 1][1] = int(parentMbx);

		if(int(childMbx) > 1){ // ANGLE ========================
			const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
			const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

			const SimTK::Transform& A_X_Fa = parentMobod.getInboardFrame(advState); // A_X_Fa
			const SimTK::Transform& B_X_Ma = parentMobod.getOutboardFrame(advState); // B_X_Ma
			const SimTK::Transform& Fa_X_Ma = parentMobod.getMobilizerTransform(advState); // Fa_X_Ma
			Transform A_X_B = A_X_Fa * Fa_X_Ma * (~B_X_Ma);
			Transform B_X_A = ~A_X_B;

			// STATS
			//Transform B_X_C = B_X_Fb * Fb_X_Mb * (~C_X_Mb);
			// SimTK::Vec4 angAx_BXC = B_X_C.R().convertRotationToAngleAxis(); // [a vx vy vz] angle-axis; -Pi < a <= Pi and |v|=1
			// SimTK::Real ang_BXC = angAx_BXC[0]; // [a] angle
			// Vec3 ax_BXC = Vec3(angAx_BXC[1], angAx_BXC[2], angAx_BXC[3]); // [vx vy vz] axis
			// std::cout << "angle-axis: " << int(childMbx) <<" "<< ang_BXC <<" "<< SimTK::dot(ax_BXC, B_X_C.p()) << std::endl;

			// WORK ==========================================
			Vec3 xAXB_A = A_X_B.R()(0);
			Vec3 xBXFb_B = (B_X_Fb.R())(0);
			Vec3 xBXMa_B = (B_X_Ma.R())(0);
			Vec3 _xBXMa_B = -1 * xBXMa_B;

			ANGLEBends[int(childMbx) - 1] = std::acos(SimTK::dot( xBXFb_B, _xBXMa_B )); // correct correct
			// ================================================

			ZMatrix[int(childMbx) - 1][2] = int(grandMbx);

			if(int(childMbx) > 2){ // TORSION =======================
				const SimTK::MobilizedBody& grandGrandMobod = parentMobod.getParentMobilizedBody().getParentMobilizedBody();

				const SimTK::Transform& T_X_Ft = grandMobod.getInboardFrame(advState); // T_X_Ft
				const SimTK::Transform& Ft_X_Mt = grandMobod.getMobilizerTransform(advState); // Ft_X_Mt
				const SimTK::Transform& A_X_Mt = grandMobod.getOutboardFrame(advState); // A_X_Mt
				const SimTK::Transform& G_X_T = grandGrandMobod.getBodyTransform(advState); // G_X_T
				Transform T_X_A = T_X_Ft * Ft_X_Mt * (~A_X_Mt);

				// WORK ==========================================

				Vec3 xAXFa_A = (A_X_Fa.R())(0);
				Vec3 xAXMt_A = (A_X_Mt.R())(0);
				
				Vec3 v2_B = SimTK::cross(xBXFb_B, xBXMa_B);
				Vec3 v1_A = SimTK::cross(xAXFa_A, xAXMt_A);

				Vec3 v2_B_hat = v2_B.normalize();
				Vec3 v1_A_hat = v1_A.normalize();

				Vec3 v1_B_hat = (~(A_X_B.R())) * v1_A_hat;

				Vec3 v3_B = SimTK::cross(v2_B_hat, v1_B_hat);
				Vec3 v3_B_hat = v3_B.normalize();

				SimTK::Real tors_cos = SimTK::dot(v1_B_hat, v2_B_hat);
				SimTK::Real tors_sin = SimTK::dot(v3_B, _xBXMa_B);

				TORSIONAngles[int(childMbx) - 1] = std::atan2(tors_sin, tors_cos);
				// ==============================================

				ZMatrix[int(childMbx) - 1][3] = int(grandGrandIx);

			}
		}

		childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;
	
	} // _end_ for mbx

	// // Print
	// for (int BOIx = 0; BOIx < BONDLengths.size(); BOIx++){
	// 	std::cout << "ZMatrixBATSimbody:"
	// 		<<" " << ZMatrix[BOIx][0] << " " << ZMatrix[BOIx][1] << " " << ZMatrix[BOIx][2] << " " << ZMatrix[BOIx][3]
	// 		<<" "<< BONDLengths[BOIx] << " " << ANGLEBends[BOIx] << " " << TORSIONAngles[BOIx] << std::endl;
	// }
	*/
}


/*! <!--  --> */
void World::PrintAllTransforms() const
{
	const SimTK::State& advState = integrator->getAdvancedState();

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){

		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		const SimTK::MobilizedBody& parentMobod = mobod.getParentMobilizedBody();
		const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

		const SimTK::Transform& X_GpP = parentMobod.getBodyTransform(advState);
		const SimTK::Transform& X_GP = mobod.getBodyTransform(advState);

		const SimTK::Transform X_PP = (~X_GpP) * X_GP;
		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(advState);

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(advState);

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(advState);

		SimTK::Test::PrintTransform(X_PP, 6, "X_PP", "X_PP:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
		SimTK::Test::PrintTransform(X_PF, 6, "X_PF", "X_PF:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
		SimTK::Test::PrintTransform(X_FM, 6, "X_FM", "X_FM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
		SimTK::Test::PrintTransform(X_BM, 6, "X_BM", "X_BM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
	}
}

/*! <!--  --> */
void World::PrintDefaultTransforms() const
{

	const SimTK::State& advState = integrator->getAdvancedState();

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		std::cout << "mobod " << mbx << std::endl;

		// Get mobod inboard frame X_PF
		const SimTK::Transform& X_PF = mobod.getInboardFrame(advState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(advState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;

		SimTK::Test::PrintTransform(X_PF, 6, "X_PF", "X_PF:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));
		SimTK::Test::PrintTransform(X_BM, 6, "X_BM", "X_BM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));

	}
}

/*! <!--  --> */
void World::PrintXFMs() const
{

	const SimTK::State& advState = integrator->getAdvancedState();

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		std::cout << "mobod " << mbx << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(advState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		SimTK::Test::PrintTransform(X_FM, 6, "X_FM", "X_FM:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(mbx)));

	}
}

/*!
 * <!--  -->
*/
void World::PrintXBMps() const
{
	const SimTK::State& advState = integrator->getAdvancedState();

	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		//std::cout << mbx;

		// Get mobod inboard frame X_BM
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(advState);
		std::cout <<" " << X_BM.p()[0];
	}
}

/*! <!--  --> */
const SimTK::Vector & World::getBMps()
{

	const SimTK::State& advState = integrator->getAdvancedState();

	if(BMps.size() == 0){
		//BMps.resize(matter->getNQ(advState));
		BMps.resize(matter->getNumBodies());
	}

	int bIx = -1;
	for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		//for(int moQIx = 0; moQIx < mobod.getNumQ(advState); moQIx++){
			bIx++;

			const SimTK::Transform& X_BM = mobod.getOutboardFrame(advState);
			BMps[int(mbx)] = X_BM.p()[0];
		//}
		
	}

	return BMps;
}

/*! <!--  --> */
const SimTK::Vector & World::getPFrs()
{

	const SimTK::State& advState = integrator->getAdvancedState();

	if(PFrs.size() == 0){
		//PFrs.resize(matter->getNQ(advState));
		PFrs.resize(matter->getNumBodies());
	}

	int bIx = -1;
	for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		//for(int moQIx = 0; moQIx < mobod.getNumQ(advState); moQIx++){
			bIx++;

			const SimTK::Transform& X_PF = mobod.getInboardFrame(advState);
			PFrs[int(mbx)] = std::acos(X_PF.R()(0)(0));
		//}
		
	}

	return PFrs;
}

/*!
 * <!--  -->
*/
const SimTK::Vector & World::getAdvancedQs()
{
	return matter->getQ(integrator->updAdvancedState());
}

/*!
 * <!--  -->
*/
const void World::PrintAdvancedQs() const
{
	const SimTK::Vector &Qs = matter->getQ(integrator->getAdvancedState());
	std::cout << ownWorldIndex;
	for(int i = 0; i < Qs.size(); i++){
		std::cout << " " << Qs[i];
	}
	std::cout << std::endl;
}



/*!
 * <!--  -->
*/
const SimTK::Vector & World::getAdvancedUs()
{
	return matter->getU(integrator->updAdvancedState());
}

/*!
 * <!--  -->
*/
int World::getNQs(void)
{
	return matter->getNQ(integrator->updAdvancedState());
}

/*!
 * <!--  -->
*/
int World::getNUs(void)
{
	return matter->getNU(integrator->updAdvancedState());
}

// RANDOM_WALK functions
void World::setTopologyIXs(std::vector<int> argTopologyIXs){
	topologyIXs = argTopologyIXs;
}

void World::setAmberAtomIXs(std::vector<std::vector<int>> argAmberAtomIXs){
	amberAtomIXs = argAmberAtomIXs;
}

SimTK::Vec3 
World::getGeometricCenterOfSelection(const SimTK::State & state 
									 //const std::vector<int>& topologyIx, 
									 //const std::vector<std::vector<int>>& amberAtomList
									 )
{
	SimTK_ASSERT_ALWAYS(false, "World::getGeometricCenterOfSelection not implemented yet");
	return SimTK::Vec3(0,0,0);

// 	// return Vec3
// 	SimTK::Vec3 geometricCenter={0,0,0};
// 	// We could just divide by the size of amberAtomList
// 	// but this works even if the user *mistakenly* repeats
// 	// indices
// 	int nOfPoints=0;


// 	// Just a quick check, to skip unnecessary computation in case of 
// 	// user error.
// 	if (amberAtomIXs.size() == 0) {
// 		std::cerr << "Warning: getGeometricCenterOfSelection called with amberAtomList of size 0" << std::endl;
// 		return geometricCenter;
// 	}
	
// 	std::cout << "topologies atoms size " << topologyIXs.size() << " " << topologyIXs.size() << std::endl;

// 	for (int i = 0; i < topologyIXs.size(); i++) {
// 		const auto& topology = topologies[topologyIXs[i]];
// 		const auto& atoms = amberAtomIXs[i];

// 		int amberIx=0;

// 		// Iterate through atoms in said topology and check 	
// 		// if they are in the list
// 		for (auto& atom : topology.subAtomList) {
// 			if (std::find(atoms.begin(), atoms.end(), amberIx) != atoms.end()){
// 				// found
// 				// Get Compound atom index
// 				auto compoundAtomIndex = atom.identity.compoundAtomIndex;
// 				// Get DuMM atom index
// 				const SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(compoundAtomIndex);
// 				// Get Mobilized Body index
// 				const MobilizedBodyIndex mobilizedBodyIndex = forceField->getAtomBody(dAIx);
// 				// Get DuMM Atom Station on its body.
// 				const Vec3 dAS_B = forceField->getAtomStationOnBody(dAIx);
// 				// Re-Express in G
// 				const SimTK::MobilizedBody& mobod_A = matter->getMobilizedBody(mobilizedBodyIndex);
// 				const SimTK::Vec3 dAS_G = mobod_A.findStationLocationInGround(state, dAS_B);
// 				/* const SimTK::Transform& X_GP = mobod_A.getBodyTransform(state);
// 				const SimTK::Vec3 dAS_G = X_GP*dAS_B; */
// 				geometricCenter += dAS_G;
// 				nOfPoints += 1;

// /* 				std::cout << "amberIx: " << amberIx << " dAIx: " << dAIx
// 				<< " MobilizedBodyIndex: " << mobilizedBodyIndex 
// 				<< " dAS_G: " << dAS_G << " nOfPoints: " << nOfPoints
// 				<< std::endl; */
// 			}
			
// 			amberIx += 1;
// 		}
// 	}

// 	// This can probably be done better, but is it clearer?
// 	for(int i=0;i<3;++i)
// 		geometricCenter[i] = geometricCenter[i] / nOfPoints;
// 	std::cout << "geometricCenter : " << geometricCenter << std::endl;

// 	return geometricCenter;

}

/** Nice print helper for get/setAtomsLocations */
void World::PrintAtomsLocations(const std::vector<std::vector<
	std::pair<RoboAtom *, SimTK::Vec3> > >& someAtomsLocations)
{	
	std::cout << "myAtomsLocations[0]" << std::endl;
	for(std::size_t j = 0; j < someAtomsLocations[0].size(); j++){
		int compoundAtomIndex = someAtomsLocations[0][j].first->identity.compoundAtomIndex;
		auto loc = someAtomsLocations[0][j].second;

		printf("%d %.10f %.10f %.10f\n", compoundAtomIndex, loc[0], loc[1], loc[2]);
	}
}

/**
 * Write coordinates to a rst7 file
 * Very costly
*/
void World::WriteRst7FromTopology(std::string FN)

{
	SimTK_ASSERT_ALWAYS(false, "World::WriteRst7FromTopology not implemented yet");

	// updateAtomListsFromSimbody(integrator->updAdvancedState());

	// FILE *File = fopen(FN.c_str(), "w+");

	// int Natoms = 0;
	// for(auto& topology : topologies){
	// 	Natoms += topology.getNumAtoms();
	// }

	// fprintf(File, "TITLE: Created by Robosample with %d atoms\n", Natoms);
	// fprintf(File, "%6d\n", Natoms);

	// int atomCnt = -1;
	// for(auto& topology : topologies){
	// 	for (auto& atom : topology.subAtomList) {			++atomCnt;
	// 		fprintf(File, "%12.7f%12.7f%12.7f", 
	// 			atom.getX() * 10.0, atom.getY() * 10.0, atom.getZ() * 10.0);

	// 		if(atomCnt % 2 == 1){
	// 			fprintf(File, "\n");
	// 		}
	// 	}
	// }

	// if(atomCnt % 2 == 0){
	// 	fprintf(File, "\n");
	// }

	// fflush(File);
	// fclose(File);
}

/** Print transformation geometries */
void World::PrintFullTransformationGeometry(std::string indS, const SimTK::State& someState,
		bool x_pf_r, bool x_fm_r, bool x_bm_r,
		bool x_pf_p, bool x_fm_p, bool x_bm_p)
{
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "X_PF X_FM X_MB\n";
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		const SimTK::Transform& X_PF = mobod.getDefaultInboardFrame();
		const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
		const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
		const SimTK::Transform& X_MB = ~X_BM;

		if( x_pf_r && x_fm_r && x_bm_r && x_pf_p && x_fm_p && x_bm_p){
			std::cout << "mobod " << int(mbx) << std::endl
				<< std::fixed << std::setprecision(3);

			std::stringstream ss;
			ss << indS << int(mbx);
			std::string pref = ss.str();

			printf("%s %9.6f %9.6f %9.6f %9.6f ", pref.c_str(),  X_PF.toMat44()[0][0], X_PF.toMat44()[0][1], X_PF.toMat44()[0][2], X_PF.toMat44()[0][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f ",                X_FM.toMat44()[0][0], X_FM.toMat44()[0][1], X_FM.toMat44()[0][2], X_FM.toMat44()[0][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f\n",               X_MB.toMat44()[0][0], X_MB.toMat44()[0][1], X_MB.toMat44()[0][2], X_MB.toMat44()[0][3]);

			printf("%s %9.6f %9.6f %9.6f %9.6f ", pref.c_str(),  X_PF.toMat44()[1][0], X_PF.toMat44()[1][1], X_PF.toMat44()[1][2], X_PF.toMat44()[1][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f ",                X_FM.toMat44()[1][0], X_FM.toMat44()[1][1], X_FM.toMat44()[1][2], X_FM.toMat44()[1][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f\n",               X_MB.toMat44()[1][0], X_MB.toMat44()[1][1], X_MB.toMat44()[1][2], X_MB.toMat44()[1][3]);

			printf("%s %9.6f %9.6f %9.6f %9.6f ", pref.c_str(),  X_PF.toMat44()[2][0], X_PF.toMat44()[2][1], X_PF.toMat44()[2][2], X_PF.toMat44()[2][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f ",                X_FM.toMat44()[2][0], X_FM.toMat44()[2][1], X_FM.toMat44()[2][2], X_FM.toMat44()[2][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f\n",               X_MB.toMat44()[2][0], X_MB.toMat44()[2][1], X_MB.toMat44()[2][2], X_MB.toMat44()[2][3]);

			printf("%s %9.6f %9.6f %9.6f %9.6f ", pref.c_str(),  X_PF.toMat44()[3][0], X_PF.toMat44()[3][1], X_PF.toMat44()[3][2], X_PF.toMat44()[3][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f ",                X_FM.toMat44()[3][0], X_FM.toMat44()[3][1], X_FM.toMat44()[3][2], X_FM.toMat44()[3][3]);
			printf("   %9.6f %9.6f %9.6f %9.6f\n",               X_MB.toMat44()[3][0], X_MB.toMat44()[3][1], X_MB.toMat44()[3][2], X_MB.toMat44()[3][3]);
			
		}else if(x_pf_p && x_fm_p && x_bm_p){
			std::cout << X_PF.p()[0] << " " << X_PF.p()[1] << " " << X_PF.p()[2] << " ";
			std::cout << X_FM.p()[0] << " " << X_FM.p()[1] << " " << X_FM.p()[2] << " ";
			std::cout << X_BM.p()[0] << " " << X_BM.p()[1] << " " << X_BM.p()[2] << " ";
	
		}else if(x_pf_r){
			PrintMat33(X_PF.R().toMat33(), 3, "X_PF.R");
		}
		else if(x_fm_r){
			PrintMat33(X_FM.R().toMat33(), 3, "X_FM.R");
		}
		else if(x_bm_r){
			PrintMat33(X_BM.R().toMat33(), 3, "X_BM.R");
		}

		else if(x_pf_p){
			std::cout << X_PF.p()[0] << " " << X_PF.p()[1] << " " << X_PF.p()[2] << " ";
		}
		else if(x_fm_p){
			std::cout << X_FM.p()[0] << " " << X_FM.p()[1] << " " << X_FM.p()[2] << " ";
		}
		else if(x_bm_p){
			std::cout << X_BM.p()[0] << " " << X_BM.p()[1] << " " << X_BM.p()[2] << " ";
		}

		std::cout << std::endl;

	}
}

/** Put coordinates into bAtomLists of Topologies.
 * When provided with a State, calcAtomLocationInGroundFrame
 * realizes Position and uses matter to calculate locations **/
void World::updateAtomListsFromSimbody(const SimTK::State &state)
{
	// Iterate through topologies
	for (auto& topology : topologies){

		// Iterate through atoms
		for (auto& atom : topology.updAtoms()) {
			const auto compoundAtomIndex = atom.identity.compoundAtomIndex;
			atom.position = topology.calcAtomLocationInGroundFrameThroughSimbody(compoundAtomIndex, *forceField, *matter, state);
			//std::cout << "updateAtomListsFromCompound (after f_x_m, ix= " << compoundAtomIndex << ") " << atom.getX() << ", " << atom.getY() << ", " << atom.getZ() << std::endl;
		}
	}
}



/**
 * RMSD function
*/
SimTK::Real World::RMSD(
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations
		) const
{

	assert( srcWorldsAtomsLocations.size() == destWorldsAtomsLocations.size() && 
		!("RMSD different size") );

	SimTK::Real rmsdVal = 0;
	int localNofAtoms = 0;

	for(std::size_t i = 0; i < srcWorldsAtomsLocations.size(); i++){
		for(std::size_t j = 0; j < srcWorldsAtomsLocations[i].size(); j++){

			auto srcLoc = srcWorldsAtomsLocations[i][j].second;
			auto destLoc = destWorldsAtomsLocations[i][j].second;

			SimTK::Real sqX = (srcLoc[0] - destLoc[0]);
			sqX = sqX * sqX;
			SimTK::Real sqY = (srcLoc[1] - destLoc[1]);
			sqY = sqY * sqY;
			SimTK::Real sqZ = (srcLoc[2] - destLoc[2]);
			sqZ = sqZ * sqZ;

			rmsdVal += (sqX + sqY + sqZ);

			localNofAtoms++;

		}
	}

	rmsdVal /= SimTK::Real(localNofAtoms);
	rmsdVal = std::sqrt(rmsdVal);

	return rmsdVal;

}

/**
 * Maximum distance between two corresponding atoms
*/
std::pair<int, SimTK::Real> World::maxAtomDeviation(
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<RoboAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations
		) const
{

	assert( srcWorldsAtomsLocations.size() == destWorldsAtomsLocations.size() && 
		!("RMSD different size") );

	SimTK::Real distVal = 0;
	SimTK::Real maxVal = 0;
	int atomCnt = -1;
	int pairIndex = -1;

	for(std::size_t i = 0; i < srcWorldsAtomsLocations.size(); i++){
		for(std::size_t j = 0; j < srcWorldsAtomsLocations[i].size(); j++){

			atomCnt++;

			auto srcLoc = srcWorldsAtomsLocations[i][j].second;
			auto destLoc = destWorldsAtomsLocations[i][j].second;

			SimTK::Real sqX = (srcLoc[0] - destLoc[0]);
			sqX = sqX * sqX;
			SimTK::Real sqY = (srcLoc[1] - destLoc[1]);
			sqY = sqY * sqY;
			SimTK::Real sqZ = (srcLoc[2] - destLoc[2]);
			sqZ = sqZ * sqZ;

			SimTK::Real dist = std::sqrt(sqX + sqY + sqZ);
			
			if(dist > distVal){
				distVal = dist;
				pairIndex = atomCnt;
			}

		}
	}

	return std::pair<int, SimTK::Real>(pairIndex, distVal);

}


/**
 * This function is only intended for root atoms
*/
std::vector<SimTK::Transform>
World::calcMobodToMobodTransforms(
	Topology& topology,
	SimTK::Compound::AtomIndex rootAIx,
	const SimTK::State& someState)
{

	// Get body and parentBody
	// don't know why it works
	//SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(rootAIx);
 	// this should be the correct version
	SimTK::MobilizedBodyIndex mbx =
		topology.getAtomMobilizedBodyIndexThroughDumm(rootAIx, *forceField);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Get the neighbor atom in the parent mobilized body
	SimTK::Compound::AtomIndex chemParentAIx = topology.getChemicalParentOfMobodRootAtom(rootAIx, *matter, *forceField);

	//SimTK::Compound::AtomIndex chemParentAIx =
	//	topology.getNeighbourWithSmallerAIx(rootAIx, forceField);

	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC = topology.getDefaultBondCenterFrameInOtherBondCenterFrame(rootAIx, chemParentAIx);  

	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform(rootAIx);

	// Get Top to parent frame
	const std::pair<int, SimTK::Compound::AtomIndex>& topoAtomPair = getMobodRootAtomIndex(parentMbx);
	SimTK::Compound::AtomIndex parentMobodAIx = topoAtomPair.second;

	//SimTK::Compound::AtomIndex parentRootAIx = mbx2aIx[parentMbx];
	SimTK::Compound::AtomIndex parentRootAIx = parentMobodAIx;

	// Origin of the parent mobod
	SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;
	
	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle =
		topology.bgetDefaultInboardDihedralAngle(rootAIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get inboard bond length
	SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(rootAIx);
	SimTK::Transform InboardLength_mZAxis
		= SimTK::Transform(SimTK::Rotation(), SimTK::Vec3(0, 0, -inboardBondlength));

	// Samuel Flores' terminology
	SimTK::Transform M_X_pin =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

	// Get the old PxFxMxB transform
	SimTK::Transform oldX_PB =
		Proot_X_T * T_X_root
		//* InboardDihedral_XAxis * X_to_Z
		//* InboardDihedral_ZAxis * Z_to_X
		;

	// B_X_Ms
	SimTK::Transform B_X_M = X_to_Z; // aka M_X_pin
	SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
	SimTK::Transform B_X_M_pin 		= X_parentBC_childBC * X_to_Z;
	SimTK::Transform B_X_M_univ 	= X_parentBC_childBC * Y_to_Z;

	// P_X_Fs = old P_X_B * B_X_M
	SimTK::Transform P_X_F 			= oldX_PB * B_X_M;
	SimTK::Transform P_X_F_anglePin = oldX_PB * B_X_M_anglePin;
	SimTK::Transform P_X_F_pin 		= oldX_PB * B_X_M_pin;
	SimTK::Transform P_X_F_univ 	= oldX_PB * B_X_M;

	//SimTK::Transform B_X_M_spheric  = X_parentBC_childBC * X_to_Z * InboardLength_mZAxis;
	//SimTK::Transform B_X_M_spheric  = Transform();
	//SimTK::Transform B_X_M_spheric = B_X_M_pin;
	SimTK::Transform B_X_M_spheric  = X_parentBC_childBC * X_to_Z;
	SimTK::Transform B_X_M_orthospheric = X_parentBC_childBC * X_to_Z;
	
	//SimTK::Transform P_X_F_spheric  = oldX_PB * B_X_M_pin;
	//SimTK::Transform P_X_F_spheric = Transform();
	//SimTK::Transform P_X_F_spheric = P_X_F_pin;
	SimTK::Transform P_X_F_spheric = oldX_PB * B_X_M_pin;
	SimTK::Transform P_X_F_orthospheric = oldX_PB * B_X_M_pin;
	
	// Get mobility (joint type)
	const auto& atom = topology.getAtom(rootAIx);
	SimTK::BondMobility::Mobility mobility;
	RoboBond bond = topology.getBondByGlobalAtomIndex(topology.getGlobalAtomIndex(rootAIx), topology.getGlobalAtomIndex(chemParentAIx));
	mobility = bond.getBondMobility(ownWorldIndex);

	bool anglePin_OR = mobility == SimTK::BondMobility::Mobility::AnglePin ||
					   mobility == SimTK::BondMobility::Mobility::Slider ||
					   mobility == SimTK::BondMobility::Mobility::BendStretch;
	
	if (anglePin_OR && atom.connectivity.neighborsGlobalIndices.size() == 1) {
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};
	} else if (anglePin_OR && atom.connectivity.neighborsGlobalIndices.size() != 1) {
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};
	} else if (mobility == SimTK::BondMobility::Mobility::Torsion 
			|| mobility == SimTK::BondMobility::Mobility::Cylinder) {
		return std::vector<SimTK::Transform> {P_X_F_pin, B_X_M_pin};
	} else if (mobility == SimTK::BondMobility::Mobility::BallM
			|| mobility == SimTK::BondMobility::Mobility::Rigid
			|| mobility == SimTK::BondMobility::Mobility::Translation) { // Cartesian
		return std::vector<SimTK::Transform> {P_X_F, B_X_M};
	} else if (mobility == SimTK::BondMobility::Mobility::Spherical) { // Spherical
		return std::vector<SimTK::Transform> {P_X_F_spheric, B_X_M_spheric};
	} else if (mobility == SimTK::BondMobility::Mobility::OrthoSpherical) { // OrthoSpherical
		return std::vector<SimTK::Transform> {P_X_F_orthospheric, B_X_M_orthospheric};		
	} else {
		std::cout << "Warning: unknown mobility\n";
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};
	}

	//// The Molmodel notation
	// M0 and Mr are actually used for F
	//SimTK::Transform root_X_M0 = InboardDihedral_XAxis; // name used in Molmodel
	//SimTK::Transform T_X_M0 = T_X_root[int(mbx)] * root_X_M0;
	//SimTK::Transform Proot_X_M0 = Proot_X_T * T_X_M0;
	//Transform oldX_PF = Proot_X_M0 * M_X_pin;
	//Transform oldX_BM = M_X_pin;
	//Transform oldX_MB = ~oldX_BM;
	//Transform oldX_FM = InboardDihedral_ZAxis;
	//Transform oldX_PB = oldX_PF * oldX_FM * oldX_MB;
}


/*!
 * <!--  -->
*/
SimTK::BondMobility::Mobility World::determineMobilityFrom_H(SimTK::MobilizedBodyIndex mbx, SimTK::State& someState)
{

	warn("Incomplete.");

	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

	// Use HCol to determine degrees of freedom
	for(SimTK::MobilizerUIndex uIx = SimTK::MobilizerUIndex(0); uIx < mobod.getNumU(someState); uIx++){
		SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(uIx));
		std::cout <<"mbx uIx " << mbx ;
		std::cout <<" " << uIx <<" H_FMCol" 
			<<" " << H_FMCol[0][0] <<" " << H_FMCol[0][1] <<" " << H_FMCol[0][2]
			<<" " << H_FMCol[1][0] <<" " << H_FMCol[1][1] <<" " << H_FMCol[1][2] << std::endl;
	}

	return SimTK::BondMobility::Mobility::Rigid;
}

/*!
 * <!--	 -->
*/
SimTK::Real World::getRootAngle(
	Topology& topology,
	SimTK::Compound::AtomIndex rootAIx,
	const SimTK::State& someState
){
	SimTK::Real bondAngle = SimTK::NaN;

	// Get body and parentBody
	SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(rootAIx, *forceField);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Get the neighbor atom in the parent mobilized body
	SimTK::Compound::AtomIndex chemParentAIx = topology.getChemicalParentOfMobodRootAtom(rootAIx, *matter, *forceField);

	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC = topology.getDefaultBondCenterFrameInOtherBondCenterFrame(rootAIx, chemParentAIx);

	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform(rootAIx);

	// Get Top to parent frame
	const std::pair<int, SimTK::Compound::AtomIndex>& topoAtomPair = getMobodRootAtomIndex(parentMbx);
	SimTK::Compound::AtomIndex parentMobodAIx = topoAtomPair.second;
	SimTK::Compound::AtomIndex parentRootAIx = parentMobodAIx;
	
	SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;

	// chemical parent atom
	SimTK::Transform T_X_chemProot = topology.getTopTransform(chemParentAIx);

	// BEGIN GET ANGLE
	SimTK::Compound::AtomIndex chemGrandParentIx;
	SimTK::Transform T_X_grand;

	if(chemParentAIx > 0){

		chemGrandParentIx = topology.getInboardAtomIndex(chemParentAIx);

		T_X_grand = topology.getTopTransform(chemGrandParentIx);

		SimTK::Vec3 V1 = (~(T_X_root.R())) * T_X_grand.p();
		SimTK::Vec3 V2 = (~(T_X_root.R())) * T_X_chemProot.p();
		SimTK::Vec3 V3 = (~(T_X_root.R())) * T_X_root.p();

		bondAngle = bAngle(V2, V1, V3);

		return bondAngle;
	}

	return bondAngle;

}



/*!
 * <!--	Calc X_FM transforms for reconstruction for root atoms -->
*/
SimTK::Transform World::calcX_FMTransforms(
	Topology& topology,
	SimTK::Compound::AtomIndex rootAIx,
	const SimTK::State& someState)
{

	// X_FM return value
	SimTK::Transform X_FM;

	// Get body and parentBody
	SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(rootAIx, *forceField);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
	// Get the neighbor atom in the parent mobilized body
	SimTK::Compound::AtomIndex chemParentAIx = topology.getChemicalParentOfMobodRootAtom(rootAIx, *matter, *forceField);
	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC = topology.getDefaultBondCenterFrameInOtherBondCenterFrame(rootAIx, chemParentAIx);
	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform(rootAIx);
	// Get Top to parent frame
	const std::pair<int, SimTK::Compound::AtomIndex>& topoAtomPair = getMobodRootAtomIndex(parentMbx);
	SimTK::Compound::AtomIndex parentMobodAIx = topoAtomPair.second;
	SimTK::Compound::AtomIndex parentRootAIx = parentMobodAIx;
	SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;
	// chemical parent atom
	SimTK::Transform T_X_chemProot = topology.getTopTransform(chemParentAIx);
	// BEGIN GET ANGLE
	SimTK::Compound::AtomIndex chemGrandParentIx;
	SimTK::Transform T_X_grand;
	if(chemParentAIx > 0){
		chemGrandParentIx = topology.getInboardAtomIndex(chemParentAIx);
		T_X_grand = topology.getTopTransform(chemGrandParentIx);
		SimTK::Vec3 v1 = (~(T_X_root.R())) * T_X_grand.p();
		SimTK::Vec3 v2 = (~(T_X_root.R())) * T_X_chemProot.p();
		SimTK::Vec3 v3 = (~(T_X_root.R()))  * T_X_root.p();
		SimTK::Real bondAngle = bAngle(v2, v1, v3);
		//if(ownWorldIndex == 1){
			std::cout << "World::calcMobodToMobodTransforms chemGrandParentIx chemParentAIx rootAIx angle "
				<< chemGrandParentIx << " " << chemParentAIx << " " << rootAIx << " " 
				//<< std::endl << T_X_grand << T_X_chemProot << T_X_root << std::endl
				//<< std::endl << topology.getTopTransform(chemGrandParentIx) << topology.getTopTransform(chemParentAIx) << topology.getTopTransform(rootAIx) << std::endl
				//<< "============================="
				//<< std::endl << topology.getTopTransform(Compound::AtomIndex(2)) << topology.getTopTransform(Compound::AtomIndex(4)) << topology.getTopTransform(Compound::AtomIndex(7)) << std::endl
				//<< v1 << " " << v2 << " " << v3 << " "
				<< bondAngle << std::endl;
		//}
	}
	// END GET ANGLE

	// Get the angle of root->parent->grand parent
	//SimTK::Real bondAngle = getRootAngle(topology, rootAIx, someState);

	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle =
		topology.bgetDefaultInboardDihedralAngle(rootAIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get inboard bond length
	SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(rootAIx);
	SimTK::Transform InboardLength_mZAxis
		= SimTK::Transform(SimTK::Rotation(), SimTK::Vec3(0, 0, -inboardBondlength));

	// Samuel Flores' terminology
	SimTK::Transform M_X_pin =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

	// Get the old PxFxMxB transform
	SimTK::Transform oldX_PB =
		Proot_X_T * T_X_root
		* InboardDihedral_XAxis * X_to_Z
		* InboardDihedral_ZAxis * Z_to_X;

	// B_X_Ms
	SimTK::Transform B_X_M = X_to_Z; // aka M_X_pin
	SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
	SimTK::Transform B_X_M_pin 		= X_parentBC_childBC * X_to_Z;
	SimTK::Transform B_X_M_univ 	= X_parentBC_childBC * Y_to_Z;

	// P_X_Fs = old P_X_B * B_X_M
	SimTK::Transform P_X_F 			= oldX_PB * B_X_M;
	SimTK::Transform P_X_F_anglePin = oldX_PB * B_X_M_anglePin;
	SimTK::Transform P_X_F_pin 		= oldX_PB * B_X_M_pin;
	SimTK::Transform P_X_F_univ 	= oldX_PB * B_X_M;

	// Get mobility (joint type)
	const RoboAtom& atom = topology.getAtom(rootAIx);
	SimTK::BondMobility::Mobility mobility;
	RoboBond bond = topology.getBondByGlobalAtomIndex(topology.getGlobalAtomIndex(rootAIx), topology.getGlobalAtomIndex(chemParentAIx));
	mobility = bond.getBondMobility(ownWorldIndex);

	// Convenient bool
		bool anglePin_OR = 
		   mobility == SimTK::BondMobility::Mobility::AnglePin
		|| mobility == SimTK::BondMobility::Mobility::Slider
		|| mobility == SimTK::BondMobility::Mobility::BendStretch;
	
	// Set X_FM value
	if (anglePin_OR && atom.connectivity.neighborsGlobalIndices.size() == 1) {
		X_FM = SimTK::Transform();
		return X_FM;
	} else if (anglePin_OR && atom.connectivity.neighborsGlobalIndices.size() != 1) {
		X_FM = SimTK::Transform();
		return X_FM;
	} else if (mobility == SimTK::BondMobility::Mobility::Torsion 
			|| mobility == SimTK::BondMobility::Mobility::Cylinder) {
		X_FM = SimTK::Transform();
		return X_FM;
	} else if (mobility == SimTK::BondMobility::Mobility::BallM
			|| mobility == SimTK::BondMobility::Mobility::Rigid
			|| mobility == SimTK::BondMobility::Mobility::Translation) { // Cartesian
		X_FM = SimTK::Transform();
		return X_FM;
	} else if (mobility == SimTK::BondMobility::Mobility::Spherical) { // Spherical
		X_FM = SimTK::Transform();
		// X_FM = InboardLength_mZAxis;
		return X_FM;
	} else if (mobility == SimTK::BondMobility::Mobility::OrthoSpherical) { // Spherical
		X_FM = SimTK::Transform();
		// X_FM = InboardLength_mZAxis;
		return X_FM;		
	} else {
		X_FM = SimTK::Transform();
		std::cout << "Warning: unknown mobility\n";
		return X_FM;
	}
}

/** Set up Fixman torque **/
void World::addFixmanTorque()
{
	// Set flag
	assert(!isUsingFixmanTorque());
	useFixmanTorque = true;

	// Alloc memory for FixmanTorque implementation and add to forces
	FixmanTorqueImpl = new FixmanTorque(matter.get());
	FixmanTorqueForce = std::make_unique<SimTK::Force::Custom>(*forces, FixmanTorqueImpl);

	FixmanTorqueExtImpl = new FixmanTorqueExt(matter.get());
	FixmanTorqueExtForce = std::make_unique<SimTK::Force::Custom>(*forces, FixmanTorqueExtImpl);

	// FixmanTorqueImpl = new FixmanTorque(matter->get());			//	
	// FixmanTorqueForce = new Force::Custom(forces, FixmanTorqueImpl);			//
	// FixmanTorqueExtImpl = new FixmanTorqueExt(matter->get());	//
	// FixmanTorqueExtForce = new Force::Custom(forces, FixmanTorqueExtImpl);		//

	// for (int i = 0; i < 10; i++) {
	// 	controller.push_back(std::make_unique<SimTK::ConformationalController>(forces, matter, SimTK::MobilizedBodyIndex(i), SimTK::Vec3(0,0,0)));
	// 	controlForce.push_back(std::make_unique<SimTK::Force::Custom>(forces, controller.back().get()));
	// }
}

/** Check if the Fixman torque flag is set **/
bool World::isUsingFixmanTorque() const
{
	return useFixmanTorque;
}

/** Get writble pointer to FixmanTorque implementation **/
FixmanTorque * World::updFixmanTorque()
{
	assert(isUsingFixmanTorque());
	// return FixmanTorqueImpl.get();
	return FixmanTorqueImpl;
}

/** Get pointer to FixmanTorque implementation **/
FixmanTorque * World::getFixmanTorque() const
{
	assert(isUsingFixmanTorque());
	// return FixmanTorqueImpl.get();
	return FixmanTorqueImpl;
}

// ----------------------
// --- Thermodynamics ---
// ----------------------

/** Get the World temperature **/
SimTK::Real World::getTemperature()
{
	return this->temperature;
}

/** Set this World temperature but also ths samplers and
Fixman torque temperature. **/
void World::setTemperature(SimTK::Real argTemperature)
{
	// Set the temperature for this World
	this->temperature = argTemperature;

	// Set the boost temperature for the samplers
	for (auto& sampler: samplers) {
		sampler->setTemperature(argTemperature);
	}

	// Set the temperature for the Fixman torque also
	if(useFixmanTorque){
		FixmanTorqueImpl->setTemperature(this->temperature);
		FixmanTorqueExtImpl->setTemperature(this->temperature);
	}
}
//...............

/** Set this World temperature but also ths samplers and
Fixman torque temperature. **/
void World::setBoostTemperature(SimTK::Real argTemperature)
{
	// Set the boost temperature for the samplers
	for (auto& sampler: samplers) {
		sampler->setBoostTemperature(argTemperature);
	}
}
//...............

//...................
// --- Simulation ---
//...................

/** Get/Set seed for reproducibility. **/
void World::setSeed(uint32_t argSeed)
{
	randomEngine = buildRandom32(argSeed);
	// forceField->setOpenMMseed(randomEngine());
}


/** Amber like scale factors. **/
void World::setAmberForceFieldScaleFactors()
{
	forceField->setVdw12ScaleFactor(0.0);
	forceField->setVdw13ScaleFactor(0.0);
	forceField->setVdw14ScaleFactor(0.5); // RESTORE from OpenMM
	//forceField->setVdw14ScaleFactor(0.0); // for OpenMM
	forceField->setVdw15ScaleFactor(1.0);

	//* RESTORE SAFETY
	forceField->setCoulomb12ScaleFactor(0.0);
	forceField->setCoulomb13ScaleFactor(0.0);
	forceField->setCoulomb14ScaleFactor(0.8333333333); // RESTORE from OpenMM
	//forceField->setCoulomb14ScaleFactor(0.0); // for OpenMM
	forceField->setCoulomb15ScaleFactor(1.0);
	//forceField->setVdwMixingRule(
	//       SimTK::DuMMForceFieldSubsystem::LorentzBerthelot); */

 	/* DANGER ! no electrostatics
	forceField->setCoulomb12ScaleFactor(0.0);
	forceField->setCoulomb13ScaleFactor(0.0);
	forceField->setCoulomb14ScaleFactor(0.0);
	forceField->setCoulomb15ScaleFactor(0.0);
	// */
}

/** Set a global scaling factor for all the terms in the forcefield **/
void World::setGlobalForceFieldScaleFactor(SimTK::Real scaleFactor)
{
	forceField->setBondStretchGlobalScaleFactor(scaleFactor);
	forceField->setBondBendGlobalScaleFactor(scaleFactor);
	forceField->setBondTorsionGlobalScaleFactor(scaleFactor);
	forceField->setAmberImproperTorsionGlobalScaleFactor(scaleFactor);

	forceField->setVdw12ScaleFactor(scaleFactor);
	forceField->setVdw13ScaleFactor(scaleFactor);
	forceField->setVdw14ScaleFactor(scaleFactor);
	forceField->setVdw15ScaleFactor(scaleFactor);
	forceField->setVdwGlobalScaleFactor(scaleFactor);

	forceField->setCoulomb12ScaleFactor(scaleFactor);
	forceField->setCoulomb13ScaleFactor(scaleFactor);
	forceField->setCoulomb14ScaleFactor(scaleFactor);
	forceField->setCoulomb15ScaleFactor(scaleFactor);
	forceField->setCoulombGlobalScaleFactor(scaleFactor);


//	std::cout << "GLOBAL SCALE FACTORS SET TO 0.0\n";
//	forceField->setBondTorsionGlobalScaleFactor(0);
//	forceField->setAmberImproperTorsionGlobalScaleFactor(0);
//
//	forceField->setVdw12ScaleFactor(0);
//	forceField->setVdw13ScaleFactor(0);
//	forceField->setVdw14ScaleFactor(0);
//	forceField->setVdw15ScaleFactor(0);
//	forceField->setVdwGlobalScaleFactor(0);
//
//	forceField->setCoulomb12ScaleFactor(0);
//	forceField->setCoulomb13ScaleFactor(0);
//	forceField->setCoulomb14ScaleFactor(0);
//	forceField->setCoulomb15ScaleFactor(0);
//	forceField->setCoulombGlobalScaleFactor(0);

}

/** Set GBSA implicit solvent scale factor. **/
void World::setGbsaGlobalScaleFactor(SimTK::Real scaleFactor)
{
	forceField->setGbsaGlobalScaleFactor(scaleFactor);
}

//...................
// --- Statistics ---
//...................

/** How many samples do we have so far **/
std::size_t World::getNofSamples() const
{
	// Zero it every time the user asks
	std::size_t nofSamples = 0;

	// Gather samples from all the samplers
	for(size_t i = 0; i < samplers.size(); i++){
		nofSamples += (samplers[i])->getNofSamples();
	}

	return nofSamples;
}

/** How many Samplers does this World have. **/
std::size_t World::getNofSamplers() const {
	return samplers.size();
}

/*! <!-- Add a sampler to this World using the specialized struct
 * for samplers names. -->
*/
bool World::addSampler(SamplerName samplerName,
	IntegratorType integratorType,
	ThermostatName thermostatName,
	bool useFixmanPotential)
{
	// // Non-bonded forces will always be calculated with OpenMM regardless of the integrator type
	// // However, if the integrator is OMMVV, we want more that that so we set it to false
	// // forceField->setUseOpenMMCalcOnlyNonBonded(integratorType != IntegratorType::OMMVV);
	// forceField->setUseOpenMMCalcOnlyNonBonded(false);

	if(samplerName == SamplerName::HMC) {

		// Construct a new sampler
		samplers.emplace_back(std::make_unique<HMCSampler>(*this, *compoundSystem, *matter, topologies, *forceField, *forces, *timeStepper));

		// Set sampler parameters
		samplers.back()->setIntegratorType(integratorType);
		samplers.back()->setThermostat(thermostatName);
		samplers.back()->setSeed(randomEngine);

		// Initialize the sampler
		samplers.back()->initialize();

		// TODO should this be inherited from parent world?
		if (useFixmanPotential) {
			samplers.back()->useFixmanPotential();
		}
	} else {
		SimTK_ASSERT_ALWAYS(false, "World::addSampler(): sampler name not recognized.");
		return false;
	}

	return true;
}

void World::useOpenMM(bool ommvv, SimTK::Real boostTemp, SimTK::Real timestep) {
	forceField->setUseOpenMMAcceleration(true);

	if (ommvv) {
		forceField->setUseOpenMMIntegration(true);
		forceField->setUseOpenMMCalcOnlyNonBonded(false);
		forceField->setDuMMTemperature(boostTemp);
		forceField->setDuMMTimestep(timestep);
	} else {
		forceField->setUseOpenMMCalcOnlyNonBonded(false);
	}
}

// Get a sampler based on its position in the samplers vector
// TODO Use ampler polymorphism (was const BaseSampler *)
BaseSampler * World::getSampler(std::size_t which) const
{
	return samplers[which].get();
}

// Get a writable sampler based on its position in the samplers vector
// TODO Use Sampler polymorphism
BaseSampler * World::updSampler(std::size_t which)
{
	return samplers[which].get();
}

/** Get a const reference to a molecule **/
const Topology& World::getTopology(std::size_t moleculeNumber) const{
	return topologies[moleculeNumber];
}

/** Get a writble reference to the last molecule. **/
Topology& World::updTopology(std::size_t moleculeNumber){
	//return topologies.back();
	return topologies[moleculeNumber];
}

// Calculate Fixman potential
SimTK::Real World::calcFixman(void)
{
    SimTK::State& currentAdvancedState = integrator->updAdvancedState();
    updateAtomListsFromSimbody(currentAdvancedState); // for det(MBAT)
	SimTK::Real Fixman = updSampler(0)->calcFixman(currentAdvancedState);
	return Fixman;
}

SimTK::Real World::findDecorrelationTime(const SimTK::State& state, int equilSteps, int tuneSteps, SimTK::Real timestep) {

	std::vector<SimTK::Vec3> initialPositions(numAtoms), positions(numAtoms);
	
	// Set initial positions
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (const auto& atom : topologies[topoIx].getAtoms()) {
			const SimTK::Vec3 location = atomTargetLocaltionsCache[topoIx][atom.identity.compoundAtomIndex];
			initialPositions[atom.identity.globalIndex] = location;
		}
	}

	// Run some steps to equilibrate and get initial positions
	if (equilSteps < 1) {
		throw std::invalid_argument("World::tuneOpenMM(): equilSteps must be at least 1");
	}
	OPENMM::get().integrateTrajectory(initialPositions, positions, true, equilSteps, timestep);

	// Allocate memory
	const int numRows = zMatrix.get().size() - 1;
	const int numAngles = zMatrix.get().size() - 2;
	const int numCoords = numRows + numAngles;

	std::vector<std::vector<SimTK::Real>> coordinates(numCoords);
	coordinates.resize(numCoords);
	for (auto& ts : coordinates) {
		ts.reserve(tuneSteps);
	}

	// Run the trajectory and collect coordinate time series
	for (int i = 0; i < tuneSteps; i++) {
		OPENMM::get().integrateTrajectory(positions, positions, false, 1, timestep);

		int coordinateIndex = 0;

		for (const auto& row : zMatrix.get()) {
			if (row.globalIndices[0] != -1 && row.globalIndices[1] != -1) {
				const SimTK::Vec3& pos1 = positions[row.globalIndices[0]];
				const SimTK::Vec3& pos2 = positions[row.globalIndices[1]];
				SimTK::Real bondLength = (pos1 - pos2).norm();

				coordinates[coordinateIndex].push_back(bondLength);
				coordinateIndex++;
			}

			if (row.globalIndices[0] != -1 && row.globalIndices[1] != -1 && row.globalIndices[2] != -1) {
				const SimTK::Vec3& pos1 = positions[row.globalIndices[0]];
				const SimTK::Vec3& pos2 = positions[row.globalIndices[1]];
				const SimTK::Vec3& pos3 = positions[row.globalIndices[2]];
				SimTK::Real angleInRad = bAngle(pos1, pos2, pos3);

				coordinates[coordinateIndex].push_back(angleInRad);
				coordinateIndex++;
			}
		}
	}

	// Save scale
	std::vector<SimTK::Real> scale;
	scale.reserve(numCoords);

	for (const auto& row : zMatrix.get()) {
		if (row.globalIndices[0] != -1 && row.globalIndices[1] != -1) {
			const auto& bond = topologies[row.moleculeIndex].getBondByCompoundAtomIndex(row.compoundAtomIndices[0], row.compoundAtomIndices[1]);
			scale.push_back(std::sqrt(bond.stiffnessInKJPerNmSq));
		}

		if (row.globalIndices[0] != -1 && row.globalIndices[1] != -1 && row.globalIndices[2] != -1) {
			const auto& angle = topologies[row.moleculeIndex].getAngleByCompoundAtomIndex(row.compoundAtomIndices[0], row.compoundAtomIndices[1], row.compoundAtomIndices[2]);
			scale.push_back(std::sqrt(angle.stiffnessInKJPerRadSq));
		}
	}

	// Standardize each coordinate time series to have mean 0 and stddev 1 for better numerical stability in ACF calculation
	optimizeCoordinates(coordinates, scale);

	// Find maximum decorrelation time
	SimTK::Real tauMax = 0.0;
	for (const auto& ts : coordinates) {
		SimTK::Real tau = integratedAutocorrelation(ts);
		if (tau > tauMax) tauMax = tau;
	}

	return tauMax;
}

void World::optimizeCoordinates(std::vector<std::vector<SimTK::Real>>& coordinates, const std::vector<SimTK::Real>& scale) const {
    for (size_t coordIdx = 0; coordIdx < coordinates.size(); ++coordIdx) {
        auto& ts = coordinates[coordIdx];
        const size_t N = ts.size();
        if (N < 2) continue;

        const SimTK::Real n_f = static_cast<SimTK::Real>(N);

        // Recompute index sums using closed-form equations
        const SimTK::Real sumX = n_f * (n_f - 1.0) * 0.5;
        const SimTK::Real sumXX = n_f * (n_f - 1.0) * (2.0 * n_f - 1.0) / 6.0;
        const SimTK::Real det = n_f * sumXX - sumX * sumX;

        // Pass 1: Compute sumY and sumXY for linear regression
        SimTK::Real sumY = 0.0;
        SimTK::Real sumXY = 0.0;
        
        #pragma omp simd reduction(+:sumY, sumXY)
        for (size_t i = 0; i < N; ++i) {
            sumY += ts[i];
            sumXY += static_cast<SimTK::Real>(i) * ts[i];
        }

        // Calculate regression coefficients
        const SimTK::Real slope = (n_f * sumXY - sumX * sumY) / det;
        const SimTK::Real intercept = (sumY - slope * sumX) / n_f;
        const SimTK::Real k = scale[coordIdx];

        // Pass 2: Apply detrending + scaling AND compute new variance in one go
        // Mathematically: ts_new = (ts_old - (intercept + slope*i)) * k
        // We need the mean of ts_new to standardize. 
        // Note: Detrending technically makes the mean 0 automatically.
        
        SimTK::Real m2 = 0.0; // sum of squares for variance
        #pragma omp simd reduction(+:m2)
        for (size_t i = 0; i < N; ++i) {
            SimTK::Real val = (ts[i] - (intercept + slope * i)) * k;
            ts[i] = val;
            m2 += val * val; 
        }

        // Final Pass: Standardize (Mean is 0 due to detrending)
        SimTK::Real stddev = std::sqrt(m2 / n_f);
        if (stddev > 1e-12) {
            SimTK::Real invStd = 1.0 / stddev;
            #pragma omp simd
            for (size_t i = 0; i < N; ++i) {
                ts[i] *= invStd;
            }
        }
    }
}

SimTK::Real World::integratedAutocorrelation(const std::vector<SimTK::Real>& x, int maxLag) const {
    const size_t N = x.size();
    if (N < 2) return 0.0;

    if (maxLag < 0) maxLag = std::min<int>(1000, static_cast<int>(N / 2));

    SimTK::Real mean = std::accumulate(x.begin(), x.end(), 0.0) / N;

    // Pre-center data for SIMD
    // This allows the inner loop to be a simple multiply-accumulate
    std::vector<SimTK::Real> centered(N);
    SimTK::Real var_sum = 0.0;
    
    #pragma omp simd reduction(+:var_sum)
    for (size_t i = 0; i < N; ++i) {
        centered[i] = x[i] - mean;
        var_sum += centered[i] * centered[i];
    }

    SimTK::Real var = var_sum / N;
    if (var < 1e-15) return 0.0;

    // Compute autocorrelation
    SimTK::Real tau = 0.5;
    
    for (int lag = 1; lag <= maxLag; ++lag) {
        SimTK::Real c = 0.0;
        const SimTK::Real* p_base = centered.data();
        const SimTK::Real* p_shift = centered.data() + lag;
        const size_t limit = N - lag;

        // Hinting to the compiler to vectorize this dot product
        #pragma omp simd reduction(+:c)
        for (size_t i = 0; i < limit; ++i) {
            c += p_base[i] * p_shift[i];
        }
        
        tau += (c / N) / var;
    }

    return tau;
}

bool World::generateSamples(int howManySamplesPerRound, std::stringstream& worldOutStream, const std::string& header, bool verbose)
{
	SimTK::State& state = integrator->updAdvancedState();
	updSampler(0)->reinitialize(state, worldOutStream, verbose);

	bool validated = true;

	if (getSampler(0)->getIntegratorType() == IntegratorType::OMMVV) {
		
		// if (!tuned) {
		// 	std::cout << "[OpenMM tune]: Estimating decorrelation time across coordinates to set random step parameters for the OMMVV sampler..." << std::endl;

		// 	// const SimTK::Real tau = findDecorrelationTime(state, 10'000, 100'000, 0.002);
		// 	// int minSteps = std::round(5 * tau);
		// 	// int maxSteps = std::round(10 * tau);

		// 	int minSteps = 242;
		// 	int maxSteps = 484;

		// 	samplers[0]->setCartesianRandomSteps(minSteps, maxSteps);

		// 	// std::cout << "[OpenMM tune]: Estimated max autocorrelation time (tau) across coordinates: " << tau << " steps\n";
		// 	std::cout << "[OpenMM tune]: L_min: " << minSteps << " steps" << std::endl;
		// 	std::cout << "[OpenMM tune]: L_max: " << maxSteps << " steps" << std::endl;

		// 	// const qualifier prevents us from reusing the same positions vector, so we need to copy it here
		// 	std::vector<SimTK::Vec3> positions(numAtoms);
		// 	for (const auto& pos : OPENMM::get().getPositions()) {
		// 		positions.push_back(pos);
		// 	}
		// 	const bool resetPositions = false;

		// 	std::uniform_real_distribution<SimTK::Real> uniformRealDistribution(0.0, 1.0);

		// 	while (true) {
		// 		const int numAttempts = 50;
		// 		int numAccepted = 0;
		// 		std::uniform_int_distribution<> uniformIntDistribution(minSteps, maxSteps);

		// 		for (int i = 0; i < numAttempts; i++) {
		// 			const SimTK::Real oldE = OPENMM::get().getPotentialEnergy() + OPENMM::get().getKineticEnergy();
		// 			OPENMM::get().integrateTrajectory(positions, positions, resetPositions, uniformIntDistribution(randomEngine), 0.002);
		// 			const SimTK::Real newE = OPENMM::get().getPotentialEnergy() + OPENMM::get().getKineticEnergy();
		// 			const SimTK::Real deltaE = newE - oldE;
		// 			const SimTK::Real beta = 1.0 / (temperature * SimTK_BOLTZMANN_CONSTANT_MD);
		// 			const SimTK::Real metropolisCriterion = (deltaE < 0) ? 1.0 : std::exp(-1.0 * beta * deltaE);
					
		// 			const SimTK::Real randomReal = uniformRealDistribution(randomEngine);
		// 			const bool accepted = (randomReal < metropolisCriterion);
		// 			if (accepted) {
		// 				numAccepted++;
		// 			}
		// 		}

		// 		const SimTK::Real acceptanceRate = static_cast<SimTK::Real>(numAccepted) / static_cast<SimTK::Real>(numAttempts);
		// 		std::cout << "[OpenMM tune]: Acceptance rate: " << acceptanceRate << " | L_min: " << minSteps << " | L_max: " << maxSteps << std::endl;
				
		// 		if (acceptanceRate > 0.8) {
		// 			minSteps = std::round(minSteps * 1.2);
		// 			maxSteps = std::round(maxSteps * 1.2);
		// 		}
		// 		if (acceptanceRate < 0.35) {
		// 			minSteps = std::round(minSteps * 0.8);
		// 			maxSteps = std::round(maxSteps * 0.8);
		// 		}

		// 	}

		// 	tuned = true;
		// }

		// OpenMM does cartesian integration, so no locking here
		for (int sampleIx = 0; sampleIx < howManySamplesPerRound; ++sampleIx) {
			// TODO do we reinitialize() here?
			validated = updSampler(0)->sample_iteration(state, worldOutStream, verbose) && validated;
		}
	} else {
		// // Simbody supports locking mobilizers
		// for (const auto& mobodLock : mobodLocks) {
		// 	for (const auto mbx : mobodLock) {
		// 		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		// 		mobod.lock(state);
		// 	}	

		// 	// compoundSystem->realize(state, SimTK::Stage::Position);

		// 	// Sample
		// 	for (int sampleIx = 0; sampleIx < howManySamplesPerRound; ++sampleIx) {
		// 		validated &= updSampler(0)->sample_iteration(state, worldOutStream, verbose);
		// 		if (!validated) {
		// 			continue;
		// 		}
		// 	}

		// 	if (!validated) {
		// 		std::cout << "\tWorld " << ownWorldIndex << " sample rejected during sampling mobilizer " << std::endl;
		// 		break;
		// 	}

		// 	for (const auto mbx : mobodLock) {
		// 		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		// 		mobod.unlock(state);
		// 	}	
		// }

		// TODO the above will lock literally everything, so it won't simulate anything
		for (int sampleIx = 0; sampleIx < howManySamplesPerRound; ++sampleIx) {
			validated &= updSampler(0)->sample_iteration(state, worldOutStream, verbose);
			if (!validated) {
				continue;
			}
		}
	}

	// Update atom target locations cache after sampling
	for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
		for (SimTK::Compound::AtomIndex cAIx = SimTK::Compound::AtomIndex(0); cAIx < topologies[topoIx].getAtoms().size(); cAIx++) {

			// Get location in ground frame (Cartesian coordinates)
			const SimTK::Vec3 location = topologies[topoIx].calcAtomLocationInGroundFrameThroughSimbody(cAIx, *forceField, *matter, state);

			// Update cache
			atomTargetLocaltionsCache[topoIx][cAIx] = location;
		}
	}

	if (testing && updSampler(0)->getIntegratorType() == IntegratorType::OMMVV) {
		const auto& positions = OPENMM::get().getPositions();
		std::vector<SimTK::Compound::AtomTargetLocations> atomTargetLocationsFromOpenMM (topologies.size());
		
		for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
			for (const auto& atom : topologies[topoIx].getAtoms()) {
				const SimTK::Vec3 locationFromOpenMM = positions[atom.identity.globalIndex];
				atomTargetLocationsFromOpenMM[topoIx][atom.identity.compoundAtomIndex] = locationFromOpenMM;
			}
		}

		const SimTK::State& state = compoundSystem->realizeTopology();
		compoundSystem->realize(state, SimTK::Stage::Position);

		// We ignore match residuals for now
		// They only work if you need to match positions to the Topology list
		// This is only needed when transfering from one world to another
		coordinateTransferErrors.push_back(checkCoordinateTransfer(atomTargetLocationsFromOpenMM));

		// Calculate RMSD after sampling if the sample is accepted
		SimTK::Real sumSquaredDist = 0.0;
		for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
			for (SimTK::Compound::AtomIndex cAIx = SimTK::Compound::AtomIndex(0); cAIx < topologies[topoIx].getAtoms().size(); cAIx++) {
				const SimTK::Vec3& newLocation = atomTargetLocationsFromOpenMM[topoIx][cAIx];
				const SimTK::Vec3& oldLocation = atomTargetLocaltionsCacheOld[topoIx][cAIx];

				sumSquaredDist += (newLocation - oldLocation).normSqr();
			}
		}
		SimTK::Real rmsd = std::sqrt(sumSquaredDist / numAtoms);
		acceptanceRMSD.emplace_back(std::make_pair(validated, rmsd));
	}

    return validated;
}

/** Print information about Simbody systems. For debugging purpose. **/
void World::PrintSimbodyStateCache(SimTK::State& someState){
	std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
	for(int i = 0; i < someState.getNumSubsystems(); i++){
		std::cout << " Subsystem " << i
			<< " Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
			<< " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
			<< " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i))
			<< std::endl;
	}
}


// void World::initializeTaskSpace(SimTK::CompoundSystem &compoundSystem, SimTK::GeneralForceSubsystem& force, SimTK::SimbodyMatterSubsystem& matter) {
	
// void World::getLocationsForTaskSpace() {
// }


void World::setSamplesPerRound(int samples) {
	samplesPerRound = samples;
}

int World::getSamplesPerRound() const {
	return samplesPerRound;
}

// void World::setDistortOption(int distort) {
// 	distortOption = distort;
// }

// int World::getDistortOption() const {
// 	return distortOption;
// }

void World::setRootMobility(ROOT_MOBILITY rootMobility) {
	switch (rootMobility)
	{
	case ROOT_MOBILITY::FREE:
		rootMobilizer = "Free";
		break;
	case ROOT_MOBILITY::CARTESIAN:
		rootMobilizer = "Cartesian";
		break;
	case ROOT_MOBILITY::WELD:
		rootMobilizer = "Weld";
		break;
	case ROOT_MOBILITY::FREE_LINE:
		rootMobilizer = "FreeLine";
		break;
	case ROOT_MOBILITY::BALL:
		rootMobilizer = "Ball";
		break;
	case ROOT_MOBILITY::PIN:
		rootMobilizer = "Pin";
		break;
	default:
		break;
	}
}

const SimTK::String& World::getRootMobility() const {
	return rootMobilizer;
}

/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_bon(){return forceField->getEnergies_drl_bon();}
/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_ang(){return forceField->getEnergies_drl_ang();}
/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_tor(){return forceField->getEnergies_drl_tor();}
/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_n14(){return forceField->getEnergies_drl_n14();}
/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_vdw(){return forceField->getEnergies_drl_vdw();}
/*!
 * <!-- Drill -->
*/
const std::vector<std::vector<SimTK::Real>>& World::getEnergies_drl_cou(){return forceField->getEnergies_drl_cou();}
/*!
 * <!-- Drill -->
*/
const std::vector<OpenMM::Vec3>& World::getForces_drl_bon(){return forceField->getForces_drl_bon();}
/*!
 * <!-- Drill -->
*/
const std::vector<OpenMM::Vec3>& World::getForces_drl_ang(){return forceField->getForces_drl_ang();}
/*!
 * <!-- Drill -->
*/
const std::vector<OpenMM::Vec3>& World::getForces_drl_tor(){return forceField->getForces_drl_tor();}
/*!
 * <!-- Drill -->
*/
const std::vector<OpenMM::Vec3>& World::getForces_drl_n14(){return forceField->getForces_drl_n14();}

/*!
 * <!--  -->
*/
void World::printDrilling(void)
{

#ifdef __DRILLING__

	// for (DuMM::NonbondAtomIndex nax(0); nax < forceField->getNumNonbondAtoms(); ++nax) {
	// 	//const DuMM::DuMMAtom&        dummAtom = forceField->getAtom(forceField->getAtomIndexOfNonbondAtom(nax));
	// 	const SimTK::DuMM::AtomIndex dax = forceField->getAtomIndexOfNonbondAtom(nax);
	// 	//const DuMM::IncludedAtomIndex& iax = dummAtom.getIncludedAtomIndex();
	// 	std::cout << "drl World::newFunction dax nax"
	// 		<< " " << dax << " " << nax //<< " " << iax 
	// 		<< std::endl;
	// }

	const std::vector<std::vector<SimTK::Real>>& drl_bon_Energies = forceField->getEnergies_drl_bon();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World bonE");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_bon_Energies[fIx][fJx]);
		}
		printf("\n");
	}
	const std::vector<std::vector<SimTK::Real>>& drl_ang_Energies = forceField->getEnergies_drl_ang();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World angE");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_ang_Energies[fIx][fJx]);
		}
		printf("\n");
	}        
	const std::vector<std::vector<SimTK::Real>>& drl_tor_Energies = forceField->getEnergies_drl_tor();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World torE");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_tor_Energies[fIx][fJx]);
		}
		printf("\n");
	}
	const std::vector<std::vector<SimTK::Real>>& drl_n14_Energies = forceField->getEnergies_drl_n14();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World n14E");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_n14_Energies[fIx][fJx]);
		}
		printf("\n");
	}             
	const std::vector<std::vector<SimTK::Real>>& drl_vdw_Energies = forceField->getEnergies_drl_vdw();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World vdwE");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_vdw_Energies[fIx][fJx]);
		}
		printf("\n");
	}             
	const std::vector<std::vector<SimTK::Real>>& drl_cou_Energies = forceField->getEnergies_drl_cou();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		printf("drl World couE");
		for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx){
			printf(" %f", drl_cou_Energies[fIx][fJx]);
		}
		printf("\n");
	}             

	const std::vector<OpenMM::Vec3>& drl_bon_Forces = forceField->getForces_drl_bon();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		const OpenMM::Vec3& ommForce = drl_bon_Forces[fIx];
		const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
		printf("drl World bonF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
	}
	const std::vector<OpenMM::Vec3>& drl_ang_Forces = forceField->getForces_drl_ang();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		const OpenMM::Vec3& ommForce = drl_ang_Forces[fIx];
		const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
		printf("drl World angF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
	}
	const std::vector<OpenMM::Vec3>& drl_tor_Forces = forceField->getForces_drl_tor();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		const OpenMM::Vec3& ommForce = drl_tor_Forces[fIx];
		const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
		printf("drl World torF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
	}
	const std::vector<OpenMM::Vec3>& drl_n14_Forces = forceField->getForces_drl_n14();
	printf("drl World::newFunction\n");
	for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx){
		const OpenMM::Vec3& ommForce = drl_n14_Forces[fIx];
		const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
		printf("drl OMMPlug n14F %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
	}

#endif // __DRILLING__ 

}

//////////////////////////////////
/////      Z Matrix BAT      /////
//////////////////////////////////
/*!
 * <!--	zmatrixbat_ -->
*/
// void World::setZMatrixBATValue(size_t rowIndex, size_t colIndex, SimTK::Real value) {
// 	// Set the value at the specified position
// 	zMatrixBAT[rowIndex][colIndex] = value;
// }


/*!
 * <!-- zmatrixbat_ -->
*/
// void World::calcZMatrixBAT(SimTK::State& someState)
// {
// 	assert(!"Not implemented");
// }

//////////////////////////////////
/////      Z Matrix BAT      /////
//////////////////////////////////
