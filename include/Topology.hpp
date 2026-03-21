#pragma once

#include "TopologyElements.hpp"
#include "OpenMM.hpp"
#include <cstddef>
#include <unordered_map>

using CompoundAtomIndexPair = std::pair<SimTK::Compound::AtomIndex, SimTK::Compound::AtomIndex>;

inline CompoundAtomIndexPair canonical(SimTK::Compound::AtomIndex a, SimTK::Compound::AtomIndex b) {
	return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

/** Topological information (bonds graph) for one molecule.
It maps to one compound in Molmodel thus it is derived 
from Molmodel Compound class.
It does the following things:
   - loads information from input files such as Amber input prmtop / inpcrd
   - adds parameters to a DuMM force field which belongs to the World class
	 because one DuMM class should be used for multiple molecules
   - contructs the graph based on a list of Atom objects each of 
	 which already contains bonding information from the input files
   - defines the rigid bodies based on imput files provided by the users.
Contains a list of atoms bAtomList which consists of Atom 
objects **/
class Topology : public SimTK::Compound {
public:

	Topology(const SimTK::Compound::Name& name, SimTK::CompoundSystem::CompoundIndex compoundIndex, int rootGlobalAtomIx);

	~Topology() override = default;

	void setAtoms(Span<RoboAtom> atoms) { subAtomList = atoms; }
	const Span<RoboAtom> getAtoms() const { return subAtomList; }
	Span<RoboAtom> updAtoms() { return subAtomList; }

	void setBonds(Span<RoboBond> bonds) {
		subBondList = bonds;

		for (std::size_t i = 0; i < subBondList.size(); ++i) {
			const auto& bond = subBondList[i];
			const auto canonicalBond = canonicalizeBond(bond.compoundAtomIndices[0], bond.compoundAtomIndices[1]);
			const auto cAIx0 = canonicalBond.first;
			const auto cAIx1 = canonicalBond.second;

			const std::string name = subAtomList[cAIx0].identity.uniqueAtomName + "-" + subAtomList[cAIx1].identity.uniqueAtomName;
			atomName2bond[name] = i;
		}
	}
	const RoboBond& getBondByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0, SimTK::Compound::AtomIndex cAIx1) const {
		const auto canonicalBond = canonicalizeBond(cAIx0, cAIx1);
		const auto name = subAtomList[canonicalBond.first].identity.uniqueAtomName + "-" + subAtomList[canonicalBond.second].identity.uniqueAtomName;
		const auto bondIt = atomName2bond.find(name);
		if (bondIt == atomName2bond.end()) {
			const auto& atom1Name = subAtomList[cAIx0].identity.uniqueAtomName;
			const auto& atom2Name = subAtomList[cAIx1].identity.uniqueAtomName;
			throw std::runtime_error("No bond found between " + atom1Name + " and " + atom2Name);
		}
		return subBondList[bondIt->second];
	}
	const Span<RoboBond> getBonds() const { return subBondList; }
	Span<RoboBond> updBonds() { return subBondList; }

	void setAngles(Span<RoboAngle> angles) {
		subAngleList = angles;

		for (std::size_t i = 0; i < subAngleList.size(); ++i) {
			const auto& angle = subAngleList[i];
			const auto canonicalAngle = canonicalizeAngle(angle.compoundAtomIndices[0], angle.compoundAtomIndices[1], angle.compoundAtomIndices[2]);
			const auto cAIx0 = canonicalAngle[0];
			const auto cAIx1 = canonicalAngle[1];
			const auto cAIx2 = canonicalAngle[2];

			const std::string name = subAtomList[cAIx0].identity.uniqueAtomName + "-" + subAtomList[cAIx1].identity.uniqueAtomName + "-" + subAtomList[cAIx2].identity.uniqueAtomName;
			atomName2angle[name] = i;
		}
	}
	const RoboAngle& getAngleByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0, SimTK::Compound::AtomIndex cAIx1, SimTK::Compound::AtomIndex cAIx2) const {
		const auto canonicalAngle = canonicalizeAngle(cAIx0, cAIx1, cAIx2);
		const auto name = subAtomList[canonicalAngle[0]].identity.uniqueAtomName + "-" + subAtomList[canonicalAngle[1]].identity.uniqueAtomName + "-" + subAtomList[canonicalAngle[2]].identity.uniqueAtomName;
		const auto angleIt = atomName2angle.find(name);
		if (angleIt == atomName2angle.end()) {
			const auto& atom1Name = subAtomList[cAIx0].identity.uniqueAtomName;
			const auto& atom2Name = subAtomList[cAIx1].identity.uniqueAtomName;
			const auto& atom3Name = subAtomList[cAIx2].identity.uniqueAtomName;
			throw std::runtime_error("No angle found between " + atom1Name + ", " + atom2Name + " and " + atom3Name);
		}
		return subAngleList[angleIt->second];
	}
	const Span<RoboAngle> getAngles() const { return subAngleList; }
	Span<RoboAngle> updAngles() { return subAngleList; }

	void setPeriodicTorsions(Span<RoboPeriodicTorsion> periodicTorsions) { subPeriodicTorsions = periodicTorsions; }
	const Span<RoboPeriodicTorsion> getPeriodicTorsions() const { return subPeriodicTorsions; }
	Span<RoboPeriodicTorsion> updPeriodicTorsions() { return subPeriodicTorsions; }

	void setImproperHarmonicTorsions(Span<RoboHarmonicImproperTorsion> improperHarmonicTorsions) { subImproperHarmonicTorsions = improperHarmonicTorsions; }
	const Span<RoboHarmonicImproperTorsion> getImproperHarmonicTorsions() const { return subImproperHarmonicTorsions; }
	Span<RoboHarmonicImproperTorsion> updImproperHarmonicTorsions() { return subImproperHarmonicTorsions; }

	/**	
	* @brief Get the name of this molecule
	* @return name of the molecule
	*/
	const std::string getName() const {return this->name;}

	/**	
	* @brief Get own CompoundIndex in CompoundSystem
	* 
	* This is equivalent to the index of the molecule in the CompoundSystem's vector of Compounds.
	* 
	* @return CompoundIndex
	*/
	inline SimTK::CompoundSystem::CompoundIndex getCompoundIndex() const {
		return compoundIndex;
	}

	/**
	* @brief Computes log(sin^2(pitch)) for the root atom's orientation in ground frame.
	* 
	* Extracts the pitch angle from the atom’s quaternion orientation.
	* The result is numerically stabilized near pitch ≈ 0 or ±π,
	* where sin(pitch) → 0 and log(sin²(pitch)) would diverge.
	* 
	* @note This quantity may represent an orientation regularization term.
	* @return log(sin²(pitch)), computed safely with analytic limits near singularities.
	*/
	SimTK::Real calcLogSineSqrGamma2(const SimTK::State &quatState) const;

	SimTK::Real calcLogDetMBATGamma2Contribution(const SimTK::State& quatState) const;

	/**
	 * @brief Get a reference to the atom object in the atom list of this Compound.
	 * 
	 * @param cAIx Compound Atom Index. This is in range [0, num_atoms-1] for this Compound. Not to confuse with the global atom index.
	 * @return Reference to the Atom object.
	 */
	const RoboAtom& getAtom(SimTK::Compound::AtomIndex cAIx) const;

	/**
	 * @brief Get a reference to the bond object in the bond list of this Compound.
	 * 
	 * This is not the Compound Atom Index but the global atom index.
	 * 
	 * @param aIx0 Global Atom Index of one atom in the bond.
	 * @param aIx1 Global Atom Index of the other atom in the bond.
	 * @return Reference to the BondLink object.
	 */
	const RoboBond& getBondByGlobalAtomIndex(int aIx0, int aIx1) const;

	/**	
	* @brief Get the bonded neighbor atom in the parent mobilized body.
	* @param aIx Compound Atom Index
	* @return Compound atom index of the root
	*/
	SimTK::Compound::AtomIndex getChemicalParentOfMobodRootAtom(SimTK::Compound::AtomIndex aIx, const SimTK::SimbodyMatterSubsystem& matter, const SimTK::DuMMForceFieldSubsystem& dumm) const;

	/**	
	* @brief Calculate all atom frames in top frame. It avoids calling 
	* calcDefaultAtomFrameInCompoundFrame multiple times. This has to be called
	* every time the coordinates change though.
	* @param : 
	* @return
	*/
	void calcAtomsTopTransforms();
	
	/**	
	* @brief 
	* @return
	*/
	void printTopTransforms();

	/**	
	* @brief Get atom Top level transform from the existing Topology map
	* @param cAIx: atom Compound AtomIndex
	* @return Atom's Top level transform
	*/
	const SimTK::Transform& getTopTransform(SimTK::Compound::AtomIndex cAIx) const;

	// Interface to access the maps

	// Retunr mbx by calling DuMM functions
	SimTK::MobilizedBodyIndex getAtomMobilizedBodyIndexThroughDumm(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm) const;

	// Get atom location on mobod through DuMM functions
	SimTK::Vec3 getAtomLocationInMobilizedBodyFrameThroughDumm(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm) const;

	SimTK::Vec3 calcAtomLocationInGroundFrameThroughSimbody(SimTK::Compound::AtomIndex aIx, const SimTK::DuMMForceFieldSubsystem& dumm, const SimTK::SimbodyMatterSubsystem& matter, const SimTK::State& someState) const;

	SimTK::Transform matchAtomTargetLocations(const SimTK::Compound::AtomTargetLocations& atomTargets);
	SimTK::Real getMatchError(const SimTK::Compound::AtomTargetLocations& atomTargets);

	void writeAtomListPdb(std::string dirname,
			              std::string prefix,
			              std::string sufix,
			              int maxNofDigits,
			              int index) const;

    /**
	 * @brief Create a mapping between the local compound atom indices and the global atom indices.
	 */
	void loadIndicesMaps(void);
	
	/**
	 * @brief Get the global atom index from the local compound atom index.
	 * @param cAIx Compound Atom Index
	 * @return Global Atom Index
	 */
	int getGlobalAtomIndex(SimTK::Compound::AtomIndex cAIx);
        
	/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index
	 * maps **/
	void printMaps();

	const std::vector<SimTK::Transform>& getAtomFrameCache() const { return atomFrameCache; }
	std::vector<SimTK::Transform>& updAtomFrameCache() { return atomFrameCache; }

private:
	Span<RoboAtom> subAtomList;
	Span<RoboBond> subBondList;
	Span<RoboAngle> subAngleList;
	Span<RoboPeriodicTorsion> subPeriodicTorsions;
	Span<RoboHarmonicImproperTorsion> subImproperHarmonicTorsions;

	std::unordered_map<std::string, std::size_t> atomName2bond, atomName2angle;

	// Map aIx to its Transform Default top transform
	std::vector<SimTK::Transform> aIx2TopTransform;

	// Map Atom number to aIx
	std::vector<int> compound2GlobalAtomIndex;
	std::vector<std::pair<int, SimTK::Compound::AtomIndex>> global2CompoundAtomIndex;


	std::vector<std::pair<CompoundAtomIndexPair, int>> aIxPair2Bonds;

	// Atom frames in Top frame
	std::vector<SimTK::Transform> atomFrameCache;

	std::string name;

	/** Every Compound has an index which is the position in the dvector
	 * of Compounds in CompoundSystem
	 */
	SimTK::CompoundSystem::CompoundIndex compoundIndex;

	std::size_t rootGlobalAtomIx = 0; // in global atom index
	SimTK::Compound::AtomIndex rootCompoundAtomIx; // in subAtomList index

	bool flipAllChirality = false;
};
