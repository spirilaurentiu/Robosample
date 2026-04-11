#include "Topology.hpp"

#include "Compound.h"

Topology::Topology(const SimTK::Compound::Name& name,
                   SimTK::CompoundSystem::CompoundIndex compoundIndex,
                   int rootGlobalAtomIx)
    : SimTK::Compound(name)
    , compoundIndex(compoundIndex)
    , rootGlobalAtomIx(rootGlobalAtomIx) {
    setCompoundName(name);
}

SimTK::Real Topology::calcLogSineSqrGamma2(const SimTK::State& quatState) const {
    // Get atom transform and convert to quaternion (w,x,y,z)
    const SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, rootCompoundAtomIx);
    const SimTK::Quaternion quat = X.R().convertRotationToQuaternion();

    const SimTK::Real w = quat[0];
    const SimTK::Real x = quat[1];
    const SimTK::Real y = quat[2];
    const SimTK::Real z = quat[3];

    // Compute sin(pitch) = 2(wy - zx)
    SimTK::Real sinPitch = 2.0 * ((w * y) - (z * x));

    // Clamp to account for floating-point drift outside [-1, 1]
    sinPitch = std::clamp(sinPitch, SimTK::Real(-1.0), SimTK::Real(1.0));

    // Compute pitch and evaluate the stable log(sin²)
    const SimTK::Real pitch = std::asin(sinPitch);
    return safeLogSineSqr(pitch);
}

SimTK::Real Topology::calcLogDetMBATGamma2Contribution(const SimTK::State& quatState) const {
    SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, rootCompoundAtomIx);
    SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();

    SimTK::Real w = quat[0];
    SimTK::Real x = quat[1];
    SimTK::Real y = quat[2];
    SimTK::Real z = quat[3];
    SimTK::Real sinPitch = 2 * (w * y - z * x);

    std::cout << std::setprecision(20) << std::fixed;
    std::cout << "sinpitch " << sinPitch << std::endl;
    std::cout << "sinpitchsq" << sinPitch * sinPitch << std::endl;

    SimTK::Real pitch = std::asin(sinPitch);
    std::cout << "pitch " << pitch << std::endl;
    // if(pitch < 0){
    //	pitch = pitch + (2*SimTK_PI);
    //	std::cout << "sin converted pitch " << std::sin(pitch) << std::endl;
    // }

    if (sinPitch < SimTK::Eps) { // consider using SimTK::Eps
        return -SimTK::Infinity;
    }
    SimTK::Real result = std::log(sinPitch * sinPitch);
    return result;
}

/** Get a pointer to an atom object in the atom list inquiring
by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
// TODO: Optimize use CompoundAtomIx2GmolAtomIx instead
const RoboAtom& Topology::getAtom(SimTK::Compound::AtomIndex cAIx) const {
    for (const auto& atom : subAtomList) {
        if (atom.identity.compoundAtomIndex == cAIx) {
            return atom;
        }
    }

    // This should never trigger, but just in case
    SimTK_ASSERT_ALWAYS(false, "Topology::getAtom(): Atom with specified Compound::AtomIndex not found.");
}

/** **/
const RoboBond& Topology::getBondByGlobalAtomIndex(int aIx0, int aIx1) const {
    // for (const auto& b : subBondList) {
    // 	if (b.bond.globalIndices[0] == aIx0 && b.bond.globalIndices[1] == aIx1) {
    // 		return b;
    // 	} else if (b.bond.globalIndices[0] == aIx1 && b.bond.globalIndices[1] == aIx0) {
    // 		return b;
    // 	}
    // }

    // This should never trigger, but just in case
    SimTK_ASSERT_ALWAYS(false,
                        "Topology::getBondByGlobalAtomIndex(): No bond with specified atom indices found.");
}

void Topology::loadIndicesMaps() {
    // // print spans
    // for (const auto& a : subAtomList) {
    // 	std::cout << "Topology::loadIndicesMaps(): Atom " << a.identity.uniqueAtomName << std::endl;
    // }
    // for (const auto& b : subBondList) {
    // 	std::cout << "Topology::loadIndicesMaps(): Bond between global atom indices "
    // 		<< b.bond.globalIndices[0] << " and " << b.bond.globalIndices[1] << std::endl;
    // }
    // for (const auto& angle : subAngleList) {
    // 	std::cout << "Topology::loadIndicesMaps(): Angle between global atom indices "
    // 		<< angle.getGlobalIndex1() << ", " << angle.getGlobalIndex2() << ", " << angle.getGlobalIndex3()
    // << std::endl;
    // }
    // for (const auto& torsion : subTorsionList) {
    // 	std::cout << "Topology::loadIndicesMaps(): Torsion between global atom indices "
    // 		<< torsion.getGlobalIndex1() << ", " << torsion.getGlobalIndex2() << ", "
    // 		<< torsion.getGlobalIndex3() << ", " << torsion.getGlobalIndex4() << std::endl;
    // }

    aIx2TopTransform = std::vector<SimTK::Transform>(subAtomList.size(), SimTK::Transform());

    // Find the root atom index in the subAtomList
    for (const auto& a : subAtomList) {
        if (a.connectivity.root) {
            rootCompoundAtomIx = a.identity.compoundAtomIndex;
            break;
        }
    }

    // compound2GlobalAtomIndex = std::vector<int>(subAtomList.size(), -1);
    // global2CompoundAtomIndex = std::vector<std::pair<int, SimTK::Compound::AtomIndex>>(subAtomList.size(),
    // std::make_pair(1, SimTK::Compound::AtomIndex(1))); // needs to be initialized to a valid (positive)
    // Compound Atom Index for (const auto& a : subAtomList) { 	SimTK::Compound::AtomIndex aIx =
    // a.identity.compoundAtomIndex; 	int globalAtomIndex = a.getGlobalIndex();

    // 	compound2GlobalAtomIndex[aIx] = globalAtomIndex;
    // 	global2CompoundAtomIndex[globalAtomIndex] = std::make_pair(globalAtomIndex, aIx);
    // }

    // // TODO Doesn't work for two molecules and i don't understand why this code exists in the first place
    // // Traverse all bonds and save their indices
    // for (std::size_t i = 0; i < subAtomList.size(); ++i) {
    // 	const auto& b = subBondList[i];

    // 	// Highly inefficient way to get the Compound::AtomIndex from the global atom index
    // 	SimTK::Compound::AtomIndex cAIx0;
    // 	int aIx0 = b.bond.globalIndices[1];
    // 	for (const auto& a : subAtomList) {
    // 		if (a.getGlobalIndex() == aIx0) {
    // 			cAIx0 = a.identity.compoundAtomIndex;
    // 			break;
    // 		}
    // 	}
    // 	SimTK_ASSERT_ALWAYS(cAIx0.isValid() && cAIx0.isValidExtended(), "Topology::loadIndicesMaps(): Invalid
    // Compound::AtomIndex found.");

    // 	// Another inefficient way to get the Compound::AtomIndex from the global atom index
    // 	SimTK::Compound::AtomIndex cAIx1;
    // 	int aIx1 = b.bond.globalIndices[0];
    // 	for (const auto& a : subAtomList) {
    // 		if (a.getGlobalIndex() == aIx1) {
    // 			cAIx1 = a.identity.compoundAtomIndex;
    // 			break;
    // 		}
    // 	}
    // 	SimTK_ASSERT_ALWAYS(cAIx1.isValid() && cAIx1.isValidExtended(), "Topology::loadIndicesMaps(): Invalid
    // Compound::AtomIndex found.");

    // 	// We sort pairs to make searching easier
    // 	CompoundAtomIndexPair pair = canonical(cAIx0, cAIx1);
    // 	aIxPair2Bonds.push_back(std::make_pair(pair, i));
    // }

    // // Sort the vector of pairs
    // std::sort(aIxPair2Bonds.begin(), aIxPair2Bonds.end());

    atomFrameCache = std::vector<SimTK::Transform>(subAtomList.size(), SimTK::Transform());
}

int Topology::getGlobalAtomIndex(SimTK::Compound::AtomIndex cAIx) {
    SimTK_ASSERT_ALWAYS(false, "Topology::getGlobalAtomIndex(): Not implemented yet.");
    // return compound2GlobalAtomIndex[cAIx];
}

/*!
 * <!-- Calculate all atom frames in top frame. It avoids calling
 * calcDefaultAtomFrameInCompoundFrame multiple times. This has
 * to be called every time the coordinates change though. -->
 */
void Topology::calcAtomsTopTransforms() {
    for (const auto& a : subAtomList) {
        SimTK::Compound::AtomIndex aIx = a.identity.compoundAtomIndex;
        aIx2TopTransform[aIx] = calcDefaultAtomFrameInCompoundFrame(aIx);
    }
}

/*!
 * <!--  -->
 */
void Topology::printTopTransforms() {
    std::cout << "Topology TopTransforms " << std::endl;
    for (unsigned int i = 0; i < getNumAtoms(); ++i) {
        SimTK::Compound::AtomIndex aIx = (subAtomList[i]).identity.compoundAtomIndex;
        std::cout << aIx << " " << aIx2TopTransform[aIx] << std::endl;
    }
}

/*!
 * <!-- Get atom Top level transform from the existing Topology map -->
 */
const SimTK::Transform& Topology::getTopTransform(SimTK::Compound::AtomIndex aIx) const {
    return aIx2TopTransform[aIx];
}

// Return mbx by calling DuMM functions
SimTK::MobilizedBodyIndex
Topology::getAtomMobilizedBodyIndexThroughDumm(SimTK::Compound::AtomIndex aIx,
                                               const SimTK::DuMMForceFieldSubsystem& dumm) const {
    SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
    return dumm.getAtomBody(dAIx);
}

// Get atom location on mobod through DuMM functions
SimTK::Vec3
Topology::getAtomLocationInMobilizedBodyFrameThroughDumm(SimTK::Compound::AtomIndex aIx,
                                                         const SimTK::DuMMForceFieldSubsystem& dumm) const {
    SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
    return dumm.getAtomStationOnBody(dAIx);
}

SimTK::Vec3 Topology::calcAtomLocationInGroundFrameThroughSimbody(SimTK::Compound::AtomIndex aIx,
                                                                  const SimTK::DuMMForceFieldSubsystem& dumm,
                                                                  const SimTK::SimbodyMatterSubsystem& matter,
                                                                  const SimTK::State& someState) const {
    const SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm);
    const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx); // vector

    const SimTK::Transform& X_GB = mobod.getBodyTransform(someState);
    const SimTK::Rotation& R_GB = X_GB.R();
    const SimTK::Vec3& p_GB = X_GB.p();

    // ROBOFACTOR this is part of a vector and should accessed in sequence
    SimTK::Vec3 station = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);

    const SimTK::Vec3 p_BS_G = R_GB * station;
    return p_GB + p_BS_G;
}

SimTK::Transform Topology::matchAtomTargetLocations(const SimTK::Compound::AtomTargetLocations& atomTargets) {
    matchDefaultBondLengths(atomTargets);
    matchDefaultAtomChirality(atomTargets, 0.01, flipAllChirality);
    matchDefaultBondAngles(atomTargets);
    matchDefaultDirections(atomTargets);
    matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
    matchDefaultTopLevelTransform(atomTargets);

    // Get the Ground to Top Transform
    const SimTK::Transform G_X_T = getTopLevelTransform();

    // Recalculate atom frames in top compound frame
    calcAtomsTopTransforms();

    // // Ensure residual is low enough only on debug builds
    // assert(getTransformAndResidual(atomTargets).residual < 1e-5 && "After setAtoms_Compound_Match, residual
    // too high");

    return G_X_T;
}

SimTK::Real Topology::getMatchError(const SimTK::Compound::AtomTargetLocations& atomTargets) {
    // std::vector<SimTK::Transform> atomSourceFrames(getNumAtoms());
    // invalidateAtomFrameCache(atomSourceFrames, getNumAtoms());
    // calcDefaultAtomFramesInCompoundFrame(atomSourceFrames);

    // SimTK::Real cumulativeError = 0.0;
    // for (const auto& target : atomTargets)
    // {
    // 	const Compound::AtomIndex atomIndex = target.first;
    // 	const SimTK::Vec3& targetVec = target.second;
    // 	const SimTK::Vec3 source = getTopLevelTransform() * atomSourceFrames[atomIndex].T();

    // 	cumulativeError += (source - targetVec).norm();
    // }

    // return cumulativeError;

    return 0;
}

/** Print maps **/
void Topology::printMaps() {
    // std::cout << "Topology map CompoundAtomIx2GmolAtomIx:" << std::endl;
    // std::map< SimTK::Compound::AtomIndex, int >::const_iterator aIx2gmolaIxIt;
    // for(aIx2gmolaIxIt = compound2GlobalAtomIndex.begin();
    //    aIx2gmolaIxIt != compound2GlobalAtomIndex.end(); ++aIx2gmolaIxIt)
    // {
    // 	std::cout << "atomIndex " << aIx2gmolaIxIt->first
    // 		<< " gmolaIx " << aIx2gmolaIxIt->second
    // 		<< std::endl << std::flush;
    // }
}

/** Write a pdb with bAtomList coordinates and inNames **/
void Topology::writeAtomListPdb(std::string dirname,
                                std::string prefix,
                                std::string sufix,
                                int maxNofDigits,
                                int index) const {
    // // Using floor here is no buneo because the index can be zero
    // std::string zeros("");
    // int nofDigits = static_cast<int>(std::to_string(index).size());
    // if(maxNofDigits > nofDigits){
    // 	zeros = std::string(maxNofDigits - nofDigits, '0');
    // }

    // std::stringstream sstream;
    // sstream << dirname << "/"
    // 	<< prefix << zeros << std::to_string(index) << sufix;
    // std::string ofilename = sstream.str();

    // FILE *oF = fopen (ofilename.c_str(),"w");
    // if (oF) {
    // 	// Pdb lines
    // 	for(int i = 0; i < getNumAtoms(); i++){
    // 		fprintf(oF, "%-6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  %4.2f%6.2f          %2s\n"
    // 			, "ATOM"                 // record
    // 			, i                      // index
    // 			, subAtomList[i].getInName().c_str()  // name
    // 			, "UNK"                  // residue name
    // 			, 'A'                    // chain
    // 			, 1                      // residue index
    // 			, 10.0*subAtomList[i].getX()    // x in A
    // 			, 10.0*subAtomList[i].getY()    // y in A
    // 			, 10.0*subAtomList[i].getZ()    // z in A
    // 			, 1.0                    // occupancy
    // 			, 0.0                    // beta factor
    // 			, "  ");                 // element
    // 	}

    // 	fclose(oF);
    // 	//std::cout << "\tTopology written to '" << ofilename << "'\n";
    // } else {
    // 	std::cout << "FAILED TO OPEN '" << ofilename << "' TO WRIE!\n";
    // }
}

/** Get the neighbor atom bonded to aIx atom in the parent mobilized body.
TODO: No chemical parent for satellite subAtomList or first atom. **/
SimTK::Compound::AtomIndex
Topology::getChemicalParentOfMobodRootAtom(SimTK::Compound::AtomIndex aIx,
                                           const SimTK::SimbodyMatterSubsystem& matter,
                                           const SimTK::DuMMForceFieldSubsystem& dumm) const {
    // Check if this atom is at the origin of the mobilized body
    SimTK::Vec3 v = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);
    SimTK_ASSERT_ALWAYS(v.norm() < 1e-12, "Atom is not at the origin of the mobilized body");

    // Get body, parentBody, parentAtom
    SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm);
    const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);
    const SimTK::MobilizedBody& parentMobod = mobod.getParentMobilizedBody();

    SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
    SimTK_ASSERT_ALWAYS(parentMbx != 0, "Parent mobilized body is Ground"); // what is tho

    SimTK::Compound::AtomIndex chemParentAIx; // Invalid by default
    const RoboAtom& origin = getAtom(aIx);

    // Traverse all neighbors and check if both bonded aatoms are in the parent mobod
    for (auto neighborGlobalIx : origin.connectivity.neighborsGlobalIndices) {
        // Get this bond
        for (const auto& b : subBondList) {
            if ((b.globalIndices[1] == origin.identity.globalIndex && b.globalIndices[0] == neighborGlobalIx)
                || (b.globalIndices[0] == origin.identity.globalIndex
                    && b.globalIndices[1] == neighborGlobalIx)) {
                // We found the bond, now get the neighbor atom
                for (const auto& neighbor : subAtomList) {
                    if (neighbor.identity.globalIndex == neighborGlobalIx) {
                        // We found the neighbor atom
                        Compound::AtomIndex candidateChemParentAIx = neighbor.identity.compoundAtomIndex;

                        // Check if neighbor atom's mobod is a parent mobod
                        if (getAtomMobilizedBodyIndexThroughDumm(candidateChemParentAIx, dumm) == parentMbx) {
                            if (!b.ringClosing) { // No ring subAtomList are allowed
                                chemParentAIx = candidateChemParentAIx;
                                return chemParentAIx;
                            }
                        }
                    }
                }
            }
        }
    }

    return chemParentAIx;
}
