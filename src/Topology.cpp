#include "Topology.hpp"

#include "Compound.h"

inline auto canonical(SimTK::Compound::AtomIndex a, SimTK::Compound::AtomIndex b) -> CompoundAtomIndexPair {
    return (a < b) ? std::make_pair(a, b) : std::make_pair(b, a);
}

Topology::Topology(const SimTK::Compound::Name& name,
                   SimTK::CompoundSystem::CompoundIndex compoundIndex,
                   int rootGlobalAtomIx,
                   SimTK::RootMobility rootMobility)
    : SimTK::Compound(name)
    , compoundIndex(compoundIndex)
    , rootGlobalAtomIx(rootGlobalAtomIx)
    , rootMobility(rootMobility) {
    setCompoundName(name);
}

void Topology::setAtoms(Span<RoboAtom> atoms) {
    subAtomList = atoms;
}

void Topology::setBonds(Span<RoboBond> bonds) {
    subBondList = bonds;

    for (std::size_t i = 0; i < subBondList.size(); ++i) {
        const auto& bond = subBondList[i];
        const auto canonicalBond = canonicalizeBond(bond.compoundAtomIndices[0], bond.compoundAtomIndices[1]);
        const auto cAIx0 = canonicalBond.first;
        const auto cAIx1 = canonicalBond.second;

        const std::string name =
            subAtomList[cAIx0].identity.uniqueAtomName + "-" + subAtomList[cAIx1].identity.uniqueAtomName;
        atomName2bond[name] = i;
    }
}

auto Topology::getBondByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0,
                                          SimTK::Compound::AtomIndex cAIx1) const -> const RoboBond& {
    const auto canonicalBond = canonicalizeBond(cAIx0, cAIx1);
    const auto name = subAtomList[canonicalBond.first].identity.uniqueAtomName + "-"
                      + subAtomList[canonicalBond.second].identity.uniqueAtomName;
    const auto bondIt = atomName2bond.find(name);
    if (bondIt == atomName2bond.end()) {
        const auto& atom1Name = subAtomList[cAIx0].identity.uniqueAtomName;
        const auto& atom2Name = subAtomList[cAIx1].identity.uniqueAtomName;
        throw std::runtime_error("No bond found between " + atom1Name + " and " + atom2Name);
    }
    return subBondList[bondIt->second];
}

void Topology::setAngles(Span<RoboAngle> angles) {
    subAngleList = angles;

    for (std::size_t i = 0; i < subAngleList.size(); ++i) {
        const auto& angle = subAngleList[i];
        const auto canonicalAngle = canonicalizeAngle(angle.compoundAtomIndices[0],
                                                      angle.compoundAtomIndices[1],
                                                      angle.compoundAtomIndices[2]);
        const auto cAIx0 = canonicalAngle[0];
        const auto cAIx1 = canonicalAngle[1];
        const auto cAIx2 = canonicalAngle[2];

        const std::string name = subAtomList[cAIx0].identity.uniqueAtomName + "-"
                                 + subAtomList[cAIx1].identity.uniqueAtomName + "-"
                                 + subAtomList[cAIx2].identity.uniqueAtomName;
        atomName2angle[name] = i;
    }
}

auto Topology::getAngleByCompoundAtomIndex(SimTK::Compound::AtomIndex cAIx0,
                                           SimTK::Compound::AtomIndex cAIx1,
                                           SimTK::Compound::AtomIndex cAIx2) const -> const RoboAngle& {
    const auto canonicalAngle = canonicalizeAngle(cAIx0, cAIx1, cAIx2);
    const auto name = subAtomList[canonicalAngle[0]].identity.uniqueAtomName + "-"
                      + subAtomList[canonicalAngle[1]].identity.uniqueAtomName + "-"
                      + subAtomList[canonicalAngle[2]].identity.uniqueAtomName;
    const auto angleIt = atomName2angle.find(name);
    if (angleIt == atomName2angle.end()) {
        const auto& atom1Name = subAtomList[cAIx0].identity.uniqueAtomName;
        const auto& atom2Name = subAtomList[cAIx1].identity.uniqueAtomName;
        const auto& atom3Name = subAtomList[cAIx2].identity.uniqueAtomName;
        throw std::runtime_error("No angle found between " + atom1Name + ", " + atom2Name + " and "
                                 + atom3Name);
    }
    return subAngleList[angleIt->second];
}

void Topology::setPeriodicTorsions(Span<RoboPeriodicTorsion> periodicTorsions) {
    subPeriodicTorsions = periodicTorsions;
}

void Topology::setImproperHarmonicTorsions(Span<RoboHarmonicImproperTorsion> improperHarmonicTorsions) {
    subImproperHarmonicTorsions = improperHarmonicTorsions;
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

void Topology::loadIndicesMaps(const SimTK::Compound::AtomTargetLocations& atomTargets) {
    // Find the root atom index in the subAtomList
    for (const auto& a : subAtomList) {
        if (a.connectivity.root) {
            rootCompoundAtomIx = a.identity.compoundAtomIndex;
            break;
        }
    }

    atomFrameCache = std::vector<SimTK::Transform>(subAtomList.size(), SimTK::Transform());
    buildCache(atomTargets);
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
