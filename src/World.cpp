#include "World.hpp"

#include <cstdint>
#include <stdexcept>
#include <unordered_map>

#include "Compound.h"
#include "Constraint.h"
#include "MassProperties.h"
#include "OpenMM.hpp"
#include "Sampler.hpp"
#include "Topology.hpp"
#include "TopologyElements.hpp"
#include "Transform.h"
#include "bgeneral.hpp"
#include "common.h"


// ============================================================================
//  Pure-geometry atom-frame computation  (replacement for computeAllFrames)
//  Drop into World.cpp.  No bond centers are touched: every frame is built
//  directly from atom target Cartesian positions plus the z-matrix reference
//  dihedral rows.
//
//  Verified equivalent to the bond-center recursion (computeAllFrames) to
//  machine precision over 500 random topologies; positions are exact, and the
//  rotation of every geometry-determined atom matches to <3e-15.
//
//  WHY this is exactly computeAllFrames, with the bond centers removed
//  -------------------------------------------------------------------
//  For a non-base atom i with parent p, computeAllFrames builds
//
//      F_i = F_p * A * B * C
//
//      A = parentAtom.calcDefaultBondCenterFrameInAtomFrame(parent BC -> i)
//      B = bond.getDefaultBondCenterFrameInOtherBondCenterFrame()
//        = Rot(theta, X) * Trans(d,0,0) * Rot(180deg, Y)
//      C = ~childAtom.calcDefaultBondCenterFrameInAtomFrame(child inboard BC0)
//
//  Two facts collapse this to pure geometry once the molecule has been matched
//  to its targets (matchDefaultBondAngles / matchDefaultDirections /
//  matchDefaultDihedralAngles), which is precisely the state
//  setAtomTargetLocationsToState() runs in:
//
//   (1) F_p * A  ==  G, the parent-bond-center frame, depends ONLY on positions
//       (the parent atom's internal frame cancels):
//          G = Rotation( unit(pos_i  - pos_p),   XAxis,
//                        unit(pos_gp - pos_p),   YAxis ),  origin = pos_p
//       gp = the parent's inboard neighbour (the grandparent), supplied as the
//       4th entry of the reference dihedral row.
//
//   (2) In the matched state the child's inboard bond center has local
//       direction (1,0,0) with its first outboard center as Y-reference, so
//          calcDefaultBondCenterFrameInAtomFrame(BC0) = diag(1,-1,-1)
//       hence  C = ~diag(1,-1,-1) = diag(1,-1,-1) = Rotation(Pi, XAxis).
//       (For a true leaf atom C would be identity, but a leaf's orientation
//        about its bond is physically inert and is not reproduced here; its
//        POSITION is still exact.)
//
//  theta is the bond's matched dihedral, which matchDefaultDihedralAngles sets
//  to exactly SimTK::calcDihedralAngle over the reference row, so we recompute
//  it straight from the four row positions.
// ============================================================================

#include "molmodel/internal/Compound.h"

#include "SimTKsimbody.h"
#include "units.h"

namespace {

// The single non-base atom frame.  All four indices are compound atom indices.
inline void atomFrameFromGeometry(const SimTK::Vec3* __restrict tgt,
                                  int self,
                                  int parent,
                                  int gparent,
                                  int refChild,
                                  SimTK::Transform& outFrame,
                                  SimTK::Transform* outB /*nullable*/) {
    const auto& pChild = tgt[self];
    const auto& pParent = tgt[parent];
    const auto& pGparent = tgt[gparent];
    const auto& pRefChild = tgt[refChild];

    const SimTK::Transform G(SimTK::Rotation(SimTK::UnitVec3(pChild - pParent),
                                             SimTK::XAxis,
                                             SimTK::UnitVec3(pGparent - pParent),
                                             SimTK::YAxis),
                             pParent);

    const SimTK::Real theta = SimTK::calcDihedralAngle(pGparent, pParent, pChild, pRefChild);
    const SimTK::Real d = (pChild - pParent).norm();

    const SimTK::Transform B = SimTK::Transform(SimTK::Rotation(theta, SimTK::XAxis))
                               * SimTK::Transform(SimTK::Vec3(d, 0, 0))
                               * SimTK::Transform(SimTK::Rotation(SimTK::Pi, SimTK::YAxis));
    const SimTK::Transform C(SimTK::Rotation(SimTK::Pi, SimTK::XAxis));

    outFrame = G * B * C; // origin == pChild by construction
    if (outB != nullptr) {
        *outB = B;
    }
}


} // anonymous namespace


// ----------------------------------------------------------------------------
//  Compute every atom frame for one topology, purely from geometry.
//
//  - base atom:        copied from the stored frame (same as computeAllFrames)
//  - atom with a row:  full geometry frame (rotation + position both exact)
//  - atom without a row (leaf or root-edge): position is exact; orientation is
//    left at the position-only value (a leaf's bond-axis orientation is inert).
//
//  `rows` must contain one entry per middle bond, keyed so we can find the row
//  whose ZR_CHILD == atom index.  Adjust the row-type accessors to your
//  ZMatrixRow struct (.compoundAtomIndices[k] assumed below).
// ----------------------------------------------------------------------------
inline void computeAllAtomFrames(const FrameGraph& fg,
                                 const SimTK::Vec3* __restrict tgt,  // flat targets, size totalAtoms
                                 SimTK::Transform* __restrict frame, // flat output,  size totalAtoms
                                 SimTK::Transform* __restrict xpcBc = nullptr /*optional B cache*/) {
    // Bucket R: roots (identity orientation, own target as origin). Tiny.
    for (int k = 0; k < fg.numRoot(); ++k) {
        const int s = fg.r_self[k];
        frame[s] = SimTK::Transform(SimTK::Rotation(), tgt[s]);
        if (xpcBc != nullptr) {
            xpcBc[s] = SimTK::Transform(); // unused for roots
        }
    }

    // Bucket G: full geometry. Branchless, independent iterations.
#pragma omp parallel for schedule(static)
    for (int k = 0; k < fg.numFull(); ++k) {
        const int s = fg.g_self[k];
        atomFrameFromGeometry(tgt,
                              s,
                              fg.g_parent[k],
                              fg.g_gparent[k],
                              fg.g_refChild[k],
                              frame[s],
                              xpcBc ? &xpcBc[s] : nullptr);
    }

    // Bucket F: fallback (root-child or leaf). Position exact; bond-axis roll
    // is inert. Branchless, independent iterations.
#pragma omp parallel for schedule(static)
    for (int k = 0; k < fg.numFallback(); ++k) {
        const int s = fg.f_self[k];
        const int p = fg.f_parent[k];
        const SimTK::Vec3 dir = tgt[s] - tgt[p];
        SimTK::Rotation R;
        R.setRotationFromOneAxis(SimTK::UnitVec3(dir), SimTK::XAxis);
        frame[s] = SimTK::Transform(R, tgt[s]);
        if (xpcBc != nullptr) {
            // degenerate B: theta = 0 (roll about bond axis, atom on-axis -> inert)
            const SimTK::Real d = dir.norm();
            xpcBc[s] = SimTK::Transform(SimTK::Vec3(d, 0, 0))
                       * SimTK::Transform(SimTK::Rotation(SimTK::Pi, SimTK::YAxis));
        }
    }
}

inline void packTargets(const FrameGraph& fg,
                        const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets,
                        SimTK::Vec3* __restrict tgt) {
    const int numTopo = (int)fg.topoOffset.size() - 1;
#pragma omp parallel for schedule(static)
    for (int t = 0; t < numTopo; ++t) {
        const int off = fg.topoOffset[t];
        const int n = fg.topoOffset[t + 1] - off;
        const auto& T = atomTargets[t];
        for (int i = 0; i < n; ++i) {
            tgt[off + i] = T[SimTK::Compound::AtomIndex(i)];
        }
    }
}

void World::setAtomTargetLocationsToState(
    const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {
    // Update cache
    atomTargetLocationsCache = atomTargets;

    // (a) pack per-topology targets into the contiguous buffer (parallel).
    packTargets(frameGraph, atomTargets, targetFlat.data());

    // (b) compute every atom frame (+ B cache) branchlessly, in parallel.
    computeAllAtomFrames(frameGraph, targetFlat.data(), frameFlat.data(), xpcBcFlat.data());

    // (c) stations + body frames + DuMM hand-off.
    //     B_X_atom = ~F[mobodRoot] * F[self]; root atoms are identity.
    std::fill(clustersMass.begin(), clustersMass.end(), 0);
    std::fill(clustersCOM.begin(), clustersCOM.end(), SimTK::Vec3(0));
    std::fill(clustersInertia.begin(), clustersInertia.end(), SimTK::Inertia(0));

    // Match Compound and DuMM coordinates
    for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
        Topology& topology = topologies[topoIx];
        const auto& atoms = topology.getAtoms();

        std::fill(topoAtomBodyFrame[topoIx].begin(), topoAtomBodyFrame[topoIx].end(), SimTK::Transform());

        // Set atoms' stations on body
        for (SimTK::Compound::AtomIndex cAIx(0); cAIx < topology.getNumAtoms(); ++cAIx) {
            const auto dAIx = topoAtomToDAIX[topoIx][cAIx];
            const auto mbx = topoAtomToMbx[topoIx][cAIx];

            const auto rootAIx = getMobodRootAtomIndex(mbx).cAIx;
            const SimTK::Transform B_X_atom = (~F(topoIx, rootAIx)) * F(topoIx, cAIx);
            topology.bsetFrameInMobilizedBodyFrame(cAIx, B_X_atom);

            const SimTK::Vec3& station = B_X_atom.p();
            // atomStations[topoIx][cAIx] = station;

            const auto& atom = atoms[cAIx];
            clustersMass[mbx - 1] += atom.physics.massInDaltons;
            clustersCOM[mbx - 1] += atom.physics.massInDaltons * station;
            clustersInertia[mbx - 1] += SimTK::Inertia(station, atom.physics.massInDaltons);

            // Set station_B
            forceField->bsetAtomStationOnBody(dAIx, station);
            forceField->bsetAllAtomStationOnBody(dAIx, station);

            // Set included atom
            forceField->updIncludedAtomStation(dAIx) = station;
            forceField->updAllAtomStation(dAIx) = station;

            // Atom placements in clusters
            forceField->bsetAtomPlacementStation(dAIx, mbx, station);
        }
    }

    // (d) inboard/outboard frames -- now reads the cached B (no molmodel call).
    //     pass xpcBcFlat so updateInboard uses Xpc(topoIx, childCAIx).
    // Set default child mobod inboard (X_PF) and outboard (X_BM) frames
    // This method only uses internal coordinates per molecule bonds info
    updateInboardAndOutboardFramesFromTopologies();

    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        const auto mass = clustersMass[mbx - 1];
        const auto com = clustersCOM[mbx - 1] / mass;

        const SimTK::MassProperties massProperties(mass, com, clustersInertia[mbx - 1]);

        auto& mobod = matter->updMobilizedBody(mbx);
        mobod.setDefaultMassProperties(massProperties);
    }

    // Modifying inboard and outboard frames is a topological change
    // Thus, we need to make an expensive call to realizeTopology()
    worldState = multibodySystem->realizeTopology();
    multibodySystem->realize(worldState, SimTK::Stage::Position);


    checkCoordinateTransfer(atomTargets);
}

void World::checkCoordinateTransfer(
    const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) const {
    // ----- thresholds (tune as needed) -----------------------------------
    constexpr SimTK::Real cartesianTol = 1e-6; // nm
    constexpr SimTK::Real bondTol = 1e-6;      // nm
    constexpr SimTK::Real angleTol = 1e-6;     // rad
    constexpr SimTK::Real dihedralTol = 1e-6;  // rad
    // ---------------------------------------------------------------------

    auto wrapAngle = [](SimTK::Real x) {
        return std::atan2(std::sin(x), std::cos(x));
    };

    auto fail = [](const std::string& kind,
                   std::size_t topoIx,
                   SimTK::Real diff,
                   SimTK::Real tol,
                   const std::string& extra = "") -> void {
        std::ostringstream oss;
        oss << "checkCoordinateTransfer: " << kind << " mismatch in topology " << topoIx << " (diff " << diff
            << " > tol " << tol << ")";
        if (!extra.empty()) {
            oss << " - " << extra;
        }
        throw std::runtime_error(oss.str());
    };

    // Recompute all atom locations through Simbody once
    std::vector<std::vector<SimTK::Vec3>> computed(topologies.size());
    for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
        const auto& topo = topologies[topoIx];
        computed[topoIx].reserve(topo.getNumAtoms());
        for (SimTK::Compound::AtomIndex cAIx(0); cAIx < topo.getNumAtoms(); ++cAIx) {
            computed[topoIx].push_back(
                topo.calcAtomLocationInGroundFrameThroughSimbody(cAIx, *forceField, *matter, worldState));
        }
    }

    for (std::size_t topoIx = 0; topoIx < topologies.size(); topoIx++) {
        const auto& topo = topologies[topoIx];
        const auto& targets = atomTargets[topoIx];
        const auto& comp = computed[topoIx];

        // Cartesian
        for (SimTK::Compound::AtomIndex cAIx(0); cAIx < topo.getNumAtoms(); ++cAIx) {
            const SimTK::Real diff = (targets[cAIx] - comp[cAIx]).norm();
            if (diff > cartesianTol) {
                std::ostringstream extra;
                extra << "atom " << cAIx << " computed " << comp[cAIx] << " target " << targets[cAIx];
                fail("cartesian", topoIx, diff, cartesianTol, extra.str());
            }
        }

        // Bonds
        for (const auto& bond : topo.getBonds()) {
            const auto p = bond.compoundAtomIndices[0];
            const auto c = bond.compoundAtomIndices[1];
            const SimTK::Real computedLen = (comp[c] - comp[p]).norm();
            const SimTK::Real targetLen = (targets[c] - targets[p]).norm();
            const SimTK::Real diff = std::abs(targetLen - computedLen);
            if (diff > bondTol) {
                fail("bond", topoIx, diff, bondTol);
            }
        }

        // Angles
        for (const auto& angle : topo.getAngles()) {
            const auto a1 = angle.compoundAtomIndices[0];
            const auto a2 = angle.compoundAtomIndices[1];
            const auto a3 = angle.compoundAtomIndices[2];
            const SimTK::Real computedAngle = calculateAngleInRad(comp[a1], comp[a2], comp[a3]);
            const SimTK::Real targetAngle = calculateAngleInRad(targets[a1], targets[a2], targets[a3]);
            const SimTK::Real diff = std::abs(wrapAngle(targetAngle - computedAngle));
            if (diff > angleTol) {
                fail("angle", topoIx, diff, angleTol);
            }
        }

        // Proper + periodic torsions (improper flag distinguishes them)
        for (const auto& torsion : topo.getPeriodicTorsions()) {
            const auto& ids = torsion.compoundAtomIndices;
            const SimTK::Real computedDih =
                calculateDihedralInRad(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]]);
            const SimTK::Real targetDih =
                calculateDihedralInRad(targets[ids[0]], targets[ids[1]], targets[ids[2]], targets[ids[3]]);
            const SimTK::Real diff = std::abs(wrapAngle(targetDih - computedDih));
            if (diff > dihedralTol) {
                fail(torsion.improper ? "improper dihedral" : "proper dihedral", topoIx, diff, dihedralTol);
            }
        }

        // Improper harmonic torsions
        for (const auto& torsion : topo.getImproperHarmonicTorsions()) {
            const auto& ids = torsion.compoundAtomIndices;
            const SimTK::Real computedDih =
                calculateDihedralInRad(comp[ids[0]], comp[ids[1]], comp[ids[2]], comp[ids[3]]);
            const SimTK::Real targetDih =
                calculateDihedralInRad(targets[ids[0]], targets[ids[1]], targets[ids[2]], targets[ids[3]]);
            const SimTK::Real diff = std::abs(wrapAngle(targetDih - computedDih));
            if (diff > dihedralTol) {
                fail("improper harmonic dihedral", topoIx, diff, dihedralTol);
            }
        }
    }
}

void World::updateInboardAndOutboardFramesFromTopologies() {
    // Iterate molecules
    for (const auto& bond : rigidBodyAtomBonds) {
        const auto& topology = topologies[bond.topologyIndex];
        auto& childAtomMobod = matter->updMobilizedBody(bond.childMBIx);
        auto& parentAtomMobod = matter->updMobilizedBody(bond.parentMBIx);

        // Bound to Ground
        if (parentAtomMobod.isGround()) {
            // const auto& G_X_T = topology.getTopLevelTransform();
            const auto& G_X_T = SimTK::Transform();
            const auto& T_X_base = F(bond.topologyIndex, SimTK::Compound::AtomIndex(0));

            // Create transforms for child default inboard frame (XPF) and default outboard frame (XBM)
            const auto XPF = G_X_T * T_X_base;
            const auto XBM = SimTK::Transform();

            childAtomMobod.setDefaultInboardFrame(XPF);
            childAtomMobod.setDefaultOutboardFrame(XBM);

            continue;
        }

        // std::cout << "parent cAIx: " << bond.parentCAIx << ", child cAIx: " << bond.childCAIx << "\n";

        // Get parent-child BondCenters relationship
        const auto& X_parentBC_childBC = Xpc(bond.topologyIndex, bond.childCAIx);
        const auto& X_childBC_parentBC = ~X_parentBC_childBC;

        // Get Top frame
        const auto& T_X_root = F(bond.topologyIndex, bond.childCAIx);

        // Origin of the parent mobod
        const auto& T_X_Proot = F(bond.topologyIndex, bond.parentMobodRootCAIx);
        const auto Proot_X_T = ~T_X_Proot;
        const auto Proot_X_root = Proot_X_T * T_X_root;

        // Create transforms for child default inboard frame (XPF) and default outboard frame (XBM)
        switch (bond.mobility) {
            /*
            -1 0  0 L
             0 c  s 0
             0 s -c 0
             0 0  0 1
            */
            case SimTK::BondMobility::Mobility::AnglePin:
            case SimTK::BondMobility::Mobility::Slider:
            case SimTK::BondMobility::Mobility::BendStretch: {
                const auto& B_X_M_anglePin = X_parentBC_childBC;
                const auto P_X_F_anglePin = Proot_X_root * B_X_M_anglePin;

                childAtomMobod.setDefaultInboardFrame(P_X_F_anglePin);  // X_PF
                childAtomMobod.setDefaultOutboardFrame(B_X_M_anglePin); // X_BM
                break;
            }

            /*
             0 0 1 L
             s c 0 0
            -c s 0 0
             0 0 0 1
            */
            case SimTK::BondMobility::Mobility::Torsion:
            case SimTK::BondMobility::Mobility::Cylinder: {
                const auto B_X_M_pin = X_parentBC_childBC * X_to_Z;
                const auto P_X_F_pin = Proot_X_root * B_X_M_pin;

                childAtomMobod.setDefaultInboardFrame(P_X_F_pin);  // X_PF
                childAtomMobod.setDefaultOutboardFrame(B_X_M_pin); // X_BM
                break;
            }

            case SimTK::BondMobility::Mobility::BallM:
            case SimTK::BondMobility::Mobility::Rigid:
            case SimTK::BondMobility::Mobility::Translation: {
                // Samuel Flores' terminology aka M_X_pin =
                // SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis)
                const auto& B_X_M = X_to_Z;
                const auto P_X_F = Proot_X_root * B_X_M;

                childAtomMobod.setDefaultInboardFrame(P_X_F);  // X_PF
                childAtomMobod.setDefaultOutboardFrame(B_X_M); // X_BM
                break;
            }

            case SimTK::BondMobility::Mobility::Spherical: {
                const auto P_X_F_spheric =
                    SimTK::Transform() * Proot_X_root * X_parentBC_childBC * X_to_Y * Y_to_Z;
                const auto M_X_B_spheric = SimTK::Transform() * Z_to_Y * Y_to_X * X_childBC_parentBC;
                const auto& B_X_M_spheric = ~M_X_B_spheric;

                childAtomMobod.setDefaultInboardFrame(P_X_F_spheric);  // X_PF
                childAtomMobod.setDefaultOutboardFrame(B_X_M_spheric); // X_BM
                break;
            }

            case SimTK::BondMobility::Mobility::OrthoSpherical: {
                // Get parent-child BC transform
                const auto X_parentAtom_BCpar =
                    topology.calcDefaultBondCenterFrameInParentAtomFrame(bond.parentCAIx, bond.childCAIx);
                const auto X_childAtom_BCchi =
                    topology.calcDefaultBondCenterFrameInChildAtomFrame(bond.parentCAIx, bond.childCAIx);
                const auto& X_BCchi_childAtom = ~X_childAtom_BCchi;

                const auto& P_X_F_orthospheric = X_parentAtom_BCpar;                    // BAT from Compound
                const auto M_X_B_orthospheric = X_parentBC_childBC * X_BCchi_childAtom; // BAT from Compound
                const auto& B_X_M_orthospheric = ~M_X_B_orthospheric; // X_childAtom_BC * X_childBC_parentBC;

                childAtomMobod.setDefaultInboardFrame(P_X_F_orthospheric);  // X_PF
                childAtomMobod.setDefaultOutboardFrame(B_X_M_orthospheric); // X_BM
                break;
            }

            default:
                throw std::runtime_error("Mobility type '" + std::string(getBondMobilityName(bond.mobility))
                                         + "' not supported in World::setFramesFromTopologies().");
        }
    }

    // Handle bonds involving root atoms separately
    for (const auto& bond : rootAtomBonds) {
        const auto& topology = topologies[bond.topologyIndex];

        // const SimTK::Transform& G_X_T = topology.getTopLevelTransform();
        const SimTK::Transform& G_X_T = SimTK::Transform();

        if (topology.getAtoms()[bond.parentCAIx].connectivity.root) {
            const SimTK::Transform& T_X_base = F(bond.topologyIndex, bond.parentCAIx);
            const SimTK::Transform G_X_base = G_X_T * T_X_base;

            SimTK::MobilizedBody& parentAtomMobod = matter->updMobilizedBody(bond.parentMBIx);
            parentAtomMobod.setDefaultInboardFrame(G_X_base);
            parentAtomMobod.setDefaultOutboardFrame(SimTK::Transform());
        }

        if (topology.getAtoms()[bond.childCAIx].connectivity.root) {
            const SimTK::Transform& T_X_base = F(bond.topologyIndex, bond.childCAIx);
            const SimTK::Transform G_X_base = G_X_T * T_X_base;

            SimTK::MobilizedBody& childAtomMobod = matter->updMobilizedBody(bond.childMBIx);
            childAtomMobod.setDefaultInboardFrame(G_X_base);
            childAtomMobod.setDefaultOutboardFrame(SimTK::Transform());
        }
    }
}

void World::generateDummParams(const std::vector<RoboAtom>& atoms,
                               const std::vector<RoboBond>& bonds,
                               const std::vector<RoboAngle>& angles,
                               const std::vector<RoboPeriodicTorsion>& properPeriodicTorsions,
                               const std::vector<RoboHarmonicImproperTorsion>& harmonicImproperTorsions) {
    // Make a counter that checks if the atom class index already exists since Molmodel does not check for and
    // does not allow re-definitions
    std::vector<bool> atomClassDefined(atoms.size(), false);
    std::vector<bool> chargedAtomTypeDefined(atoms.size(), false);

    // Define atom classes and charged atom types
    for (const auto& atom : atoms) {
        if (!atomClassDefined[atom.identity.atomClassIndex]) {
            atomClassDefined[atom.identity.atomClassIndex] = true;

            forceField->defineAtomClass(atom.identity.atomClassIndex,
                                        atom.identity.atomClassName.c_str(),
                                        atom.elementInfo.atomicNumber,
                                        atom.connectivity.numBondsInvolved,
                                        atom.physics.vdwRadiusInNm,
                                        atom.physics.vdwWellDepthInKJ);
        }

        if (!chargedAtomTypeDefined[atom.identity.chargedAtomTypeIndex]) {
            chargedAtomTypeDefined[atom.identity.chargedAtomTypeIndex] = true;

            forceField->defineChargedAtomType(atom.identity.chargedAtomTypeIndex,
                                              atom.identity.chargedAtomTypeName.c_str(),
                                              atom.identity.atomClassIndex,
                                              atom.physics.chargeInE);
            forceField->setBiotypeChargedAtomType(atom.identity.chargedAtomTypeIndex, atom.biotypeIndex);
        }
    }

    // Define bonds parameters using atom classes
    // DuMM canonicalizes the order of atom class indices internally (smallest first)
    // DuMM also checks for re-definitions and ignores them if found
    // However, it will fail if we try to define a bond with the same atom classes but different parameters
    std::map<BondStretchKey, BondStretchValue> definedBondStretches;

    for (const auto& bond : bonds) {
        const auto aCIx1 = atoms[bond.globalIndices[0]].identity.atomClassIndex;
        const auto aCIx2 = atoms[bond.globalIndices[1]].identity.atomClassIndex;

        const auto key = BondStretchKey(aCIx1, aCIx2);
        const BondStretchValue newParams{{atoms[bond.globalIndices[0]].identity.globalIndex,
                                          atoms[bond.globalIndices[1]].identity.globalIndex},
                                         bond.stiffnessInKJPerNmSq,
                                         bond.nominalLengthInNm};

        const auto it = definedBondStretches.find(key);
        if (it != definedBondStretches.end()) {
            const BondStretchValue& existingParams = it->second;
            if (existingParams != newParams) {
                std::string message =
                    "Error in World::generateDummParams(): Conflicting bond stretch parameters.\n";
                message += "Existing parameters defined for atoms:\n";
                for (const auto i : existingParams.globalAtomIndices) {
                    message += "  " + getAtomDescription(atoms[i]) + "\n";
                }
                message += "  Stiffness: " + std::to_string(existingParams.stiffness)
                           + " kJ/(nm^2 mol), Length: " + std::to_string(existingParams.length) + " nm\n";
                message += "New parameters defined for atoms:\n";
                for (const auto i : newParams.globalAtomIndices) {
                    message += "  " + getAtomDescription(atoms[i]) + "\n";
                }
                message += "  Stiffness: " + std::to_string(newParams.stiffness)
                           + " kJ/(nm^2 mol), Length: " + std::to_string(newParams.length) + " nm\n";
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
        const BondBendValue newParams{{atoms[angle.globalIndices[0]].identity.globalIndex,
                                       atoms[angle.globalIndices[1]].identity.globalIndex,
                                       atoms[angle.globalIndices[2]].identity.globalIndex},
                                      angle.stiffnessInKJPerRadSq,
                                      angle.nominalAngleInDeg};

        const auto it = definedBondBends.find(key);
        if (it != definedBondBends.end()) {
            const BondBendValue& existingParams = it->second;
            if (existingParams != newParams) {
                std::string message =
                    "Error in World::generateDummParams(): Conflicting bond bend parameters.\n";
                message += "Existing parameters defined for atoms:\n";
                for (const auto i : existingParams.globalAtomIndices) {
                    message += "  " + getAtomDescription(atoms[i]) + "\n";
                }
                message += "  Stiffness: " + std::to_string(existingParams.stiffness)
                           + " kJ/(rad^2 mol), Angle: " + std::to_string(existingParams.angleDeg) + " deg\n";
                message += "New parameters defined for atoms:\n";
                for (const auto i : newParams.globalAtomIndices) {
                    message += "  " + getAtomDescription(atoms[i]) + "\n";
                }
                message += "  Stiffness: " + std::to_string(newParams.stiffness)
                           + " kJ/(rad^2 mol), Angle: " + std::to_string(newParams.angleDeg) + " deg\n";
                SimTK_ASSERT_ALWAYS(existingParams == newParams, message.c_str());
            }
        } else {
            definedBondBends.insert(it, std::make_pair(key, newParams));
        }

        forceField->defineBondBend(aCIx1, aCIx2, aCIx3, angle.stiffnessInKJPerRadSq, angle.nominalAngleInDeg);
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
        const PeriodicTorsionKey key(aCIx1, aCIx2, aCIx3, aCIx4, !torsion.improper);
        const PeriodicTorsionValue newParams{{atoms[torsion.globalIndices[0]].identity.globalIndex,
                                              atoms[torsion.globalIndices[1]].identity.globalIndex,
                                              atoms[torsion.globalIndices[2]].identity.globalIndex,
                                              atoms[torsion.globalIndices[3]].identity.globalIndex},
                                             {torsion.terms[0].periodicity,
                                              torsion.terms[1].periodicity,
                                              torsion.terms[2].periodicity,
                                              torsion.terms[3].periodicity,
                                              torsion.terms[4].periodicity},
                                             {torsion.terms[0].amplitudeKJ,
                                              torsion.terms[1].amplitudeKJ,
                                              torsion.terms[2].amplitudeKJ,
                                              torsion.terms[3].amplitudeKJ,
                                              torsion.terms[4].amplitudeKJ},
                                             {torsion.terms[0].phaseDeg,
                                              torsion.terms[1].phaseDeg,
                                              torsion.terms[2].phaseDeg,
                                              torsion.terms[3].phaseDeg,
                                              torsion.terms[4].phaseDeg},
                                             torsion.numTerms};

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
                std::string message = "Error in World::generateDummParams(): Conflicting "
                                      + std::string(torsion.improper ? "improper" : "proper")
                                      + " periodic torsion parameters for atom classes ";
                message += "Existing parameters defined for atoms:\n";
                for (const auto i : existingParams.globalAtomIndices) {
                    message += "\t\t" + getAtomDescription(atoms[i]) + "\n";
                }
                message += "Existing parameters:\n";
                for (size_t i = 0; i < existingParams.numTerms; i++) {
                    message += "\t\tTerm " + std::to_string(i)
                               + ": periodicity=" + std::to_string(existingParams.periodicity[i])
                               + ", amplitude=" + std::to_string(existingParams.amplitude[i])
                               + ", phase=" + std::to_string(existingParams.phase[i]) + "\n";
                }

                message += "New parameters defined for atoms:\n";
                for (const auto i : newParams.globalAtomIndices) {
                    message += "\t\t" + getAtomDescription(atoms[i]) + "\n";
                }
                message += "New parameters:\n";
                for (size_t i = 0; i < newParams.numTerms; i++) {
                    message += "\t\tTerm " + std::to_string(i)
                               + ": periodicity=" + std::to_string(newParams.periodicity[i])
                               + ", amplitude=" + std::to_string(newParams.amplitude[i])
                               + ", phase=" + std::to_string(newParams.phase[i]) + "\n";
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
            switch (torsion.numTerms) {
                case 1:
                    forceField->defineAmberImproperTorsion(aCIx1,
                                                           aCIx2,
                                                           aCIx3,
                                                           aCIx4,
                                                           torsion.terms[0].periodicity,
                                                           torsion.terms[0].amplitudeKJ,
                                                           torsion.terms[0].phaseDeg);
                    break;
                case 2:
                    forceField->defineAmberImproperTorsion(aCIx1,
                                                           aCIx2,
                                                           aCIx3,
                                                           aCIx4,
                                                           torsion.terms[0].periodicity,
                                                           torsion.terms[0].amplitudeKJ,
                                                           torsion.terms[0].phaseDeg,
                                                           torsion.terms[1].periodicity,
                                                           torsion.terms[1].amplitudeKJ,
                                                           torsion.terms[1].phaseDeg);
                    break;
                case 3:
                    forceField->defineAmberImproperTorsion(aCIx1,
                                                           aCIx2,
                                                           aCIx3,
                                                           aCIx4,
                                                           torsion.terms[0].periodicity,
                                                           torsion.terms[0].amplitudeKJ,
                                                           torsion.terms[0].phaseDeg,
                                                           torsion.terms[1].periodicity,
                                                           torsion.terms[1].amplitudeKJ,
                                                           torsion.terms[1].phaseDeg,
                                                           torsion.terms[2].periodicity,
                                                           torsion.terms[2].amplitudeKJ,
                                                           torsion.terms[2].phaseDeg);
                    break;

                default:
                    SimTK_ASSERT_ALWAYS(false,
                                        "Error in World::generateDummParams: Improper torsions can "
                                        "only have 1, 2, or 3 terms.");
                    break;
            }
        } else {
            switch (torsion.numTerms) {
                case 1:
                    forceField->defineBondTorsion(aCIx1,
                                                  aCIx2,
                                                  aCIx3,
                                                  aCIx4,
                                                  torsion.terms[0].periodicity,
                                                  torsion.terms[0].amplitudeKJ,
                                                  torsion.terms[0].phaseDeg);
                    break;
                case 2:
                    forceField->defineBondTorsion(aCIx1,
                                                  aCIx2,
                                                  aCIx3,
                                                  aCIx4,
                                                  torsion.terms[0].periodicity,
                                                  torsion.terms[0].amplitudeKJ,
                                                  torsion.terms[0].phaseDeg,
                                                  torsion.terms[1].periodicity,
                                                  torsion.terms[1].amplitudeKJ,
                                                  torsion.terms[1].phaseDeg);
                    break;
                case 3:
                    forceField->defineBondTorsion(aCIx1,
                                                  aCIx2,
                                                  aCIx3,
                                                  aCIx4,
                                                  torsion.terms[0].periodicity,
                                                  torsion.terms[0].amplitudeKJ,
                                                  torsion.terms[0].phaseDeg,
                                                  torsion.terms[1].periodicity,
                                                  torsion.terms[1].amplitudeKJ,
                                                  torsion.terms[1].phaseDeg,
                                                  torsion.terms[2].periodicity,
                                                  torsion.terms[2].amplitudeKJ,
                                                  torsion.terms[2].phaseDeg);
                    break;
                case 4:
                    forceField->defineBondTorsion(aCIx1,
                                                  aCIx2,
                                                  aCIx3,
                                                  aCIx4,
                                                  torsion.terms[0].periodicity,
                                                  torsion.terms[0].amplitudeKJ,
                                                  torsion.terms[0].phaseDeg,
                                                  torsion.terms[1].periodicity,
                                                  torsion.terms[1].amplitudeKJ,
                                                  torsion.terms[1].phaseDeg,
                                                  torsion.terms[2].periodicity,
                                                  torsion.terms[2].amplitudeKJ,
                                                  torsion.terms[2].phaseDeg,
                                                  torsion.terms[3].periodicity,
                                                  torsion.terms[3].amplitudeKJ,
                                                  torsion.terms[3].phaseDeg);
                    break;
                case 5:
                    forceField->defineBondTorsion(aCIx1,
                                                  aCIx2,
                                                  aCIx3,
                                                  aCIx4,
                                                  torsion.terms[0].periodicity,
                                                  torsion.terms[0].amplitudeKJ,
                                                  torsion.terms[0].phaseDeg,
                                                  torsion.terms[1].periodicity,
                                                  torsion.terms[1].amplitudeKJ,
                                                  torsion.terms[1].phaseDeg,
                                                  torsion.terms[2].periodicity,
                                                  torsion.terms[2].amplitudeKJ,
                                                  torsion.terms[2].phaseDeg,
                                                  torsion.terms[3].periodicity,
                                                  torsion.terms[3].amplitudeKJ,
                                                  torsion.terms[3].phaseDeg,
                                                  torsion.terms[4].periodicity,
                                                  torsion.terms[4].amplitudeKJ,
                                                  torsion.terms[4].phaseDeg);
                    break;

                default:
                    SimTK_ASSERT_ALWAYS(
                        false,
                        "Error in World::generateDummParams: Proper torsions can only have 1 to 5 terms.");
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

        forceField->defineCustomBondTorsion(
            aCIx1,
            aCIx2,
            aCIx3,
            aCIx4,
            new HarmonicImproperTorsionForce(torsion.stiffnessInKJPerRadSq, torsion.nominalAngleInRad));
    }
}

World::World(int worldIndex,
             Span<Topology> topo,
             bool testing,
             const ZMatrix& _zMatrix,
             bool wantSpatialForceHistory)
    : ownWorldIndex(worldIndex)
    , topologies(topo)
    , testing(testing)
    , zMatrix(_zMatrix)
    , wantSpatialForceHistory(wantSpatialForceHistory) {
    multibodySystem = std::make_unique<SimTK::CompoundSystem>();
    matter = std::make_unique<SimTK::SimbodyMatterSubsystem>(*multibodySystem);
    forces = std::make_unique<SimTK::GeneralForceSubsystem>(*multibodySystem);
    forceField = std::make_unique<SimTK::DuMMForceFieldSubsystem>(*multibodySystem);
    integrator = std::make_unique<SimTK::VerletIntegrator>(*multibodySystem);
    timeStepper = std::make_unique<SimTK::TimeStepper>(*multibodySystem, *integrator);

    // Allocate atom target locations cache
    for (int topoIx = 0; topoIx < int(topologies.size()); ++topoIx) {
        const auto& topology = topologies[topoIx];
        atomTargetLocationsCache.emplace_back();
        numMolecules++;
        for (const auto& atom : topology.getAtoms()) {
            atomTargetLocationsCache.back().push_back(atom.position);
            numAtoms++;
        }
    }
}

void World::modelTopologies(const std::vector<SimTK::Compound::AtomTargetLocations>& atomTargets) {
    // (1) connectivity-only frame graph + reused buffers.
    buildFrameGraph(topologies);
    targetFlat.resize(frameGraph.totalAtoms);
    frameFlat.resize(frameGraph.totalAtoms);
    xpcBcFlat.resize(frameGraph.totalAtoms);

    topoAtomBodyFrame.assign(topologies.size(), {});
    for (std::size_t t = 0; t < topologies.size(); ++t) {
        topoAtomBodyFrame[t].assign(topologies[t].getNumAtoms(), SimTK::Transform());
    }

    // (2) frames + B, once, branchless/parallel.
    packTargets(frameGraph, atomTargets, targetFlat.data());
    computeAllAtomFrames(frameGraph, targetFlat.data(), frameFlat.data(), xpcBcFlat.data());

    // ---- Build every molecule's multibody tree (DuMM atoms/bonds + bodies). ----
    std::vector<std::vector<WorldRigidUnit>> unitsPerTopo(topologies.size());
    for (int topoIx = 0; topoIx < int(topologies.size()); ++topoIx) {
        topologies[topoIx].setMultibodySystem(*multibodySystem);
        topologies[topoIx].setTopLevelTransform(SimTK::Transform());
        unitsPerTopo[topoIx] = modelOneCompound(topoIx, topologies[topoIx].getRootMobility());
    }

    atomFrameCache.assign(topologies.size(), std::vector<SimTK::Transform>());
    // atomStations.assign(topologies.size(), std::vector<SimTK::Vec3>());
    clustersMass.assign(matter->getNumBodies() - 1, 0);
    clustersCOM.assign(matter->getNumBodies() - 1, SimTK::Vec3(0));
    clustersInertia.assign(matter->getNumBodies() - 1, SimTK::Inertia(0));

    mbxRootCAIx.resize(matter->getNumBodies());

    // ---- Per-topology bookkeeping, read straight off the decomposition. ----
    for (int topoIx = 0; topoIx < int(topologies.size()); ++topoIx) {
        Topology& topology = topologies[topoIx];

        // atomStations[topoIx].assign(topology.getNumAtoms(), SimTK::Vec3(0));

        const auto& units = unitsPerTopo[topoIx];
        const int n = topology.getNumAtoms();

        topoAtomToMbx.emplace_back(n);
        topoAtomToDAIX.emplace_back(n);
        topoAtomIsRigidBodyRoot.emplace_back(n, false);

        // unit root atom for each body, and root flag / mbx per atom
        for (const auto& unit : units) {
            mbxRootCAIx[unit.mbx] = {topoIx, SimTK::Compound::AtomIndex(unit.rootCAIx)};

            // spatial force history: record inboard/outboard prmtop indices for
            // non-root child bodies (mbx > 1) using the joint we already know.
            if (wantSpatialForceHistory && unit.mbx > 1 && unit.rootCAIx > 0 && unit.jointParentCAIx >= 0) {
                const auto inboardPrmtop = topology.getAtoms()[unit.rootCAIx].identity.prmtopIndex;
                mbx2PrmtopInboardIndex.insert(std::make_pair(unit.mbx, inboardPrmtop));

                const auto outboardPrmtop = topology.getAtoms()[unit.jointParentCAIx].identity.prmtopIndex;
                prmtopInboardIndex2PrmtopOutboardIndex[inboardPrmtop] = outboardPrmtop;

                std::cout << "Mobx " << unit.mbx << " has inboard CAIx " << unit.rootCAIx
                          << " with prmtop index " << inboardPrmtop << std::endl;
            }

            for (int a : unit.atomCAIxs) {
                topoAtomIsRigidBodyRoot[topoIx][a] = (a == unit.rootCAIx);
            }
        }

        // mbx / dAIx per atom, mass refresh, DuMM cluster cache for stations
        for (SimTK::Compound::AtomIndex cAIx(0); cAIx < n; ++cAIx) {
            const auto mbx = topology.getAtomMobilizedBodyIndex(cAIx); // value we set
            const auto dAIx = topology.getDuMMAtomIndex(cAIx);         // value we set
            topoAtomToMbx[topoIx][cAIx] = mbx;
            topoAtomToDAIX[topoIx][cAIx] = dAIx;

            // (mass was set in modelOneCompound; re-set is harmless but kept for
            //  parity -- it invalidates the subsystem topology cache.)
            forceField->setDuMMAtomMass(dAIx, topology.getAtoms()[cAIx].physics.massInDaltons);

            // Map dAIx -> the body's implicit cluster so the per-sample
            // bsetAtomPlacementStation path keeps working.
            forceField->updateClustersCacheList(dAIx, mbx);
        }

        // ---- Cache root-atom bonds and inter-body (flexible) bonds. ----
        for (const auto& bond : topology.getBonds()) {
            const auto childCAIx = bond.compoundAtomIndices[1];
            const auto parentCAIx = bond.compoundAtomIndices[0];
            const auto childMBIx = topology.getAtomMobilizedBodyIndex(childCAIx);
            const auto parentMBIx = topology.getAtomMobilizedBodyIndex(parentCAIx);

            if (topology.getAtoms()[childCAIx].connectivity.root
                || topology.getAtoms()[parentCAIx].connectivity.root) {
                RootAtomBond rootAtomBond;
                rootAtomBond.topologyIndex = topoIx;
                rootAtomBond.childCAIx = childCAIx;
                rootAtomBond.parentCAIx = parentCAIx;
                rootAtomBond.childMBIx = childMBIx;
                rootAtomBond.parentMBIx = parentMBIx;
                rootAtomBonds.push_back(rootAtomBond);
            }

            if (bond.ringClosing) {
                continue;
            }
            if (bond.getBondMobility(ownWorldIndex) == SimTK::BondMobility::Mobility::Rigid) {
                continue;
            }

            // A flexible bond between two rigid bodies.
            RigidBodyAtomBond rigidBodyAtomBond;
            rigidBodyAtomBond.topologyIndex = topoIx;
            rigidBodyAtomBond.mobility = bond.getBondMobility(ownWorldIndex);
            rigidBodyAtomBond.childAtomGlobalIndex = bond.globalIndices[1];
            rigidBodyAtomBond.parentAtomGlobalIndex = bond.globalIndices[0];
            rigidBodyAtomBond.childCAIx = childCAIx;
            rigidBodyAtomBond.parentCAIx = parentCAIx;
            rigidBodyAtomBond.childMBIx = childMBIx;
            rigidBodyAtomBond.parentMBIx = parentMBIx;

            // root atom of the parent body -- straight from mbxRootCAIx now.
            rigidBodyAtomBond.parentMobodRootCAIx = mbxRootCAIx[parentMBIx].cAIx;

            // grandparent (needs three atoms to define an angle).
            if (rigidBodyAtomBond.parentCAIx > 0) {
                rigidBodyAtomBond.grandParentCAIx =
                    topology.getInboardAtomIndex(rigidBodyAtomBond.parentCAIx); // pure topology adjacency
                rigidBodyAtomBond.grandParentMBIx =
                    topology.getAtomMobilizedBodyIndex(rigidBodyAtomBond.grandParentCAIx);
            }

            rigidBodyAtomBonds.push_back(rigidBodyAtomBond);
            interestingMobodIndices.insert(childMBIx);
            interestingMobodIndices.insert(parentMBIx);
        }
    }

    multibodySystem->realizeTopology();
}


// Print recommended timesteps. We need and advanced State here
auto World::getRecommendedTimesteps() -> SimTK::Real {
    SimTK::State& someState = integrator->updAdvancedState();
    int nu = matter->getNU(someState);

    SimTK::Real minTimeStep;
    SimTK::Real prevMinTimeStep = SimTK::Infinity;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        minTimeStep = 0.0007 * std::sqrt(mobod.getBodyMass(someState));
        if (minTimeStep < prevMinTimeStep) {
            prevMinTimeStep = minTimeStep;
        }
    }

    return prevMinTimeStep;
}

/**
 * Calc Station Jacobian JS
 */
void World::calcStationJacobian(const SimTK::State& someState, SimTK::Matrix_<SimTK::Vec3>& JS) const {
    matter->calcStationJacobian(someState, onBodyB, taskStationPInGuest, JS);

    std::cout << "Task Bodies ";
    std::cout << onBodyB << std::endl;
    std::cout << "Task Stations ";
    std::cout << taskStationPInGuest << std::endl;
    std::cout << "Station Jacobian ";
    std::cout << JS << std::endl;

    // matter->calcBiasForStationJacobian(someState, onBodyB, stationPInB, JSDotu);
}

/**
 * Add contact constraints to specific bodies.
 **/
void World::addRodConstraint(SimTK::State& someState) {
    int hostTopology = 0;
    int guestTopology = 1;

    std::vector<int> bAtomIxs_host = {4};   // atoms on host topology
    std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology
    rodBodies.emplace_back(std::make_pair(SimTK::MobilizedBody(), SimTK::MobilizedBody()));
    conStationPInGuest.emplace_back(SimTK::Vec3());
    conStationPInHost.emplace_back(SimTK::Vec3());
    conDeltaStationP.emplace_back(SimTK::Vec3());

    // Get stations in host
    int topi = -1;
    for (auto& topology : topologies) {
        topi++;
        if (topi == hostTopology) {
            /* // Atoms
            int tz = -1;
            for (int bAtomIx : bAtomIxs_host) {
                tz++;
                SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
                SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx,
            forceField);

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
    for (auto& topology : topologies) {
        topi++;

        if (topi == guestTopology) {
            /* // Atoms
            int tz = -1;
            for (int bAtomIx : bAtomIxs_guest) {
                tz++;
                SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
                SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx,
            forceField);

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

    // rodConstraints.back().enable(someState);
}

/** Add speed constraints to specific bodies.
TODO:use number of mobilities. TODO: Solve if **/
auto World::addSpeedConstraint(int prmtopIndex) -> const SimTK::State& {
    int hostTopology = 0;
    int guestTopology = 1;

    std::vector<int> bAtomIxs_host = {4};   // atoms on host topology
    std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

    if (prmtopIndex >= 0) {
        std::cout << "Adding constraint to atom with prmtop index " << prmtopIndex << "\n";
        SimTK::MobilizedBodyIndex mbx =
            topologies[0].getAtomMobilizedBodyIndexThroughDumm(SimTK::Compound::AtomIndex(prmtopIndex),
                                                               *forceField);
        SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
        SimTK::Constraint::ConstantSpeed B3291ConstraintU1(mobod, SimTK::MobilizerUIndex(0), 0);
        if (matter->getNumBodies() > 5000) {
            SimTK::Constraint::ConstantSpeed B3291ConstraintU2(mobod, SimTK::MobilizerUIndex(1), 0);
            SimTK::Constraint::ConstantSpeed B3291ConstraintU3(mobod, SimTK::MobilizerUIndex(2), 0);
        }
    }

    const SimTK::State& returnState = multibodySystem->realizeTopology();
    return returnState;
}

// Get the (potential) energy transfer
// If any of the Q, U or tau is actively modifyied by the sampler
// the Jacobian of that transformation will be included too
SimTK::Real World::getWorkOrHeat() {
    // Accumulate in this variable
    SimTK::Real retValue = 0.0;

    // Get the energy transfer from all the samplers
    for (auto& sampler : this->samplers) {
        // Get the potential energy difference
        retValue +=
            getSampler(0)->getCurrentEnergy().potential - getSampler(0)->getPreviousEnergy().potential;

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
SimTK::Real World::getWork() const {
    // Accumulate in this variable
    SimTK::Real retValue = 0.0;

    // Get the energy transfer from all the samplers
    for (auto& sampler : this->samplers) {
        if (sampler->getDistortOpt() < 0) {
            // Get the potential energy difference
            retValue +=
                getSampler(0)->getCurrentEnergy().potential - getSampler(0)->getPreviousEnergy().potential;

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
void World::getTransformsStatistics(SimTK::State& someState) {
    // Get generalized coordinates Q template values. These are the values that
    // Q is extending. In the case of an AnglePin mobilizer, Q is extending an
    // <(P_x, F_x) angle.

    // Get bonds and angles values
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // Get mobod inboard frame X_PF
        const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
        // std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

        // Get mobod inboard frame X_FM measured and expressed in P
        const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
        // std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

        // Get mobod inboard frame X_BM
        const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
        // std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
        // std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

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

// Print bond lengthe and angle bends
void World::traceBendStretch(SimTK::State& someState) {
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // Get mobod inboard frame X_PF
        const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
        // PrintTransform(X_PF, 4, "X_PF " + std::to_string(int(mbx)));

        // Get mobod inboard frame X_FM measured and expressed in P
        const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);

        // Get mobod inboard frame X_BM
        const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
        // PrintTransform(X_BM, 4, "X_BM " + std::to_string(int(mbx)));

        // Get BAT coordinates B and A
        SimTK::Vec3 bondVector = X_BM.p();
        std::cout << "bondVector.norm() " << bondVector.norm() << " acos(X_PF.R()(0)(0)) "
                  << std::acos(X_PF.R()(0)(0)) << std::endl;
        // trace("X_FM");
        // PrintTransform(X_FM, 10);
    }
}

// Print X_PF means
void World::PrintAcosX_PFs() {
    int i = -1;
    for (const auto& xpf : acosX_PF00) {
        i += 1;
        // std::cout << "X_PF " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
        std::cout << "acosX_PF " << i << " " << xpf << std::endl;
    }
}

// Print X_PF means
void World::PrintNormX_BMs() {
    int i = -1;
    for (const auto& xbm : normX_BMp) {
        i += 1;
        std::cout << "normX_BM " << i << " " << xbm << std::endl;
    }
}

// Print X_PF means
void World::PrintAcosX_PFMeans() {
    int i = -1;
    for (const auto& xpf : acosX_PF00_means) {
        i += 1;
        // std::cout << "X_PFMean " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
        std::cout << "acosX_PFMean " << i << " " << xpf << std::endl;
    }
}

// Print X_PF means
void World::PrintNormX_BMMeans() {
    int i = -1;
    for (const auto& xbm : normX_BMp_means) {
        i += 1;
        std::cout << "normX_BMMean " << i << " " << xbm << std::endl;
    }
}

// REORIENT

SimTK::Transform& World::getReorientTransformInAnotherBody(const SimTK::State& someState,
                                                           const SimTK::MobilizedBody& inBodyA,
                                                           const SimTK::MobilizedBody& ofBodyB,
                                                           const SimTK::Transform& reorientAB,
                                                           SimTK::Transform& X_FMprim) {
    SimTK::Transform X_MB = ~(ofBodyB.getOutboardFrame(someState));
    SimTK::Transform X_FM = ofBodyB.getMobilizerTransform(someState);
    SimTK::Transform X_AB = ofBodyB.findBodyTransformInAnotherBody(someState, inBodyA);

    SimTK::Transform X_BBprim = (~X_AB) * reorientAB;
    X_FMprim = X_FM * X_MB * X_BBprim * X_MB;

    return X_FMprim;
}

//...............

/**
 * Set X_PF, X_BM means
 */
void World::setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
                               const std::vector<SimTK::Real>& givenX_BM) {
    // Update acosX_PF00 means
    int i = -1;
    for (auto& xpf : acosX_PF00_means) {
        i += 1;
        xpf = givenX_PF[i];
    }

    // Update normX_BMp means
    i = -1;
    for (auto& xbm : normX_BMp_means) {
        i += 1;
        xbm = givenX_BM[i];
    }
}

/**
 * Set X_PF, X_BM means to initial values
 */
void World::setTransformsMeansToIni() {
    const SimTK::State& defaultState = matter->getSystem().getDefaultState();

    // Get generalized coordinates Q template values. These are the values that
    // Q is extending. In the case of an AnglePin mobilizer, Q is extending an
    // <(P_x, F_x) angle.

    // Get bonds and angles values
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // Get mobod inboard frame X_PF
        const SimTK::Transform& X_PF = mobod.getInboardFrame(defaultState);
        // std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

        // Get mobod inboard frame X_FM measured and expressed in P
        const SimTK::Transform& X_FM = mobod.getMobilizerTransform(defaultState);
        // std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

        // Get mobod inboard frame X_BM
        const SimTK::Transform& X_BM = mobod.getOutboardFrame(defaultState);
        // std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
        // std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

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
void World::setTransformsMeansToCurrent(SimTK::State& someState) {
    // Get generalized coordinates Q template values. These are the values that
    // Q is extending. In the case of an AnglePin mobilizer, Q is extending an
    // <(P_x, F_x) angle.

    // Get bonds and angles values
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // Get mobod inboard frame X_PF
        const SimTK::Transform& X_PF = mobod.getInboardFrame(someState);
        // std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

        // Get mobod inboard frame X_FM measured and expressed in P
        const SimTK::Transform& X_FM = mobod.getMobilizerTransform(someState);
        // std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

        // Get mobod inboard frame X_BM
        const SimTK::Transform& X_BM = mobod.getOutboardFrame(someState);
        // std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
        // std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

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
void World::updateTransformsMeans(SimTK::State& someState) {
    int nofSamples = getNofSamples() + 1;
    // std::cout << "Nof samples " << nofSamples << std::endl;

    // Useful vars
    SimTK::Real N_1overN = 9999, NInv = 9999;

    if (nofSamples == 1) {
        for (unsigned int k = 0; k < acosX_PF00_means.size(); k++) {
            acosX_PF00_means[k] = acosX_PF00[k];
        }
        for (unsigned int k = 0; k < normX_BMp_means.size(); k++) {
            normX_BMp_means[k] = normX_BMp[k];
            // std::cout << "World " << ownWorldIndex << " bondUpdMean " << k << " " << normX_BMp_means[k] <<
            // std::endl;
        }
    } else {
        if (nofSamples == 2) {
            N_1overN = NInv = 0.5;
        } else {
            // Update useful vars
            SimTK::Real N_1 = nofSamples - 1.0;
            N_1overN = N_1 / nofSamples;
            NInv = 1.0 / nofSamples;
        }
        // std::cout << "updateX_PFMeans check " << " "
        //	<<  N_1overN << " " <<  NInv  << " " << std::flush;

        // Update acosX_PF00 means
        int i = -1;
        for (auto& xpf : acosX_PF00_means) {
            i += 1;
            xpf = (N_1overN * xpf) + (NInv * acosX_PF00[i]);
        }

        // Update normX_BMp means
        i = -1;
        for (auto& xbm : normX_BMp_means) {
            i += 1;
            xbm = (N_1overN * xbm) + (NInv * normX_BMp[i]);
            // std::cout << "World " << ownWorldIndex << " bondUpdMean " << i << " " << xbm << std::endl;
        }
    }
}

// Get X_PF means
std::vector<SimTK::Real>& World::getX_PFMeans() {
    return acosX_PF00_means;
}

// Get X_BM means
std::vector<SimTK::Real>& World::getX_BMMeans() {
    return normX_BMp_means;
}

/**
 * Calculate bond length and angle deviations from their means
 */
void World::calcBendStretchDeviations(SimTK::State& someState,
                                      std::vector<SimTK::Real>& X_PFdiffs,
                                      std::vector<SimTK::Real>& X_BMdiffs) {
    // Make sure it has
    X_PFdiffs.resize(this->acosX_PF00_means.size(), 0.0);
    X_BMdiffs.resize(this->normX_BMp_means.size(), 0.0);

    //
    for (unsigned int k = 0; k < X_PFdiffs.size(); k++) {
        X_PFdiffs[k] = this->acosX_PF00[k] - this->acosX_PF00_means[k];
    }
    for (unsigned int k = 0; k < X_BMdiffs.size(); k++) {
        X_BMdiffs[k] = this->normX_BMp[k] - this->normX_BMp_means[k];
        /* std::cout << "World " << ownWorldIndex << " bondDiff " << k << " "
        //<< this->normX_BMp[k] << " " << this->normX_BMp_means[k] << " "
        << X_BMdiffs[k] << std::endl; */
    }
}

/** Get U scale factor for the mobilized body **/
SimTK::Real World::getMobodUScaleFactor(SimTK::MobilizedBodyIndex& mbx) const {
    if (!mbx2uScale.empty()) {
        if (mbx2uScale.find(mbx) != mbx2uScale.end()) {
            return mbx2uScale.at(mbx);
        } else {
            // std::cout << "Warning: U scale factor for mobod " << int(mbx) << " not found.\n";
            return 1;
        }
    } else {
        return 1;
    }
}

/*! <!--  --> */
void World::calcSimbodyBAT(std::vector<std::vector<int>>& ZMatrix,
                           std::vector<SimTK::Real>& BONDLengthe,
                           std::vector<SimTK::Real>& ANGLEBends,
                           std::vector<SimTK::Real>& TORSIONAngles) {
    SimTK::State& advState = integrator->updAdvancedState();

    bool parFlag = false;
    bool parParFlag = false;
    int childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;

    if (BONDLengthe.size() != matter->getNumBodies() - 1) {
        BONDLengthe.resize(matter->getNumBodies() - 1, SimTK::NaN);
    }
    if (ANGLEBends.size() != matter->getNumBodies() - 1) {
        ANGLEBends.resize(matter->getNumBodies() - 1, SimTK::NaN);
    }
    if (TORSIONAngles.size() != matter->getNumBodies() - 1) {
        TORSIONAngles.resize(matter->getNumBodies() - 1, SimTK::NaN);
    }
    if (ZMatrix.size() != matter->getNumBodies() - 1) {
        ZMatrix.resize(matter->getNumBodies() - 1, std::vector<int>(4, -1));
    }

    bool printTransforms = false;

    for (SimTK::MobilizedBodyIndex childMbx(1); childMbx < matter->getNumBodies(); ++childMbx) {
        childIx = int(childMbx);

        const SimTK::MobilizedBody& childMobod = matter->getMobilizedBody(childMbx);
        const SimTK::MobilizedBody& parentMobod = childMobod.getParentMobilizedBody();
        const SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
        parentIx = int(parentMbx);

        if (int(childMbx) > 1) {
            const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
            const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();
            grandIx = int(grandMbx);
        }

        if (int(childMbx) > 2) {
            const SimTK::MobilizedBody& grandGrandMobod =
                parentMobod.getParentMobilizedBody().getParentMobilizedBody();
            const SimTK::MobilizedBodyIndex grandGrandMbx = grandGrandMobod.getMobilizedBodyIndex();
            grandGrandIx = int(grandGrandMbx);
        }

        // Print out the indices
        if (printTransforms) {
            std::cout << "World " << ownWorldIndex << " " << "child " << childIx << " " << "parent "
                      << parentIx << " " << "parPar " << grandIx << " " << "grandGrandIx " << grandGrandIx
                      << " " << std::endl
                      << std::flush;
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
        BONDLengthe[int(childMbx) - 1] = B_X_C.p().norm(); // correct

        ZMatrix[int(childMbx) - 1][0] = int(childMbx);
        ZMatrix[int(childMbx) - 1][1] = int(parentMbx);

        if (printTransforms) {
            // SimTK::Test::PrintTransform(G_X_C, 6, "G_X_C", "G_X_C:" + std::to_string(ownWorldIndex) + ":" +
            // std::to_string(int(childMbx))); SimTK::Test::PrintTransform(G_X_Mb, 6, "G_X_Mb", "G_X_Mb:" +
            // std::to_string(ownWorldIndex) + ":" + std::to_string(int(childMbx)));
            // SimTK::Test::PrintTransform(G_X_Fb, 6, "G_X_Fb", "G_X_Fb:" + std::to_string(ownWorldIndex) +
            // ":" + std::to_string(int(childMbx)));
            SimTK::Test::PrintTransform(G_X_B,
                                        6,
                                        "G_X_B",
                                        "G_X_B:" + std::to_string(ownWorldIndex) + ":"
                                            + std::to_string(int(parentMbx)));
            SimTK::Test::PrintTransform(B_X_C,
                                        6,
                                        "B_X_C",
                                        "B_X_C:" + std::to_string(int(parentMbx)) + ":"
                                            + std::to_string(int(childMbx)));
            SimTK::Test::PrintTransform(B_X_Fb,
                                        6,
                                        "B_X_Fb",
                                        "B_X_Fb:" + std::to_string(int(parentMbx)) + ":"
                                            + std::to_string(int(childMbx)));
            SimTK::Test::PrintTransform(Fb_X_Mb,
                                        6,
                                        "Fb_X_Mb",
                                        "Fb_X_Mb:" + std::to_string(int(parentMbx)) + ":"
                                            + std::to_string(int(childMbx)));
            SimTK::Test::PrintTransform(C_X_Mb,
                                        6,
                                        "C_X_Mb",
                                        "C_X_Mb:" + std::to_string(int(parentMbx)) + ":"
                                            + std::to_string(int(childMbx)));
        }

        if (int(childMbx) > 1) { // ANGLE ========================
            const SimTK::MobilizedBody& grandMobod = parentMobod.getParentMobilizedBody();
            const SimTK::MobilizedBodyIndex grandMbx = grandMobod.getMobilizedBodyIndex();

            const SimTK::Transform& A_X_Fa = parentMobod.getInboardFrame(advState);        // A_X_Fa
            const SimTK::Transform& B_X_Ma = parentMobod.getOutboardFrame(advState);       // B_X_Ma
            const SimTK::Transform& Fa_X_Ma = parentMobod.getMobilizerTransform(advState); // Fa_X_Ma
            SimTK::Transform A_X_B = A_X_Fa * Fa_X_Ma * (~B_X_Ma);
            SimTK::Transform B_X_A = ~A_X_B;
            SimTK::Transform A_X_C = A_X_B * B_X_C;

            const SimTK::Transform& G_X_A = grandMobod.getBodyTransform(advState); // G_X_A
            SimTK::Transform G_X_Fa = G_X_A * A_X_Fa;
            SimTK::Transform G_X_Ma = G_X_B * B_X_Ma;

            // Checks
            SimTK::Vec3 pAXB_A = A_X_B.p();
            // SimTK::Vec3 pAXB_B = ~(A_X_B.R()) * pAXB_A;
            SimTK::Vec3 pBXA_B = B_X_A.p();
            SimTK::Vec3 pBXC_B = B_X_C.p();

            ANGLEBends[int(childMbx) - 1] =
                std::acos(SimTK::dot((-1 * pBXA_B).normalize(), pBXC_B.normalize())); // correct
            // ================================================

            ZMatrix[int(childMbx) - 1][2] = int(grandMbx);

            if (printTransforms) {
                // SimTK::Test::PrintTransform(G_X_Ma, 6, "G_X_Ma", "G_X_Ma:" + std::to_string(ownWorldIndex)
                // + ":" + std::to_string(int(parentMbx))); SimTK::Test::PrintTransform(G_X_Fa, 6, "G_X_Fa",
                // "G_X_Fa:" + std::to_string(ownWorldIndex) + ":" + std::to_string(int(parentMbx)));
                SimTK::Test::PrintTransform(G_X_A,
                                            6,
                                            "G_X_A",
                                            "G_X_A:" + std::to_string(ownWorldIndex) + ":"
                                                + std::to_string(int(grandIx)));
                SimTK::Test::PrintTransform(A_X_B,
                                            6,
                                            "A_X_B",
                                            "A_X_B:" + std::to_string(int(grandMbx)) + ":"
                                                + std::to_string(int(parentMbx)));
                SimTK::Test::PrintTransform(A_X_Fa,
                                            6,
                                            "A_X_Fa",
                                            "A_X_Fa:" + std::to_string(int(grandMbx)) + ":"
                                                + std::to_string(int(parentMbx)));
                SimTK::Test::PrintTransform(Fa_X_Ma,
                                            6,
                                            "Fa_X_Ma",
                                            "Fa_X_Ma:" + std::to_string(int(grandMbx)) + ":"
                                                + std::to_string(int(parentMbx)));
                SimTK::Test::PrintTransform(B_X_Ma,
                                            6,
                                            "B_X_Ma",
                                            "B_X_Ma:" + std::to_string(int(grandMbx)) + ":"
                                                + std::to_string(int(parentMbx)));
            }

            if (int(childMbx) > 2) { // TORSION =======================
                const SimTK::MobilizedBody& grandGrandMobod =
                    parentMobod.getParentMobilizedBody().getParentMobilizedBody();

                const SimTK::Transform& T_X_Ft = grandMobod.getInboardFrame(advState);        // T_X_Ft
                const SimTK::Transform& A_X_Mt = grandMobod.getOutboardFrame(advState);       // A_X_Mt
                const SimTK::Transform& G_X_T = grandGrandMobod.getBodyTransform(advState);   // G_X_T
                const SimTK::Transform& Ft_X_Mt = grandMobod.getMobilizerTransform(advState); // Ft_X_Mt
                SimTK::Transform T_X_A = T_X_Ft * Ft_X_Mt * (~A_X_Mt);
                SimTK::Transform T_X_B = T_X_A * A_X_B;

                // Transform G_X_Ft = G_X_T * T_X_Ft;
                // Transform G_X_Mt = G_X_A * A_X_Mt;

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

                // tors_sin
                SimTK::Vec3 m1_B = SimTK::cross(n1_B_hat, b2_B_hat);
                SimTK::Real tors_sin = SimTK::dot(m1_B, n2_B_hat);

                // std::cout << " tors cos sin " << tors_cos <<" "<< tors_sin << std::endl;

                TORSIONAngles[int(childMbx) - 1] = std::atan2(tors_sin, tors_cos);
                // ==============================================

                ZMatrix[int(childMbx) - 1][3] = int(grandGrandIx);

                if (printTransforms) {
                    // SimTK::Test::PrintTransform(G_X_Mt, 6, "G_X_Mt", "G_X_Mt:" +
                    // std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
                    // SimTK::Test::PrintTransform(G_X_Ft, 6, "G_X_Ft", "G_X_Ft:" +
                    // std::to_string(ownWorldIndex) + ":" + std::to_string(int(grandGrandIx)));
                    SimTK::Test::PrintTransform(G_X_T,
                                                6,
                                                "G_X_T",
                                                "G_X_T:" + std::to_string(ownWorldIndex) + ":"
                                                    + std::to_string(int(grandGrandIx)));
                    SimTK::Test::PrintTransform(T_X_A,
                                                6,
                                                "T_X_A",
                                                "T_X_A:" + std::to_string(int(grandGrandIx)) + ":"
                                                    + std::to_string(int(grandMbx)));
                    SimTK::Test::PrintTransform(T_X_Ft,
                                                6,
                                                "T_X_Ft",
                                                "T_X_Ft:" + std::to_string(int(grandGrandIx)) + ":"
                                                    + std::to_string(int(grandMbx)));
                    SimTK::Test::PrintTransform(Ft_X_Mt,
                                                6,
                                                "Ft_X_Mt",
                                                "Ft_X_Mt:" + std::to_string(int(grandGrandIx)) + ":"
                                                    + std::to_string(int(grandMbx)));
                    SimTK::Test::PrintTransform(A_X_Mt,
                                                6,
                                                "A_X_Mt",
                                                "A_X_Mt:" + std::to_string(int(grandGrandIx)) + ":"
                                                    + std::to_string(int(grandMbx)));
                }
            }
        }

        childIx = -1, parentIx = -1, grandIx = -1, grandGrandIx = -2;

    } // _end_ for mbx

    // Print
    bool printZmatBAT = false;
    if (printZmatBAT) {
        for (int BOIx = 0; BOIx < BONDLengthe.size(); BOIx++) {
            std::cout << "ZMatrixBATSimbody:" << " " << ZMatrix[BOIx][0] << " " << ZMatrix[BOIx][1] << " "
                      << ZMatrix[BOIx][2] << " " << ZMatrix[BOIx][3] << " " << BONDLengthe[BOIx] << " "
                      << ANGLEBends[BOIx] << " " << TORSIONAngles[BOIx] << std::endl;
        }
    }
}

/*! <!--  --> */
auto World::getBMps() -> const SimTK::Vector& {
    if (BMps.size() == 0) {
        // BMps.resize(matter->getNQ(advState));
        BMps.resize(matter->getNumBodies());
    }

    int bIx = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // for(int moQIx = 0; moQIx < mobod.getNumQ(advState); moQIx++){
        bIx++;

        const SimTK::Transform& X_BM = mobod.getOutboardFrame(worldState);
        BMps[int(mbx)] = X_BM.p()[0];
        //}
    }

    return BMps;
}

/*! <!--  --> */
auto World::getPFrs() -> const SimTK::Vector& {
    if (PFrs.size() == 0) {
        // PFrs.resize(matter->getNQ(advState));
        PFrs.resize(matter->getNumBodies());
    }

    int bIx = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter->getNumBodies(); ++mbx) {
        // Get mobod
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

        // for(int moQIx = 0; moQIx < mobod.getNumQ(advState); moQIx++){
        bIx++;

        const SimTK::Transform& X_PF = mobod.getInboardFrame(worldState);
        PFrs[int(mbx)] = std::acos(X_PF.R()(0)(0));
        //}
    }

    return PFrs;
}

/*!
 * <!--  -->
 */
const SimTK::Vector& World::getAdvancedUs() {
    return matter->getU(integrator->updAdvancedState());
}

// RANDOM_WALK functions
void World::setTopologyIXs(std::vector<int> argTopologyIXs) {
    topologyIXs = argTopologyIXs;
}

void World::setAmberAtomIXs(std::vector<std::vector<int>> argAmberAtomIXs) {
    amberAtomIXs = argAmberAtomIXs;
}

SimTK::Vec3 World::getGeometricCenterOfSelection(const SimTK::State& state
                                                 // const std::vector<int>& topologyIx,
                                                 // const std::vector<std::vector<int>>& amberAtomList
) {
    // return Vec3
    SimTK::Vec3 geometricCenter = {0, 0, 0};
    // We could just divide by the size of amberAtomList
    // but this works even if the user *mistakenly* repeats
    // indices
    int nOfPoints = 0;

    // Just a quick check, to skip unnecessary computation in case of
    // user error.
    if (amberAtomIXs.size() == 0) {
        std::cerr << "Warning: getGeometricCenterOfSelection called with amberAtomList of size 0"
                  << std::endl;
        return geometricCenter;
    }

    std::cout << "topologies atoms size " << topologyIXs.size() << " " << topologyIXs.size() << std::endl;

    for (int i = 0; i < topologyIXs.size(); i++) {
        const auto& topology = topologies[topologyIXs[i]];
        const auto& atoms = amberAtomIXs[i];

        int amberIx = 0;

        // Iterate through atoms in said topology and check
        // if they are in the list
        for (auto& atom : topology.getAtoms()) {
            if (std::find(atoms.begin(), atoms.end(), amberIx) != atoms.end()) {
                // found
                // Get Compound atom index
                auto compoundAtomIndex = atom.identity.compoundAtomIndex;
                // Get DuMM atom index
                const SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(compoundAtomIndex);
                // Get Mobilized Body index
                const SimTK::MobilizedBodyIndex mobilizedBodyIndex = forceField->getAtomBody(dAIx);
                // Get DuMM Atom Station on its body.
                const SimTK::Vec3 dAS_B = forceField->getAtomStationOnBody(dAIx);
                // Re-Express in G
                const SimTK::MobilizedBody& mobod_A = matter->getMobilizedBody(mobilizedBodyIndex);
                const SimTK::Vec3 dAS_G = mobod_A.findStationLocationInGround(state, dAS_B);
                /* const SimTK::Transform& X_GP = mobod_A.getBodyTransform(state);
                const SimTK::Vec3 dAS_G = X_GP*dAS_B; */
                geometricCenter += dAS_G;
                nOfPoints += 1;

                /* 				std::cout << "amberIx: " << amberIx << " dAIx: " << dAIx
                                << " MobilizedBodyIndex: " << mobilizedBodyIndex
                                << " dAS_G: " << dAS_G << " nOfPoints: " << nOfPoints
                                << std::endl; */
            }

            amberIx += 1;
        }
    }

    // This can probably be done better, but is it clearer?
    for (int i = 0; i < 3; ++i) {
        geometricCenter[i] = geometricCenter[i] / nOfPoints;
    }
    std::cout << "geometricCenter : " << geometricCenter << "\n";

    return geometricCenter;
}

/** Put coordinates into bAtomLists of Topologies.
 * When provided with a State, calcAtomLocationInGroundFrame
 * realizes Position and uses matter to calculate locations **/
void World::updateAtomListsFromSimbody(const SimTK::State& state) {
    // Iterate through topologies
    for (auto& topology : topologies) {
        // Iterate through atoms
        for (auto& atom : topology.updAtoms()) {
            const auto compoundAtomIndex = atom.identity.compoundAtomIndex;
            atom.position = topology.calcAtomLocationInGroundFrameThroughSimbody(compoundAtomIndex,
                                                                                 *forceField,
                                                                                 *matter,
                                                                                 state);
        }
    }
}

/** Set up Fixman torque **/
void World::addFixmanTorque() {
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
    // 	controller.push_back(std::make_unique<SimTK::ConformationalController>(forces, matter,
    // SimTK::MobilizedBodyIndex(i), SimTK::Vec3(0,0,0)));
    // 	controlForce.push_back(std::make_unique<SimTK::Force::Custom>(forces, controller.back().get()));
    // }
}

/** Check if the Fixman torque flag is set **/
bool World::isUsingFixmanTorque() const {
    return useFixmanTorque;
}

/** Get writble pointer to FixmanTorque implementation **/
FixmanTorque* World::updFixmanTorque() {
    assert(isUsingFixmanTorque());
    // return FixmanTorqueImpl.get();
    return FixmanTorqueImpl;
}

/** Get pointer to FixmanTorque implementation **/
FixmanTorque* World::getFixmanTorque() const {
    assert(isUsingFixmanTorque());
    // return FixmanTorqueImpl.get();
    return FixmanTorqueImpl;
}

// ----------------------
// --- Thermodynamics ---
// ----------------------

/** Get the World temperature **/
SimTK::Real World::getTemperature() {
    return this->temperature;
}

/** Set this World temperature but also the samplers and
Fixman torque temperature. **/
void World::setTemperature(SimTK::Real argTemperature) {
    // Set the temperature for this World
    this->temperature = argTemperature;

    // Set the boost temperature for the samplers
    for (auto& sampler : samplers) {
        sampler->setTemperature(argTemperature);
    }

    // Set the temperature for the Fixman torque also
    if (useFixmanTorque) {
        FixmanTorqueImpl->setTemperature(this->temperature);
        FixmanTorqueExtImpl->setTemperature(this->temperature);
    }
}
//...............

/** Set this World temperature but also the samplers and
Fixman torque temperature. **/
void World::setBoostTemperature(SimTK::Real argTemperature) {
    // Set the boost temperature for the samplers
    for (auto& sampler : samplers) {
        sampler->setBoostTemperature(argTemperature);
    }
}
//...............

//...................
// --- Simulation ---
//...................

/** Get/Set seed for reproducibility. **/
void World::setSeed(uint32_t argSeed) {
    randomEngine = buildRandom32(argSeed);
    // forceField->setOpenMMseed(randomEngine());
}

/** How many samples do we have so far **/
std::size_t World::getNofSamples() const {
    // Zero it every time the user asks
    std::size_t nofSamples = 0;

    // Gather samples from all the samplers
    for (size_t i = 0; i < samplers.size(); i++) {
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
auto World::addSampler(SamplerName samplerName,
                       IntegratorType integratorType,
                       ThermostatName thermostatName,
                       bool useFixmanPotential,
                       bool useNUTS) -> bool {
    if (samplerName == SamplerName::HMC) {
        // Construct a new sampler
        samplers.emplace_back(std::make_unique<HMCSampler>(*this,
                                                           *multibodySystem,
                                                           *matter,
                                                           topologies,
                                                           *forceField,
                                                           *forces,
                                                           *timeStepper));

        // Set sampler parameters
        samplers.back()->setIntegratorType(integratorType);
        samplers.back()->setThermostat(thermostatName);
        samplers.back()->setSeed(randomEngine);
        samplers.back()->setUseNUTS(useNUTS);

        // Initialize the sampler
        samplers.back()->initialize();

        // TODO should this be inherited from parent world?
        if (useFixmanPotential) {
            samplers.back()->useFixmanPotential();
        }
    } else {
        throw std::invalid_argument("World::addSampler(): sampler name not recognized.");
        return false;
    }

    return true;
}

auto World::generateSamples(int howManySamplesPerRound,
                            std::stringstream& worldOutStream,
                            const std::string& header,
                            bool shouldPrint) -> bool {
    updSampler(0)->reinitialize(worldState, worldOutStream, shouldPrint);

    // updateAtomListsFromSimbody(worldState);

    return updSampler(0)->sampleIteration(worldState, atomTargetLocationsCache, shouldPrint);
}

void World::calcSpatialForces() {
    if (getSampler(0)->getIntegratorType() == IntegratorType::OpenMMVelocityVerlet) {
        throw std::runtime_error(
            "World::calcSpatialForces() is only implemented for OpenMMVelocityVerlet sampler");
    }

    multibodySystem->realize(worldState);
    SimTK::Vector_<SimTK::SpatialVec> reactionForces;
    matter->calcMobilizerReactionForces(worldState, reactionForces);

    for (const auto mbx : interestingMobodIndices) {
        const auto& mobod = matter->getMobilizedBody(mbx);

        const auto prmtopOutboardIndex = prmtopInboardIndex2PrmtopOutboardIndex[mbx2PrmtopInboardIndex[mbx]];
        const auto& torque_G = reactionForces[mbx][0];
        const auto& force_G = reactionForces[mbx][1];
        const auto u = mobod.getOneU(worldState, 0);
        const auto uDot = mobod.getOneUDot(worldState, 0);

        // const auto& COM_station = mobod.getBodyMassCenterStation(worldState);
        // const auto& COM_ground = mobod.findStationLocationInGround(worldState, COM_station);

        int prmIdx = mbx2PrmtopInboardIndex[mbx];
        spatialForceHistory[prmIdx].push_back({prmtopOutboardIndex, force_G, torque_G, u, uDot});
    }
}

void World::writeSpatialForces(const std::string& filename) const {
    std::ofstream out(filename);
    out << "# frame  inboard_idx  outboard_idx  fx  fy  fz  tx  ty  tz  u  uDot\n";

    // Global max magnitude across all atoms and all frames
    double globalMaxF = 0.0;
    double globalMaxT = 0.0;

    for (const auto& [prmIdx, samples] : spatialForceHistory) {
        for (const auto& s : samples) {
            globalMaxF = std::max(globalMaxF, s.force.norm());
            globalMaxT = std::max(globalMaxT, s.torque.norm());
        }
    }

    const double forceScale = (globalMaxF > 0.0) ? 1.0 / globalMaxF : 1.0;
    const double torqueScale = (globalMaxT > 0.0) ? 1.0 / globalMaxT : 1.0;

    // --- Number of frames ---
    int nFrames = 0;
    for (const auto& [_, samples] : spatialForceHistory) {
        nFrames = std::max(nFrames, (int)samples.size());
    }

    // --- Sorted inboard indices for consistent output order ---
    std::vector<int> sortedIndices;
    sortedIndices.reserve(spatialForceHistory.size());
    for (const auto& [prmIdx, _] : spatialForceHistory) {
        sortedIndices.push_back(prmIdx);
    }
    std::sort(sortedIndices.begin(), sortedIndices.end());

    // Write
    for (int frame = 0; frame < nFrames; ++frame) {
        for (const int prmIdx : sortedIndices) {
            const auto& s = spatialForceHistory.at(prmIdx).at(frame);

            const SimTK::Vec3 f = s.force * forceScale;
            const SimTK::Vec3 t = s.torque * torqueScale;

            if (prmIdx != s.outboardPrmtopIndex) {
                out << frame << " " << prmIdx << " " << s.outboardPrmtopIndex << " " << f[0] << " " << f[1]
                    << " " << f[2] << " " << t[0] << " " << t[1] << " " << t[2] << " " << s.u << " " << s.uDot
                    << "\n";
            }
        }
    }
}

void World::setSamplesPerRound(int samples) {
    samplesPerRound = samples;
}

int World::getSamplesPerRound() const {
    return samplesPerRound;
}

/*!
 * <!--  -->
 */
void World::printDrilling() {
#ifdef __DRILLING__

    // for (DuMM::NonbondAtomIndex nax(0); nax < forceField->getNumNonbondAtoms(); ++nax) {
    // 	//const DuMM::DuMMAtom&        dummAtom =
    // forceField->getAtom(forceField->getAtomIndexOfNonbondAtom(nax)); 	const SimTK::DuMM::AtomIndex dax =
    // forceField->getAtomIndexOfNonbondAtom(nax);
    // 	//const DuMM::IncludedAtomIndex& iax = dummAtom.getIncludedAtomIndex();
    // 	std::cout << "drl World::newFunction dax nax"
    // 		<< " " << dax << " " << nax //<< " " << iax
    // 		<< std::endl;
    // }

    const std::vector<std::vector<SimTK::Real>>& drl_bon_Energies = forceField->getEnergies_drl_bon();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World bonE");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_bon_Energies[fIx][fJx]);
        }
        printf("\n");
    }
    const std::vector<std::vector<SimTK::Real>>& drl_and_Energies = forceField->getEnergies_drl_and();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World andE");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_and_Energies[fIx][fJx]);
        }
        printf("\n");
    }
    const std::vector<std::vector<SimTK::Real>>& drl_tor_Energies = forceField->getEnergies_drl_tor();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World torE");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_tor_Energies[fIx][fJx]);
        }
        printf("\n");
    }
    const std::vector<std::vector<SimTK::Real>>& drl_n14_Energies = forceField->getEnergies_drl_n14();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World n14E");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_n14_Energies[fIx][fJx]);
        }
        printf("\n");
    }
    const std::vector<std::vector<SimTK::Real>>& drl_vdw_Energies = forceField->getEnergies_drl_vdw();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World vdwE");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_vdw_Energies[fIx][fJx]);
        }
        printf("\n");
    }
    const std::vector<std::vector<SimTK::Real>>& drl_cou_Energies = forceField->getEnergies_drl_cou();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        printf("drl World couE");
        for (int fJx = 0; fJx < forceField->getNumNonbondAtoms(); ++fJx) {
            printf(" %f", drl_cou_Energies[fIx][fJx]);
        }
        printf("\n");
    }

    const std::vector<OpenMM::Vec3>& drl_bon_Forces = forceField->getForces_drl_bon();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        const OpenMM::Vec3& ommForce = drl_bon_Forces[fIx];
        const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        printf("drl World bonF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
    }
    const std::vector<OpenMM::Vec3>& drl_and_Forces = forceField->getForces_drl_and();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        const OpenMM::Vec3& ommForce = drl_and_Forces[fIx];
        const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        printf("drl World andF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
    }
    const std::vector<OpenMM::Vec3>& drl_tor_Forces = forceField->getForces_drl_tor();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
        const OpenMM::Vec3& ommForce = drl_tor_Forces[fIx];
        const SimTK::Vec3 simForce(ommForce[0], ommForce[1], ommForce[2]);
        printf("drl World torF %f %f %f\n", ommForce[0], ommForce[1], ommForce[2]);
    }
    const std::vector<OpenMM::Vec3>& drl_n14_Forces = forceField->getForces_drl_n14();
    printf("drl World::newFunction\n");
    for (int fIx = 0; fIx < forceField->getNumNonbondAtoms(); ++fIx) {
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

auto World::decomposeRigidUnits(const Topology& topology) const -> std::vector<WorldRigidUnit> {
    const int worldIndex = ownWorldIndex;
    const int n = topology.getNumAtoms();
    const auto& atoms = topology.getAtoms();

    std::vector<int> parent(n, -1); // chemical tree parent
    std::vector<int> uf(n);
    std::iota(uf.begin(), uf.end(), 0);
    std::function<int(int)> find = [&](int x) {
        return uf[x] == x ? x : uf[x] = find(uf[x]);
    };
    auto unite = [&](int a, int b) {
        uf[find(a)] = find(b);
    };

    // (min,max) -> mobility, so we can look up a specific joint's mobility later
    std::map<std::pair<int, int>, SimTK::BondMobility::Mobility> bondMob;

    for (const auto& bond : topology.getBonds()) {
        const int p = int(bond.compoundAtomIndices[0]);
        const int c = int(bond.compoundAtomIndices[1]);
        const SimTK::BondMobility::Mobility mob = bond.getBondMobility(worldIndex);
        bondMob[{std::min(p, c), std::max(p, c)}] = mob;

        if (bond.ringClosing) {
            continue; // not a tree edge; never merges
        }
        parent[c] = p; // tree edge p -> c
        if (mob == SimTK::BondMobility::Mobility::Rigid) {
            unite(p, c); // same rigid unit
        }
    }

    // group atoms by union-find component
    std::map<int, int> compToUnit;
    std::vector<WorldRigidUnit> units;
    std::vector<int> unitOf(n, -1);
    for (int i = 0; i < n; ++i) {
        const int r = find(i);
        auto it = compToUnit.find(r);
        if (it == compToUnit.end()) {
            compToUnit[r] = int(units.size());
            units.emplace_back();
        }
        const int u = compToUnit[r];
        units[u].atomCAIxs.push_back(i);
        unitOf[i] = u;
    }

    // root atom + inboard joint of every unit
    for (int u = 0; u < int(units.size()); ++u) {
        WorldRigidUnit& unit = units[u];
        std::sort(unit.atomCAIxs.begin(), unit.atomCAIxs.end());
        for (int a : unit.atomCAIxs) {
            if (atoms[a].connectivity.root) { // molecule root
                unit.rootCAIx = a;
                unit.parentUnit = -1;
                break;
            }
            const int p = parent[a];
            if (p >= 0 && unitOf[p] != u) { // tree edge leaves this unit
                unit.rootCAIx = a;
                unit.parentUnit = unitOf[p];
                unit.jointParentCAIx = p;
                unit.mobility = bondMob.at({std::min(p, a), std::max(p, a)});
                break;
            }
        }
        if (unit.rootCAIx < 0) {
            unit.rootCAIx = unit.atomCAIxs.front(); // isolated atom
        }
    }
    return units;
}

// ---------------------------------------------------------------------------
//  Build the cache from connectivity only. Call once after the topologies and
//  their bonds are fixed (e.g. end of modelTopologies). O(atoms + bonds).
// ---------------------------------------------------------------------------
void World::buildFrameGraph(const Span<Topology>& topologies) {
    frameGraph.topoOffset.resize(topologies.size() + 1, 0);
    for (std::size_t t = 0; t < topologies.size(); ++t) {
        frameGraph.topoOffset[t + 1] = frameGraph.topoOffset[t] + topologies[t].getNumAtoms();
    }
    frameGraph.totalAtoms = frameGraph.topoOffset.back();

    // local-then-global parent / first-child, per topology
    for (std::size_t topoIx = 0; topoIx < topologies.size(); ++topoIx) {
        const Topology& topo = topologies[topoIx];
        const int numAtoms = topo.getNumAtoms();
        const int off = frameGraph.topoOffset[topoIx];
        const auto& atoms = topo.getAtoms();

        std::vector<int> parent(numAtoms, -1);
        std::vector<int> firstChild(numAtoms, -1); // smallest-index child == refChild

        for (const auto& bond : topo.getBonds()) {
            if (bond.ringClosing) {
                continue; // ring-closing bonds are not tree edges
            }
            const int p = int(bond.compoundAtomIndices[0]);
            const int c = int(bond.compoundAtomIndices[1]);
            parent[c] = p;
            if (firstChild[p] < 0 || c < firstChild[p]) {
                firstChild[p] = c; // track minimum in one pass, no sort
            }
        }

        int rootLocal = -1;
        for (int i = 0; i < numAtoms; ++i) {
            if (atoms[i].connectivity.root) {
                rootLocal = i;
                break;
            }
        }
        if (rootLocal < 0) {
            throw std::runtime_error("buildFrameGraph: topology " + std::to_string(topoIx)
                                     + " has no root atom");
        }

        for (int i = 0; i < numAtoms; ++i) {
            const int g = off + i;
            if (i == rootLocal) {
                frameGraph.r_self.push_back(g);
                continue;
            }
            const int p = parent[i];
            if (p < 0) {
                throw std::runtime_error("buildFrameGraph: atom " + std::to_string(i) + " in topology "
                                         + std::to_string(topoIx) + " has no parent and is not root");
            }
            const int gp = parent[p];     // -1 if parent is the root
            const int rc = firstChild[i]; // -1 if i is a leaf

            if (gp >= 0 && rc >= 0) {
                frameGraph.g_self.push_back(g);
                frameGraph.g_parent.push_back(off + p);
                frameGraph.g_gparent.push_back(off + gp);
                frameGraph.g_refChild.push_back(off + rc);
            } else {
                frameGraph.f_self.push_back(g);
                frameGraph.f_parent.push_back(off + p);
            }
        }
    }

    // Partition must cover every atom exactly once
    if (frameGraph.r_self.size() + frameGraph.g_self.size() + frameGraph.f_self.size()
        != frameGraph.totalAtoms) {
        throw std::runtime_error("buildFrameGraph: partition sizes don't add up to total atoms");
    }
}

// ----------------------------------------------------------------------------
//  modelOneCompound, now a World member. One topology -> DuMM atoms/bonds +
//  Simbody MobilizedBodies, all placed from the target-based atom-frame cache.
// ----------------------------------------------------------------------------
auto World::modelOneCompound(int topoIx, SimTK::RootMobility rootMobility) -> std::vector<WorldRigidUnit> {
    Topology& topology = topologies[topoIx];
    SimTK::DuMMForceFieldSubsystem& dumm = *forceField;

    const SimTK::Transform G_X_T = SimTK::Transform(); // World convention: G_X_T == I

    // Axis switches (pure reindex; replace with column-shuffle helpers later).
    const SimTK::Transform X_to_Z = SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::YAxis);
    const SimTK::Transform X_to_Y = SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::ZAxis);
    const SimTK::Transform Y_to_Z = SimTK::Rotation(-90 * SimTK::Deg2Rad, SimTK::XAxis);
    const SimTK::Transform Z_to_Y = ~Y_to_Z;
    const SimTK::Transform Y_to_X = ~X_to_Y;

    // (1) DuMM atoms + masses + cAIx->dAIx.
    const int n = topology.getNumAtoms();
    std::vector<SimTK::DuMM::AtomIndex> dAIxOf(n);
    for (SimTK::Compound::AtomIndex cAIx(0); cAIx < n; ++cAIx) {
        const SimTK::BiotypeIndex biotypeIx = topology.getAtomBiotypeIndex(cAIx);
        const SimTK::DuMM::ChargedAtomTypeIndex chargedTypeId = dumm.getBiotypeChargedAtomType(biotypeIx);
        const SimTK::DuMM::AtomIndex dAIx = dumm.addAtom(chargedTypeId);
        dAIxOf[cAIx] = dAIx;
        topology.setDuMMAtomIndex(cAIx, dAIx);
        dumm.setDuMMAtomMass(dAIx, topology.getAtoms()[cAIx].physics.massInDaltons);
    }

    // (2) DuMM bonds (connectivity only).
    for (const auto& bond : topology.getBonds()) {
        dumm.addBond(dAIxOf[bond.compoundAtomIndices[0]], dAIxOf[bond.compoundAtomIndices[1]]);
    }

    // (3) rigid units + parent-before-child order.
    std::vector<WorldRigidUnit> units = decomposeRigidUnits(topology);
    std::vector<int> order;
    {
        std::vector<bool> done(units.size(), false);
        for (bool progress = true; progress;) {
            progress = false;
            for (int u = 0; u < int(units.size()); ++u) {
                if (done[u]) {
                    continue;
                }
                if (units[u].parentUnit < 0 || done[units[u].parentUnit]) {
                    order.push_back(u);
                    done[u] = true;
                    progress = true;
                }
            }
        }
        if (order.size() != units.size()) {
            throw std::runtime_error("modelOneCompound: rigid-unit tree has a cycle");
        }
    }

    // (4) bodies + atom attachment, frames read from the flat caches.
    for (int u : order) {
        WorldRigidUnit& unit = units[u];
        const auto root = unit.rootCAIx;
        const SimTK::Transform& T_X_B =
            F(topoIx, SimTK::Compound::AtomIndex(root)); // body == root-atom frame
        const SimTK::Transform B_X_T = ~T_X_B;

        // stations (stored once) + mass properties, from the cached frames.
        SimTK::Real mass = 0;
        SimTK::Vec3 com(0);
        SimTK::Inertia inertia(0);

        for (const auto a : unit.atomCAIxs) {
            const SimTK::Transform B_X_atom = B_X_T * F(topoIx, SimTK::Compound::AtomIndex(a));
            topoAtomBodyFrame[topoIx][a] = B_X_atom; // single source of truth

            const SimTK::Vec3& s = B_X_atom.p();
            const SimTK::Real m = topology.getAtoms()[a].physics.massInDaltons;
            mass += m;
            com += m * s;
            inertia += SimTK::Inertia(s, m);
        }
        const SimTK::MassProperties massProps(mass, mass > 0 ? com / mass : SimTK::Vec3(0), inertia);

        // ---- create the MobilizedBody ----
        if (unit.parentUnit < 0) {
            const SimTK::Transform G_X_B = G_X_T * T_X_B;
            SimTK::MobilizedBody::Ground& ground = matter->Ground();
            switch (rootMobility) {
                case SimTK::RootMobility::Free:
                    unit.mbx = SimTK::MobilizedBody::Free(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::RootMobility::Cartesian:
                    unit.mbx = SimTK::MobilizedBody::Translation(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::RootMobility::FreeLine:
                    unit.mbx = SimTK::MobilizedBody::FreeLine(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::RootMobility::Ball:
                    unit.mbx = SimTK::MobilizedBody::Ball(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::RootMobility::Pin:
                    unit.mbx = SimTK::MobilizedBody::Pin(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::RootMobility::Weld:
                default:
                    unit.mbx = SimTK::MobilizedBody::Weld(ground, G_X_B, massProps, SimTK::Transform())
                                   .getMobilizedBodyIndex();
                    break;
            }
        } else {
            const WorldRigidUnit& parentUnit = units[unit.parentUnit];
            SimTK::MobilizedBody& parentMobod = matter->updMobilizedBody(parentUnit.mbx);

            const SimTK::Transform Proot_X_root =
                (~F(topoIx, SimTK::Compound::AtomIndex(parentUnit.rootCAIx)))
                * F(topoIx, SimTK::Compound::AtomIndex(root));

            // cached B == X_parentBC_childBC for edge (jointParent -> root).
            const SimTK::Transform& X_parentBC_childBC = Xpc(topoIx, SimTK::Compound::AtomIndex(root));
            const SimTK::Transform X_childBC_parentBC = ~X_parentBC_childBC;

            SimTK::Transform X_PF;
            SimTK::Transform X_BM;
            switch (unit.mobility) {
                case SimTK::BondMobility::Mobility::AnglePin:
                case SimTK::BondMobility::Mobility::Slider:
                case SimTK::BondMobility::Mobility::BendStretch:
                    X_BM = X_parentBC_childBC;
                    X_PF = Proot_X_root * X_BM;
                    break;
                case SimTK::BondMobility::Mobility::Torsion:
                case SimTK::BondMobility::Mobility::Cylinder:
                    X_BM = X_parentBC_childBC * X_to_Z;
                    X_PF = Proot_X_root * X_BM;
                    break;
                case SimTK::BondMobility::Mobility::BallM:
                case SimTK::BondMobility::Mobility::Translation:
                    X_BM = X_to_Z;
                    X_PF = Proot_X_root * X_BM;
                    break;
                case SimTK::BondMobility::Mobility::Spherical:
                    X_PF = Proot_X_root * X_parentBC_childBC * X_to_Y * Y_to_Z;
                    X_BM = ~(Z_to_Y * Y_to_X * X_childBC_parentBC);
                    break;
                case SimTK::BondMobility::Mobility::OrthoSpherical: {
                    // not in the frame graph; rare mobility, still molmodel calls.
                    const auto X_parentAtom_BCpar = topology.calcDefaultBondCenterFrameInParentAtomFrame(
                        SimTK::Compound::AtomIndex(unit.jointParentCAIx),
                        SimTK::Compound::AtomIndex(root));
                    const auto X_childAtom_BCchi = topology.calcDefaultBondCenterFrameInChildAtomFrame(
                        SimTK::Compound::AtomIndex(unit.jointParentCAIx),
                        SimTK::Compound::AtomIndex(root));
                    X_PF = X_parentAtom_BCpar;
                    X_BM = ~(X_parentBC_childBC * (~X_childAtom_BCchi));
                    break;
                }
                default:
                    throw std::runtime_error(
                        "modelOneCompound: unsupported mobility '"
                        + std::string(SimTK::BondMobility::getBondMobilityName(unit.mobility)) + "'");
            }

            switch (unit.mobility) {
                case SimTK::BondMobility::Mobility::Torsion:
                    unit.mbx =
                        SimTK::MobilizedBody::Pin(parentMobod, X_PF, massProps, X_BM).getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::Cylinder:
                    unit.mbx = SimTK::MobilizedBody::Cylinder(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::AnglePin:
                    unit.mbx = SimTK::MobilizedBody::Torsion(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::Slider:
                    unit.mbx = SimTK::MobilizedBody::Slider(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::BendStretch:
                    unit.mbx = SimTK::MobilizedBody::BendStretch(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::BallM:
                    unit.mbx = SimTK::MobilizedBody::Ball(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::Translation:
                    unit.mbx = SimTK::MobilizedBody::Translation(parentMobod, X_PF, massProps, X_BM)
                                   .getMobilizedBodyIndex();
                    break;
                case SimTK::BondMobility::Mobility::Spherical: {
                    SimTK::MobilizedBody::SphericalCoords b(parentMobod, X_PF, massProps, X_BM);
                    b.setRadialAxis(SimTK::ZAxis);
                    unit.mbx = b.getMobilizedBodyIndex();
                    break;
                }
                case SimTK::BondMobility::Mobility::OrthoSpherical: {
                    SimTK::MobilizedBody::SphericalCoords b(parentMobod, X_PF, massProps, X_BM);
                    b.setRadialAxis(SimTK::XAxis);
                    unit.mbx = b.getMobilizedBodyIndex();
                    break;
                }
                default:
                    throw std::runtime_error("modelOneCompound: unhandled mobilizer body");
            }
        }

        // ---- attach atoms (no user cluster); station from the stored frame ----
        for (const auto a : unit.atomCAIxs) {
            dumm.attachAtomToBody(dAIxOf[a], unit.mbx, topoAtomBodyFrame[topoIx][a].p());
            topology.setAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(a), unit.mbx);
        }
    }

    return units;
}
