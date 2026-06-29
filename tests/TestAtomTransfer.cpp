// ============================================================================
//  TestAtomTransfer.cpp -- the ATOM/TRANSFER layer the engine suite never
//  exercised. RobotBuilders.hpp deliberately builds models with numAtoms == 0
//  ("no atoms ... those belong to the geometry/transfer layer"), so every test
//  that ships today validates the multibody MECHANICS on a bare body tree and
//  leaves the boundary that the original molmodel/DuMM "problematic calls" were
//  all about completely untested:
//
//    original molmodel/DuMM call                  ported responsibility (current)
//    ------------------------------------------   --------------------------------
//    DuMM::bsetAtomStationOnBody / station_B      RobotModel::atomStation_B
//    CompoundAtom::setFrameInMobilizedBodyFrame   atomStation_B (+ body frame)
//    getAtomStationOnBody / getLocationInMBF      (~X_GB[b]) * atomPosG  (invariant)
//    calcAtomLocationInGroundFrameThroughSimbody  RobotEngine::realizePosition tail
//      / DuMM::realizeSubsystemPositionImpl         + fillAtomPositionsFromBodies
//    updIncludedAtomStation/getIncludedAtomStn    RobotState::atomStationG (R_GB*station)
//    calcAtomVelocityInGroundFrame                V_GB rigid transport to the atom
//    calcAtomAccelerationInGroundFrame            A_GB + V_GB rigid point accel
//    DuMM::calcMassProperties ("we do it ourselves")  bodyMass/Com_B/UnitInertia_B
//                                                     -> Mk_G in realizePosition
//    HMCSampler::sampleIteration accept-path       transfer round trip
//      target locations <-> internal state
//
//  These are all engine/model-level facts (no World, no OpenMM, no frame graph),
//  so they are tested here the same way the rest of the suite tests mechanics:
//  hand-built RobotModels, deterministic RNG, finite-difference cross-checks.
//  The CONSTRAINTS chosen below are the ones the original code silently relied
//  on; see each test's banner for the precise invariant and its tolerance.
// ============================================================================
#include <array>
#include <gtest/gtest.h>
#include <vector>

#include "RobotBuilders.hpp"
#include "RobotEngine.hpp"
#include "RobotModel.hpp"
#include "RobotState.hpp"
#include "TestHelpers.hpp"

using namespace robo;
using rtest::BodySpec;
using rtest::buildForest;
using rtest::randomizeState;
using rtest::Rng;

namespace {

// ---------------------------------------------------------------------------
//  Attach `perBody` random atoms to every non-Ground body of an already-built
//  model, then DERIVE the body mass properties from those atom point masses --
//  exactly as World::recomputeGeometry does (bodyMass = sum m, com = mass-
//  weighted station mean, unit inertia = sum Inertia(station,m) / mass). This
//  makes the model's mass properties and its atoms self-consistent, which is the
//  precondition the engine's Mk_G build assumes.
// ---------------------------------------------------------------------------
void attachAtoms(RobotModel& m, Rng& rng, int perBody = 3) {
    const int B = m.numBodies;
    m.atomBody.clear();
    m.atomStation_B.clear();
    m.atomMass.clear();

    std::vector<std::vector<int>> bodyAtoms(B);
    for (int b = 1; b < B; ++b) {
        for (int k = 0; k < perBody; ++k) {
            const int a = static_cast<int>(m.atomBody.size());
            m.atomBody.push_back(b);
            m.atomStation_B.push_back(rng.vec3(Real(-0.3), Real(0.3)));
            m.atomMass.push_back(rng.uniform(Real(1.0), Real(16.0))); // ~H..O daltons
            bodyAtoms[b].push_back(a);
        }
    }
    m.numAtoms = static_cast<int>(m.atomBody.size());

    // body-sorted atom CSR (optional for the engine, used by our own sums)
    m.bodyAtomsBeg.assign(B, 0);
    m.bodyAtomsEnd.assign(B, 0);
    m.bodyAtoms.clear();
    for (int b = 0; b < B; ++b) {
        m.bodyAtomsBeg[b] = static_cast<int>(m.bodyAtoms.size());
        for (int a : bodyAtoms[b]) {
            m.bodyAtoms.push_back(a);
        }
        m.bodyAtomsEnd[b] = static_cast<int>(m.bodyAtoms.size());
    }

    // mass properties from the atom point masses (the calcMassProperties port).
    m.bodyMass.assign(B, Real(0));
    m.bodyCom_B.assign(B, Vec3(0));
    m.bodyUnitInertia_B.assign(B, UnitInertia(Real(0), Real(0), Real(0)));
    for (int b = 1; b < B; ++b) {
        Real mass = 0;
        Vec3 com(0);
        Inertia inertia(0);
        for (int ci = m.bodyAtomsBeg[b]; ci < m.bodyAtomsEnd[b]; ++ci) {
            const int a = m.bodyAtoms[ci];
            const Vec3& st = m.atomStation_B[a];
            const Real ma = m.atomMass[a];
            mass += ma;
            com = com + (st * ma);
            inertia += Inertia(st, ma);
        }
        const MassProperties mp(mass, mass > 0 ? Vec3(com * (Real(1) / mass)) : Vec3(0), inertia);
        m.bodyMass[b] = mass;
        m.bodyCom_B[b] = mp.getMassCenter();
        m.bodyUnitInertia_B[b] = mp.getUnitInertia();
    }
}

// A forest of several robots on the shared Ground, one of every joint type, so
// the atom transfer is exercised across the full mobilizer spread at once.
RobotModel makeAtomicForest(Rng& rng, int perBody = 3) {
    auto F = [&](int parent, JointType jt) {
        BodySpec s;
        s.parent = parent;
        s.joint = jt;
        s.X_PF = Transform(rng.rotation(), rng.vec3());
        s.X_BM = Transform(rng.rotation(), rng.vec3());
        return s;
    };
    std::vector<BodySpec> specs;
    specs.push_back(F(0, JointType::Free));        // body 1
    specs.push_back(F(1, JointType::Torsion));     // body 2
    specs.push_back(F(2, JointType::BendStretch)); // body 3 (q-dependent H_FM)
    specs.push_back(F(0, JointType::Ball));        // body 4
    specs.push_back(F(0, JointType::Cartesian));   // body 5
    specs.push_back(F(5, JointType::Slider));      // body 6
    specs.push_back(F(0, JointType::FreeLine));    // body 7
    specs.push_back(F(0, JointType::Rigid));       // body 8 (frozen)
    RobotModel m = buildForest(specs);
    attachAtoms(m, rng, perBody);
    return m;
}

// atom velocity in Ground = rigid transport of the body spatial velocity to the
// atom: v_a = v_o + w x (r_a - o). (calcAtomVelocityInGroundFrame.)
Vec3 atomVelG(const RobotState& s, const RobotModel& m, int a) {
    const int b = m.atomBody[a];
    const Vec3 r = s.atomPosG()[a] - s.X_GB()[b].p();
    return s.V_GB()[b][1] + (s.V_GB()[b][0] % r);
}

} // namespace

// ---------------------------------------------------------------------------
//  T1. station_B is BODY-FIXED: for ANY configuration, an atom's Ground position
//      re-expressed in its body frame returns exactly its station. This is the
//      property getAtomStationOnBody / getLocationInMobilizedBodyFrame silently
//      depended on (the station is a topological constant; only X_GB moves).
//      Tolerance kTight (machine algebra): the map is one rigid inverse.
// ---------------------------------------------------------------------------
TEST(AtomTransfer, StationIsBodyFixedUnderAnyConfiguration) {
    Rng rng(0xA70505);
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);
    for (int rep = 0; rep < 30; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = m.atomBody[a];
            const Vec3 stationBack = (~s.X_GB()[b]) * s.atomPosG()[a];
            EXPECT_TRUE(rtest::NearVec3(stationBack, m.atomStation_B[a], rtest::kTight))
                << "atom " << a << " body " << b << " rep " << rep;
        }
    }
}

// ---------------------------------------------------------------------------
//  T2. Ground-position currency is internally consistent: the per-atom positions
//      written by realizePosition equal both (a) the explicit station map
//      X_GB.p + R_GB*station (= fillAtomPositionsFromBodies, the accept-path) and
//      (b) origin + atomStationG, where atomStationG = R_GB*station is the
//      getIncludedAtomStation referent fed to the force bridge. All three must
//      agree to the last bit -- they are the SAME quantity computed two ways, so
//      any drift is a real bug, not FD noise. Tolerance kTight.
// ---------------------------------------------------------------------------
TEST(AtomTransfer, RealizePositionFillAndStationGAllAgree) {
    Rng rng(0xB22222);
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);

    RobotEngine::realizePosition(m, s);
    std::vector<Vec3> posFromRealize(s.atomPosG(), s.atomPosG() + m.numAtoms);
    std::vector<Vec3> stationG(s.atomStationG(), s.atomStationG() + m.numAtoms);

    // scribble over atomPosG, then recompute via the accept-path entry point.
    for (int a = 0; a < m.numAtoms; ++a) {
        s.atomPosG()[a] = Vec3(0);
    }
    RobotEngine::fillAtomPositionsFromBodies(m, s);

    for (int a = 0; a < m.numAtoms; ++a) {
        const int b = m.atomBody[a];
        const Vec3 viaStationMap = s.X_GB()[b].p() + s.X_GB()[b].R() * m.atomStation_B[a];
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], posFromRealize[a], rtest::kTight))
            << "fill vs realize " << a;
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], viaStationMap, rtest::kTight)) << "fill vs map " << a;
        // atomStationG == R_GB * station, and origin + stationG == position.
        EXPECT_TRUE(rtest::NearVec3(stationG[a], s.X_GB()[b].R() * m.atomStation_B[a], rtest::kTight))
            << "stationG " << a;
        EXPECT_TRUE(rtest::NearVec3(s.X_GB()[b].p() + stationG[a], posFromRealize[a], rtest::kTight))
            << "origin+stationG " << a;
    }
}

// ---------------------------------------------------------------------------
//  T3. Atom velocity in Ground (calcAtomVelocityInGroundFrame) is the rigid
//      transport of the body spatial velocity, and must equal d/dt of the atom's
//      Ground position along the true trajectory q(t)=q0+t*qdot. Central FD,
//      tolerance 1e-4 (~ the kinematic FD tolerance used elsewhere in the suite).
// ---------------------------------------------------------------------------
TEST(AtomTransfer, AtomVelocityIsDerivativeOfAtomPosition) {
    Rng rng(0xC33033);
    const Real h = 1e-6;
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 20; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);

        std::vector<Real> q0(s.q(), s.q() + m.nq), qdot(m.nq);
        RobotEngine::calcQDot(m, s, qdot.data());
        std::vector<Vec3> vAnalytic(m.numAtoms);
        for (int a = 0; a < m.numAtoms; ++a) {
            vAnalytic[a] = atomVelG(s, m, a);
        }

        auto posAt = [&](Real t, std::vector<Vec3>& P) {
            for (int i = 0; i < m.nq; ++i) {
                s.q()[i] = q0[i] + t * qdot[i];
            }
            RobotEngine::normalizeQuaternions(m, s);
            RobotEngine::realizePosition(m, s);
            P.assign(s.atomPosG(), s.atomPosG() + m.numAtoms);
        };
        std::vector<Vec3> Pp, Pm;
        posAt(+h, Pp);
        posAt(-h, Pm);
        for (int a = 0; a < m.numAtoms; ++a) {
            const Vec3 vFD = (Pp[a] - Pm[a]) * (Real(1) / (Real(2) * h));
            EXPECT_TRUE(rtest::NearVec3(vAnalytic[a], vFD, 1e-4)) << "atom " << a << " rep " << rep;
        }
        std::copy(q0.begin(), q0.end(), s.q());
    }
}

// ---------------------------------------------------------------------------
//  T4. Atom acceleration in Ground (calcAtomAccelerationInGroundFrame) is the
//      rigid point acceleration a_a = a_o + alpha x r + w x (w x r), and must
//      equal d/dt of the atom velocity along the true (q,u) trajectory. Drives a
//      random generalized force through full forward dynamics first. Central FD
//      of the atom velocity, tolerance 1e-3 (second-kinematic-derivative FD, as
//      in TestMobilizer's A_GB leg).
// ---------------------------------------------------------------------------
TEST(AtomTransfer, AtomAccelerationIsDerivativeOfAtomVelocity) {
    Rng rng(0xD44044);
    const Real h = 1e-6;
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 12; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);
        RobotEngine::realizeVelocity(m, s);
        RobotEngine::realizeArticulatedBodyInertias(m, s);
        std::fill(s.bodyForceG(), s.bodyForceG() + m.numBodies, SpatialVec(Vec3(0), Vec3(0)));
        for (int i = 0; i < m.nu; ++i) {
            s.mobilityForce()[i] = rng.gaussian(0, Real(0.5));
        }
        RobotEngine::calcUDot(m, s);

        std::vector<Real> q0(s.q(), s.q() + m.nq), u0(s.u(), s.u() + m.nu);
        std::vector<Real> qdot0(m.nq), udot0(s.udot(), s.udot() + m.nu);
        RobotEngine::calcQDot(m, s, qdot0.data());

        // analytic atom acceleration from body spatial accel + velocity.
        std::vector<Vec3> aAnalytic(m.numAtoms);
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = m.atomBody[a];
            const Vec3 r = s.atomPosG()[a] - s.X_GB()[b].p();
            const Vec3 w = s.V_GB()[b][0];
            const Vec3 alpha = s.A_GB()[b][0];
            const Vec3 a_o = s.A_GB()[b][1];
            aAnalytic[a] = a_o + (alpha % r) + (w % (w % r));
        }

        auto atomVelAt = [&](Real t, std::vector<Vec3>& V) {
            for (int i = 0; i < m.nq; ++i) {
                s.q()[i] = q0[i] + t * qdot0[i];
            }
            for (int i = 0; i < m.nu; ++i) {
                s.u()[i] = u0[i] + t * udot0[i];
            }
            RobotEngine::normalizeQuaternions(m, s);
            RobotEngine::realizePosition(m, s);
            RobotEngine::realizeVelocity(m, s);
            V.resize(m.numAtoms);
            for (int a = 0; a < m.numAtoms; ++a) {
                V[a] = atomVelG(s, m, a);
            }
        };
        std::vector<Vec3> Vp, Vm;
        atomVelAt(+h, Vp);
        atomVelAt(-h, Vm);
        for (int a = 0; a < m.numAtoms; ++a) {
            const Vec3 aFD = (Vp[a] - Vm[a]) * (Real(1) / (Real(2) * h));
            EXPECT_TRUE(rtest::NearVec3(aAnalytic[a], aFD, 1e-3)) << "atom " << a << " rep " << rep;
        }
        std::copy(q0.begin(), q0.end(), s.q());
        std::copy(u0.begin(), u0.end(), s.u());
    }
}

// ---------------------------------------------------------------------------
//  T5. The calcMassProperties REPLACEMENT. The body spatial inertia the engine
//      builds (Mk_G, about the body origin, in Ground) must equal the spatial
//      inertia assembled INDEPENDENTLY from the atom point masses at their Ground
//      positions. Concretely:
//        * bodyMass[b]            == sum_a m_a                       (kTight)
//        * comG[b]                == (sum_a m_a posG[a]) / M          (kAlg)
//        * Mk_G unit inertia * M  == sum_a Inertia(posG[a]-origin,m_a)(kAlg)
//      This is the whole point of "we handle mass properties ourselves": the
//      reduction from atoms to a rigid body must be exact, in Ground, every
//      transfer. (Validates recomputeGeometry's sums AND realizePosition's
//      reexpression bodyUnitInertia_B -> G_Bo_G in one shot.)
// ---------------------------------------------------------------------------
TEST(AtomTransfer, BodySpatialInertiaEqualsAtomPointMassSum) {
    Rng rng(0xE55055);
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    for (int b = 1; b < m.numBodies; ++b) {
        if (m.bodyAtomsBeg[b] == m.bodyAtomsEnd[b]) {
            continue;
        }
        Real mass = 0;
        Vec3 comG(0);
        Inertia inertiaAboutOriginG(0);
        const Vec3 origin = s.X_GB()[b].p();
        for (int ci = m.bodyAtomsBeg[b]; ci < m.bodyAtomsEnd[b]; ++ci) {
            const int a = m.bodyAtoms[ci];
            const Real ma = m.atomMass[a];
            const Vec3 r = s.atomPosG()[a] - origin; // atom offset from origin, in Ground
            mass += ma;
            comG = comG + (s.atomPosG()[a] * ma);
            inertiaAboutOriginG += Inertia(r, ma);
        }
        comG = comG * (Real(1) / mass);

        // (a) total mass
        EXPECT_NEAR(m.bodyMass[b], mass, rtest::kTight) << "mass body " << b;
        // (b) COM in Ground: engine comG[b] is origin + R_GB*com_B.
        EXPECT_TRUE(rtest::NearVec3(s.comG()[b], comG, rtest::kAlg)) << "comG body " << b;
        // (c) inertia about the body origin, expressed in Ground.
        const Mat33 fromEngine = s.Mk_G()[b].getUnitInertia().full() * m.bodyMass[b];
        const Mat33 fromAtoms = inertiaAboutOriginG.asSymMat33().full();
        EXPECT_TRUE(rtest::NearMat33(fromEngine, fromAtoms, rtest::kAlg)) << "inertia body " << b;
    }
}

// ---------------------------------------------------------------------------
//  T6. TRANSFER ROUND TRIP (HMCSampler::sampleIteration accept-path <-> the
//      getAtomLocationInMobilizedBodyFrameThroughDumm inverse). Given ANY set of
//      per-atom Cartesian "targets", the rigid decomposition recovers them
//      EXACTLY: station' = (~X_GB[b]) * target, then refill posG = X_GB * station'
//      reproduces target to machine precision -- for every body configuration.
//
//      This is the analogue of the original TestAtomTargetLocations residual
//      check, but where molmodel's iterative matchDefault*/getTransformAndResidual
//      only reached residual < 0.02 nm, the SoA decomposition is exact (the
//      station absorbs the within-body offset), so the bar here is kTight (~1e-13)
//      rather than 2e-2. Tightening the constraint is the correctness gain.
// ---------------------------------------------------------------------------
TEST(AtomTransfer, CartesianTargetsRoundTripExactly) {
    Rng rng(0xF66066);
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);

    for (int rep = 0; rep < 15; ++rep) {
        randomizeState(m, s, rng);
        RobotEngine::realizePosition(m, s);

        // arbitrary targets (need not be a rigid placement; the station absorbs
        // whatever offset, which is exactly why a transfer is lossless).
        std::vector<Vec3> targets(m.numAtoms);
        for (int a = 0; a < m.numAtoms; ++a) {
            targets[a] = rng.vec3(Real(-5), Real(5));
        }

        // inverse: derive a station per atom from the target and current frames.
        std::vector<Vec3> savedStation = m.atomStation_B;
        for (int a = 0; a < m.numAtoms; ++a) {
            const int b = m.atomBody[a];
            m.atomStation_B[a] = (~s.X_GB()[b]) * targets[a];
        }
        // forward again: recover the targets.
        RobotEngine::fillAtomPositionsFromBodies(m, s);
        for (int a = 0; a < m.numAtoms; ++a) {
            EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], targets[a], rtest::kTight))
                << "atom " << a << " rep " << rep;
        }
        m.atomStation_B = savedStation;
    }
}

// ---------------------------------------------------------------------------
//  T7. Atoms follow forest indexing: an atom is moved ONLY by its own robot's
//      coordinates. Perturbing robot R1's q must not move any atom whose body
//      belongs to robot R3 -- the atom<->body map and the children CSR carry the
//      molecule="robot" independence all the way down to the per-atom positions.
//      Tolerance kTight.
// ---------------------------------------------------------------------------
TEST(AtomTransfer, AtomsMoveOnlyWithTheirOwnRobot) {
    Rng rng(0x121212);
    RobotModel m = makeAtomicForest(rng);
    RobotState s;
    s.allocateFull(m);
    randomizeState(m, s, rng);
    RobotEngine::realizePosition(m, s);

    // robot rooted at the Cartesian body 5 occupies bodies {5,6}; snapshot its atoms.
    const std::array<int, 2> otherBodies = {5, 6};
    std::vector<std::pair<int, Vec3>> before;
    for (int a = 0; a < m.numAtoms; ++a) {
        for (int b : otherBodies) {
            if (m.atomBody[a] == b) {
                before.emplace_back(a, s.atomPosG()[a]);
            }
        }
    }
    ASSERT_TRUE(!before.empty());

    // perturb robot 1 (bodies {1,2,3}).
    for (int b : {1, 2, 3}) {
        for (int j = 0; j < m.bodyNQ[b]; ++j) {
            s.q()[m.bodyQIndex[b] + j] += Real(0.4);
        }
    }
    RobotEngine::normalizeQuaternions(m, s);
    RobotEngine::realizePosition(m, s);

    for (auto& [a, p] : before) {
        EXPECT_TRUE(rtest::NearVec3(s.atomPosG()[a], p, rtest::kTight)) << "atom " << a << " moved";
    }
}