// ============================================================================
//  TestForceReducer.cpp -- pins the CUDA reduceForces kernel
//  (src/OpenMMContext.cpp) against the shared host reduction
//  reduceAtomForcesToBodies (include/bridge/ForceReducer.hpp), per
//  docs/specs/refactor/SPLIT-DEDUP-FORCEREDUCER.md.
//
//  INV-1 (force->wrench convention): for body b with Ground origin
//  X_GB[b].p(), bodyForceG[b].linear = sum_a f_a and bodyForceG[b].angular =
//  sum_a (r_a - origin_b) x f_a, summed over b's REAL atoms. Both the host
//  reducer and the CUDA reduceForces kernel MUST produce this identically --
//  this test drives the SAME per-atom Ground forces (an actual OpenMM
//  NonbondedForce evaluation, so the input is deterministic and
//  reproducible) through both and asserts the per-body wrenches agree.
//
//  INV-2 (virtual-site skip): the fixture includes a massless ThreeParticle-
//  AverageSite (an M-site-like extra point). OpenMM's force pipeline
//  redistributes the site's force onto its parents AND leaves the original
//  force in the site's own slot (ReferenceVirtualSites::distributeForces /
//  CudaIntegrationUtilities::distributeForcesFromVirtualSites both read-then-
//  add without zeroing the source slot); re-reducing that leftover slot would
//  double-count it. This test asserts the reduced body wrench equals an
//  INDEPENDENT sum over real atoms only -- computed without calling
//  reduceAtomForcesToBodies -- so a missing skip in either implementation is
//  caught.
//
//  Two bodies (nontrivial rotation + translation each), one real-atom-only
//  body and one body carrying a 3-parent virtual site, driven by an ordinary
//  NoCutoff NonbondedForce (LJ + Coulomb) so the per-atom forces are a real
//  OpenMM force-field evaluation, not a hand-typed array.
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <gtest/gtest.h>
#include <vector>

#include "OpenMMContext.hpp"
#include "bridge/ForceReducer.hpp"
#include "robot_math.hpp"

namespace {

using robo::Real;
using robo::Rotation;
using robo::SpatialVec;
using robo::Transform;
using robo::Vec3;
using robo::XAxis;
using robo::YAxis;

// Absolute-floor + relative tolerance: forces here span ~1-1e3 kJ/mol/nm, and
// the CUDA half accumulates in a mixed-precision fixed-point buffer (see
// SPLIT-DEDUP-FORCEREDUCER.md "tight band, not exact-bit"), so a single fixed
// absolute epsilon is either too loose for small components or too tight for
// large ones.
constexpr double kAbsFloor = 1e-4;
constexpr double kRel = 1e-5;

void expectNear(double a, double b, const char* what) {
    const double tol = std::max(kAbsFloor, kRel * std::max(std::fabs(a), std::fabs(b)));
    EXPECT_NEAR(a, b, tol) << what;
}

void expectVec3Near(const Vec3& a, const Vec3& b, const char* what) {
    expectNear(a[0], b[0], what);
    expectNear(a[1], b[1], what);
    expectNear(a[2], b[2], what);
}

// 6 atoms / 2 bodies: body 0 = {0,1} (plain LJ+charge, no vsite); body 1 =
// {2,3,4} (parents) + {5} (massless ThreeParticleAverageSite EP, M-site-like
// charge, no LJ). Positions come from FIXED body-frame stations pushed
// through a nontrivial per-body rotation + translation (X_GB), so the host
// and device paths consume the identical Ground positions.
struct Fixture {
    static constexpr int kNumAtoms = 6;
    static constexpr int kNumBodies = 2;

    std::vector<int> atomBody{0, 0, 1, 1, 1, 1};
    std::vector<Real> atomMass{12.0, 12.0, 12.0, 1.0, 1.0, 0.0}; // atom 5 = EP
    std::vector<Vec3> stationB{
        Vec3(0.00, 0.00, 0.00),  // atom 0 (body 0)
        Vec3(0.32, 0.03, -0.02), // atom 1 (body 0)
        Vec3(0.00, 0.00, 0.00),  // atom 2 (body 1, vsite parent 1)
        Vec3(0.35, 0.00, 0.00),  // atom 3 (body 1, vsite parent 2)
        Vec3(0.10, 0.35, 0.00),  // atom 4 (body 1, vsite parent 3)
        Vec3(0.00, 0.00, 0.00),  // atom 5 (EP; station unused -- OpenMM
                                 // recomputes its position from parents 2/3/4)
    };
    std::vector<Transform> X_GB{
        Transform(Rotation(0.6, XAxis), Vec3(0.0, 0.0, 0.0)),
        Transform(Rotation(-0.4, YAxis), Vec3(0.6, 0.1, -0.05)),
    };
    // ThreeParticleAverageSite weights (sum to 1), parents = atoms 2,3,4.
    static constexpr double kW1 = 0.5, kW2 = 0.3, kW3 = 0.2;

    [[nodiscard]] auto atomPosG() const -> std::vector<Vec3> {
        std::vector<Vec3> pos(static_cast<std::size_t>(kNumAtoms));
        for (int a = 0; a < kNumAtoms; ++a) {
            pos[static_cast<std::size_t>(a)] =
                X_GB[static_cast<std::size_t>(atomBody[static_cast<std::size_t>(a)])]
                * stationB[static_cast<std::size_t>(a)];
        }
        return pos;
    }

    [[nodiscard]] auto buildSystemTopology() const -> SystemTopology {
        SystemTopology s;
        s.numAtoms = kNumAtoms;
        s.numMolecules = kNumAtoms; // no intramolecular bonds/angles/torsions needed
        s.nonbondedMethod = NonbondedMethod::NoCutoff;
        s.nonbondedCutoff = 2.0;

        const std::vector<double> charge{+0.30, -0.25, +0.20, +0.15, +0.15, -0.80};
        const std::vector<double> sigma{0.30, 0.32, 0.30, 0.24, 0.24, 0.00};
        const std::vector<double> epsilon{0.30, 0.28, 0.30, 0.10, 0.10, 0.00};

        const std::vector<Vec3> pos = atomPosG();
        for (int a = 0; a < kNumAtoms; ++a) {
            s.atomsMass.push_back(atomMass[static_cast<std::size_t>(a)]);
            s.atomsCharge.push_back(charge[static_cast<std::size_t>(a)]);
            s.atomsSigma.push_back(sigma[static_cast<std::size_t>(a)]);
            s.atomsEpsilon.push_back(epsilon[static_cast<std::size_t>(a)]);
            s.atomsRadius.push_back(0.0);
            s.atomsScreen.push_back(0.0);
            s.atomsX.push_back(pos[static_cast<std::size_t>(a)][0]);
            s.atomsY.push_back(pos[static_cast<std::size_t>(a)][1]);
            s.atomsZ.push_back(pos[static_cast<std::size_t>(a)][2]);
        }

        // The massless EP (atom 5) as a 3-particle average site of atoms 2,3,4.
        s.numVirtualSites = 1;
        s.vsSite = {5};
        s.vsAtom1 = {2};
        s.vsAtom2 = {3};
        s.vsAtom3 = {4};
        s.vsWeight1 = {kW1};
        s.vsWeight2 = {kW2};
        s.vsWeight3 = {kW3};

        // Exclude the EP from its own parents (same "molecule" -- otherwise the
        // near-coincident EP/parent Coulomb term dominates and is physically
        // meaningless for a bonded extra point).
        s.numExclusions = 3;
        s.exclusionI = {2, 3, 4};
        s.exclusionJ = {5, 5, 5};

        return s;
    }
};

} // namespace

// Host half: always runs. Drives an OpenMM NonbondedForce evaluation for a
// FIXED, reproducible set of per-atom Ground positions, feeds the resulting
// per-atom forces through reduceAtomForcesToBodies, and checks the result
// against an INDEPENDENT sum over real atoms only (INV-2).
TEST(ForceReducer, HostReductionMatchesRealAtomOracle) {
    Fixture fx;
    auto& omm = OpenMMContext::get();
    omm.shutdown();
    omm.initialize(fx.buildSystemTopology());

    const std::vector<Vec3> posG = fx.atomPosG();
    std::vector<OpenMM::Vec3> posCache(Fixture::kNumAtoms);
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        posCache[static_cast<std::size_t>(a)] =
            OpenMM::Vec3(posG[static_cast<std::size_t>(a)][0], posG[static_cast<std::size_t>(a)][1],
                        posG[static_cast<std::size_t>(a)][2]);
    }

    std::vector<OpenMM::Vec3> rawForces;
    omm.evaluateForcesFromPositionsCache(posCache, rawForces);
    ASSERT_EQ(rawForces.size(), static_cast<std::size_t>(Fixture::kNumAtoms));

    std::vector<Vec3> forceG(Fixture::kNumAtoms);
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        forceG[static_cast<std::size_t>(a)] =
            Vec3(rawForces[static_cast<std::size_t>(a)][0], rawForces[static_cast<std::size_t>(a)][1],
                rawForces[static_cast<std::size_t>(a)][2]);
    }

    // Sanity: the EP's own slot carries a nonzero leftover force (the thing
    // INV-2 exists to skip). If this is ever zero the fixture stopped
    // exercising the skip and the test below would pass vacuously.
    const Vec3& epForce = forceG[5];
    EXPECT_GT(epForce.norm(), 1.0) << "fixture no longer gives the EP a nonzero own-slot force";

    std::vector<SpatialVec> hostBF(Fixture::kNumBodies, SpatialVec(Vec3(0), Vec3(0)));
    reduceAtomForcesToBodies(forceG.data(), posG.data(), fx.atomMass.data(), fx.atomBody.data(), fx.X_GB.data(),
                             Fixture::kNumAtoms, Fixture::kNumBodies, hostBF.data());

    // Independent oracle: sum over REAL atoms only (mass != 0), NOT via
    // reduceAtomForcesToBodies -- exercises the same INV-2 skip from scratch.
    std::vector<SpatialVec> oracleBF(Fixture::kNumBodies, SpatialVec(Vec3(0), Vec3(0)));
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        if (fx.atomMass[static_cast<std::size_t>(a)] == Real(0)) {
            continue;
        }
        const int b = fx.atomBody[static_cast<std::size_t>(a)];
        const Vec3 r = posG[static_cast<std::size_t>(a)] - fx.X_GB[static_cast<std::size_t>(b)].p();
        oracleBF[static_cast<std::size_t>(b)].linear += forceG[static_cast<std::size_t>(a)];
        oracleBF[static_cast<std::size_t>(b)].angular += r % forceG[static_cast<std::size_t>(a)];
    }

    for (int b = 0; b < Fixture::kNumBodies; ++b) {
        expectVec3Near(hostBF[static_cast<std::size_t>(b)].linear, oracleBF[static_cast<std::size_t>(b)].linear,
                       "linear (body, host vs real-atom-only oracle)");
        expectVec3Near(hostBF[static_cast<std::size_t>(b)].angular,
                       oracleBF[static_cast<std::size_t>(b)].angular,
                       "angular (body, host vs real-atom-only oracle)");
    }

    omm.shutdown();
}

// CUDA half: drives the SAME fixture through OpenMMContext's real device
// path (ensureKinematicsConstants -> pushBodyTransforms ->
// computeForcesAndEnergyOnDevice -> reduceForcesToBodies, i.e. the
// production reduceForces kernel in src/OpenMMContext.cpp) and asserts the
// per-body wrenches agree with the host reduction to a tight (not exact-bit)
// tolerance. GTEST_SKIP()-gated on non-CUDA builds; the host half above
// still runs there.
TEST(ForceReducer, CudaKernelMatchesHostReduction) {
#if !USE_CUDA
    GTEST_SKIP() << "CUDA-only: reduceForces kernel differential test";
#else
    Fixture fx;
    auto& omm = OpenMMContext::get();
    omm.shutdown();
    omm.initialize(fx.buildSystemTopology());

    const std::vector<Vec3> posG = fx.atomPosG();
    std::vector<OpenMM::Vec3> posCache(Fixture::kNumAtoms);
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        posCache[static_cast<std::size_t>(a)] =
            OpenMM::Vec3(posG[static_cast<std::size_t>(a)][0], posG[static_cast<std::size_t>(a)][1],
                        posG[static_cast<std::size_t>(a)][2]);
    }
    std::vector<OpenMM::Vec3> rawForces;
    omm.evaluateForcesFromPositionsCache(posCache, rawForces);
    std::vector<Vec3> forceG(Fixture::kNumAtoms);
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        forceG[static_cast<std::size_t>(a)] =
            Vec3(rawForces[static_cast<std::size_t>(a)][0], rawForces[static_cast<std::size_t>(a)][1],
                rawForces[static_cast<std::size_t>(a)][2]);
    }
    std::vector<SpatialVec> hostBF(Fixture::kNumBodies, SpatialVec(Vec3(0), Vec3(0)));
    reduceAtomForcesToBodies(forceG.data(), posG.data(), fx.atomMass.data(), fx.atomBody.data(), fx.X_GB.data(),
                             Fixture::kNumAtoms, Fixture::kNumBodies, hostBF.data());

    // Drive the real device pipeline: K1 (push X_GB*station into posq), the
    // system's registered forces on device, then K2 (reduceForces).
    static int worldTokenTag = 0;
    std::vector<double> stationFlat(3 * Fixture::kNumAtoms, 0.0);
    std::vector<int> isVirtual(Fixture::kNumAtoms, 0);
    isVirtual[5] = 1;
    for (int a = 0; a < Fixture::kNumAtoms; ++a) {
        const Vec3& st = fx.stationB[static_cast<std::size_t>(a)];
        stationFlat[static_cast<std::size_t>(3 * a)] = st[0];
        stationFlat[static_cast<std::size_t>(3 * a + 1)] = st[1];
        stationFlat[static_cast<std::size_t>(3 * a + 2)] = st[2];
    }
    const std::vector<int> bodyAtomsBeg{0, 2};
    const std::vector<int> bodyAtomsEnd{2, 6};
    const std::vector<int> bodyAtoms{0, 1, 2, 3, 4, 5};

    omm.ensureKinematicsConstants(static_cast<const void*>(&worldTokenTag), Fixture::kNumAtoms,
                                  Fixture::kNumBodies, stationFlat.data(), fx.atomBody.data(), isVirtual.data(),
                                  bodyAtomsBeg.data(), bodyAtomsEnd.data(), bodyAtoms.data(),
                                  /*stationsChanged=*/true);

    std::vector<double> xgbFlat(12 * Fixture::kNumBodies, 0.0);
    for (int b = 0; b < Fixture::kNumBodies; ++b) {
        const std::array<Real, 9>& e = fx.X_GB[static_cast<std::size_t>(b)].R().elems;
        double* d = xgbFlat.data() + 12 * b;
        for (int k = 0; k < 9; ++k) {
            d[k] = e[static_cast<std::size_t>(k)];
        }
        const Vec3& p = fx.X_GB[static_cast<std::size_t>(b)].p();
        d[9] = p[0];
        d[10] = p[1];
        d[11] = p[2];
    }
    omm.pushBodyTransforms(xgbFlat.data());
    omm.computeForcesAndEnergyOnDevice();

    std::vector<double> deviceBodyForce(6 * Fixture::kNumBodies, 0.0);
    omm.reduceForcesToBodies(deviceBodyForce.data());

    for (int b = 0; b < Fixture::kNumBodies; ++b) {
        const double* f = deviceBodyForce.data() + 6 * b;
        const Vec3 devAngular(f[0], f[1], f[2]);
        const Vec3 devLinear(f[3], f[4], f[5]);
        expectVec3Near(devLinear, hostBF[static_cast<std::size_t>(b)].linear, "linear (body, CUDA vs host)");
        expectVec3Near(devAngular, hostBF[static_cast<std::size_t>(b)].angular, "angular (body, CUDA vs host)");
    }

    omm.shutdown();
#endif
}
