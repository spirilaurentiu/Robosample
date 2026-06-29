// ============================================================================
//  TestPeriodicBoundary.cpp -- explicit-solvent / periodic-boundary correctness.
//
//  Three OpenMM-free pieces govern PBC behaviour, and all three were untested:
//
//   1. The minimum-image (triclinic wrap) used by the clash scan to measure
//      cross-boundary distances under explicit solvent. This is the primitive
//      that decides whether two atoms on opposite faces of the box are treated as
//      neighbours (they are) or as box-length apart (the "phantom bond" bug the
//      Context.cpp / OpenMMContext.hpp comments warn about). Hoisted into
//      robo::pbc::minimumImage (PeriodicBox.hpp); the production clash scan now
//      calls the same function these tests exercise.
//
//   2. The reduced box-vector construction (lengths+angles -> lower-triangular
//      a,b,c), the pure math of OpenMMContext::computePeriodicBoxVectors_Context.
//
//   3. The DCDWriter periodic-box record: the CHARMM cosine convention written
//      into each frame's EXTRA_BLOCK, which trajectory tools read to reconstruct
//      the cell.
//
//  CORRECTNESS NOTE on minimum image. A single-pass c->b->a reduction is GLOBALLY
//  minimal for an orthorhombic cell, but NOT in general for a skewed triclinic
//  cell (the true nearest image can sit in a corner the single pass does not
//  visit). This matches OpenMM's reduced-box behaviour. So the tests assert:
//    - orthorhombic  -> the wrap equals the brute-force global minimum image;
//    - triclinic     -> the wrap is a VALID canonical image: lattice-translation
//                       invariant and reproducing the production algorithm,
//                       rather than (incorrectly) demanding global minimality.
// ============================================================================
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <gtest/gtest.h>
#include <vector>

#include "DCDWriter.hpp"
#include "PeriodicBox.hpp"
#include "TestHelpers.hpp"

using robo::pbc::BoxVectors;
using robo::pbc::minimumImage;
using robo::pbc::minimumImageDistance;
using robo::pbc::orthorhombic;
using robo::pbc::reducedBoxVectors;
using rtest::Rng;

namespace {

constexpr double kPi = 3.14159265358979323846;
inline double deg2rad(double d) {
    return d * kPi / 180.0;
}

// Brute-force global minimum-image distance: scan lattice images i,j,k in
// [-R,R]. With both points confined to the primary cell, R=2 is more than ample.
double bruteMinImageDistance(const std::array<double, 3>& pa,
                             const std::array<double, 3>& pb,
                             const double* bv,
                             int R) {
    double best = 1e30;
    for (int i = -R; i <= R; ++i) {
        for (int j = -R; j <= R; ++j) {
            for (int k = -R; k <= R; ++k) {
                const double ix = pb[0] + i * bv[0] + j * bv[3] + k * bv[6];
                const double iy = pb[1] + i * bv[1] + j * bv[4] + k * bv[7];
                const double iz = pb[2] + i * bv[2] + j * bv[5] + k * bv[8];
                const double dx = pa[0] - ix, dy = pa[1] - iy, dz = pa[2] - iz;
                best = std::min(best, std::sqrt((dx * dx) + (dy * dy) + (dz * dz)));
            }
        }
    }
    return best;
}

} // namespace

// ---------------------------------------------------------------------------
//  1. Orthorhombic wrap is the GLOBAL minimum image. Points confined to the
//     primary box; the single-pass result must match the brute-force search.
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, OrthorhombicMinimumImageIsGlobal) {
    Rng rng(0x9101);
    const double Lx = 3.0, Ly = 4.0, Lz = 5.0;
    const BoxVectors B = orthorhombic(Lx, Ly, Lz);
    const double* bv = B.data();

    double maxErr = 0;
    for (int t = 0; t < 20000; ++t) {
        const std::array<double, 3> pa{rng.uniform(0, Lx), rng.uniform(0, Ly), rng.uniform(0, Lz)};
        const std::array<double, 3> pb{rng.uniform(0, Lx), rng.uniform(0, Ly), rng.uniform(0, Lz)};
        const double dWrap = minimumImageDistance(pa.data(), pb.data(), bv);
        const double dBrute = bruteMinImageDistance(pa, pb, bv, 2);
        maxErr = std::max(maxErr, std::abs(dWrap - dBrute));
    }
    EXPECT_LT(maxErr, rtest::kTight) << "orthorhombic single-pass is not the global minimum image";
}

// ---------------------------------------------------------------------------
//  2. Canonical range: each wrapped component of an orthorhombic displacement
//     lands in [-L/2, L/2].
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, OrthorhombicWrapInCanonicalRange) {
    Rng rng(0x9102);
    const double Lx = 3.0, Ly = 4.0, Lz = 5.0;
    const BoxVectors B = orthorhombic(Lx, Ly, Lz);
    const double* bv = B.data();

    for (int t = 0; t < 20000; ++t) {
        double dx = rng.uniform(-12, 12), dy = rng.uniform(-12, 12), dz = rng.uniform(-12, 12);
        minimumImage(dx, dy, dz, bv);
        // allow a rounding sliver beyond exactly L/2.
        EXPECT_LE(std::abs(dx), Lx / 2 + 1e-9) << "x out of canonical range";
        EXPECT_LE(std::abs(dy), Ly / 2 + 1e-9) << "y out of canonical range";
        EXPECT_LE(std::abs(dz), Lz / 2 + 1e-9) << "z out of canonical range";
    }
}

// ---------------------------------------------------------------------------
//  3. Lattice-translation invariance (orthorhombic AND triclinic): adding any
//     whole combination of box vectors to a displacement does not change its
//     minimum image. This is the defining property of a correct wrap and holds
//     for triclinic cells even though global minimality does not.
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, WrapIsLatticeTranslationInvariant) {
    Rng rng(0x9103);
    const BoxVectors boxes[] = {
        orthorhombic(3.0, 4.0, 5.0),
        reducedBoxVectors(3.0, 3.2, 3.5, deg2rad(80), deg2rad(95), deg2rad(70)),
        reducedBoxVectors(4.0, 4.0, 4.0, deg2rad(88), deg2rad(92), deg2rad(89)),
    };
    for (const BoxVectors& B : boxes) {
        const double* bv = B.data();
        double maxViol = 0;
        for (int t = 0; t < 20000; ++t) {
            const double dx = rng.uniform(-6, 6), dy = rng.uniform(-6, 6), dz = rng.uniform(-6, 6);
            double wx = dx, wy = dy, wz = dz;
            minimumImage(wx, wy, wz, bv);

            // add a random integer lattice translation, then re-wrap.
            const int na = static_cast<int>(std::lround(rng.uniform(-3, 3)));
            const int nb = static_cast<int>(std::lround(rng.uniform(-3, 3)));
            const int nc = static_cast<int>(std::lround(rng.uniform(-3, 3)));
            double tx = dx + na * bv[0] + nb * bv[3] + nc * bv[6];
            double ty = dy + na * bv[1] + nb * bv[4] + nc * bv[7];
            double tz = dz + na * bv[2] + nb * bv[5] + nc * bv[8];
            minimumImage(tx, ty, tz, bv);

            maxViol = std::max({maxViol, std::abs(tx - wx), std::abs(ty - wy), std::abs(tz - wz)});
        }
        EXPECT_LT(maxViol, rtest::kAlg) << "wrap is not lattice-translation invariant";
    }
}

// ---------------------------------------------------------------------------
//  4. Explicit-solvent straddle: two atoms of one molecule sit just inside
//     opposite faces of the box. Their NAIVE separation is ~box, but the true
//     (minimum-image) separation is the small real gap across the boundary. This
//     is exactly the case the clash scan must get right -- otherwise a solvated
//     molecule straddling a face reports a phantom ~box-length bond (or misses a
//     real cross-boundary clash). The bug the Context.cpp comment warns about.
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, ExplicitSolventStraddleIsNeighbour) {
    const double L = 3.0;
    const BoxVectors B = orthorhombic(L, L, L);
    const double* bv = B.data();

    // gap of 0.08 across the x-boundary: atom A near x=0, atom B near x=L.
    const double gap = 0.08;
    const std::array<double, 3> atomA{0.02, 1.5, 1.5};
    const std::array<double, 3> atomB{0.02 + L - gap, 1.5, 1.5};

    const double naive = std::abs(atomA[0] - atomB[0]); // ~2.92
    const double mimg = minimumImageDistance(atomA.data(), atomB.data(), bv);

    EXPECT_GT(naive, L / 2) << "test setup: atoms should be naively far apart";
    EXPECT_NEAR(mimg, gap, 1e-12) << "minimum image failed to identify the cross-boundary neighbour";
}

// ---------------------------------------------------------------------------
//  5. Triclinic wrap is a valid canonical image that matches the production
//     algorithm. We do NOT assert global minimality (single-pass triclinic is
//     not globally minimal, by construction); we assert the wrapped displacement
//     equals an explicit re-implementation of the c->b->a reduction and is itself
//     a fixed point (re-wrapping does nothing).
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, TriclinicWrapMatchesReductionAndIsFixedPoint) {
    Rng rng(0x9105);
    const BoxVectors B = reducedBoxVectors(3.0, 3.2, 3.5, deg2rad(80), deg2rad(95), deg2rad(70));
    const double* bv = B.data();

    double maxAlgErr = 0, maxFixErr = 0;
    for (int t = 0; t < 20000; ++t) {
        double dx = rng.uniform(-6, 6), dy = rng.uniform(-6, 6), dz = rng.uniform(-6, 6);

        // explicit c->b->a reduction (independent re-derivation)
        double ex = dx, ey = dy, ez = dz;
        double n = std::round(ez / bv[8]);
        ex -= n * bv[6];
        ey -= n * bv[7];
        ez -= n * bv[8];
        n = std::round(ey / bv[4]);
        ex -= n * bv[3];
        ey -= n * bv[4];
        n = std::round(ex / bv[0]);
        ex -= n * bv[0];

        double wx = dx, wy = dy, wz = dz;
        minimumImage(wx, wy, wz, bv);
        maxAlgErr = std::max({maxAlgErr, std::abs(wx - ex), std::abs(wy - ey), std::abs(wz - ez)});

        // fixed point: wrapping an already-wrapped displacement changes nothing.
        double fx = wx, fy = wy, fz = wz;
        minimumImage(fx, fy, fz, bv);
        maxFixErr = std::max({maxFixErr, std::abs(fx - wx), std::abs(fy - wy), std::abs(fz - wz)});
    }
    EXPECT_LT(maxAlgErr, rtest::kTight) << "wrap disagrees with the c->b->a reduction";
    EXPECT_LT(maxFixErr, rtest::kTight) << "wrap is not idempotent on an already-canonical image";
}

// ---------------------------------------------------------------------------
//  6. Reduced box vectors: orthogonal angles produce a diagonal (orthorhombic)
//     cell, and the reduced form is lower-triangular with the off-diagonals
//     reduced to the smallest images (|bx|<=ax/2, |cx|<=ax/2, |cy|<=by/2).
//     Also: the cell VOLUME (det) matches the analytic triclinic volume.
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, ReducedBoxVectorsAreLowerTriangularAndReduced) {
    // orthogonal angles -> diagonal
    {
        const BoxVectors B = reducedBoxVectors(2.0, 3.0, 4.0, deg2rad(90), deg2rad(90), deg2rad(90));
        const double* bv = B.data();
        EXPECT_NEAR(bv[0], 2.0, 1e-12);
        EXPECT_NEAR(bv[4], 3.0, 1e-12);
        EXPECT_NEAR(bv[8], 4.0, 1e-12);
        for (int idx : {1, 2, 3, 5, 6, 7}) {
            EXPECT_NEAR(bv[idx], 0.0, 1e-9) << "off-diagonal " << idx << " nonzero";
        }
    }

    // general triclinic: lower-triangular + reduced + correct volume. Angles
    // chosen so the UNREDUCED c_y exceeds b_y/2 -- i.e. the c -= round(c_y/b_y) b
    // reduction is actually required (a no-skew box would not exercise it).
    {
        const double aL = 3.0, bL = 3.4, cL = 3.8;
        const double al = deg2rad(120), be = deg2rad(90), ga = deg2rad(90);
        const BoxVectors B = reducedBoxVectors(aL, bL, cL, al, be, ga);
        const double* bv = B.data();

        // lower-triangular: a has only x; b has only x,y.
        EXPECT_NEAR(bv[1], 0.0, 1e-12) << "a_y != 0";
        EXPECT_NEAR(bv[2], 0.0, 1e-12) << "a_z != 0";
        EXPECT_NEAR(bv[5], 0.0, 1e-12) << "b_z != 0";

        // reduced: off-diagonals are the smallest images.
        EXPECT_LE(std::abs(bv[3]), bv[0] / 2 + 1e-9) << "b_x not reduced";
        EXPECT_LE(std::abs(bv[6]), bv[0] / 2 + 1e-9) << "c_x not reduced";
        EXPECT_LE(std::abs(bv[7]), bv[4] / 2 + 1e-9) << "c_y not reduced";

        // volume = |det| = ax * by * cz (lower-triangular) == analytic triclinic V.
        const double detV = bv[0] * bv[4] * bv[8];
        const double ca = std::cos(al), cb = std::cos(be), cg = std::cos(ga);
        const double analytic = aL * bL * cL * std::sqrt(1 - ca * ca - cb * cb - cg * cg + 2 * ca * cb * cg);
        EXPECT_NEAR(detV, analytic, 1e-9) << "reduced-cell volume wrong";
    }
}

// ---------------------------------------------------------------------------
//  7. DCD periodic-box record round-trips through the CHARMM cosine convention.
//     Write one frame with a known box, read the 6-double unit-cell record back
//     out of the file, and check {a, cosAB, b, cosAC, cosBC, c}: orthogonal ->
//     cosines 0; gamma=60 deg -> cosAB = cos(60) = 0.5.
// ---------------------------------------------------------------------------
TEST(PeriodicBoundary, DcdBoxRecordUsesCharmmCosineConvention) {
    const char* path = "/tmp/test_pbc_box.dcd";
    const int kHeaderBytes = 276;               // CORD + title blocks
    const int kUcDataOffset = kHeaderBytes + 4; // skip the leading Fortran record marker (=48)

    auto writeAndReadBox = [&](const dcd::Box& box) -> std::array<double, 6> {
        {
            dcd::Writer w;
            w.initialize(path, /*numAtoms=*/2, /*withBox=*/true);
            const std::vector<double> coords = {0, 0, 0, 1, 1, 1};
            w.append(coords, box);
            w.close();
        }
        std::ifstream f(path, std::ios::binary);
        f.seekg(kUcDataOffset);
        std::array<double, 6> uc{};
        f.read(reinterpret_cast<char*>(uc.data()), 48);
        EXPECT_TRUE(f.good()) << "failed to read box record back";
        return uc;
    };

    // orthogonal 30x40x50: all cosines zero, sides preserved.
    {
        dcd::Box box;
        box.sideA = 30;
        box.sideB = 40;
        box.sideC = 50;
        box.angleAlpha = 90;
        box.angleBeta = 90;
        box.angleGamma = 90;
        const std::array<double, 6> uc = writeAndReadBox(box);
        EXPECT_NEAR(uc[0], 30.0, 1e-9); // a
        EXPECT_NEAR(uc[2], 40.0, 1e-9); // b
        EXPECT_NEAR(uc[5], 50.0, 1e-9); // c
        EXPECT_NEAR(uc[1], 0.0, 1e-12); // cosAB
        EXPECT_NEAR(uc[3], 0.0, 1e-12); // cosAC
        EXPECT_NEAR(uc[4], 0.0, 1e-12); // cosBC
    }

    // gamma = 60 deg: cosAB stored = cos(60) = 0.5; other angles still 90.
    {
        dcd::Box box;
        box.sideA = 30;
        box.sideB = 40;
        box.sideC = 50;
        box.angleAlpha = 90;
        box.angleBeta = 90;
        box.angleGamma = 60;
        const std::array<double, 6> uc = writeAndReadBox(box);
        EXPECT_NEAR(uc[1], 0.5, 1e-9) << "cosAB should be cos(60 deg) = 0.5";
        EXPECT_NEAR(uc[3], 0.0, 1e-12) << "cosAC should remain 0";
        EXPECT_NEAR(uc[4], 0.0, 1e-12) << "cosBC should remain 0";
    }

    // alpha = 120 deg: cosBC stored = cos(120) = -0.5 (sign carried through).
    {
        dcd::Box box;
        box.sideA = 30;
        box.sideB = 40;
        box.sideC = 50;
        box.angleAlpha = 120;
        box.angleBeta = 90;
        box.angleGamma = 90;
        const std::array<double, 6> uc = writeAndReadBox(box);
        EXPECT_NEAR(uc[4], -0.5, 1e-9) << "cosBC should be cos(120 deg) = -0.5";
    }
}