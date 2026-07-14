#pragma once
// ============================================================================
//  PeriodicBox.hpp -- the OpenMM-free periodic-boundary primitives for explicit
//  solvent: the minimum-image (triclinic wrap) and the reduced box-vector
//  construction.
//
//  SINGLE SOURCE OF TRUTH. The minimum-image wrap was a file-local lambda inside
//  src/Context.cpp's clash scan; the reduced-box-vector construction was the pure
//  body of OpenMMContext::computePeriodicBoxVectors_Context. Both are extracted
//  here verbatim (only wrapped as free functions) so the production code and the
//  unit tests bind the SAME logic -- the hoisting pattern already used for
//  robo_linalg / EngineHelpers. Nothing here touches OpenMM, so it can be unit
//  tested with zero forcefield dependency.
//
//  CONVENTION. Box vectors are REDUCED (lower-triangular):
//      a = (ax, 0,  0 )   b = (bx, by, 0 )   c = (cx, cy, cz)
//  stored as a flat 9-array  bv = [ax,ay,az, bx,by,bz, cx,cy,cz]  (a=[0..2],
//  b=[3..5], c=[6..8]) -- the same layout SystemTopology::boxVectors uses and the
//  order OpenMM reduces them. The diagonal is positive; off-diagonals satisfy the
//  reduction |bx| <= ax/2, |cx| <= ax/2, |cy| <= by/2.
// ============================================================================

#include <array>
#include <cmath>

namespace robo { namespace pbc {

/**
 * @brief Three reduced (lower-triangular) periodic-box lattice vectors, stored
 *        flat as [a.x a.y a.z  b.x b.y b.z  c.x c.y c.z] in nm.
 * @note Reduced form: a=(ax,0,0), b=(bx,by,0), c=(cx,cy,cz), positive diagonal,
 *       |bx|<=ax/2, |cx|<=ax/2, |cy|<=by/2. Same layout and order as
 *       SystemTopology::boxVectors and OpenMM's reduction.
 */
struct BoxVectors {
    /** @brief The nine components a,b,c stacked (nm). */
    std::array<double, 9> v{};

    /** @brief Borrowed read pointer to the nine components. */
    [[nodiscard]] auto data() const -> const double* {
        return v.data();
    }
    /** @brief Borrowed mutable pointer to the nine components. */
    [[nodiscard]] auto data() -> double* {
        return v.data();
    }
};

/**
 * @brief Wrap a displacement to its minimum image in a reduced cell, in place.
 * @param[in,out] dx,dy,dz displacement components (nm), overwritten with the
 *        minimum-image displacement.
 * @param[in] bv borrowed reduced box vectors, 9 doubles (BoxVectors layout).
 * @pre @p bv is reduced (lower-triangular); the single c->b->a subtraction pass
 *      yields the canonical minimum image only for a reduced cell.
 */
inline void minimumImage(double& dx, double& dy, double& dz, const double* bv) {
    // c then b then a  (bv layout: a=[0..2], b=[3..5], c=[6..8]).
    double n = std::round(dz / bv[8]);
    dx -= n * bv[6];
    dy -= n * bv[7];
    dz -= n * bv[8];
    n = std::round(dy / bv[4]);
    dx -= n * bv[3];
    dy -= n * bv[4];
    n = std::round(dx / bv[0]);
    dx -= n * bv[0];
}

/**
 * @brief Minimum-image distance between two points under a reduced cell.
 * @param[in] pa,pb borrowed 3-component point coordinates (nm).
 * @param[in] bv borrowed reduced box vectors, 9 doubles.
 * @return the minimum-image separation (nm).
 * @pre @p bv is reduced (see minimumImage()).
 */
inline double minimumImageDistance(const double* pa, const double* pb, const double* bv) {
    double dx = pa[0] - pb[0];
    double dy = pa[1] - pb[1];
    double dz = pa[2] - pb[2];
    minimumImage(dx, dy, dz, bv);
    return std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
}

/**
 * @brief Build reduced box vectors from cell lengths and angles.
 * @param[in] aLen,bLen,cLen cell edge lengths (nm).
 * @param[in] alpha,beta,gamma cell angles (radians): alpha=angle(b,c),
 *        beta=angle(a,c), gamma=angle(a,b).
 * @return reduced lower-triangular BoxVectors (nm): a along x, b in the xy-plane,
 *         c fixed by the angles, then off-diagonals reduced to smallest images.
 * @note This is the triclinic construction + reduction OpenMM applies. In the
 *       current tree it has no production caller: production consumes ParmEd's
 *       already-reduced box vectors (SystemTopology::boxVectors) directly. It is
 *       exercised by TestPeriodicBoundary as the reduction contract (see findings).
 */
inline BoxVectors
reducedBoxVectors(double aLen, double bLen, double cLen, double alpha, double beta, double gamma) {
    constexpr double TOL = 1e-6;
    std::array<double, 3> a{aLen, 0.0, 0.0};
    std::array<double, 3> b{bLen * std::cos(gamma), bLen * std::sin(gamma), 0.0};
    const double cx = cLen * std::cos(beta);
    const double cy = cLen * (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / std::sin(gamma);
    const double cz = std::sqrt((cLen * cLen) - (cx * cx) - (cy * cy));
    std::array<double, 3> c{cx, cy, cz};

    auto clampSmall = [&](std::array<double, 3>& v) {
        for (int i = 0; i < 3; ++i) {
            if (std::abs(v[static_cast<std::size_t>(i)]) < TOL) {
                v[static_cast<std::size_t>(i)] = 0.0;
            }
        }
    };
    clampSmall(a);
    clampSmall(b);
    clampSmall(c);

    auto axpy = [](std::array<double, 3>& dst, const std::array<double, 3>& src, double s) {
        dst[0] -= s * src[0];
        dst[1] -= s * src[1];
        dst[2] -= s * src[2];
    };
    if (b[1] != 0.0) {
        axpy(c, b, std::round(c[1] / b[1]));
    }
    if (a[0] != 0.0) {
        axpy(c, a, std::round(c[0] / a[0]));
    }
    if (a[0] != 0.0) {
        axpy(b, a, std::round(b[0] / a[0]));
    }

    BoxVectors out;
    out.v = {a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]};
    return out;
}

/**
 * @brief Diagonal (orthorhombic) reduced box of the given side lengths.
 * @param[in] lx,ly,lz box edge lengths along x, y, z (nm).
 * @return BoxVectors with a=(lx,0,0), b=(0,ly,0), c=(0,0,lz).
 */
inline BoxVectors orthorhombic(double lx, double ly, double lz) {
    BoxVectors out;
    out.v = {lx, 0, 0, 0, ly, 0, 0, 0, lz};
    return out;
}

}} // namespace robo::pbc