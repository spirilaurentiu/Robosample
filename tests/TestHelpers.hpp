// ============================================================================
//  TestHelpers.hpp -- shared utilities for the Robosample primitive test suite
//
//  Deterministic RNG (fixed seed, Rule 9), random generators for the robo::
//  value types, double-cover-aware comparators, and finite-difference helpers.
//  Everything here is pure and header-only so each test translation unit can
//  include it without a link dependency.
// ============================================================================
#pragma once

#include <array>
#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>
#include <random>

#include "robot_math.hpp"

namespace rtest {

using robo::Mat33;
using robo::Quat;
using robo::Real;
using robo::Rotation;
using robo::SpatialVec;
using robo::SymMat33;
using robo::Transform;
using robo::UnitVec3;
using robo::Vec3;
using robo::Vec4;

// ---- tolerances (named so intent is explicit at each call site) ------------
inline constexpr Real kTight = 1e-13; // machine-precision algebra
inline constexpr Real kAlg = 1e-12;   // analytic identities
inline constexpr Real kFD = 1e-6;     // central finite difference (~sqrt(eps)*scale)
inline constexpr Real kLoose = 1e-9;  // accumulated / eigen tolerances

// ---- deterministic RNG -----------------------------------------------------
// A single fixed seed makes every run bit-reproducible. Tests that need an
// independent stream pass their own seed.
class Rng {
    public:
    explicit Rng(std::uint64_t seed = 0x9E3779B97F4A7C15ULL)
        : gen_(seed) {
    }

    auto uniform(Real lo, Real hi) -> Real {
        return std::uniform_real_distribution<Real>(lo, hi)(gen_);
    }
    auto gaussian(Real mean = 0, Real sd = 1) -> Real {
        return std::normal_distribution<Real>(mean, sd)(gen_);
    }
    auto vec3(Real lo = -2, Real hi = 2) -> Vec3 {
        return Vec3(uniform(lo, hi), uniform(lo, hi), uniform(lo, hi));
    }
    // uniform direction on S^2 (Marsaglia)
    auto unitVec3() -> Vec3 {
        Real x, y, s;
        do {
            x = uniform(-1, 1);
            y = uniform(-1, 1);
            s = x * x + y * y;
        } while (s >= Real(1) || s < Real(1e-12));
        const Real f = Real(2) * std::sqrt(Real(1) - s);
        return Vec3(x * f, y * f, Real(1) - Real(2) * s);
    }
    // uniform on SO(3): unit quaternion via Shoemake -> rotation
    auto unitQuat() -> Quat {
        const Real u1 = uniform(0, 1), u2 = uniform(0, 1), u3 = uniform(0, 1);
        const Real s1 = std::sqrt(Real(1) - u1), s2 = std::sqrt(u1);
        const Real t2 = Real(2) * M_PI * u2, t3 = Real(2) * M_PI * u3;
        Quat q(s2 * std::cos(t3), s1 * std::sin(t2), s1 * std::cos(t2), s2 * std::sin(t3));
        q.normalize();
        return q;
    }
    auto rotation() -> Rotation {
        return Rotation::fromQuaternion(unitQuat());
    }

    // random SPD nxn (row-major) with eigenvalues in [lo,hi]; A = Q diag Q^T.
    void spd(int n, Real* out, Real lo = 0.3, Real hi = 5.0) {
        std::array<Real, 36> Q{};
        // random orthonormal Q via Gram-Schmidt on a Gaussian matrix
        for (int i = 0; i < n * n; ++i) {
            Q[i] = gaussian();
        }
        for (int c = 0; c < n; ++c) {
            for (int p = 0; p < c; ++p) {
                Real d = 0;
                for (int r = 0; r < n; ++r) {
                    d += Q[r * n + c] * Q[r * n + p];
                }
                for (int r = 0; r < n; ++r) {
                    Q[r * n + c] -= d * Q[r * n + p];
                }
            }
            Real nrm = 0;
            for (int r = 0; r < n; ++r) {
                nrm += Q[r * n + c] * Q[r * n + c];
            }
            nrm = std::sqrt(nrm);
            for (int r = 0; r < n; ++r) {
                Q[r * n + c] /= nrm;
            }
        }
        std::array<Real, 6> lam{};
        for (int k = 0; k < n; ++k) {
            lam[k] = uniform(lo, hi);
        }
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                Real acc = 0;
                for (int k = 0; k < n; ++k) {
                    acc += Q[i * n + k] * lam[k] * Q[j * n + k];
                }
                out[i * n + j] = acc;
            }
        }
    }

    private:
    std::mt19937_64 gen_;
};

// ---- small vector ops not on the value types ------------------------------
inline auto cross(const Vec3& a, const Vec3& b) -> Vec3 {
    return a % b;
}
inline auto vnorm(const Vec3& v) -> Real {
    return v.norm();
}

// ---- comparators -----------------------------------------------------------
inline ::testing::AssertionResult NearVec3(const Vec3& a, const Vec3& b, Real tol) {
    const Real d = (a - b).norm();
    if (d <= tol) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure()
           << "Vec3 differ by " << d << " (tol " << tol << "): a=(" << a[0] << "," << a[1] << "," << a[2]
           << ") b=(" << b[0] << "," << b[1] << "," << b[2] << ")";
}

inline ::testing::AssertionResult NearMat33(const Mat33& a, const Mat33& b, Real tol) {
    Real d = 0;
    for (int i = 0; i < 9; ++i) {
        const Real e = a.elems[static_cast<std::size_t>(i)] - b.elems[static_cast<std::size_t>(i)];
        d = std::max(d, std::abs(e));
    }
    if (d <= tol) {
        return ::testing::AssertionSuccess();
    }
    return ::testing::AssertionFailure() << "Mat33 max elem diff " << d << " (tol " << tol << ")";
}

// double-cover aware quaternion distance: min(|q-q'|, |q+q'|)
inline auto quatDist(const Quat& a, const Quat& b) -> Real {
    Real dp = 0, dm = 0;
    for (int i = 0; i < 4; ++i) {
        dp += (a.elems[static_cast<std::size_t>(i)] - b.elems[static_cast<std::size_t>(i)])
              * (a.elems[static_cast<std::size_t>(i)] - b.elems[static_cast<std::size_t>(i)]);
        dm += (a.elems[static_cast<std::size_t>(i)] + b.elems[static_cast<std::size_t>(i)])
              * (a.elems[static_cast<std::size_t>(i)] + b.elems[static_cast<std::size_t>(i)]);
    }
    return std::sqrt(std::min(dp, dm));
}

inline auto quatDist(const Vec4& a, const Vec4& b) -> Real {
    return quatDist(Quat(a[0], a[1], a[2], a[3]), Quat(b[0], b[1], b[2], b[3]));
}

// Frobenius diff of two row-major nxn
inline auto matDiff(const Real* a, const Real* b, int n) -> Real {
    Real d = 0;
    for (int i = 0; i < n * n; ++i) {
        d = std::max(d, std::abs(a[i] - b[i]));
    }
    return d;
}

// dense matrix * vector / matrix (row-major)
inline void matVec(const Real* A, const Real* x, int n, Real* y) {
    for (int i = 0; i < n; ++i) {
        Real s = 0;
        for (int j = 0; j < n; ++j) {
            s += A[i * n + j] * x[j];
        }
        y[i] = s;
    }
}
inline void matMat(const Real* A, const Real* B, int n, Real* C) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Real s = 0;
            for (int k = 0; k < n; ++k) {
                s += A[i * n + k] * B[k * n + j];
            }
            C[i * n + j] = s;
        }
    }
}

inline void identity(int n, Real* I) {
    for (int i = 0; i < n * n; ++i) {
        I[i] = 0;
    }
    for (int i = 0; i < n; ++i) {
        I[i * n + i] = 1;
    }
}

// ---- finite difference of a Vec3-valued function of a scalar --------------
template <class F>
inline auto fdVec3(F f, Real t, Real h = 1e-6) -> Vec3 {
    const Vec3 fp = f(t + h);
    const Vec3 fm = f(t - h);
    return (fp - fm) * (Real(1) / (Real(2) * h));
}

// ---------------------------------------------------------------------------
//  Slow statistical-tier tests GTEST_SKIP unless ROBOSAMPLE_SLOW_TESTS is set.
//  A bare `ctest` run then reports SKIPPED, which is easy to miss -- so also
//  print one prominent stderr line naming the test and how to enable it. The
//  authoritative gate (nox -s tests) always sets ROBOSAMPLE_SLOW_TESTS=1, so
//  this only fires on an incomplete/bare dev run; call it right next to
//  GTEST_SKIP(), never in place of it.
// ---------------------------------------------------------------------------
inline void warnSlowTierSkipped() {
    const ::testing::TestInfo* ti = ::testing::UnitTest::GetInstance()->current_test_info();
    std::fprintf(stderr,
                 "[SLOW TEST SKIPPED] %s.%s did not run -- set ROBOSAMPLE_SLOW_TESTS=1 to enable "
                 "(the authoritative `nox -s tests` gate always sets it; a bare run is INCOMPLETE)\n",
                 ti != nullptr ? ti->test_suite_name() : "?", ti != nullptr ? ti->name() : "?");
}

} // namespace rtest
