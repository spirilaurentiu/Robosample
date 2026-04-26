/* -------------------------------------------------------------------------- *
 *                         SimTK Simbody: SimTKmath                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */
#include <gtest/gtest.h>
#include <iostream>

#include "SimTKlapack.h"
#include "SimTKmath.h"

using SimTK::Differentiator;
using SimTK::Matrix;
using SimTK::Real;
using SimTK::Vector;

class MyVectorFunc : public Differentiator::JacobianFunction {
    public:
    MyVectorFunc(int nf, int ny)
        : Differentiator::JacobianFunction(nf, ny)
        , time(0) {
    }

    void setTime(Real t) {
        time = t;
    }
    Real getTime() const {
        return time;
    }

    int f(const Vector& y, Vector& fy) const override;

    private:
    Real time;
};

class MyObjectiveFunc : public Differentiator::GradientFunction {
    public:
    MyObjectiveFunc(int ny)
        : Differentiator::GradientFunction(ny)
        , time(0) {
    }

    void setTime(Real t) {
        time = t;
    }
    Real getTime() const {
        return time;
    }

    int f(const Vector& y, Real& fy) const override;

    private:
    Real time;
};

class GenericScalarFunc : public Differentiator::ScalarFunction {
    typedef Real (*CFunc)(Real);

    public:
    GenericScalarFunc(CFunc cf)
        : Differentiator::ScalarFunction()
        , cp(cf) {
    }

    int f(Real x, Real& fx) const override {
        fx = cp(x);
        return 0;
    }

    CFunc cp;
};

class SinOmegaX : public Differentiator::ScalarFunction {
    public:
    SinOmegaX(Real omega, Real acc)
        : ScalarFunction(acc)
        , w(omega) {
    }

    Real calc(Real x) const {
        return std::sin(w * x);
    }
    Real calcD1(Real x) const {
        return w * std::cos(w * x);
    }
    Real calcD2(Real x) const {
        return -w * w * std::sin(w * x);
    }
    // Must provide this virtual function.
    int f(Real x, Real& fx) const override {
        Real ffx = (Real)calc(x);
        fx = ffx;
        return 0; // success
    }

    private:
    const Real w;
};

class Cubic : public Differentiator::ScalarFunction {
    public:
    // ax^3+bx^2+cx+d
    Cubic(Real aa, Real bb, Real cc, Real dd, Real acc)
        : ScalarFunction(acc)
        , a(aa)
        , b(bb)
        , c(cc)
        , d(dd) {
    }

    Real calc(Real x) const {
        return (a * x * x * x + b * x * x + c * x + d) * std::exp(a * x);
    }

    Real calcD1(Real x) const {
        return (3 * a * x * x + 2 * b * x + c) * std::exp(a * x) + a * calc(x);
    }

    Real calcD2(Real x) const {
        return (6 * a * x + 2 * b) * std::exp(a * x) + (3 * a * x * x + 2 * b * x + c) * a * std::exp(a * x)
               + a * calcD1(x);
    }

    // Must provide this virtual function.
    int f(Real x, Real& fx) const override {
        Real ffx = (Real)calc(x);
        fx = ffx;
        return 0; // success
    }

    private:
    const Real a, b, c, d;
};

// ===================== Pendulum =====================

static int pendODE(Real, const Vector& yy, Vector& fy) {
    const Real g = 13.7503716373294544;

    const Real x = yy[0];
    const Real y = yy[1];
    const Real xd = yy[2];
    const Real yd = yy[3];

    const Real tmp = xd * xd + yd * yd - g * y;

    fy[0] = xd;
    fy[1] = yd;
    fy[2] = -x * tmp;
    fy[3] = -y * tmp - g;

    return 0;
}

int MyVectorFunc::f(const Vector& yy, Vector& fy) const {
    return pendODE(getTime(), yy, fy);
}

int MyObjectiveFunc::f(const Vector& yy, Real& fy) const {
    Vector tmp(4);
    const int res = pendODE(getTime(), yy, tmp);
    fy = tmp.norm();
    return res;
}

// ===================== Tests =====================

TEST(DifferentiatorTest, ScalarDerivative) {
    auto mysin = [](Real x) {
        return std::sin(x);
    };
    GenericScalarFunc gf(mysin);

    Differentiator dsin(gf, Differentiator::CentralDifference);

    const Real x = SimTK_PI / 6; // 30 deg
    const Real approx = dsin.calcDerivative(x);
    const Real exact = std::cos(x);

    EXPECT_NEAR(approx, exact, 1e-6) << "sin' approx=" << approx << " exact=" << exact << "\n";
}

TEST(DifferentiatorTest, GradientConsistency) {
    MyObjectiveFunc sf(4);
    Differentiator gradf(sf);

    const Real rp[] = {.01, .02, .03, -.14};
    Vector y0(4, rp);
    Vector dy(4, rp);

    Real f0;
    sf.f(y0, f0);

    Vector grad;
    gradf.calcGradient(y0, f0, grad);

    Real f1;
    sf.f(y0 + dy, f1);

    const Real linearApprox = f0 + ~grad * dy;

    EXPECT_NEAR(f1, linearApprox, 1e-2) << "f(y0+dy)=" << f1 << " linear approx=" << linearApprox << "\n";
}

TEST(DifferentiatorTest, JacobianConsistency) {
    MyVectorFunc vf(4, 4);
    Differentiator df(vf);

    const Real rp[] = {.01, .02, .03, -.14};
    Vector y0(4, rp);
    Vector dy(4, rp);

    Vector f0(4);
    vf.f(y0, f0);

    Matrix J;
    df.calcJacobian(y0, f0, J);

    Vector f1(4);
    vf.f(y0 + dy, f1);

    Vector linearApprox = f0 + J * dy;

    EXPECT_NEAR((f1 - linearApprox).norm(), 0.0, 1e-2) << "||error||=" << (f1 - linearApprox).norm() << "\n";
}

TEST(DifferentiatorTest, SinOmegaX_Derivatives) {
    const Real w = 0.5;
    const Real acc = 1e-8;

    SinOmegaX func(w, acc);
    Differentiator d(func, Differentiator::CentralDifference);

    const Real x = 0.7;

    const Real d1_exact = func.calcD1(x);
    const Real d1_num = d.calcDerivative(x);

    EXPECT_NEAR(d1_num, d1_exact, 1e-6);

    // crude second derivative via differentiating first derivative numerically
    const Real h = 1e-5;
    const Real d2_num = (d.calcDerivative(x + h) - d.calcDerivative(x - h)) / (2 * h);
    const Real d2_exact = func.calcD2(x);

    EXPECT_NEAR(d2_num, d2_exact, 1e-4) << "SinOmegaX d2 approx=" << d2_num << " exact=" << d2_exact << "\n";
}

TEST(DifferentiatorTest, SinOmegaX_ErrorSweep) {
    for (int digits = 0; digits <= 10; ++digits) {
        Real acc;
        if (digits < 40) {
            acc = std::pow(10., -(digits / (1.5 * sizeof(double) / sizeof(Real))));
        } else if (digits == 40) {
            acc = SimTK::NTraits<Real>::getSignificant();
        } else {
            acc = SimTK::NTraits<Real>::getEps();
        }

        const Real w = 0.01;
        SinOmegaX func(w, acc);
        Differentiator d(func);

        const int N = 200; // reduced from 1000 for unit test
        const Real offs = 0.1;
        const Real inc = (Real)SimTK_PI / N;

        Real err1rms = 0, err1max = 0;
        Real err2rms = 0, err2max = 0;

        for (int i = 0; i < N; ++i) {
            const Real x = offs + i * inc;

            const Real exact = func.calcD1(x);
            const Real d1 = d.calcDerivative(x); // default (usually forward)
            const Real d2 = d.calcDerivative(x, Differentiator::CentralDifference);

            const Real e1 = std::abs((d1 - exact) / exact);
            const Real e2 = std::abs((d2 - exact) / exact);

            err1rms += e1 * e1;
            err2rms += e2 * e2;

            err1max = std::max(err1max, e1);
            err2max = std::max(err2max, e2);
        }

        err1rms = std::sqrt(err1rms / N);
        err2rms = std::sqrt(err2rms / N);

        // Central difference should be at least as good as forward
        EXPECT_LE(err2rms, err1rms * 1.2)
            << " acc=" << acc << " err1 rms=" << err1rms << " err2 rms=" << err2rms << "\n";

        // Errors should remain bounded (sanity regression)
        EXPECT_LT(err2rms, 1e-2) << " acc=" << acc << " err2 rms=" << err2rms << "\n";
        EXPECT_LT(err1rms, 1e-1) << " acc=" << acc << " err1 rms=" << err1rms << "\n";
    }
}

TEST(DifferentiatorTest, Cubic_Derivatives) {
    const Real a = -1.0;
    const Real b = -2.0;
    const Real c = 3.0;
    const Real d0 = 4.0;
    const Real acc = 1e-8;

    Cubic func(a, b, c, d0, acc);
    Differentiator d(func, Differentiator::CentralDifference);

    const Real x = 0.3;

    const Real d1_exact = func.calcD1(x);
    const Real d1_num = d.calcDerivative(x);

    EXPECT_NEAR(d1_num, d1_exact, 1e-6);

    const Real h = 1e-5;
    const Real d2_num = (d.calcDerivative(x + h) - d.calcDerivative(x - h)) / (2 * h);
    const Real d2_exact = func.calcD2(x);

    EXPECT_NEAR(d2_num, d2_exact, 1e-4) << "Cubic d2 approx=" << d2_num << " exact=" << d2_exact << "\n";
}
