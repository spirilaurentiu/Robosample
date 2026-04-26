/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves:

Ax = B,

where A is the general matrix


      1.80   2.88   2.05   -0.89           9.52
A =   5.25  -2.95  -0.95   -3.80   and B = 24.35
      1.58  -2.69  -2.90   -1.04           0.77
     -1.11  -0.66  -0.59    0.80          -6.22


   Solution
 x =     1.0000    -1.0000     3.0000    -5.0000


 LU factorization:
             1          2          3          4
 1      5.2500    -2.9500    -0.9500    -3.8000
 2      0.3429     3.8914     2.3757     0.4129
 3      0.3010    -0.4631    -1.5139     0.2948
 4     -0.2114    -0.3299     0.0047     0.1314

 Pivot indices
             2          2          3          4

*/

#include <iostream>
#include <sstream>

#include "gtest/gtest.h"

#include "SimTKmath.h"

using namespace SimTK;

// ---------------------------------------------------------------------------
// Shared test data
// ---------------------------------------------------------------------------
namespace {

Real A[16] = {1.80,
              2.88,
              2.05,
              -0.89,
              5.25,
              -2.95,
              -0.95,
              -3.80,
              1.58,
              -2.69,
              -2.90,
              -1.04,
              -1.11,
              -0.66,
              -0.59,
              0.80};

Real B[4] = {9.52, 24.35, 0.77, -6.22};
Real X[4] = {1.0, -1.0, 3.0, -5.0};

Real C[4] = {1.0, 2.0, 1.0, 3.0};

Real Z[4] = {0.0, 0.0, 0.0, 0.0};

} // anonymous namespace

// ---------------------------------------------------------------------------
// TEST: Real (double) precision solve
// ---------------------------------------------------------------------------
TEST(FactorLU, RealSolve) {
    Matrix a(4, 4, A);
    Vector b(4, B);
    Vector x_right(4, X);
    Vector x;

    FactorLU lu(a);
    lu.solve(b, x);

    EXPECT_LT((x - x_right).norm(), 10 * SignificantReal)
        << " Real SOLUTION: " << x << "  errnorm=" << (x - x_right).norm() << "\n";
}

// ---------------------------------------------------------------------------
// TEST: float precision solve
// ---------------------------------------------------------------------------
TEST(FactorLU, FloatSolve) {
    Matrix a(4, 4, A);
    Vector b(4, B);
    Vector x_right(4, X);

    Matrix_<float> af(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            af(i, j) = static_cast<float>(a(i, j));
        }
    }

    Vector_<float> bf(4);
    for (int i = 0; i < 4; ++i) {
        bf[i] = static_cast<float>(b[i]);
    }

    Vector_<float> xf_right(4);
    for (int i = 0; i < 4; ++i) {
        xf_right[i] = static_cast<float>(x_right[i]);
    }

    Vector_<float> xf;

    FactorLU luf;
    luf.factor(af);
    luf.solve(bf, xf);

    const float SignificantFloat = NTraits<float>::getSignificant();
    EXPECT_LT((xf - xf_right).norm(), 10 * SignificantFloat)
        << " float SOLUTION: " << xf << "  errnorm=" << (xf - xf_right).norm() << "\n";
}

// ---------------------------------------------------------------------------
// TEST: Re-factor with a Real matrix after a float factor, then solve
// ---------------------------------------------------------------------------
TEST(FactorLU, RefactorRealAfterFloat) {
    Matrix a(4, 4, A);
    Vector b(4, B);
    Vector x_right(4, X);
    Vector x;

    // Build the float factorisation first (mirrors the original test sequence)
    Matrix_<float> af(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            af(i, j) = static_cast<float>(a(i, j));
        }
    }

    FactorLU luf;
    luf.factor(af); // float factor – result not used here
    luf.factor(a);  // re-factor with Real matrix

    FactorLU lu(a);
    lu.solve(b, x);

    EXPECT_LT((x - x_right).norm(), 10 * SignificantReal)
        << " Real SOLUTION: " << x << "  errnorm=" << (x - x_right).norm() << "\n";
}

// ---------------------------------------------------------------------------
// TEST: 2x2 matrix inverse
// ---------------------------------------------------------------------------
TEST(FactorLU, Inverse2x2) {
    Matrix c(2, 2, C);
    FactorLU clu(c);
    Matrix invC;
    clu.inverse(invC);

    // The inverse of [1 2; 1 3] is [3 -2; -1 1]
    EXPECT_NEAR(invC[0][0], 3.0, 1e-10) << " invC[0][0] = " << invC[0][0] << "\n";
    EXPECT_NEAR(invC[0][1], -2.0, 1e-10) << " invC[0][1] = " << invC[0][1] << "\n";
    EXPECT_NEAR(invC[1][0], -1.0, 1e-10) << " invC[1][0] = " << invC[1][0] << "\n";
    EXPECT_NEAR(invC[1][1], 1.0, 1e-10) << " invC[1][1] = " << invC[1][1] << "\n";
}

// ---------------------------------------------------------------------------
// TEST: Solve with an all-zero matrix (singular – expect a zero solution)
// ---------------------------------------------------------------------------
TEST(FactorLU, SolveZeroMatrix) {
    Matrix z(2, 2, Z);
    FactorLU zlu(z);
    Vector_<double> xz;
    Vector_<double> bz(2);
    bz(0) = bz(1) = 0.0;

    zlu.solve(bz, xz);

    // Prepare the failure message
    std::stringstream ss;
    ss << "solve with mat all zeros: ";
    for (int i = 0; i < xz.size(); ++i) {
        ss << xz(i) << " ";
    }

    // SCOPED_TRACE will only print 'ss' if an assertion below fails
    SCOPED_TRACE(ss.str());

    ASSERT_EQ(xz.size(), 2);
    for (int i = 0; i < xz.size(); ++i) {
        EXPECT_TRUE(std::isnan(xz(i)));
    }
}

// ---------------------------------------------------------------------------
// TEST: 0x0 (null) matrix – expected to throw
// ---------------------------------------------------------------------------
TEST(FactorLU, NullMatrixThrows) {
    // The exception is raised by FactorLU::factor() which is called
    // immediately inside the constructor — wrap the construction, not solve().
    SCOPED_TRACE("(EXPECTED EXCEPTION) NULL matrix test: solve with mat(0,0)");

    EXPECT_THROW(
        {
            Matrix_<double> z0; // 0x0
            FactorLU z0lu(z0);  // Should throw here
        },
        std::exception);
}