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

/**@file
 * This is a test program which uses the FactorSVD  class to compute
 * eigen values and eigen vectors
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves for the singular values and vectors for the
following matrix



A = 2.27   0.28  -0.48   1.07  -2.35   0.62
   -1.54  -1.67  -3.09   1.22   2.93  -7.39
    1.15   0.94   0.99   0.79  -1.45   1.03
   -1.94  -0.78  -0.21   0.63   2.30  -2.57



SOLUTION =
Singular values
     9.9966  3.6831  1.3569  0.5000
 Left singular vectors
          1       2       3       4
 1  -0.1921  0.8030  0.0041 -0.5642
 2   0.8794  0.3926 -0.0752  0.2587
 3  -0.2140  0.2980  0.7827  0.5027
 4   0.3795 -0.3351  0.6178 -0.6017

 Right singular vectors by row (first m rows of V**T)
          1       2       3       4       5       6
 1  -0.2774 -0.2020 -0.2918  0.0938  0.4213 -0.7816
 2   0.6003  0.0301 -0.3348  0.3699 -0.5266 -0.3353
 3  -0.1277  0.2805  0.6453  0.6781  0.0413 -0.1645
 4   0.1323  0.7034  0.1906 -0.5399 -0.0575 -0.3957

 Error estimate for the singular values
        1.1E-15

 Error estimates for the left singular vectors
        1.8E-16    4.8E-16    1.3E-15    1.3E-15

 Error estimates for the right singular vectors
        1.8E-16    4.8E-16    1.3E-15    2.2E-15
*/

#include <gtest/gtest.h>
#include <sstream>
#include <string>

#include "SimTKmath.h"

using namespace SimTK;

namespace {

// Input matrix A (4x6), stored row-major.
static const Real A[24] = {2.27, 0.28, -0.48, 1.07, -2.35, 0.62, -1.54, -1.67, -3.09, 1.22, 2.93, -7.39,
                           1.15, 0.94, 0.99,  0.79, -1.45, 1.03, -1.94, -0.78, -0.21, 0.63, 2.30, -2.57};

// Expected singular values.
static const Real X[4] = {9.9966, 3.6831, 1.3569, 0.5000};

static const double EPS = 0.001;

// ---------------------------------------------------------------------------
// Formatting helpers – attached to assertion failure messages only.
// ---------------------------------------------------------------------------

std::string FormatVector(const Vector& v, const char* label = "vector") {
    std::ostringstream oss;
    oss << label << " (size=" << v.size() << "):\n";
    for (int i = 0; i < v.size(); ++i) {
        oss << "  [" << i << "] = " << v[i] << "\n";
    }
    return oss.str();
}

std::string FormatMatrix(const Matrix& m, const char* label = "matrix") {
    std::ostringstream oss;
    oss << label << " (" << m.nrow() << "x" << m.ncol() << "):\n";
    for (int i = 0; i < m.nrow(); ++i) {
        oss << " ";
        for (int j = 0; j < m.ncol(); ++j) {
            oss << "  " << m(i, j);
        }
        oss << "\n";
    }
    return oss.str();
}

} // anonymous namespace


// ---------------------------------------------------------------------------
// Test: singular values with a custom rcond threshold
// ---------------------------------------------------------------------------
TEST(FactorSVD, SingularValuesWithRcond) {
    SCOPED_TRACE("getSingularValues with rcond = 0.01");

    Matrix a(4, 6, A);
    Vector singularValues(4);
    Vector expectedValues(4, X);
    FactorSVD svd(a, 0.01);

    svd.getSingularValues(singularValues);

    EXPECT_LT((singularValues - expectedValues).norm(), EPS)
        << FormatVector(singularValues, "computed") << FormatVector(expectedValues, "expected");
}


// ---------------------------------------------------------------------------
// Test: singular values and vectors with default rcond
// ---------------------------------------------------------------------------
TEST(FactorSVD, SingularValuesAndVectors) {
    SCOPED_TRACE("getSingularValuesAndVectors with default rcond");

    Matrix a(4, 6, A);
    Vector singularValues(4);
    Vector expectedValues(4, X);
    Matrix leftVectors;
    Matrix rightVectors;
    FactorSVD svd(a);

    svd.getSingularValuesAndVectors(singularValues, leftVectors, rightVectors);

    EXPECT_LT((singularValues - expectedValues).norm(), EPS)
        << FormatVector(singularValues, "computed singular values")
        << FormatVector(expectedValues, "expected singular values");

    // Left vectors U must be column-orthonormal: U^T * U ≈ I.
    {
        SCOPED_TRACE("left vectors column-orthonormality (U^T * U ≈ I)");
        Matrix UtU = leftVectors.transpose() * leftVectors;
        for (int i = 0; i < UtU.nrow(); ++i) {
            for (int j = 0; j < UtU.ncol(); ++j) {
                EXPECT_NEAR(UtU(i, j), (i == j) ? 1.0 : 0.0, EPS)
                    << "U^T*U(" << i << "," << j << ") out of tolerance\n"
                    << FormatMatrix(leftVectors, "leftVectors");
            }
        }
    }

    // Right vectors V must be row-orthonormal: V * V^T ≈ I.
    {
        SCOPED_TRACE("right vectors row-orthonormality (V * V^T ≈ I)");
        Matrix VVt = rightVectors * rightVectors.transpose();
        for (int i = 0; i < VVt.nrow(); ++i) {
            for (int j = 0; j < VVt.ncol(); ++j) {
                EXPECT_NEAR(VVt(i, j), (i == j) ? 1.0 : 0.0, EPS)
                    << "V*V^T(" << i << "," << j << ") out of tolerance\n"
                    << FormatMatrix(rightVectors, "rightVectors");
            }
        }
    }
}


// ---------------------------------------------------------------------------
// Test: matrix inverse via SVD
// ---------------------------------------------------------------------------
TEST(FactorSVD, Inverse2x2) {
    SCOPED_TRACE("inverse of 2x2 matrix via SVD");

    Real C[4] = {1.0, 2.0, 1.0, 3.0};
    Matrix c(2, 2, C);
    FactorSVD csvd(c);
    Matrix invSVD;

    csvd.inverse(invSVD);

    // Verify A * A^-1 ≈ I.
    Matrix identity = c * invSVD;
    for (int i = 0; i < identity.nrow(); ++i) {
        for (int j = 0; j < identity.ncol(); ++j) {
            EXPECT_NEAR(identity(i, j), (i == j) ? 1.0 : 0.0, EPS)
                << "A * inv(A) (" << i << "," << j << ") out of tolerance\n"
                << FormatMatrix(c, "A") << FormatMatrix(invSVD, "inv(A)")
                << FormatMatrix(identity, "A * inv(A)");
        }
    }
}


// ---------------------------------------------------------------------------
// Test: solve with an all-zeros matrix
// ---------------------------------------------------------------------------
TEST(FactorSVD, SolveAllZerosMatrix) {
    SCOPED_TRACE("solve with all-zeros 2x2 matrix");

    Real Z[4] = {0.0, 0.0, 0.0, 0.0};
    Matrix z(2, 2, Z);
    FactorSVD zsvd(z);

    Vector_<double> xz;
    Vector_<double> bz(2, 0.0);

    zsvd.solve(bz, xz);

    EXPECT_EQ(xz.size(), 2) << "solution vector has unexpected size";
}


// ---------------------------------------------------------------------------
// Test: solve with a 0x0 (empty) matrix
// ---------------------------------------------------------------------------
TEST(FactorSVD, SolveEmptyMatrix) {
    SCOPED_TRACE("solve with empty (0x0) matrix");

    Matrix_<double> z0;
    FactorSVD z0svd(z0);

    Vector_<double> xz;
    Vector_<double> bz0(0);

    z0svd.solve(bz0, xz);

    EXPECT_EQ(xz.size(), 0) << "solution vector of empty system should be empty";
}


// ---------------------------------------------------------------------------
// Test: SVD factorisation of a 0x0 (empty) matrix
// ---------------------------------------------------------------------------
TEST(FactorSVD, SVDEmptyMatrix) {
    SCOPED_TRACE("SVD factorisation of empty (0x0) matrix");

    Matrix_<double> z0;
    FactorSVD z0fsvd(z0);

    Vector singularValues;
    Matrix leftVectors;
    Matrix rightVectors;

    z0fsvd.getSingularValuesAndVectors(singularValues, leftVectors, rightVectors);

    EXPECT_EQ(singularValues.size(), 0) << "singular values of empty matrix should be empty";
    EXPECT_EQ(leftVectors.nrow() * leftVectors.ncol(), 0)
        << FormatMatrix(leftVectors, "leftVectors (should be empty)");
    EXPECT_EQ(rightVectors.nrow() * rightVectors.ncol(), 0)
        << FormatMatrix(rightVectors, "rightVectors (should be empty)");
}