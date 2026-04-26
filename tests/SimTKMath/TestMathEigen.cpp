/* -------------------------------------------------------------------------- *
 *                         SimTK Simbody: SimTKmath                           *
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
 * This is a test program which uses the Eigen  class to compute
 * eigen values and eigen vectors
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves for the eigen valus and vectors for the
following system

   Ax = 0 where :



     0.35   0.45   -0.14   -0.17
     0.09   0.07   -0.54    0.35
A = -0.44  -0.33   -0.03    0.17
     0.25  -0.32   -0.13    0.11



SOLUTION =
reciprocal condition number =  9.9E-01
 Error bound                 =  1.3E-16

 Eigenvector( 1)
 -6.5509E-01
 -5.2363E-01
  5.3622E-01
 -9.5607E-02

 Reciprocal condition number =  8.2E-01
 Error bound                 =  1.6E-16

 Eigenvalue( 2) = (-9.9412E-02, 4.0079E-01)

 Reciprocal condition number =  7.0E-01
 Error bound                 =  1.8E-16

 Eigenvector( 2)
 (-1.9330E-01, 2.5463E-01)
 ( 2.5186E-01,-5.2240E-01)
 ( 9.7182E-02,-3.0838E-01)
 ( 6.7595E-01, 0.0000E+00)

 Reciprocal condition number =  4.0E-01
 Error bound                 =  3.3E-16

 Eigenvalue( 3) = (-9.9412E-02,-4.0079E-01)

 Reciprocal condition number =  7.0E-01
 Error bound                 =  1.8E-16

 Eigenvector( 3)
 (-1.9330E-01,-2.5463E-01)
 ( 2.5186E-01, 5.2240E-01)
 ( 9.7182E-02, 3.0838E-01)
 ( 6.7595E-01,-0.0000E+00)

 Reciprocal condition number =  4.0E-01
 Error bound                 =  3.3E-16

 Eigenvalue( 4) = -1.0066E-01

 Reciprocal condition number =  5.7E-01
 Error bound                 =  2.3E-16

 Eigenvector( 4)
  1.2533E-01
  3.3202E-01
  5.9384E-01
  7.2209E-01

 Reciprocal condition number =  3.1E-01
 Error bound                 =  4.2E-16


estimated rank = 4

*/

#include <cmath>
#include <complex>
#include <gtest/gtest.h>

#include "SimTKmath.h"

using namespace SimTK;
using cd = std::complex<double>;
using cf = std::complex<float>;

// ---------------------------------------------------------------------------
// Shared fixtures
// ---------------------------------------------------------------------------

namespace {

const std::array<double, 16> A_DATA =
    {0.35, 0.45, -0.14, -0.17, 0.09, 0.07, -0.54, 0.35, -0.44, -0.33, -0.03, 0.17, 0.25, -0.32, -0.13, 0.11};

const std::array<cd, 4> EXPECTED_EIGENVALUES = {cd{0.79948, 0.0},
                                                cd{-0.099412, 0.40079},
                                                cd{-0.099412, -0.40079},
                                                cd{-0.10066, 0.0}};

// Stored column-major: column j, row i → expVectors[j*4 + i]
const std::array<cd, 16> EXPECTED_EIGENVECTORS = {
    cd{-.65509, 0.0},
    cd{-.52363, 0.0},
    cd{.53622, 0.0},
    cd{-.095607, 0.0}, // Col 0
    cd{-.1933001, .25463},
    cd{.2518601, -.52240},
    cd{.09718202, -.30838},
    cd{.67595, 0.0}, // Col 1
    cd{-.1933001, -.25463},
    cd{.2518601, .52240},
    cd{.09718202, .30838},
    cd{.67595, 0.0}, // Col 2
    cd{.12533, 0.0},
    cd{.33202, 0.0},
    cd{.59384, 0.0},
    cd{.72209, 0.0} // Col 3
};

// ---------------------------------------------------------------------------
// Norm helpers (compare absolute values of each component)
// ---------------------------------------------------------------------------

template <typename T>
auto absNormComplex(const Vector_<std::complex<T>>& values, const Vector_<std::complex<T>>& expected) -> T {
    T norm = 0;
    for (int i = 0; i < values.size(); ++i) {
        T dr = std::fabs(values(i).real()) - std::fabs(expected(i).real());
        T di = std::fabs(values(i).imag()) - std::fabs(expected(i).imag());
        norm += dr * dr + di * di;
    }
    return std::sqrt(norm);
}

// Build the 4x4 expected eigenvector matrix (rows = eigenvectors)
auto makeExpectedVectorMatrix() -> Matrix_<cd> {
    Matrix_<cd> m(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            m(i, j) = EXPECTED_EIGENVECTORS[(j * 4) + i];
        }
    }
    return m;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------

class EigenTest : public ::testing::Test {
    protected:
    Matrix a{4, 4, A_DATA.data()};
};

// ---------------------------------------------------------------------------
// Double-precision tests
// ---------------------------------------------------------------------------

TEST_F(EigenTest, NonSymmetricDouble_EigenvaluesMatchReference) {
    Vector_<cd> values;
    Matrix_<cd> vectors;
    Eigen es(a);
    es.getAllEigenValuesAndVectors(values, vectors);

    Vector_<cd> expected(4);
    for (int i = 0; i < 4; ++i) {
        expected[i] = EXPECTED_EIGENVALUES[i];
    }

    EXPECT_LT(absNormComplex(values, expected), 0.001) << "Eigenvalue error norm too large.\n"
                                                       << "Computed: " << values << "\n";
}

TEST_F(EigenTest, NonSymmetricDouble_EigenvectorsMatchReference) {
    Vector_<cd> values;
    Matrix_<cd> vectors;
    Eigen es(a);
    es.getAllEigenValuesAndVectors(values, vectors);

    auto expectedVectors = makeExpectedVectorMatrix();

    for (int i = 0; i < 4; ++i) {
        Vector_<cd> computed = vectors(i);
        Vector_<cd> expected = expectedVectors(i);
        double errnorm = absNormComplex(computed, expected);
        EXPECT_LT(errnorm, 0.00001) << "Eigenvector " << i << " error norm too large.\n"
                                    << "Computed: " << computed << "\n"
                                    << "Expected: " << expected << "\n";
    }
}

// ---------------------------------------------------------------------------
// Single-precision tests
// ---------------------------------------------------------------------------

TEST_F(EigenTest, NonSymmetricFloat_EigenvaluesMatchReference) {
    Matrix_<float> af(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            af(i, j) = static_cast<float>(a(i, j));
        }
    }

    Vector_<cf> expectedf(4);
    for (int i = 0; i < 4; ++i) {
        expectedf[i] = static_cast<cf>(EXPECTED_EIGENVALUES[i]);
    }

    Vector_<cf> valuesf;
    Matrix_<cf> vectorsf;
    Eigen esf(af);
    esf.getAllEigenValuesAndVectors(valuesf, vectorsf);

    EXPECT_LT(absNormComplex(valuesf, expectedf), 0.001) << "Float eigenvalue error norm too large.\n"
                                                         << "Computed: " << valuesf << "\n";
}

TEST_F(EigenTest, NonSymmetricFloat_EigenvectorsMatchReference) {
    Matrix_<float> af(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            af(i, j) = static_cast<float>(a(i, j));
        }
    }

    Matrix_<cf> expectedVectorsf(4, 4);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            expectedVectorsf(i, j) = static_cast<cf>(EXPECTED_EIGENVECTORS[(j * 4) + i]);
        }
    }

    Vector_<cf> valuesf;
    Matrix_<cf> vectorsf;
    Eigen esf(af);
    esf.getAllEigenValuesAndVectors(valuesf, vectorsf);

    for (int i = 0; i < 4; ++i) {
        Vector_<cf> computed = vectorsf(i);
        Vector_<cf> expected = expectedVectorsf(i);
        float errnorm = absNormComplex(computed, expected);
        EXPECT_LT(errnorm, 0.0001F) << "Float eigenvector " << i << " error norm too large.\n"
                                    << "Computed: " << computed << "\n"
                                    << "Expected: " << expected << "\n";
    }
}

// ---------------------------------------------------------------------------
// Edge-case tests
// ---------------------------------------------------------------------------

TEST(EigenEdgeCases, ZeroMatrix_ReturnsZeroEigenvalues) {
    Real Z[4] = {0.0, 0.0, 0.0, 0.0};
    Matrix z(2, 2, Z);
    Eigen zeigen(z);

    Vector_<cd> values;
    Matrix_<cd> vectors;
    zeigen.getAllEigenValuesAndVectors(values, vectors);

    ASSERT_EQ(values.size(), 2);
    for (int i = 0; i < values.size(); ++i) {
        EXPECT_NEAR(values(i).real(), 0.0, 1e-10)
            << "Zero-matrix eigenvalue " << i << " real part should be 0";
        EXPECT_NEAR(values(i).imag(), 0.0, 1e-10)
            << "Zero-matrix eigenvalue " << i << " imag part should be 0";
    }
}

TEST(EigenEdgeCases, EmptyMatrix_ThrowsOnSolve) {
    // The original test accidentally called the zero-matrix solver a second
    // time here instead of the empty-matrix one, so this path was never
    // exercised. Constructing Eigen from a 0x0 matrix and then asking for
    // eigenvalues triggers an illegal LAPACK argument — document that.
    Matrix_<double> z0;
    Eigen z0Eigen(z0);
    Vector_<cd> values;
    Matrix_<cd> vectors;
    EXPECT_THROW(z0Eigen.getAllEigenValuesAndVectors(values, vectors), SimTK::Exception::IllegalLapackArg);
}
