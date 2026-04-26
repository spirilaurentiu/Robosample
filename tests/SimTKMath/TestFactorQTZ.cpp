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
 * This is a test program which uses the FactorQTZ  class to do an QTZ
 * factorization on a system of linear equations and then use the
 * factored QTZ matrix to solve a find a least squares solution
 * for a particular right hand side
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves the least squares problem:

Ax = B,

where A is the general matrix


     -0.09   0.14  -0.46    0.68   1.29        7.4
     -1.56   0.20   0.29    1.09   0.51        4.2
A =  -1.48  -0.43   0.89   -0.71  -0.96    B= -8.3
     -1.09   0.84   0.77    2.11  -1.27        1.8
      0.08   0.55  -1.13    0.14   1.74        8.6
     -1.59  -0.72   1.06    1.24   0.34        2.1

The default tolerance of 0.01 is used to determine the effective rank of A


SOLUTION =
0.6344     0.9699    -1.4402     3.3678     3.3992

estimated rank = 4

*/

#include <gtest/gtest.h>
#include <sstream>

#include "SimTKmath.h"

using namespace SimTK;

namespace {

Real A[30] = {-0.09, 0.14,  -0.46, 0.68,  1.29,  -1.56, 0.20,  0.29, 1.09, 0.51,
              -1.48, -0.43, 0.89,  -0.71, -0.96, -1.09, 0.84,  0.77, 2.11, -1.27,
              0.08,  0.55,  -1.13, 0.14,  1.74,  -1.59, -0.72, 1.06, 1.24, 0.34};

Real B[6] = {7.4, 4.2, -8.3, 1.8, 8.6, 2.1};
Real X[5] = {0.6344, 0.9699, -1.4402, 3.3678, 3.3992};

} // namespace

TEST(FactorQTZ, OverdeterminedDoubleDefaultRcond) {
    Matrix a(6, 5, A);
    FactorQTZ qtz;

    qtz.factor(a);

    SCOPED_TRACE("Checking estimated rank with default rcond");
    EXPECT_EQ(qtz.getRank(), 5) << "Estimated rank with default rcond: " << qtz.getRank() << "\n";
}

TEST(FactorQTZ, OverdeterminedDouble) {
    Matrix a(6, 5, A);
    Vector b(6, B);
    Vector x_right(5, X);
    Vector x;

    FactorQTZ qtz;
    qtz.factor(a, 0.01);
    qtz.solve(b, x);

    {
        SCOPED_TRACE("Checking estimated rank with rcond = 0.01");
        EXPECT_EQ(qtz.getRank(), 4) << "Estimated rank with rcond = 0.01: " << qtz.getRank() << "\n";
    }

    double errnorm = (x - x_right).norm();
    {
        SCOPED_TRACE("Overdetermined Double SOLUTION");
        EXPECT_LT(errnorm, 0.001) << "Overdetermined Double SOLUTION: " << x << "  errnorm=" << errnorm
                                  << "\n";
    }
}

TEST(FactorQTZ, OverdeterminedDoubleCopyConstructor) {
    Matrix a(6, 5, A);
    Vector b(6, B);
    Vector x_right(5, X);
    Vector xc;

    FactorQTZ qtz;
    qtz.factor(a, 0.01);

    FactorQTZ qtzCopy(qtz);
    qtzCopy.solve(b, xc);

    double errnorm = (xc - x_right).norm();
    {
        SCOPED_TRACE("Copy constructor SOLUTION");
        EXPECT_LT(errnorm, 0.001) << "Copy constructor SOLUTION: " << xc << "  errnorm=" << errnorm << "\n";
    }
}

TEST(FactorQTZ, OverdeterminedDoubleCopyAssign) {
    Matrix a(6, 5, A);
    Vector b(6, B);
    Vector x_right(5, X);
    Vector xa;

    FactorQTZ qtz;
    qtz.factor(a, 0.01);

    FactorQTZ qtzAssign = qtz;
    qtzAssign.solve(b, xa);

    double errnorm = (xa - x_right).norm();
    {
        SCOPED_TRACE("Copy assign SOLUTION");
        EXPECT_LT(errnorm, 0.001) << "Copy assign SOLUTION: " << xa << "  errnorm=" << errnorm << "\n";
    }
}

TEST(FactorQTZ, OverdeterminedFloat) {
    Matrix a(6, 5, A);
    Vector b(6, B);
    Vector x_right(5, X);

    Matrix_<float> af(6, 5);
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < 5; ++j) {
            af(i, j) = (float)a(i, j);
        }
    }

    Vector_<float> bf(6);
    for (int i = 0; i < 6; ++i) {
        bf[i] = (float)b[i];
    }

    Vector_<float> xf_right(5);
    for (int i = 0; i < 5; ++i) {
        xf_right[i] = (float)x_right[i];
    }

    Vector_<float> xf;

    FactorQTZ qtz;
    qtz.factor(af, (float)0.01);
    qtz.solve(bf, xf);

    float errnorm = (xf - xf_right).norm();
    {
        SCOPED_TRACE("Overdetermined Float SOLUTION");
        EXPECT_LT(errnorm, 0.001f) << "Overdetermined Float SOLUTION: " << xf << "  errnorm=" << errnorm
                                   << "\n";
    }
}

TEST(FactorQTZ, UnderdeterminedDouble) {
    Real Au[12] = {2, 5, 3, 4, 7, 1, 3, 5, 4, 3, 6, 2};
    Real Bu[3] = {3, 1, 6};
    Real Xu[4] = {-0.0376844, 0.350628, 0.986164, -0.409066};

    Matrix au(3, 4, Au);
    Vector bu(3, Bu);
    Vector xu_right(4, Xu);
    Vector xu;

    FactorQTZ qtzu(au);
    qtzu.solve(bu, xu);

    double errnorm = (xu - xu_right).norm();
    {
        SCOPED_TRACE("Underdetermined Double SOLUTION");
        EXPECT_LT(errnorm, 0.001) << "Underdetermined Double SOLUTION: " << xu << "  errnorm=" << errnorm
                                  << "\n";
    }
}

TEST(FactorQTZ, UnderdeterminedDoubleMultipleRhs) {
    Real Au[12] = {2, 5, 3, 4, 7, 1, 3, 5, 4, 3, 6, 2};
    Real Bu[3] = {3, 1, 6};

    Matrix au(3, 4, Au);
    Vector bu(3, Bu);

    Matrix bu2(3, 2);
    bu2(0) = bu;
    bu2(1) = 2 * bu;
    Matrix xu2;

    FactorQTZ qtzu(au);
    qtzu.solve(bu2, xu2);

    SCOPED_TRACE("Underdetermined Double multiple RHS");
    EXPECT_EQ(xu2.nrow(), 4) << "Multiple RHS solution (double):\n" << xu2 << "\n";
    EXPECT_EQ(xu2.ncol(), 2) << "Multiple RHS solution (double):\n" << xu2 << "\n";
}

TEST(FactorQTZ, UnderdeterminedFloat) {
    Real Au[12] = {2, 5, 3, 4, 7, 1, 3, 5, 4, 3, 6, 2};
    Real Bu[3] = {3, 1, 6};
    Real Xu[4] = {-0.0376844, 0.350628, 0.986164, -0.409066};

    Matrix au(3, 4, Au);
    Vector bu(3, Bu);
    Vector xu_right(4, Xu);

    Matrix_<float> afu(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            afu(i, j) = (float)au(i, j);
        }
    }

    Vector_<float> bfu(3);
    for (int i = 0; i < 3; ++i) {
        bfu[i] = (float)bu[i];
    }

    Vector_<float> xfu_right(4);
    for (int i = 0; i < 4; ++i) {
        xfu_right[i] = (float)xu_right[i];
    }

    Vector_<float> xfu;

    FactorQTZ qtzfu(afu);
    qtzfu.solve(bfu, xfu);

    float errnorm = (xfu - xfu_right).norm();
    {
        SCOPED_TRACE("Underdetermined Float SOLUTION");
        EXPECT_LT(errnorm, 0.001f) << "Underdetermined Float SOLUTION: " << xfu << "  errnorm=" << errnorm
                                   << "\n";
    }
}

TEST(FactorQTZ, UnderdeterminedFloatMultipleRhs) {
    Real Au[12] = {2, 5, 3, 4, 7, 1, 3, 5, 4, 3, 6, 2};
    Real Bu[3] = {3, 1, 6};

    Matrix au(3, 4, Au);
    Vector bu(3, Bu);

    Matrix_<float> afu(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            afu(i, j) = (float)au(i, j);
        }
    }

    Vector_<float> bfu(3);
    for (int i = 0; i < 3; ++i) {
        bfu[i] = (float)bu[i];
    }

    Matrix_<float> bfu2(3, 2);
    bfu2(0) = bfu;
    bfu2(1) = 2 * bfu;
    Matrix_<float> xfu2;

    FactorQTZ qtzfu(afu);
    qtzfu.solve(bfu2, xfu2);

    SCOPED_TRACE("Underdetermined Float multiple RHS");
    EXPECT_EQ(xfu2.nrow(), 4) << "Multiple RHS solution (float):\n" << xfu2 << "\n";
    EXPECT_EQ(xfu2.ncol(), 2) << "Multiple RHS solution (float):\n" << xfu2 << "\n";
}

TEST(FactorQTZ, Inverse2x2) {
    Real C[4] = {1.0, 2.0, 1.0, 3.0};
    Matrix c(2, 2, C);

    FactorQTZ cqtz(c);
    Matrix invQTZ;
    cqtz.inverse(invQTZ);

    SCOPED_TRACE("FactorQTZ.inverse 2x2");
    ASSERT_EQ(invQTZ.nrow(), 2) << "FactorQTZ.inverse result:\n" << invQTZ[0] << "\n" << invQTZ[1] << "\n";
    ASSERT_EQ(invQTZ.ncol(), 2) << "FactorQTZ.inverse result:\n" << invQTZ[0] << "\n" << invQTZ[1] << "\n";

    // Verify A * A^-1 = I
    Matrix identity = c * invQTZ;
    EXPECT_NEAR(identity(0, 0), 1.0, 1e-10) << "FactorQTZ.inverse:\n"
                                            << invQTZ[0] << "\n"
                                            << invQTZ[1] << "\n";
    EXPECT_NEAR(identity(0, 1), 0.0, 1e-10) << "FactorQTZ.inverse:\n"
                                            << invQTZ[0] << "\n"
                                            << invQTZ[1] << "\n";
    EXPECT_NEAR(identity(1, 0), 0.0, 1e-10) << "FactorQTZ.inverse:\n"
                                            << invQTZ[0] << "\n"
                                            << invQTZ[1] << "\n";
    EXPECT_NEAR(identity(1, 1), 1.0, 1e-10) << "FactorQTZ.inverse:\n"
                                            << invQTZ[0] << "\n"
                                            << invQTZ[1] << "\n";
}

TEST(FactorQTZ, SolveAllZeroMatrix) {
    Real Z[4] = {0.0, 0.0, 0.0, 0.0};
    Matrix z(2, 2, Z);

    FactorQTZ zqtz(z);
    Vector_<double> xz;
    Vector_<double> bz(2);
    bz(0) = bz(1) = 0.0;

    SCOPED_TRACE("Solve with all-zero matrix");
    EXPECT_NO_THROW(zqtz.solve(bz, xz)) << "Expected solve with zero matrix to succeed\n";

    for (int i = 0; i < xz.size(); ++i) {
        EXPECT_NEAR(xz(i), 0.0, 1e-10) << "xz(" << i << ") = " << xz(i) << "\n";
    }
}

TEST(FactorQTZ, SolveNullMatrixThrows) {
    SCOPED_TRACE("Solve with null (0x0) matrix should throw");
    EXPECT_THROW(
        {
            Matrix_<double> z0;
            FactorQTZ z0qtz(z0);
            Vector_<double> bz0(0);
            Vector_<double> xz;
            z0qtz.solve(bz0, xz);
        },
        std::exception);
}
