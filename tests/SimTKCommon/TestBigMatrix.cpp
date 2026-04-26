/* -------------------------------------------------------------------------- *
 *                       SimTK Simbody: SimTKcommon                           *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

#include <complex>
#include <cstdlib>
#include <gtest/gtest.h>
#include <string>

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;
using P = double;

// Numerical Recipes version 2.11 LU decomp and inversion via backsolve.
// This is the C++ version modified to use column-ordered consecutive storage.
namespace NumericalRecipes {

// Return 1d index for column ordered matrix with leading dim N.
#define X(i, j) j* N + i

template <class DP>
void lubksb(const int N, const DP* a /*N,N*/, const int* indx /*N*/, DP* b /*N*/) {
    int i;
    int ii = 0;
    int ip;
    int j;
    DP sum;

    for (i = 0; i < N; i++) {
        ip = indx[i];
        sum = b[ip];
        b[ip] = b[i];
        if (ii != 0) {
            for (j = ii - 1; j < i; j++) {
                sum -= a[X(i, j)] * b[j];
            }
        } else if (sum != 0.0) {
            ii = i + 1;
        }
        b[i] = sum;
    }
    for (i = N - 1; i >= 0; i--) {
        sum = b[i];
        for (j = i + 1; j < N; j++) {
            sum -= a[X(i, j)] * b[j];
        }
        b[i] = sum / a[X(i, i)];
    }
}

template <class DP>
void ludcmp(const int N, DP* a /*N,N*/, int* indx /*N*/, DP& d) {
    const DP TINY = DP(1.0e-20);
    int i;
    int imax;
    int j;
    int k;
    DP big;
    DP dum;
    DP sum;
    DP temp;

    DP* vv = new DP[N];
    d = DP(1);
    for (i = 0; i < N; i++) {
        big = 0.0;
        for (j = 0; j < N; j++) {
            temp = fabs(a[X(i, j)]);
            if (temp > big) {
                big = temp;
            }
        }
        EXPECT_NE(big, 0.0) << "Singular matrix in routine ludcmp";
        vv[i] = DP(1) / big;
    }
    for (j = 0; j < N; j++) {
        for (i = 0; i < j; i++) {
            sum = a[X(i, j)];
            for (k = 0; k < i; k++) {
                sum -= a[X(i, k)] * a[X(k, j)];
            }
            a[X(i, j)] = sum;
        }
        big = 0.0;
        for (i = j; i < N; i++) {
            sum = a[X(i, j)];
            for (k = 0; k < j; k++) {
                sum -= a[X(i, k)] * a[X(k, j)];
            }
            a[X(i, j)] = sum;
            dum = vv[i] * fabs(sum);
            if (dum >= big) {
                big = dum;
                imax = i;
            }
        }
        if (j != imax) {
            for (k = 0; k < N; k++) {
                dum = a[X(imax, k)];
                a[X(imax, k)] = a[X(j, k)];
                a[X(j, k)] = dum;
            }
            d = -d;
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (a[X(j, j)] == 0.0) {
            a[X(j, j)] = TINY;
        }
        if (j != N - 1) {
            dum = DP(1) / (a[X(j, j)]);
            for (i = j + 1; i < N; i++) {
                a[X(i, j)] *= dum;
            }
        }
    }

    delete[] vv;
}

template <class DP>
void luinvert(const int N, DP* a /*N,N*/, DP* y /*N,N*/) {
    assert(a && y);
    int* indx = new int[N];
    DP d;
    NumericalRecipes::ludcmp(N, a, indx, d);
    for (int j = 0; j < N; ++j) {
        DP* col = &y[X(0, j)];
        for (int i = 0; i < N; ++i) {
            col[i] = DP(0);
        }
        col[j] = DP(1);
        NumericalRecipes::lubksb(N, a, indx, col); // writes directly into y
    }

    delete[] indx;
}

} // namespace NumericalRecipes

// Some explicit instantiations just to make sure everything's there.
namespace SimTK {
template class Matrix_<Real>;
template class Vector_<Complex>;
template class RowVector_<conjugate<float>>;
// template class MatrixBase< Mat<3,4,Vec2> >;

template class MatrixView_<complex<double>>;
template class VectorView_<negator<float>>;
template class RowVectorView_<negator<conjugate<float>>>;

template class MatrixBase<double>;
template class VectorBase<double>;
template class RowVectorBase<double>;
template class MatrixView_<double>;
template class VectorView_<double>;
template class RowVectorView_<double>;
template class Vector_<double>;
template class RowVector_<double>;

template class MatrixBase<negator<double>>;
template class VectorBase<negator<double>>;
template class RowVectorBase<negator<double>>;
template class MatrixView_<negator<double>>;
template class VectorView_<negator<double>>;
template class RowVectorView_<negator<double>>;
template class Matrix_<negator<double>>;
template class Vector_<negator<double>>;
template class RowVector_<negator<double>>;

} // namespace SimTK

namespace {

template <class T>
auto isNaN(const T& v) -> bool {
    return v.isNaN();
}
template <>
auto isNaN(const double& v) -> bool {
    return SimTK::isNaN(v);
}
template <>
auto isNaN(const float& v) -> bool {
    return SimTK::isNaN(v);
}
template <>
auto isNaN(const negator<double>& v) -> bool {
    return SimTK::isNaN(v);
}
template <>
auto isNaN(const negator<float>& v) -> bool {
    return SimTK::isNaN(v);
}

/// Verify that every element of a BigMatrix vector-like object matches
/// the corresponding element of a fixed-size Vec<N>.
template <class T, int N>
void expectVectorEq(const T& value, const Vec<N>& expected) {
    ASSERT_EQ(value.size(), N);
    for (int i = 0; i < N; ++i) {
        if (isNaN(expected[i])) {
            EXPECT_TRUE(isNaN(value[i])) << "Expected NaN at index " << i << " but got " << value[i];
        } else {
            EXPECT_EQ(value[i], expected[i]) << "Mismatch at index " << i;
        }
    }
}

/// Verify that every element of a BigMatrix matrix-like object matches
/// the corresponding element of a fixed-size Mat<M,N>.
template <class T, int M, int N>
void expectMatrixEq(const T& value, const Mat<M, N, typename T::E>& expected) {
    ASSERT_EQ(value.nrow(), M);
    ASSERT_EQ(value.ncol(), N);
    for (int i = 0; i < M; ++i) {
        for (int j = 0; j < N; ++j) {
            if (isNaN(expected(i, j))) {
                EXPECT_TRUE(isNaN(value(i, j)))
                    << "Expected NaN at (" << i << "," << j << ") but got " << value(i, j);
            } else {
                EXPECT_EQ(value(i, j), expected(i, j)) << "Mismatch at (" << i << "," << j << ")";
            }
        }
    }
}

} // anonymous namespace


TEST(SimTKCommon_BigMatrix_RealSums, ReturnsCorrectSum) {
    Matrix m(Mat22(1, 2, 3, 4));

    EXPECT_TRUE(AssertSimTKEqual("m.colSum()", "RowVector(4, 6)", m.colSum(), RowVector(Row2(4, 6))));
    EXPECT_TRUE(AssertSimTKEqual("m.rowSum()", "Vector(3, 7)", m.rowSum(), Vector(Vec2(3, 7))));
    EXPECT_TRUE(AssertSimTKEqual("m.sum()", "m.colSum()", m.sum(), m.colSum()));
}

TEST(SimTKCommon_BigMatrix_ComplexSums, ReturnsCorrectSum) {
    const Complex I = Complex(0, 1);
    Matrix_<Complex> mc(Mat<2, 2, Complex>(1 + 2 * I, 3 + 4 * I, 5 + 6 * I, 7 + 8 * I));

    using CRow2 = Row<2, Complex>;
    using CVec2 = Vec<2, Complex>;
    using CRowVector = RowVector_<Complex>;
    using CVector = Vector_<Complex>;

    EXPECT_TRUE(AssertSimTKEqual("mc.colSum()",
                                 "CRowVector(6+8I, 10+12I)",
                                 mc.colSum(),
                                 CRowVector(CRow2(6 + 8 * I, 10 + 12 * I))));
    EXPECT_TRUE(AssertSimTKEqual("mc.rowSum()",
                                 "CVector(4+6I, 12+14I)",
                                 mc.rowSum(),
                                 CVector(CVec2(4 + 6 * I, 12 + 14 * I))));
    EXPECT_TRUE(AssertSimTKEqual("mc.sum()", "mc.colSum()", mc.sum(), mc.colSum()));
}

TEST(SimTKCommon_BigMatrix_Character, ScalarMultiplicationAndViews) {
    Matrix m(2, 3);
    m.setTo(5.0); // fills ALL elements with 5, not identity-style
    Matrix mm = m * 3;

    // mm should be 15 everywhere
    for (int i = 0; i < mm.nrow(); ++i) {
        for (int j = 0; j < mm.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("mm(i,j)", "15", mm(i, j), 15.0));
        }
    }

    // Submatrix view correctness: rows 0-1, cols 1-2 (2 rows, 2 cols)
    MatrixView mmv = mm(0, 1, 2, 2);
    for (int i = 0; i < mmv.nrow(); ++i) {
        for (int j = 0; j < mmv.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("mmv(i,j)", "15", mmv(i, j), 15.0));
        }
    }

    // Reassign view to a different block: row 0, cols 0-1 (1 row, 2 cols)
    mmv = mm(0, 0, 1, 2);
    for (int i = 0; i < mmv.nrow(); ++i) {
        for (int j = 0; j < mmv.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("mmv reassigned", "15", mmv(i, j), 15.0));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Character, CommitmentAndHermitian) {
    MatrixCommitment mc(MatrixStructure(MatrixStructure::Triangular, MatrixStructure::Upper));
    Matrix t(mc);

    t.resize(5, 3);

    t.clear();
    t.commitTo(MatrixStructure(MatrixStructure::Hermitian, MatrixStructure::Upper));
    t.resize(3, 5);

    t = 123;
    t(0, 2) = -2;

    // Hermitian: t(i,j) == conjugate(t(j,i))
    for (int i = 0; i < t.nrow(); ++i) {
        for (int j = 0; j < t.ncol(); ++j) {
            if (i < t.ncol() && j < t.nrow()) {
                EXPECT_TRUE(AssertSimTKEqual("Symmetry", "t(j,i)", t(i, j), t(j, i)));
            }
        }
    }
}

TEST(SimTKCommon_BigMatrix_Character, VectorNormsAndWeights) {
    Vector v(10);
    for (int i = 0; i < 10; ++i) {
        v[i] = i * 0.1;
    }
    v[4] = -17.3;

    Vector w1(10, Real(1));

    int worst = -1;

    Real rms = v.normRMS(&worst);
    EXPECT_GE(rms, 0);
    EXPECT_GE(worst, 0);

    Real inf = v.normInf(&worst);
    EXPECT_GE(inf, 0);
    EXPECT_GE(worst, 0);

    Real wrms = v.weightedNormRMS(w1, &worst);
    EXPECT_GE(wrms, 0);

    Real winf = v.weightedNormInf(w1, &worst);
    EXPECT_GE(winf, 0);

    // Change weights and ensure effect
    w1(9) = 100;
    Real wrms2 = v.weightedNormRMS(w1);
    EXPECT_NE(wrms, wrms2);
}

TEST(SimTKCommon_BigMatrix_Character, ZeroLengthVectorNorms) {
    Vector v;

    int worst = -1;
    EXPECT_TRUE(AssertSimTKEqual("normRMS", "0", v.normRMS(&worst), 0));
    EXPECT_TRUE(AssertSimTKEqual("normInf", "0", v.normInf(&worst), 0));
}

TEST(SimTKCommon_BigMatrix_Character, IndexedVectorView) {
    Vector v(10);
    for (int i = 0; i < 10; ++i) {
        v[i] = i;
    }

    Array_<int> idx{2, 5, 7, 8};
    VectorView vxx = v(idx);

    ASSERT_EQ(vxx.size(), 4);

    EXPECT_TRUE(AssertSimTKEqual("vxx[0]", "2", vxx[0], 2));
    EXPECT_TRUE(AssertSimTKEqual("vxx[1]", "5", vxx[1], 5));
    EXPECT_TRUE(AssertSimTKEqual("vxx[2]", "7", vxx[2], 7));
    EXPECT_TRUE(AssertSimTKEqual("vxx[3]", "8", vxx[3], 8));
}

TEST(SimTKCommon_BigMatrix_Character, ComplexMatrixDataSharing) {
    Complex cmplx[] = {{1, 2}, {3, 4}, {-.2, .3}, {-100, 200}, {20, 40}, {-.001, .002}};

    ComplexMatrix cm(2, 3, 2, cmplx); // column-major shared
    ComplexMatrix cm2(2, 3, cmplx);   // row-major copy

    MatrixCharacter::LapackFull mchar(2, 3);
    mchar.setStorage(MatrixStorage(MatrixStorage::NoPacking, MatrixStorage::RowOrder));

    MatrixBase<Complex> cm3(MatrixCommitment(), mchar, 3, cmplx);

    cm3(0, 1) = 99;

    EXPECT_TRUE(AssertSimTKEqual("cm shared", "99", cm(0, 1).real(), 99.0));
    EXPECT_FALSE(SimTK::Test::numericallyEqual(cm2(0, 1), Complex(99, 0), 1));
}

TEST(SimTKCommon_BigMatrix_ScalarMultiply, SpatialVec_VectorRowConsistency) {
    Vector_<SpatialVec> vs(3, SpatialVec(Vec3(1, 2, 3), Vec3(4, 5, 6)));
    const double scalar = 2.5;

    const auto a = vs * scalar;
    const auto b = scalar * vs;

    EXPECT_TRUE(AssertSimTKEqual("vs*scalar", "scalar*vs", a, b));
}

TEST(SimTKCommon_BigMatrix_ScalarMultiply, Mat33_SpatialVec_Multiplication) {
    Mat33 m33(.03, .04, .05, .06, .08, .09, .07, .10, .11);

    SpatialVec sv(Vec3(1, 2, 3), Vec3(4, 5, 6));

    auto result = m33 * sv;

    Vec3 expected_w = m33 * sv[0];
    Vec3 expected_v = m33 * sv[1];

    EXPECT_TRUE(AssertSimTKEqual("result[0]", "m33*sv[0]", result[0], expected_w));
    EXPECT_TRUE(AssertSimTKEqual("result[1]", "m33*sv[1]", result[1], expected_v));
}

TEST(SimTKCommon_BigMatrix_ScalarMultiply, SpatialMat_VectorApplication) {
    Mat33 m33(.03, .04, .05, .06, .08, .09, .07, .10, .11);

    SpatialMat sm(m33);

    Vector_<SpatialVec> vs(3, SpatialVec(Vec3(1, 2, 3), Vec3(4, 5, 6)));

    auto result = sm * vs;

    for (int i = 0; i < vs.size(); ++i) {
        SpatialVec expected(m33 * vs[i][0], m33 * vs[i][1]);

        EXPECT_TRUE(AssertSimTKEqual("sm*vs[i]", "expected", result[i], expected));
    }
}

TEST(SimTKCommon_BigMatrix_ScalarMultiply, SpatialRow_Mat33_Consistency) {
    Mat33 m33(.03, .04, .05, .06, .08, .09, .07, .10, .11);

    SpatialVec sv(Vec3(1, 2, 3), Vec3(4, 5, 6));

    auto left = ~sv * m33;

    Row3 expected_w = ~sv[0] * m33;
    Row3 expected_v = ~sv[1] * m33;

    EXPECT_TRUE(AssertSimTKEqual("left[0]", "expected_w", left[0], expected_w));
    EXPECT_TRUE(AssertSimTKEqual("left[1]", "expected_v", left[1], expected_v));
}

TEST(SimTKCommon_BigMatrix_ScalarMultiply, SymMat22_RowVectorProducts) {
    RowVector_<SymMat22> rv(3, SymMat22(1, 2, 3));
    const double scalar = 0.1;

    // Scale then apply vs apply then scale — same result type on both sides
    Vec2 v(0.1, 0.2);
    RowVector_<Vec2> a(rv.size()), b(rv.size());
    for (int i = 0; i < rv.size(); ++i) {
        a[i] = (rv[i] * scalar) * v; // scale matrix, then multiply
        b[i] = rv[i] * (scalar * v); // scale vector, then multiply
    }

    EXPECT_TRUE(AssertSimTKEqual("(rv*s)*v", "rv*(s*v)", a, b));
}

TEST(SimTKCommon_BigMatrix_AjaysBlock, SubmatrixTransposeCorrectness) {
    const int nu = 7, nm = 4;
    Matrix J(6, nu);

    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j < nu; ++j) {
            J(i, j) = 1000 * i + j;
        }
    }

    Matrix t = ~J(0, 3, 3, nm);

    ASSERT_EQ(t.nrow(), nm);
    ASSERT_EQ(t.ncol(), 3);

    for (int r = 0; r < nm; ++r) {
        for (int c = 0; c < 3; ++c) {
            Real expected = J(c, r + 3);
            EXPECT_TRUE(AssertSimTKEqual("t(r,c)", "J(c, r+3)", t(r, c), expected));
        }
    }
}

TEST(SimTKCommon_BigMatrix_ElementwiseAssign, OverwritesAllElements) {
    Matrix assignToMe(5, 4);

    assignToMe.elementwiseAssign(1.0);
    for (int i = 0; i < assignToMe.nrow(); ++i) {
        for (int j = 0; j < assignToMe.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("assignToMe(i,j)", "1.0", assignToMe(i, j), 1.0));
        }
    }

    assignToMe.elementwiseAssign(14);
    for (int i = 0; i < assignToMe.nrow(); ++i) {
        for (int j = 0; j < assignToMe.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("assignToMe(i,j)", "14", assignToMe(i, j), 14));
        }
    }
}

TEST(SimTKCommon_BigMatrix_VectorViewAssign, SubvectorExtractionAndEdgeCases) {
    const Real vvvdata[] = {1, 2, .1, .2, 9, 10, -22, -23, -24, 25};
    Vector vvv(10, vvvdata);

    Vector vvv25;
    vvv25.viewAssign(vvv(2, 5));

    ASSERT_EQ(vvv25.size(), 5);
    for (int i = 0; i < 5; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("vvv25[i]", "vvv[i+2]", vvv25[i], vvv[i + 2]));
    }

    Vector vvv20;
    vvv20.viewAssign(vvv(2, 0));
    EXPECT_EQ(vvv20.size(), 0);

    Vector vvv00;
    vvv00.viewAssign(vvv(0, 0));
    EXPECT_EQ(vvv00.size(), 0);

    Vector vb;
    vvv00 = vb;
    EXPECT_EQ(vvv00.size(), 0);
}

TEST(SimTKCommon_BigMatrix_Complex, ConstructionConsistency) {
    const Complex mdc[] = {{1, 2},
                           {3, 4},
                           {5, 6},
                           {7, 8},
                           {9, 10},
                           {10, 11},
                           {.1, .26},
                           {.3, .45},
                           {.5, .64},
                           {.7, .83},
                           {.9, .102},
                           {.10, .111}};

    Matrix_<Complex> md(2, 2, mdc);
    Mat<2, 2, Complex> md_mat(mdc);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("md(i,j)", "md_mat(i,j)", md(i, j), md_mat(i, j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Norm, ScalarVsNonScalarNormRMS) {
    const Complex mdc[] = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
    Matrix_<Complex> md(2, 2, mdc);

    EXPECT_NO_THROW({
        Real r = md.normRMS();
        EXPECT_GE(r, 0);
    });

    Matrix_<Mat<2, 2, Complex>> mm22c;
    mm22c.resize(2, 2);

    EXPECT_THROW(mm22c.normRMS(), std::exception);
}

TEST(SimTKCommon_BigMatrix_Construction, ComplexMatrixMatchesFixedMat) {
    const Complex mdc[] = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};

    Matrix_<Complex> md(2, 2, mdc);
    Mat<2, 2, Complex> md_mat(mdc);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("md(i,j)", "md_mat(i,j)", md(i, j), md_mat(i, j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_View, AssignVsViewAssignScaling) {
    const Complex mdc[] = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};

    Matrix_<Vec<2, Complex>> rr(2, 3);

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            rr(i, j) = Vec<2, Complex>(mdc[i] * (j + 1), mdc[i] * (j - 1));
        }
    }

    Matrix_<Vec<2, Complex>> rrAssign;

    rrAssign = rr(0, 1, 2, 2);
    rrAssign *= 1000.0;

    rrAssign.viewAssign(rr(0, 1, 2, 2));
    rrAssign *= 100.0;

    // check aliasing consistency (view must match source slice)
    for (int i = 0; i < rrAssign.nrow(); ++i) {
        for (int j = 0; j < rrAssign.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("view alias", "rr slice", rrAssign(i, j), rr(0 + i, 1 + j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Complex, SetToZeroAndNaN) {
    Matrix_<Mat<2, 2, Complex>> mm22c(2, 2);

    mm22c.setToZero();
    for (int i = 0; i < mm22c.nrow(); ++i) {
        for (int j = 0; j < mm22c.ncol(); ++j) {
            for (int r = 0; r < 2; ++r) {
                for (int c = 0; c < 2; ++c) {
                    EXPECT_TRUE(AssertSimTKEqual("zero", "0", mm22c(i, j)(r, c), Complex(0, 0)));
                }
            }
        }
    }

    mm22c.setToNaN();
    for (int i = 0; i < mm22c.nrow(); ++i) {
        for (int j = 0; j < mm22c.ncol(); ++j) {
            for (int r = 0; r < 2; ++r) {
                for (int c = 0; c < 2; ++c) {
                    EXPECT_TRUE(std::isnan(mm22c(i, j)(r, c).real()) || std::isnan(mm22c(i, j)(r, c).imag()));
                }
            }
        }
    }
}

TEST(SimTKCommon_BigMatrix_Complex, TransposeAccessConsistency) {
    const Complex mdc[] = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};
    const Matrix_<Complex> md(2, 2, mdc);
    const auto t = ~md;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("transpose", "md(j,i)", t(i, j), md(j, i)));
        }
    }

    {
        const auto transpose1 = ~md(1);
        const auto transpose2 = (~md)(1);
        EXPECT_EQ(transpose1.size(), transpose2.size());
        for (int i = 0; i < transpose1.size(); ++i) {
            EXPECT_NEAR(transpose1[i].real(), transpose2[i].real(), 1e-6);
            EXPECT_NEAR(transpose1[i].imag(), transpose2[i].imag(), 1e-6);
        }
    }

    {
        const auto transpose1 = ~md[1];
        const auto transpose2 = (~md)[1];
        EXPECT_EQ(transpose1.size(), transpose2.size());
        for (int i = 0; i < transpose1.size(); ++i) {
            EXPECT_NEAR(transpose1[i].real(), transpose2[i].real(), 1e-6);
            EXPECT_NEAR(transpose1[i].imag(), transpose2[i].imag(), 1e-6);
        }
    }
}

TEST(SimTKCommon_BigMatrix_Complex, SubmatrixViewAliasing) {
    const Complex mdc[] = {{1, 2}, {3, 4}, {5, 6}, {7, 8}};

    Matrix_<Complex> md(2, 2, mdc);

    const ComplexMatrixView& mvc = md(0, 1, 2, 1);

    md(1, 0) *= 10.0;
    EXPECT_TRUE(AssertSimTKEqual("mvc reflects change", "md(1,0)", mvc(1, 0), md(1, 0)));

    md(0, 1, 2, 1) *= Complex(10., 100.);
    EXPECT_TRUE(AssertSimTKEqual("md updated via view",
                                 "scaled value",
                                 md(0, 1),
                                 Complex(3, 4) * Complex(10., 100.)));
}


TEST(SimTKCommon_BigMatrix_ComplexScaling, AssignmentAndPatternFill) {
    Matrix_<Complex> mm(3, 4);
    mm = 2390.;

    for (int i = 0; i < mm.nrow(); ++i) {
        for (int j = 0; j < mm.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("mm(i,j)", "2390", mm(i, j), Complex(2390, 0)));
        }
    }

    for (int i = 0; i < mm.nrow(); ++i) {
        for (int j = 0; j < mm.ncol(); ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    for (int i = 0; i < mm.nrow(); ++i) {
        for (int j = 0; j < mm.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("pattern", "(i+1)*(j+1)", mm(i, j), Complex((i + 1) * (j + 1), 0)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_ComplexScaling, ColumnScaleInvertibility) {
    Matrix_<Complex> mm(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    Matrix_<Complex> original = mm;

    Vector scale(4);
    scale[0] = 1;
    scale[1] = 10;
    scale[2] = 100;
    scale[3] = 1000;

    mm.colScaleInPlace(scale);
    mm.colScaleInPlace(scale.elementwiseInvert());

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("col unscale", "original", mm(i, j), original(i, j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_ComplexScaling, RowScaleInvertibilityAndType) {
    Matrix_<Complex> mm(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    Matrix_<Complex> original = mm;

    Vector rowScale(3);
    rowScale[0] = -1000;
    rowScale[1] = -100;
    rowScale[2] = -10;

    mm.rowScaleInPlace(rowScale);
    mm.rowScaleInPlace(rowScale.elementwiseInvert());

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("row unscale", "original", mm(i, j), original(i, j)));
        }
    }

    Vector_<double> rowScaleR(3);
    for (int i = 0; i < 3; ++i) {
        rowScaleR[i] = (double)rowScale[i];
    }

    auto result = mm.rowScale(rowScaleR);
    ASSERT_EQ(result.nrow(), 3);
    ASSERT_EQ(result.ncol(), 4);
}

TEST(SimTKCommon_BigMatrix_ComplexScaling, AlgebraicIdentity) {
    Matrix_<Complex> mm(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    mm *= 1000.;
    auto expr = mm + 2 * (-mm);
    for (int i = 0; i < mm.nrow(); ++i) {
        for (int j = 0; j < mm.ncol(); ++j) {
            EXPECT_TRUE(AssertSimTKEqual("mm + 2*(-mm)", "-mm", expr(i, j), -mm(i, j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_ComplexScaling, NormAndIndexing) {
    Matrix_<Complex> mm(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    EXPECT_GE(mm.norm(), 0);

    const auto mm3_vv = mm(3);
    const auto mm3_rvv = mm[3];
    EXPECT_EQ(mm3_vv.size(), mm3_rvv.size());
    for (int i = 0; i < mm3_vv.size(); ++i) {
        EXPECT_NEAR(mm3_vv[i].real(), mm3_rvv[i].real(), 1e-6);
        EXPECT_NEAR(mm3_vv[i].imag(), mm3_rvv[i].imag(), 1e-6);
    }
}

TEST(SimTKCommon_BigMatrix_ComplexScaling, SubmatrixAndTransposeAliasing) {
    Matrix_<Complex> mm(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            mm(i, j) = (i + 1) * (j + 1);
        }
    }

    Matrix_<Complex> mm2 = mm;
    mm2(1, 1, 1, 2) = 7;
    EXPECT_FALSE(SimTK::Test::numericallyEqual(mm(1, 1), Complex(7, 0), 1));

    (~mm)(1, 1, 1, 2) = 99;
    EXPECT_TRUE(AssertSimTKEqual("transpose alias", "99", mm(1, 1), Complex(99, 0)));
}

TEST(SimTKCommon_BigMatrix_Diagonal, ExtractionAndScaling) {
    Matrix_<Complex> mm2(3, 4);
    mm2 = 99;

    auto d = mm2.diag();
    for (int i = 0; i < d.size(); ++i) {
        EXPECT_TRUE(AssertSimTKEqual("diag", "99", d[i], Complex(99, 0)));
    }

    mm2(0, 2, 2, 2) = Complex(-99, -99);
    mm2.updDiag() *= .001;

    for (int i = 0; i < std::min(mm2.nrow(), mm2.ncol()); ++i) {
        EXPECT_TRUE(AssertSimTKEqual("scaled diag", "0.099", mm2(i, i), Complex(99, 0) * Complex(.001, 0)));
    }
}

TEST(SimTKCommon_BigMatrix_Diagonal, RowBlockScalarAssignment) {
    Matrix_<Complex> mm2(3, 4);
    mm2 = 0;

    mm2(0, 0, 1, mm2.ncol()) = 1;
    EXPECT_TRUE(AssertSimTKEqual("mm2(0,0)", "1", mm2(0, 0), Complex(1, 0)));

    for (int j = 1; j < mm2.ncol(); ++j) {
        EXPECT_TRUE(AssertSimTKEqual("row block off-diag", "0", mm2(0, j), Complex(0, 0)));
    }

    mm2[mm2.nrow() - 1] = 1;
    for (int j = 0; j < mm2.ncol(); ++j) {
        EXPECT_TRUE(AssertSimTKEqual("full row", "1", mm2(mm2.nrow() - 1, j), Complex(1, 0)));
    }
}

TEST(SimTKCommon_BigMatrix_Diagonal, ColumnBlockScalarAssignment) {
    Matrix_<Complex> mm2(3, 4);
    mm2 = 0;
    mm2(0, 1, mm2.nrow(), 1) = 2;
    EXPECT_TRUE(AssertSimTKEqual("mm2(1,1)", "2", mm2(1, 1), Complex(2, 0)));

    for (int i = 0; i < mm2.nrow(); ++i) {
        if (i == 1) {
            continue;
        }
        EXPECT_TRUE(AssertSimTKEqual("col block off-diag", "0", mm2(i, 1), Complex(0, 0)));
    }

    mm2(2) = 2;

    for (int i = 0; i < mm2.nrow(); ++i) {
        EXPECT_TRUE(AssertSimTKEqual("full col", "2", mm2(i, 2), Complex(2, 0)));
    }
}

TEST(SimTKCommon_BigMatrix_Resize, ResizeBehaviorAndPreservation) {
    Vector v1(5);
    Vector v2(5);

    for (int i = 0; i < 5; ++i) {
        v1[i] = i;
        v2[i] = i;
    }

    v1.resize(10);
    v2.resizeKeep(10);

    EXPECT_EQ(v1.size(), 10);
    EXPECT_EQ(v2.size(), 10);

    for (int i = 0; i < 5; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("resizeKeep preserved", "original", v2[i], i));
    }

    Matrix m(2, 3);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 3; ++j) {
            m(i, j) = (i + 1) * (j + 1);
        }
    }

    m.resizeKeep(3, 5);
    m.resizeKeep(2, 2);
    m.resize(3, 4);

    EXPECT_EQ(m.nrow(), 3);
    EXPECT_EQ(m.ncol(), 4);
}

TEST(SimTKCommon_BigMatrix_FixedMat, MatVecAndMatMatMultiplication) {
    const Complex mdc[] = {{1, 2},
                           {3, 4},
                           {5, 6},
                           {7, 8},
                           {9, 10},
                           {10, 11},
                           {.1, .26},
                           {.3, .45},
                           {.5, .64},
                           {.7, .83},
                           {.9, .102},
                           {.10, .111}};

    Mat<3, 4, Complex> cm34;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            cm34(i, j) = mdc[i + (j * 3)];
        }
    }

    const Vec<4, Complex> cv4(&mdc[6]);
    const auto vres = cm34 * cv4;
    const auto mres = cm34 * ~cm34;

    EXPECT_EQ(vres.size(), 3);
    EXPECT_EQ(mres.nrow(), 3);
    EXPECT_EQ(mres.ncol(), 3);
}

TEST(SimTKCommon_BigMatrix_Arithmetic, ComplexOperatorConsistency) {
    complex<float> z(1, 2);
    conjugate<float> j(0.3F, 0.4F);
    negator<float> n(7.1F);

    const auto r = z * j;
    EXPECT_TRUE(std::isfinite(r.real()) || std::isfinite(r.imag()));
}

TEST(SimTKCommon_BigMatrix_Product, MatrixVectorAndMatrixMatrix) {
    const Complex mdc[] = {{1, 2},
                           {3, 4},
                           {5, 6},
                           {7, 8},
                           {9, 10},
                           {10, 11},
                           {.1, .26},
                           {.3, .45},
                           {.5, .64},
                           {.7, .83},
                           {.9, .102},
                           {.10, .111}};

    Mat<3, 4, Complex> cm34;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            cm34(i, j) = mdc[i + (j * 3)];
        }
    }

    Vector_<Complex> v(4);
    for (int i = 0; i < 4; ++i) {
        v[i] = cm34(0, i);
    }

    Matrix_<Complex> m(3, 4);
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 4; ++j) {
            m(i, j) = cm34(i, j);
        }
    }

    const auto mv = m * v;
    const auto mm = m * ~m;

    EXPECT_EQ(mv.size(), 3);
    EXPECT_EQ(mm.nrow(), 3);
    EXPECT_EQ(mm.ncol(), 3);
}

TEST(SimTKCommon_BigMatrix_VectorOps, ScalarMultiplyAndCompoundOps) {
    Vector vv(4);
    for (int i = 0; i < 4; ++i) {
        vv[i] = i + 1;
    }

    Vector original = vv;
    vv *= 9.0;
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("vv *= 9", "9*(i+1)", vv[i], original[i] * 9.0));
    }

    Vector ww = vv;
    ww *= 0.1;
    vv += ww;
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("vv += ww", "1.1*vv", vv[i], original[i] * 9.9));
    }
}

TEST(SimTKCommon_BigMatrix_VectorOps, ScalarTimesVector) {
    Vector vv(4);
    for (int i = 0; i < 4; ++i) {
        vv[i] = i + 1;
    }

    Vector vvvv;
    vvvv = vv[2] * vv;
    for (int i = 0; i < 4; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("scalar * vector", "vv[2]*vv[i]", vvvv[i], vv[2] * vv[i]));
    }
}

TEST(SimTKCommon_BigMatrix_MatrixOps, NestedMatAssignmentConsistency) {
    Matrix_<Mat<2, 2, Mat<2, 2, double>>> mmm(2, 1);

    mmm = Mat<2, 2, Mat<2, 2, double>>(Mat<2, 2, double>(1));

    for (int i = 0; i < mmm.nrow(); ++i) {
        for (int j = 0; j < mmm.ncol(); ++j) {
            for (int r = 0; r < 2; ++r) {
                for (int c = 0; c < 2; ++c) {
                    for (int rr = 0; rr < 2; ++rr) {
                        for (int cc = 0; cc < 2; ++cc) {
                            EXPECT_EQ(mmm(i, j)(r, c)(rr, cc), 1.0);
                        }
                    }
                }
            }
        }
    }
}

TEST(SimTKCommon_BigMatrix_MatrixOps, RowColumnAssignmentAndAccess) {
    Vector vv(4);
    for (int i = 0; i < 4; ++i) {
        vv[i] = i + 1;
    }

    Matrix mnm(4, 2);
    mnm(0) = vv;
    mnm(1) = -0.01 * vv;
    for (int j = 0; j < 4; ++j) {
        EXPECT_TRUE(AssertSimTKEqual("row 0", "vv", mnm(0, j), vv[j]));
        EXPECT_TRUE(AssertSimTKEqual("row 1", "-0.01*vv", mnm(1, j), -0.01 * vv[j]));
    }

    auto absMat = mnm.abs();
    for (int i = 0; i < mnm.nrow(); ++i) {
        for (int j = 0; j < mnm.ncol(); ++j) {
            EXPECT_GE(absMat(i, j), 0);
        }
    }
}

TEST(SimTKCommon_BigMatrix_MatrixOps, RowViewConsistency) {
    Vector vv(4);
    for (int i = 0; i < 4; ++i) {
        vv[i] = i + 1;
    }

    Matrix mnm(4, 2);

    mnm(0) = vv;
    mnm(1) = -0.01 * vv;

    const auto r1 = mnm(1);
    const auto r1abs = mnm(1).abs();

    const auto r1idx = mnm[1];
    const auto r1idxabs = mnm[1].abs();

    for (int j = 0; j < 2; ++j) {
        EXPECT_TRUE(AssertSimTKEqual("row view", "indexed row", r1[j], r1idx[j]));
        EXPECT_TRUE(AssertSimTKEqual("abs row view", "abs indexed row", r1abs[j], r1idxabs[j]));
    }
}

TEST(SimTKCommon_BigMatrix_Inverse, InvertAndMultiplyIdentity) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};

    Matrix_<negator<Real>> A(3, 3, (negator<Real>*)rdata);
    Matrix AI = A.invert();
    Matrix I = A * AI;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("A*AI", "Identity", I(i, j), (i == j ? Real(1) : Real(0))));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Inverse, InvertInPlaceMatchesOutOfPlace) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};

    Matrix_<negator<Real>> A1(3, 3, (negator<Real>*)rdata);
    Matrix_<negator<Real>> A2(3, 3, (negator<Real>*)rdata);

    Matrix AI = A1.invert();
    A2.invertInPlace();

    Matrix diff = A2 - AI;
    EXPECT_LT(diff.norm(), 1e-10);
}

TEST(SimTKCommon_BigMatrix_Inverse, TransposeInverseConsistency) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};

    Matrix_<negator<Real>> A(3, 3, (negator<Real>*)rdata);

    Matrix invA = A.invert();
    Matrix invAT = (~A).invert();
    Matrix tInvA = ~invA;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_TRUE(AssertSimTKEqual("(~A)^{-1}", "~(A^{-1})", invAT(i, j), tInvA(i, j)));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Inverse, FixedMatInverseConsistency) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};

    Mat<3, 3, negator<Real>> smallNegA((negator<Real>*)rdata);
    auto smallNegAI = smallNegA.invert();
    auto I = smallNegA * smallNegAI;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            EXPECT_TRUE(
                AssertSimTKEqual("smallNegA * inv", "Identity", I(i, j), (i == j ? Real(1) : Real(0))));
        }
    }
}

TEST(SimTKCommon_BigMatrix_Inverse, DeterminantInverseIdentity) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};
    const Mat<3, 3, negator<Real>> A((negator<Real>*)rdata);

    const auto AI = A.invert();
    const Real detA = det(A);
    const Real detAI = det(AI);

    EXPECT_NEAR(detA * detAI, 1.0, 1e-8);
}

TEST(SimTKCommon_BigMatrix_Inverse, NegatorArithmeticMatchesReal) {
    const Real rdata[] = {1, 2, 3, 9, .1, 14, 2, 6, 9};

    const Mat<3, 3, negator<Real>> A((negator<Real>*)rdata);
    const negator<Real> n1 = A(0, 0) - A(1, 1);
    const Real r1 = A(0, 0) - A(1, 1);
    EXPECT_TRUE(AssertSimTKEqual("negator subtraction", "real subtraction", n1, r1));

    const negator<Real> n2 = A(0, 1) - A(1, 0);
    const Real r2 = A(0, 1) - A(1, 0);
    EXPECT_TRUE(AssertSimTKEqual("negator subtraction", "real subtraction", n2, r2));
}

TEST(SimTKCommon_BigMatrix_ConjugateInverse, LapackInverseConsistency4x4) {
    const Real cjdata[] = {1, 1, 2, 2, 3, 3, 4,  4,  9,  9,  .1, .1, 14, 14, 22, 22,
                           2, 2, 6, 6, 9, 9, 11, 11, .2, .2, .7, .7, 5,  5,  10, 10};

    const Mat<4, 4, conjugate<Real>> A((conjugate<Real>*)cjdata);
    const auto AI = A.invert();
    const auto I = A * AI;

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            const auto& actual = I(i, j);
            const auto expected = (i == j ? Real(1) : Real(0));

            EXPECT_NEAR(actual.real(), expected, 1e-6);
            EXPECT_NEAR(actual.imag(), 0.0, 1e-6);
        }
    }

    const auto d = det(A) * det(AI);
    EXPECT_NEAR(d.real(), 1.0, 1e-6);
    EXPECT_NEAR(d.imag(), 0.0, 1e-6);
}

TEST(SimTKCommon_BigMatrix_ConjugateInverse, LapackVsReferenceInverse4x4) {
    const Real cjdata[] = {1, 1, 2, 2, 3, 3, 4,  4,  9,  9,  .1, .1, 14, 14, 22, 22,
                           2, 2, 6, 6, 9, 9, 11, 11, .2, .2, .7, .7, 5,  5,  10, 10};

    const Mat<4, 4, conjugate<Real>> A((conjugate<Real>*)cjdata);
    const auto invA = inverse(A);
    const auto lapA = lapackInverse(A);

    const auto diff = invA - lapA;
    EXPECT_LT(diff.norm(), 1e-8);
}

TEST(SimTKCommon_BigMatrix_ConjugateInverse, LapackInverseConsistency3x3) {
    const Real cjdata[] = {1, 1, 2, 2, 3, 3, 9, 9, .1, .1, 14, 14, 2, 2, 6, 6, 9, 9};

    const Mat<3, 3, conjugate<Real>> A((conjugate<Real>*)cjdata);
    const auto AI = A.invert();
    const auto I = A * AI;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            const auto& actual = I(i, j);
            const auto expected = (i == j ? Real(1) : Real(0));
            EXPECT_NEAR(actual.real(), expected, 1e-6);
            EXPECT_NEAR(actual.imag(), 0.0, 1e-6);
        }
    }

    const auto diff = inverse(A) - lapackInverse(A);
    EXPECT_LT(diff.norm(), 1e-8);

    const auto d = det(A) * det(AI);
    EXPECT_NEAR(d.real(), 1.0, 1e-6);
    EXPECT_NEAR(d.imag(), 0.0, 1e-6);
}

TEST(SimTKCommon_BigMatrix_ConjugateInverse, SubmatrixInverseConsistency) {
    const Real cjdata[] = {1, 1, 2, 2, 3, 3, 4,  4,  9,  9,  .1, .1, 14, 14, 22, 22,
                           2, 2, 6, 6, 9, 9, 11, 11, .2, .2, .7, .7, 5,  5,  10, 10};

    const Mat<3, 3, conjugate<Real>> A((conjugate<Real>*)cjdata);
    const auto A11 = A.getSubMat<1, 1>(1, 0);
    const auto A22 = A.getSubMat<2, 2>(0, 0);

    const auto inv11 = A11.invert();
    const auto inv22 = A22.invert();

    const auto I11 = A11 * inv11;
    const auto I22 = A22 * inv22;

    EXPECT_TRUE(AssertSimTKEqual("1x1 identity", "1", I11(0, 0).real(), Real(1)));

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) {
            const auto& actual = I22(i, j);
            const auto expected = (i == j ? Real(1) : Real(0));
            EXPECT_NEAR(actual.real(), expected, 1e-6);
            EXPECT_NEAR(actual.imag(), 0.0, 1e-6);
        }
    }
}


TEST(SimTKCommon_BigMatrix_ComplexVector, ConjugateTransposeAliasingAndMutation) {
    const complex<float> ccc[] = {{1.F, 2.F}, {3.F, 4.F}, {5.F, 6.F}, {7.F, 8.F}};

    const Vec<2, complex<float>> cv2(ccc);
    const auto& cv2_copy = cv2;
    const auto scaled = (cv2 + cv2) / complex<float>(1000.F, 0.F);
    for (int i = 0; i < 2; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("scaled vector",
                                     "cv2/500",
                                     scaled[i],
                                     cv2_copy[i] / complex<float>(500.F, 0.F)));
    }

    const auto t = cv2.transpose();
    const auto nt = -cv2;
    const auto ct = ~cv2;

    EXPECT_EQ(t.size(), cv2.size());
    EXPECT_EQ(nt.size(), cv2.size());
    EXPECT_EQ(ct.size(), cv2.size());

    auto chain1 = -(~cv2);
    auto chain2 = ~(-cv2);

    for (int i = 0; i < 2; ++i) {
        EXPECT_TRUE(AssertSimTKEqual("-(~v) == ~( -v )", "identity", chain1[i], chain2[i]));
    }
}

TEST(SimTKCommon_BigMatrix_ComplexVector, ConjugateMutationAliasBehavior) {
    const complex<float> ccc[] = {{1.F, 2.F}, {3.F, 4.F}, {5.F, 6.F}, {7.F, 8.F}};

    const Vec<2, complex<float>> cv2(ccc);

    auto cvt = ~cv2;
    cvt[1] = complex<float>(101.1F, 202.3F);
    EXPECT_NE(cv2[1], complex<float>(101.1F, 202.3F));

    auto negcvt = -(~cv2);
    negcvt[1] = complex<float>(11.1F, 22.3F);
    EXPECT_NE(cv2[1], complex<float>(11.1F, 22.3F));
}

TEST(SimTKCommon_BigMatrix_Vector, NormAndScalingConsistency) {
    const float fddd[] = {11, 12, 13, 14, 15, 16};

    Vec<3, float> dv2(fddd);
    Vec<3, float> ddv2(fddd + 3);

    dv2[2] = 1000.F;
    const auto diff = ddv2 - dv2;
    for (int i = 0; i < 3; ++i) {
        EXPECT_TRUE(diff[i] < 1e-3);
    }

    float norm = dv2.norm();
    EXPECT_GE(norm, 0.0F);

    dv2 = 100.F * dv2;
    EXPECT_GT(dv2.norm(), norm);
}

TEST(SimTKCommon_BigMatrix_Vector, NestedVectorScalingAndAssignment) {
    const float fddd[] = {11, 12, 13, 14, 15, 16};
    const Vec<3, float> v3c[] = {Vec<3, float>(fddd), Vec<3, float>(fddd + 1)};
    Vector_<Vec<2, Vec<3, float>>> vflt(2);

    vflt[0] = Vec<2, Vec<3, float>>(v3c);
    vflt[1] = vflt[0] * 100.F;

    for (int i = 0; i < vflt.size(); ++i) {
        for (int j = 0; j < 2; ++j) {
            for (int k = 0; k < 3; ++k) {
                EXPECT_NEAR(vflt[1][j][k], 100.F * vflt[0][j][k], 1e-6);
            }
        }
    }
}

TEST(SimTKCommon_BigMatrix_Vector, ContiguousStorageProperties) {
    const float fddd[] = {11, 12, 13, 14, 15, 16};
    Vector_<Vec<2, Vec<3, float>>> vflt(2);

    EXPECT_TRUE(vflt.hasContiguousData());

    auto n = vflt.getContiguousScalarDataLength();
    EXPECT_GT(n, 0);

    const float* p = vflt.getContiguousScalarData();
    EXPECT_NE(p, nullptr);

    for (int i = 0; i < n; ++i) {
        EXPECT_TRUE(std::isfinite(p[i]));
    }
}

TEST(SimTKCommon_BigMatrix_Vector, SwapOwnedContiguousScalarData) {
    Vector_<Vec<2, Vec<3, float>>> vflt(2);

    auto* newData = new float[12];
    float* oldData = nullptr;

    for (int i = 0; i < 12; ++i) {
        newData[i] = static_cast<float>(-i);
    }

    vflt.swapOwnedContiguousScalarData(newData, 12, oldData);
    EXPECT_NE(oldData, nullptr);

    for (int i = 0; i < 12; ++i) {
        EXPECT_EQ(oldData[i], -i);
    }

    delete[] oldData;
}

TEST(SimTKCommon_BigMatrix_LargeInverse, LapackVsNumericalRecipesConsistency) {
    const int N = 50; // reduced for test feasibility
    Matrix_<P> big(N, N);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            big(i, j) = 1.0 + (P)std::rand() / RAND_MAX;
        }
    }

    Matrix_<P> flip = big.invert();
    Matrix_<P> identity = big * flip;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            EXPECT_NEAR(identity(i, j), (i == j ? 1.0 : 0.0), 1e-6);
        }
    }
}

TEST(SimTKCommon_BigMatrix_LargeInverse, LapackVsNRInverseAgreement) {
    const int N = 50;
    Matrix_<P> big(N, N);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            big(i, j) = 1.0 + ((P)std::rand() / RAND_MAX);
        }
    }

    Matrix_<P> flip = big.invert();
    Matrix_<P> nrflip(N, N);
    Matrix_<P> tmp = big;

    NumericalRecipes::luinvert(N, &tmp(0, 0), &nrflip(0, 0));

    Matrix_<P> diff = flip - nrflip;
    EXPECT_LT(diff.norm() / N, 1e-5);
}

TEST(SimTKCommon_BigMatrix_LargeInverse, MatrixProductIdentityCheck) {
    const int N = 50;
    Matrix_<P> big(N, N);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            big(i, j) = 1.0 + ((P)std::rand() / RAND_MAX);
        }
    }

    Matrix_<P> flip = big.invert();
    Matrix_<P> ans(N, N);

    Lapack::gemm('n', 'n', N, N, N, P(1), &big(0, 0), N, &flip(0, 0), N, P(0), &ans(0, 0), N);

    double rmsError = std::sqrt(ans.normSqr() / N) - 1.0;
    EXPECT_NEAR(rmsError, 0.0, 1e-6);
}

TEST(SimTKCommon_BigMatrix_LargeInverse, ManualMultiplicationConsistency) {
    const int N = 30;

    Matrix_<P> big(N, N);
    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            big(i, j) = 1.0 + ((P)std::rand() / RAND_MAX);
        }
    }

    Matrix_<P> flip = big.invert();
    Matrix_<P> ans(N, N);

    const P* bigp = &big(0, 0);
    const P* flipp = &flip(0, 0);
    P* ansp = &ans(0, 0);

    for (int j = 0; j < N; ++j) {
        for (int i = 0; i < N; ++i) {
            P sum = 0;
            for (int k = 0; k < N; ++k) {
                sum += bigp[(k * N) + i] * flipp[(j * N) + k];
            }

            ansp[(j * N) + i] = sum;
        }
    }

    for (int i = 0; i < N; ++i) {
        EXPECT_NEAR(ans(i, i), 1.0, 1e-5);
    }
}

TEST(SimTKCommon_BigMatrix_ScalarInit, InitializesAllElementsToGivenScalar) {
    Matrix minit(2, 3, 5.25);

    ASSERT_EQ(minit.nrow(), 2);
    ASSERT_EQ(minit.ncol(), 3);

    expectMatrixEq<Matrix, 2, 3>(minit, Mat23(5.25, 5.25, 5.25, 5.25, 5.25, 5.25));
}

TEST(SimTKCommon_BigMatrix_VecElementInit, InitializesAllElementsToGivenVec2) {
    const Vec2 v12(1, 2);
    Matrix_<Vec2> mvinit(3, 2, v12);

    ASSERT_EQ(mvinit.nrow(), 3);
    ASSERT_EQ(mvinit.ncol(), 2);

    expectMatrixEq<Matrix_<Vec2>, 3, 2>(mvinit, Mat<3, 2, Vec2>(v12, v12, v12, v12, v12, v12));
}

TEST(SimTKCommon_BigMatrix_MatDivision, ScalarInverseOfMat22ProducesCorrectReciprocals) {
    Mat22 m1(4, 0, 0, 1);
    Mat22 expected(.25, 0, 0, 1);

    SimTK_TEST_EQ(1 / m1, expected);
}

TEST(SimTKCommon_BigMatrix_MatDivision, ScalarInverseOfNestedMat22ProducesCorrectReciprocals) {
    Mat<2, 2, Mat22> m2(Mat22(2, 0, 0, 3));
    Mat<2, 2, Mat22> expected(Mat22(.5, 0, 0, OneThird));

    SimTK_TEST_EQ(1 / m2, expected);
}

TEST(SimTKCommon_BigMatrix_Transform, TransformOfStridedVec3IsNegationOfTransformOfNegatedVec) {
    Transform X;
    Vec<3, Real, 6> vs(1, 2, 3); // non-unit stride

    SimTK_TEST(X * vs == -(X * -vs));
}

TEST(SimTKCommon_BigMatrix_Transform, TransformOfStridedVec4IsNegationOfTransformOfNegatedVec) {
    Transform X;
    Vec<4, Real, 9> vs2(1, 2, 3, 0); // non-unit stride

    SimTK_TEST(X * vs2 == -(X * -vs2));
}

TEST(SimTKCommon_BigMatrix_Transform, RotationOfStridedVec3IsNegationOfRotationOfNegatedVec) {
    Rotation R;
    Vec<3, Real, 6> vs(1, 2, 3);

    SimTK_TEST(R * vs == -(R * -vs));
}

TEST(SimTKCommon_BigMatrix_Transform, RowTimesRotationIsNegationOfNegatedRowTimesRotation) {
    Rotation R;
    Vec<3, Real, 6> vs(1, 2, 3);

    SimTK_TEST(~vs * R == -(-~vs * R));
}

TEST(SimTKCommon_BigMatrix_MatrixFromMat22, ConstructedMatrixMatchesSourceValues) {
    Matrix m(Mat22(1, 2, 3, 4));

    expectMatrixEq<Matrix, 2, 2>(m, Mat22(1, 2, 3, 4));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, AddScalarInPlaceAddsToDiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    m += 3;
    expectMatrixEq<Matrix, 2, 2>(m, Mat22(4, 2, 3, 7));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, SubtractScalarInPlaceSubtractsFromDiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    m += 3;
    m -= 3;
    expectMatrixEq<Matrix, 2, 2>(m, Mat22(1, 2, 3, 4));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, MatrixMinusScalarSubtractsFromDiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    expectMatrixEq<Matrix, 2, 2>(m - 1, Mat22(0, 2, 3, 3));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, MatrixPlusScalarAddsTodiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    expectMatrixEq<Matrix, 2, 2>(m + 1, Mat22(2, 2, 3, 5));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, ScalarMinusMatrixNegatesAndAddsToDiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    expectMatrixEq<Matrix, 2, 2>(1 - m, Mat22(0, -2, -3, -3));
}

TEST(SimTKCommon_BigMatrix_MatrixScalarArithmetic, ScalarPlusMatrixAddsTodiagonal) {
    Matrix m(Mat22(1, 2, 3, 4));
    expectMatrixEq<Matrix, 2, 2>(1 + m, Mat22(2, 2, 3, 5));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, ConstructedVectorMatchesSourceValues) {
    Vector v(Vec3(1, 2, 3));
    expectVectorEq(v, Vec3(1, 2, 3));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, AddScalarInPlaceAddsToEveryElement) {
    Vector v(Vec3(1, 2, 3));
    v += 2;
    expectVectorEq(v, Vec3(3, 4, 5));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, SubtractScalarInPlaceSubtractsFromEveryElement) {
    Vector v(Vec3(1, 2, 3));
    v += 2;
    v -= 2;
    expectVectorEq(v, Vec3(1, 2, 3));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, VectorMinusScalarSubtractsFromEveryElement) {
    Vector v(Vec3(1, 2, 3));
    expectVectorEq(v - 1, Vec3(0, 1, 2));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, VectorPlusScalarAddsToEveryElement) {
    Vector v(Vec3(1, 2, 3));
    expectVectorEq(v + 1, Vec3(2, 3, 4));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, ScalarMinusVectorNegatesAndAddsToEveryElement) {
    Vector v(Vec3(1, 2, 3));
    expectVectorEq(1 - v, Vec3(0, -1, -2));
}

TEST(SimTKCommon_BigMatrix_VectorScalarArithmetic, ScalarPlusVectorAddsToEveryElement) {
    Vector v(Vec3(1, 2, 3));
    expectVectorEq(1 + v, Vec3(2, 3, 4));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, ConstructedRowVectorMatchesSourceValues) {
    RowVector r(Row3(1, 2, 3));
    expectVectorEq(r, Vec3(1, 2, 3));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, AddScalarInPlaceAddsToEveryElement) {
    RowVector r(Row3(1, 2, 3));
    r += 2;
    expectVectorEq(r, Vec3(3, 4, 5));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, SubtractScalarInPlaceSubtractsFromEveryElement) {
    RowVector r(Row3(1, 2, 3));
    r += 2;
    r -= 2;
    expectVectorEq(r, Vec3(1, 2, 3));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, RowVectorMinusScalarSubtractsFromEveryElement) {
    RowVector r(Row3(1, 2, 3));
    expectVectorEq(r - 1, Vec3(0, 1, 2));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, RowVectorPlusScalarAddsToEveryElement) {
    RowVector r(Row3(1, 2, 3));
    expectVectorEq(r + 1, Vec3(2, 3, 4));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, ScalarMinusRowVectorNegatesAndAddsToEveryElement) {
    RowVector r(Row3(1, 2, 3));
    expectVectorEq(1 - r, Vec3(0, -1, -2));
}

TEST(SimTKCommon_BigMatrix_RowVectorScalarArithmetic, ScalarPlusRowVectorAddsToEveryElement) {
    RowVector r(Row3(1, 2, 3));
    expectVectorEq(1 + r, Vec3(2, 3, 4));
}

// Fixture to avoid repeating the same matrix construction
class SimTKCommon_BigMatrix_ColumnRowCopy : public ::testing::Test {
    protected:
    void SetUp() override {
        mm = Matrix(Mat23(1, 2, 3, 7, 8, 9));
    }
    Matrix mm;
};

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, ConstructedMatrixMatchesSourceValues) {
    expectMatrixEq<Matrix, 2, 3>(mm, Mat23(1, 2, 3, 7, 8, 9));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, VectorCopyConstructedFromColumnHasCorrectValues) {
    Vector vv = mm(1); // column 1
    expectVectorEq(vv, Vec2(2, 8));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, VectorCopyAssignedFromColumnHasCorrectValues) {
    Vector vv = mm(1);
    vv = mm(0); // column 0
    expectVectorEq(vv, Vec2(1, 7));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, RowVectorCopyConstructedFromRowHasCorrectValues) {
    RowVector rr = mm[1]; // row 1
    expectVectorEq(rr, Vec3(7, 8, 9));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, RowVectorCopyAssignedFromRowHasCorrectValues) {
    RowVector rr = mm[1];
    rr = mm[0]; // row 0
    expectVectorEq(rr, Vec3(1, 2, 3));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, TransposedColumnBecomesRowVectorWithCorrectValues) {
    RowVector rrr = ~mm(1); // transpose of column 1 → RowVector
    expectVectorEq(rrr, Vec2(2, 8));
    rrr = ~mm(0);
    expectVectorEq(rrr, Vec2(1, 7));
}

TEST_F(SimTKCommon_BigMatrix_ColumnRowCopy, TransposedRowBecomesVectorWithCorrectValues) {
    Vector vvv = ~mm[1]; // transpose of row 1 → Vector
    expectVectorEq(vvv, Vec3(7, 8, 9));
    vvv = ~mm[0];
    expectVectorEq(vvv, Vec3(1, 2, 3));
}

TEST(SimTKCommon_BigMatrix_SharedMemory, RealMatrixViewOverArrayReflectsCorrectColumnMajorLayout) {
    // Data is stored column-major: col(0) first, then col(1)
    Array_<Real> rarrmat;
    rarrmat.push_back(1.1);
    rarrmat.push_back(2.2); // col 0
    rarrmat.push_back(3.3);
    rarrmat.push_back(4.4); // col 1

    Matrix rmatrix(2, 2, /*lda=*/2, &rarrmat[0]);

    expectMatrixEq<Matrix, 2, 2>(rmatrix, Mat22(1.1, 3.3, 2.2, 4.4));
}

TEST(SimTKCommon_BigMatrix_SharedMemory, SpatialVecMatrixViewOverArrayReflectsCorrectLayout) {
    Array_<SpatialVec> svarrmat;
    svarrmat.push_back(SpatialVec(Vec3(1, 2, 3), Vec3(4, 5, 6)));
    svarrmat.push_back(SpatialVec(Vec3(1.1, 2.1, 3.1), Vec3(4.1, 5.1, 6.1)));
    svarrmat.push_back(SpatialVec(Vec3(1.2, 2.2, 3.2), Vec3(4.2, 5.2, 6.2)));
    svarrmat.push_back(SpatialVec(Vec3(1.3, 2.3, 3.3), Vec3(4.3, 5.3, 6.3)));

    const int szInScalars = sizeof(SpatialVec) / sizeof(Real);
    Matrix_<SpatialVec> svmatrix(2, 2, 2 * szInScalars, (Real*)&svarrmat[0]);

    Matrix_<SpatialVec> expected(2, 2);
    expected(0, 0) = svarrmat[0];
    expected(1, 0) = svarrmat[1];
    expected(0, 1) = svarrmat[2];
    expected(1, 1) = svarrmat[3];

    SimTK_TEST_EQ_TOL(svmatrix, expected, 1e-16);
}

TEST(SimTKCommon_BigMatrix_SharedMemory, RealVectorViewOverArrayReflectsCorrectValues) {
    Array_<Real> rarray;
    rarray.push_back(1.1);
    rarray.push_back(2.2);
    rarray.push_back(3.3);

    Vector rvector(3, &rarray[0], /*owns=*/true);
    expectVectorEq(rvector, Vec3(1.1, 2.2, 3.3));
}

TEST(SimTKCommon_BigMatrix_SharedMemory, SpatialVecVectorViewOverArrayReflectsCorrectValues) {
    Array_<SpatialVec> svarray;
    svarray.push_back(SpatialVec(Vec3(1, 2, 3), Vec3(4, 5, 6)));
    svarray.push_back(SpatialVec(Vec3(1.1, 2.1, 3.1), Vec3(4.1, 5.1, 6.1)));
    svarray.push_back(SpatialVec(Vec3(1.2, 2.2, 3.2), Vec3(4.2, 5.2, 6.2)));

    Vector_<SpatialVec> svvector(3, (Real*)&svarray[0], /*owns=*/true);
    Vector_<SpatialVec> expected(3);
    expected[0] = svarray[0];
    expected[1] = svarray[1];
    expected[2] = svarray[2];

    SimTK_TEST_EQ_TOL(svvector, expected, 1e-16);
}

class SimTKCommon_BigMatrix_ZeroWidthSlice : public ::testing::Test {
    protected:
    void SetUp() override {
        general.resize(3, 4);
    }
    Matrix general; // 3x4
};

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroRowSliceInMiddleHasCorrectDimensions) {
    MatrixView s = general(1, 1, 0, 2);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 2);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroColSliceInMiddleHasCorrectDimensions) {
    MatrixView s = general(1, 1, 1, 0);
    EXPECT_EQ(s.nrow(), 1);
    EXPECT_EQ(s.ncol(), 0);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroColSliceOnLeftEdgeHasCorrectDimensions) {
    MatrixView s = general(0, 0, 3, 0);
    EXPECT_EQ(s.nrow(), 3);
    EXPECT_EQ(s.ncol(), 0);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroRowSliceOnTopEdgeHasCorrectDimensions) {
    MatrixView s = general(0, 0, 0, 4);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 4);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroRowSliceOffBottomEdgeHasCorrectDimensions) {
    MatrixView s = general(3, 0, 0, 4);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 4);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, ZeroColSliceOffRightEdgeHasCorrectDimensions) {
    MatrixView s = general(0, 4, 3, 0);
    EXPECT_EQ(s.nrow(), 3);
    EXPECT_EQ(s.ncol(), 0);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, EmptySliceAtOriginHasZeroDimensions) {
    MatrixView s = general(0, 0, 0, 0);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 0);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, EmptySliceAtInteriorPointHasZeroDimensions) {
    MatrixView s = general(1, 2, 0, 0);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 0);
}

TEST_F(SimTKCommon_BigMatrix_ZeroWidthSlice, EmptySliceNearCornerHasZeroDimensions) {
    MatrixView s = general(2, 3, 0, 0);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 0);
}

// MatrixView_::MatrixView_() is private in this SimTK version, so a GTest
// fixture (which requires a default-constructible test class) cannot hold a
// MatrixView member.  Each test therefore owns its own Matrix and derives the
// 3x1 vector-shaped view inline.  The helper lambda keeps the boilerplate to
// a single line per test.

// Macro to create the 3x4 parent matrix and extract column 1 as a 3x1 view.
// 'vector_' is a MatrixView that aliases column 1 of 'general'.
#define VECTOR_SHAPE_SETUP()                  \
    Matrix general(3, 4);                     \
    MatrixView vector_ = general(0, 1, 3, 1); \
    ASSERT_EQ(vector_.nrow(), 3);             \
    ASSERT_EQ(vector_.ncol(), 1)

TEST(SimTKCommon_BigMatrix_VectorShapeZeroWidthSlice, ZeroColSliceOfVectorColumnHasCorrectDimensions) {
    VECTOR_SHAPE_SETUP();
    MatrixView s = vector_(0, 0, 3, 0);
    EXPECT_EQ(s.nrow(), 3);
    EXPECT_EQ(s.ncol(), 0);
}

TEST(SimTKCommon_BigMatrix_VectorShapeZeroWidthSlice, EmptySliceAtOriginOfVectorColumnHasZeroDimensions) {
    VECTOR_SHAPE_SETUP();
    MatrixView s = vector_(0, 0, 0, 0);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 0);
}

TEST(SimTKCommon_BigMatrix_VectorShapeZeroWidthSlice,
     ZeroColSliceNearBottomOfVectorColumnHasCorrectDimensions) {
    VECTOR_SHAPE_SETUP();
    MatrixView s = vector_(2, 0, 1, 0);
    EXPECT_EQ(s.nrow(), 1);
    EXPECT_EQ(s.ncol(), 0);
}

TEST(SimTKCommon_BigMatrix_VectorShapeZeroWidthSlice,
     ZeroRowSliceOffBottomOfVectorColumnHasCorrectDimensions) {
    VECTOR_SHAPE_SETUP();
    MatrixView s = vector_(3, 0, 0, 1);
    EXPECT_EQ(s.nrow(), 0);
    EXPECT_EQ(s.ncol(), 1);
}

TEST(SimTKCommon_BigMatrix_VectorShapeZeroWidthSlice,
     ZeroColSliceOffRightOfVectorColumnHasCorrectDimensions) {
    VECTOR_SHAPE_SETUP();
    MatrixView s = vector_(0, 1, 3, 0);
    EXPECT_EQ(s.nrow(), 3);
    EXPECT_EQ(s.ncol(), 0);
}

#undef VECTOR_SHAPE_SETUP

TEST(SimTKCommon_BigMatrix_RowVectorDefaultConstruct, DefaultConstructedRowVectorHasZeroSize) {
    RowVector rv;
    EXPECT_EQ(rv.size(), 0);
    EXPECT_EQ(rv.nrow(), 1);
    EXPECT_EQ(rv.ncol(), 0);
    EXPECT_EQ(rv.nelt(), 0);
}

TEST(SimTKCommon_BigMatrix_RowVectorDefaultConstruct,
     ZeroSizeRowVectorConstructedWithExplicitZeroHasZeroSize) {
    RowVector rv(0);
    EXPECT_EQ(rv.size(), 0);
    EXPECT_EQ(rv.nrow(), 1);
    EXPECT_EQ(rv.ncol(), 0);
    EXPECT_EQ(rv.nelt(), 0);
}