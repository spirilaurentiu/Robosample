/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Matthew Millard (the testNaturalCubicSpline code)            *
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
#include <vector>

#include "SimTKmath.h"
#include "Util.hpp"

using namespace SimTK;
using namespace std;

const Real TESTTOL = 1e-4;

/**
 * This function will compute the value and first two
 * derivatives of an analytic function at the point x.
 *
 * @params x       the input value
 * @params fcnType the function to compute. There are
 *                 currently 5 choices (see below)
 * @returns Vector: a 3x1 vector of the value,
 *                 first derivative and second derivative
 */
auto getAnalyticFunction(double x, int fcnType) -> Vector {
    Vector fdF(3);
    fdF = -1;

    switch (fcnType) {
        case 0: // f(x) = 0;
            fdF = 0;
            break;
        case 1:             // f(x) = 2*x
            fdF(0) = 2 * x; // f
            fdF(1) = 2;     // fx
            fdF(2) = 0;
            break;
        case 2:             // f(x) = x^2
            fdF(0) = x * x; // f
            fdF(1) = 2 * x; // fx
            fdF(2) = 2;
            break;
        case 3: // f(x) = 2*x + x*x;
            fdF(0) = (2 * x) + (x * x);
            fdF(1) = 2 + (2 * x);
            fdF(2) = 2;
            break;
        case 4: // f(x)  =2*x + x*x + 5*x*x*x
            fdF(0) = (2 * x) + (x * x) + (5 * x * x * x);
            fdF(1) = 2 + (2 * x) + (15 * x * x);
            fdF(2) = 2 + (30 * x);
            break;
        case 5: // fx(x) = sin(x)
            fdF(0) = sin(x);
            fdF(1) = cos(x);
            fdF(2) = -sin(x);
            break;
        default:
            SCOPED_TRACE("Invalid fcnType in getAnalyticFunction: " + to_string(fcnType));
            ADD_FAILURE() << "Invalid fcnType: " << fcnType << "\n";
    }

    return fdF;
}

/**
 * This function tests the accuracy of the natural cubic spline sp.
 * The accuracy of the spline is tested in the following manner:
 *
 *    a.    Spline must pass through the knots given
 *             -Error between spline and input data at the knots
 *             (should be zero)
 *    b.   The first derivatives are continuous at the knot points
 *             -Error between the value of the first derivative at
 *             the knot point, and what a linear extrapolation would
 *             predict just to the left and right ofthe knot point.
 *             (should be zero, within a tolerace affected by the
 *             step size in xD)
 *    c.   The second derivatives are continuous at the knots points
 *             -Error between the value of the numerically calculated
 *             derivative at the knot point, and what a linear
 *             extrapolation would predict just to the left and
 *             right of the knot point. (should be zero, within a
 *             tolerace affected by the step size in xD)
 *    d.  The second derivative is zero at the end points.
 *             -Numerically calculated extrapolation of the 2nd
 *             derivative should be zero at the end points within
 *             some tolerance
 *
 */
auto benchmarkNaturalCubicSpline(Function* sp, Vector xK, Vector yK, Vector xM, Vector xD) -> Vector {
    int size = xK.size();
    int sizeD = xD.size();
    int sizeDK = xD.size() / (xK.size() - 1);
    double deltaD = (xK(xK.size() - 1) - xK(0)) / xD.size();

    Matrix ysp_K(size, 2);
    Matrix ysp_M(size - 1, 2);
    Matrix ysp_D(sizeD, 4);
    Vector errVec(4);
    errVec = 1;
    ysp_K = 0;
    ysp_M = 0;
    ysp_D = 0;

    vector<int> derOrder(1);
    derOrder[0] = 0;

    // Evaluate the spline at the knots, the mid points and then a dense sample
    Vector tmpV1(1);
    double xVal = 0;
    for (int i = 0; i < size; i++) {
        xVal = xK(i);
        tmpV1(0) = xK(i);
        ysp_K(i, 0) = sp->calcValue(tmpV1);
        ysp_K(i, 1) = sp->calcDerivative(derOrder, tmpV1);
    }
    for (int i = 0; i < size - 1; i++) {
        xVal = xM(i);
        tmpV1(0) = xM(i);
        ysp_M(i, 0) = sp->calcValue(tmpV1);
        ysp_M(i, 1) = sp->calcDerivative(derOrder, tmpV1);
    }
    for (int i = 0; i < sizeD; i++) {
        xVal = xD(i);
        tmpV1(0) = xD(i);
        ysp_D(i, 0) = sp->calcValue(tmpV1);
        ysp_D(i, 1) = sp->calcDerivative(derOrder, tmpV1);
    }

    // Compute the second derivative of the spline (using central differences), and linearly interpolate to
    // get the end points. The end points should go to exactly zero because the second derivative is linear in
    // a cubic spline, as is the linear extrapolation Also compute the 3rd derivative using the same method.
    // The 3rd derivative is required in the test to determine if the second derivative is continuous at the
    // knots or not.
    ysp_D(2) = getCentralDifference(xD, ysp_D(1), true);
    ysp_D(3) = getCentralDifference(xD, ysp_D(2), true);

    // Now check to see if the splines meet the conditions of a natural cubic spline:
    Vector tmpK(size, size);
    Vector tmpM(size - 1, size - 1);

    // a. Spline passes through all knot points given
    tmpK = yK - ysp_K(0);
    errVec(0) = tmpK.norm();

    // b. The first derivative is continuous at the knot points.
    // Apply a continuity test to the data points that defines the second derivative
    // Continuity test: a linear extrapolation of first derivative of the curve in interest on either side of
    // the point in interest should equal the point in interest
    double ykL;
    double ykR;
    double y0L;
    double dydxL;
    double y0R;
    double dydxR = 0;

    for (int i = 1; i < size - 1; i++) {
        y0L = ysp_D((i * sizeDK) - 1, 1);
        y0R = ysp_D((i * sizeDK) + 1, 1);
        dydxL = ysp_D((i * sizeDK) - 1, 2); // Found using central differences
        dydxR = ysp_D((i * sizeDK) + 1, 2); // Found using central differences
        ykL = y0L + (dydxL * deltaD);
        ykR = y0R - (dydxR * deltaD);
        errVec(1) = (ysp_D(i * sizeDK, 1) - ykL) + (ysp_D(i * sizeDK, 1) - ykR);
    }

    // c. The second derivative is continuous at the knot points.
    // Apply a continuity test to the data points that define the second derivative. This also tests if the
    // first derivative is smooth.
    // Continuity test: a linear extrapolation of first derivative of the curve in interest on either side of
    // the point in interest should equal the point in interest;
    for (int i = 1; i < size - 1; i++) {
        y0L = ysp_D((i * sizeDK) - 1, 2);
        y0R = ysp_D((i * sizeDK) + 1, 2);
        dydxL = ysp_D((i * sizeDK) - 1, 3); // Found using central differences
        dydxR = ysp_D((i * sizeDK) + 1, 3); // Found using central differences
        ykL = y0L + (dydxL * deltaD);
        ykR = y0R - (dydxR * deltaD);
        errVec(2) = (ysp_D(i * sizeDK, 2) - ykL) + (ysp_D(i * sizeDK, 2) - ykR);
    }

    // d. The second derivative is zero at the end points
    errVec(3) = abs(ysp_D(0, 2)) + abs(ysp_D(sizeD - 1, 2));

    return errVec;
}

void testNaturalCubicSpline(int fcnType) {
    std::string fcnName;
    switch (fcnType) {
        case 0:
            fcnName = "f(x) = 0";
            break;
        case 1:
            fcnName = "f(x) = 2*x";
            break;
        case 2:
            fcnName = "f(x) = x^2";
            break;
        case 3:
            fcnName = "f(x) = 2*x + x^2";
            break;
        case 4:
            fcnName = "f(x) = 2*x + x^2 + 5x^3";
            break;
        case 5:
            fcnName = "f(x) = sin(x)";
            break;
        default:
            SCOPED_TRACE("Invalid fcnType in testNaturalCubicSpline: " + to_string(fcnType));
            ADD_FAILURE() << "Invalid fcnType: " << fcnType << "\n";
            break;
    }
    SCOPED_TRACE("Testing natural cubic spline with function: " + fcnName);

    // Number of knot points
    const int size = 6;

    // Number of points per knot in the densely sampled vector
    const int sizeDK = 100;

    // Number of points in a densely sampled interpolation
    int sizeD = sizeDK * (size - 1);

    // Domain vector variables
    double xmin = Pi / 4; // Value of first knot
    double xmax = Pi / 2; // Value of the final knot
    double deltaX = (xmax - xmin) / (size - 1);
    double deltaD = (xmax - xmin) / (sizeD - 1);
    double etime = 0;

    // This matrix stores the results of the 5 tests in each row entry, for each of the 2 spline classes
    // tested. SimTK SplineFitter results are stored in column 0
    // OpenSim::NaturalCubicSpline results are stored in column 1 testResults.elementwiseAssign(0.0);
    Matrix testResults(4, 1);

    testResults = -1;
    Vector tmpV1(1);

    // Generate initialization knot points (denoted by a 'K') and the mid points (denoted by a 'M') and for
    // the densely sampled interpolation vector (denoted by a 'D')
    Vector xK(size);
    Vector xM(size - 1);
    Vector xD(sizeD);
    Matrix yK(size, 3);
    Matrix yM(size - 1, 3);
    Matrix yD(sizeD, 3);

    for (int i = 0; i < size; i++) {
        xK(i) = xmin + (((double)i) * deltaX);
        if (i < size - 1) {
            xM(i) = xmin + (deltaX / (double)2) + (((double)i) * deltaX);
        }
    }
    for (int i = 0; i < sizeD; i++) {
        xD(i) = xmin + (deltaD * (double)i);
    }

    // Get the function values at the knot points
    Vector tmp(3);
    tmp = 0;
    for (int i = 0; i < size; i++) {
        tmp = getAnalyticFunction(xK(i), fcnType);
        for (int k = 0; k < 3; k++) {
            yK(i, k) = tmp(k);
        }
    }
    Vector yKVal = yK(0);

    // Get the function y, dy, ddy at the mid points
    for (int i = 0; i < size - 1; i++) {
        tmp = getAnalyticFunction(xM(i), fcnType);
        for (int k = 0; k < 3; k++) {
            yM(i, k) = tmp(k);
        }
    }
    // Get the function y, dy, ddy at the dense points
    for (int i = 0; i < sizeD; i++) {
        tmp = getAnalyticFunction(xD(i), fcnType);
        for (int k = 0; k < 3; k++) {
            yD(i, k) = tmp(k);
        }
    }

    // SplineFitter
    Vector sfDerivs1(xK.size());
    Spline_<Real> sTK = SplineFitter<Real>::fitForSmoothingParameter(3, xK, yKVal, 0.0).getSpline();

    testResults(0) = benchmarkNaturalCubicSpline(&sTK, xK, yK(0), xM, xD);

    // Run numerical assertions on each test
    double tol = 0;
    for (int k = 0; k < testResults.ncol(); k++) {
        for (int i = 0; i < testResults.nrow(); i++) {
            switch (i) {
                case 0: // Equal at knots
                    tol = 1e-14;
                    break;
                case 1: // Continuous 1st derivative
                    tol = deltaD;
                    break;
                case 2: // Continuous 2nd derivative
                    tol = 10 * deltaD;
                    break;
                case 3: // 2nd derivative zero at end points
                    tol = deltaD / 10;
                    break;
                case 4: // 3rd derivative zero at end points
                case 5: // 4th derivative zero at end points
                    tol = 100000000 * deltaD;
                    break;
                default:
                    SCOPED_TRACE("Invalid fcnType in testNaturalCubicSpline: " + to_string(fcnType));
                    ADD_FAILURE() << "Invalid fcnType: " << fcnType << "\n";
                    break;
            }
            EXPECT_NEAR(testResults(i, k), 0, tol)
                << "Test " << i << " failed for spline class " << k << " with function type " << fcnType
                << " with value " << testResults(i, k) << " and tolerance " << tol;
        }
    }
}

TEST(Spline, Spline) {
    Vector_<Vec3> coeff(5);
    coeff[0] = Vec3(0, 1, 2);
    coeff[1] = Vec3(1, 4, 1);
    coeff[2] = Vec3(2, 2, 20);
    coeff[3] = Vec3(1, -1, 2);
    coeff[4] = Vec3(0, 0, 1);
    Vector x(Vec5(0, 1, 2, 5, 10));

    // Create a linear spline, and verify that it interpolates linearly between the control points.
    Spline_<Vec3> spline(1, x, coeff);
    for (int i = 0; i < x.size(); ++i) {
        SimTK_TEST_EQ(coeff[i], spline.calcValue(Vector(1, x[i])));
    }
    std::vector<int> deriv;
    deriv.push_back(0);
    for (int i = 0; i < x.size() - 1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i + 1.0) / 12.0;
            Real t = x[i] + (fract * (x[i + 1] - x[i]));

            const auto expectedValue = coeff[i] + fract * (coeff[i + 1] - coeff[i]);
            const auto actualValue = spline.calcValue(Vector(1, t));

            const auto expectedDeriv = (coeff[i + 1] - coeff[i]) / (x[i + 1] - x[i]);
            const auto actualDeriv = spline.calcDerivative(deriv, Vector(1, t));

            EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualValue, expectedValue, 1, TESTTOL))
                << " at t=" << t << " expected " << expectedValue << " actual " << actualValue;
            EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualDeriv, expectedDeriv, 1, TESTTOL))
                << " at t=" << t << " expected " << expectedDeriv << " actual " << actualDeriv;
        }
    }

    // Create a cubic spline and verify the derivative calculations.
    spline = Spline_<Vec3>(3, x, coeff);
    const Real delta = 1e-10;

    for (int i = 0; i < x.size() - 1; ++i) {
        for (int j = 0; j < 10; ++j) {
            const Real fract = (i + 1.0) / 12.0;
            const Real t = x[i] + (fract * (x[i + 1] - x[i]));
            const Vec3 value1 = spline.calcValue(Vector(1, t - delta));
            const Vec3 value2 = spline.calcValue(Vector(1, t + delta));

            const auto expectedDeriv = 0.5 * (value1 + value2);
            const auto actualDeriv = spline.calcValue(Vector(1, t));
            EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualDeriv, expectedDeriv, 1, TESTTOL))
                << " at t=" << t << " expected " << expectedDeriv << " actual " << actualDeriv;
        }
    }
}

TEST(Spline, SplineFitter) {
    Real stddev = 0.5;
    const int n = 100;
    Random::Gaussian random(0.0, stddev);
    Vector x(n);
    Vector_<Vec3> truey(n);
    Vector_<Vec3> y(n);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i * 0.1;
        truey[i] = Vec3(sin(x[i]), 3.0 * sin(2 * x[i]), cos(x[i]));
        y[i] = truey[i] + Vec3(random.getValue(), random.getValue(), random.getValue());
    }
    SplineFitter<Vec3> fitter = SplineFitter<Vec3>::fitFromGCV(3, x, y);
    Spline_<Vec3> spline1 = fitter.getSpline();

    // The fitting should have reduced the error.
    Vec3 originalError = mean(abs(y - truey));
    Vec3 fittedError = mean(abs(spline1.getControlPointValues() - truey));
    EXPECT_LE(fittedError[0], originalError[0]);
    EXPECT_LE(fittedError[1], originalError[1]);
    EXPECT_LE(fittedError[2], originalError[2]);

    // If we perform the fitting again, explicitly specifying the same value for the smoothing parameter, it
    // should produce identical results.
    const auto& expectedControlPoints = spline1.getControlPointValues();
    const auto actualControlPoints =
        SplineFitter<Vec3>::fitForSmoothingParameter(3, x, y, fitter.getSmoothingParameter())
            .getSpline()
            .getControlPointValues();
    EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualControlPoints, expectedControlPoints, 1, TESTTOL))
        << " expected " << expectedControlPoints << " actual " << actualControlPoints;

    // Likewise, specifying the same number of degrees of freedom should produce identical results
    const auto& expectedControlPoints2 = spline1.getControlPointValues();
    const auto actualControlPoints2 = SplineFitter<Vec3>::fitFromDOF(3, x, y, fitter.getDegreesOfFreedom())
                                          .getSpline()
                                          .getControlPointValues();
    EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualControlPoints2, expectedControlPoints2, 1, TESTTOL))
        << " expected " << expectedControlPoints2 << " actual " << actualControlPoints2;

    // If we specify a smoothing parameter of 0, it should exactly reproduce the original data.
    const Spline_<Vec3> nosmoothing = SplineFitter<Vec3>::fitForSmoothingParameter(3, x, y, 0.0).getSpline();
    for (int i = 0; i < x.size(); ++i) {
        const auto& expectedValue = y[i];
        const auto actualValue = nosmoothing.calcValue(Vector(1, x[i]));
        EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualValue, expectedValue, 1, TESTTOL))
            << " at x=" << x[i] << " expected " << expectedValue << " actual " << actualValue;
    }
}

TEST(Spline, RealSpline) {
    Vector coeff(5);
    coeff[0] = 0;
    coeff[1] = 1;
    coeff[2] = 2;
    coeff[3] = 1;
    coeff[4] = 0;
    Vector x(Vec5(0, 1, 2, 5, 10));

    // Create a linear spline, and verify that it interpolates linearly between the control points.
    Spline spline(1, x, coeff);
    for (int i = 0; i < x.size(); ++i) {
        SimTK_TEST_EQ_TOL(coeff[i], spline.calcValue(Vector(1, x[i])), TESTTOL);
    }
    Array_<int> deriv;
    deriv.push_back(0);
    for (int i = 0; i < x.size() - 1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i + 1.0) / 12.0;
            Real t = x[i] + (fract * (x[i + 1] - x[i]));

            const auto expectedDeriv = coeff[i] + (fract * (coeff[i + 1] - coeff[i]));
            const auto actualDeriv = spline.calcValue(Vector(1, t));
            EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualDeriv, expectedDeriv, 1, TESTTOL))
                << " at t=" << t << " expected " << expectedDeriv << " actual " << actualDeriv;

            const auto expectedDeriv2 = (coeff[i + 1] - coeff[i]) / (x[i + 1] - x[i]);
            const auto actualDeriv2 = spline.calcDerivative(deriv, Vector(1, t));
            EXPECT_TRUE(SimTK ::Test ::numericallyEqual(actualDeriv2, expectedDeriv2, 1, TESTTOL))
                << " at t=" << t << " expected " << expectedDeriv2 << " actual " << actualDeriv2;
        }
    }
    EXPECT_NEAR(1, spline.getControlPointValues()[1], TESTTOL);

    // Try using a SplineFitter.
    SplineFitter<Real> fitter = SplineFitter<Real>::fitFromGCV(3, x, coeff);
    const Spline spline2 = fitter.getSpline();
    EXPECT_NEAR(3, spline2.getSplineDegree(), TESTTOL);
}

/**
 * The test works by seeing if the tested splines have the properties of
 * a natural cubic spline. To do so, this test file has several steps
 *
 * User Steps: Configure the script
 *        a. Choose the function to be interpolated
 *        b. Choose the location and number of knot points
 *        c. Choose the density of a high resolution interpolation
 *
 * Test Script Steps:
 * 0. Initialize the input vectors xK, xM, and xD for the knot locations,
 *    mid knot location and high resolution step locations respectively
 * 1. Initialize the analytically computed output yK, yM and yD
 * 2. Create each of the spline objects.
 * 3. Evaluate the numerical accuracy of the splines by calling
 *    testNaturalCubicSpline
 */
TEST(Spline, NaturalCubicSpline) {
    testNaturalCubicSpline(0);
    testNaturalCubicSpline(1);
    testNaturalCubicSpline(2);
    testNaturalCubicSpline(3);
    // testNaturalCubicSpline(4);
    testNaturalCubicSpline(5);
}