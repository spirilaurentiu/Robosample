#include <gtest/gtest.h>
#include <vector>

#include "SimTKcommon.h"
#include "Util.hpp"

using namespace SimTK;

TEST(SimTKCommon_Function_FunctionConstant, ValueAndDerivative) {
    Function_<Vec3>::Constant f(Vec3(1, 2, 3), 2);

    EXPECT_EQ(f.getArgumentSize(), 2);

    Vector x(2);
    EXPECT_TRUE(AssertSimTKEqual("f(x)", "Vec3(1,2,3)", f.calcValue(x), Vec3(1, 2, 3)));

    Array_<int> derivComponents(1);
    const Vec3 df = f.calcDerivative(derivComponents, x);

    EXPECT_TRUE(AssertSimTKEqual("df", "Vec3(0)", df, Vec3(0)));
}

TEST(SimTKCommon_Function_FunctionLinearVec3, EvaluationAndDerivatives) {
    Vector_<Vec3> coeff(3);
    coeff[0] = Vec3(1, 2, 3);
    coeff[1] = Vec3(4, 3, 2);
    coeff[2] = Vec3(-1, -2, -3);

    Function_<Vec3>::Linear f(coeff);

    EXPECT_EQ(f.getArgumentSize(), 2);

    EXPECT_TRUE(AssertSimTKEqual("f(0,0)", "", f.calcValue(Vector(Vec2(0, 0))), Vec3(-1, -2, -3)));
    EXPECT_TRUE(AssertSimTKEqual("f(1,0)", "", f.calcValue(Vector(Vec2(1, 0))), Vec3(0, 0, 0)));
    EXPECT_TRUE(
        AssertSimTKEqual("f(0.5,-0.5)", "", f.calcValue(Vector(Vec2(0.5, -0.5))), Vec3(-2.5, -2.5, -2.5)));

    std::vector<int> d1{1};
    EXPECT_TRUE(AssertSimTKEqual("df/dx1", "", f.calcDerivative(d1, Vector(Vec2(1, 0))), Vec3(4, 3, 2)));

    std::vector<int> d2(2);
    EXPECT_TRUE(AssertSimTKEqual("d2f", "", f.calcDerivative(d2, Vector(Vec2(1, 0))), Vec3(0)));
}

TEST(SimTKCommon_Function_FunctionPolynomial, EvaluationAndDerivatives) {
    Vector_<Vec3> coeff(3);
    coeff[0] = Vec3(1, 2, 3);
    coeff[1] = Vec3(4, 3, 2);
    coeff[2] = Vec3(-1, -2, -3);

    Function_<Vec3>::Polynomial f(coeff);

    EXPECT_EQ(f.getArgumentSize(), 1);

    EXPECT_TRUE(AssertSimTKEqual("f(0)", "", f.calcValue(Vector(Vec1(0))), Vec3(-1, -2, -3)));
    EXPECT_TRUE(AssertSimTKEqual("f(1)", "", f.calcValue(Vector(Vec1(1))), Vec3(4, 3, 2)));
    EXPECT_TRUE(AssertSimTKEqual("f(2)", "", f.calcValue(Vector(Vec1(2))), Vec3(11, 12, 13)));

    std::vector<int> d1{0};
    EXPECT_TRUE(AssertSimTKEqual("f'(0)", "", f.calcDerivative(d1, Vector(Vec1(0))), Vec3(4, 3, 2)));
    EXPECT_TRUE(AssertSimTKEqual("f'(1)", "", f.calcDerivative(d1, Vector(Vec1(1))), Vec3(6, 7, 8)));

    std::vector<int> d2(2);
    EXPECT_TRUE(AssertSimTKEqual("f''", "", f.calcDerivative(d2, Vector(Vec1(0))), Vec3(2, 4, 6)));

    std::vector<int> d3(3);
    EXPECT_TRUE(AssertSimTKEqual("f'''", "", f.calcDerivative(d3, Vector(Vec1(1))), Vec3(0)));
}


TEST(SimTKCommon_Function_FunctionLinearReal, EvaluationAndDerivatives) {
    Vector coeff(3);
    coeff[0] = 1.0;
    coeff[1] = 4.0;
    coeff[2] = -1.0;

    Function::Linear f(coeff);

    EXPECT_EQ(f.getArgumentSize(), 2);

    EXPECT_DOUBLE_EQ(f.calcValue(Vector(Vec2(0, 0))), -1.0);
    EXPECT_DOUBLE_EQ(f.calcValue(Vector(Vec2(1, 0))), 0.0);
    EXPECT_DOUBLE_EQ(f.calcValue(Vector(Vec2(0.5, -0.5))), -2.5);

    Array_<int> d1(1);
    d1[0] = 1;
    EXPECT_DOUBLE_EQ(f.calcDerivative(d1, Vector(Vec2(1, 0))), 4.0);

    Array_<int> d2(2);
    EXPECT_DOUBLE_EQ(f.calcDerivative(d2, Vector(Vec2(1, 0))), 0.0);
}

TEST(SimTKCommon_Function_FunctionSinusoid, HigherOrderDerivatives) {
    const Real a = 11.23;
    const Real w = 1.1;
    const Real p = Pi / 4;

    Function::Sinusoid s(a, w, p);

    Vector t(1, 0.23);

    EXPECT_DOUBLE_EQ(s.calcValue(Vector(1, 0.0)), a * std::sin(p));
    EXPECT_DOUBLE_EQ(s.calcValue(t), a * std::sin((w * t[0]) + p));

    Array_<int> deriv;

    // 0th derivative
    EXPECT_DOUBLE_EQ(s.calcDerivative(deriv, t), s.calcValue(t));

    // 1st derivative
    deriv.push_back(0);
    EXPECT_DOUBLE_EQ(s.calcDerivative(deriv, t), a * w * std::cos((w * t[0]) + p));

    // 2nd derivative
    deriv.push_back(0);
    EXPECT_DOUBLE_EQ(s.calcDerivative(deriv, t), -a * w * w * std::sin((w * t[0]) + p));
}

// TEST(SimTKCommon_Function_FunctionStep, BoundaryAndDerivatives) {
//     Function::Step s(-1, 1, 0, 1);

//     EXPECT_DOUBLE_EQ(s.calcValue(Vector(1, 0)), -1);
//     EXPECT_DOUBLE_EQ(s.calcValue(Vector(1, 1)), 1);
//     EXPECT_DOUBLE_EQ(s.calcValue(Vector(1, 0.5)), 0);

//     Array_<int> d1(1);
//     EXPECT_DOUBLE_EQ(s.calcDerivative(d1, Vector(1, 0)), 0);
//     EXPECT_DOUBLE_EQ(s.calcDerivative(d1, Vector(1, 1)), 0);

//     Array_<int> d2(2);
//     EXPECT_DOUBLE_EQ(s.calcDerivative(d2, Vector(1, 0)), 0);

//     // Vec3 interpolation
//     Function_<Vec3>::Step sv(Vec3(1, 2, 3), Vec3(4, 5, 6), 0, 1);

//     EXPECT_TRUE(AssertSimTKEqual("sv(0.5)", "", sv.calcValue(Vector(1, 0.5)), Vec3(2.5, 3.5, 4.5)));

//     EXPECT_TRUE(AssertSimTKEqual("sv''", "", sv.calcDerivative(d2, Vector(1, -29.3)), Vec3(0)));
// }