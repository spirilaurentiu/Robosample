// Test the functioning of Simbody operators which involve the mass matrix,
// and other system matrices like the Jacobian (partial velocity matrix) that
// maps between generalized and spatial coordinates.
// The O(N) operators like multiplyByM() and multiplyByMInv() are supposed to
// behave *as though* they used the mass matrix, without actually forming it.

#include <exception>
#include <gtest/gtest.h>
#include <iostream>

#include "SimTKcommon/Testing.h"

#include "SimTKsimbody.h"

using namespace SimTK;

#define EXPECT_SIMTK_SIZE(expected, actual, n)                            \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), (n))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nSize:      " << (n) << "\n"

#define EXPECT_NEAR_CUSTOM_TOL_SIMTK(expected, actual, tol)                              \
    EXPECT_TRUE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1, (tol))) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\nTolerance: " << (tol) << "\n"

#define EXPECT_NEAR_DEFAULT_TOL_SIMTK(expected, actual) \
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((expected), (actual), 1e-6)

#define EXPECT_NOT_NEAR_DEFAULT_TOL_SIMTK(expected, actual)                        \
    EXPECT_FALSE(SimTK::Test::numericallyEqual((expected), (actual), /*scale*/ 1)) \
        << "Expected:  " << (expected) << "\nActual:    " << (actual) << "\n"

// This will apply a constant set of mobility forces that can be set
// externally.
class MyForceImpl : public Force::Custom::Implementation {
    public:
    MyForceImpl() = default;

    void calcForce(const State& state,
                   Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces,
                   Vector& mobilityForces) const override {
        EXPECT_TRUE(f.size() == 0 || f.size() == mobilityForces.size());
        EXPECT_TRUE(F.size() == 0 || F.size() == bodyForces.size());
        if (f.size()) {
            mobilityForces += f;
        }
        if (F.size()) {
            bodyForces += F;
        }
    }
    [[nodiscard]] auto calcPotentialEnergy(const State& state) const -> Real override {
        return 0;
    }

    void setMobilityForces(const Vector& mobFrc) {
        f = mobFrc;
    }
    void setBodyForces(const Vector_<SpatialVec>& bodFrc) {
        F = bodFrc;
    }

    private:
    Vector f;
    Vector_<SpatialVec> F;
};

// This is an imitation of SD/FAST's sdrel2cart() subroutine. We
// are given a station point S fixed to a body B. S is given by
// the constant vector p_BS from B's origin Bo to point S, expressed
// in B's frame. Denote the position of S in the ground frame G
// p_GS = p_GB + p_BS_G, where p_BS_G=R_GB*p_BS is the vector p_BS
// reexpressed in G. The velocity of S in G is v_GS = d/dt p_GS, taken
// in G. So v_GS = v_GB + w_GB X p_BS_G = v_GB - p_BS_G % w_GB.
//
// We would like to obtain the partial velocity of S with respect to each of
// the generalized speeds u, taken in the Ground frame, that is,
// JS=d v_GS / du. (JS is a 3xnu matrix, or a single row of Vec3s.) We have
// a method that can calculate J=d V_GB / du where V_GB=[w_GB;v_GB] is the
// spatial velocity of B at its origin. So we need to calculate
//        d v_GS   d v_GS   d V_GB
//   JS = ------ = ------ * ------ =
//          du     d V_GB     du
//
//               = [ -px | eye(3) ] * J
// where px is the cross product matrix of p_BS_G and eye(3) is a 3x3
// identity matrix.
//
// This function should produce the same result as the SimbodyMatterSubsystem
// method calcStationJacobian().
void sbrel2cart(const State& state,
                const SimbodyMatterSubsystem& matter,
                MobilizedBodyIndex bodyIx,
                const Vec3& p_BS,       // point in body frame
                RowVector_<Vec3>& dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of body origin.
    Vector_<SpatialVec> J(nu);
    J = SpatialVec(Vec3(0), Vec3(0)); // or J.setToZero();

    Vector u(nu);
    u = 0;
    Vector_<SpatialVec> Ju(nu); // d allV / d ui
    for (int i = 0; i < nu; ++i) {
        u[i] = 1;
        matter.multiplyBySystemJacobian(state, u, Ju);
        u[i] = 0;
        J[i] = Ju[bodyIx]; // pick out the body of interest
    }

    Row<2, Mat33> dvdV(-crossMat(p_BS_G), Mat33(1));
    dvdu.resize(nu);
    for (int i = 0; i < nu; ++i) {
        dvdu[i] = dvdV * J[i]; // or J[i][0] % p_BS_G + J[i][1]
    }
}

// Another way to calculate exactly what sbrel2cart() does -- can we do it
// faster using f=J^T*F rather than V=J*u?
void sbrel2cart2(const State& state,
                 const SimbodyMatterSubsystem& matter,
                 MobilizedBodyIndex bodyIx,
                 const Vec3& p_BS,       // point in body frame
                 RowVector_<Vec3>& dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of body origin.
    // (This is one row of J.)
    Matrix Jt(nu, 6); // a column of Jt but with scalar elements

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx];               // the only one we'll change
    for (int which = 0; which < 2; ++which) { // moment, force
        for (int i = 0; i < 3; ++i) {
            Fb[which][i] = 1;
            VectorView col = Jt(3 * which + i);
            matter.multiplyBySystemJacobianTranspose(state, F, col);
            Fb[which][i] = 0;
        }
    }

    Row<2, Mat33> dvdV(-crossMat(p_BS_G), Mat33(1));
    dvdu.resize(nu);
    for (int i = 0; i < nu; ++i) {
        const RowVectorView r = Jt[i];
        SpatialVec V(Vec3::getAs(&r[0]), Vec3::getAs(&r[3]));
        dvdu[i] = dvdV * V; // or J[i][0] % p_BS_G + J[i][1]
    }
}

// This is a further refinement that still calculates exactly what sbrel2cart()
// does. But now try it without the intermediate storage for a row of J^T and
// using only 3 J*v multiplies since only the translational (station) Jacobian
// is wanted.
void sbrel2cart3(const State& state,
                 const SimbodyMatterSubsystem& matter,
                 MobilizedBodyIndex bodyIx,
                 const Vec3& p_BS,       // point in body frame
                 RowVector_<Vec3>& dvdu) // v is dS/dt in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dvdu where v is linear velocity of p_BS.
    // (This is three rows of J.)
    dvdu.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx]; // the only one we'll change
    Vector col(nu);             // temporary to hold column of J^T
    for (int i = 0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F
        matter.multiplyBySystemJacobianTranspose(state, F, col);
        for (int r = 0; r < nu; ++r) {
            dvdu[r][i] = col[r];
        }
        Fb[1][i] = 0;
    }
}

// Using the method of sbrel2cart3() but with 6 J*v multiplies, this gives
// the full 6xnu "Frame Jacobian" for one body for a specified frame on that
// body (only the origin, not the orientation, matters).
// This should produce the same result as the built-in calcFrameJacobian()
// method.
void sbrel2cart4(const State& state,
                 const SimbodyMatterSubsystem& matter,
                 MobilizedBodyIndex bodyIx,
                 const Vec3& p_BS,             // point in body frame
                 RowVector_<SpatialVec>& dVdu) // V is [w,v] in G
{
    const int nu = state.getNU();
    const int nb = matter.getNumBodies(); // includes ground

    const MobilizedBody& mobod = matter.getMobilizedBody(bodyIx);
    const Vec3 p_BS_G = mobod.expressVectorInGroundFrame(state, p_BS);

    // Calculate J=dVdu where V is spatial velocity of p_BS.
    // (This is six rows of J.)
    dVdu.resize(nu);

    Vector_<SpatialVec> F(nb, SpatialVec(Vec3(0)));
    SpatialVec& Fb = F[bodyIx]; // the only one we'll change
    Vector col(nu);             // temporary to hold column of J^T
    // Rotational part.
    for (int i = 0; i < 3; ++i) {
        Fb[0][i] = 1;
        matter.multiplyBySystemJacobianTranspose(state, F, col);
        for (int r = 0; r < nu; ++r) {
            dVdu[r][0][i] = col[r];
        }
        Fb[0][i] = 0;
    }
    // Translational part.
    for (int i = 0; i < 3; ++i) {
        Fb[1][i] = 1;
        Fb[0] = p_BS_G % Fb[1]; // r X F
        matter.multiplyBySystemJacobianTranspose(state, F, col);
        for (int r = 0; r < nu; ++r) {
            dVdu[r][1][i] = col[r];
        }
        Fb[1][i] = 0;
    }
}

// Compare two representations of the same matrix: one as an mXn matrix
// of SpatialVecs, the other as a 6mXn matrix of scalars. Note that this
// will also work if the first actual parameter is a Vector_<SpatialVec>
// or RowVector_<SpatialVec> since those have implicit conversions to mX1
// or 1Xn Matrix_<SpatialVec>, resp.
[[nodiscard]] static auto compareElementwise(const Matrix_<SpatialVec>& J, const Matrix& Jf) -> bool {
    const int m = J.nrow();
    const int n = J.ncol();
    EXPECT_EQ(Jf.nrow(), 6 * m);
    EXPECT_EQ(Jf.ncol(), n);

    for (int b = 0; b < m; ++b) {
        const int r = 6 * b; // row start for Jf
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < n; ++j) {
                if (J(b, j)[i / 3][i % 3] != Jf(r + i, j)) {
                    return false;
                }
            }
        }
    }

    return true;
}

// Same thing but for comparing matrices where one has Vec3 elements.
[[nodiscard]] static auto compareElementwise(const Matrix_<Vec3>& JS, const Matrix& JSf) -> bool {
    const int m = JS.nrow();
    const int n = JS.ncol();
    EXPECT_EQ(JSf.nrow(), 3 * m);
    EXPECT_EQ(JSf.ncol(), n);

    for (int b = 0; b < m; ++b) {
        const int r = 3 * b; // row start for JSf
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < n; ++j) {
                if (JS(b, j)[i] != JSf(r + i, j)) {
                    return false;
                }
            }
        }
    }

    return true;
}

TEST(Simbody_MobilizedBody_RelativeToCartesianJacobianConsistency,
     StationAndFrameJacobianAgreementAcrossConfigurations) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);

    // Pendulum of length 1, initially along -x like this:
    //     B ----------- * O
    //    -1,0,0          0,0,0
    // At q=0, partial(B)/partial(u) = 0,-1,0.
    // At q=pi/2, partial(B)/partial(u) = 1,0,0.
    // Then try this with station S=(0,1,0)_B:
    //      S
    //      |
    //      |
    //      B ------ * O
    // Now |OS| = sqrt(2). At q=Pi/4, S will be horizontal
    // so partial(S)/partial(u) = (0, -sqrt(2)/2, 0).
    //
    // In all cases the partial angular velocity is (0,0,1).
    //
    MobilizedBody::Pin pinBody(matter.Ground(),
                               Transform(),
                               MassProperties(1, Vec3(0), Inertia(1)),
                               Vec3(1, 0, 0));
    State state = system.realizeTopology();

    RowVector_<Vec3> dvdu;
    RowVector_<Vec3> dvdu2;
    RowVector_<Vec3> dvdu3;
    RowVector_<Vec3> JS;
    RowVector_<SpatialVec> dvdu4;
    RowVector_<SpatialVec> JF;
    Matrix_<SpatialVec> J;
    Matrix_<SpatialVec> Jn;
    Matrix JSf;
    Matrix JFf;
    Matrix Jf; // flat

    // We'll compute Jacobians for the body origin Bo in 2 configurations,
    // q==0,pi/2 then for station S (0,1,0) at q==pi/4. Not
    // much of a test, I know, but at least we know the right answer.

    // q == 0; answer is JB == (0,-1,0)

    pinBody.setQ(state, 0);
    system.realize(state, Stage::Position);
    sbrel2cart(state, matter, pinBody, Vec3(0), dvdu);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu[0], Vec3(0, -1, 0));

    sbrel2cart2(state, matter, pinBody, Vec3(0), dvdu2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu2[0], Vec3(0, -1, 0));

    sbrel2cart3(state, matter, pinBody, Vec3(0), dvdu3);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu3[0], Vec3(0, -1, 0));

    matter.calcStationJacobian(state, pinBody, Vec3(0), JS);
    matter.calcStationJacobian(state, pinBody, Vec3(0), JSf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JS, dvdu3); // == dvdu2 == dvdu
    EXPECT_TRUE(compareElementwise(JS, JSf));

    sbrel2cart4(state, matter, pinBody, Vec3(0), dvdu4);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][0], Vec3(0, 0, 1));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][1], Vec3(0, -1, 0));

    matter.calcFrameJacobian(state, pinBody, Vec3(0), JF);
    matter.calcFrameJacobian(state, pinBody, Vec3(0), JFf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4, JF);
    EXPECT_TRUE(compareElementwise(JF, JFf));

    // Calculate the whole system Jacobian at q==0 in two different
    // representations and make sure they are the same.
    matter.calcSystemJacobian(state, J);
    matter.calcSystemJacobian(state, Jf);
    EXPECT_TRUE(compareElementwise(J, Jf));

    // q == 90 degrees; answer is JB == (1,0,0)

    pinBody.setQ(state, Pi / 2);
    system.realize(state, Stage::Position);
    sbrel2cart(state, matter, pinBody, Vec3(0), dvdu);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu[0], Vec3(1, 0, 0));

    sbrel2cart2(state, matter, pinBody, Vec3(0), dvdu2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu2[0], Vec3(1, 0, 0));

    sbrel2cart3(state, matter, pinBody, Vec3(0), dvdu3);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu3[0], Vec3(1, 0, 0));

    matter.calcStationJacobian(state, pinBody, Vec3(0), JS);
    matter.calcStationJacobian(state, pinBody, Vec3(0), JSf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JS, dvdu3); // == dvdu2 == dvdu
    EXPECT_TRUE(compareElementwise(JS, JSf));

    sbrel2cart4(state, matter, pinBody, Vec3(0), dvdu4);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][0], Vec3(0, 0, 1));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][1], Vec3(1, 0, 0));

    matter.calcFrameJacobian(state, pinBody, Vec3(0), JF);
    matter.calcFrameJacobian(state, pinBody, Vec3(0), JFf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4, JF);
    EXPECT_TRUE(compareElementwise(JF, JFf));

    // Calculate the whole system Jacobian at q==pi/2 in two different
    // representations and make sure they are the same.
    matter.calcSystemJacobian(state, J);
    matter.calcSystemJacobian(state, Jf);
    EXPECT_TRUE(compareElementwise(J, Jf));

    // now station S, q == 45 degrees; answer is JS == (0,-sqrt(2),0)

    pinBody.setQ(state, Pi / 4);
    system.realize(state, Stage::Position);
    sbrel2cart(state, matter, pinBody, Vec3(0, 1, 0), dvdu);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu[0], Vec3(0, -Sqrt2, 0));

    sbrel2cart2(state, matter, pinBody, Vec3(0, 1, 0), dvdu2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu2[0], Vec3(0, -Sqrt2, 0));

    sbrel2cart3(state, matter, pinBody, Vec3(0, 1, 0), dvdu3);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu3[0], Vec3(0, -Sqrt2, 0));

    matter.calcStationJacobian(state, pinBody, Vec3(0, 1, 0), JS);
    matter.calcStationJacobian(state, pinBody, Vec3(0, 1, 0), JSf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JS, dvdu3); // == dvdu2 == dvdu
    EXPECT_TRUE(compareElementwise(JS, JSf));

    // Calculate station Jacobian JS by multiplication to test that the
    // multiplyByStationJacobian[Transpose] methods are working.
    Vec3 JSn; // 3xnu
    Vector u(1);
    u[0] = 1.;
    JSn = matter.multiplyByStationJacobian(state, pinBody, Vec3(0, 1, 0), u);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JSn, dvdu3[0]);

    Vec3 FS(0);
    Row3 JSnt;
    for (int i = 0; i < 3; ++i) {
        FS[i] = 1;
        matter.multiplyByStationJacobianTranspose(state, pinBody, Vec3(0, 1, 0), FS, u);
        FS[i] = 0;
        JSnt[i] = u[0];
    }
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JSnt, ~dvdu3[0]);

    sbrel2cart4(state, matter, pinBody, Vec3(0, 1, 0), dvdu4);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][0], Vec3(0, 0, 1));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4[0][1], Vec3(0, -Sqrt2, 0));

    matter.calcFrameJacobian(state, pinBody, Vec3(0, 1, 0), JF);
    matter.calcFrameJacobian(state, pinBody, Vec3(0, 1, 0), JFf);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(dvdu4, JF);
    EXPECT_TRUE(compareElementwise(JF, JFf));

    // Calculate frame Jacobian JF by multiplication to test that the
    // multiplyByFrameJacobian[Transpose] methods are working.
    SpatialVec JFn; // 6xnu
    u[0] = 1.;
    JFn = matter.multiplyByFrameJacobian(state, pinBody, Vec3(0, 1, 0), u);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JFn, dvdu4[0]);

    SpatialVec FF(Vec3(0));
    SpatialRow JFnt;
    for (int i = 0; i < 6; ++i) {
        FF[i / 3][i % 3] = 1;
        matter.multiplyByFrameJacobianTranspose(state, pinBody, Vec3(0, 1, 0), FF, u);
        FF[i / 3][i % 3] = 0;
        JFnt[i / 3][i % 3] = u[0];
    }
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JFnt, ~dvdu4[0]);

    // Calculate the whole system Jacobian at q==pi/4 in two different
    // representations and make sure they are the same.
    matter.calcSystemJacobian(state, J);
    matter.calcSystemJacobian(state, Jf);
    EXPECT_TRUE(compareElementwise(J, Jf));

    // Generate the whole system Jacobian one body at a time by multiplication.
    Jn.resize(matter.getNumBodies(), matter.getNumMobilities());
    u[0] = 1;
    for (MobodIndex i(0); i < matter.getNumBodies(); ++i) {
        Jn[i] = matter.multiplyByFrameJacobian(state, i, Vec3(0), u);
    }
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(Jn, J);
}

void makeSystem(bool constrained, MultibodySystem& mbs, MyForceImpl*& frcp) {
    SimbodyMatterSubsystem pend(mbs);
    GeneralForceSubsystem forces(mbs);
    frcp = new MyForceImpl();
    Force::Custom(forces, frcp);

    const Real randomAngle1 = (Pi / 2) * Test::randReal();
    const Real randomAngle2 = (Pi / 2) * Test::randReal();
    Vector_<Vec3> randomVecs(10);
    for (int i = 0; i < 10; ++i) {
        randomVecs[i] = Test::randVec3();
    }

    const Real mass = 2.3;
    const Vec3 com = randomVecs[5];
    const Inertia inertia = Inertia(3, 4, 5, .01, -.02, .04).shiftFromMassCenter(com, mass);
    Body::Rigid pendulumBody = Body::Rigid(MassProperties(mass, com, inertia));

    MobilizedBody::Ball pendBody1(pend.Ground(),
                                  Transform(Rotation(randomAngle1, randomVecs[0]), randomVecs[1]),
                                  pendulumBody,
                                  Transform(Rotation(randomAngle2, randomVecs[2]), randomVecs[3]));
    MobilizedBody::Weld pendBody2(pendBody1,
                                  Transform(Rotation(randomAngle1, randomVecs[4]), randomVecs[5]),
                                  pendulumBody,
                                  Transform(Rotation(randomAngle2, randomVecs[6]), randomVecs[7]));

    MobilizedBody::Pin pendBody3(pendBody2,
                                 Transform(Rotation(randomAngle1, randomVecs[8]), randomVecs[9]),
                                 pendulumBody,
                                 Transform(Rotation(randomAngle2, randomVecs[8]), randomVecs[7]));
    MobilizedBody::Screw pendBody4(pendBody3,
                                   Transform(Rotation(randomAngle2, randomVecs[6]), randomVecs[5]),
                                   pendulumBody,
                                   Transform(Rotation(randomAngle1, randomVecs[4]), randomVecs[3]),
                                   3); // pitch
    MobilizedBody::Translation pendBody5(pendBody4,
                                         Transform(Rotation(randomAngle2, randomVecs[2]), randomVecs[1]),
                                         pendulumBody,
                                         Transform(Rotation(randomAngle1, randomVecs[0]), randomVecs[1]));

    // Now add some side branches.
    MobilizedBody::BendStretch pendBody1a(pendBody1,
                                          Transform(Rotation(randomAngle2, randomVecs[2]), randomVecs[3]),
                                          pendulumBody,
                                          Transform(Rotation(randomAngle1, randomVecs[4]), randomVecs[5]));

    MobilizedBody::Slider pendBody2a(pendBody2,
                                     Transform(Rotation(randomAngle1, randomVecs[6]), randomVecs[7]),
                                     pendulumBody,
                                     Transform(Rotation(randomAngle2, randomVecs[8]), randomVecs[9]));

    MobilizedBody::Universal pendBody2b(pendBody2a,
                                        Transform(Rotation(randomAngle1, randomVecs[8]), randomVecs[7]),
                                        pendulumBody,
                                        Transform(Rotation(randomAngle2, randomVecs[6]), randomVecs[5]));
    MobilizedBody::Slider pendBody2x(pendBody2b,
                                     Transform(Rotation(randomAngle1, randomVecs[6]), randomVecs[7]),
                                     pendulumBody,
                                     Transform(Rotation(randomAngle2, randomVecs[8]), randomVecs[9]));

    MobilizedBody::Planar pendBody4a(pendBody4,
                                     Transform(Rotation(randomAngle1, randomVecs[4]), randomVecs[3]),
                                     pendulumBody,
                                     Transform(Rotation(randomAngle2, randomVecs[2]), randomVecs[1]));

    // Probably can't be satisfied, but doesn't matter
    if (constrained) {
        // Holonomic
        Constraint::Rod(pendBody4, pendBody2b, 1.);

        // Nonholonomic
        Constraint::ConstantSpeed(pendBody2a, MobilizerUIndex(0), -3.);

        // Acceleration only
        Constraint::ConstantAcceleration(pendBody5, MobilizerUIndex(2), 0.01);

        // Weld
        Constraint::Weld(pendBody4a, Test::randTransform(), pendBody4, Test::randTransform());
    }
}


// Test calculations of Jacobian "bias" terms, where bias=JDot*u.
// We can estimate JDot using a numerical directional derivative
// since JDot = (DJ/Dq)*qdot ~= (J(q+h*qdot)-J(q-h*qdot))/2h.
// Then we multiply JDot*u and compare with the bias calculations.
// Or, we can estimate JDot*u directly with
//       JDotu ~= (J(q+h*qdot)*u - J(q-h*qdot)*u)/2h
// using the fast "multiply by Jacobian" methods.
// We use both methods below.
TEST(Simbody_MobilizedBody_JacobianBiasTerms, JDotUConsistencyWithAnalyticBiasAndFiniteDifference) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    system.realizeModel(state);

    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);

    const MobilizedBodyIndex whichBod(8);
    const Vec3 whichPt(1, 2, 3);
    system.realize(state, Stage::Velocity);
    const Vector& q = state.getQ();
    const Vector& u = state.getU();
    const Vector& qdot = state.getQDot();

    // sbias, fbias, sysbias are the JDot*u quantities we want to check.
    const Vec3 sbias = matter.calcBiasForStationJacobian(state, whichBod, whichPt);
    const SpatialVec fbias = matter.calcBiasForFrameJacobian(state, whichBod, whichPt);
    Vector_<SpatialVec> sysbias;
    matter.calcBiasForSystemJacobian(state, sysbias);

    // These are for computing JDot first.
    RowVector_<Vec3> JS_P;
    RowVector_<Vec3> JS1_P;
    RowVector_<Vec3> JS2_P;
    RowVector_<Vec3> JSDot_P;
    RowVector_<SpatialVec> JF_P;
    RowVector_<SpatialVec> JF1_P;
    RowVector_<SpatialVec> JF2_P;
    RowVector_<SpatialVec> JFDot_P;
    Matrix_<SpatialVec> J;
    Matrix_<SpatialVec> J1;
    Matrix_<SpatialVec> J2;
    Matrix_<SpatialVec> JDot;

    // These are for computing JDot*u directly.
    Vec3 JS_Pu;
    Vec3 JS1_Pu;
    Vec3 JS2_Pu;
    Vec3 JSDot_Pu;
    SpatialVec JF_Pu;
    SpatialVec JF1_Pu;
    SpatialVec JF2_Pu;
    SpatialVec JFDot_Pu;
    Vector_<SpatialVec> Ju;
    Vector_<SpatialVec> J1u;
    Vector_<SpatialVec> J2u;
    Vector_<SpatialVec> JDotu;

    // Unperturbed:
    matter.calcStationJacobian(state, whichBod, whichPt, JS_P);
    matter.calcFrameJacobian(state, whichBod, whichPt, JF_P);
    matter.calcSystemJacobian(state, J);

    JS_Pu = matter.multiplyByStationJacobian(state, whichBod, whichPt, u);
    JF_Pu = matter.multiplyByFrameJacobian(state, whichBod, whichPt, u);
    matter.multiplyBySystemJacobian(state, u, Ju);

    const Real Delta = 5e-6; // we'll use central difference
    State perturbq = state;

    // Perturbed +:
    perturbq.updQ() = q + Delta * qdot;
    system.realize(perturbq, Stage::Position);
    matter.calcStationJacobian(perturbq, whichBod, whichPt, JS2_P);
    matter.calcFrameJacobian(perturbq, whichBod, whichPt, JF2_P);
    matter.calcSystemJacobian(perturbq, J2);

    JS2_Pu = matter.multiplyByStationJacobian(perturbq, whichBod, whichPt, u);
    JF2_Pu = matter.multiplyByFrameJacobian(perturbq, whichBod, whichPt, u);
    matter.multiplyBySystemJacobian(perturbq, u, J2u);

    // Perturbed -:
    perturbq.updQ() = q - Delta * qdot;
    system.realize(perturbq, Stage::Position);
    matter.calcStationJacobian(perturbq, whichBod, whichPt, JS1_P);
    matter.calcFrameJacobian(perturbq, whichBod, whichPt, JF1_P);
    matter.calcSystemJacobian(perturbq, J1);

    JS1_Pu = matter.multiplyByStationJacobian(perturbq, whichBod, whichPt, u);
    JF1_Pu = matter.multiplyByFrameJacobian(perturbq, whichBod, whichPt, u);
    matter.multiplyBySystemJacobian(perturbq, u, J1u);

    // Estimate JDots:
    JSDot_P = (JS2_P - JS1_P) / Delta / 2;
    JFDot_P = (JF2_P - JF1_P) / Delta / 2;
    JDot = (J2 - J1) / Delta / 2;

    // Estimate JDotus:
    JSDot_Pu = (JS2_Pu - JS1_Pu) / Delta / 2;
    JFDot_Pu = (JF2_Pu - JF1_Pu) / Delta / 2;
    JDotu = (J2u - J1u) / Delta / 2;

    // Calculate errors in JDot*u:
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JSDot_P * u - sbias).norm(), 0, SqrtEps);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JFDot_P * u - fbias).norm(), 0, SqrtEps);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JDot * u - sysbias).norm(), 0, SqrtEps);

    // Calculate errors in JDotu:
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JSDot_Pu - sbias).norm(), 0, SqrtEps);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JFDot_Pu - fbias).norm(), 0, SqrtEps);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK((JDotu - sysbias).norm(), 0, SqrtEps);
}

TEST(Simbody_MultibodySystem_UnconstrainedDynamicsOperators,
     MassMatrixInverseDynamicsAndResidualForceConsistency) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    // Attainable accuracy drops with problem size.
    const Real Slop = nu * SignificantReal;

    system.realizeModel(state);
    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);

    Vector randVec = 100 * SimTK::Test::randVector(nu);
    Vector result1;
    Vector result2;

    // result1 = M*v
    system.realize(state, Stage::Position);
    matter.multiplyByM(state, randVec, result1);
    EXPECT_EQ(result1.size(), nu);

    // result2 = M^-1 * result1 == M^-1 * M * v == v
    system.realize(state, Stage::Dynamics);
    matter.multiplyByMInv(state, result1, result2);
    EXPECT_EQ(result2.size(), nu);

    EXPECT_NEAR_CUSTOM_TOL_SIMTK(result2, randVec, Slop);

    Matrix M(nu, nu);
    Matrix MInv(nu, nu);

    Vector v(nu, Real(0));
    for (int j = 0; j < nu; ++j) {
        v[j] = 1;
        matter.multiplyByM(state, v, M(j));
        matter.multiplyByMInv(state, v, MInv(j));
        v[j] = 0;
    }

    Matrix MInvCalc(M);
    MInvCalc.invertInPlace();
    EXPECT_SIMTK_SIZE(MInv, MInvCalc, nu);

    Matrix identity(nu, nu);
    identity = 1;
    EXPECT_SIMTK_SIZE(M * MInv, identity, nu);
    EXPECT_SIMTK_SIZE(MInv * M, identity, nu);

    // Compare above-calculated values with values returned by the calcM() and calcMInv() methods.
    Matrix MM;
    Matrix MMInv;

    matter.calcM(state, MM);
    matter.calcMInv(state, MMInv);

    EXPECT_SIMTK_SIZE(MM, M, nu);
    EXPECT_SIMTK_SIZE(MMInv, MInv, nu);

    frcp->setMobilityForces(randVec);
    system.realize(state, Stage::Acceleration);
    Vector accel = state.getUDot();

    matter.multiplyByMInv(state, randVec, result1);
    EXPECT_NOT_NEAR_DEFAULT_TOL_SIMTK(accel, result1);
    EXPECT_GT((accel - result1).norm(), SignificantReal);

    // With no velocities M^-1*f should match calculated acceleration.
    state.updU() = 0;
    system.realize(state, Stage::Acceleration);
    accel = state.getUDot();
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(accel, result1);

    // And then M*a should = f.
    matter.multiplyByM(state, accel, result2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(result2, randVec, Slop);

    // Test forward and inverse dynamics operators.
    // Apply random forces and a random prescribed acceleration to
    // get back the residual generalized forces. Then applying those
    // should result in zero residual, and applying them.

    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);

    // Inverse dynamics should require realization only to Velocity stage.
    system.realize(state, Stage::Velocity);

    // Randomize body forces.
    Vector_<SpatialVec> bodyForces(nb);
    for (int i = 0; i < nb; ++i) {
        bodyForces[i] = SimTK::Test::randSpatialVec();
    }

    // Random mobility forces and known udots.
    const Vector mobilityForces = SimTK::Test::randVector(nu);
    const Vector knownUdots = SimTK::Test::randVector(nu);

    // Check self consistency: compute residual, apply it, should be no remaining residual.
    Vector residualForces;
    Vector shouldBeZeroResidualForces;

    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                bodyForces,
                                                knownUdots,
                                                residualForces);
    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces + residualForces,
                                                bodyForces,
                                                knownUdots,
                                                shouldBeZeroResidualForces);
    EXPECT_LE(shouldBeZeroResidualForces.norm(), Slop);

    // Now apply these forces in forward dynamics and see if we get the desired
    // acceleration. State must be realized to Dynamics stage.
    system.realize(state, Stage::Dynamics);
    Vector udots;
    Vector_<SpatialVec> bodyAccels;
    matter.calcAccelerationIgnoringConstraints(state,
                                               mobilityForces + residualForces,
                                               bodyForces,
                                               udots,
                                               bodyAccels);

    EXPECT_NEAR_CUSTOM_TOL_SIMTK(udots, knownUdots, Slop);

    // See if we get back the same body accelerations by feeding in
    // these udots.
    Vector_<SpatialVec> A_GB;
    Vector_<SpatialVec> AC_GB;
    matter.calcBodyAccelerationFromUDot(state, udots, A_GB);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(A_GB, bodyAccels, Slop);

    // Collect coriolis accelerations.
    AC_GB.resize(matter.getNumBodies());
    for (MobodIndex i(0); i < nb; ++i) {
        AC_GB[i] = matter.getTotalCoriolisAcceleration(state, i);
    }

    // Verify that either a zero-length or all-zero udot gives just
    // coriolis accelerations.
    matter.calcBodyAccelerationFromUDot(state, Vector(), A_GB);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(A_GB, AC_GB, Slop);

    Vector allZeroUdot(matter.getNumMobilities(), Real(0));
    matter.calcBodyAccelerationFromUDot(state, allZeroUdot, A_GB);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(A_GB, AC_GB, Slop);

    // Now let's test noncontiguous input and output vectors.
    Matrix MatUdot(3, nu); // use middle row
    MatUdot.setToNaN();
    MatUdot[1] = ~udots;
    Matrix_<SpatialRow> MatA_GB(3, nb); // use middle row
    MatA_GB.setToNaN();
    matter.calcBodyAccelerationFromUDot(state, ~MatUdot[1], ~MatA_GB[1]);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(MatA_GB[1], ~bodyAccels, Slop);

    // Verify that leaving out arguments makes them act like zeroes.
    Vector residualForces1;
    Vector residualForces2;
    matter.calcResidualForceIgnoringConstraints(state,
                                                0 * mobilityForces,
                                                0 * bodyForces,
                                                0 * knownUdots,
                                                residualForces1);
    // no, the residual is not zero here because of the angular velocities
    matter.calcResidualForceIgnoringConstraints(state,
                                                Vector(),
                                                Vector_<SpatialVec>(),
                                                Vector(),
                                                residualForces2);

    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    // We just calculated f_residual = M udot + f_inertial - f_applied, with
    // both udot and f_applied zero, i.e. f_residual=f_inertial. That should
    // be the same as what is returned by getTotalCentrifugalForces().
    Vector_<SpatialVec> F_inertial(nb);
    Vector f_inertial;
    for (MobodIndex i(0); i < nb; ++i) {
        F_inertial[i] = matter.getTotalCentrifugalForces(state, i);
    }
    matter.multiplyBySystemJacobianTranspose(state, F_inertial, f_inertial);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(f_inertial, residualForces1, Slop);

    // This should also match total Mass*Coriolis acceleration + gyro force.
    Vector_<SpatialVec> F_coriolis(nb);
    Vector_<SpatialVec> F_gyro(nb);
    Vector_<SpatialVec> F_total(nb);
    Vector f_total;

    for (MobodIndex i(0); i < nb; ++i) {
        if (i == 0) {
            F_coriolis[i] = SpatialVec(Vec3(0), Vec3(0));
        } else {
            F_coriolis[i] = matter.getMobilizedBody(i).getBodySpatialInertiaInGround(state) * AC_GB[i];
        }
        F_gyro[i] = matter.getGyroscopicForce(state, i);
    }

    F_total = F_coriolis + F_gyro;
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(F_inertial, F_total, Slop);

    // Same, but leave out combinations of arguments.
    matter.calcResidualForceIgnoringConstraints(state,
                                                0 * mobilityForces,
                                                bodyForces,
                                                knownUdots,
                                                residualForces1);
    matter.calcResidualForceIgnoringConstraints(state, Vector(), bodyForces, knownUdots, residualForces2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                0 * bodyForces,
                                                knownUdots,
                                                residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                Vector_<SpatialVec>(),
                                                knownUdots,
                                                residualForces2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                bodyForces,
                                                0 * knownUdots,
                                                residualForces1);
    matter.calcResidualForceIgnoringConstraints(state, mobilityForces, bodyForces, Vector(), residualForces2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    matter.calcResidualForceIgnoringConstraints(state,
                                                0 * mobilityForces,
                                                bodyForces,
                                                0 * knownUdots,
                                                residualForces1);
    matter.calcResidualForceIgnoringConstraints(state, Vector(), bodyForces, Vector(), residualForces2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                0 * bodyForces,
                                                0 * knownUdots,
                                                residualForces1);
    matter.calcResidualForceIgnoringConstraints(state,
                                                mobilityForces,
                                                Vector_<SpatialVec>(),
                                                Vector(),
                                                residualForces2);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residualForces2, residualForces1, Slop);

    // Check that we object to wrong-length arguments.
    EXPECT_THROW(matter.calcResidualForceIgnoringConstraints(state,
                                                             Vector(3, Zero),
                                                             bodyForces,
                                                             knownUdots,
                                                             residualForces2),
                 std::exception);
    EXPECT_THROW(matter.calcResidualForceIgnoringConstraints(state,
                                                             mobilityForces,
                                                             Vector_<SpatialVec>(5),
                                                             knownUdots,
                                                             residualForces2),
                 std::exception);
    EXPECT_THROW(matter.calcResidualForceIgnoringConstraints(state,
                                                             mobilityForces,
                                                             bodyForces,
                                                             Vector(2),
                                                             residualForces2),
                 std::exception);
}

TEST(Simbody_MultibodySystem_ConstrainedDynamicsOperators,
     ResidualForceAndConstraintConsistencyAcrossRandomStatesAndForces) {
    MultibodySystem mbs;
    MyForceImpl* frcp;
    makeSystem(true, mbs, frcp);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Instance); // allocate multipliers, etc.

    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int m = state.getNMultipliers();
    const int nb = matter.getNumBodies();

    // Attainable accuracy drops with problem size.
    const Real Slop = nu * SignificantReal;

    mbs.realizeModel(state);

    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);

    Vector randMobFrc = 100 * SimTK::Test::randVector(nu);
    Vector_<SpatialVec> randBodyFrc(nb);
    for (int i = 0; i < nb; ++i) {
        randBodyFrc[i] = SimTK::Test::randSpatialVec();
    }

    // Apply random mobility forces
    frcp->setMobilityForces(randMobFrc);

    mbs.realize(state); // calculate accelerations and multipliers
    Vector udot = state.getUDot();
    Vector lambda = state.getMultipliers();
    Vector residual;
    matter.calcResidualForce(state, randMobFrc, Vector_<SpatialVec>(), udot, lambda, residual);

    // Residual should be zero since we accounted for everything.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residual, 0 * randMobFrc, Slop);

    Vector abias;
    Vector mgbias;

    // These are the acceleration error bias terms.
    matter.calcBiasForAccelerationConstraints(state, abias);

    // These use pverr (velocity-level errors) for holonomic constraints.
    matter.calcBiasForMultiplyByG(state, mgbias);

    Vector mgGudot;
    matter.multiplyByG(state, udot, mgbias, mgGudot);

    Matrix G;
    matter.calcG(state, G);

    Vector Gudot = G * udot;
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(mgGudot, Gudot, Slop);

    // Won't be zero because bad constraints
    Vector aerr = state.getUDotErr();
    Vector GudotPlusBias = Gudot + abias;
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(GudotPlusBias, aerr, Slop);

    // Add in some body forces
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    frcp->setBodyForces(randBodyFrc);
    mbs.realize(state);
    udot = state.getUDot();
    lambda = state.getMultipliers();
    matter.calcResidualForce(state, randMobFrc, randBodyFrc, udot, lambda, residual);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residual, 0 * randMobFrc, Slop);

    // Try body forces only.
    state.invalidateAllCacheAtOrAbove(Stage::Dynamics);
    frcp->setMobilityForces(0 * randMobFrc);
    mbs.realize(state);
    udot = state.getUDot();
    lambda = state.getMultipliers();
    matter.calcResidualForce(state, Vector(), randBodyFrc, udot, lambda, residual);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residual, 0 * randMobFrc, Slop);

    // Put vectors in noncontiguous storage.
    Matrix udotmat(3, nu); // rows are noncontig
    Matrix mobFrcMat(11, nu);
    Matrix lambdamat(5, m);
    Matrix_<SpatialRow> bodyFrcMat(3, nb);
    udotmat[2] = ~udot;
    lambdamat[3] = ~lambda;
    mobFrcMat[8] = ~randMobFrc;
    bodyFrcMat[2] = ~randBodyFrc;
    Matrix residmat(4, nu);

    // We last computed udot,lambda with no mobility forces. This time
    // will throw some in and then make sure the residual tries to cancel them.
    matter.calcResidualForce(state, ~mobFrcMat[8], ~bodyFrcMat[2], ~udotmat[2], ~lambdamat[3], ~residmat[2]);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(residmat[2], -1 * mobFrcMat[8], Slop);
}


TEST(Simbody_SimbodyMatterSubsystem_CompositeBodyInertia, InertiaProjectionAndCacheInvalidationConsistency) {
    MultibodySystem mbs;
    SimbodyMatterSubsystem pend(mbs);

    // This will be a point mass but the body origin is not at the COM.
    // That's important to make sure the rigid body shifts are working right.
    Body::Rigid pointMass(MassProperties(3, Vec3(3, 0, 0), UnitInertia(0, 9, 9)));

    // Point mass at x=4.5 with origin at x=1.5 rotating about (0,0,0).
    MobilizedBody::Pin body1(pend.Ground(), Transform(), pointMass, Vec3(-1.5, 0, 0));
    const MobilizedBodyIndex body1x = body1.getMobilizedBodyIndex();

    // A second body 2 units further along x, rotating about the
    // first point mass origin. So this one's origin is at x=3.5 and its
    // mass is at x=6.5.
    MobilizedBody::Pin body2(body1, Transform(), pointMass, Vec3(-2, 0, 0));
    const MobilizedBodyIndex body2x = body2.getMobilizedBodyIndex();

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Position);

    Array_<SpatialInertia, MobilizedBodyIndex> R(pend.getNumBodies());
    pend.calcCompositeBodyInertias(state, R);

    // Calculate expected inertias about the joint axes.
    const Real expInertia2 = body2.getBodyMassProperties(state).getMass() * square(2 + 3);
    const Real expInertia1 = (body1.getBodyMassProperties(state).getMass() * square(1.5 + 3))
                             + (body2.getBodyMassProperties(state).getMass() * square(3.5 + 3));

    // Should be able to recover these inertias by projecting the composite
    // body inertias onto the joint axes using H matrices.
    const SpatialVec H1 = body1.getHCol(state, MobilizerUIndex(0));
    const SpatialVec H2 = body2.getHCol(state, MobilizerUIndex(0));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(~H2 * (R[body2x] * H2), expInertia2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(~H1 * (R[body1x] * H1), expInertia1);

    // This should force realization of the composite body inertias.
    SpatialInertia cbi = pend.getCompositeBodyInertia(state, body1);
    body2.setAngle(state, Pi / 4);

    // This is not allowed until PositionKinematics stage.
    EXPECT_THROW(pend.getCompositeBodyInertia(state, body1), std::exception);
    mbs.realize(state, Stage::Position);

    // Now it should be OK.
    cbi = pend.getCompositeBodyInertia(state, body1);
    body2.setAngle(state, Pi / 5);
    pend.realizePositionKinematics(state);
    EXPECT_EQ(state.getSystemStage(), Stage::Time);
    EXPECT_FALSE(pend.isCompositeBodyInertiasRealized(state));

    pend.realizeCompositeBodyInertias(state);
    EXPECT_TRUE(pend.isCompositeBodyInertiasRealized(state));

    pend.invalidateCompositeBodyInertias(state);
    EXPECT_FALSE(pend.isCompositeBodyInertiasRealized(state));
}

// Currently just testing validity/invalidation, not correctness.
TEST(Simbody_SimbodyMatterSubsystem_ArticulatedBodyInertia, CacheRealizationInvalidationAndStageDependence) {
    MultibodySystem mbs;
    MyForceImpl* frcp;
    makeSystem(true, mbs, frcp);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    const MobodIndex body1(1);
    const MobodIndex body2(2);

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Position);

    // This should force realization of the articulated body inertias.
    ArticulatedInertia abi = matter.getArticulatedBodyInertia(state, body1);
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));

    matter.invalidateArticulatedBodyInertias(state);
    EXPECT_FALSE(matter.isArticulatedBodyInertiasRealized(state));

    matter.realizeArticulatedBodyInertias(state);
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));

    state.updQ() = 0.1;
    EXPECT_EQ(state.getSystemStage(), Stage::Time);
    EXPECT_FALSE(matter.isArticulatedBodyInertiasRealized(state));

    // This is not allowed until PositionKinematics stage.
    EXPECT_THROW(matter.getArticulatedBodyInertia(state, body2), std::exception);
    mbs.realize(state, Stage::Position);

    // Now it should be OK.
    abi = matter.getArticulatedBodyInertia(state, body1);

    matter.invalidatePositionKinematics(state); // a prerequisite
    EXPECT_FALSE(matter.isArticulatedBodyInertiasRealized(state));

    matter.realizePositionKinematics(state);
    matter.getArticulatedBodyInertia(state, body2); // ok; implicit realization
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));
    EXPECT_EQ(state.getSystemStage(), Stage::Time);

    state.updQ() = 0.2;
    EXPECT_FALSE(matter.isArticulatedBodyInertiasRealized(state));

    mbs.realize(state, Stage::Dynamics); // not high enough
    EXPECT_FALSE(matter.isArticulatedBodyInertiasRealized(state));

    mbs.realize(state, Stage::Acceleration); // implicit realization of abis
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));

    state.updTime() = 1.; // shouldn't affect time-independent stuff
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);
    EXPECT_TRUE(matter.isPositionKinematicsRealized(state));
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));
}

// Currently just testing validity/invalidation, not correctness.
TEST(Simbody_SimbodyMatterSubsystem_ArticulatedBodyVelocity,
     CacheRealizationDependenciesAndInvalidationBehavior) {
    MultibodySystem mbs;
    MyForceImpl* frcp;
    makeSystem(true, mbs, frcp);
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    const MobodIndex body1(1);
    const MobodIndex body2(2);

    State state = mbs.realizeTopology();
    mbs.realize(state, Stage::Position);

    // Neither ABIs nor VelocityKinematics valid.
    EXPECT_THROW(matter.realizeArticulatedBodyVelocity(state), std::exception);

    mbs.realize(state, Stage::Velocity);
    // ABIs still not valid.
    EXPECT_THROW(matter.realizeArticulatedBodyVelocity(state), std::exception);

    matter.realizeArticulatedBodyInertias(state);
    matter.realizeArticulatedBodyVelocity(state); // OK now

    state.updU() = 0.1; // invalidates VelocityKinematics
    EXPECT_EQ(state.getSystemStage(), Stage::Position);
    EXPECT_FALSE(matter.isArticulatedBodyVelocityRealized(state));

    mbs.realize(state, Stage::Dynamics); // not enough
    EXPECT_FALSE(matter.isArticulatedBodyVelocityRealized(state));

    mbs.realize(state, Stage::Acceleration); // that should do it!
    EXPECT_TRUE(matter.isArticulatedBodyVelocityRealized(state));

    matter.invalidateArticulatedBodyInertias(state); // a prerequisite
    EXPECT_FALSE(matter.isArticulatedBodyVelocityRealized(state));

    mbs.realize(state, Stage::Acceleration);
    EXPECT_TRUE(matter.isArticulatedBodyVelocityRealized(state));

    matter.invalidateVelocityKinematics(state); // another prerequisite
    EXPECT_FALSE(matter.isArticulatedBodyVelocityRealized(state));

    mbs.realize(state, Stage::Acceleration);
    EXPECT_TRUE(matter.isArticulatedBodyVelocityRealized(state));

    state.updTime() = 1.; // shouldn't affect time-independent stuff
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);
    EXPECT_TRUE(matter.isPositionKinematicsRealized(state));
    EXPECT_TRUE(matter.isVelocityKinematicsRealized(state));
    EXPECT_TRUE(matter.isArticulatedBodyInertiasRealized(state));
    EXPECT_TRUE(matter.isArticulatedBodyVelocityRealized(state));
}

TEST(Simbody_SimbodyMatterSubsystem_TaskJacobians_BlockAndScalarConsistency,
     CorrectnessAcrossRepresentations) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    // Attainable accuracy drops with problem size.
    const Real Slop = nu * SignificantReal;

    system.realizeModel(state);
    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);

    system.realize(state, Stage::Position);

    Matrix_<SpatialVec> J;
    Matrix Jmat, Jmat2;
    matter.calcSystemJacobian(state, J);
    EXPECT_EQ(J.nrow(), nb);
    EXPECT_EQ(J.ncol(), nu);

    matter.calcSystemJacobian(state, Jmat);
    EXPECT_EQ(Jmat.nrow(), 6 * nb);
    EXPECT_EQ(Jmat.ncol(), nu);

    // Unpack J into Jmat2 and compare with Jmat.
    Jmat2.resize(6 * nb, nu);
    for (int row = 0; row < nb; ++row) {
        const int nxtr = 6 * row; // row index into scalar matrix
        for (int col = 0; col < nu; ++col) {
            for (int k = 0; k < 3; ++k) {
                Jmat2(nxtr + k, col) = J(row, col)[0][k];
                Jmat2(nxtr + 3 + k, col) = J(row, col)[1][k];
            }
        }
    }
    // These should be exactly the same.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(Jmat2, Jmat, SignificantReal);

    Vector randU = 100. * SimTK::Test::randVector(nu);
    Vector resultU1;
    Vector resultU2;
    Vector_<SpatialVec> randF(nb);
    Vector_<SpatialVec> resultF1;
    Vector_<SpatialVec> resultF2;
    for (int i = 0; i < nb; ++i) {
        randF[i] = 100. * SimTK::Test::randSpatialVec();
    }

    matter.multiplyBySystemJacobian(state, randU, resultF1);
    resultF2 = J * randU;
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(resultF1, resultF2, Slop);

    matter.multiplyBySystemJacobianTranspose(state, randF, resultU1);
    resultU2 = ~J * randF;
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(resultU1, resultU2, Slop);

    // See if Station Jacobian can be used to duplicate the translation
    // rows of the System Jacobian, and if Frame Jacobian can be used to
    // duplicate the whole thing.
    Array_<MobilizedBodyIndex> allBodies(nb);
    for (int i = 0; i < nb; ++i) {
        allBodies[i] = MobilizedBodyIndex(i);
    }
    Array_<Vec3> allOrigins(nb, Vec3(0));

    Matrix_<Vec3> JS;
    Matrix_<Vec3> JS2;
    Matrix_<Vec3> JSbyrow;
    Matrix_<SpatialVec> JF;
    Matrix_<SpatialVec> JF2;
    Matrix_<SpatialVec> JFbyrow;

    matter.calcStationJacobian(state, allBodies, allOrigins, JS);
    matter.calcFrameJacobian(state, allBodies, allOrigins, JF);
    for (int i = 0; i < nb; ++i) {
        for (int j = 0; j < nu; ++j) {
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(JS(i, j), J(i, j)[1]);
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(JF(i, j), J(i, j));
        }
    }

    // Now use random stations to calculate JS & JF.
    Array_<Vec3> randS(nb);
    for (int i = 0; i < nb; ++i) {
        randS[i] = 10. * SimTK::Test::randVec3();
    }
    matter.calcStationJacobian(state, allBodies, randS, JS);
    matter.calcFrameJacobian(state, allBodies, randS, JF);

    // Recalculate one row at a time to test non-contiguous memory handling.
    // Do it backwards just to show off.
    JSbyrow.resize(nb, nu);
    JFbyrow.resize(nb, nu);
    for (int i = nb - 1; i >= 0; --i) {
        matter.calcStationJacobian(state, allBodies[i], randS[i], JSbyrow[i]);
        matter.calcFrameJacobian(state, allBodies[i], randS[i], JFbyrow[i]);
    }
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JS, JSbyrow);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(JF, JFbyrow);

    // Calculate JS2=JS and JF2=JF again using multiplication by mobility-space
    // unit vectors.
    JS2.resize(nb, nu);
    JF2.resize(nb, nu);
    Vector zeroU(nu, 0.);
    for (int i = 0; i < nu; ++i) {
        zeroU[i] = 1;
        matter.multiplyByStationJacobian(state, allBodies, randS, zeroU, JS2(i));
        matter.multiplyByFrameJacobian(state, allBodies, randS, zeroU, JF2(i));
        zeroU[i] = 0;
    }
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JS2, JS, Slop);
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JF2, JF, Slop);

    // Calculate JS2t=~JS using multiplication by force-space unit vectors.
    Matrix_<Row3> JS2t(nu, nb);
    Vector_<Vec3> zeroF(nb, Vec3(0));
    // While we're at it, let's test non-contiguous vectors by filling in
    // this scalar version and using its non-contig rows as column temps.
    Matrix JS3mat(3 * nb, nu);
    for (int b = 0; b < nb; ++b) {
        for (int k = 0; k < 3; ++k) {
            zeroF[b][k] = 1;
            RowVectorView JS3matr = JS3mat[(3 * b) + k];
            matter.multiplyByStationJacobianTranspose(state, allBodies, randS, zeroF, ~JS3matr);
            zeroF[b][k] = 0;
            for (int u = 0; u < nu; ++u) {
                JS2t(u, b)[k] = JS3matr[u];
            }
        }
    }
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JS2, ~JS2t, Slop); // we'll check JS3mat below

    // Calculate JF2t=~JF using multiplication by force-space unit vectors.
    Matrix_<SpatialRow> JF2t(nu, nb);
    Vector_<SpatialVec> zeroSF(nb, SpatialVec(Vec3(0)));
    // While we're at it, let's test non-contiguous vectors by filling in
    // this scalar version and using its non-contig rows as column temps.
    Matrix JF3mat(6 * nb, nu);
    for (int b = 0; b < nb; ++b) {
        for (int k = 0; k < 6; ++k) {
            zeroSF[b][k / 3][k % 3] = 1;
            RowVectorView JF3matr = JF3mat[(6 * b) + k];
            matter.multiplyByFrameJacobianTranspose(state, allBodies, randS, zeroSF, ~JF3matr);
            zeroSF[b][k / 3][k % 3] = 0;
            for (int u = 0; u < nu; ++u) {
                JF2t(u, b)[k / 3][k % 3] = JF3matr[u];
            }
        }
    }
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JF2, ~JF2t, Slop); // we'll check JS3mat below

    // All three methods match. Now let's see if they are right by shifting
    // the System Jacobian to the new stations.

    for (int i = 0; i < nb; ++i) {
        const MobilizedBody& mobod = matter.getMobilizedBody(allBodies[i]);
        const Rotation& R_GB = mobod.getBodyRotation(state);
        const Vec3 S_G = R_GB * randS[i];
        for (int j = 0; j < nu; ++j) {
            const Vec3 w = J(i, j)[0];
            const Vec3 v = J(i, j)[1];
            const Vec3 vJ = v + w % S_G; // Shift
            const Vec3 vS = JS2(i, j);
            const SpatialVec vF = JF2(i, j);
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(vS, vJ);
            EXPECT_NEAR_DEFAULT_TOL_SIMTK(vF, SpatialVec(w, vJ));
        }
    }

    // Now create a scalar version of JS and make sure it matches the Vec3 one.
    Matrix JSmat;
    Matrix JSmat2;
    Matrix JFmat;
    Matrix JFmat2;

    matter.calcStationJacobian(state, allBodies, randS, JSmat);
    matter.calcFrameJacobian(state, allBodies, randS, JFmat);
    EXPECT_EQ(JSmat.nrow(), 3 * nb);
    EXPECT_EQ(JSmat.ncol(), nu);
    EXPECT_EQ(JFmat.nrow(), 6 * nb);
    EXPECT_EQ(JFmat.ncol(), nu);

    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JSmat, JS3mat, Slop); // same as above?
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JFmat, JF3mat, Slop); // same as above?

    // Unpack JS into JSmat2 and compare with JSmat.
    JSmat2.resize(3 * nb, nu);
    for (int row = 0; row < nb; ++row) {
        const int nxtr = 3 * row; // row index into scalar matrix
        for (int col = 0; col < nu; ++col) {
            for (int k = 0; k < 3; ++k) {
                JSmat2(nxtr + k, col) = JS(row, col)[k];
            }
        }
    }
    // These should be exactly the same.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JSmat2, JSmat, SignificantReal);

    // Unpack JF into JFmat2 and compare with JFmat.
    JFmat2.resize(6 * nb, nu);
    for (int row = 0; row < nb; ++row) {
        const int nxtr = 6 * row; // row index into scalar matrix
        for (int col = 0; col < nu; ++col) {
            for (int k = 0; k < 6; ++k) {
                JFmat2(nxtr + k, col) = JF(row, col)[k / 3][k % 3];
            }
        }
    }

    // These should be exactly the same.
    EXPECT_NEAR_CUSTOM_TOL_SIMTK(JFmat2, JFmat, SignificantReal);
}

/*
 * Position kinematics should be valid if:
 * - realize(Position) has been done
 * - or, realize(Instance) + realizePositionKinematics()
 * It should be invalidated when:
 * - any q changes
 * - Instance stage changes
 * It should *not* be invalidated when:
 * - time changes
 */
TEST(Simbody_SimbodyMatterSubsystem_PositionKinematics_CacheValidityAndInvalidationRules,
     InvalidatesOnQAndInstanceButNotTime) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(true, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    system.realizeModel(state);

    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    EXPECT_EQ(state.getSystemStage(), Stage::Model);

    Vec3 syscom;
    EXPECT_THROW // stage is too low
        (syscom = matter.calcSystemMassCenterLocationInGround(state), std::exception);

    system.realize(state, Stage::Instance);
    EXPECT_THROW // Instance alone is not enough
        (syscom = matter.calcSystemMassCenterLocationInGround(state), std::exception);

    matter.realizePositionKinematics(state);
    syscom = matter.calcSystemMassCenterLocationInGround(state); // OK
    matter.invalidatePositionKinematics(state);
    EXPECT_THROW // No good again
        (syscom = matter.calcSystemMassCenterLocationInGround(state), std::exception);

    system.realize(state, Stage::Position);
    syscom = matter.calcSystemMassCenterLocationInGround(state); // OK

    state.setTime(1.); // should not invalidate position kinematics
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);

    syscom = matter.calcSystemMassCenterLocationInGround(state); // OK

    state.updQ() = SimTK::Test::randVector(nq); // should invalidate position kinematics
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);
    EXPECT_THROW // Not up to date with q's
        (syscom = matter.calcSystemMassCenterLocationInGround(state), std::exception);

    matter.realizePositionKinematics(state);
    syscom = matter.calcSystemMassCenterLocationInGround(state); // OK

    state.invalidateAllCacheAtOrAbove(Stage::Instance);
    EXPECT_EQ(state.getSystemStage(), Stage::Model);
    EXPECT_THROW // stage too low again
        (syscom = matter.calcSystemMassCenterLocationInGround(state), std::exception);
}

/*
 * Velocity kinematics should be valid if:
 * - realize(Velocity) has been done
 * - or, position kinematics is valid and realizeVelocityKinematics() has
 *   been called (and that implies that Instance stage has been realized)
 * It should be invalidated when:
 * - any q or u changes
 * - position kinematics changes
 * - Instance stage changes
 * It should *not* be invalidated when:
 * - time changes
 */
TEST(Simbody_SimbodyMatterSubsystem_VelocityKinematics_CacheValidityAndDependencyOrdering,
     InvalidatesOnQ_U_PositionAndInstanceButNotTime) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(true, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int nb = matter.getNumBodies();

    system.realizeModel(state);

    // Randomize state.
    state.updQ() = SimTK::Test::randVector(nq);
    state.updU() = SimTK::Test::randVector(nu);
    EXPECT_EQ(state.getSystemStage(), Stage::Model);

    Vec3 syscomv;
    EXPECT_THROW // stage is too low
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    system.realize(state, Stage::Instance);
    EXPECT_THROW // Instance alone is not enough
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizePositionKinematics(state);
    EXPECT_THROW // Instance + pos kinematics is not enough
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizeVelocityKinematics(state);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK

    matter.invalidateVelocityKinematics(state);
    EXPECT_THROW // No good again
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizeVelocityKinematics(state);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK

    matter.invalidatePositionKinematics(state);
    EXPECT_THROW // No good again
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    system.realize(state, Stage::Velocity);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK

    state.setTime(1.); // should not invalidate velocity kinematics
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);

    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK
    state.updU() = SimTK::Test::randVector(nu);                   // should invalidate velocity kinematics
    EXPECT_EQ(state.getSystemStage(), Stage::Instance);
    EXPECT_THROW // Not up to date with u's
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizeVelocityKinematics(state);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK
    state.updQ() = SimTK::Test::randVector(nq);                   // invalidates position kinematics
    EXPECT_THROW                                                  // Pos kinematics no good
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizePositionKinematics(state);
    EXPECT_THROW // Still not enough; need to recalc vel kinematics
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    matter.realizeVelocityKinematics(state);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK

    state.invalidateAllCacheAtOrAbove(Stage::Instance);
    EXPECT_EQ(state.getSystemStage(), Stage::Model);
    EXPECT_THROW // stage too low again
        (syscomv = matter.calcSystemMassCenterVelocityInGround(state), std::exception);

    // Fails because realizePositionKinematics() or realize(Position) needed.
    EXPECT_THROW(matter.realizeVelocityKinematics(state), std::exception);
    system.realize(state, Stage::Position);
    matter.realizeVelocityKinematics(state);
    syscomv = matter.calcSystemMassCenterVelocityInGround(state); // OK
}
