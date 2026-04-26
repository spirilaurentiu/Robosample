/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

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


/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Translation.
 */
class CustomTranslation : public MobilizedBody::Custom::Implementation {
    public:
    CustomTranslation(SimbodyMatterSubsystem& matter)
        : Implementation(matter, 3, 3, 0) {
    }
    Implementation* clone() const override {
        return new CustomTranslation(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const override {
        EXPECT_EQ(nq, 3);
        return Transform(Vec3(q[0], q[1], q[2]));
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 3);
        return SpatialVec(Vec3(0), Vec3(u[0], u[1], u[2]));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(f) = F[1];
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 3);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(f) = Vec3(0);
    }
    void setQToFitTransform(const State&, const Transform& X_FM, int nq, Real* q) const override {
        EXPECT_EQ(nq, 3);
        Vec3::updAs(q) = X_FM.p();
    }
    void setUToFitVelocity(const State&, const SpatialVec& V_FM, int nu, Real* u) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(u) = V_FM[1];
    }
};

/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Ball.
 */

class CustomBall : public MobilizedBody::Custom::Implementation {
    public:
    CustomBall(SimbodyMatterSubsystem& matter)
        : Implementation(matter, 3, 4, 4) {
    }
    Implementation* clone() const override {
        return new CustomBall(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const override {
        Transform t(Vec3(0));
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nq, 3);
            t.updR().setRotationToBodyFixedXYZ(Vec3::getAs(q));
        } else {
            EXPECT_EQ(nq, 4);
            t.updR().setRotationFromQuaternion(Quaternion(Vec4::getAs(q)));
        }
        return t;
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 3);
        return SpatialVec(Vec3(u[0], u[1], u[2]), Vec3(0));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(f) = F[0];
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 3);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(f) = Vec3(0);
    }
    void multiplyByN(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 3);
            EXPECT_EQ(nOut, 3);

            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 N = Rotation::calcNForBodyXYZInBodyFrame(Vec3::getAs(&q[0])) * ~R_FM;
            if (transposeMatrix) {
                Row3::updAs(out) = Row3::getAs(in) * N;
            } else {
                Vec3::updAs(out) = N * Vec3::getAs(in);
            }
        } else {
            if (transposeMatrix) {
                EXPECT_EQ(nIn, 4);
                EXPECT_EQ(nOut, 3);
            } else {
                EXPECT_EQ(nIn, 3);
                EXPECT_EQ(nOut, 4);
            }
            const Mat43 N = Rotation::calcUnnormalizedNForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) {
                Row3::updAs(out) = Row4::getAs(in) * N;
            } else {
                Vec4::updAs(out) = N * Vec3::getAs(in);
            }
        }
    }
    void multiplyByNInv(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 3);
            EXPECT_EQ(nOut, 3);

            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 NInv = R_FM * Rotation::calcNInvForBodyXYZInBodyFrame(Vec3::getAs(&q[0]));
            if (transposeMatrix) {
                Row3::updAs(out) = Row3::getAs(in) * NInv;
            } else {
                Vec3::updAs(out) = NInv * Vec3::getAs(in);
            }
        } else {
            if (transposeMatrix) {
                EXPECT_EQ(nIn, 3);
                EXPECT_EQ(nOut, 4);
            } else {
                EXPECT_EQ(nIn, 4);
                EXPECT_EQ(nOut, 3);
            }
            const Mat34 NInv = Rotation::calcUnnormalizedNInvForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) {
                Row4::updAs(out) = Row3::getAs(in) * NInv;
            } else {
                Vec3::updAs(out) = NInv * Vec4::getAs(in);
            }
        }
    }
    void multiplyByNDot(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 3);
            EXPECT_EQ(nOut, 3);
            const Rotation& R_FM = getMobilizerTransform(s).R();
            Vec3::updAs(out) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot(Vec3::getAs(&q[0]),
                                                                                    ~R_FM * Vec3::getAs(in),
                                                                                    Vec3(0));
        } else if (transposeMatrix) {
            EXPECT_EQ(nIn, 4);
            EXPECT_EQ(nOut, 3);
            EXPECT_FALSE(true); // I didn't bother to implement this case, since it currently is never used.
        } else {
            EXPECT_EQ(nIn, 3);
            EXPECT_EQ(nOut, 4);
            Vec4::updAs(out) =
                Rotation::convertAngVelDotToQuaternionDotDot(Vec4::getAs(&q[0]), Vec3::getAs(in), Vec3(0));
        }
    }
    void setQToFitTransform(const State& s, const Transform& X_FM, int nq, Real* q) const override {
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nq, 3);
            Vec3::updAs(q) = X_FM.R().convertRotationToBodyFixedXYZ();
        } else {
            EXPECT_EQ(nq, 4);
            Vec4::updAs(q) = X_FM.R().convertRotationToQuaternion().asVec4();
        }
    }
    void setUToFitVelocity(const State& s, const SpatialVec& V_FM, int nu, Real* u) const override {
        EXPECT_EQ(nu, 3);
        Vec3::updAs(u) = V_FM[0];
    }
};

/**
 * This is a custom MobilizedBody that is identical to MobilizedBody::Free.
 */
class CustomFree : public MobilizedBody::Custom::Implementation {
    public:
    CustomFree(SimbodyMatterSubsystem& matter)
        : Implementation(matter, 6, 7, 4) {
    }
    Implementation* clone() const override {
        return new CustomFree(*this);
    }
    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const override {
        Transform t(Vec3::getAs(&q[nq - 3]));
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nq, 6);
            t.updR().setRotationToBodyFixedXYZ(Vec3::getAs(q));
        } else {
            EXPECT_EQ(nq, 7);
            t.updR().setRotationFromQuaternion(Quaternion(Vec4::getAs(q)));
        }
        return t;
    }
    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 6);
        return SpatialVec(Vec3(u[0], u[1], u[2]), Vec3(u[3], u[4], u[5]));
    }
    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(f)) = F;
    }
    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const override {
        EXPECT_EQ(nu, 6);
        return SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        EXPECT_EQ(nu, 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(f)) = SpatialVec(Vec3(0), Vec3(0));
    }
    void multiplyByN(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 6);
            EXPECT_EQ(nOut, 6);

            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 N = Rotation::calcNForBodyXYZInBodyFrame(Vec3::getAs(&q[0])) * ~R_FM;
            if (transposeMatrix) {
                Row3::updAs(out) = Row3::getAs(in) * N;
            } else {
                Vec3::updAs(out) = N * Vec3::getAs(in);
            }
        } else {
            if (transposeMatrix) {
                EXPECT_EQ(nIn, 7);
                EXPECT_EQ(nOut, 6);
            } else {
                EXPECT_EQ(nIn, 6);
                EXPECT_EQ(nOut, 7);
            }
            const Mat43 N = Rotation::calcUnnormalizedNForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) {
                Row3::updAs(out) = Row4::getAs(in) * N;
            } else {
                Vec4::updAs(out) = N * Vec3::getAs(in);
            }
        }
        Vec3::updAs(&out[nOut - 3]) = Vec3::getAs(&in[nIn - 3]);
    }
    void multiplyByNInv(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 6);
            EXPECT_EQ(nOut, 6);

            Rotation R_FM;
            R_FM.setRotationToBodyFixedXYZ(Vec3::getAs(&q[0]));
            const Mat33 NInv = R_FM * Rotation::calcNInvForBodyXYZInBodyFrame(Vec3::getAs(&q[0]));
            if (transposeMatrix) {
                Row3::updAs(out) = Row3::getAs(in) * NInv;
            } else {
                Vec3::updAs(out) = NInv * Vec3::getAs(in);
            }
        } else {
            if (transposeMatrix) {
                EXPECT_EQ(nIn, 6);
                EXPECT_EQ(nOut, 7);
            } else {
                EXPECT_EQ(nIn, 7);
                EXPECT_EQ(nOut, 6);
            }
            const Mat34 NInv = Rotation::calcUnnormalizedNInvForQuaternion(Vec4::getAs(&q[0]));
            if (transposeMatrix) {
                Row4::updAs(out) = Row3::getAs(in) * NInv;
            } else {
                Vec3::updAs(out) = NInv * Vec4::getAs(in);
            }
        }
        Vec3::updAs(&out[nOut - 3]) = Vec3::getAs(&in[nIn - 3]);
    }
    void multiplyByNDot(const State& s, bool transposeMatrix, int nIn, const Real* in, int nOut, Real* out)
        const override {
        const Vector q = getQ(s);
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nIn, 6);
            EXPECT_EQ(nOut, 6);
            const Rotation& R_FM = getMobilizerTransform(s).R();
            Vec3::updAs(out) = Rotation::convertAngVelDotInBodyFrameToBodyXYZDotDot(Vec3::getAs(&q[0]),
                                                                                    ~R_FM * Vec3::getAs(in),
                                                                                    Vec3(0));
        } else if (transposeMatrix) {
            EXPECT_EQ(nIn, 7);
            EXPECT_EQ(nOut, 6);
            EXPECT_FALSE(true); // I didn't bother to implement this case, since it currently is never used.
        } else {
            EXPECT_EQ(nIn, 6);
            EXPECT_EQ(nOut, 7);
            Vec4::updAs(out) =
                Rotation::convertAngVelDotToQuaternionDotDot(Vec4::getAs(&q[0]), Vec3::getAs(in), Vec3(0));
        }
        Vec3::updAs(&out[nOut - 3]) = Vec3(0);
    }
    void setQToFitTransform(const State& s, const Transform& X_FM, int nq, Real* q) const override {
        if (getUseEulerAngles(s)) {
            EXPECT_EQ(nq, 6);
            Vec3::updAs(q) = X_FM.R().convertRotationToBodyFixedXYZ();
        } else {
            EXPECT_EQ(nq, 7);
            Vec4::updAs(q) = X_FM.R().convertRotationToQuaternion().asVec4();
        }
        Vec3::updAs(&q[nq - 3]) = X_FM.p();
    }
    void setUToFitVelocity(const State& s, const SpatialVec& V_FM, int nu, Real* u) const override {
        EXPECT_EQ(nu, 6);
        SpatialVec::updAs(reinterpret_cast<Vec3*>(u)) = V_FM;
    }
};

void compareMobilizedBodies(const MobilizedBody& b1,
                            const MobilizedBody& b2,
                            bool eulerAngles,
                            int expectedQ,
                            int expectedU) {
    const SimbodyMatterSubsystem& matter = b1.getMatterSubsystem();
    const System& system = matter.getSystem();

    // Set whether to use Euler angles.
    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, eulerAngles);
    system.realizeModel(state);

    // Make sure the number of state variables is correct.
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getNumQ(state), expectedQ);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getNumU(state), expectedU);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getNumQ(state), expectedQ);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getNumU(state), expectedU);

    // Set all the state variables to random values.
    Random::Gaussian random;
    int nq = state.getNQ() / 2;
    for (int i = 0; i < nq; ++i) {
        state.updQ()[i] = state.updQ()[i + nq] = random.getValue();
    }
    int nu = state.getNU() / 2;
    for (int i = 0; i < nu; ++i) {
        state.updU()[i] = state.updU()[i + nu] = random.getValue();
    }
    system.realize(state, Stage::Acceleration);

    // Compare state variables and their derivatives.
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneQ(state, i), b2.getOneQ(state, i));
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneQDot(state, i), b2.getOneQDot(state, i));
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneQDotDot(state, i), b2.getOneQDotDot(state, i));
    }
    for (int i = 0; i < b1.getNumU(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneU(state, i), b2.getOneU(state, i));
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneUDot(state, i), b2.getOneUDot(state, i));
    }

    // Compare lots of properties of the two bodies.
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyTransform(state), b2.getBodyTransform(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyVelocity(state), b2.getBodyVelocity(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyAcceleration(state), b2.getBodyAcceleration(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyOriginVelocity(state), b2.getBodyOriginVelocity(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyOriginAcceleration(state), b2.getBodyOriginAcceleration(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getMobilizerTransform(state), b2.getMobilizerTransform(state));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getMobilizerVelocity(state), b2.getMobilizerVelocity(state));

    // Test methods that multiply by various matrices.
    Vector tempq(state.getNQ());
    Vector tempu(state.getNU());
    matter.multiplyByN(state, false, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneFromQPartition(state, i, tempq),
                                      b2.getOneFromQPartition(state, i, tempq));
    }
    matter.multiplyByN(state, true, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneFromUPartition(state, i, tempu),
                                      b2.getOneFromUPartition(state, i, tempu));
    }
    matter.multiplyByNInv(state, false, state.getQ(), tempu);
    for (int i = 0; i < b1.getNumU(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneFromUPartition(state, i, tempu),
                                      b2.getOneFromUPartition(state, i, tempu));
    }
    matter.multiplyByNInv(state, true, state.getU(), tempq);
    for (int i = 0; i < b1.getNumQ(state); ++i) {
        EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getOneFromQPartition(state, i, tempq),
                                      b2.getOneFromQPartition(state, i, tempq));
    }

    // Have them calculate q and u, and see if they agree.
    Transform t(Rotation(random.getValue(), Vec3(random.getValue(), random.getValue(), random.getValue())),
                Vec3(random.getValue(), random.getValue(), random.getValue()));
    b1.setQToFitTransform(state, t);
    b2.setQToFitTransform(state, t);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getQAsVector(state), b2.getQAsVector(state));

    SpatialVec v(Vec3(random.getValue(), random.getValue(), random.getValue()),
                 Vec3(random.getValue(), random.getValue(), random.getValue()));
    b1.setUToFitVelocity(state, v);
    b2.setUToFitVelocity(state, v);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getUAsVector(state), b2.getUAsVector(state));

    // Simulate the system, and see if the two bodies remain identical.
    VerletIntegrator integ(system);
    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getQAsVector(integ.getState()), b2.getQAsVector(integ.getState()));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getUAsVector(integ.getState()), b2.getUAsVector(integ.getState()));
}

TEST(Simbody_CustomMobilizedBodies, CustomTranslation) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Translation t1(matter.Ground(), body);
    MobilizedBody::Translation t2(t1, body);
    MobilizedBody::Custom c1(matter.Ground(), new CustomTranslation(matter), body);
    MobilizedBody::Custom c2(c1, new CustomTranslation(matter), body);
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 3, 3);
    compareMobilizedBodies(t2, c2, true, 3, 3);
}

TEST(Simbody_CustomMobilizedBodies, CustomBall) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Ball t1(matter.Ground(), body);
    MobilizedBody::Ball t2(t1, Vec3(1, 0, 0), body, Vec3(0, 2, 0));
    MobilizedBody::Custom c1(matter.Ground(), new CustomBall(matter), body);
    MobilizedBody::Custom c2(c1, new CustomBall(matter), Vec3(1, 0, 0), body, Vec3(0, 2, 0));
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 4, 3);
    compareMobilizedBodies(t2, c2, true, 3, 3);
}

TEST(Simbody_CustomMobilizedBodies, CustomFree) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free t1(matter.Ground(), body);
    MobilizedBody::Free t2(t1, body);
    MobilizedBody::Custom c1(matter.Ground(), new CustomFree(matter), body);
    MobilizedBody::Custom c2(c1, new CustomFree(matter), body);
    system.realizeTopology();
    compareMobilizedBodies(t2, c2, false, 7, 6);
    compareMobilizedBodies(t2, c2, true, 6, 6);
}
