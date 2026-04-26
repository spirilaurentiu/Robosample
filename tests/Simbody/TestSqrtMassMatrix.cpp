#include <gtest/gtest.h>

#include "SimTKcommon/Testing.h"

#include "SimTKsimbody.h"

using namespace SimTK;

// This will apply a constant set of mobility forces that can be set externally.
class MyForceImpl : public Force::Custom::Implementation {
    public:
    MyForceImpl() = default;
    void calcForce(const State& state,
                   Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces,
                   Vector& mobilityForces) const override {
        SimTK_TEST(f.size() == 0 || f.size() == mobilityForces.size());
        SimTK_TEST(F.size() == 0 || F.size() == bodyForces.size());
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

TEST(Simbody_SqrtMassMatrix_SqrtMinv_ReconstructsMinv, ConsistencyWithInverse) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    system.realizeModel(state);

    const int nu = state.getNU();
    const Real Slop = nu * SignificantReal;

    state.updQ() = SimTK::Test::randVector(state.getNQ());
    state.updU() = SimTK::Test::randVector(nu);

    system.realize(state, Stage::Position);

    Matrix SqrtMInv(nu, nu, 0.0);

    Vector v(nu), e(nu);

    for (int i = 0; i < nu; ++i) {
        e.setToZero();
        e[i] = 1;

        matter.multiplyBySqrtMInv(state, e, v);
        SqrtMInv(i) = v;
    }

    Matrix MInv;
    matter.calcMInv(state, MInv);

    EXPECT_TRUE(SimTK::Test::numericallyEqual(SqrtMInv * ~SqrtMInv, MInv, nu));
    EXPECT_TRUE(SimTK::Test::numericallyEqual(SqrtMInv * ~SqrtMInv, MInv, 1, Slop));
}

TEST(Simbody_SqrtMassMatrix_MultiplyByM, ForwardMapConsistency) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    system.realizeModel(state);

    const int nu = state.getNU();
    const Real Slop = nu * SignificantReal;

    state.updQ() = SimTK::Test::randVector(state.getNQ());
    state.updU() = SimTK::Test::randVector(nu);

    system.realize(state, Stage::Velocity);

    Vector v = SimTK::Test::randVector(nu);
    Vector Mv(nu), Mv_direct(nu);

    matter.multiplyByM(state, v, Mv);

    Matrix M;
    matter.calcM(state, M);
    Mv_direct = M * v;

    EXPECT_TRUE(SimTK::Test::numericallyEqual(Mv, Mv_direct, 1, Slop));
}

TEST(Simbody_SqrtMassMatrix_MultiplyByMInv, InverseConsistency) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    system.realizeModel(state);

    const int nu = state.getNU();
    const Real Slop = nu * SignificantReal;

    state.updQ() = SimTK::Test::randVector(state.getNQ());
    state.updU() = SimTK::Test::randVector(nu);

    system.realize(state, Stage::Velocity);

    Vector v = SimTK::Test::randVector(nu);
    Vector Mv(nu), back(nu);

    matter.multiplyByM(state, v, Mv);
    matter.multiplyByMInv(state, Mv, back);

    EXPECT_TRUE(SimTK::Test::numericallyEqual(v, back, 1, Slop));
}

TEST(Simbody_SqrtMassMatrix_MaxwellBoltzmann_Invariance, KineticEnergyConsistency) {
    MultibodySystem system;
    MyForceImpl* frcp;
    makeSystem(false, system, frcp);
    const SimbodyMatterSubsystem& matter = system.getMatterSubsystem();

    State state = system.realizeTopology();
    system.realizeModel(state);

    const int nu = state.getNU();
    const Real Slop = nu * SignificantReal;

    state.updQ() = SimTK::Test::randVector(state.getNQ());
    state.updU() = SimTK::Test::randVector(nu);

    system.realize(state, Stage::Velocity);

    Vector xi = SimTK::Test::randVector(nu);
    Vector v(nu);

    matter.multiplyBySqrtMInv(state, xi, v);

    Vector Mv(nu);
    matter.multiplyByM(state, v, Mv);

    EXPECT_EQ(v.size(), Mv.size());
    Real KE_physical = 0;
    for (std::size_t i = 0; i < v.size(); ++i) {
        KE_physical += 0.5 * v[i] * Mv[i];
    }

    Real KE_white = 0;
    for (double i : xi) {
        KE_white += 0.5 * i * i;
    }

    EXPECT_NEAR(KE_physical, KE_white, Slop);
}
