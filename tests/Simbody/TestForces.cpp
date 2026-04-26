#include <array>
#include <gtest/gtest.h>

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

const int NUM_BODIES = 10;
const Real BOND_LENGTH = 0.5;

void verifyForces(const Force& force,
                  const State& state,
                  Vector_<SpatialVec> bodyForces,
                  Vector_<Vec3> particleForces,
                  Vector mobilityForces) {
    Vector_<SpatialVec> actualBodyForces(bodyForces.size());
    Vector_<Vec3> actualParticleForces(particleForces.size());
    Vector actualMobilityForces(mobilityForces.size());
    force.calcForceContribution(state, actualBodyForces, actualParticleForces, actualMobilityForces);

    for (int i = 0; i < bodyForces.size(); ++i) {
        EXPECT_NEAR(0.0, (bodyForces[i] - actualBodyForces[i]).norm(), 1e-10);
    }
    for (int i = 0; i < particleForces.size(); ++i) {
        EXPECT_NEAR(0.0, (particleForces[i] - actualParticleForces[i]).norm(), 1e-10);
    }
    for (int i = 0; i < mobilityForces.size(); ++i) {
        EXPECT_NEAR(0.0, std::abs(mobilityForces[i] - actualMobilityForces[i]), 1e-10);
    }
}

class MyForceImpl : public Force::Custom::Implementation {
    public:
    mutable std::array<bool, Stage::Report + 1> hasRealized;

    MyForceImpl() {
        for (int i = 0; i < Stage::NValid; i++) {
            hasRealized[i] = false;
        }
    }

    void calcForce(const State& state,
                   Vector_<SpatialVec>& bodyForces,
                   Vector_<Vec3>& particleForces,
                   Vector& mobilityForces) const override {
        for (int i = 0; i < mobilityForces.size(); ++i) {
            mobilityForces[i] += i;
        }
    }
    auto calcPotentialEnergy(const State& state) const -> Real override {
        return 0.0;
    }

    void realizeTopology(State& state) const override {
        hasRealized[Stage::Topology] = true;
    }

    void realizeModel(State& state) const override {
        hasRealized[Stage::Model] = true;
    }

    void realizeInstance(const State& state) const override {
        hasRealized[Stage::Instance] = true;
    }

    void realizeTime(const State& state) const override {
        hasRealized[Stage::Time] = true;
    }

    void realizePosition(const State& state) const override {
        hasRealized[Stage::Position] = true;
    }

    void realizeVelocity(const State& state) const override {
        hasRealized[Stage::Velocity] = true;
    }

    void realizeDynamics(const State& state) const override {
        hasRealized[Stage::Dynamics] = true;
    }

    void realizeAcceleration(const State& state) const override {
        hasRealized[Stage::Acceleration] = true;
    }

    void realizeReport(const State& state) const override {
        hasRealized[Stage::Report] = true;
    }
};

/**
 * Test all of the standard Force subclasses, and make sure they generate correct forces.
 */
TEST(Simbody_Forces_StandardForceSubclasses, CorrectBodyParticleAndMobilityForceEvaluation) {
    // Create a system consisting of a chain of bodies.

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNumBodies() - 1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }

    // Add one of each type of force.
    MobilizedBody& body1 = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& body9 = matter.updMobilizedBody(MobilizedBodyIndex(9));
    Force::ConstantForce constantForce(forces, body1, Vec3(0), Vec3(1, 2, 3));
    Force::ConstantTorque constantTorque(forces, body1, Vec3(1, 2, 3));
    Force::GlobalDamper globalDamper(forces, matter, 2.0);
    Force::MobilityConstantForce mobilityConstantForce(forces, body1, 1, 2.0);
    Force::MobilityLinearDamper mobilityLinearDamper(forces, body1, 1, 2.0);
    Force::MobilityLinearSpring mobilityLinearSpring(forces, body1, 1, 2.0, 1.0);
    Force::TwoPointConstantForce twoPointConstantForce(forces, body1, Vec3(0), body9, Vec3(0), 2.0);
    Force::TwoPointLinearDamper twoPointLinearDamper(forces, body1, Vec3(0), body9, Vec3(0), 2.0);
    Force::TwoPointLinearSpring twoPointLinearSpring(forces, body1, Vec3(0), body9, Vec3(0), 2.0, 0.5);
    Force::UniformGravity uniformGravity(forces, matter, Vec3(0, -2.0, 0));
    Force::Custom custom(forces, new MyForceImpl());

    // Create a random state for it.

    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Uniform random;
    for (int i = 0; i < state.getNY(); ++i) {
        state.updY()[i] = random.getValue();
    }
    system.realize(state, Stage::Velocity);
    Vec3 pos1 = body1.getBodyOriginLocation(state);
    Vec3 pos9 = body9.getBodyOriginLocation(state);
    Vec3 delta19 = pos9 - pos1;

    // Calculate each force component and see if it is correct.
    Vector_<SpatialVec> bodyForces(matter.getNumBodies());
    Vector_<Vec3> particleForces(0);
    Vector mobilityForces(state.getNU());
    Real pe = 0;

    // Check ConstantForce
    {
        SCOPED_TRACE("ConstantForce");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        bodyForces[1][1] = Vec3(1, 2, 3);
        verifyForces(constantForce, state, bodyForces, particleForces, mobilityForces);
    }

    // Check ConstantTorque
    {
        SCOPED_TRACE("ConstantTorque");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        bodyForces[1][0] = Vec3(1, 2, 3);
        verifyForces(constantTorque, state, bodyForces, particleForces, mobilityForces);
    }

    // Check GlobalDamper
    {
        SCOPED_TRACE("GlobalDamper");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = -2.0 * state.getU();
        verifyForces(globalDamper, state, bodyForces, particleForces, mobilityForces);
    }

    // Check MobilityConstantForce
    {
        SCOPED_TRACE("MobilityConstantForce");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        body1.updOneFromUPartition(state, 1, mobilityForces) = 2.0;
        verifyForces(mobilityConstantForce, state, bodyForces, particleForces, mobilityForces);
    }

    // Check MobilityLinearDamper
    {
        SCOPED_TRACE("MobilityLinearDamper");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        body1.updOneFromUPartition(state, 1, mobilityForces) = -2.0 * body1.getOneU(state, 1);
        verifyForces(mobilityLinearDamper, state, bodyForces, particleForces, mobilityForces);
    }

    // Check MobilityLinearSpring
    {
        SCOPED_TRACE("MobilityLinearSpring");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        body1.updOneFromUPartition(state, 1, mobilityForces) = -2.0 * (body1.getOneQ(state, 1) - 1.0);
        verifyForces(mobilityLinearSpring, state, bodyForces, particleForces, mobilityForces);
    }

    // Check TwoPointConstantForce
    {
        SCOPED_TRACE("TwoPointConstantForce");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        bodyForces[1][1] = -2.0 * delta19.normalize();
        bodyForces[9][1] = 2.0 * delta19.normalize();
        verifyForces(twoPointConstantForce, state, bodyForces, particleForces, mobilityForces);
    }

    // Check TwoPointLinearDamper
    {
        SCOPED_TRACE("TwoPointLinearDamper");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        Vec3 v19 = body9.getBodyOriginVelocity(state) - body1.getBodyOriginVelocity(state);
        Real twoPointLinearDamperForce = 2.0 * dot(v19, delta19.normalize());
        bodyForces[1][1] = twoPointLinearDamperForce * delta19.normalize();
        bodyForces[9][1] = -twoPointLinearDamperForce * delta19.normalize();
        verifyForces(twoPointLinearDamper, state, bodyForces, particleForces, mobilityForces);
    }

    // Check TwoPointLinearSpring
    {
        SCOPED_TRACE("TwoPointLinearSpring");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        Real twoPointLinearSpringForce = 2.0 * (delta19.norm() - 0.5);
        bodyForces[1][1] = twoPointLinearSpringForce * delta19.normalize();
        bodyForces[9][1] = -twoPointLinearSpringForce * delta19.normalize();
        verifyForces(twoPointLinearSpring, state, bodyForces, particleForces, mobilityForces);
    }

    // Check UniformGravity
    {
        SCOPED_TRACE("UniformGravity");
        bodyForces = SpatialVec(Vec3(0), Vec3(0, -2.0, 0));
        bodyForces[0] = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        mobilityForces = 0;
        verifyForces(uniformGravity, state, bodyForces, particleForces, mobilityForces);
    }

    // Check Custom
    {
        SCOPED_TRACE("Custom");
        bodyForces = SpatialVec(Vec3(0), Vec3(0));
        particleForces = Vec3(0);
        for (int i = 0; i < mobilityForces.size(); ++i) {
            mobilityForces[i] = i;
        }
        verifyForces(custom, state, bodyForces, particleForces, mobilityForces);
    }
}

/**
 * Test the standard conservative forces to make sure they really conserve energy.
 */
TEST(Simbody_Forces_ConservativeForces_ConserveTotalEnergyOverTime, MaintainsEnergyWithinTolerance) {
    // Create a system consisting of a chain of bodies.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    for (int i = 0; i < NUM_BODIES; ++i) {
        MobilizedBody& parent = matter.updMobilizedBody(MobilizedBodyIndex(matter.getNumBodies() - 1));
        MobilizedBody::Gimbal b(parent, Transform(Vec3(0)), body, Transform(Vec3(BOND_LENGTH, 0, 0)));
    }

    // Add one of each type of conservative force.
    MobilizedBody& body1 = matter.updMobilizedBody(MobilizedBodyIndex(1));
    MobilizedBody& body9 = matter.updMobilizedBody(MobilizedBodyIndex(9));
    Force::MobilityLinearSpring mobilityLinearSpring(forces, body1, 1, 0.1, 1.0);
    Force::TwoPointLinearSpring twoPointLinearSpring(forces, body1, Vec3(0), body9, Vec3(0), 1.0, 4.0);
    Force::UniformGravity uniformGravity(forces, matter, Vec3(0, -1.0, 0));

    // Create a random initial state for it.
    system.realizeTopology();
    State state = system.getDefaultState();
    Random::Uniform random;
    for (int i = 0; i < state.getNY(); ++i) {
        state.updY()[i] = random.getValue();
    }

    // Simulate it for a while and see if the energy changes.
    system.realize(state, Stage::Dynamics);
    Real initialEnergy = system.calcEnergy(state);

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-4);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(10.0));

    system.realize(state, Stage::Dynamics);
    Real finalEnergy = system.calcEnergy(timeStepper.getState());
    EXPECT_LT(std::abs((initialEnergy / finalEnergy) - 1.0), 0.005);
}

/**
 * Make sure that all the "realize" methods on a custom force actually get called.
 */
TEST(Simbody_Forces_CustomForce_RealizeMethods, AllStagesAreCalled) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    auto* impl = new MyForceImpl();
    Force::Custom custom(forces, impl);
    State state = system.realizeTopology();
    for (Stage j = Stage::Model; j <= Stage::Report; j++) {
        system.realize(state, j);
        for (Stage i = Stage::Topology; i <= Stage::Report; i++) {
            EXPECT_EQ(impl->hasRealized[i], (i <= j));
        }
    }
}

/**
 * Test enabling and disabling forces.
 */
TEST(Simbody_Forces_DisablingForces, CorrectlyEnablesAndDisablesForces) {
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Free body1(matter.updGround(), Vec3(0), body, Vec3(0));
    MobilizedBody::Free body2(matter.updGround(), Vec3(0), body, Vec3(0));
    Force::TwoPointLinearSpring spring(forces, body1, Vec3(0), body2, Vec3(0), 2.0, 0.5);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -2.0, 0));

    // Create an initial state.
    State state = system.realizeTopology();
    body1.setQToFitTranslation(state, Vec3(0, 1, 0));
    body2.setQToFitTranslation(state, Vec3(1, 1, 0));

    // These are the contribution of each force to the energy and to the force on body1.
    const Real springEnergy = 0.5 * 2.0 * 0.5 * 0.5;
    const SpatialVec springForce(Vec3(0), Vec3(2.0 * 0.5, 0, 0));
    const Real gravityEnergy = 2 * 2.0;
    const SpatialVec gravityForce(Vec3(0), Vec3(0, -2.0, 0));

    // Verify the force and energy for each combination of the forces being enabled or disabled.
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(springEnergy + gravityEnergy, system.calcEnergy(state));
    EXPECT_LT((springForce + gravityForce - system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm(),
              1e-10);
    EXPECT_FALSE(forces.isForceDisabled(state, gravity.getForceIndex()));
    EXPECT_FALSE(forces.isForceDisabled(state, spring.getForceIndex()));

    forces.setForceIsDisabled(state, spring.getForceIndex(), true);
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(gravityEnergy, system.calcEnergy(state));
    EXPECT_LT((gravityForce - system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm(), 1e-10);
    EXPECT_FALSE(forces.isForceDisabled(state, gravity.getForceIndex()));
    EXPECT_TRUE(forces.isForceDisabled(state, spring.getForceIndex()));

    forces.setForceIsDisabled(state, gravity.getForceIndex(), true);
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(0, system.calcEnergy(state));
    EXPECT_LT((system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm(), 1e-10);
    EXPECT_TRUE(forces.isForceDisabled(state, gravity.getForceIndex()));
    EXPECT_TRUE(forces.isForceDisabled(state, spring.getForceIndex()));

    forces.setForceIsDisabled(state, spring.getForceIndex(), false);
    system.realize(state, Stage::Dynamics);
    EXPECT_EQ(springEnergy, system.calcEnergy(state));
    EXPECT_LT((springForce - system.getRigidBodyForces(state, Stage::Dynamics)[1]).norm(), 1e-10);
    EXPECT_TRUE(forces.isForceDisabled(state, gravity.getForceIndex()));
    EXPECT_FALSE(forces.isForceDisabled(state, spring.getForceIndex()));
}
