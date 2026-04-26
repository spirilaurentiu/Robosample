#include <gtest/gtest.h>
#include <vector>

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
        state.updU()[i] = state.updU()[i + nu] = (eulerAngles ? 0.0 : random.getValue());
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
    if (!eulerAngles) { // The optimizer does not work reliably for Euler angles, since it can hit a
                        // singularity
        Transform t = b1.getBodyTransform(state);
        b1.setQFromVector(state, Vector(b1.getNumQ(state), 0.0));
        b2.setQFromVector(state, Vector(b2.getNumQ(state), 0.0));
        b1.setQToFitTransform(state, t);
        b2.setQToFitTransform(state, t);
        system.realize(state, Stage::Velocity);
        EXPECT_NEAR_CUSTOM_TOL_SIMTK(b1.getBodyOriginLocation(state), b2.getBodyOriginLocation(state), 1e-2);
        EXPECT_NEAR_CUSTOM_TOL_SIMTK(
            (~b1.getBodyRotation(state) * b2.getBodyRotation(state)).convertRotationToAngleAxis()[0],
            0.0,
            1e-2);
        SpatialVec v = b1.getBodyVelocity(state);
        b1.setUFromVector(state, Vector(b1.getNumU(state), 0.0));
        b2.setUFromVector(state, Vector(b2.getNumU(state), 0.0));
        b1.setUToFitVelocity(state, v);
        b2.setUToFitVelocity(state, v);
        EXPECT_NEAR_CUSTOM_TOL_SIMTK(b1.getUAsVector(state), b2.getUAsVector(state), 1e-2);
    }

    // Simulate the system, and see if the two bodies remain identical.
    b2.setQFromVector(state, b1.getQAsVector(state));
    b2.setUFromVector(state, b1.getUAsVector(state));

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getQAsVector(integ.getState()), b2.getQAsVector(integ.getState()));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getQDotAsVector(integ.getState()), b2.getQDotAsVector(integ.getState()));
}

class ConstantFunction : public Function {
    // Implements a simple constant function, y = C
    private:
    Real C;

    public:
    // Default constructor
    ConstantFunction() {
        C = 0.0;
    }

    // Convenience constructor to specify constant value
    ConstantFunction(Real constant) {
        C = constant;
    }

    Real calcValue(const Vector& x) const override {
        return C;
    }

    // This is the pure virtual signature.
    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const override {
        return 0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const override {
        // constant has no arguments
        return 0;
    }

    int getMaxDerivativeOrder() const override {
        return 10;
    }
};

class LinearFunction : public Function {
    // Implements a simple linear functional relationship, y = m*x + b
    private:
    Real m;
    Real b;

    public:
    // Default constructor
    LinearFunction() {
        m = 1.0;
        b = 0.0;
    }

    // Convenience constructor to specify the slope and Y-intercept of the linear r
    LinearFunction(Real slope, Real intercept) {
        m = slope;
        b = intercept;
    }

    Real calcValue(const Vector& x) const override {
        return m * x[0] + b;
    }

    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const override {
        if (derivComponents.size() == 1) {
            return m;
        }
        return 0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const override {
        return 1;
    }

    int getMaxDerivativeOrder() const override {
        return 10;
    }
};

class NonlinearFunction : public Function {
    public:
    NonlinearFunction() {
    }
    Real calcValue(const Vector& x) const override {
        return x[0] + x[1] * x[1];
    }
    Real calcDerivative(const Array_<int>& derivComponents, const Vector& x) const override {
        switch (derivComponents.size()) {
            case 1:
                return (derivComponents[0] == 0 ? 1.0 : x[1]);
            case 2:
                return (derivComponents[0] == 1 && derivComponents[1] == 1 ? 1.0 : 0.0);
        }
        return 0.0;
    }
    // This is just a local method providing std::vector compatibility without copying.
    Real calcDerivative(const std::vector<int>& derivComponents, const Vector& x) const {
        return calcDerivative(ArrayViewConst_<int>(derivComponents), x);
    }

    int getArgumentSize() const override {
        return 2;
    }
    int getMaxDerivativeOrder() const override {
        return std::numeric_limits<int>::max();
    }
};

int defineMobilizerFunctions(const std::vector<bool>& isdof,
                             std::vector<std::vector<int>>& coordIndices,
                             std::vector<const Function*>& functions1,
                             std::vector<const Function*>& functions2) {
    int nm = 0;
    for (int i = 0; i < 6; i++) {
        if (isdof[i]) {
            std::vector<int> findex(1);
            findex[0] = nm++;
            functions1.push_back(new LinearFunction());
            functions2.push_back(new LinearFunction());
            coordIndices.push_back(findex);
        } else {
            std::vector<int> findex(0);
            functions1.push_back(new ConstantFunction());
            functions2.push_back(new ConstantFunction());
            coordIndices.push_back(findex);
        }
    }
    return nm;
}
TEST(Simbody_FunctionBasedMobilizedBodies_PinEquivalenceWithBuiltIn,
     MatchesKinematicsVelocitiesAndAccelerations) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Pin");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set the 1 spatial rotation about Z to be mobility
    isdof[2] = true; // rot Z
    const int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0.5), Inertia(1)));
    MobilizedBody::Pin p1(matter.Ground(), body);
    MobilizedBody::Pin p2(p1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);

    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_PinEquivalenceSkewedAxes,
     MatchesBuiltInPinKinematicsUnderSkewedAxisDefinition) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Pin with skewed axes");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);
    std::vector<Vec3> axes(6);

    // Set the 1 spatial rotation about first axis
    isdof[0] = true; // rot 1
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    double angle = 0;

    axes[0] = Vec3(0, 0, 1);
    axes[1] = Vec3(0, 1, 0);
    axes[2] = Vec3(1, 0, 0);
    axes[3] = Vec3(1, 0, 0);
    axes[4] = Vec3(0, 1, 0);
    axes[5] = Vec3(0, 0, 1);

    Transform inParentPin = Transform(Rotation(angle, YAxis), Vec3(0));
    Transform inChildPin = Transform(Rotation(angle, YAxis), Vec3(0, 1, 0));

    Transform inParentFB = Transform(Vec3(0));
    Transform inChildFB = Transform(Vec3(0, 1, 0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));

    // Built-in
    MobilizedBody::Pin p1(matter.Ground(), inParentPin, body, inChildPin);
    MobilizedBody::Pin p2(p1, inParentPin, body, inChildPin);
    // Function-based
    MobilizedBody::FunctionBased
        fb1(matter.Ground(), inParentFB, body, inChildFB, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParentFB, body, inChildFB, nm, functions2, coordIndices, axes);

    system.realizeTopology();
    compareMobilizedBodies(p2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_SliderEquivalence, MatchesBuiltInSliderKinematics) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Slider");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set the 1 spatial translation along X to be mobility
    isdof[3] = true; // trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), body);
    MobilizedBody::Slider s2(s1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);

    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_SliderEquivalenceSkewedAxes,
     MatchesBuiltInSliderKinematicsUnderSkewedAxisDefinition) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Slider with skewed axes");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);
    std::vector<Vec3> axes(6);

    axes[0] = Vec3(1, 0, 0);
    axes[1] = Vec3(0, 1, 0);
    axes[2] = Vec3(0, 0, 1);
    axes[3] = Vec3(0, 0, 1);
    axes[4] = Vec3(0, 1, 0);
    axes[5] = Vec3(1, 0, 0);

    // Set the 1 spatial translation along X to be mobility
    isdof[5] = true; // trans X
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    Transform inParent = Transform(Vec3(0)); // Transform(Rotation(-Pi/2, YAxis));
    Transform inChild = Transform(Vec3(0));

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Slider s1(matter.Ground(), inParent, body, inChild);
    MobilizedBody::Slider s2(s1, inParent, body, inChild);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices, axes);

    system.realizeTopology();
    compareMobilizedBodies(s2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_CylinderEquivalence, MatchesBuiltInCylinderKinematics) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Cylinder");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set 2 mobilities: rotation about and translation along Z
    isdof[2] = true; // rot Z
    isdof[5] = true; // trans Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Cylinder c1(matter.Ground(), body);
    MobilizedBody::Cylinder c2(c1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);

    system.realizeTopology();
    compareMobilizedBodies(c2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_UniversalEquivalence, MatchesBuiltInUniversalKinematics) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Universal");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set 2 rotation mobilities about body's X then Y
    isdof[0] = true; // rot X
    isdof[1] = true; // rot Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0), Inertia(1)));
    MobilizedBody::Universal u1(matter.Ground(), body);
    MobilizedBody::Universal u2(u1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(u1, body, nm, functions2, coordIndices);

    system.realizeTopology();
    compareMobilizedBodies(u2, fb2, true, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_PlanarEquivalence, MatchesBuiltInPlanarKinematics) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Planar");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[2] = true; // rot Z
    isdof[3] = true; // trans X
    isdof[4] = true; // trans Y
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, 0, 0), Inertia(1)));
    MobilizedBody::Planar u1(matter.Ground(), body);
    MobilizedBody::Planar u2(u1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);

    system.realizeTopology();
    compareMobilizedBodies(u2, fb2, false, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_GimbalEquivalence, MatchesBuiltInGimbalKinematics) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Gimbal");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<bool> isdof(6, false);

    // Set 3 mobilities: Z rotation and translation along body's X then Y
    isdof[0] = true; // rot X
    isdof[1] = true; // rot Y
    isdof[2] = true; // rot Z
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::Gimbal b1(matter.Ground(), body);
    MobilizedBody::Gimbal b2(b1, body);
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();
    compareMobilizedBodies(b2, fb2, true, nm, nm);
}

TEST(Simbody_FunctionBasedMobilizedBodies_GimbalUserAxes_EquivalenceWithBuiltIn,
     MatchesKinematicsVelocitiesAndAccelerationsUnderSkewedAxes) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Gimbal with user-defined axes");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;
    std::vector<Vec3> axes(6);
    std::vector<bool> isdof(6, false);

    isdof[0] = true; // rot 1
    isdof[1] = true; // rot 2
    isdof[2] = true; // rot 3
    int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    // Sherm 20130213: I replaced the random number generator with some firm numbers to prevent singularities
    // from occurring on some platforms based on different random number output.
    axes[0] = Vec3(0.05, 1.4, 0);
    axes[1] = Vec3(0.6, 0, -1.2);
    axes[2] = Vec3(0, 2, -0.055);
    axes[3] = Vec3(1, 0, 0);
    axes[4] = Vec3(0, 1, 0);
    axes[5] = Vec3(0, 0, 1);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    const Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    const Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));

    // Use massless bodies for generationg skewed-axes
    const Body::Massless massLessBody;

    const Transform inParent = Transform(Vec3(0));
    const Transform inChild = Transform(Vec3(0, 1, 0));

    // Compared to standard built-in pin mobilizers with skewed axes
    // Pin rotates about Z-axis and need to align with first axis
    const Transform parentPinAxis0 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0, 0, 0));
    const Transform childPinAxis0 = Transform(Rotation(UnitVec3(axes[0]), ZAxis), Vec3(0, 0, 0));
    const Transform parentPinAxis1 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0, 0, 0));
    const Transform childPinAxis1 = Transform(Rotation(UnitVec3(axes[1]), ZAxis), Vec3(0, 0, 0));
    const Transform parentPinAxis2 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0, 0, 0));
    const Transform childPinAxis2 = Transform(Rotation(UnitVec3(axes[2]), ZAxis), Vec3(0, 1, 0));

    MobilizedBody::Pin masslessPin0(matter.Ground(), parentPinAxis0, massLessBody, childPinAxis0);
    MobilizedBody::Pin masslessPin1(masslessPin0, parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin b1(masslessPin1, parentPinAxis2, body, childPinAxis2);

    MobilizedBody::Pin masslessPin00(b1, parentPinAxis0, massLessBody, childPinAxis0);
    MobilizedBody::Pin masslessPin01(masslessPin00, parentPinAxis1, massLessBody, childPinAxis1);
    MobilizedBody::Pin b2(masslessPin01, parentPinAxis2, body, childPinAxis2);

    MobilizedBody::FunctionBased
        fb1(matter.Ground(), inParent, body, inChild, nm, functions1, coordIndices, axes);
    MobilizedBody::FunctionBased fb2(fb1, inParent, body, inChild, nm, functions2, coordIndices, axes);
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    // These were generated randomly but we want repeatability across machines so we'll use the same numbers
    // every time. Note that we'll re-use each of these twice, once for the pin joint system and once for the
    // function based mobilizers.
    Real initq[] = {1.41292, 0.048025, -1.19474, 0.618909, -0.0552235, 2.043930};
    Real initu[] = {1.53485, 0.546119, -1.55779, -1.872230, 0.0982929, 0.118798};

    const int nq = state.getNQ() / 2;
    EXPECT_LE(nq, 6); // make sure we have enough random numbers!
    for (int i = 0; i < nq; ++i) {
        state.updQ()[i] = state.updQ()[i + nq] = initq[i];
    }

    int nu = state.getNU() / 2;
    for (int i = 0; i < nu; ++i) {
        state.updU()[i] = state.updU()[i + nu] = initu[i];
    }

    system.realize(state, Stage::Acceleration);

    const Transform Xb2 = b2.getBodyTransform(state);
    const Transform Xfb2 = fb2.getBodyTransform(state);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(Xb2, Xfb2);

    const SpatialVec A_b2 = b2.getBodyAcceleration(state);
    const SpatialVec A_fb2 = fb2.getBodyAcceleration(state);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(A_b2, A_fb2);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    const State& result = timeStepper.getState();
    const Vec3& com_bin = b2.getBodyOriginLocation(result);
    const Vec3& com_fb = fb2.getBodyOriginLocation(result);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(com_bin, com_fb);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getBodyVelocity(result), fb2.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}

TEST(Simbody_FunctionBasedMobilizedBodies_SliderEquivalenceWithBuiltIn, MatchesKinematicsUnderSkewedAxes) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Translation");

    // Test against built-in Translation mobilizer
    // for a total of 3 coordinates and 3 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof(6, true);

    // No rotations
    isdof[0] = false; // rot X
    isdof[1] = false; // rot Y
    isdof[2] = false; // rot Z

    const int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    // Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));

    // Built-in mobilized bodies
    MobilizedBody::Translation b1(matter.Ground(), body);
    MobilizedBody::Translation b2(b1, body);

    // Function-based
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);
    MobilizedBody::FunctionBased fb2(fb1, body, nm, functions2, coordIndices);
    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    const int nq = state.getNQ() / 2;
    for (int i = 0; i < nq; ++i) {
        state.updQ()[i] = state.updQ()[i + nq] = random.getValue();
    }

    const int nu = state.getNU() / 2;
    for (int i = 0; i < nu; ++i) {
        state.updU()[i] = state.updU()[i + nu] = random.getValue(); // 0.0; //
    }

    system.realize(state, Stage::Acceleration);

    const Transform& Xb2 = b2.getBodyTransform(state);
    const Transform& Xfb2 = fb2.getBodyTransform(state);

    const SpatialVec& A_b2 = b2.getBodyAcceleration(state);
    const SpatialVec& A_fb2 = fb2.getBodyAcceleration(state);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(A_b2, A_fb2);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    const State& result = timeStepper.getState();

    const Vec3& com_bin = b2.getBodyOriginLocation(result);
    const Vec3& com_fb = fb2.getBodyOriginLocation(result);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(com_bin, com_fb);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getBodyVelocity(result), fb2.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b2.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}

TEST(Simbody_FunctionBasedMobilizedBodies_FreeEquivalenceWithBuiltIn,
     MatchesFull6DOFKinematicsVelocitiesAndAccelerations) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Free");

    // Test against free joint using Euler angles for orientation (q)
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices;
    std::vector<const Function*> functions1;
    std::vector<const Function*> functions2;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof(6, true);

    const int nm = defineMobilizerFunctions(isdof, coordIndices, functions1, functions2);

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0.2, -0.5, 0.1), Inertia(1.2)));

    // Built-in free
    MobilizedBody::Free b1(matter.Ground(), body);

    // Function-based equivalent?
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm, functions1, coordIndices);

    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    int nq = state.getNQ() - nm;
    EXPECT_EQ(nm, state.getNU() / 2);

    // Get random q's and u's and set equivalent on both bodies.
    // (Not really random so we can get repeatability on all platforms.)
    Real initq[] = {0.455189, -0.383271, 1.21353, -0.510623, -1.71438, 0.968387};
    EXPECT_LE(nm, 6);

    // Free has slots for 4 rot q's and fb only has 3
    for (int i = 0; i < nm; ++i) {
        state.updQ()[i] = state.updQ()[i + nq] = initq[i];
    }

    system.realize(state, Stage::Position);
    SpatialVec inputVelocity(Vec3(-0.962157, 0.523767, 1.94993), Vec3(-1.15752, 0.436991, -0.787116));

    b1.setUToFitVelocity(state, inputVelocity);
    fb1.setUToFitVelocity(state, inputVelocity);

    system.realize(state, Stage::Acceleration);

    Transform Xb1 = b1.getBodyTransform(state);
    Transform Xfb1 = fb1.getBodyTransform(state);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(Xb1, Xfb1);

    const SpatialVec& A_b1 = b1.getBodyAcceleration(state);
    const SpatialVec& A_fb1 = fb1.getBodyAcceleration(state);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(A_b1, A_fb1);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    const State& result = timeStepper.getState();

    Xb1 = b1.getBodyTransform(result);
    Xfb1 = fb1.getBodyTransform(result);

    const Vec3& com_bin = b1.getBodyOriginLocation(result);
    const Vec3& com_fb = fb1.getBodyOriginLocation(result);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(com_bin, com_fb);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyVelocity(result), fb1.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(b1.getBodyAcceleration(result), fb1.getBodyAcceleration(result));
}

TEST(Simbody_FunctionBasedMobilizedBodies_FreeVsTranslationGimbal_Equivalence,
     MatchesKinematicsVelocitiesAndAccelerations) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against built-in Translation + Gimbal");

    // Test function-based free against a combination of Translation and Gimbal mobilizer
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices1;
    std::vector<std::vector<int>> coordIndices2a;
    std::vector<std::vector<int>> coordIndices2b;
    std::vector<const Function*> functions1;
    std::vector<const Function*> temp;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof1(6, true);
    std::vector<bool> isdof2a(6, true);
    std::vector<bool> isdof2b(6, true);

    const int nm1 = defineMobilizerFunctions(isdof1, coordIndices1, functions1, temp);
    const int nm2a = 3;
    const int nm2b = 3;

    // Check that we constructed the correct number of functions
    EXPECT_EQ(nm1, nm2a + nm2b);

    // Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(0.1, Vec3(0.25, -0.5, 0.1), Inertia(0.5)));

    // One Free-like function-based mmobilizer
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm1, functions1, coordIndices1);

    // Two function-based mmobilizers: 2a for translation and 2b for rotation
    MobilizedBody::Translation massLessTrans(matter.Ground(), massLessBody);
    MobilizedBody::Gimbal b1(massLessTrans, body);

    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    int nq = state.getNQ() / 2;
    EXPECT_EQ(nq, nm1);

    // Set rotation states first
    for (int i = 0; i < nm2b; ++i) {
        state.updQ()[i] = state.updQ()[i + nq + nm2a] = 0.0; // random.getValue();
        state.updU()[i] = state.updU()[i + nq + nm2a] = random.getValue();
    }

    // Set translations states second
    for (int i = 0; i < nm2a; ++i) {
        state.updQ()[i + nm2a] = state.updQ()[i + nq] = random.getValue();
        state.updU()[i + nm2a] = state.updU()[i + nq] = random.getValue();
    }

    system.realize(state, Stage::Acceleration);

    Transform Xfb1 = fb1.getBodyTransform(state);
    Transform Xb1 = b1.getBodyTransform(state);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(Xfb1, Xb1);

    SpatialVec A_fb1 = fb1.getBodyAcceleration(state);
    SpatialVec A_b1 = b1.getBodyAcceleration(state);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(A_fb1, A_b1);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    const State& result = timeStepper.getState();

    Xfb1 = fb1.getBodyTransform(result);
    Xb1 = b1.getBodyTransform(result);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(Xfb1, Xb1);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb1.getBodyVelocity(result), b1.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb1.getBodyAcceleration(result), b1.getBodyAcceleration(result));
}

TEST(Simbody_FunctionBasedMobilizedBodies_FreeVsTwoFunctionBased_Equivalence,
     MatchesKinematicsVelocitiesAndAccelerations) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody against combination of two FunctionBased mobilizers");

    // Test against free joint that is a combination of Translation and Gimbal mobilizer
    // for a total of 6 coordinates and 6 mobilities

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices1;
    std::vector<std::vector<int>> coordIndices2a;
    std::vector<std::vector<int>> coordIndices2b;
    std::vector<const Function*> functions1;
    std::vector<const Function*> temp;
    std::vector<const Function*> functions2a;
    std::vector<const Function*> functions2b;

    // Set 6 mobilities: rotation and translation about body's X, Y, and then Z axes
    std::vector<bool> isdof1(6, true);
    std::vector<bool> isdof2a(6, true);
    std::vector<bool> isdof2b(6, true);

    // Just translation
    isdof2a[0] = false; // rot X
    isdof2a[1] = false; // rot Y
    isdof2a[2] = false; // rot Z

    // Just rotation
    isdof2b[3] = false; // trans X
    isdof2b[4] = false; // trans Y
    isdof2b[5] = false; // trans Z

    const int nm1 = defineMobilizerFunctions(isdof1, coordIndices1, functions1, temp);
    const int nm2a = defineMobilizerFunctions(isdof2a, coordIndices2a, functions2a, temp);
    const int nm2b = defineMobilizerFunctions(isdof2b, coordIndices2b, functions2b, temp);

    // Check that we constructed the correct number of functions
    EXPECT_EQ(nm1, nm2a + nm2b);

    // Use massless body for translation
    Body::Massless massLessBody;

    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(0.1, Vec3(0.25, -0.5, 0.1), Inertia(0.5)));

    // One Free-like function-based mmobilizer
    MobilizedBody::FunctionBased fb1(matter.Ground(), body, nm1, functions1, coordIndices1);

    // Two function-based mmobilizers: 2a for translation and 2b for rotation
    MobilizedBody::FunctionBased massLessfb(matter.Ground(), massLessBody, nm2a, functions2a, coordIndices2a);
    MobilizedBody::FunctionBased fb2(massLessfb, body, nm2b, functions2b, coordIndices2b);

    system.realizeTopology();

    State state = system.getDefaultState();
    matter.setUseEulerAngles(state, true);
    system.realizeModel(state);

    Random::Gaussian random;

    const int nq = state.getNQ() / 2;
    EXPECT_EQ(nq, nm1);

    // Set rotation states first
    for (int i = 0; i < nm2b; ++i) {
        state.updQ()[i] = state.updQ()[i + nq + nm2a] = random.getValue();
        state.updU()[i] = state.updU()[i + nq + nm2a] = random.getValue();
    }

    // Set translations states second
    for (int i = 0; i < nm2a; ++i) {
        state.updQ()[i + nm2a] = state.updQ()[i + nq] = random.getValue();
        state.updU()[i + nm2a] = state.updU()[i + nq] = random.getValue();
    }

    system.realize(state, Stage::Acceleration);

    const SpatialVec& A_fb1 = fb1.getBodyAcceleration(state);
    const SpatialVec& A_fb2 = fb2.getBodyAcceleration(state);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(A_fb1, A_fb2);

    // Simulate it.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(1.0));

    const State& result = timeStepper.getState();

    const Vec3& com_fb1 = fb1.getBodyOriginLocation(result);
    const Vec3& com_fb2 = fb2.getBodyOriginLocation(result);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(com_fb1, com_fb2);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb1.getBodyVelocity(result), fb2.getBodyVelocity(result));

    // stepTo() only guarantees realization through velocity stage.
    system.realize(result, Stage::Acceleration);
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb1.getBodyAcceleration(result), fb2.getBodyAcceleration(result));
}


/**
 * Test a mobilized body based on functions that take multiple arguments.
 */
TEST(Simbody_FunctionBasedMobilizedBodies_MultipleArgumentFunctions_EvaluatesCorrectKinematics,
     ComputesCorrectTransformAndVelocity) {
    SCOPED_TRACE("Testing FunctionBased MobilizedBody with functions that take multiple arguments");

    // Define the functions that specify the FunctionBased Mobilized Body.
    std::vector<std::vector<int>> coordIndices(6);
    std::vector<const Function*> functions(6);
    Vector coeff(3);
    coeff[0] = 0.5;
    coeff[1] = -0.5;
    coeff[2] = 1.0;
    functions[0] = new Function::Constant(0.0, 0);
    functions[1] = new Function::Constant(0.0, 0);
    functions[2] = new Function::Constant(0.0, 0);
    functions[3] = new NonlinearFunction();
    functions[4] = new Function::Linear(coeff);
    functions[5] = new Function::Constant(0.0, 0);
    coordIndices[3].push_back(0);
    coordIndices[3].push_back(1);
    coordIndices[4].push_back(0);
    coordIndices[4].push_back(1);

    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity gravity(forces, matter, Vec3(0, -9.8, 0));
    Body::Rigid body(MassProperties(1.0, Vec3(0, -0.5, 0), Inertia(0.5)));
    MobilizedBody::FunctionBased fb(matter.Ground(), body, 2, functions, coordIndices);
    State state = system.realizeTopology();

    // See if coordinates and velocities are calculated correctly.
    ASSERT_EQ(state.getNQ(), 2);
    ASSERT_EQ(state.getNU(), 2);

    state.updQ()[0] = 2.0;
    state.updQ()[1] = -3.0;
    state.updU()[0] = 0.1;
    state.updU()[1] = -0.4;
    system.realize(state, Stage::Acceleration);

    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb.getBodyTransform(state), Transform(Vec3(11.0, 3.5, 0.0)));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(fb.getBodyVelocity(state),
                                  SpatialVec(Vec3(0.0), Vec3(0.1 + 3.0 * 0.4, 0.25, 0.0)));

    // Simulate it.
    Real energy = system.calcEnergy(state);
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));
    EXPECT_NO_THROW(timeStepper.stepTo(5.0));
    EXPECT_NEAR_DEFAULT_TOL_SIMTK(energy, system.calcEnergy(timeStepper.getState()));
}
