#include <gtest/gtest.h>

#include "SimTKsimbody.h"

using namespace SimTK;
using namespace std;

class HarmonicOscillator {
    public:
    class OscillatorReporter : public PeriodicEventReporter {
        public:
        mutable int eventCount{};
        mutable Real sumEnergy{};
        mutable Real sumEnergySquared{};
        mutable Real sumVelocity{};
        mutable Real sumAbsVelocity{};
        mutable Real sumVelocitySquared{};
        mutable Real sumPosition{};
        mutable Real sumRMSVelPos{};

        OscillatorReporter(HarmonicOscillator& oscillator, Real reportInterval)
            : PeriodicEventReporter(reportInterval)
            , oscillator(oscillator) {
        }

        void handleEvent(const State& state) const override {
            // Equilibrate a bit before collecting data
            if (state.getTime() <= 1) {
                return;
            }

            ++eventCount;
            oscillator.updSystem().realize(state, Stage::Dynamics);

            Real energy = oscillator.getSystem().calcKineticEnergy(state);
            sumEnergy += energy;
            sumEnergySquared += energy * energy;

            Real position = oscillator.getPosition(state);
            sumPosition += position;

            Real velocity = oscillator.getVelocity(state);
            sumVelocity += velocity;
            sumAbsVelocity += std::abs(velocity);
            sumVelocitySquared += velocity * velocity;

            sumRMSVelPos += std::sqrt(position * position * velocity * velocity);
        }

        private:
        HarmonicOscillator& oscillator;
    };

    HarmonicOscillator()
        : matter(system)
        , forces(system)
        , mass(1.0) {
        Vec3 station(0.0);

        // Use a slider to constrain the oscillator to one dimension of motion
        MobilizedBody::Slider body(matter.updGround(),
                                   Body::Rigid(MassProperties(mass, station, Inertia(1))));
        body.setDefaultLength(-2.0); // initial position at -2
        sliderIndex = body.getMobilizedBodyIndex();

        // Unit spring constant to match example in Frenkel and Smit
        Force::TwoPointLinearSpring(forces, matter.getGround(), Vec3(0), body, Vec3(0), 1.0, 0.0);
    }

    void simulate() {
        reporter = new OscillatorReporter(*this, 0.1);
        system.addEventReporter(reporter);

        State state = system.realizeTopology();
        Random::Uniform rand(-1, 1);

        // Simulate it.
        VerletIntegrator integ(system);
        // RungeKuttaMersonIntegrator integ(system);
        // integ.setAccuracy(0.01);

        TimeStepper timeStepper(system, integ);
        EXPECT_NO_THROW(timeStepper.initialize(state));
        EXPECT_NO_THROW(timeStepper.stepTo(150.0));
    }

    void assertTemperature(Real temperature) const {
        // ensure we collected some data
        EXPECT_GT(reporter->eventCount, 100);

        int degreesOfFreedom = 1;

        // Mean position should be zero
        const Real expectedMeanPosition = 0.0;
        const Real measuredMeanPosition = reporter->sumPosition / reporter->eventCount;
        EXPECT_LT(std::abs(expectedMeanPosition - measuredMeanPosition), 0.2);

        // Mean velocity should be zero
        const Real expectedMeanVelocity = 0.0;
        const Real measuredMeanVelocity = reporter->sumVelocity / reporter->eventCount;
        EXPECT_LT(std::abs(expectedMeanVelocity - measuredMeanVelocity), 0.2);

        // Check temperature
        const Real measuredMeanEnergy = reporter->sumEnergy / reporter->eventCount;
        const Real expectedMeanEnergy =
            degreesOfFreedom * 0.5 * SimTK_BOLTZMANN_CONSTANT_MD * temperature; // kT/2 per degree of freedom
        EXPECT_LT(std::abs(1.0 - (measuredMeanEnergy / expectedMeanEnergy)), 0.2);

        // Boltzmann distribution stuff

        // Mean squared velocity should be dof*kT/mass
        // Boltzmann distribution
        const Real expectedMeanVelocitySquared =
            degreesOfFreedom * SimTK_BOLTZMANN_CONSTANT_MD * temperature / mass;
        const Real measuredMeanVelocitySquared = reporter->sumVelocitySquared / reporter->eventCount;
        EXPECT_LT(std::abs(1.0 - measuredMeanVelocitySquared / expectedMeanVelocitySquared), 0.2);

        // TODO: check this formula
        // Mean absolute velocity should be (8*v2bar/3PI)^1/2
        const Real expectedMeanAbsVelocity =
            std::sqrt(8.0 * expectedMeanVelocitySquared / (degreesOfFreedom * SimTK_PI));
        const Real measuredMeanAbsVelocity = reporter->sumAbsVelocity / reporter->eventCount;
        // ASSERT(std::abs(1.0 - measuredMeanAbsVelocity/expectedMeanAbsVelocity) < 0.2);
    }

    auto updSystem() -> MultibodySystem& {
        return system;
    }

    [[nodiscard]] auto getSystem() const -> const MultibodySystem& {
        return system;
    }

    SimbodyMatterSubsystem& updMatterSubsystem() {
        return matter;
    }

    [[nodiscard]] auto getMatterSubsystem() const -> const SimbodyMatterSubsystem& {
        return matter;
    }

    auto updForceSubsystem() -> GeneralForceSubsystem& {
        return forces;
    }

    [[nodiscard]] auto getForceSubsystem() const -> const GeneralForceSubsystem& {
        return forces;
    }

    // slider coordinate
    [[nodiscard]] auto getPosition(const State& state) const -> Real {
        const MobilizedBody::Slider& slider =
            MobilizedBody::Slider::downcast(matter.getMobilizedBody(sliderIndex));

        return slider.getLength(state);
    }

    [[nodiscard]] auto getVelocity(const State& state) const -> Real {
        const MobilizedBody::Slider& slider =
            MobilizedBody::Slider::downcast(matter.getMobilizedBody(sliderIndex));

        return slider.getRate(state);
    }

    static auto getTime(const State& state) -> Real {
        return state.getTime();
    }

    private:
    MultibodySystem system;
    SimbodyMatterSubsystem matter;
    GeneralForceSubsystem forces;

    Real mass;
    MobilizedBodyIndex sliderIndex;

    OscillatorReporter* reporter;
};


// Case study 12, page 155 in
// Understanding Molecular Simulation: From Algorithms to Applications
// Frenkel and Smit
TEST(Simbody_NoseHooverThermostat_HarmonicOscillatorNoThermostat, RunsWithoutErrors) {
    SCOPED_TRACE("testHarmonicOscillatorNoThermostat");

    HarmonicOscillator oscillator;
    oscillator.simulate();
}

TEST(Simbody_NoseHooverThermostat_ConstructorSmoke, ConstructsAndRuns) {
    SCOPED_TRACE("testNoseHooverConstructorSmoke");

    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat(forces, oscillator.getMatterSubsystem(), SimTK_BOLTZMANN_CONSTANT_MD, 300, .1);
    oscillator.simulate();
}

TEST(Simbody_NoseHooverThermostat_Temperature100K, MaintainsTargetTemperature) {
    SCOPED_TRACE("oscillator 100K");

    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat nhc(forces, oscillator.getMatterSubsystem(), SimTK_BOLTZMANN_CONSTANT_MD, 100, 0.1);
    oscillator.simulate();
    oscillator.assertTemperature(100);
}

TEST(Simbody_NoseHooverThermostat_Temperature300K, MaintainsTargetTemperature) {
    SCOPED_TRACE("oscillator 300K");

    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat nhc(forces, oscillator.getMatterSubsystem(), SimTK_BOLTZMANN_CONSTANT_MD, 300.0, 0.1);
    oscillator.simulate();
    oscillator.assertTemperature(300.0);
}

TEST(Simbody_NoseHooverThermostat_Temperature5000K, MaintainsTargetTemperature) {
    SCOPED_TRACE("oscillator 5000K");

    HarmonicOscillator oscillator;
    GeneralForceSubsystem& forces = oscillator.updForceSubsystem();
    Force::Thermostat nhc(forces, oscillator.getMatterSubsystem(), SimTK_BOLTZMANN_CONSTANT_MD, 5000.0, 0.1);
    oscillator.simulate();
    oscillator.assertTemperature(5000.0);
}
