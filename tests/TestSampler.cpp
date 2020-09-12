/* -------------------------------------------------------------------------- *
 *                               TestSampler                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include <iostream>
#include <numeric>
#include "Sampler.hpp"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#define ASSERT_EQUAL(val1, val2) {ASSERT(std::abs(val1-val2) < 1e-10);}

class DerivedSampler : virtual public Sampler{
public:
    // Constructor
    DerivedSampler (SimTK::CompoundSystem *argCompoundSystem,
                    SimTK::SimbodyMatterSubsystem *argMatter,
                    SimTK::Compound *argResidue,
                    SimTK::DuMMForceFieldSubsystem *argDumm,
                    SimTK::GeneralForceSubsystem *argForces,
                    SimTK::TimeStepper *argTimeStepper)
    : Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper){}

    ~DerivedSampler(){}
    void propose(SimTK::State& someState){}
    void update(SimTK::State& someState){}
};

void CalcAverageAndVariance(const std::vector<SimTK::Real>& v, SimTK::Real& average, SimTK::Real& variance) {
    // Average
    average = static_cast<SimTK::Real>(std::accumulate(v.begin(), v.end(), 0.0)) / (v.size());

    // Variance
    variance = 0.0;
    for(const auto it : v) {
        variance += ((it - average) * (it - average));
    }
    
    variance /= v.size();
}

SimTK::Real KolmogorovSmirnov(const std::vector<SimTK::Real>& v) {
    SimTK::Real res = std::numeric_limits<SimTK::Real>::min();
    for(size_t i = 0; i < v.size(); i++) {
        auto DPlus = (i + 1) / static_cast<SimTK::Real>(v.size()) - v[i],
            DMinus = v[i] - (i - 1) / static_cast<SimTK::Real>(v.size());

        res = std::max({DPlus, DMinus, res});
    }

    return res;
}

void testSampler(){

    // Molmodel System derived from Simbody System
    SimTK::CompoundSystem compoundSystem;

    // Simbody subsystems (Minimum requirements)
    SimTK::SimbodyMatterSubsystem matter(compoundSystem);
    SimTK::GeneralForceSubsystem forces(compoundSystem);

    // Initialize Molmodel default ForceSubsystem (DuMM)
    SimTK::DuMMForceFieldSubsystem forceField(compoundSystem);

    // Intialize an integrator and a TimeStepper to manage it
    SimTK::VerletIntegrator integ(compoundSystem);
    SimTK::TimeStepper ts(compoundSystem, integ);

    // Empty Compound
    SimTK::Compound compound;

    DerivedSampler sampler(&compoundSystem, &matter, &compound,
                    &forceField, &forces, &ts);

    // Units check
    sampler.setTemperature(10);
    ASSERT_EQUAL(sampler.getRT(), 83.14472e-3L);

    // Check random number generator
    sampler.setSeed(std::time(0));
    std::vector<SimTK::Real> v;

    for (unsigned int experiment = 0; experiment < 5; experiment++) {
        // Check uniform distribution
        for(unsigned int i = 0; i < 10000; i++) {
            v.emplace_back(sampler.generateRandomNumber(GmolRandDistributionType::UNIFORM));
        }

        SimTK::Real average, variance;
        CalcAverageAndVariance(v, average, variance);

        ASSERT( std::abs(average - 0.5) < 0.01 );
        ASSERT( std::abs(variance - 0.0833) < 0.01 );
        ASSERT( KolmogorovSmirnov(v) > 0.95 );

        // Check normal distribution
        v.clear();
        for(unsigned int i = 0; i < 100000; i++){
            v.emplace_back(sampler.generateRandomNumber(GmolRandDistributionType::NORMAL));
        }

        CalcAverageAndVariance(v, average, variance);

        // Standard deviation
        SimTK::Real stdev = std::sqrt(variance);

        // Skewness
        SimTK::Real skew = 0;
        for(const auto it : v ) {
            skew += (it * it * it);
        }
        skew /= v.size();

        skew = skew - (3 * average * variance) - (average * average * average);
        skew /= variance * stdev;

        //std::cout << "Average: " << average << std::endl;
        //std::cout << "Variance: " << variance << std::endl;
        //std::cout << "Skewness: " << skew << std::endl;

        ASSERT( std::abs(average - 0.0) < 0.05 );
        ASSERT( std::abs(variance - 1.0) < 0.05 );
        ASSERT( std::abs(skew - 0.0) < 0.05 );

        v.clear();
    }
}

int main() {
    try {
        testSampler();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
