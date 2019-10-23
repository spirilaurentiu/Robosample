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
    bool update(SimTK::State& someState){}
};

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

    // TODO Kolmogorov-Smirnov test
    // Check random number generator
    sampler.setSeed(std::time(0));
    std::vector<SimTK::Real> X;
    std::vector<SimTK::Real>::iterator XIt;
    unsigned int sizeOfSample = 0;
    SimTK::Real variance;
    SimTK::Real average;
    SimTK::Real skew;

    for (unsigned int experiment = 0; experiment < 5; experiment++){

        ///////////////////////////
        // Check uniform distribution
        sizeOfSample = 10000;
        for(unsigned int i = 0; i < sizeOfSample; i++){
            X.emplace_back(sampler.generateRandomNumber(UNIFORM));
        }

        // Average
        average = 0;
        SimTK::Real average = static_cast<SimTK::Real>(std::accumulate
                (X.begin(), X.end(), 0.0)) / (X.size());

        // Variance
        SimTK::Real variance = 0.0;
        for(XIt = X.begin(); XIt != X.end(); ++XIt){
            //std::cout << *XIt << std::endl;
            variance += ((*XIt - average) * (*XIt - average));
        }
        variance /= X.size();

        ASSERT( std::abs(average - 0.5) < 0.01 );
        ASSERT( std::abs(variance - 0.0833) < 0.01 );

        ///////////////////////////
        // Check normal distribution
        sizeOfSample = 100000;
        X.clear();

        // Average
        average = 0;
        for(unsigned int i = 0; i < sizeOfSample; i++){
            X.emplace_back(sampler.generateRandomNumber(NORMAL));
        }

        average = static_cast<SimTK::Real>(std::accumulate
                (X.begin(), X.end(), 0.0)) / (X.size());

        // Variance
        variance = 0.0;
        for(XIt = X.begin(); XIt != X.end(); ++XIt){
            variance += ((*XIt - average) * (*XIt - average));
        }
        variance /= X.size();
        SimTK::Real stdev = std::sqrt(variance);

        // Skewness
        skew = 0;
        for(XIt = X.begin(); XIt != X.end(); ++XIt){
            skew += (*XIt * *XIt * *XIt);
        }
        skew /= X.size();

        skew = skew - (3 * average * variance) - (average * average * average);
        skew /= variance * stdev;

        //std::cout << "Average: " << average << std::endl;
        //std::cout << "Variance: " << variance << std::endl;
        //std::cout << "Skewness: " << skew << std::endl;

        ASSERT( std::abs(average - 0.0) < 0.05 );
        ASSERT( std::abs(variance - 1.0) < 0.05 );
        ASSERT( std::abs(skew - 0.0) < 0.05 );

        X.clear();

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
