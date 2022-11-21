/* -------------------------------------------------------------------------- *
 *                               TestSampler                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include <iostream>
#include <numeric>
#include "Sampler.hpp"
#include "Topology.hpp"
#include "HMCSampler.hpp"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#define ASSERT_EQUAL(val1, val2) {ASSERT(std::abs(val1-val2) < 1e-10);}


class DerivedSampler : virtual public Sampler {
public:
    // Constructor
    DerivedSampler (World *argWorld,
            SimTK::CompoundSystem *argCompoundSystem,
            SimTK::SimbodyMatterSubsystem *argMatter,
            std::vector<Topology> &argTopologies,
            SimTK::DuMMForceFieldSubsystem *argDumm,
            SimTK::GeneralForceSubsystem *argForces,
            SimTK::TimeStepper *argTimeStepper)
        : Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper) {}

    ~DerivedSampler(){}
    bool propose(SimTK::State& someState) { return false; }
    void shiftQ(SimTK::State& someState, SimTK::Real scalingFactor, int numIgnoredQs) { }
    void update(SimTK::State& someState){}
};

void CalcAverageAndVariance(const std::vector<SimTK::Real>& v, SimTK::Real& average, SimTK::Real& variance, SimTK::Real& stdev) {
    // Number of elements
    const SimTK::Real N = static_cast<SimTK::Real>(v.size());

    // Average
    const SimTK::Real sum = std::accumulate(v.begin(), v.end(), SimTK::Zero);
    average = sum / N;

    // Variance
    auto variance_func = [&average, &N] (SimTK::Real accumulator, const SimTK::Real& val) {
        const auto diff = val - average;

        // We divide and then multiply in order to mitigate any overflows
        return accumulator + ((diff / (N - 1)) * diff);
    };

    variance =  std::accumulate(v.begin(), v.end(), 0.0, variance_func);
    
    // Standard deviation
    auto stdev_func = [&average, &N] (SimTK::Real accumulator, const SimTK::Real& val) {
        const auto diff = val - average;

        // We divide and then multiply in order to mitigate any overflows
        return accumulator + ((diff / N) * diff);
    };

    stdev =  std::accumulate(v.begin(), v.end(), 0.0, stdev_func);
    stdev = sqrt(stdev);
}

void CalcSkewness(const std::vector<SimTK::Real>& v, SimTK::Real& average, SimTK::Real& b1, SimTK::Real& g1) {
    // See https://en.wikipedia.org/wiki/Skewness#Sample_skewness for more details

    // Number of elements
    const SimTK::Real N = static_cast<SimTK::Real>(v.size());

    // m3
    auto m3_func = [&average, &N] (SimTK::Real accumulator, const SimTK::Real& val) {
        const auto diff = val - average;

        // We divide and then multiply in order to mitigate any overflows
        return accumulator + ((diff * diff * diff) / N);
    };

    const SimTK::Real m3 =  std::accumulate(v.begin(), v.end(), 0.0, m3_func);

    // s3
    auto s3_func = [&average, &N] (SimTK::Real accumulator, const SimTK::Real& val) {
        const auto diff = val - average;

        // We divide and then multiply in order to mitigate any overflows
        return accumulator + ((diff * diff) / (N - 1));
    };

    const SimTK::Real s3 =  std::accumulate(v.begin(), v.end(), 0.0, s3_func);

    // m2
    auto m2_func = [&average, &N] (SimTK::Real accumulator, const SimTK::Real& val) {
        const auto diff = val - average;

        // We divide and then multiply in order to mitigate any overflows
        return accumulator + ((diff * diff) / N);
    };

    const SimTK::Real m2 =  std::accumulate(v.begin(), v.end(), 0.0, m2_func);

    // b1
    b1 = m3 / pow(s3, 1.5);

    // g1
    g1 = m3 / pow(m2, 1.5);
}

void CalcKurtosis() {
    // Kurotsis
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

SimTK::Real percentage(SimTK::Real first, SimTK::Real second) {
    if (first < second) {
        return std::abs(first / first);
    } else {
        return std::abs(second / first);
    }
}

void testSampler(bool Verbose) {

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

    // Create derived sampler. We don't need a world or topologies for testing.
    std::vector<Topology> Topologies;
    Topologies.push_back({});

    DerivedSampler sampler(nullptr, &compoundSystem, &matter, Topologies, &forceField, &forces, &ts);

    // Units check
    sampler.setTemperature(10);
    ASSERT_EQUAL(sampler.getRT(), 83.14472e-3L);

    // Check random number generator
    sampler.setSeed(0);

    constexpr int EXPERIMENTS = 100;
    constexpr int SAMPLES = 15'000 * 5; // samples per molecule * rounds

    constexpr auto B = std::numeric_limits<Sampler::RANDOM_ENGINE_INIT_RESULT_TYPE>::max(),
        A = std::numeric_limits<Sampler::RANDOM_ENGINE_INIT_RESULT_TYPE>::min();
    constexpr __float128 RNG_INIT_AVG = A / 2.0 + B / 2.0;
    constexpr __float128 RNG_INIT_VAR = ((B - A) * (B - A + 2.0)) / 12.0; // TODO is this cross platform?

    // This is where we will store our random numbers
    std::vector<SimTK::Real> v(SAMPLES, 0.0);

    // We store our stats here
    SimTK::Real average, variance, stdev, b1, g1, kurt;

    for (int experiment = 0; experiment < EXPERIMENTS; experiment++) {
        if (Verbose) {
            std::cout << "Running experiment " << experiment << "/" << EXPERIMENTS << "...\n";
        }

        ////////////////////////////////////
        // Check initial RNG - must be on 624, not 10'000
        ////////////////////////////////////
        std::generate(v.begin(), v.end(), [&sampler]() {
            // TODO we calculate on double, but this generates integers
            // In theory, this should not be a problem: https://stackoverflow.com/questions/1848700/biggest-integer-that-can-be-stored-in-a-double
            return static_cast<SimTK::Real>(sampler.randomEngineInit());
        });

        CalcAverageAndVariance(v, average, variance, stdev);

        if (Verbose) {
            std::cout << " - Generated " << SAMPLES << " random uniform integers:\n";
            std::cout << "\tAverage:  " << average << std::endl;
            std::cout << "\tVariance: " << variance << std::endl;
        }
        
        ASSERT( percentage(average, RNG_INIT_AVG) > 0.95 );
        ASSERT( percentage(variance, static_cast<SimTK::Real>(RNG_INIT_VAR)) > 0.95 );
        ASSERT( KolmogorovSmirnov(v) > 0.95 );

        ////////////////////////////////////
        // Check uniform distribution
        ////////////////////////////////////
        std::generate(v.begin(), v.end(), [&sampler]() {
            return sampler.generateRandomNumber(GmolRandDistributionType::UNIFORM);
        });

        CalcAverageAndVariance(v, average, variance, stdev);

        if (Verbose) {
            std::cout << " - Generated " << SAMPLES << " random uniform reals:\n";
            std::cout << "\tAverage:  " << average << std::endl;
            std::cout << "\tVariance: " << variance << std::endl;
        }

        ASSERT( percentage(average, 0.5) > 0.99 );
        ASSERT( percentage(variance, 1.0 / 12.0) > 0.99 );
        ASSERT( KolmogorovSmirnov(v) > 0.95 );

        ////////////////////////////////////
        // Check normal distribution
        ////////////////////////////////////
        std::generate(v.begin(), v.end(), [&sampler]() {
            return sampler.generateRandomNumber(GmolRandDistributionType::NORMAL);
        });

        CalcAverageAndVariance(v, average, variance, stdev);
        CalcSkewness(v, average, b1, g1);

        if (Verbose) {
            std::cout << " - Generated " << SAMPLES << " random Gaussian numbers:\n";
            std::cout << "\tAverage:  " << average << std::endl;
            std::cout << "\tVariance: " << variance << std::endl;
            std::cout << "\tSkewness:\n";
            std::cout << "\t - b1: " << b1 << "\n\t - g1: " << g1 << std::endl;
        }

        ASSERT( std::abs(average - 0.0) < 0.05 );
        ASSERT( std::abs(variance - 1.0) < 0.05 );
        ASSERT( std::abs(b1) < 0.05 );
        ASSERT( std::abs(g1) < 0.05 );
        // ASSERT(b1 < g1);
    }
}

int main(int argc, char* argv[]) {
    try {
        bool Verbose = false;

        if (argc == 2) {
            auto arg = std::string(argv[1]);
            if (arg == "--verbose") {
                Verbose = true;
            }
        }

        testSampler(Verbose);
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        std::cout << "Run with --verbose to see what is wrong." << std::endl;
        std::cout << "Keep in mind that tests like this can fail occasionally. ";
        std::cout << "Run this a couple more times before saying concluding that there is something wrong here." << std::endl;

        return 1;
    }
    std::cout << "Done" << std::endl;

    return 0;
}
