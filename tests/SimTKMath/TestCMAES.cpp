/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-14 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
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
#include <iostream>

#include "OptimizerSystems.h"
#include "SimTKmath.h"

using SimTK::Optimizer;
using SimTK::OptimizerSystem;
using SimTK::Real;
using SimTK::Vector;
using std::cout;
using std::endl;

auto vectorsAreEqual(const Vector& actual, const Vector& expected, double tol, bool printWhenNotEqual = false)
    -> bool {
    bool isEqual = true;
    for (unsigned int i = 0; i < actual.size(); ++i) {
        if (!SimTK::Test::numericallyEqual(actual[i], expected[i], 1, tol)) {
            if (printWhenNotEqual) {
                printf("error actual[%d] = %f  expected[%d] = %f \n", i, actual[i], i, expected[i]);
            }
            isEqual = false;
        }
    }
    return isEqual;
}

#define EXPECT_OPT(opt, results, tol)                                            \
    do {                                                                         \
        Real funval = opt.optimize(results);                                     \
        const TestOptimizerSystem& sys =                                         \
            *static_cast<const TestOptimizerSystem*>(&opt.getOptimizerSystem()); \
        bool passed = vectorsAreEqual(results, sys.optimalParameters(), tol);    \
        if (!SimTK::Test::numericallyEqual(funval, sys.optimalValue(), 1, tol))  \
            passed = false;                                                      \
        EXPECT_TRUE(passed);                                                     \
    } while (false)

/* --- Tests --- */

TEST(CMAES, Available) {
    EXPECT_TRUE(Optimizer::isAlgorithmAvailable(SimTK::CMAES));
}

TEST(CMAES, RequiresAtLeastTwoParameters) {
    EXPECT_THROW(Optimizer opt(Cigtab(1), SimTK::CMAES), SimTK::Exception::ValueOutOfRange);
}

TEST(CMAES, MaxIterations) {
    Cigtab sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(0.5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setAdvancedRealOption("init_stepsize", 0.3);
    opt.setMaxIterations(500);

    opt.optimize(results);

    EXPECT_FALSE(vectorsAreEqual(results, sys.optimalParameters(), 1e-5));
}

TEST(CMAES, CigtabOptimum) {
    Cigtab sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(0.5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("init_stepsize", 0.3);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    EXPECT_OPT(opt, results, 1e-5);
}

TEST(CMAES, ParameterLimits) {
    Easom sys;
    Vector results(sys.getNumParameters());
    results.setTo(100);

    Optimizer opt(sys, SimTK::CMAES);
    opt.optimize(results);

    results.setTo(100.01);
    EXPECT_THROW(opt.optimize(results), SimTK::Exception::APIArgcheckFailed);
}

TEST(CMAES, SigmaAndAckley) {
    Ackley sys(2);
    Vector results(sys.getNumParameters());
    results.setTo(25);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    opt.setAdvancedIntOption("popsize", 50);
    opt.setAdvancedIntOption("seed", 30);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    Real f1 = opt.optimize(results);

    Vector expected(2, 24.999749);
    EXPECT_TRUE(vectorsAreEqual(results, expected, 1e-5));

    opt.setAdvancedRealOption("init_stepsize", 0.5 * 64);
    results.setTo(25);

    EXPECT_OPT(opt, results, 1e-5);
}

TEST(CMAES, DropWave) {
    DropWave sys;
    Vector results(sys.getNumParameters());
    results.setTo(2);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-5);
    opt.setMaxIterations(5000);
    opt.setAdvancedRealOption("init_stepsize", 3.5);
    opt.setAdvancedIntOption("popsize", 1000);
    opt.setAdvancedIntOption("stopMaxFunEvals", 100000);
    opt.setAdvancedIntOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    EXPECT_OPT(opt, results, 1e-2);
}

TEST(CMAES, MaxFunEvals) {
    Cigtab sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setAdvancedRealOption("init_stepsize", 0.3);
    opt.setAdvancedIntOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    opt.setAdvancedIntOption("stopMaxFunEvals", 1);
    opt.optimize(results);

    EXPECT_FALSE(vectorsAreEqual(results, sys.optimalParameters(), 1e-4, false));

    opt.setAdvancedIntOption("stopMaxFunEvals", 100000);
    results.setTo(5);

    EXPECT_OPT(opt, results, 1e-4);
}

TEST(CMAES, Seed) {
    Ackley sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(25);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setAdvancedRealOption("init_stepsize", 1);

    opt.setAdvancedIntOption("seed", -10);
    EXPECT_THROW(opt.optimize(results), SimTK::Exception::ValueWasNegative);

    opt.setMaxIterations(100);
    opt.setAdvancedIntOption("seed", 42);

    Real f1 = opt.optimize(results);
    Vector r1 = results;

    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    results.setTo(25);
    Real f2 = opt.optimize(results);
    Vector r2 = results;

    results.setTo(25);
    Real f3 = opt.optimize(results);
    Vector r3 = results;

    EXPECT_NEAR(f2, f3, 1e-10);
    EXPECT_TRUE(vectorsAreEqual(r2, r3, 1e-10));

    opt.setAdvancedIntOption("seed", 50);
    results.setTo(25);
    Real f4 = opt.optimize(results);
    Vector r4 = results;

    EXPECT_FALSE(vectorsAreEqual(r2, r4, 1e-10, false));
}

TEST(CMAES, ConvergenceTolerance) {
    Cigtab sys(2);
    Vector results(sys.getNumParameters());
    results.setTo(5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setAdvancedIntOption("seed", 10);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    opt.setConvergenceTolerance(0.001);
    Real f = opt.optimize(results);
    EXPECT_NEAR(f, sys.optimalValue(), 0.001);

    opt.setConvergenceTolerance(1e-10);
    results.setTo(5);
    f = opt.optimize(results);
    EXPECT_NEAR(f, sys.optimalValue(), 1e-10);
}

TEST(CMAES, Rosenbrock) {
    Rosenbrock sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(0.5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(100000);
    opt.setAdvancedRealOption("init_stepsize", 0.3);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    EXPECT_OPT(opt, results, 1e-6);
}

TEST(CMAES, Schwefel) {
    Schwefel sys(4);
    Vector results(sys.getNumParameters());
    results.setTo(200);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-4);
    opt.setAdvancedIntOption("popsize", 200);
    opt.setAdvancedRealOption("init_stepsize", 300);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    EXPECT_OPT(opt, results, 1e-4);
}

TEST(CMAES, Easom) {
    Easom sys;
    Vector results(sys.getNumParameters());
    results.setTo(-10);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setAdvancedIntOption("popsize", 500);
    opt.setAdvancedRealOption("init_stepsize", 25);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    EXPECT_OPT(opt, results, 1e-5);
}

TEST(CMAES, StopFitness) {
    Ackley sys(2);
    Vector results(sys.getNumParameters());
    results.setTo(25);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    opt.setAdvancedIntOption("popsize", 50);
    Vector step(2, 0.5 * 64);
    opt.setAdvancedVectorOption("init_stepsize", step);
    opt.setAdvancedIntOption("seed", 30);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);

    opt.setAdvancedRealOption("stopFitness", 5);
    Real f = opt.optimize(results);

    EXPECT_GT(f, 0.01);
}

TEST(CMAES, Multithreading) {
    Cigtab sys(22);
    Vector results(sys.getNumParameters());
    results.setTo(0.5);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setConvergenceTolerance(1e-12);
    opt.setMaxIterations(5000);
    Vector step(sys.getNumParameters(), 0.3);
    opt.setAdvancedVectorOption("init_stepsize", step);
    opt.setAdvancedIntOption("seed", 42);
    opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
    opt.setAdvancedStrOption("parallel", "multithreading");

    EXPECT_OPT(opt, results, 1e-5);

    opt.setAdvancedIntOption("nthreads", 2);
    EXPECT_OPT(opt, results, 1e-5);
}

TEST(CMAES, InitStepSizeConflict) {
    Easom sys;
    Vector results(sys.getNumParameters());
    results.setTo(-10);

    Optimizer opt(sys, SimTK::CMAES);
    opt.setAdvancedRealOption("init_stepsize", 25);
    Vector step(sys.getNumParameters(), 25);
    opt.setAdvancedVectorOption("init_stepsize", step);

    EXPECT_THROW(opt.optimize(results), std::logic_error);
}