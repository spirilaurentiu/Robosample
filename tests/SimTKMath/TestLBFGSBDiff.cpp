/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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

#include "SimTKmath.h"

using SimTK::Optimizer;
using SimTK::OptimizerSystem;
using SimTK::Real;
using SimTK::Vector;

const static int NUMBER_OF_PARAMETERS = 25;

class ProblemSystem : public OptimizerSystem {
    public:
    ProblemSystem(int numParameters)
        : OptimizerSystem(numParameters) {
    }

    int objectiveFunc(const Vector& coefficients, bool, Real& f) const override {
        const Real* x = &coefficients[0];

        f = .25 * (x[0] - 1.0) * (x[0] - 1.0);
        for (int i = 1; i < getNumParameters(); i++) {
            f += std::pow(x[i] - x[i - 1] * x[i - 1], 2.0);
        }
        f = 4.0 * f;
        return 0;
    }

    int gradientFunc(const Vector& coefficients, bool, Vector& gradient) const override {
        const Real* x = &coefficients[0];
        Real t1, t2;

        t1 = x[1] - (x[0] * x[0]);
        gradient[0] = 2.0 * (x[0] - 1.0) - 16.0 * x[0] * t1;

        for (int i = 1; i < getNumParameters() - 1; i++) {
            t2 = t1;
            t1 = x[i + 1] - (x[i] * x[i]);
            gradient[i] = 8.0 * t2 - 16.0 * x[i] * t1;
        }

        gradient[getNumParameters() - 1] = 8.0 * t1;
        return 0;
    }
};

static bool equalToTol(Real v1, Real v2, Real tol) {
    const Real scale = std::max(std::max(std::abs(v1), std::abs(v2)), Real(1));
    return std::abs(v1 - v2) < scale * tol;
}

TEST(SimbodyOptimizer, LBFGSBDiffTest) {
    SCOPED_TRACE("Testing LBFGSB with numerical gradient\n");

    const int n = NUMBER_OF_PARAMETERS;

    Vector results(n);
    Vector lower_bounds(n);
    Vector upper_bounds(n);

    ProblemSystem sys(n);

    for (int i = 0; i < n; i++) {
        results[i] = 3.0;
    }

    for (int i = 0; i < n; i += 2) {
        lower_bounds[i] = 1.0;
        upper_bounds[i] = 100.0;
    }
    for (int i = 1; i < n; i += 2) {
        lower_bounds[i] = -100.0;
        upper_bounds[i] = 100.0;
    }

    sys.setParameterLimits(lower_bounds, upper_bounds);

    Real f = 0;

    EXPECT_NO_THROW({
        Optimizer opt(sys);
        opt.setConvergenceTolerance(.0001);
        opt.useNumericalGradient(true);
        f = opt.optimize(results);
    });

    static const Real TOL = 1e-3;
    Real expected[] = {1.000000, 0.999998, 1.000000, 1.000001, 1.000003, 1.000006, 1.000007,
                       1.000012, 1.000022, 1.000040, 1.000081, 1.000161, 1.000325, 1.000650,
                       1.001302, 1.002603, 1.005214, 1.010450, 1.021013, 1.042466, 1.086736,
                       1.180997, 1.394759, 1.945352, 3.784388};

    for (int i = 0; i < n; i++) {
        SCOPED_TRACE("Checking parameter index " + std::to_string(i) + "\n");
        EXPECT_TRUE(equalToTol(results[i], expected[i], TOL))
            << "result=" << results[i] << " expected=" << expected[i] << "\n";
    }
}
