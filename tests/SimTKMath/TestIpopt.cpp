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

using SimTK::Matrix;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;
using SimTK::Real;
using SimTK::Vector;

static int NUMBER_OF_PARAMETERS = 4;
static int NUMBER_OF_EQUALITY_CONSTRAINTS = 1;
static int NUMBER_OF_INEQUALITY_CONSTRAINTS = 1;

class ProblemSystem : public OptimizerSystem {
    public:
    int objectiveFunc(const Vector& coefficients, bool, Real& f) const override {
        const Real* x = &coefficients[0];
        f = x[0] * x[3] * (x[0] + x[1] + x[2]) + x[2];
        return 0;
    }

    int gradientFunc(const Vector& coefficients, bool, Vector& gradient) const override {
        const Real* x = &coefficients[0];
        gradient[0] = x[0] * x[3] + x[3] * (x[0] + x[1] + x[2]);
        gradient[1] = x[0] * x[3];
        gradient[2] = x[0] * x[3] + 1;
        gradient[3] = x[0] * (x[0] + x[1] + x[2]);
        return 0;
    }

    int constraintFunc(const Vector& coefficients, bool, Vector& constraints) const override {
        const Real* x = &coefficients[0];
        constraints[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] - 40.0;
        constraints[1] = x[0] * x[1] * x[2] * x[3] - 25.0;
        return 0;
    }

    int constraintJacobian(const Vector& coefficients, bool, Matrix& jac) const override {
        const Real* x = &coefficients[0];
        jac(0, 0) = 2 * x[0];
        jac(0, 1) = 2 * x[1];
        jac(0, 2) = 2 * x[2];
        jac(0, 3) = 2 * x[3];
        jac(1, 0) = x[1] * x[2] * x[3];
        jac(1, 1) = x[0] * x[2] * x[3];
        jac(1, 2) = x[0] * x[1] * x[3];
        jac(1, 3) = x[0] * x[1] * x[2];
        return 0;
    }

    ProblemSystem(int numParams, int numEq, int numIneq)
        : OptimizerSystem(numParams) {
        setNumEqualityConstraints(numEq);
        setNumInequalityConstraints(numIneq);
    }
};

TEST(SimbodyOptimizer, IpoptHS071) {
    SCOPED_TRACE("Testing Ipopt HS071 equivalent problem\n");

    ProblemSystem sys(NUMBER_OF_PARAMETERS, NUMBER_OF_EQUALITY_CONSTRAINTS, NUMBER_OF_INEQUALITY_CONSTRAINTS);

    Vector results(NUMBER_OF_PARAMETERS);
    Vector lower_bounds(NUMBER_OF_PARAMETERS);
    Vector upper_bounds(NUMBER_OF_PARAMETERS);

    results[0] = 1.0;
    results[1] = 5.0;
    results[2] = 5.0;
    results[3] = 1.0;

    for (int i = 0; i < NUMBER_OF_PARAMETERS; ++i) {
        lower_bounds[i] = 1.0;
        upper_bounds[i] = 5.0;
    }

    sys.setParameterLimits(lower_bounds, upper_bounds);

    Real f = 0;

    EXPECT_NO_THROW({
        Optimizer opt(sys);

        opt.setConvergenceTolerance(1e-4);
        opt.setDiagnosticsLevel(7);
        opt.setLimitedMemoryHistory(500);

        opt.setAdvancedBoolOption("warm_start", true);
        opt.setAdvancedRealOption("obj_scaling_factor", 1);
        opt.setAdvancedRealOption("nlp_scaling_max_gradient", 1);

        f = opt.optimize(results);
    });

    static const Real TOL = 1e-4;
    Real expected[] = {1.00000000, 4.74299963, 3.82114998, 1.37940829};

    for (int i = 0; i < NUMBER_OF_PARAMETERS; ++i) {
        SCOPED_TRACE("Checking parameter index " + std::to_string(i) + "\n");
        EXPECT_NEAR(results[i], expected[i], TOL)
            << "result=" << results[i] << " expected=" << expected[i] << "\n";
    }
}
