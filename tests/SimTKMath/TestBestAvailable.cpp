/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
 * Authors: Chris Dembia                                                      *
 * Contributors: Jack Middleton                                               *
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

#include "SimTKmath.h"
using SimTK::Matrix;
using SimTK::Optimizer;
using SimTK::OptimizerSystem;
using SimTK::Real;
using SimTK::Vector;
using std::cout;
using std::endl;

/* This test doesn't actually run any optimizations. It simply creates an
 * Optimizer, using the same OptimizerSystem's (not the EXACT same) used in the
 * tests for each of the algorithms. NOTE: Since BestAvailable never selects
 * CFSQP, we do not test for CFSQP.
 */

/* See IpoptTest. 'BestAvailable' will select Ipopt because this sytem has
 * constraints.
 */

class IpoptSystem : public OptimizerSystem {
    public:
    auto objectiveFunc(const Vector& coefficients, bool, Real& f) const -> int override {
        const Real* x = &coefficients[0];
        f = (x[0] * x[3] * (x[0] + x[1] + x[2])) + x[2];
        return 0;
    }

    auto gradientFunc(const Vector& coefficients, bool, Vector& gradient) const -> int override {
        const Real* x = &coefficients[0];
        gradient[0] = (x[0] * x[3]) + (x[3] * (x[0] + x[1] + x[2]));
        gradient[1] = x[0] * x[3];
        gradient[2] = (x[0] * x[3]) + 1;
        gradient[3] = x[0] * (x[0] + x[1] + x[2]);
        return 0;
    }

    auto constraintFunc(const Vector& coefficients, bool, Vector& constraints) const -> int override {
        const Real* x = &coefficients[0];
        constraints[0] = (x[0] * x[0]) + (x[1] * x[1]) + (x[2] * x[2]) + (x[3] * x[3]) - 40.0;
        constraints[1] = (x[0] * x[1] * x[2] * x[3]) - 25.0;
        return 0;
    }

    auto constraintJacobian(const Vector& coefficients, bool, Matrix& jac) const -> int override {
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

    IpoptSystem()
        : OptimizerSystem(4) {
        setNumEqualityConstraints(1);
        setNumInequalityConstraints(1);
    }
};


/* See LBFGSBTest.cpp. 'BestAvailable' will select LBFGSB because this system
 * does NOT have constraints and does have parameter limits.
 */
class LBFGSBSystem : public OptimizerSystem {
    public:
    LBFGSBSystem()
        : OptimizerSystem(25) {
    }

    auto objectiveFunc(const Vector& coefficients, bool, Real& f) const -> int override {
        const Real* x = &coefficients[0];
        f = .25 * (x[0] - 1.0) * (x[0] - 1.0);
        for (int i = 1; i < getNumParameters(); i++) {
            f += pow(x[i] - (x[i - 1] * x[i - 1]), 2.0);
        }
        f *= 4.0;
        return 0;
    }

    auto gradientFunc(const Vector& coefficients, bool, Vector& gradient) const -> int override {
        const Real* x = &coefficients[0];
        Real t1 = x[1] - (x[0] * x[0]);
        gradient[0] = (2.0 * (x[0] - 1.0)) - (16.0 * x[0] * t1);

        for (int i = 1; i < getNumParameters() - 1; i++) {
            const Real t2 = t1;
            t1 = x[i + 1] - (x[i] * x[i]);
            gradient[i] = (8.0 * t2) - (16.0 * x[i] * t1);
        }
        gradient[getNumParameters() - 1] = 8.0 * t1;
        return 0;
    }
};


/* See LBFGSTest.cpp. 'BestAvailable' will select LBFGS because this system
 * does NOT have constraints and does NOT have parameter limits.
 */
class LBFGSSystem : public OptimizerSystem {
    public:
    LBFGSSystem()
        : OptimizerSystem(2) {
    }

    auto objectiveFunc(const Vector& coefficients, bool, Real& f) const -> int override {
        const Real x = coefficients[0];
        const Real y = coefficients[1];
        f = (0.5 * (3 * x * x + 4 * x * y + 6 * y * y)) - (2 * x) + (8 * y);
        return 0;
    }

    auto gradientFunc(const Vector& coefficients, bool, Vector& gradient) const -> int override {
        const Real x = coefficients[0];
        const Real y = coefficients[1];
        gradient[0] = (3 * x) + (2 * y) - 2;
        gradient[1] = (2 * x) + (6 * y) + 8;
        return 0;
    }
};


TEST(OptimizerTest, SelectsInteriorPointForConstrainedSystem) {
    IpoptSystem sys;
    Vector lb(sys.getNumParameters());
    Vector ub(sys.getNumParameters());

    for (int i = 0; i < sys.getNumParameters(); i++) {
        lb[i] = 1.0;
        ub[i] = 5.0;
    }
    sys.setParameterLimits(lb, ub);

    Optimizer opt(sys);
    EXPECT_EQ(opt.getAlgorithm(), SimTK::InteriorPoint);
}

TEST(OptimizerTest, SelectsLBFGSBForBoundedUnconstrainedSystem) {
    LBFGSBSystem sys;
    Vector lb(sys.getNumParameters());
    Vector ub(sys.getNumParameters());

    for (int i = 0; i < sys.getNumParameters(); i += 2) {
        lb[i] = 1.0;
        ub[i] = 100.0;
    }
    for (int i = 1; i < sys.getNumParameters(); i += 2) {
        lb[i] = -100.0;
        ub[i] = 100.0;
    }

    sys.setParameterLimits(lb, ub);

    Optimizer opt(sys);
    EXPECT_EQ(opt.getAlgorithm(), SimTK::LBFGSB);
}

TEST(OptimizerTest, SelectsLBFGSForUnconstrainedUnboundedSystem) {
    LBFGSSystem sys;
    Optimizer opt(sys);
    EXPECT_EQ(opt.getAlgorithm(), SimTK::LBFGS);
}
