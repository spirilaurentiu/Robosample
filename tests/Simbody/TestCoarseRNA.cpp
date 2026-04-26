/* -------------------------------------------------------------------------- *
 *                           SimTK Simbody(tm)                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

#include <cmath>
#include <cstdio>
#include <gtest/gtest.h>
#include <vector>

#include "SimTKsimbody.h"

using namespace SimTK;

static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
static const int GroundBodyNum = 0; // ground is always body 0

static const Real g = 9.8; // meters/s^2; apply in �y direction

static const Real DuplexRadius = 3; // A
static const Real HalfHeight = 10;  // A
static const Real CylinderSlop = 1; // A

static const int NAtoms = 20;
static const Real AtomMass = 12;  // Daltons
static const Real AtomRadius = 1; // A

static const Real ConnectorRadius = 1;     // A
static const Real ConnectorHalfHeight = 3; // A
static const Real ConnectorEndSlop = 0.2;  // A
static const Real ConnectorDensity = 10;   // Dalton/A^3

static int NSegments = 3;
bool shouldFlop = false;

class MyRNAExample : public SimbodyMatterSubsystem {
    struct PerBodyInfo {
        PerBodyInfo(MobilizedBodyIndex b, bool d)
            : bnum(b)
            , isDuplex(d) {
        }
        MobilizedBodyIndex bnum;
        bool isDuplex;
    };
    std::vector<PerBodyInfo> bodyInfo;
    MobilizedBodyIndex end1, end2;

    public:
    MyRNAExample(MultibodySystem& mbs, int nsegs, bool shouldFlop)
        : SimbodyMatterSubsystem(mbs) {
        bodyInfo.emplace_back(GroundIndex, false); // placeholder for ground
        end1 = makeChain(GroundIndex, Vec3(0), nsegs, shouldFlop);
        end2 = makeChain(GroundIndex, Vec3(20, 0, 0), nsegs, shouldFlop);

        Constraint::Rod theConstraint2(updMobilizedBody(end1),
                                       Vec3(0, -HalfHeight, 0),
                                       updMobilizedBody(end2),
                                       Vec3(0, -HalfHeight, 0),
                                       10);
    }

    private:
    auto makeChain(MobilizedBodyIndex startBodyId, const Vec3& startOrigin, int nSegs, bool shouldFlop)
        -> MobilizedBodyIndex {
        MobilizedBody baseBody = updMobilizedBody(startBodyId);
        Vec3 origin = startOrigin;
        MobilizedBody lastDup;
        for (int seg = 0; seg < nSegs; ++seg) {
            MobilizedBody::Ball left1(
                baseBody,
                Transform(origin + Vec3(-DuplexRadius, -HalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            left1.setDefaultRadius(1.5);
            bodyInfo.emplace_back(left1, false);

            MobilizedBody::Ball left2(
                left1,
                Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            left2.setDefaultRadius(1.5);
            bodyInfo.emplace_back(left2, false);

            MobilizedBody::Ball rt1(
                baseBody,
                Transform(origin + Vec3(DuplexRadius, -HalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            rt1.setDefaultRadius(1.5);
            bodyInfo.emplace_back(rt1, false);

            MobilizedBody::Ball rt2(
                rt1,
                Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            rt2.setDefaultRadius(1.5);
            bodyInfo.emplace_back(rt2, false);

            MobilizedBody::Ball dup(
                rt2,
                Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcDuplexMassProps(DuplexRadius, HalfHeight, NAtoms, AtomMass)),
                Transform(Vec3(-DuplexRadius, HalfHeight, 0)));
            dup.setDefaultRadius(1.5);
            bodyInfo.emplace_back(dup, true);

            if (!shouldFlop) {
                Constraint::Ball theConstraint(left2,
                                               Vec3(0, -ConnectorHalfHeight, 0),
                                               dup,
                                               Vec3(DuplexRadius, HalfHeight, 0));
                theConstraint.setDefaultRadius(1.5);
            }

            baseBody = lastDup = dup;
            origin = Vec3(0);
        }
        return lastDup;
    }

    static auto calcDuplexMassProps(Real halfHeight, Real r, int nAtoms, Real atomMass) -> MassProperties {
        const Real pitch = 2 * Pi / halfHeight;
        const Real trans = (2 * halfHeight) / (nAtoms - 1);
        const Real rot = pitch * trans;
        Inertia iner(0);
        Vec3 com(0);
        Real mass = 0;
        for (int i = 0; i < nAtoms; ++i) {
            const Real h = halfHeight - (i * trans);
            const Real th = i * rot;
            const Vec3 p1(-r * cos(th), h, r * sin(th));
            const Vec3 p2(r * cos(th), h, -r * sin(th));
            mass += 2 * atomMass;
            iner += Inertia(p1, atomMass) + Inertia(p2, atomMass);
            com += atomMass * p1 + atomMass * p2;
        }
        return {mass, com / mass, iner};
    }

    static auto calcConnectorMassProps(Real r, Real halfHeight, Real density) -> MassProperties {
        const Real volume = Pi * r * r * halfHeight;
        const Real mass = volume * density;
        const Vec3 com = Vec3(0);
        const Inertia iner = mass * UnitInertia::cylinderAlongY(r, halfHeight);

        return {mass, com, iner};
    }
};

TEST(Simbody_CoarseRNA, Basic) {
    // Create a multibody system using Simbody.
    MultibodySystem mbs;
    MyRNAExample myRNA(mbs, NSegments, shouldFlop);
    GeneralForceSubsystem forces(mbs);
    Force::UniformGravity ugs(forces, myRNA, Vec3(0, -g, 0), -0.8);

    const Vec3 attachPt(150, -40, -50);
    Force::TwoPointLinearSpring(forces,
                                myRNA.Ground(),
                                attachPt,
                                myRNA.getMobilizedBody(MobilizedBodyIndex(myRNA.getNumBodies() - 1)),
                                Vec3(0),
                                1000., // stiffness
                                1.);   // natural length
    Force::GlobalDamper(forces, myRNA, 1000);

    State state = mbs.realizeTopology();

    mbs.realizeModel(state);
    mbs.realize(state, Stage::Position);

    EXPECT_EQ(myRNA.getNumBodies(), 31)
        << "Expected Ground + 2 chains x 3 segments x 5 bodies per segment = 31 total bodies. "
           "Check makeChain() or NSegments if this fails.";
    EXPECT_EQ(myRNA.getNumConstraints(), 7)
        << "Expected 6 Ball constraints (2 chains x 3 segments, shouldFlop=false) "
           "+ 1 Rod constraint between chain endpoints = 7 total. "
           "If shouldFlop is true, expect 1 (only the Rod).";
    EXPECT_EQ(myRNA.getNumQuaternionsInUse(state), 30)
        << "Each Ball mobilizer requires a quaternion for orientation. "
           "2 chains x 3 segments x 5 Ball joints = 30 quaternions expected.";
    EXPECT_LT(myRNA.getQErr(state).normRMS(), 1e-10)
        << "Position constraint error at the initial assembled configuration should be "
           "machine-precision zero before any integration step is taken.";

    for (ConstraintIndex cid(0); cid < myRNA.getNumConstraints(); ++cid) {
        const Constraint& constraint = myRNA.getConstraint(cid);
        int mp{};
        int mv{};
        int ma{};
        constraint.getNumConstraintEquationsInUse(state, mp, mv, ma);

        std::cout << "CONSTRAINT " << cid << " constrained bodies=" << constraint.getNumConstrainedBodies()
                  << " ancestor=" << constraint.getAncestorMobilizedBody().getMobilizedBodyIndex()
                  << " constrained mobilizers/nq/nu=" << constraint.getNumConstrainedMobilizers() << "/"
                  << constraint.getNumConstrainedQ(state) << "/" << constraint.getNumConstrainedU(state)
                  << " mp,mv,ma=" << mp << "," << mv << "," << ma << std::endl;

        for (ConstrainedBodyIndex cid(0); cid < constraint.getNumConstrainedBodies(); ++cid) {
            std::cout << "  constrained body: "
                      << constraint.getMobilizedBodyFromConstrainedBody(cid).getMobilizedBodyIndex();
            std::cout << std::endl;
        }

        for (ConstrainedMobilizerIndex cmx(0); cmx < constraint.getNumConstrainedMobilizers(); ++cmx) {
            std::cout << "  constrained mobilizer "
                      << constraint.getMobilizedBodyFromConstrainedMobilizer(cmx).getMobilizedBodyIndex()
                      << ", q(" << constraint.getNumConstrainedQ(state, cmx) << ")=";

            for (MobilizerQIndex i(0); i < constraint.getNumConstrainedQ(state, cmx); ++i) {
                std::cout << " " << constraint.getConstrainedQIndex(state, cmx, i);
            }

            std::cout << ", u(" << constraint.getNumConstrainedU(state, cmx) << ")=";

            for (MobilizerUIndex i(0); i < constraint.getNumConstrainedU(state, cmx); ++i) {
                std::cout << " " << constraint.getConstrainedUIndex(state, cmx, i);
            }

            std::cout << std::endl;
        }

        std::cout << constraint.getSubtree();

        std::cout << "   d(perrdot)/du=" << constraint.calcPositionConstraintMatrixP(state);
        std::cout << "   d(perrdot)/du=" << ~constraint.calcPositionConstraintMatrixPt(state);

        std::cout << "   d(perr)/dq=" << constraint.calcPositionConstraintMatrixPNInv(state);
    }


    SimbodyMatterSubtree sub(myRNA);
    sub.addTerminalBody(myRNA.getMobilizedBody(MobilizedBodyIndex(7)));
    sub.addTerminalBody(myRNA.getMobilizedBody(MobilizedBodyIndex(10)));
    EXPECT_NO_THROW(sub.realizeTopology());

    EXPECT_EQ((int)sub.getAllBodies().size(), 6)
        << "Subtree from common ancestor (body 5) to terminals 7 and 10 should span exactly 6 bodies: "
           "{5(ancestor), 6(left1), 7(left2/terminal), 8(rt1), 9(rt2), 10(dup/terminal)}.";
    EXPECT_EQ((int)sub.getAncestorMobilizedBodyIndex(), 5)
        << "Lowest common ancestor of bodies 7 and 10 should be body 5 "
           "(dup of segment 0, chain 1). Both paths converge at this body.";
    EXPECT_EQ((int)sub.getParentSubtreeBodyIndex(SubtreeBodyIndex(0)), -1)
        << "SubtreeBodyIndex 0 is always the ancestor; its parent index must be -1 (no parent).";

    std::cout << "sub.ancestor=" << sub.getAncestorMobilizedBodyIndex();
    //    std::cout << "  sub.terminalBodies=" << sub.getTerminalBodies() << std::endl;
    //    std::cout << "sub.allBodies=" << sub.getAllBodies() << std::endl;
    for (SubtreeBodyIndex i(0); i < (int)sub.getAllBodies().size(); ++i) {
        std::cout << "sub.parent[" << i << "]=" << sub.getParentSubtreeBodyIndex(i);
        //       std::cout << "  sub.children[" << i << "]=" << sub.getChildSubtreeBodyIndexs(i) << std::endl;
    }

    printf("# quaternions in use = %d\n", myRNA.getNumQuaternionsInUse(state));
    for (MobilizedBodyIndex i(0); i < myRNA.getNumBodies(); ++i) {
        printf("body %2d: using quat? %s; quat index=%d\n",
               (int)i,
               myRNA.isUsingQuaternion(state, i) ? "true" : "false",
               (int)myRNA.getQuaternionPoolIndex(state, i));
    }

    // And a study using the Runge Kutta Merson integrator
    bool suppressProject = false;

    RungeKuttaMersonIntegrator myStudy(mbs);
    // RungeKuttaFeldbergIntegrator myStudy(mbs);
    // RungeKutta3Integrator myStudy(mbs);
    // CPodesIntegrator  myStudy(mbs);
    // VerletIntegrator myStudy(mbs);
    // ExplicitEulerIntegrator myStudy(mbs);

    myStudy.setAccuracy(1e-2);
    myStudy.setConstraintTolerance(1e-3);
    myStudy.setProjectEveryStep(false);

    const Real dt = 1. / 30; // output intervals

    state.updTime() = 0;
    myStudy.initialize(state);
    std::cout << "Using Integrator " << std::string(myStudy.getMethodName()) << ":\n";
    std::cout << "ACCURACY IN USE=" << myStudy.getAccuracyInUse() << std::endl;
    std::cout << "CTOL IN USE=" << myStudy.getConstraintToleranceInUse() << std::endl;
    std::cout << "TIMESCALE=" << mbs.getDefaultTimeScale() << std::endl;
    std::cout << "U WEIGHTS=" << state.getUWeights() << std::endl;
    std::cout << "Z WEIGHTS=" << state.getZWeights() << std::endl;
    std::cout << "1/QTOLS=" << state.getQErrWeights() << std::endl;
    std::cout << "1/UTOLS=" << state.getUErrWeights() << std::endl;

    EXPECT_NEAR(myStudy.getAccuracyInUse(), 1e-2, 1e-12)
        << "Integrator accuracy in use must match the requested value of 1e-2. "
           "A mismatch indicates the integrator silently clamped or overrode the setting.";
    EXPECT_NEAR(myStudy.getConstraintToleranceInUse(), 1e-3, 1e-12)
        << "Constraint tolerance in use must match the requested value of 1e-3. "
           "Used to validate QErr/UErr acceptance thresholds during projection.";
    EXPECT_EQ(myStudy.getState().getTime(), 0.0)
        << "Integrator state time must be exactly 0 immediately after initialize(). "
           "Non-zero would indicate state was not reset properly before integration.";

    const double startReal = realTime();
    const double startCPU = cpuTime();
    int stepNum = 0;

    for (;;) {
        const State& ss = myStudy.getState();

        mbs.realize(ss);

        if ((stepNum++ % 100) == 0) {
            printf("%5g qerr=%10.4g uerr=%10.4g hNext=%g\n",
                   ss.getTime(),
                   myRNA.getQErr(ss).normRMS(),
                   myRNA.getUErr(ss).normRMS(),
                   myStudy.getPredictedNextStepSize());
            printf("      E=%14.8g (pe=%10.4g ke=%10.4g)\n",
                   mbs.calcEnergy(ss),
                   mbs.calcPotentialEnergy(ss),
                   mbs.calcKineticEnergy(ss));

            std::cout << "QERR=" << ss.getQErr() << std::endl;
            std::cout << "UERR=" << ss.getUErr() << std::endl;
        }

        if (ss.getTime() >= 10) {
            break;
        }

        // TODO: should check for errors or have or teach RKM to throw.
        myStudy.stepTo(ss.getTime() + dt, Infinity);
    }

    const State& finalState = myStudy.getState();
    mbs.realize(finalState, Stage::Dynamics);

    printf("CPU time=%gs, REAL time=%gs\n", cpuTime() - startCPU, realTime() - startReal);
    printf("Using Integrator %s:\n", myStudy.getMethodName());
    printf("# STEPS/ATTEMPTS = %d/%d\n", myStudy.getNumStepsTaken(), myStudy.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", myStudy.getNumErrorTestFailures());
    printf("# CONVERGENCE FAILS = %d\n", myStudy.getNumConvergenceTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", myStudy.getNumRealizations(), myStudy.getNumProjections());
    printf("# PROJECTION FAILS = %d\n", myStudy.getNumProjectionFailures());
}
