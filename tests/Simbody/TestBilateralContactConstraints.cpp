/* -------------------------------------------------------------------------- *
 *           Simbody(tm) Adhoc Test: Bilateral Contact Constraints            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

/* This uses a variety of bilateral constraints that are intended as the
underpinnings for unilateral constraints. The most important check here is
that energy should be conserved perfectly (to integration accuracy) since
these are all non-working constraints. (Look at the total energy in the
visualizer.) Even the friction constraints are
non-working because the underlying constraints represent rolling (a.k.a.
"stiction"); sliding is imposed elsewhere by disabling the rolling constraints
and replacing them with different conditions. */

#include <gtest/gtest.h>

#include "Simbody.h"

using namespace SimTK;

TEST(Simbody_BilateralContactConstraints, Basic) {
    // Define the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::Gravity gravity(forces, matter, -YAxis, 9.8 / 10);

    // Describe mass and visualization properties for a generic body.
    Real mass = 2;
    Vec3 hdim(1, .5, .25);
    Body::Rigid bodyInfo(MassProperties(mass, Vec3(0), UnitInertia::brick(hdim)));
    bodyInfo.addDecoration(Transform(), DecorativeBrick(hdim).setColor(Orange).setOpacity(.3));

    Real pmass = .1;
    Vec3 phdim(5, .5, 2);
    Body::Rigid platformBody(MassProperties(10 * mass, Vec3(0), UnitInertia::ellipsoid(phdim)));
    platformBody.addDecoration(Transform(),
                               DecorativeEllipsoid(phdim).setColor(Cyan).setOpacity(.1).setResolution(5));

    MobilizedBody::Ball platform(matter.Ground(), Vec3(0), platformBody, phdim / 2);

    // Create the moving (mobilized) bodies of the pendulum.
    MobilizedBody::Free brick(matter.Ground(), Transform(Vec3(0)), bodyInfo, Transform(Vec3(0)));

    const Rotation ZtoY(-Pi / 2, XAxis);

    Constraint::SphereOnSphereContact
        sphereOnSphere(brick, hdim, 0.5, platform, Vec3(-3, 1, -.5), 1.2, false);

    Constraint::Rod rod1(brick, Vec3(0, hdim[1], hdim[2]), platform, Vec3(0, 3, -.5), 1.5 * 1.2);

    // Try edge/edge contact.
    Constraint::LineOnLineContact lineOneLine(
        platform,
        Transform(Rotation(UnitVec3(1, 1, 1), XAxis, UnitVec3(-XAxis), ZAxis), Vec3(1, 1, 1)),
        2, // hlen
        brick,
        Transform(Rotation(UnitVec3(ZAxis), XAxis, Vec3(-1, -1, 0), ZAxis), Vec3(-hdim[0], -hdim[1], 0)),
        2, // hlen
        true);

    // Initialize the system and acquire default state.
    State state = system.realizeTopology();

    // -----------------------------------------------------------------------
    // Topology checks: confirm the expected number of bodies and constraints
    // were registered. Ground counts as body 0, so total = Ground + platform
    // + brick = 3.  Constraints are ss, rod1, and ll.
    // -----------------------------------------------------------------------
    EXPECT_EQ(matter.getNumBodies(), 3) << "Expected Ground + platform + brick";
    EXPECT_EQ(matter.getNumConstraints(), 3) << "Expected SphereOnSphere + Rod + LineOnLine";

    // The state must have reached at least the Topology stage so that
    // subsequent realize() calls have a valid starting point.
    EXPECT_GE(state.getSystemStage(), Stage::Topology)
        << "System stage should be at least Topology after realizeTopology()";

    brick.setQToFitTransform(state, Vec3(0, 5, 0));
    brick.setUToFitAngularVelocity(state, Vec3(10, 10, 10));

    // -----------------------------------------------------------------------
    // Initial-condition checks: verify the prescribed position and angular
    // velocity were actually applied to the brick body.
    // -----------------------------------------------------------------------
    system.realize(state, Stage::Velocity);

    // The brick was placed at Y = 5; allow a modest tolerance because
    // setQToFitTransform works in body coordinates and may introduce small
    // numeric adjustments to satisfy joint constraints.
    Vec3 brickPos = brick.getBodyOriginLocation(state);
    EXPECT_NEAR(brickPos[1], 5.0, 0.5)
        << "Brick Y position should be approximately 5 after setQToFitTransform";

    // Angular velocity components should match the prescribed (10,10,10) rad/s.
    Vec3 brickAngVel = brick.getBodyAngularVelocity(state);
    EXPECT_NEAR(brickAngVel[0], 10.0, 1e-10) << "Brick angular velocity X component";
    EXPECT_NEAR(brickAngVel[1], 10.0, 1e-10) << "Brick angular velocity Y component";
    EXPECT_NEAR(brickAngVel[2], 10.0, 1e-10) << "Brick angular velocity Z component";

    Assembler asmb(system);
    asmb.assemble(state);

    // -----------------------------------------------------------------------
    // Post-assembly constraint-satisfaction checks: all position-level and
    // velocity-level constraint errors must be near zero after the Assembler
    // has run.  The Assembler's own convergence tolerance is tight, so 1e-6
    // is a reasonable upper bound on each scalar error.
    // -----------------------------------------------------------------------
    system.realize(state, Stage::Position);
    const Vector& qerr = state.getQErr();
    for (int i = 0; i < qerr.size(); ++i) {
        EXPECT_NEAR(qerr[i], 0.0, 1e-6)
            << "Position-level constraint error at index " << i << " should be near zero after assembly";
    }

    // Note: the Assembler only satisfies position-level constraints (qerr).
    // Velocity-level constraints (uerr) are projected by the integrator during
    // ts.initialize(), so uerr is checked there, not here.

    // Choose integrator and simulate for 10 seconds.
    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-8);

    TimeStepper timeStepper(system, integ);
    timeStepper.initialize(state);

    // -----------------------------------------------------------------------
    // Baseline energy capture: record total mechanical energy (KE + PE) at
    // the start of the simulation.  All constraints in this test are
    // non-working (rolling / stiction), so energy must be an invariant of
    // the motion.  We also verify the initial energy is a finite number.
    //
    // Velocity-level constraint errors are also checked here: ts.initialize()
    // projects velocities onto the constraint manifold, so uerr must be near
    // zero before the first step.
    // -----------------------------------------------------------------------
    {
        const State& initialState = integ.getState();

        system.realize(initialState, Stage::Velocity);
        const Vector& uerr0 = initialState.getUErr();
        for (int i = 0; i < uerr0.size(); ++i) {
            EXPECT_NEAR(uerr0[i], 0.0, 1e-6)
                << "Velocity-level constraint error at index " << i
                << " should be near zero after ts.initialize() projects velocities";
        }
        system.realize(initialState, Stage::Dynamics);
        const Real initialEnergy = system.calcEnergy(initialState);

        EXPECT_FALSE(std::isnan(initialEnergy)) << "Initial total energy must not be NaN";
        EXPECT_FALSE(std::isinf(initialEnergy)) << "Initial total energy must not be infinite";

        EXPECT_NO_THROW(timeStepper.stepTo(100.0));

        // -----------------------------------------------------------------------
        // Post-simulation checks
        // -----------------------------------------------------------------------
        const State& finalState = integ.getState();

        // --- Energy conservation ---
        // Because every constraint is non-working and the only external force
        // (gravity) is conservative, total mechanical energy must be preserved
        // to the integrator's accuracy.  With RKMerson at 1e-8 over 100 s, a
        // relative drift of 1e-3 is a generous but meaningful bound.
        system.realize(finalState, Stage::Dynamics);
        const Real finalEnergy = system.calcEnergy(finalState);

        EXPECT_FALSE(std::isnan(finalEnergy)) << "Final total energy must not be NaN";
        EXPECT_FALSE(std::isinf(finalEnergy)) << "Final total energy must not be infinite";
        EXPECT_NEAR(finalEnergy, initialEnergy, std::abs(initialEnergy) * 1e-3)
            << "Total mechanical energy should be conserved to within 0.1% "
               "over the 100-second simulation";

        // --- Constraint satisfaction at end of simulation ---
        // The integrator enforces constraints at each step; check that
        // accumulated drift has not pushed errors outside a coarse tolerance.
        system.realize(finalState, Stage::Position);
        const Vector& finalQerr = finalState.getQErr();
        for (int i = 0; i < finalQerr.size(); ++i) {
            EXPECT_NEAR(finalQerr[i], 0.0, 1e-4) << "Final position-level constraint error at index " << i;
        }

        system.realize(finalState, Stage::Velocity);
        const Vector& finalUerr = finalState.getUErr();
        for (int i = 0; i < finalUerr.size(); ++i) {
            EXPECT_NEAR(finalUerr[i], 0.0, 1e-4) << "Final velocity-level constraint error at index " << i;
        }
    }
}
