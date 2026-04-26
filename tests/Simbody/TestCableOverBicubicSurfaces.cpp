/* -------------------------------------------------------------------------- *
 *            Simbody(tm) Adhoc Test: Cable Over Bicubic Surfaces             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012 Stanford University and the Authors.           *
 * Authors: Michael Sherman, Andreas Scholz                                   *
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

/* Simbody CableOverBicubicSurfaces
This example shows how to use a CableTrackerSubsystem to follow the motion of
a cable that crosses bicubic surfaces. We'll then
create a force element that generates spring forces that result from the
stretching and stretching rate of the cable. */

#include <gtest/gtest.h>
#include <iostream>

#include "Simbody.h"

using namespace SimTK;

TEST(Simbody_CableOverBicubicSurfaces, Basic) {
    // Create the system.
    MultibodySystem system;
    SimbodyMatterSubsystem matter(system);
    CableTrackerSubsystem cables(system);
    GeneralForceSubsystem forces(system);

    Force::Gravity gravity(forces, matter, -YAxis, 9.81);
    // Force::GlobalDamper(forces, matter, 5);

    system.setUseUniformBackground(true);   // no ground plane in display
    MobilizedBody Ground = matter.Ground(); // convenient abbreviation

    // Read in some bones.
    PolygonalMesh femur;
    PolygonalMesh tibia;
    femur.loadVtpFile("CableOverBicubicSurfaces-femur.vtp");
    tibia.loadVtpFile("CableOverBicubicSurfaces-tibia.vtp");
    femur.scaleMesh(20);
    tibia.scaleMesh(20);

    // Create some bicubic surfaces.
    constexpr int Nx = 4;
    constexpr int Ny = 5;

    constexpr std::array<Real, Nx> xData = {0.1, 1.0, 2.0, 4.0};
    constexpr std::array<Real, Ny> yData = {-3.0, -2.0, 0.0, 1.0, 3.0};
    constexpr std::array<Real, Nx * Ny> fData = {1.0, 2.0, 3.0, 3.0, 2.0, 1.1, 2.1, 3.1, 3.1, 2.1,
                                                 1.0, 2.0, 7.0, 3.0, 2.0, 1.2, 2.2, 3.2, 3.2, 2.2};

    const Vector x(Nx, xData.data());
    const Vector y(Ny, yData.data());
    const Matrix f(Nx, Ny, fData.data());

    const BicubicSurface rough(x, y, f, 0);  // raw
    const BicubicSurface smooth(x, y, f, 1); // smoothed

    const Vector xp(Vec2(.25, 3.25));
    const Vector yp(Vec2(.75, 5.75));

    // One-hump patch:
    const Matrix fp(Mat22(1, 1, 1, 1));
    const Matrix fxp(Mat22(1, 1, -1, -1));
    const Matrix fyp(Mat22(1, -1, 1, -1));
    const Matrix fxyp(0.5 * Mat22(1, 3, -3, 4));
    const BicubicSurface patch(xp, yp, fp, fxp, fyp, fxyp);
    const Rotation xm90(-Pi / 2, XAxis);
    const Transform patchPose(xm90, Vec3(4, 2, 0));

    // Ask the bicubic surfaces for some meshes we can use for display.
    Real resolution = 31;
    PolygonalMesh patchMesh = patch.createPolygonalMesh(resolution);
    PolygonalMesh roughMesh = rough.createPolygonalMesh(resolution);
    PolygonalMesh smoothMesh = smooth.createPolygonalMesh(resolution);

    const Vec3 SmoothOrigin(-3, -3, -3);

    Body::Rigid someBody(MassProperties(2.0,
                                        Vec3(0, -4, 0),
                                        UnitInertia::cylinderAlongY(1, 4).shiftFromCentroid(Vec3(0, 4, 0))));
    MobilizedBody::Free body1(Ground, Transform(Vec3(0)), someBody, Transform(Vec3(0, 0, 0)));

    CablePath path1(cables,
                    Ground,
                    Vec3(.5, -.5, 0), // origin
                    body1,
                    Vec3(0, 0, 0)); // termination

    CableObstacle::Surface obstacle1(path1, Ground, SmoothOrigin, ContactGeometry::SmoothHeightMap(smooth));

    // Provide an initial guess for P and Q (in frame of "smooth").
    Vec3 P1(1.5, 1, 3.75);
    Vec3 Q1(1.5, -1, 3.75);
    obstacle1.setContactPointHints(P1, Q1);

    CableSpring cable1(forces, path1, 50., 8., 0.1);

    // Initialize the system and state.
    EXPECT_NO_THROW(system.realizeTopology());
    State state = system.getDefaultState();

    // Random::Gaussian random;
    // for (int i = 0; i < state.getNQ(); ++i)
    //     state.updQ()[i] = random.getValue();
    // for (int i = 0; i < state.getNU(); ++i)
    //     state.updU()[i] = 0.1*random.getValue();

    body1.setQToFitTranslation(state, Vec3(4, -10, -3));
    body1.setQToFitRotation(state, Rotation(-Pi, ZAxis));

    EXPECT_NO_THROW(system.realize(state, Stage::Position));

    path1.setIntegratedCableLengthDot(state, path1.getCableLength(state));

    RungeKuttaMersonIntegrator integ(system);
    integ.setAccuracy(1e-3);

    TimeStepper timeStepper(system, integ);
    EXPECT_NO_THROW(timeStepper.initialize(state));

    const Real finalTime = 10;
    const double startTime = realTime();
    const double startCPU = cpuTime();
    EXPECT_NO_THROW(timeStepper.stepTo(finalTime));
}
