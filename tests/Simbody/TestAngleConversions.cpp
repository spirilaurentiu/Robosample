/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKsimbody.h"

using namespace SimTK;

TEST(Simbody_AngleConversions, QuaternionEulerRoundTrip) {
    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    Body::Rigid body = Body::Rigid(MassProperties(1, Vec3(0), Inertia(1)));
    Random::Uniform random(0.0, 2.0);

    MobilizedBody lastBody =
        MobilizedBody::Pin(matter.Ground(),
                           Transform(Vec3(0, 0, 0)),
                           body,
                           Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Slider(lastBody,
                              Transform(Vec3(0, 0, 0)),
                              body,
                              Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Universal(lastBody,
                                 Transform(Vec3(0, 0, 0)),
                                 body,
                                 Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Cylinder(lastBody,
                                Transform(Vec3(0, 0, 0)),
                                body,
                                Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::BendStretch(lastBody,
                                   Transform(Vec3(0, 0, 0)),
                                   body,
                                   Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Planar(lastBody,
                              Transform(Vec3(0, 0, 0)),
                              body,
                              Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Gimbal(lastBody,
                              Transform(Vec3(0, 0, 0)),
                              body,
                              Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Ball(lastBody,
                                   Transform(Vec3(0, 0, 0)),
                                   body,
                                   Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::Translation(lastBody,
                                   Transform(Vec3(0, 0, 0)),
                                   body,
                                   Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Free(lastBody,
                                   Transform(Vec3(0, 0, 0)),
                                   body,
                                   Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::LineOrientation(
        lastBody,
        Transform(Vec3(0, 0, 0)),
        body,
        Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody =
        MobilizedBody::FreeLine(lastBody,
                                Transform(Vec3(0, 0, 0)),
                                body,
                                Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Weld(lastBody,
                                   Transform(Vec3(0, 0, 0)),
                                   body,
                                   Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));
    lastBody = MobilizedBody::Screw(lastBody,
                                    Transform(Vec3(0, 0, 0)),
                                    body,
                                    Transform(Vec3(random.getValue(), random.getValue(), random.getValue())),
                                    0.5);
    lastBody =
        MobilizedBody::Ellipsoid(lastBody,
                                 Transform(Vec3(0, 0, 0)),
                                 body,
                                 Transform(Vec3(random.getValue(), random.getValue(), random.getValue())));

    mbs.realizeTopology();
    State& state = mbs.updDefaultState();
    mbs.realizeModel(state);

    for (int i = 0; i < state.getNQ(); ++i) {
        state.updQ()[i] = random.getValue();
    }
    mbs.realize(state, Stage::Instance);
    mbs.project(state, 0.01);
    mbs.realize(state, Stage::Position);

    // Convert to Euler angles and verify positions are preserved.
    State euler = state;
    matter.convertToEulerAngles(state, euler);
    mbs.realize(euler, Stage::Position);
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& mob = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (mob.getBodyOriginLocation(euler) - mob.getBodyOriginLocation(state)).norm();
        EXPECT_LT(dist, 1e-5) << "Position mismatch after quaternion->Euler conversion at body " << i;
    }

    // Convert back to quaternions and verify positions are still preserved.
    State quaternions = state;
    matter.convertToQuaternions(euler, quaternions);
    mbs.realize(quaternions, Stage::Position);
    for (int i = 0; i < matter.getNumBodies(); ++i) {
        const MobilizedBody& mob = matter.getMobilizedBody(MobilizedBodyIndex(i));
        Real dist = (mob.getBodyOriginLocation(quaternions) - mob.getBodyOriginLocation(state)).norm();
        EXPECT_LT(dist, 1e-5) << "Position mismatch after Euler->quaternion conversion at body " << i;
    }

    // Verify state variables are accurately reproduced after the round trip.
    mbs.project(state, 0.01);
    Real diff = std::sqrt((state.getQ() - quaternions.getQ()).normSqr() / state.getNQ());
    EXPECT_LT(diff, 1e-5) << "Q vector RMS difference after round trip: " << diff;
}
