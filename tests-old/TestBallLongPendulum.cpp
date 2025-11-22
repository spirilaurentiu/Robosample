/* -------------------------------------------------------------------------- *
 *                      Simbody(tm) Example: Long Pendulum                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
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

/*                      Simbody ExampleLongPendulum
This example shows how to build a linked chain of bodies programmatically,
simulate it, and produce a simple animation while it is simulating. */

#include "Simbody.h"
#include <iostream>

using namespace SimTK;

int main() {
  try {    
    // Create the system.
    MultibodySystem system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem matter(system);
    GeneralForceSubsystem forces(system);
    Force::UniformGravity noGravity(forces, matter, Vec3(0, 0, 0));
    Random::Gaussian random;
    random.setSeed(time(NULL));

    // First body
    Real mass = 0;
    Inertia inertia(0);
    Vec3 com(0);
/*    Real mass = 17.007975980000001;
    com = Vec3(0, 0.059264899079425913, 0);
    inertia += Inertia(Vec3(0,0,0), 16);
    inertia += Inertia(Vec3(0,1,0), 1);*/
/*     Real mass = 32;
     com = Vec3(0, 0.5, 0);
     inertia += Inertia(Vec3(0,0,0), 16);
     inertia += Inertia(Vec3(0,1,0), 16);
    //std::cout << "Inertia tensor " << inertia << std::endl;

    //inertia = Inertia(1);

    Body::Rigid freeBody(MassProperties(mass, com, inertia));
    MobilizedBody::Free first(matter.Ground(),
                                Transform(Vec3(0, 0, 0)),      // X_PF
                                freeBody,            // Mass properties
                                Transform(Vec3(0, 0, 0)) // X_BM
    );
    //freeBody.addDecoration(Transform(), DecorativeSphere(0.1));
    freeBody.addDecoration(Transform(), DecorativeFrame());*/

    // Second body
    //inertia = Inertia(0);
    //mass = 17.007975980000001;
    //com = Vec3(0, 0.059264899079425913, 0);
    //inertia += Inertia(Vec3(0,0,0), 16);
    //inertia += Inertia(Vec3(0,1,0), 1);

    mass = 16 ; //+ 16; //+ 16; //+ 16;

    com += 16 * Vec3(0,0,0);
    //com += 16 * Vec3(1,0,0);
    //com += 16 * Vec3(0,1,0);
    //com += 16 * Vec3(0,0,1);
    com /= mass;

    inertia += Inertia(Vec3(0, 0, 0), 16);
    //inertia += Inertia(Vec3(1, 0, 0), 16);
    //inertia += Inertia(Vec3(0, 1, 0), 16);
    //inertia += Inertia(Vec3(0, 0, 1), 16);

    std::cout << "Ball body initial mass " << mass << std::endl;
    std::cout << "Ball body initial com " << com << std::endl;
    std::cout << "Ball body initial inertia " << inertia << std::endl;

    Transform X_PF_P(Vec3(random.getValue(), random.getValue(), random.getValue()));
    //Transform P_X_M(Rotation(35*SimTK_DEGREE_TO_RADIAN, XAxis), Vec3(1, 0, 0));
    //Transform M_X_pin(Rotation(-90*SimTK_DEGREE_TO_RADIAN, YAxis));
    //std::cout << "Ball body initial X_PF " << P_X_M << std::endl;
    Body::Rigid ballBody(MassProperties(mass, com, inertia));
    MobilizedBody::Ball second(matter.Ground(),
            //P_X_M * M_X_pin,
            Transform(Vec3(0, 1, 0)),      // X_PF
            ballBody,            // Mass properties
            Transform(Vec3(0, 0, 0)) // X_BM
    );
    ballBody.addDecoration(Transform(), DecorativeFrame());
    ballBody.addDecoration(Transform(Vec3(0, 0, 0)), DecorativeBrick());
    //ballBody.addDecoration(Transform(Vec3(1, 0, 0)), DecorativeBrick());
    //ballBody.addDecoration(Transform(Vec3(0, 1, 0)), DecorativeBrick());
    //ballBody.addDecoration(Transform(Vec3(0, 0, 1)), DecorativeBrick());

      // Visualizer
    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 1./30));
    
    // Initialize the system and state.
    system.realizeTopology();
    State state = system.getDefaultState();
/*    for (int i = 0; i < state.getNQ(); ++i)
        state.updQ()[i] = random.getValue();*/

      //for (int i = 0; i < state.getNU(); ++i)
       //   state.updU()[i] = random.getValue();

    // Simulate it.
    VerletIntegrator integ(system, 0.2);
    //RungeKuttaMersonIntegrator integ(system);
    TimeStepper ts(system, integ);
    ts.initialize(state);

    // Run
    int nofRounds = 10;
    State advState = integ.getAdvancedState();
    for(int i = 0; i < nofRounds; i++){
        ts.initialize(advState);
        for (int ui = 0; ui < advState.getNU(); ++ui){
            advState.updU()[ui] = random.getValue() / 0.5;
        }
        system.realize(advState, Stage::Velocity);
        ts.stepTo(i * 3.0);
        //ts.stepTo(0.000000000001);
    }
    std::cout << "Finished" << std::endl;

  } catch(const std::exception& e) {
    std::cout << "EXCEPTION: " << e.what() << std::endl;
    return 1;
  } catch (...) {
      std::cout << "UNKNOWN EXCEPTION\n";
      return 1;
  }
    return 0;
}
