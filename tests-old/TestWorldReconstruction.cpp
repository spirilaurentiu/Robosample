/* -------------------------------------------------------------------------- *
 *                               TestSampler                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include <iostream>
#include "Robo.hpp"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#define ASSERT_EQUAL(val1, val2) {ASSERT(std::abs(val1-val2) < 1e-10);}

void testsReconstruct(){

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

    // Empty Compound
    SimTK::Compound compound;



}

int main() {
    try {
        testsReconstruct();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
