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


    /*            for(unsigned int j = 0; j < otherWorldsAtomsLocations[i].size(); j++){
                SimTK::Compound::AtomIndex aIx = ((otherWorldsAtomsLocations[i][j]).first)->atomIndex;
                SimTK::Vec3 location = ((otherWorldsAtomsLocations[i][j]).second);
                std::cout << "setAtomsLoc atomTargets from previous World i aIx "
                    << j << " " << aIx << " " << location << std::endl;
            }

            for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
                SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
                SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);

                // Check station_B
                std::cout << "setAtomsLoc aIx dumm.station_B " << aIx
                    << " " << forceField->getAtomStationOnBody(dAIx) << std::endl;

            }*/

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
