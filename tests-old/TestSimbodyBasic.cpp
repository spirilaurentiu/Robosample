/* -------------------------------------------------------------------------- *
 *                 SimTK Molmodel Example: Simple Protein                     *
 * -------------------------------------------------------------------------- *
 * This is the first example from the Molmodel User's Guide. It creates a     *
 * small protein (a five-residue peptide), simulates it and generates a live  *
 * animation while it is running.                                             *
 *                                                                            *
 * Authors: Christopher Bruns, Michael Sherman                                *
 * -------------------------------------------------------------------------- */

#include "Molmodel.h"
#include <iostream>
#include <exception>
using namespace SimTK;

int main() {
try {
    
    // SIMBODY CONSTRUCTORS WORK:
    System someSystem;
    MultibodySystem someMultibodySystem; // Extends System (Simbody) // adoptSystemGuts(new MultibodySystemRep())
    SimbodyMatterSubsystem matter(someMultibodySystem); // Extends Subsystem (Simbody)
    DecorationSubsystem decorations(someMultibodySystem); // Extends Subsystem (Simbody)
    // =====================

    return 0;
} 
catch(const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    return 1;
}
catch(...) {
    std::cerr << "ERROR: An unknown exception was raised" 
              << std::endl;
    return 1;
}

}


