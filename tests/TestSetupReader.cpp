/* -------------------------------------------------------------------------- *
 *                               TestSetupReader                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include <iostream>
#include "SetupReader.hpp"

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

#define ASSERT_EQUAL(val1, val2) {ASSERT(std::abs(val1-val2) < 1e-10);}

/**  **/
class DerivedSetupReader : virtual public SetupReader{
public:
    // Constructor
    DerivedSetupReader (std::string& FN) : SetupReader(FN){}

    ~DerivedSetupReader(){}
};

/**  **/
void testSetupReader(){

    // Setup Reader
    DerivedSetupReader setupReader();
}

/**  **/
int main() {
    try {
        testSetupReader();
    }
    catch(const std::exception& e) {
        std::cout << "exception: " << e.what() << std::endl;
        return 1;
    }
    std::cout << "Done" << std::endl;
    return 0;
}
