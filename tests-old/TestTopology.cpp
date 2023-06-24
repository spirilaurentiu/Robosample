/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <sys/stat.h>

#include "Robo.hpp"

#include "Topology.hpp"

int main(int argc, char **argv)
{
    /** Constructors name assignment **/
    Topology topology1;
    assert(topology1.getName() == "no_name");

    Topology topology2("topology2");
    assert(topology2.getName() == "topology2");



    return 0;
}


