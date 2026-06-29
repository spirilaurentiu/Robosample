#pragma once
// TEST-ONLY stub of ForceBridge: the kinematics/dynamics sweeps under test never
// invoke the forcefield. The real ForceBridge.hpp pulls OpenMM; this stub lets the
// engine TU compile+link OpenMM-free for the unit tests. verletStep (the only
// caller of evaluate) is not exercised by these tests.
#include "RobotState.hpp"
class ForceBridge {
  public:
    void evaluate(RobotState&) {}
};
