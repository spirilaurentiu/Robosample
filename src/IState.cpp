#include "IState.hpp"

IState::~IState() {}

SimTK::Vec3 IState::getCartesian(int atom_no){}

void IState::updCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z){}

void IState::setCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z){}


