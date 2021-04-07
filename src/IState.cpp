#include "IState.hpp"

IState::~IState() {}

SimTK::Vec3 IState::getCartesian(int)
{
    // args were int atom_no
    return {};
}

void IState::updCartesian(int, SimTK::Real, SimTK::Real, SimTK::Real)
{
    // args were int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z
}

void IState::setCartesian(int, SimTK::Real, SimTK::Real, SimTK::Real)
{
    // args were int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z
}
