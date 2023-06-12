#ifndef __ISTATE_HPP__
#define __ISTATE_HPP__

#include "Robo.hpp"

class IState : public SimTK::State
{
public:
    virtual ~IState()=0;

    // Interface

    virtual SimTK::Vec3 getCartesian(int atom_no) = 0;

    virtual void updCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z) = 0;

    virtual void setCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z) = 0;

};

#endif // __ISTATE_HPP__
