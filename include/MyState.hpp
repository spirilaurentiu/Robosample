#ifndef __MYSTATE_HPP__
#define __MYSTATE_HPP__

#include "Robo.hpp"
#include "IState.hpp"

class MyState : public IState
{
public:

    // Constructor

    MyState();

    // Destructor

    ~MyState();

    // Interface

    SimTK::Vec3 getCartesian(int atom_no);

    void updCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z);

    void setCartesian(int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z);



};

#endif // __IMYSTATE_HPP__
