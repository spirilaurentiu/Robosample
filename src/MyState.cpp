#include "MyState.hpp"

// Constructor

MyState::MyState(){

    //;

}

// Destructor

MyState::~MyState(){
    //;
}


// int atom_no
SimTK::Vec3 MyState::getCartesian(int){
    return {};
}

// int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z
void MyState::updCartesian(int, SimTK::Real, SimTK::Real, SimTK::Real){

}

// int atom_no, SimTK::Real x, SimTK::Real y, SimTK::Real z
void MyState::setCartesian(int, SimTK::Real, SimTK::Real, SimTK::Real){

}



