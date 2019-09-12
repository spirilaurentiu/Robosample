#include "Robo.hpp"
#include "FixmanTorque.hpp"


////////////////////////////
////// GRID FORCE //////////
////////////////////////////
FixmanTorque::FixmanTorque(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter
                    ) : matter(matter){
    this->compoundSystem = compoundSystem;
    scaleFactor = 1.0;
    this->temperature = 0.0;
    this->RT = temperature * SimTK_BOLTZMANN_CONSTANT_MD;
    
}

void FixmanTorque::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
                           SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const  
{
    // Compute Fixman torque
    int nu = state.getNU();
    SimTK::Vector V3(nu);
    SimTK::Vector V4(nu);
    SimTK::Real* D0 = new SimTK::Real(1.0);
    matter.calcFixmanTorque(state, V3, V4, D0);
    delete D0;
    // end - Compute Fixman torque

    // Calculate geometric features fast
    //xstd::cout << " 1 " ;
    //xfor (SimTK::MobilizedBodyIndex mbx(2); mbx < matter.getNumBodies(); ++mbx){
    //x    const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);
    //x    const SimTK::MobilizedBody *p_mobod = &mobod;
    //x    std::cout << std::setprecision(10) << std::fixed <<  ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(state) << ' ' ;
    //x}

    //std::cout << " @ " ;
    int uslot = -1;
    for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter.getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

        for(int k = 0; k < mobod.getNumU(state); k++){
            uslot++;

            //mobod.applyOneMobilityForce(state, k, scaleFactor * V4[uslot], mobilityForces);

            mobod.applyOneMobilityForce(state, k, (-1.0) * RT * V4[uslot], mobilityForces);
            //mobod.applyOneMobilityForce(state, k, (-1.0) * V4[uslot], mobilityForces);

            //xif(int(mbx) == 2){
            //x    std::cout << std::setprecision(10) << std::fixed << (-1.0) * RT * V4[uslot] ;
            //x}
            //    << " to body " << std::setprecision(0) << int(mbx) << " slot " << uslot << " ";

            //std::cout << "Fixman torque scaleFactor " << scaleFactor << " " ;
            //std::cout << " temperature " << temperature << " -RT " << (-1.0) * RT << std::endl ;
        }
    }
    //xstd::cout << std::endl ;

    //const SimTK::Real q = knee.getOneQ(state, 0);
    //const SimTK::Real x = q < low ? q-low : (q > high ? q-high : 0);
    //knee.applyOneMobilityForce(state, 0, -k*x, mobilityForces);
}

FixmanTorque::~FixmanTorque(){}

// This should be carefully analyzed. Intended to be taken from somewhere else.
SimTK::Real FixmanTorque::calcPotentialEnergy(const SimTK::State& state) const {
    SimTK::Real energy = 0.0;
    return energy;
}

bool FixmanTorque::dependsOnlyOnPositions() const {
    return true;
}

SimTK::Real FixmanTorque::getTemperature(void)
{
    return this->temperature;
}

void FixmanTorque::setTemperature(SimTK::Real argTemperature)
{
    this->temperature = argTemperature;
    this->RT = temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

SimTK::Real FixmanTorque::getScaleFactor(void)
{
    return scaleFactor;
}

void FixmanTorque::setScaleFactor(SimTK::Real argScaleFactor)
{
    scaleFactor = argScaleFactor;
    std::cout << "FixmanTorque::setScaleFactor("<< this->scaleFactor << ")" << std::endl;
}



////////////////////////////
////// END GRID FORCE //////
////////////////////////////

