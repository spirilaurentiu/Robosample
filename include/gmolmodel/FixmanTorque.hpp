#ifndef __FIXMANTORQUE_HPP__
#define __FIXMANTORQUE_HPP__

/* -------------------------------------------------------------------------- *
 *                               Robosampling                                 *
 * -------------------------------------------------------------------------- *
 * This is part of Robosampling                                               *
 */

//==============================================================================
//                           CLASS FixmanTorque
//==============================================================================
/**
 **/

class FixmanTorque : public SimTK::Force::Custom::Implementation {
public:
    SimTK::CompoundSystem *compoundSystem;
    int *flag;

    FixmanTorque(SimTK::CompoundSystem *compoundSystem, SimTK::SimbodyMatterSubsystem& matter);
    ~FixmanTorque();

    void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
        SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const override;

    SimTK::Real calcPotentialEnergy(const SimTK::State& state) const override;

    bool dependsOnlyOnPositions() const;

    SimTK::Real getScaleFactor(void);
    void setScaleFactor(SimTK::Real);

    SimTK::Real getTemperature(void);
    void setTemperature(SimTK::Real);

private:
    SimTK::SimbodyMatterSubsystem& matter;
    SimTK::Real temperature;
    SimTK::Real RT;
    SimTK::Real scaleFactor;

};

#endif //__FIXMANTORQUE_HPP__
