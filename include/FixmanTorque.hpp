#pragma once

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
    FixmanTorque(SimTK::SimbodyMatterSubsystem* argMatter);
    ~FixmanTorque();

    void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
        SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const override;

    SimTK::Real calcPotentialEnergy(const SimTK::State& state) const override;

    bool dependsOnlyOnPositions() const override;

    SimTK::Real getScaleFactor(void);
    void setScaleFactor(SimTK::Real);

    SimTK::Real getTemperature(void);
    void setTemperature(SimTK::Real);

private:
    // std::shared_ptr<SimTK::SimbodyMatterSubsystem> matter;
    SimTK::SimbodyMatterSubsystem* matter = nullptr;
    SimTK::Real temperature;
    SimTK::Real RT;
    SimTK::Real scaleFactor;
};

//==============================================================================
//                           CLASS FixmanTorqueExt
//==============================================================================
/**
 **/

class FixmanTorqueExt : public SimTK::Force::Custom::Implementation {
public:
    FixmanTorqueExt(SimTK::SimbodyMatterSubsystem* argMatter);
    ~FixmanTorqueExt();

    void calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces,
        SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces) const override;

    SimTK::Real calcPotentialEnergy(const SimTK::State& state) const override;

    bool dependsOnlyOnPositions() const override;

    SimTK::Real getScaleFactor(void);
    void setScaleFactor(SimTK::Real);

    SimTK::Real getTemperature(void);
    void setTemperature(SimTK::Real);

private:
    // std::shared_ptr<SimTK::SimbodyMatterSubsystem> matter;
    SimTK::SimbodyMatterSubsystem* matter = nullptr;
    SimTK::Real temperature;
    SimTK::Real RT;
    SimTK::Real scaleFactor;
};
