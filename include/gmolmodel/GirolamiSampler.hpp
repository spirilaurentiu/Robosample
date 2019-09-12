#ifndef __GIROLAMISAMPLER_HPP__
#define __GIROLAMISAMPLER_HPP__

#include "HamiltonianMonteCarloSampler.hpp"

class Topology;
class IState;
void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

class GirolamiSampler : public HamiltonianMonteCarloSampler
{
public:

    // Constructor
    GirolamiSampler(SimTK::CompoundSystem *argCompoundSystem,
                                 SimTK::SimbodyMatterSubsystem *argMatter,
                                 SimTK::Compound *argResidue,
                                 SimTK::DuMMForceFieldSubsystem *argDumm,
                                 SimTK::GeneralForceSubsystem *forces,
                                 SimTK::TimeStepper *argTimeStepper);

    // Destructor
    virtual ~GirolamiSampler();

    // Initialize variables (like TVector)
    virtual void initialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman = true) ; 

    // Initialize variables (like TVector)
    virtual void reinitialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) ; 

    // Assign a random conformation. Time measured in picoseconds
    void propose(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State& someState, SimTK::Real timestep, int nosteps);

protected:
    bool useFixmanTorque;

};

#endif // __GIROLAMISAMPLER_HPP__

