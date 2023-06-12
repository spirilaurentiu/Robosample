#ifndef __TESTHMCSOA_HPP__
#define __TESTHMCSOA_HPP__

#include "MonteCarloSampler.hpp"

class Topology;
class IState;
void writePdb(      SimTK::Compound& c, SimTK::State& advanced,
         const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime);

class TestHMCSOA : public MonteCarloSampler
{
public:

    // Constructor
    TestHMCSOA(SimTK::CompoundSystem *argCompoundSystem,
                                 SimTK::SimbodyMatterSubsystem *argMatter,
                                 SimTK::Compound *argResidue,
                                 SimTK::DuMMForceFieldSubsystem *argDumm,
                                 SimTK::GeneralForceSubsystem *forces,
                                 SimTK::TimeStepper *argTimeStepper);

    // Destructor
    virtual ~TestHMCSOA();

    // Calculate O(n2) the square root of the mass matrix inverse
    void calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv);
    void calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv);
    // Calculate sqrt(M) using Eigen
    void calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper);

    // Initialize variables (like TVector)
    virtual void initialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature, bool argUseFixman = true) ; 

    // Initialize variables (like TVector)
    virtual void reinitialize(SimTK::State& advanced, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) ; 

    // Assign a random conformation. Time measured in picoseconds
    void propose(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    // Performs the acception-rejection step and sets the state of the compound
    // to the appropriate conformation
    void update(SimTK::State& someState, SimTK::Real timestep, int nosteps);

    // Get set kinetic energy
    SimTK::Real getOldKE(void) { return this->ke_o; }
    
    // Get old kinetic energy
    SimTK::Real getSetKE(void) { return this->ke_set; }
    
    // Set old kinetic energy
    void setOldKE(SimTK::Real);

    // Set set kinetic energies    
    void setSetKE(SimTK::Real);

protected:
    SimTK::Real ke_set; // The kinetic energy retained for acc-rej step
    SimTK::Real ke_o; // The kinetic energy retained for acc-rej step
    SimTK::Real etot_set; // The kinetic energy retained for acc-rej step
    SimTK::Real etot_o; // The kinetic energy retained for acc-rej step
    // TO BE DELETED
    SimTK::Matrix prevM;
    SimTK::Real prevThetaK;
    SimTK::Vector prevTheta;
    SimTK::Real prevNumDetM;
    SimTK::Real prevDetM;
    int kForTheta;
    // TO BE DELETED

};

#endif // __TESTHMCSOA_HPP__

