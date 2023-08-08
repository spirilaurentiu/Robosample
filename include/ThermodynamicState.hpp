#pragma once

#include <vector>
#include <string>
#include "Simbody.h"

constexpr SimTK::Real DEFAULT_TEMPERATURE = 300.0;

class ThermodynamicState {
public:
	ThermodynamicState(int index);

	ThermodynamicState(int index,
        SimTK::Real T,
	    std::vector<int>& argWorldIndexes,
	    std::vector<SimTK::Real>& argTimesteps,
	    std::vector<int>& argMdsteps
	);

	int getIndex() const;
	void setIndex(int someIndex);

	void setTemperature(SimTK::Real T);
	SimTK::Real getTemperature() const;
	SimTK::Real getBeta() const;

	const std::vector<int>& getWorldIndexes() const;
	std::vector<int>& updWorldIndexes();
	const std::vector<SimTK::Real>& getTimesteps() const;
	const std::vector<int>& getMdsteps() const;

	// Set the sampling method
	void setSamplers(const std::vector<std::string>& rexSamplersArg);
	const std::vector<std::string>& getSamplers() const;

	// Next functions set Q, U, tau perturbing functions options
	// for samplers
	void setDistortOptions(const std::vector<int>& rexDistortOptionsArg);
	const std::vector<int>& getDistortOptions() const;
	void setFlowOptions(const std::vector<int>& rexFlowOptionsArg);
	void setWorkOptions(const std::vector<int>& rexWorkOptionsArg);

	// Set the integrating method
	void setIntegrators(const std::vector<std::string>& rexIntegratorsArg);
	const std::vector<std::string>& getIntegrators() const;

	// If any of the distort, flow or work options are active
	int hasNonequilibriumMoves() const;
	void setNonequilibrium(int nArg);

	// Print everything about thermostate
	void Print() const;

private:

	// Index
	int myIndex = 0;

	// Temperature
	SimTK::Real temperature = DEFAULT_TEMPERATURE;
    SimTK::Real RT = DEFAULT_TEMPERATURE * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	SimTK::Real beta = 1.0 / RT;

	// Worlds related parameters 
	std::vector<int> worldIndexes; 
	std::vector<SimTK::Real> timesteps; 
	std::vector<int> mdsteps; 

	int nonequilibrium = 0;

	std::vector<std::string> rexSamplers;
	std::vector<int> rexDistortOptions;
	std::vector<int> rexFlowOptions;
	std::vector<int> rexWorkOptions;
	std::vector<std::string> rexIntegrators;
};