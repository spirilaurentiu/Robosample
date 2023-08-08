#include "ThermodynamicState.hpp"

ThermodynamicState::ThermodynamicState(int index) : myIndex(index)
{
}

ThermodynamicState::ThermodynamicState(int index,
        SimTK::Real T,
	    std::vector<int>& argWorldIndexes,
	    std::vector<SimTK::Real>& argTimesteps,
	    std::vector<int>& argMdsteps
	) : myIndex(index)
{
	// Worlds related parameters
	worldIndexes = argWorldIndexes;
	timesteps = argTimesteps;
	mdsteps = argMdsteps;

	nonequilibrium = 0;
}

int ThermodynamicState::getIndex() const {
    return myIndex;
}

void ThermodynamicState::setIndex(int someIndex) {
    myIndex = someIndex;
}

void ThermodynamicState::setTemperature(SimTK::Real T)
{
	temperature = T;
	SimTK::Real RT = temperature * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;
}

SimTK::Real ThermodynamicState::getTemperature() const
{
	return temperature;
}

SimTK::Real ThermodynamicState::getBeta() const
{
	return beta;
}

const std::vector<int>& ThermodynamicState::getWorldIndexes() const
{
	return worldIndexes;
}

std::vector<int>& ThermodynamicState::updWorldIndexes()
{
	return worldIndexes;
}

const std::vector<SimTK::Real>& ThermodynamicState::getTimesteps() const
{
	return timesteps;
}

const std::vector<int>& ThermodynamicState::getMdsteps() const
{
	return mdsteps;
}

// Set the sampling method
void ThermodynamicState::setSamplers(const std::vector<std::string>& rexSamplersArg)
{
	rexSamplers = rexSamplersArg;
}

const std::vector<std::string>& ThermodynamicState::getSamplers() const
{
	return rexSamplers;
}

// Next functions set Q, U, tau perturbing functions options for samplers
void ThermodynamicState::setDistortOptions(const std::vector<int>& rexDistortOptionsArg)
{
	rexDistortOptions = rexDistortOptionsArg;
}

const std::vector<int>& ThermodynamicState::getDistortOptions() const
{
	return rexDistortOptions;
}

void ThermodynamicState::setFlowOptions(const std::vector<int>& rexFlowOptionsArg)
{
	rexFlowOptions = rexFlowOptionsArg;

}

void ThermodynamicState::setWorkOptions(const std::vector<int>& rexWorkOptionsArg)
{
	rexWorkOptions = rexWorkOptionsArg;
}

// Set the integrating method
void ThermodynamicState::setIntegrators(const std::vector<std::string>& rexIntegratorsArg)
{
	rexIntegrators = rexIntegratorsArg;
}

const std::vector<std::string>& ThermodynamicState::getIntegrators() const
{
	return rexIntegrators;
}

int ThermodynamicState::hasNonequilibriumMoves() const {
    return nonequilibrium;
}

void ThermodynamicState::setNonequilibrium(int nArg) {
    nonequilibrium = nArg;
}

void ThermodynamicState::Print() const
{
	std::cout << "ThermodynamicState::Print index T "
		<< myIndex << " " << temperature
		<< std::endl;

	std::cout << "ThermodynamicState::Print index worldIndexes " << myIndex;
	for(auto worldIndex : worldIndexes){
		std::cout << " " << worldIndex;
	}
	std::cout << std::endl;
}