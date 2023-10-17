#include "ThermodynamicState.hpp"

ThermodynamicState::ThermodynamicState(int index)
{
	myIndex = index;
	temperature = 300;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;

	nonequilibrium = 0;
}

ThermodynamicState::ThermodynamicState(int index, const SimTK::Real& T,
		const std::vector<int>& argWorldIndexes,
		const std::vector<SimTK::Real>& argTimesteps,
		const std::vector<int>& argMdsteps
	)
{
	// Own index
	myIndex = index;

	// Temperature related
	temperature = T;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;

	// Worlds related parameters
	worldIndexes = argWorldIndexes;
	timesteps = argTimesteps;
	mdsteps = argMdsteps;

	nonequilibrium = 0;
}

void ThermodynamicState::setTemperature(SimTK::Real T)
{
	temperature = T;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
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

//
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
	this->rexSamplers = rexSamplersArg;
}

const std::vector<std::string>& ThermodynamicState::getSamplers() const
{
	return this->rexSamplers;
}

// Next functions set Q, U, tau perturbing functions options for samplers
void ThermodynamicState::setDistortOptions(const std::vector<int>& rexDistortOptionsArg)
{
	this->rexDistortOptions = rexDistortOptionsArg;
}

std::vector<int>& ThermodynamicState::getDistortOptions()
{
	return this->rexDistortOptions;
}

void ThermodynamicState::setDistortArgs(const std::vector<std::string>& rexDistortOptionsArg)
{
	this->rexDistortArgs = rexDistortOptionsArg;
}

std::vector<std::string>& ThermodynamicState::getDistortArgs()
{
	return this->rexDistortArgs;
}

void ThermodynamicState::setFlowOptions(const std::vector<int>& rexFlowOptionsArg)
{
	this->rexFlowOptions = rexFlowOptionsArg;

}

void ThermodynamicState::setWorkOptions(const std::vector<int>& rexWorkOptionsArg)
{
	this->rexWorkOptions = rexWorkOptionsArg;
}

// Set the integrating method
void ThermodynamicState::setIntegrators(const std::vector<std::string>& rexIntegratorsArg)
{
	this->rexIntegrators = rexIntegratorsArg;
}

const std::vector<std::string>& ThermodynamicState::getIntegrators() const
{
	return rexIntegrators;
}

void ThermodynamicState::Print()
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
