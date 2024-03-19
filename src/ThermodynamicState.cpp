#include "ThermodynamicState.hpp"

ThermodynamicState::ThermodynamicState(
		int index,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_,
		std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
		) :
	atoms(atoms_),
	zMatrixTable(zMatrixTable_),
	zMatrixBAT_ref(zMatrixBAT_ref_)
{
	myIndex = index;
	temperature = 300;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;

	nonequilibrium = 0;
}

ThermodynamicState::ThermodynamicState(
		int index,
		const SimTK::Real& T,
		const std::vector<int>& argWorldIndexes,
		const std::vector<SimTK::Real>& argTimesteps,
		const std::vector<int>& argMdsteps,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_,
		std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
	) :
	atoms(atoms_),
	zMatrixTable(zMatrixTable_),
	zMatrixBAT_ref(zMatrixBAT_ref_)
	
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


// BAT
/*!
 * <!--	zmatrixbat_ Getter function implementation -->
*/
int ThermodynamicState::getNofSamples() const {
    return nofSamples;
}

/*!
 * <!--	zmatrixbat_ Setter function implementation -->
*/
void ThermodynamicState::setNofSamples(int newNofSamples) {
    nofSamples = newNofSamples;
}

/*!
 * <!--	zmatrixbat_ Incrementer function implementation -->
*/
void ThermodynamicState::incrementNofSamples() {
    ++nofSamples;
}


/*!
 * <!--	zmatrixbat_ -->
*/
void ThermodynamicState::PrintZMatrixBAT() const {

	int bati = 0;
	for (const auto& row : zMatrixBAT_ref) {

		for(const auto tabValue : zMatrixTable[bati]){
			std::cout << tabValue << " ";
		}

		scout("aIxs ");

		for(const auto tabValue : zMatrixTable[bati]){
			if( tabValue >= 0){
				std::cout << atoms[tabValue].getCompoundAtomIndex() << " ";
			}else{
				std::cout << "dummy " ;
			}
		}

		for (const SimTK::Real value : row) {
			std::cout << std::setw(6) << value << " ";
		}
		std::cout << std::endl;

		bati++;
	}

}

/*!
 * <!--	zmatrixbat_ -->
*/
void
ThermodynamicState::calcZMatrixBATStats(void)
{
	SimTK::Real N = nofSamples + 1;

	if(nofSamples == 0){

        // Initialize zMatrixBATMeans, zMatrixBATDiffs, and zMatrixBATStds
        zMatrixBATMeans.resize(zMatrixBAT_ref.size(), std::vector<SimTK::Real>(3, 0));
        zMatrixBATDiffs.resize(zMatrixBAT_ref.size(), std::vector<SimTK::Real>(3, 0));
        zMatrixBATStds.resize(zMatrixBAT_ref.size(), std::vector<SimTK::Real>(3, 0));

	}else{

		// Iterate bats
	    for (size_t bati = 0; bati < zMatrixBAT_ref.size(); ++bati) {

			std::vector<SimTK::Real>& BAT = zMatrixBAT_ref[bati];
			std::vector<SimTK::Real>& BATmeans = zMatrixBATMeans[bati];
			std::vector<SimTK::Real>& BATdiffs = zMatrixBATDiffs[bati];
			std::vector<SimTK::Real>& BATstds = zMatrixBATStds[bati];			

			// Usefull vars
			SimTK::Real N_1_over_N = (N - 1.0) / N;
			SimTK::Real Ninv = 1.0 / N;
			
			// Update BAT means
			for (int i = 0; i < 3; ++i) {
				BATmeans[i] = (N_1_over_N * BATmeans[i]) + (Ninv * BAT[i]);
			}

			// Update BAT differences
			for (int i = 0; i < 3; ++i) {
				BATdiffs[i] = BAT[i] - BATmeans[i];
			}

			// Update BAT standard deviations
			for (int i = 0; i < 3; ++i) {
				BATstds[i] = std::sqrt((N_1_over_N * BATstds[i]) + (Ninv * (BATdiffs[i] * BATdiffs[i])));
			}

			// Keep track of BAT
			bati++;
		}
	}	
}


/*!
 * <!--	zmatrixbat_ -->
*/
std::vector<SimTK::Real>&
ThermodynamicState::getBATMeansRow(int rowIndex) 
{

	assert(rowIndex >= 0 && rowIndex < zMatrixBATMeans.size()
		&& "Invalid row index");

	return zMatrixBATMeans[rowIndex];
}


/*!
 * <!--	zmatrixbat_ -->
*/
std::vector<SimTK::Real>&
ThermodynamicState::getBATDiffsRow(int rowIndex) 
{

	assert(rowIndex >= 0 && rowIndex < zMatrixBATDiffs.size()
		&& "Invalid row index");

	return zMatrixBATDiffs[rowIndex];
}


/*!
 * <!--	zmatrixbat_ -->
*/
std::vector<SimTK::Real>&
ThermodynamicState::getBATStdsRow(int rowIndex) 
{

	assert(rowIndex >= 0 && rowIndex < zMatrixBATStds.size()
		&& "Invalid row index");

	return zMatrixBATStds[rowIndex];
}

