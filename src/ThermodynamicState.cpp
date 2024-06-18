#include "ThermodynamicState.hpp"

ThermodynamicState::ThermodynamicState(
		int index,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_
		//, std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
		) :
	atoms(atoms_),
	zMatrixTable(zMatrixTable_)
	//, zMatrixBAT_poi(zMatrixBAT_ref_)
{
	myIndex = index;
	temperature = 300;
	SimTK::Real RT = temperature *
		static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
	beta = 1.0 / RT;

	nonequilibrium = 0;

	allocQStatsFirstDimension();

}

ThermodynamicState::ThermodynamicState(
		int index,
		const SimTK::Real& T,
		const std::vector<int>& argWorldIndexes,
		const std::vector<SimTK::Real>& argTimesteps,
		const std::vector<int>& argMdsteps,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_
		//, std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
	) :
	atoms(atoms_),
	zMatrixTable(zMatrixTable_)
	//, zMatrixBAT_poi(zMatrixBAT_ref_)
	
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

	allocQStatsFirstDimension();

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
void ThermodynamicState::PrintZMatrixBAT(bool printBATStats) const {

	int bati = 0;

	for (const auto& row : *zMatrixBAT_poi) {
		
		std::cout << "thermoIx|" << myIndex <<"| ";
		for(const auto tabValue : zMatrixTable[bati]){
			std::cout << tabValue << " ";
		}

		scout("cAIxs ");

		for(const auto tabValue : zMatrixTable[bati]){
			if( tabValue >= 0){
				std::cout << atoms[tabValue].getCompoundAtomIndex() << " ";
			}else{
				std::cout << "dummy " ;
			}
		}

		// BAT values
		for (const SimTK::Real value : row) {
			std::cout << std::setprecision(12) << value << " ";
		}

		if(nofSamples > 0){

			if(printBATStats){
				
				// BAT means
				for(size_t dimi = 0; dimi < 3; dimi++){
					std::cout << zMatrixBATMeans[bati][dimi] <<" ";
				}
				// BAT diffs
				for(size_t dimi = 0; dimi < 3; dimi++){
					std::cout << zMatrixBATDiffs[bati][dimi] <<" ";
				}
				// BAT stds
				for(size_t dimi = 0; dimi < 3; dimi++){
					std::cout << zMatrixBATVars[bati][dimi] <<" ";
				}
			}

		}

		std::cout << std::endl;

		bati++;
	}

	std::cout << "=====" << std::endl;

}

/*!
 * <!--	zmatrixbat_ -->
*/
void
ThermodynamicState::calcZMatrixBATStats(void)
{

	// Usefull vars
	SimTK::Real N = nofSamples + 1;
	SimTK::Real N_1_over_N = (N - 1.0) / N;
	SimTK::Real Ninv = 1.0 / N;
			
	if(nofSamples == 0){

        // Initialize zMatrixBATMeans, zMatrixBATDiffs, and zMatrixBATStds
        zMatrixBATMeans.resize((*zMatrixBAT_poi).size(), std::vector<SimTK::Real>(3, 0));
        zMatrixBATDiffs.resize((*zMatrixBAT_poi).size(), std::vector<SimTK::Real>(3, 0));
        zMatrixBATVars.resize((*zMatrixBAT_poi).size(), std::vector<SimTK::Real>(3, 0));

	    for (size_t bati = 0; bati < (*zMatrixBAT_poi).size(); ++bati) {

			zMatrixBATMeans[bati][0] = (*zMatrixBAT_poi)[bati][0];
			zMatrixBATMeans[bati][1] = (*zMatrixBAT_poi)[bati][1];
			zMatrixBATMeans[bati][2] = (*zMatrixBAT_poi)[bati][2];

			zMatrixBATDiffs[bati][0] = 0;
			zMatrixBATDiffs[bati][1] = 0;
			zMatrixBATDiffs[bati][2] = 0;

			zMatrixBATVars[bati][0] = 0;
			zMatrixBATVars[bati][1] = 0;
			zMatrixBATVars[bati][2] = 0;

		}

	}else{ // nofSamples gt 2

		// Iterate bats
	    for (size_t bati = 0; bati < (*zMatrixBAT_poi).size(); ++bati) {

			std::vector<SimTK::Real>& BAT = (*zMatrixBAT_poi)[bati];

			std::vector<SimTK::Real>& BATmeans = zMatrixBATMeans[bati];
			std::vector<SimTK::Real>& BATdiffs = zMatrixBATDiffs[bati];
			std::vector<SimTK::Real>& BATstds = zMatrixBATVars[bati];			

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
				BATstds[i] = (N_1_over_N * BATstds[i]) + (Ninv * (BATdiffs[i] * BATdiffs[i]));
			}

		} // every BAT entry
	} // nofSamples gt 2

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
ThermodynamicState::getBATVarsRow(int rowIndex) 
{

	assert(rowIndex >= 0 && rowIndex < zMatrixBATVars.size()
		&& "Invalid row index");

	return zMatrixBATVars[rowIndex];
}

/*!
 * <!--  -->
*/
void ThermodynamicState::allocQStatsFirstDimension(void)
{
        // Resize
		Qmeans.resize(worldIndexes.size());
		Qdiffs.resize(worldIndexes.size());
		Qvars.resize(worldIndexes.size());
}

/*!
 * <!--  -->
*/
bool ThermodynamicState::calcQStats(const int whichWorld, const SimTK::Vector & worldQs)
{

	// Usefull vars
	SimTK::Real N = nofSamples + 1;
	SimTK::Real N_1_over_N = (N - 1.0) / N;
	SimTK::Real Ninv = 1.0 / N;

	// Search the position in cpp vector
	bool found = false;
	int vPosInVector = -1;

	for(const auto wIx : worldIndexes){
		vPosInVector++;

		if(whichWorld == wIx){

			vPosInVector = wIx;
			found = true;
			break;
		}
	}

	if(found == false){
		return false;
	}

	if(nofSamples == 0){

		// Resize
        Qmeans[vPosInVector].resize(worldQs.size());
        Qdiffs[vPosInVector].resize(worldQs.size());
        Qvars[vPosInVector].resize(worldQs.size());

		// Initialize at the first sample
		for(int qIx = 0; qIx < worldQs.size(); qIx++){

			Qmeans[vPosInVector][qIx] = worldQs[qIx];

			Qdiffs[vPosInVector][qIx] = 0;

			Qvars[vPosInVector][qIx] = 0;

		}

	}else{ // nofSamples gt 2

		// Update Q means
		for(int qIx = 0; qIx < worldQs.size(); qIx++){

			Qmeans[vPosInVector][qIx] = (N_1_over_N * Qmeans[vPosInVector][qIx]) + (Ninv * worldQs[qIx]);
		}

		// Update Q differences
		for(int qIx = 0; qIx < worldQs.size(); qIx++){
			Qdiffs[vPosInVector][qIx] = worldQs[qIx] - Qmeans[vPosInVector][qIx];
		}

		// Update Q variances
		for(int qIx = 0; qIx < worldQs.size(); qIx++){
			Qvars[vPosInVector][qIx] = (N_1_over_N * Qvars[vPosInVector][qIx]) + (Ninv * (Qdiffs[vPosInVector][qIx] * Qdiffs[vPosInVector][qIx]));
		}

	} // nofSamples gt 2

	return true;
	
}

/*!
 * <!--  -->
*/
void ThermodynamicState::printQStats(void)
{
	int wPosInVector = -1;
	for(const auto wIx : worldIndexes){
		wPosInVector++;

		std::cout << "wIx " << wIx <<": ";
		for(int qIx = 0; qIx < Qmeans[wPosInVector].size(); qIx++){
			std::cout <<" " << Qmeans[wPosInVector][qIx];
		}
		std::cout << std::endl;
	}

}

