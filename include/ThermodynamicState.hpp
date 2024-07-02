#pragma once

#include <vector>
#include <string>
#include "Simbody.h"
#include "bSpecificAtom.hpp"

constexpr SimTK::Real DEFAULT_TEMPERATURE = 300.0;

class ThermodynamicState{
  public:
	ThermodynamicState() = default;
	
	ThermodynamicState(
		int index,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_
		//,std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
	);

	ThermodynamicState(int index, const SimTK::Real& T,
		const std::vector<int>& argWorldIndexes,
		const std::vector<SimTK::Real>& argTimesteps,
		const std::vector<int>& argMdsteps,

		std::vector<bSpecificAtom>& atoms_,
		std::vector<std::vector<int>>& zMatrixTable_
		//, std::vector<std::vector<SimTK::Real>>& zMatrixBAT_ref_
	);

	const int getIndex(){return myIndex;}
	void setIndex(const int& someIndex){myIndex = someIndex;}

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
	std::vector<int>& getDistortOptions();
	void setDistortArgs(const std::vector<std::string>& rexDistortOptionsArg);
	std::vector<std::string>& getDistortArgs();
	void setFlowOptions(const std::vector<int>& rexFlowOptionsArg);
	void setWorkOptions(const std::vector<int>& rexWorkOptionsArg);

	// Set the integrating method
	void setIntegrators(const std::vector<std::string>& rexIntegratorsArg);
	const std::vector<std::string>& getIntegrators() const;

	// If any of the distort, flow or work options are active
	const int hasNonequilibriumMoves(){
		return this->nonequilibrium;}
	void setNonequilibrium(int nArg){this->nonequilibrium = nArg;}

	// Print everything about thermostate
	void Print();


    /** @name Z Matrix and BAT functions
	*/

    /**@{**/

    // Getter function for nofSamples
    int getNofSamples() const;

    // Setter function for nofSamples
    void setNofSamples(int newNofSamples);

    // Incrementer function for nofSamples
    void incrementNofSamples();
    void incrementNofSamples(int howMany);	

	void PrintZMatrixBAT(bool printBATStats = false) const ;

	//
	void calcZMatrixBATStats(void);

	// Get a row from the BAT means
	std::vector<SimTK::Real>& getBATMeansRow(int rowIndex);	

	// Get a row from the BAT diffs
	std::vector<SimTK::Real>& getBATDiffsRow(int rowIndex);	

	// Get a row from the BAT stds
	std::vector<SimTK::Real>& getBATVarsRow(int rowIndex);	

    // Setter function to update the zMatrixBAT_poi member
    void setZMatrixBATPointer(std::vector<std::vector<SimTK::Real>>& pointer) {
        zMatrixBAT_poi = &pointer;
    }

	/**@}**/

    
	/** @name Q Stats functions
	*/

    /**@{**/
	
	void allocateQStats(int nofWorlds){

		//allQs.resize(nofWorlds);
		
	}

	void setWorldQs(int whichWorld, const SimTK::Vector & worldQs){

		//allQs[whichWorld] = worldQs;

	}

	void allocQStatsFirstDimension(void);
	bool calcQStats(const int whichWorld, const SimTK::Vector & worldQs);
	void printQStats(void);

	/**@}**/	
  
  
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
	std::vector<std::string> rexDistortArgs;
	std::vector<int> rexFlowOptions;
	std::vector<int> rexWorkOptions;
	std::vector<std::string> rexIntegrators;

	int nofSamples = 0;

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	std::vector<bSpecificAtom>& atoms;

	std::vector<std::vector<int>>& zMatrixTable;
	std::vector<std::vector<SimTK::Real>>* zMatrixBAT_poi;

	std::vector<std::vector<SimTK::Real>>  zMatrixBATMeans;
	std::vector<std::vector<SimTK::Real>>  zMatrixBATDiffs;
	std::vector<std::vector<SimTK::Real>>  zMatrixBATVars;

	//////////////////////////////////
	/////      Z Matrix BAT      /////
	//////////////////////////////////

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////

	//std::vector< SimTK::Vector & > allQs;

	std::vector<std::vector<SimTK::Real>> Qdiffs; // vector[world][qindex]
	std::vector<std::vector<SimTK::Real>> Qmeans; // vector[world][qindex]
	std::vector<std::vector<SimTK::Real>> Qvars;  // vector[world][qindex]

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////	

};
