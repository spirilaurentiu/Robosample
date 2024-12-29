#pragma once

#include <vector>
#include <string>
#include "Simbody.h"
#include "bSpecificAtom.hpp"
#include "TrajectoryObject.hpp"

constexpr SimTK::Real DEFAULT_TEMPERATURE = 300.0;

class ThermodynamicState{
  public:
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
	void setTimesteps(const std::vector<SimTK::Real>& timesteps);
	const std::vector<SimTK::Real>& getTimesteps() const;
	void setMdsteps(const std::vector<int>& mdsteps);
	const std::vector<int>& getMdsteps() const;

	// Set the sampling method
	void setAcceptRejectModes(const std::vector<AcceptRejectMode>& acceptRejectModes);
	const std::vector<AcceptRejectMode>& getAcceptRejectModes() const;

	// Next functions set Q, U, tau perturbing functions options
	// for samplers
	void setDistortOptions(const std::vector<int>& rexDistortOptionsArg);
	std::vector<int>& getDistortOptions();
	void setDistortArgs(const std::vector<std::string>& rexDistortOptionsArg);
	const std::vector<std::string>& getDistortArgs() const;
	void setFlowOptions(const std::vector<int>& rexFlowOptionsArg);
	void setWorkOptions(const std::vector<int>& rexWorkOptionsArg);

	// Set the integrating method
	void setIntegrators(const std::vector<IntegratorType>& rexIntegratorsArg);
	const std::vector<IntegratorType>& getIntegrators() const;

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
	int getNofWorldsSamples() const;
    int getNofSamples() const {return nofSamples;};

    // Incrementer function for nofSamples
    void incrementWorldsNofSamples();
    void incrementWorldsNofSamples(int howMany);	
    void incrementNofSamples(){nofSamples++;}
	void incrementNofSamples(int howMany){nofSamples += howMany;}
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

	int findWorld(const int whichWorld);

	void allocQStatsFirstDimension(void);
	bool calcQStats(const int whichWorld, const SimTK::Vector & worldBMps, const SimTK::Vector & worldQs, int worldNofSamples);
	void printQStats(void);

	std::vector<SimTK::Real>& getBMps_means(const int whichWorld);
	std::vector<SimTK::Real>& get_dBMps(const int whichWorld);

	std::vector<SimTK::Real>& getCurrentQs(const int whichWorld);
	std::vector<SimTK::Real>& getQmeans(const int whichWorld);
	std::vector<SimTK::Real>& getQdiffs(const int whichWorld);
	std::vector<SimTK::Real>& getQvars(const int whichWorld);

	void appendLog(const std::string& filename);
	void appendDCDReporter(const std::string& filename, int natoms, int ntopologies);
	void writeDCD(std::vector<SimTK::Real>& x, std::vector<SimTK::Real>& y, std::vector<SimTK::Real>& z);
	
	std::ofstream logFile;
	TrajectoryObject traj;

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

	std::vector<AcceptRejectMode> rexAcceptRejectModes;
	std::vector<int> rexDistortOptions;
	std::vector<std::string> rexDistortArgs;
	std::vector<int> rexFlowOptions;
	std::vector<int> rexWorkOptions;
	std::vector<IntegratorType> rexIntegrators;

	int allWorldsNofSamples = 0;
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
	std::vector<std::vector<SimTK::Real>> BMps_means; // vector[world][qindex] 
	std::vector<std::vector<SimTK::Real>> dBMps; // vector[world][qindex] 

	std::vector<std::vector<SimTK::Real>> currQs; // vector[world][qindex] 
	std::vector<std::vector<SimTK::Real>> Qdiffs; // vector[world][qindex]
	std::vector<std::vector<SimTK::Real>> Qmeans; // vector[world][qindex]
	std::vector<std::vector<SimTK::Real>> Qvars;  // vector[world][qindex]

	//////////////////////////////////
	//---         Q Stats        -----
	//////////////////////////////////	

};
