#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "Robo.hpp"
#include <unordered_map>

// This class implements the functionality for reading arguments from a setup
// file and give back the necessary values as vectors
class SetupReader {
public:

	// Constructor
	SetupReader() {};
	SetupReader(const char *FN);
	SetupReader(const std::string& FN);

	// Read setup function
	void ReadSetup(const char *FN);
	void ReadSetup(const std::string& FN);

	// Print all the arguments
	void dump(bool PrettyPrint) const;

	// Access values by key
	bool find(const char *argKey) const;
	bool find(const std::string& argKey) const;

	const std::vector<std::string>& get(const char *argKey) const;
	const std::vector<std::string>& get(const std::string& argKey) const;
	
	// 
	const size_t getNofKeys(void) const;

	bool IsPrintable(const std::string& s) const;

	std::vector<std::string> split(const std::string& i_str, const std::string& i_delim);

	/////////////////////////
	// REX
	/////////////////////////

	// Get the number of replica requested
	int readREXConfigFile(std::string FN,
		std::vector<SimTK::Real>& temperatures,		
		
		std::vector<std::vector<std::string>>& rexSamplers,
		std::vector<std::vector<int>>& rexDistortOptions,
		std::vector<std::vector<std::string>>& rexDistortArgs,
		std::vector<std::vector<int>>& rexFlowOptions,
		std::vector<std::vector<int>>& rexWorkOptions,
		std::vector<std::vector<std::string>>& rexIntegrators,

		std::vector<std::vector<SimTK::Real>>& rexTimesteps,
		std::vector<std::vector<int>>& rexWorldIndexes,
		std::vector<std::vector<int>>& rexMdsteps,
		std::vector<std::vector<int>>& rexSamplesPerRound);

	const std::string& getRestartDirectory(void);

private:
	std::map<std::string, std::vector<std::string>> Args;
	std::vector<std::string> KeyNotFound = { "ERROR_KEY_NOT_FOUND" };

	/////////////////////////
	// REX
	/////////////////////////

	// REX config file parameter types
	const std::unordered_map<std::string, int> rexToIntKeys{
		{"TEMPERATURE", 0},

		{"SAMPLERS", 1},
		{"DISTORT_OPTIONS", 2},
		{"DISTORT_ARGS", 3},
		{"FLOW_OPTIONS", 4},
		{"WORK_OPTIONS", 5},
		{"INTEGRATORS", 6},

		{"TIMESTEPS", 7},
		{"WORLD_INDEXES", 8},
		{"MDSTEPS", 9},
		{"SAMPLES_PER_ROUND", 10}
	};

	enum RexKey{
		TEMPERATURE,

		SAMPLERS,
		DISTORT_OPTIONS,
		DISTORT_ARGS,
		FLOW_OPTIONS,
		WORK_OPTIONS,
		INTEGRATORS,
		
		TIMESTEPS,
		WORLD_INDEXES,
		MDSTEPS,
		SAMPLES_PER_ROUND
	};
	
	std::string restartDirectory = "";
};

