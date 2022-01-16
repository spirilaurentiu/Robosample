#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

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

private:
	std::map<std::string, std::vector<std::string>> Args;
	std::vector<std::string> KeyNotFound = { "ERROR_KEY_NOT_FOUND" };

	bool IsPrintable(const std::string& s) const;
};

