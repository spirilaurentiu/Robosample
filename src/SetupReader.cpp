#include "SetupReader.hpp"

// Constructor
SetupReader::SetupReader(const char *cFN)
{
	ReadSetup(cFN);
}

// Constructor
SetupReader::SetupReader(const std::string& FN)
{
	ReadSetup(FN);
}

// Read setup function
void SetupReader::ReadSetup(const char *cFN)
{
	ReadSetup(std::string(cFN));
}

void SetupReader::ReadSetup(const std::string& FN)
{
	// TODO return errors

	std::ifstream F(FN);
	std::string line;

	// read file
	while (F) {
		// read line
		std::getline(F, line);
		std::istringstream iss(line);

		// read key
		std::string key;
		iss >> key;

		// check if key has characters, all of them are printable and does not begin a comment
		if (key.length() > 0 && IsPrintable(key) && '#' != key[0]) {
			std::vector<std::string> V;
			std::string word;

			// read arguments of key
			while (iss >> word) {
				// stop on comment
				if ('#' == word[0]) {
					break;
				}

				// add word if printable
				if (IsPrintable(word)) {
					V.push_back(std::move(word));
				}
			}

			// save key and arguments
			Args.emplace(std::make_pair(std::move(key), std::move(V)));
		}
	}
}

// Print all the arguments
void SetupReader::dump(bool PrettyPrint) const
{
	if(PrettyPrint) {
		std::cout << "Dumping input file:\n";
	}

	// TODO replace cout with other streams (as param)
	for (const auto& Arg : Args) {
		if(PrettyPrint) {
			std::cout << "  ";
		}

		std::cout << Arg.first << " : ";
		for (const auto& Val : Arg.second) {
			std::cout << Val << " ";
		}

		std::cout << '\n';
	}
}

/** Check if key exists **/
bool SetupReader::find(const char *argKey) const
{
	return find(std::string(argKey));
}

bool SetupReader::find(const std::string& argKey) const
{
	// TODO std::optional
	return Args.find(argKey) != Args.end();
}

// Access values by key
const std::vector<std::string>& SetupReader::get(const char *cArgKey) const
{
	return get(std::string(cArgKey));
}

// Access values by key
const std::vector<std::string>& SetupReader::get(const std::string& argKey) const
{
	const auto it = Args.find(argKey);
	if (it != Args.cend()) {
		return it->second;
	}
	else {
        std::cout << "SetupReader::get() did not find '" << argKey << "'. Returning " << KeyNotFound[0] << ".\n";
		return KeyNotFound;
	}
}

bool SetupReader::IsPrintable(const std::string& s) const
{
	return std::all_of(s.begin(), s.end(), [](char c) { return c >= 0 && c <= 128; });
}