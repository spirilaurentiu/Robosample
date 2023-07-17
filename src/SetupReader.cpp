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
		// throw std::exception(); // TODO do we throw here? it would be ub not to, you know...
		
		return KeyNotFound;
	}
}

//
const size_t SetupReader::getNofKeys(void) const
{
	return Args.size();
}

bool SetupReader::IsPrintable(const std::string& s) const
{
	return std::all_of(s.begin(), s.end(), [](char c) { return c >= 0 && c <= 128; });
}

// TODO: Victor: Can we trust this function?
std::vector<std::string> SetupReader::split(const std::string& i_str, const std::string& i_delim)
{
    std::vector<std::string> result;
    
    size_t found = i_str.find(i_delim);
    size_t startIndex = 0;

    while(found != std::string::npos)
    {
        result.push_back(std::string(i_str.begin()+startIndex, i_str.begin()+found));
        startIndex = found + i_delim.size();
        found = i_str.find(i_delim, startIndex);
    }
    if(startIndex != i_str.size())
        result.push_back(std::string(i_str.begin()+startIndex, i_str.end()));
    return result;      
}


// Get the number of replica requested
int SetupReader::readREXConfigFile(std::string FN,
		std::vector<SimTK::Real>& temperatures,

		std::vector<std::vector<std::string>>& rexSamplers,
		std::vector<std::vector<int>>& rexDistortOptions,
		std::vector<std::vector<int>>& rexFlowOptions,
		std::vector<std::vector<int>>& rexWorkOptions,
		std::vector<std::vector<std::string>>& rexIntegrators,

		std::vector<std::vector<SimTK::Real>>& rexTimesteps,
		std::vector<std::vector<int>>& rexWorldIndexes,
		std::vector<std::vector<int>>& rexMdsteps,
		std::vector<std::vector<int>>& rexSamplesPerRound)
{
	std::ifstream F(FN);
	std::string line;
	size_t nofReplicas = 0;

	// Get the number of replicas
	while (F) {
		std::getline(F, line);
		if(line.length() > 0){
			std::vector<std::string> words = split(line, " ");
			if((words[0][0] != '#') && (words.size() >= 2) && (words[0] == "NOF_REPLICAS")){
				nofReplicas = std::stoi(words[1]);
			}
		}
	}
	std::cout << "Number of replicas requested " << nofReplicas << " \n" ;

	// Rewind the file
	F.clear();
	F.seekg(0);

	temperatures.resize(nofReplicas);

	rexSamplers.resize(nofReplicas);
	rexDistortOptions.resize(nofReplicas);
	rexFlowOptions.resize(nofReplicas);
	rexWorkOptions.resize(nofReplicas);	
	rexIntegrators.resize(nofReplicas);

	rexTimesteps.resize(nofReplicas);
	rexWorldIndexes.resize(nofReplicas);
	rexMdsteps.resize(nofReplicas);
	rexSamplesPerRound.resize(nofReplicas);

	// Load simulation parameters vectors
	while (F) {

		// Read line
		std::getline(F, line);

		if(line.length() > 0){
			std::vector<std::string> words = split(line, " ");

			if((words[0][0] != '#') && (words[0] != "NOF_REPLICAS")){
				std::cout << "REX FILE READING: ";
				for (const auto& word: words){
					std::cout << word << "|";
				}
				std::cout << std::endl;
				
				// Get replica index
				int repIx = int(std::stoi(words[0]));

				// Get keyword
				switch( rexToIntKeys.at(words[1]) ){
					case RexKey::TEMPERATURE:
						std::cout << "Loaded thermodynamic state " << repIx << " \n" ;
						temperatures[repIx] = std::stod(words[2]);
						break;

					case RexKey::SAMPLERS:
						for(int i = 2; i < words.size(); i++){
							rexSamplers[repIx].push_back(words[i]);
						}
						break;
					case RexKey::DISTORT_OPTIONS:
						for(int i = 2; i < words.size(); i++){
							rexDistortOptions[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::FLOW_OPTIONS:
						for(int i = 2; i < words.size(); i++){
							rexFlowOptions[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::WORK_OPTIONS:
						for(int i = 2; i < words.size(); i++){
							rexWorkOptions[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::INTEGRATORS:
						for(int i = 2; i < words.size(); i++){
							rexIntegrators[repIx].push_back(words[i]);
						}
						break;

					case RexKey::TIMESTEPS:
						for(int i = 2; i < words.size(); i++){
							rexTimesteps[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::WORLD_INDEXES:
						for(int i = 2; i < words.size(); i++){
							rexWorldIndexes[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::MDSTEPS:
						for(int i = 2; i < words.size(); i++){
							rexMdsteps[repIx].push_back(std::stod(words[i]));
						}
						break;
					case RexKey::SAMPLES_PER_ROUND:
						for(int i = 2; i < words.size(); i++){
							rexSamplesPerRound[repIx].push_back(std::stod(words[i]));
						}
						break;

					default:
						;
						break;
				} // keys

			} // not a commentary

		} // line is empty

	} // file

	std::cout << "All the timesteps that I got:\n" ;
	for(auto& timesteps : rexTimesteps){for(auto& ts : timesteps){std::cout << ts << " ";}std::cout << std::endl;}
	std::cout << "All the worldIndexes that I got:\n" ;
	for(auto& worldIndexes : rexWorldIndexes){for(auto& ts : worldIndexes){std::cout << ts << " ";}std::cout << std::endl;}
	std::cout << "All the mdsteps that I got:\n" ;
	for(auto& mdsteps : rexMdsteps){for(auto& ts : mdsteps){std::cout << ts << " ";}std::cout << std::endl;}
	std::cout << "All the samplesPerRound that I got:\n" ;
	for(auto& samplesPerRound : rexSamplesPerRound){for(auto& ts : samplesPerRound){std::cout << ts << " ";}std::cout << std::endl;}	
	std::cout << "All the Samplers that I got:\n" ;
	for(auto& samplers : rexSamplers){for(auto& ts : samplers){std::cout << ts << " ";}std::cout << std::endl;}	
	std::cout << "All the temperatures that I got:\n" ;
	for(auto& temperatures_ix : temperatures){{std::cout << temperatures_ix << " ";}std::cout << std::endl;}	

	return nofReplicas;
	
}


