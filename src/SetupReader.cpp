#include "SetupReader.hpp"

// Constructor
SetupReader::SetupReader(const char *cFN)
{
    ReadSetup(cFN);
}

// Constructor
SetupReader::SetupReader(std::string& FN)
{
    ReadSetup(FN);
}

// Destructor
SetupReader::~SetupReader(void)
{
}

// Read setup function
void SetupReader::ReadSetup(const char *cFN)
{
    std::string line;
    std::ifstream F(cFN);
    int line_i = -1;
    while(F.good()){
        line_i++;
        std::getline(F, line);
        std::istringstream iss(line);
        std::string word;
        std::string key;
        std::vector<std::string> V;

        int word_i = -1;
        while(iss >> word){
            if(word[0] == '#'){
                break;
            }
            word_i++;
            if(word_i == 0){
                key = word;
            }else{
                V.push_back(std::move(word));
            }
        }
        ArgsIt = Args.begin();
        Args.insert( ArgsIt, std::pair<std::string, std::vector<std::string>>( std::move(key), std::move(V) ) );
    }
}

void SetupReader::ReadSetup(std::string& FN)
{
    ReadSetup(FN.c_str());
}

// Print all the arguments
void SetupReader::dump(void)
{
    std::vector<std::string>::iterator it;
    for (ArgsIt = Args.begin(); ArgsIt != Args.end(); ++ArgsIt){
        std::cout << ArgsIt->first << " : " ;
        for (it = (ArgsIt->second).begin(); it != (ArgsIt->second).end(); ++it){
            std::cout << *it << " ";
        }
        std::cout << std::endl;
    }
}

/** Check if key exists **/
bool SetupReader::find( std::string argKey )
{

    if( Args.find( argKey ) == Args.end() ){
        return false;
    }else{
        return true;
    }

}

// Access values by key
std::vector<std::string> SetupReader::get(const char *cArgKey)
{
    return Args[std::string(cArgKey)];
}

// Access values by key
std::vector<std::string> SetupReader::get(std::string argKey)
{
    return Args[argKey];
}



