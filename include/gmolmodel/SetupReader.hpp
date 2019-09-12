#ifndef __SETUPREADER_HPP__
#define __SETUPREADER_HPP__

#include "Robo.hpp"

// This class implements the functionality for reading arguments from a setup
// file and give back the necessary values as vectors
class SetupReader{
public:

    // Constructor
    SetupReader(const char *FN);
    SetupReader(std::string& FN);

    // Destructor
    ~SetupReader();

    // Read setup function
    void ReadSetup(const char *FN);
    void ReadSetup(std::string& FN);

    // Print all the arguments
    void dump(void);

    // Access values by key
    bool find( std::string argKey );

    std::vector<std::string> get(const char *argKey);
    std::vector<std::string> get(std::string argKey);

private:
    std::map<std::string, std::vector<std::string>> Args; // arguments
    std::map<std::string, std::vector<std::string>>::iterator ArgsIt; // arguments
};


#endif //__SETUPREADER_HPP__

