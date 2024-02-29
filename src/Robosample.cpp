#include "Context.hpp"


std::string extractValueFromInputFile(
	const std::string& filename,
	const std::string& keyword)
{
    std::ifstream file(filename);
    std::string line;
    std::string value;

    if (file.is_open()) {

        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string word;

            // Iterate over words in the line
            while (iss >> word) {
                if (word == keyword) {
                    
					// Check if there is another word after the keyword
                    if (iss >> value) {
                        return value;
                    }
                }
            }
        }
        file.close();
    } else {
        std::cerr << "Error: Unable to open file " << filename << std::endl;
    }

    return value;
}


int main(int argc, char **argv)
{
	std::cout << "OS memory 0.\n" << exec("free") << std::endl;
	
	std::string helpString = "Usage: ./robosample [options]\n";
	helpString +=
		" Options:\n";
	helpString +=
		"  -h, --help for help\nUsage: ./robosample [file]\n";

	if(argc < 2) {
		std::cout << "Error: not enough parameters to run. See help below.\n";
		std::cout << helpString;

		return -1;
	}
	else {
		auto arg = std::string(argv[1]);
		if("-h" == arg || "--help" == arg) {
		std::cout << helpString;

			return -1;
		}
	}

	// Get the seed and pass it to the context
	int seed = 0;
	std::string seedStr = extractValueFromInputFile(argv[1], "SEED");
	if (!seedStr.empty()) {
		seed = std::stoi(seedStr);
	}

	Context c(300, 300, seed); // c = Context();

	bool singlePrmtop = true;

	std::string singlePrmtopOpt;
	if(argc >= 3){
		singlePrmtopOpt = argv[2];
	}

	if(singlePrmtopOpt == "singlePrmtop"){
		singlePrmtop = true;
	}

	// singlePrmtop = true;
	if (!c.initializeFromFile(argv[1], singlePrmtop)) {
		return -1;
	}

	// c.appendDCDReporter("2but.dcd");

	// c.loadAmberSystem("2but/ligand.prmtop", "2but/ligand.inpcrd");

	// c.PrintMolmodelAndDuMMTypes();
	// c.PrintSimbodyMobods();

	c.Run();
		
	//std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	// -- Run --
	
	// Write final pdbs
	//c.writeFinalPdb();

	//std::cout << "OS memory 6\n" << exec("free") << std::endl;
	//std::cout << "printStatus 1 " << std::endl;
	//context.printStatus();

	return 0;
}
