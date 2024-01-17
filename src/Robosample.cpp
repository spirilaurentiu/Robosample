#include "Context.hpp"

int main(int argc, char **argv)
{
	//std::cout << "OS memory 0.\n" << exec("free") << std::endl;
	
	std::string helpString = 
	"Usage: ./robosample [options]\n Options:\n  -h, --help for help\nUsage: ./robosample [file]\n";

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
	
	Context c;
	if (!c.initializeFromFile(argv[1], false)) {
		return -1;
	}

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
