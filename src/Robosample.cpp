#include "Context.hpp"

int main(int argc, char **argv)
{
	//std::cout << "OS memory 0.\n" << exec("free") << std::endl;
	
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
	
	Context c; // c = Context();

	// all are optional. seed is set to random in context constructor
	// c.setNumThreads(3);
	// c.setSeed(42);
	// c.setOutput("temp"); // the log file is created like log.[seed]
	// c.useVisualizer(frequency); // GLOBAL?

	/*
	c = Context()

	for i in range(10):
		c.addWorld(...)
		c.getWorld(i).useOpenMM(...)

	c.loadAmberSystem()

	for w in c.getWorlds():
		w.addSampler(...)
		w.world.getSampler(0).setNonequilibriumParameters(...)
	*/

	bool singlePrmtop = true;

	std::string singlePrmtopOpt;
	if(argc >= 3){
		singlePrmtopOpt = argv[2];
	}

	if(singlePrmtopOpt == "singlePrmtop"){
		singlePrmtop = true;
	}

	if (!c.initializeFromFile(argv[1], singlePrmtop)) {
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
