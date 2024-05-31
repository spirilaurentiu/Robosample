#include "Context.hpp"

/*!
 * <!-- Extract a value from a file with KEYWORD VAL0 VAL1 ... format -->
*/
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

/*!
 * <!-- Victor's Context test -->
*/
int testContext(int argc, char **argv)
{
	Context c(300, 300, 42);
	// c.setNonbonded(0, 1.2); // default
	c.setNumThreads(0);
	c.setPdbPrefix("2but42");
	// c.setNonbonded(0, 1.2); // default
	c.setOutput("temp"); // the log file is created like log.[seed]
	c.setNofRoundsTillReblock(10); // per world?
	c.setRequiredNofRounds(1);
	c.setPdbRestartFreq(0); // WRITE_PDBS
	c.setPrintFreq(1); // PRINT_FREQ
	c.setRunType(RUN_TYPE::DEFAULT);
	// read files and create topologies. this populates ```atoms``` and ```bonds```
	// atoms must have a dumm index which is used by bonds, angles and torsions
	// generate topologies
	// addWorld will deal with adding generating dumm parameters
	c.loadAmberSystem("2but/ligand.prmtop", "2but/ligand.rst7");
	// World 0 OPENMM
	std::vector<BOND_FLEXIBILITY> flexibilities_w0 = {
		{ 3, 1, BondMobility::Mobility::Translation },
		{ 1, 2, BondMobility::Mobility::Translation },
		{ 2, 4, BondMobility::Mobility::Translation },
		{ 1, 0, BondMobility::Mobility::Translation },
		{ 3, 8, BondMobility::Mobility::Translation },
		{ 3, 9, BondMobility::Mobility::Translation },
		{ 3, 10, BondMobility::Mobility::Translation },
		{ 0, 14, BondMobility::Mobility::Translation },
		{ 1, 5, BondMobility::Mobility::Translation },
		{ 2, 6, BondMobility::Mobility::Translation },
		{ 2, 7, BondMobility::Mobility::Translation },
		{ 4, 11, BondMobility::Mobility::Translation },
		{ 4, 12, BondMobility::Mobility::Translation },
		{ 4, 13, BondMobility::Mobility::Translation }
	};
	c.addWorld(false, 1, ROOT_MOBILITY::WELD, flexibilities_w0);
	// World 1
	std::vector<BOND_FLEXIBILITY> flexibilities_w1 = {
		{ 3, 1, BondMobility::Mobility::Torsion },
		{ 1, 2, BondMobility::Mobility::Torsion },
		{ 2, 4, BondMobility::Mobility::Torsion },
		{ 1, 0, BondMobility::Mobility::Torsion },
	};
	c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w1);
	// Does OMMVV have to be first?
	// Add samplers
	c.getWorld(0).addSampler(SamplerName::HMC, SampleGenerator::MC, IntegratorName::OMMVV, ThermostatName::ANDERSEN, 0.001, 200, 0, 300, 1, 0, 0, 0, false);
	// // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	c.getWorld(1).addSampler(SamplerName::HMC, SampleGenerator::MC, IntegratorName::VERLET, ThermostatName::ANDERSEN, 0.001, 50, 0, 300, 1, 0, 0, 0, true);
	c.appendDCDReporter("2but.dcd");
	c.Run();
	return 0;
	/////////////////////////

}

/*!
 * <!-- Main -->
*/
int main(int argc, char **argv)
{

	std::cout << "OS memory 0.\n" << exec("free") << std::endl;
	
	std::string helpString = "Usage: ./robosample [options]\n";
	helpString +=
		" Options:\n";
	helpString +=
		"  -h, --help for help\nUsage: ./robosample <file> [dcdfile]\n";

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

	Context context(300, 300, seed);

	if (!context.initializeFromFile(argv[1])) {
		return -1;
	}

	if(argc >=3){
		context.appendDCDReporter(argv[2]);
	}

	// -- Run --
	context.Run();
		
	std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	
	// Write final pdbs
	//c.writeFinalPdb();
	//c.printStatus();

	return 0;
}
