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
int testContext(int seed)
{
	// no need for tini tfin
	Context c("2but", seed, 0, 1, RUN_TYPE::REMC, 1, 0);

	// c.setNonbonded(0, 1.2); // default - set for each world
	c.setPdbPrefix("2but"); // will disappear
	c.setOutput("temp"); // the log file is created like log.[seed] - needs refactoring
	c.setRequiredNofRounds(10); // per world? what does it do?
	c.setPdbRestartFreq(1); // WRITE_PDBS
	c.setPrintFreq(1); // PRINT_FREQ

	// read files and create topologies. this populates ```atoms``` and ```bonds```
	// atoms must have a dumm index which is used by bonds, angles and torsions
	// generate topologies
	// addWorld will deal with adding generating dumm parameters
	// c.loadAmberSystem("diala_double/diala_double.prmtop", "diala_double/diala_double.rst7");
	c.loadAmberSystem("2but/ligand.prmtop", "2but/ligand.rst7");
	// c.loadAmberSystem("2kyy/2kyy.H.capped.prmtop", "2kyy/2kyy.H.capped.rst7");
	// c.loadAmberSystem("ala10/ligand.prmtop", "ala10/ligand.rst7");
	// c.loadAmberSystem("rage.fps.0.prmtop", "rage.fps.0.rst7");
	// c.loadAmberSystem("cb8/ligand.prmtop", "cb8/ligand.rst7");

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
		{ 4, 13, BondMobility::Mobility::Translation },

		// {0, 1, BondMobility::Mobility::Translation},
		// {1, 2, BondMobility::Mobility::Translation},
		// {1, 3, BondMobility::Mobility::Translation},
		// {1, 4, BondMobility::Mobility::Translation},
		// {4, 5, BondMobility::Mobility::Translation},
		// {4, 6, BondMobility::Mobility::Translation},
		// {6, 7, BondMobility::Mobility::Translation},
		// {6, 8, BondMobility::Mobility::Translation},
		// {8, 9, BondMobility::Mobility::Translation},
		// {8, 10, BondMobility::Mobility::Translation},
		// {8, 14, BondMobility::Mobility::Translation},
		// {10, 11, BondMobility::Mobility::Translation},
		// {10, 12, BondMobility::Mobility::Translation},
		// {10, 13, BondMobility::Mobility::Translation},
		// {14, 15, BondMobility::Mobility::Translation},
		// {14, 16, BondMobility::Mobility::Translation},
		// {16, 17, BondMobility::Mobility::Translation},
		// {16, 18, BondMobility::Mobility::Translation},
		// {18, 19, BondMobility::Mobility::Translation},
		// {18, 20, BondMobility::Mobility::Translation},
		// {18, 24, BondMobility::Mobility::Translation},
		// {20, 21, BondMobility::Mobility::Translation},
		// {20, 22, BondMobility::Mobility::Translation},
		// {20, 23, BondMobility::Mobility::Translation},
		// {24, 25, BondMobility::Mobility::Translation},
		// {24, 26, BondMobility::Mobility::Translation},
		// {26, 27, BondMobility::Mobility::Translation},
		// {26, 28, BondMobility::Mobility::Translation},
		// {28, 29, BondMobility::Mobility::Translation},
		// {28, 30, BondMobility::Mobility::Translation},
		// {28, 34, BondMobility::Mobility::Translation},
		// {30, 31, BondMobility::Mobility::Translation},
		// {30, 32, BondMobility::Mobility::Translation},
		// {30, 33, BondMobility::Mobility::Translation},
		// {34, 35, BondMobility::Mobility::Translation},
		// {34, 36, BondMobility::Mobility::Translation},
		// {36, 37, BondMobility::Mobility::Translation},
		// {36, 38, BondMobility::Mobility::Translation},
		// {38, 39, BondMobility::Mobility::Translation},
		// {38, 40, BondMobility::Mobility::Translation},
		// {38, 44, BondMobility::Mobility::Translation},
		// {40, 41, BondMobility::Mobility::Translation},
		// {40, 42, BondMobility::Mobility::Translation},
		// {40, 43, BondMobility::Mobility::Translation},
		// {44, 45, BondMobility::Mobility::Translation},
		// {44, 46, BondMobility::Mobility::Translation},
		// {46, 47, BondMobility::Mobility::Translation},
		// {46, 48, BondMobility::Mobility::Translation},
		// {48, 49, BondMobility::Mobility::Translation},
		// {48, 50, BondMobility::Mobility::Translation},
		// {48, 54, BondMobility::Mobility::Translation},
		// {50, 51, BondMobility::Mobility::Translation},
		// {50, 52, BondMobility::Mobility::Translation},
		// {50, 53, BondMobility::Mobility::Translation},
		// {54, 55, BondMobility::Mobility::Translation},
		// {54, 56, BondMobility::Mobility::Translation},
		// {56, 57, BondMobility::Mobility::Translation},
		// {56, 58, BondMobility::Mobility::Translation},
		// {58, 59, BondMobility::Mobility::Translation},
		// {58, 60, BondMobility::Mobility::Translation},
		// {58, 64, BondMobility::Mobility::Translation},
		// {60, 61, BondMobility::Mobility::Translation},
		// {60, 62, BondMobility::Mobility::Translation},
		// {60, 63, BondMobility::Mobility::Translation},
		// {64, 65, BondMobility::Mobility::Translation},
		// {64, 66, BondMobility::Mobility::Translation},
		// {66, 67, BondMobility::Mobility::Translation},
		// {66, 68, BondMobility::Mobility::Translation},
		// {68, 69, BondMobility::Mobility::Translation},
		// {68, 70, BondMobility::Mobility::Translation},
		// {68, 74, BondMobility::Mobility::Translation},
		// {70, 71, BondMobility::Mobility::Translation},
		// {70, 72, BondMobility::Mobility::Translation},
		// {70, 73, BondMobility::Mobility::Translation},
		// {74, 75, BondMobility::Mobility::Translation},
		// {74, 76, BondMobility::Mobility::Translation},
		// {76, 77, BondMobility::Mobility::Translation},
		// {76, 78, BondMobility::Mobility::Translation},
		// {78, 79, BondMobility::Mobility::Translation},
		// {78, 80, BondMobility::Mobility::Translation},
		// {78, 84, BondMobility::Mobility::Translation},
		// {80, 81, BondMobility::Mobility::Translation},
		// {80, 82, BondMobility::Mobility::Translation},
		// {80, 83, BondMobility::Mobility::Translation},
		// {84, 85, BondMobility::Mobility::Translation},
		// {84, 86, BondMobility::Mobility::Translation},
		// {86, 87, BondMobility::Mobility::Translation},
		// {86, 88, BondMobility::Mobility::Translation},
		// {88, 89, BondMobility::Mobility::Translation},
		// {88, 90, BondMobility::Mobility::Translation},
		// {88, 94, BondMobility::Mobility::Translation},
		// {90, 91, BondMobility::Mobility::Translation},
		// {90, 92, BondMobility::Mobility::Translation},
		// {90, 93, BondMobility::Mobility::Translation},
		// {94, 95, BondMobility::Mobility::Translation},
		// {94, 96, BondMobility::Mobility::Translation},
		// {96, 97, BondMobility::Mobility::Translation},
		// {96, 98, BondMobility::Mobility::Translation},
		// {98, 99, BondMobility::Mobility::Translation},
		// {98, 100, BondMobility::Mobility::Translation},
		// {98, 104, BondMobility::Mobility::Translation},
		// {100, 101, BondMobility::Mobility::Translation},
		// {100, 102, BondMobility::Mobility::Translation},
		// {100, 103, BondMobility::Mobility::Translation},
		// {104, 105, BondMobility::Mobility::Translation},
		// {104, 106, BondMobility::Mobility::Translation},
		// {106, 107, BondMobility::Mobility::Translation},
		// {106, 108, BondMobility::Mobility::Translation},
		// {108, 109, BondMobility::Mobility::Translation},
		// {108, 110, BondMobility::Mobility::Translation},
		// {108, 111, BondMobility::Mobility::Translation},
	};
	c.addWorld(false, 1, ROOT_MOBILITY::WELD, flexibilities_w0);

	// // World 1
	// std::vector<BOND_FLEXIBILITY> flexibilities_w1 = {
	// 	{ 3, 1, BondMobility::Mobility::Torsion },
	// 	{ 1, 2, BondMobility::Mobility::Torsion },
	// 	{ 2, 4, BondMobility::Mobility::Torsion },
	// 	{ 1, 0, BondMobility::Mobility::Torsion },

	// 	// {6, 8, BondMobility::Mobility::Torsion},
	// 	// {8, 10, BondMobility::Mobility::Torsion},
	// 	// {8, 14, BondMobility::Mobility::Torsion},
	// 	// {16, 18, BondMobility::Mobility::Torsion},
	// 	// {18, 20, BondMobility::Mobility::Torsion},
	// 	// {18, 24, BondMobility::Mobility::Torsion},
	// 	// {26, 28, BondMobility::Mobility::Torsion},
	// 	// {28, 30, BondMobility::Mobility::Torsion},
	// 	// {28, 34, BondMobility::Mobility::Torsion},
	// 	// {36, 38, BondMobility::Mobility::Torsion},
	// 	// {38, 40, BondMobility::Mobility::Torsion},
	// 	// {38, 44, BondMobility::Mobility::Torsion},
	// 	// {46, 48, BondMobility::Mobility::Torsion},
	// 	// {48, 50, BondMobility::Mobility::Torsion},
	// 	// {48, 54, BondMobility::Mobility::Torsion},
	// 	// {56, 58, BondMobility::Mobility::Torsion},
	// 	// {58, 60, BondMobility::Mobility::Torsion},
	// 	// {58, 64, BondMobility::Mobility::Torsion},
	// 	// {66, 68, BondMobility::Mobility::Torsion},
	// 	// {68, 70, BondMobility::Mobility::Torsion},
	// 	// {68, 74, BondMobility::Mobility::Torsion},
	// 	// {76, 78, BondMobility::Mobility::Torsion},
	// 	// {78, 80, BondMobility::Mobility::Torsion},
	// 	// {78, 84, BondMobility::Mobility::Torsion},
	// 	// {86, 88, BondMobility::Mobility::Torsion},
	// 	// {88, 90, BondMobility::Mobility::Torsion},
	// 	// {88, 94, BondMobility::Mobility::Torsion},
	// 	// {96, 98, BondMobility::Mobility::Torsion},
	// 	// {98, 100, BondMobility::Mobility::Torsion},
	// 	// {98, 104, BondMobility::Mobility::Torsion}
	// };
	// c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w1);

	// // World 2
	// std::vector<BOND_FLEXIBILITY> flexibilities_w2 = {
	// 	{6, 8, BondMobility::Mobility::Torsion},
	// 	{8, 14, BondMobility::Mobility::Torsion},
	// 	{16, 18, BondMobility::Mobility::Torsion},
	// 	{18, 24, BondMobility::Mobility::Torsion},
	// 	{26, 28, BondMobility::Mobility::Torsion},
	// 	{28, 34, BondMobility::Mobility::Torsion},
	// 	{36, 38, BondMobility::Mobility::Torsion},
	// 	{38, 44, BondMobility::Mobility::Torsion},
	// 	{46, 48, BondMobility::Mobility::Torsion},
	// 	{48, 54, BondMobility::Mobility::Torsion},
	// 	{56, 58, BondMobility::Mobility::Torsion},
	// 	{58, 64, BondMobility::Mobility::Torsion},
	// 	{66, 68, BondMobility::Mobility::Torsion},
	// 	{68, 74, BondMobility::Mobility::Torsion},
	// 	{76, 78, BondMobility::Mobility::Torsion},
	// 	{78, 84, BondMobility::Mobility::Torsion},
	// 	{86, 88, BondMobility::Mobility::Torsion},
	// 	{88, 94, BondMobility::Mobility::Torsion},
	// 	{96, 98, BondMobility::Mobility::Torsion},
	// 	{98, 104, BondMobility::Mobility::Torsion},
	// };
	// c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w2);

	// Does OMMVV have to be first?
	// Add samplers
	c.getWorld(0).addSampler(SamplerName::HMC, IntegratorName::OMMVV, ThermostatName::ANDERSEN, false);

	// // // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	// c.getWorld(1).addSampler(SamplerName::HMC, IntegratorName::Verlet, ThermostatName::ANDERSEN, true);

	// // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	// c.getWorld(2).addSampler(SamplerName::HMC, IntegratorName::Verlet, ThermostatName::ANDERSEN, true);

	const int num_replicas = 1; // Number of replicas
    const double T_min = 300.0;  // Starting temperature (Replica 0)
    const double T_max = 600.0;  // Ending temperature (Replica 19)

    // Calculate the common ratio for geometric progression
    double r = std::pow(T_max / T_min, 1.0 / (num_replicas - 1));

	SimTK::Real temperature = 300;
	std::vector<SimTK::Real> temperatures (num_replicas), boostTemperatures(num_replicas);
	for (int i = 0; i < num_replicas; i++) {
		temperatures[i] = T_min * std::pow(r, i);
		boostTemperatures[i] = T_min * std::pow(r, i); // used for openmm velocities
	}

	std::vector<AcceptRejectMode> acceptRejectModes = { AcceptRejectMode::MetropolisHastings, AcceptRejectMode::MetropolisHastings, AcceptRejectMode::MetropolisHastings };
	std::vector<SimTK::Real> timesteps = { 100000000000000, 0.006, 0.008 };
	std::vector<int> worldIndexes = { 0, 1, 2 };
	std::vector<int> mdsteps = { 50, 50, 50 };
	std::vector<int> boostMDSteps = { 50, 50, 50 };
	std::vector<int> samplesPerRound = { 1, 1, 1 };

	std::vector<int> distortOptions = { 0, 0, 0 };
	std::vector<std::string> distortArgs = { "0", "0", "0" };
	std::vector<int> flow = { 0, 0, 0 };
	std::vector<int> work = { 0, 0, 0 };

	acceptRejectModes.pop_back();
	timesteps.pop_back();
	worldIndexes.pop_back();
	mdsteps.pop_back();
	boostMDSteps.pop_back();
	samplesPerRound.pop_back();
	distortOptions.pop_back();
	distortArgs.pop_back();
	flow.pop_back();
	work.pop_back();

	acceptRejectModes.pop_back();
	timesteps.pop_back();
	worldIndexes.pop_back();
	mdsteps.pop_back();
	boostMDSteps.pop_back();
	samplesPerRound.pop_back();
	distortOptions.pop_back();
	distortArgs.pop_back();
	flow.pop_back();
	work.pop_back();

	for (int i = 0; i < num_replicas; i++) {
		c.addReplica(i);
		c.addThermodynamicState(i,
			temperatures[i],
			boostTemperatures[i],
			acceptRejectModes,
			distortOptions,
			distortArgs,
			flow,
			work,
			worldIndexes,
			timesteps,
			mdsteps,
			boostMDSteps
			);
	}

	c.initialize();

	// pas how many rounds to run for here
	c.RunREXNew(0, 10);

	return 0;
}

/*!
 * <!-- Helpers -->
*/
void printInertia(std::string prefix, const SimTK::Inertia inertia, int precision, bool scientific = false)
{
	const SymMat33 inertia33 = inertia.asSymMat33();

	if(scientific){
		std::cout << std::scientific;
	}

	std::cout << prefix << std::fixed << std::setprecision(15) << std::endl
		<<" " << inertia33[0][0] <<" " << inertia33[0][1] <<" " << inertia33[0][2] << std::endl
		<<" " << inertia33[1][0] <<" " << inertia33[1][1] <<" " << inertia33[1][2] << std::endl
		<<" " << inertia33[2][0] <<" " << inertia33[2][1] <<" " << inertia33[2][2] << std::endl
		<< std::endl;	
}

/*!
 * <!-- Test ground -->
*/
void testSOA(void)
{
	// Hydrogen example
	const Real HydrogenRadius = 0.0000012; // 1.2e-6;
	const Real Mass = 1;

	// Get mass properrties for a cluster of atoms
	std::vector<SimTK::Real> myMasses{1.0, 12.0};
	std::vector<SimTK::Real> myRadii{0.0000012, 0.0000027473};
	std::vector<SimTK::Vec3> myStations{SimTK::Vec3(1,1,1), SimTK::Vec3(0,1,0)};
	
	Real    mass = 0;
    Vec3    com(0);
    Inertia inertia(0);
	Inertia inertia_Spheres(0);
	int aIx = -1;
    for(auto ma : myMasses) {
		aIx++;

		// Accumulate mass
        mass    += ma;
        com     += ma * myStations[aIx];

		// Inertia as if the atoms are point masses
		Inertia pointMass = Inertia(myStations[aIx], ma);
		printInertia("inertia", pointMass , 15);
        inertia += pointMass;

		// Inertia as if the atoms are spheres of nucleus radius
		Inertia sphericalInertia = UnitInertia::sphere(myRadii[aIx]);
		printInertia("sphericalInertia unit ", sphericalInertia, 15);

		sphericalInertia *= ma;
		printInertia("sphericalInertia scaled ", sphericalInertia, 15);

		sphericalInertia += pointMass;
		printInertia("sphericalInertia shifted", sphericalInertia, 15);

		inertia_Spheres += sphericalInertia;
		std::cout << "---------" << std::endl;

    }
    com /= mass; // center of mass normalization

	std::cout << "=========" << std::endl;
	printInertia("inertia", inertia, 15, true);
	printInertia("final inertia_Spheres", inertia_Spheres, 15, true);

	//Transform IdentityT;
	//MassProperties massProps_Transformed = MassProperties(mass,com,inertia).calcTransformedMassProps(IdentityT);


}


/*!
 * <!-- Main -->
*/
int main(int argc, char **argv)
{
	int seed = std::stoi(argv[1]);
	testContext(seed);
	return 0;


	// //testSOA();
	// //exit(0);

	// std::cout << "OS memory 0.\n" << exec("free") << std::endl;
	
	// std::string helpString = "Usage: ./robosample [options]\n";
	// helpString +=
	// 	" Options:\n";
	// helpString +=
	// 	"  -h, --help for help\nUsage: ./robosample <file> [dcdfile]\n";

	// if(argc < 2) {
	// 	std::cout << "Error: not enough parameters to run. See help below.\n";
	// 	std::cout << helpString;

	// 	return -1;
	// }
	// else {
	// 	auto arg = std::string(argv[1]);
	// 	if("-h" == arg || "--help" == arg) {
	// 	std::cout << helpString;

	// 		return -1;
	// 	}
	// }

	// // Get the seed and pass it to the context
	// int seed = 0;
	// std::string seedStr = extractValueFromInputFile(argv[1], "SEED");
	// if (!seedStr.empty()) {
	// 	seed = std::stoi(seedStr);
	// }

	// Context context(300, 300, seed, 0, 10, RUN_TYPE::DEFAULT, 1, 0);

	// if (!context.initializeFromFile(argv[1])) {
	// 	return -1;
	// }

	// if(argc >=3){
	// 	context.appendDCDReporter(argv[2]);
	// }

	// // -- Run --
	// context.Run();
		
	// std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	
	// // Write final pdbs
	// //c.writeFinalPdb();
	// //c.printStatus();

	// return 0;
}
