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


void testAlanineDipeptide(int seed) {

	// auto ommvv = [](int seed) {
	// 	// no need for tini tfin
	// 	Context c("ala2", seed, 0, 1, RUN_TYPE::REMC, 1, 0);

	// 	// c.setNonbonded(0, 1.2); // default - set for each world
	// 	c.setPdbRestartFreq(0); // WRITE_PDBS
	// 	c.setPrintFreq(1); // PRINT_FREQ
	// 	c.setNonbonded(0, 1.2);
	// 	c.setGBSA(1);
	// 	c.setVerbose(false);

	// 	c.loadAmberSystem("alanine-dipeptide/alanine-dipeptide.prmtop", "alanine-dipeptide/alanine-dipeptide.rst7");

	// 	// Fully flexible for OpenMM
	// 	std::vector<BOND_FLEXIBILITY> allFlexibilities = {
	// 		{ 0, 1, BondMobility::Mobility::Translation },
	// 		{ 1, 2, BondMobility::Mobility::Translation },
	// 		{ 1, 3, BondMobility::Mobility::Translation },
	// 		{ 1, 4, BondMobility::Mobility::Translation },
	// 		{ 4, 5, BondMobility::Mobility::Translation },
	// 		{ 4, 6, BondMobility::Mobility::Translation },
	// 		{ 6, 7, BondMobility::Mobility::Translation },
	// 		{ 6, 8, BondMobility::Mobility::Translation },
	// 		{ 8, 9, BondMobility::Mobility::Translation },
	// 		{ 8, 10, BondMobility::Mobility::Translation },
	// 		{ 8, 14, BondMobility::Mobility::Translation },
	// 		{ 10, 11, BondMobility::Mobility::Translation },
	// 		{ 10, 12, BondMobility::Mobility::Translation },
	// 		{ 10, 13, BondMobility::Mobility::Translation },
	// 		{ 14, 15, BondMobility::Mobility::Translation },
	// 		{ 14, 16, BondMobility::Mobility::Translation },
	// 		{ 16, 17, BondMobility::Mobility::Translation },
	// 		{ 16, 18, BondMobility::Mobility::Translation },
	// 		{ 18, 19, BondMobility::Mobility::Translation },
	// 		{ 18, 20, BondMobility::Mobility::Translation },
	// 		{ 18, 21, BondMobility::Mobility::Translation },
	// 	};
	// 	c.addWorld(false, 1, ROOT_MOBILITY::WELD, allFlexibilities);

	// 	// Add samplers
	// 	c.getWorld(0).addSampler(SamplerName::HMC, IntegratorType::OMMVV, ThermostatName::ANDERSEN, false);

	// 	// IntegratorType::VERLET
	// 	std::vector<SimTK::Real> temperatures = { 300 }, boostTemperatures = { 300 };
	// 	std::vector<AcceptRejectMode> acceptRejectModes = { AcceptRejectMode::AlwaysAccept };
	// 	std::vector<IntegratorType> integrators = { IntegratorType::OMMVV };
	// 	std::vector<SimTK::Real> timesteps = { 0 };
	// 	std::vector<int> worldIndexes = { 0 };
	// 	std::vector<int> mdsteps = { 10 };
	// 	std::vector<int> distortOptions = { 0 };
	// 	std::vector<std::string> distortArgs = { "0" };
	// 	std::vector<int> flow = { 0 };
	// 	std::vector<int> work = { 0 };

	// 	c.addReplica(0);
	// 	c.addThermodynamicState(0, temperatures[0], acceptRejectModes, distortOptions, distortArgs, flow, work, integrators, worldIndexes, timesteps, mdsteps);

	// 	c.Initialize();

	// 	c.PrintSimbodyMobods();

	// 	// Run the simulation
	// 	c.RunREX(0, 1);
	// };

	// ommvv(42);
	// return;


	// no need for tini tfin
	Context c("ala2", seed, 0, 1, RUN_TYPE::REMC, 1, 0);

	// c.setNonbonded(0, 1.2); // default - set for each world
	c.setPdbRestartFreq(0); // WRITE_PDBS
	c.setPrintFreq(50); // PRINT_FREQ
	c.setNonbonded(0, 1.2);
	c.setGBSA(1);
	c.setVerbose(true);

	c.loadAmberSystem("alanine-dipeptide/alanine-dipeptide.prmtop", "alanine-dipeptide/alanine-dipeptide.rst7");

	// Fully flexible for OpenMM
	std::vector<BOND_FLEXIBILITY> allFlexibilities = {
		{ 0, 1, BondMobility::Mobility::Translation },
		{ 1, 2, BondMobility::Mobility::Translation },
		{ 1, 3, BondMobility::Mobility::Translation },
		{ 1, 4, BondMobility::Mobility::Translation },
		{ 4, 5, BondMobility::Mobility::Translation },
		{ 4, 6, BondMobility::Mobility::Translation },
		{ 6, 7, BondMobility::Mobility::Translation },
		{ 6, 8, BondMobility::Mobility::Translation },
		{ 8, 9, BondMobility::Mobility::Translation },
		{ 8, 10, BondMobility::Mobility::Translation },
		{ 8, 14, BondMobility::Mobility::Translation },
		{ 10, 11, BondMobility::Mobility::Translation },
		{ 10, 12, BondMobility::Mobility::Translation },
		{ 10, 13, BondMobility::Mobility::Translation },
		{ 14, 15, BondMobility::Mobility::Translation },
		{ 14, 16, BondMobility::Mobility::Translation },
		{ 16, 17, BondMobility::Mobility::Translation },
		{ 16, 18, BondMobility::Mobility::Translation },
		{ 18, 19, BondMobility::Mobility::Translation },
		{ 18, 20, BondMobility::Mobility::Translation },
		{ 18, 21, BondMobility::Mobility::Translation },
	};
	c.addWorld(false, 1, ROOT_MOBILITY::WELD, allFlexibilities);

	// Sidechain fleibilities
	std::vector<BOND_FLEXIBILITY> sidechainFlexibilities = {
		{ 8, 10, BondMobility::Mobility::Torsion },
	};
	c.addWorld(true, 1, ROOT_MOBILITY::WELD, sidechainFlexibilities);

	// Ramachandran fleibilities
	std::vector<BOND_FLEXIBILITY> ramachandranFlexibilities = {
		{ 6, 8, BondMobility::Mobility::Torsion },
		{ 8, 14, BondMobility::Mobility::Torsion },
	};
	c.addWorld(true, 1, ROOT_MOBILITY::WELD, ramachandranFlexibilities);

	// Add samplers
	c.getWorld(0).addSampler(SamplerName::HMC, IntegratorType::OMMVV, ThermostatName::ANDERSEN, false);
	c.getWorld(1).addSampler(SamplerName::HMC, IntegratorType::VERLET, ThermostatName::ANDERSEN, true);
	c.getWorld(2).addSampler(SamplerName::HMC, IntegratorType::VERLET, ThermostatName::ANDERSEN, true);

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
	std::vector<IntegratorType> integrators = { IntegratorType::OMMVV, IntegratorType::VERLET, IntegratorType::VERLET };
	std::vector<SimTK::Real> timesteps = { 0.0007, 0.035, 0.07 };
	std::vector<int> worldIndexes = { 0, 1, 2 };
	std::vector<int> mdsteps = { 10, 10, 10 };
	// std::vector<int> boostMDSteps = mdsteps;
	// std::vector<int> samplesPerRound = { 1, 1, 1 };

	std::vector<int> distortOptions = { 0, 0, 0 };
	std::vector<std::string> distortArgs = { "0", "0", "0" };
	std::vector<int> flow = { 0, 0, 0 };
	std::vector<int> work = { 0, 0, 0 };

	for (int i = 0; i < num_replicas; i++) {
		c.addReplica(i);
		c.addThermodynamicState(i,
			temperatures[i],
			acceptRejectModes,
			distortOptions,
			distortArgs,
			flow,
			work,
			integrators,
			worldIndexes,
			timesteps,
			mdsteps
		);
	}

	c.Initialize();

	c.PrintSimbodyMobods();

	// Run the simulation
	c.RunREX(0, 100);
}

/*!
 * <!-- Victor's Context test -->
*/
int testContext(int seed)
{
	// no need for tini tfin
	Context c("2but", seed, 0, 1, RUN_TYPE::REMC, 1, 0);

	// c.setNonbonded(0, 1.2); // default - set for each world
	c.setPdbPrefix("rage"); // will disappear
	c.setOutput("temp"); // the log file is created like log.[seed] - needs refactoring
	c.setRequiredNofRounds(10); // per world? what does it do?
	c.setPdbRestartFreq(0); // WRITE_PDBS
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

	// // World 2
	// std::vector<BOND_FLEXIBILITY> flexibilities_w2 = {
	// };
	// c.addWorld(true, 1, ROOT_MOBILITY::WELD, flexibilities_w2);

	// Does OMMVV have to be first?
	// Add samplers
	c.getWorld(0).addSampler(SamplerName::HMC, IntegratorType::OMMVV, ThermostatName::ANDERSEN, false);

	// // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	c.getWorld(1).addSampler(SamplerName::HMC, IntegratorType::VERLET, ThermostatName::ANDERSEN, true);

	// // // Does not work if I set OMMVV instead of VV. How do I check if it is working?
	// c.getWorld(2).addSampler(SamplerName::HMC, IntegratorType::Verlet, ThermostatName::ANDERSEN, true);

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
	std::vector<IntegratorType> integrators = { IntegratorType::OMMVV, IntegratorType::VERLET };
	std::vector<SimTK::Real> timesteps = { 0.0007, 0.002, 0.002 };
	std::vector<int> worldIndexes = { 0, 1, 2 };
	std::vector<int> mdsteps = { 25, 100, 50 };
	std::vector<int> boostMDSteps = { 25, 100, 50 };
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

	// acceptRejectModes.pop_back();
	// timesteps.pop_back();
	// worldIndexes.pop_back();
	// mdsteps.pop_back();
	// boostMDSteps.pop_back();
	// samplesPerRound.pop_back();
	// distortOptions.pop_back();
	// distortArgs.pop_back();
	// flow.pop_back();
	// work.pop_back();

	for (int i = 0; i < num_replicas; i++) {
		c.addReplica(i);
		c.addThermodynamicState(i,
			temperatures[i],
			acceptRejectModes,
			distortOptions,
			distortArgs,
			flow,
			work,
			integrators,
			worldIndexes,
			timesteps,
			mdsteps
		);
	}

	c.Initialize();

	// pas how many rounds to run for here
	// c.RunREXNew(0, 1'000'000);
	c.RunREX(10, 10);

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

/*! <!-- Main --> */
int main(int argc, char **argv)
{
    // if (argc < 10) {
    //     std::cerr << "Usage: " << argv[0] << " <name> <top> <rst7> <equilRounds> <prodRounds> <writeFreq> <baseTemperature> <runType> <seed>" << std::endl;
    //     return -1;
    // }	
	std::string name = (argv[1]); std::cout << "name: " << name << std::endl;
	std::string top = (argv[2]); std::cout << "top: " << top << std::endl;
	std::string rst7 = (argv[3]); std::cout << "rst7: " << rst7 << std::endl;
	int equilRounds = std::stoi(argv[4]); std::cout << "equilRounds: " << equilRounds << std::endl;
	int prodRounds = std::stoi(argv[5]); std::cout << "prodRounds: " << prodRounds << std::endl;
	int writeFreq = std::stoi(argv[6]); std::cout << "writeFreq: " << writeFreq << std::endl;
	float baseTemperature = std::stof(argv[7]); std::cout << "baseTemperature: " << baseTemperature << std::endl;
	
	std::string runTypeStr = argv[8];
	if (RUN_TYPE_MAP.find(runTypeStr) == RUN_TYPE_MAP.end()) {
        std::cerr << "Error: Invalid run type '" << runTypeStr 
			<< "'. Valid options are: DEFAULT, REMC, RENEMC, RENE." << std::endl;
        return -1;
    }
    RUN_TYPE runType = RUN_TYPE_MAP.at(runTypeStr); std::cout << "runType: " << runTypeStr << std::endl;

    #pragma region checkArgs //Check if the necessary files exist
	if (!fileExists(top)) {
		std::cerr << "Error: File does not exist: " << top << std::endl;
		return -1;
	}
	if (!fileExists(rst7)) {
		std::cerr << "Error: File does not exist: " << rst7 << std::endl;
		return -1;
	}
	#pragma endregion // checkArgs

	int seed = std::stoi(argv[9]); std::cout << "seed: " << seed << std::endl;

	int threads = 0;
	uint32_t nofRoundsTillReblock = -1; // not implemented yet
	int swapFreq = 1;
	int swapFixmanFreq = 0;

	Context context(name, seed, threads, nofRoundsTillReblock, runType, swapFreq, swapFixmanFreq);

	// Setup general input-output parameters
	context.setOutputDir("temp/");
	context.setPdbPrefix(name);

	context.setPdbRestartFreq(1);
	context.setPrintFreq(writeFreq);

	// Set nonbonded method and cutoff
	context.setNonbonded(0, 1.2);
	context.setGBSA(0.0);
	context.setAtomMasses();
	context.setVerbose(true);

	context.loadAmberSystem(top, rst7);

	// context.calcZMatrixTable();
	// context.PrintZMatrixTable();
	// context.reallocZMatrixBAT();

	// Fully flexible world
	std::string flexFileFN = "../tests_inputs/" + name + "/ligand.cart.flex";
	if (!fileExists(flexFileFN)) {
		std::cerr << "Error: File does not exist: " << flexFileFN << std::endl;
		return -1;
	}	
	std::vector<BOND_FLEXIBILITY> flexibilities = {};
	context.readFlexibility(flexFileFN, flexibilities);
	context.addWorld(false, 1, ROOT_MOBILITY::CARTESIAN, flexibilities, true, false, 0);

	SamplerName sampler = SamplerName::HMC;
	ThermostatName thermostat = ThermostatName::ANDERSEN;
	IntegratorType integrator = IntegratorType::OMMVV;

	context.getWorld(0).addSampler(sampler, integrator, thermostat, false);

	std::vector<SimTK::Real> temperatures = { baseTemperature };
	std::vector<SimTK::Real> boostTemperatures = { baseTemperature };
	std::vector<AcceptRejectMode> accept_reject_modes = { AcceptRejectMode::MetropolisHastings };
	std::vector<int> distort_options = { 0 };
	std::vector<std::string> distort_args = { "0" };
	std::vector<int> flow = { 0 };
	std::vector<int> work = { 0 };
	std::vector<IntegratorType> integrators = { integrator };
	std::vector<int> worldIndexes = { 0 };
	std::vector<double> timesteps = { 0.007 };
	std::vector<int> mdsteps = { 10 };
	
	// Restart directory
	std::string restartDir = "rest." + name + "." + std::to_string(seed);
	context.setRestartDir(restartDir);
	if (!fileExists(top)) {
		std::cerr << "Error: File does not exist: " << restartDir << std::endl;
		return -1;
	}

	size_t nofReplicas = 1;
	context.addReplicasAndLoadCoordinates(name, top, restartDir, nofReplicas);

	for(int replIx = 0; replIx < nofReplicas; replIx++){

		std::cout << "Adding thermo state " << replIx << std::endl;

		context.addThermodynamicState(replIx,
			300.00,
			accept_reject_modes,
			distort_options,
			distort_args,
			flow,
			work,
			integrators,
			worldIndexes,
			timesteps,
			mdsteps);
	}

	// Initialize
	context.Initialize();

	context.setSwapEvery(1);
	context.setSwapFixman(0);

	// -- Run --
	context.RunREX(equilRounds, prodRounds);
		
	std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	
	// Write final pdbs
	//c.writeFinalPdb();
	//c.printStatus();

	return 0;	

}
