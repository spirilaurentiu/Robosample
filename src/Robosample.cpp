/** @file
 * This file tests the Context class
 */

#include <string>
#include <iostream>
#include <sstream>
#include <chrono>
#include <cassert>
#include <sys/stat.h>
#include "simmain.hpp"
#include "Robo.hpp"

#include "HMCSampler.hpp"
#include "readAmberInput.hpp"
#include "SetupReader.hpp"

//#include <filesystem>

//#ifndef ROBO_DEBUG_LEVEL01
//#define ROBO_DEBUG_LEVEL01
//#endif

#include <sys/sysinfo.h>

bool LoadInputIntoSetupReader(int argc, char **argv,
	SetupReader& setupReader)
{

	std::cout << "Reading input...\n" ;

	std::string helpString =
		"Usage: Robsample [options]\n Options:\n  -h, --help for help\nUsage: Robsample file\n";

	if(argc < 2) {
		std::cout << "Error: not enough parameters to run. See help below.\n";
		std::cout << helpString;

		return false;
	}
	else {
		auto arg = std::string(argv[1]);
		if("-h" == arg || "--help" == arg) {
		std::cout << helpString;

			return false;
		}
	}

	// Read input
	setupReader.ReadSetup(argv[1]);
	setupReader.dump(true);
	std::cout << "Done.\n" ;
	return true;

}

bool CreateOutputDirectory(std::string outDir)
{
    if( !SimTK::Pathname::fileExists(outDir + "/pdbs") ){
//	if( ! std::filesystem::is_directory( std::filesystem::status ( outDir + "/pdbs") ) ){
		const int err = mkdir((outDir
			+ "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (err == -1){
			std::cout << "Error creating " << outDir + "/pdbs" << std::endl;
			return false;
		}
	}
    else{
        return true;
    }
}

std::string CreateLogfilename( std::string outDir, long long int seed )
{
	std::string logFilename = outDir + std::string("/log.") + std::to_string(seed);

	return logFilename;

}

std::string GetMoleculeDirectoryShort(std::string path)
{
	std::size_t lastSlashPos = path.find_last_of("/");
	//std::string upDir = path.substr(0, lastSlashPos);
	std::string molDir = path.substr(lastSlashPos + 1, path.length());
	return molDir;
}


vector<string> split(const string& i_str, const string& i_delim)
{
    vector<string> result;

    size_t found = i_str.find(i_delim);
    size_t startIndex = 0;

    while(found != string::npos)
    {
        result.push_back(string(i_str.begin()+startIndex, i_str.begin()+found));
        startIndex = found + i_delim.size();
        found = i_str.find(i_delim, startIndex);
    }
    if(startIndex != i_str.size())
        result.push_back(string(i_str.begin()+startIndex, i_str.end()));
    return result;
}

///////////////////////////////////////////////////////////////////////////////
/////////////////////////////   MAIN   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{	
	//std::cout << "OS memory 0.\n" << exec("free") << std::endl;

	// Read input into a SetupReader object
	SetupReader setupReader;
	if ( !LoadInputIntoSetupReader(argc, argv, setupReader) ){
		return 1;
	}

	// Declare global variables
	int currentWorldIx = 0;
	int round_mcsteps = 0;

	/////////// Create context ////////////

	// - Get molecules directory
	// TODO: move molDir in Context
	std::string molDir = GetMoleculeDirectoryShort(setupReader.get("MOLECULES")[0]);
	std::cout << "Molecule directory: " << molDir << std::endl << std::flush;

	// Create pdbs directory if necessary
	if ( ! CreateOutputDirectory(setupReader.get("OUTPUT_DIR")[0]) ){
		return 1;
	}

	// Set the log filename
	std::string logFilename = CreateLogfilename(
		setupReader.get("OUTPUT_DIR")[0],
		std::stoi(setupReader.get("SEED")[0]));

	// Instantiate a context object
	Context context(setupReader, logFilename);

	// Set the directory where the logs and the trajectories are stored
	context.setOutputDir(setupReader.get("OUTPUT_DIR")[0]);

	context.setPdbPrefix(molDir
		+ setupReader.get("SEED")[0]
		);

	// Alert user of CUDA environment variables
	if(SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT").empty()){
		std::cout << "CUDA_ROOT not set." << std::endl;
	}else{
		std::cout << "CUDA_ROOT set to "
			<< SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
	}

	// Requested nof Worlds in input. We'll try to construct them all
	// but we're not sure if we'll going to succesed.
	int requestedNofWorlds = context.getNofWorlds();
	int requestedNofMols = context.getNofMolecules();
	std::cout << "Requested " << requestedNofMols << " molecules\n" << std::flush;

	/////////// Add Worlds to context ////////////
	// Add Worlds to the Context. Every World instantiates a:
	// CompoundSystem, SimbodyMatterSubsystem, GeneralForceSubsystem,
	// DuMMForceSubsystem, Integrator, TimeStepper and optionally:
	// DecorationSubsystem, Visualizer, VisuzlizerReporter,
	//  ParaMolecularDecorator

	// Deal with visualizer before adding worlds.
	std::vector<double> visualizerFrequencies;
	int i = -1;
	for(auto ts : setupReader.get("TIMESTEPS")){
		i++;
		if (setupReader.get("VISUAL")[i] == "TRUE"){
			visualizerFrequencies.push_back(std::stod(ts));
		}else{
			visualizerFrequencies.push_back(0);
		}
	}

	// Add Worlds
	context.addEmptyWorlds(setupReader.get("WORLDS").size(),
		visualizerFrequencies);

	int nofWorlds = context.getNofWorlds();

	//std::cout << "OS memory 1.\n" << exec("free") << std::endl;
	// Request threads
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		context.setNumThreadsRequested(worldIx,
			std::stoi(setupReader.get("THREADS")[worldIx]));
	}
	context.PrintNumThreads();

	// Set force field scale factor.
	if(setupReader.get("FFSCALE")[0] == "AMBER"){
		context.setAmberForceFieldScaleFactors();
	}else{
		context.setGlobalForceFieldScaleFactor(
			std::stod(setupReader.get("FFSCALE")[0]));
	}

	// Set GBSA scale factor
	context.setGbsaGlobalScaleFactor(
		std::stod(setupReader.get("GBSA")[0]));

	// Use OpenMM if possible
	if(setupReader.get("OPENMM")[0] == "TRUE"){
		context.setUseOpenMMAcceleration(true);
	}

	if(setupReader.get("OPENMM_CalcOnlyNonbonded")[0] == "TRUE"){
		context.setUseOpenMMCalcOnlyNonBonded(true);
	}
	else context.setUseOpenMMCalcOnlyNonBonded(false);

	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		// Only NoCutoff (0) and CutoffNonPeriodic(1) methods are supported. Additional 3 methods available in
		// OpenMM could be integrated as well.
		if(setupReader.get("NONBONDED_METHOD")[worldIx] == "1"){
			context.setNonbondedMethod(worldIx, 1);
			context.setNonbondedCutoff(worldIx,
				std::stod( setupReader.get(
					"NONBONDED_CUTOFF")[worldIx] )
			);
		}

		if(setupReader.get("INTEGRATORS")[worldIx] == "OMMVV"){
			context.setUseOpenMMIntegration(worldIx, std::stod( setupReader.get("BOOST_TEMPERATURE")[worldIx]), std::stod(setupReader.get("TIMESTEPS")[worldIx]));
		}
	}


    // Set Lennard-Jones mixing rule
	context.setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Add molecules based on the setup reader
	context.AddMolecules(requestedNofMols, setupReader);

	//std::cout << "OS memory 2.\n" << exec("free") << std::endl;
	int finalNofMols = context.getNofMolecules();
	if(requestedNofMols != finalNofMols){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}
	std::cout << "Added " << finalNofMols << " molecules" << std::endl;

	// Loads parameters into DuMM
	context.addDummParams(finalNofMols, setupReader);

	// Adopts compound by the CompoundSystem
	// and loads maps of indexes
	context.model(finalNofMols, setupReader);

	// Allocate space for containers that keep statistics if we're doing any
	context.allocWorldsStatsContainers();

	// Adaptive Gibbs blocking
	context.setNofRoundsTillReblock(
		std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]));

	// Only now we can allocate memory for reblocking Q vectors
	//context.allocateReblockQsCacheQVectors();

	// Add Fixman torque (Additional ForceSubsystem) if required
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			context.addFixmanTorque(worldIx);
		}
	}

	// Is this necessary?
	context.realizeTopology();

	// Add empty samplers to the worlds.
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		BaseSampler *p = context.addSampler(worldIx,
			setupReader.get("SAMPLERS")[worldIx],
			setupReader.get("INTEGRATORS")[worldIx]);
	}

	//std::cout << "OS memory 3\n" << exec("free") << std::endl;
	//////////////////////
	// Thermodynamics
	//////////////////////
	// Set thermostats to the samplers
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {

		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			HMCSampler* sampler_p = pHMC(context.updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setThermostat(
				setupReader.get("THERMOSTAT")[worldIx]);
		}

		context.setTemperature(worldIx,
			std::stof(setupReader.get("TEMPERATURE_INI")[0]));
	}

	// Set the guidance Hamiltonian parameters
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			HMCSampler* sampler_p = pHMC(context.updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setBoostTemperature(
				std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]));
		}
	}

	//////////////////////
	// MD parameters
	//////////////////////
	// Set timesteps
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++) {

			// Set timesteps
			context.setTimestep(worldIx, samplerIx,
				std::stof(setupReader.get("TIMESTEPS")[worldIx]));

			// Activate Fixman potential if needed
			if(setupReader.get("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
				 context.useFixmanPotential(worldIx, samplerIx);
			}
		}
	}

	// Set the nunmber of MD steps
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++) {
			context.setNofMDStepsPerSample(worldIx,
				samplerIx,
				std::stoi(setupReader.get("MDSTEPS")[worldIx]));

			HMCSampler* sampler_p = pHMC(context.updWorld(worldIx)->updSampler(samplerIx));

			if(setupReader.find("MDSTEPS_STD")){
				sampler_p->setMDStepsPerSampleStd(
						std::stoi(setupReader.get("MDSTEPS_STD")[worldIx]));
			}else{
				sampler_p->setMDStepsPerSampleStd(0);
			}

		}
	}

	// Guidance Hamiltonian MD
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++) {

			HMCSampler* sampler_p = pHMC(context.updWorld(worldIx)->updSampler(samplerIx));
			sampler_p->setBoostMDSteps(
					std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]));

		}
	}

	//////////////////////
	// Simulation parameters
	//////////////////////
	// Set the seeds for reproducibility. Samplers have to be here already.
	// Let the user set one seed only. TODO: better algorithm (Victor).
	if( setupReader.find("SEED") ){
		if( !(setupReader.get("SEED").empty()) ){
			for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
				context.setSeed(worldIx, 0, std::stoi(setupReader.get("SEED")[0]) + worldIx );
			}
		}
	}

	// Initialize samplers
	/** Set simulation temperature,
	velocities to desired temperature, variables that store the configuration
	and variables that store the energies, both needed for the
	acception-rejection step. Also realize velocities and initialize
	the timestepper. **/
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		for (int samplerIx = 0;
		samplerIx < context.getWorld(worldIx)->getNofSamplers();
		samplerIx++){
			context.initializeSampler(worldIx, samplerIx);
		}
	}

	// Set the number of samples per round
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		context.setNofSamplesPerRound(worldIx,
			std::stoi(setupReader.get("SAMPLES_PER_ROUND")[worldIx]));
	}

	// Simulation parameters
	currentWorldIx = 0;
	round_mcsteps = 0;

	context.setNofRounds(std::stoi(setupReader.get("ROUNDS")[0]));

	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		round_mcsteps += context.getNofSamplesPerRound(worldIx);
	}

	// Set pdb writing frequency
	context.setPdbRestartFreq( std::stoi(setupReader.get("WRITEPDBS")[0]) );

	// Get atom indeces for geometry calculations
	if(setupReader.get("GEOMETRY")[0] == "TRUE"){
		std::vector<int> distanceIx;
		std::vector<int> dihedralIx;
		for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
			distanceIx.push_back(atoi(setupReader.get("DISTANCE")[i].c_str()));
		}
		for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
			dihedralIx.push_back(atoi(setupReader.get("DIHEDRAL")[i].c_str()));
		}
		context.addDistances(distanceIx);
		context.addDihedrals(dihedralIx);
	}

	// Write intial pdb for reference
	if(setupReader.get("WRITEPDBS")[0] == "TRUE"){
		context.writeInitialPdb();
	}

	// Get output printing frequency
	context.setPrintFreq( std::stoi(setupReader.get("PRINT_FREQ")[0]) );

	// Realize topology for all the Worlds
	context.realizeTopology();

	// U Scale Factors uses maps stored in Topology
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
		(context.updWorld(worldIx))->setUScaleFactorsToMobods();
	}

	// Realize topology for all the Worlds
	context.realizeTopology();

	//std::cout << "OS memory 4.\n" << exec("free") << std::endl;
	// Check atom stations for debug purposes
	context.checkAtomStationsThroughDumm();

	// Load/store Mobilized bodies joint types in samplers
	context.loadMbxsToMobilities();

	// -- Setup REX --
	if(setupReader.get("RUN_TYPE")[0] == "REX"){

		SetupReader rexReader;

		// Storage for thermodynamic state temperatures
		std::vector<SimTK::Real> temperatures;

		// Storage for sampling details
		std::vector<std::vector<std::string>> rexSamplers;
		std::vector<std::vector<int>> rexDistortOptions;
		std::vector<std::vector<int>> rexFlowOptions;
		std::vector<std::vector<int>> rexWorkOptions;

		// Storage for each replica simulation parameters
		std::vector<std::vector<std::string>> rexIntegrators;
		std::vector<std::vector<SimTK::Real>> rexTimesteps;
		std::vector<std::vector<int>> rexWorldIndexes;
		std::vector<std::vector<int>> rexMdsteps;
		std::vector<std::vector<int>> rexSamplesPerRound;

		// Read REX config file
		size_t nofReplicas = rexReader.readREXConfigFile(
			setupReader.get("REX_FILE")[0],
			temperatures,

			rexSamplers,
			rexDistortOptions,
			rexFlowOptions,
			rexWorkOptions,
			rexIntegrators,

			rexTimesteps,
			rexWorldIndexes,
			rexMdsteps,
			rexSamplesPerRound);

		// Add thermodynamic states
		for(int k = 0; k < temperatures.size(); k++){
			context.addThermodynamicState(k, temperatures[k],

				rexSamplers[k],
				rexDistortOptions[k],
				rexFlowOptions[k],
				rexWorkOptions[k],
				rexIntegrators[k],

				rexWorldIndexes[k], rexTimesteps[k], rexMdsteps[k]);
		}

		// Add replicas
		for(int k = 0; k < nofReplicas; k++){
			context.addReplica(k);
		}

		// Consider renaming
		context.loadReplica2ThermoIxs();

		context.PrintReplicas();

		// How many Gibbs rounds until replica swpas occurs
		context.setSwapEvery(std::stoi(setupReader.get("REX_SWAP_EVERY")[0]));

		context.setSwapFixman(std::stoi(setupReader.get("REX_SWAP_FIXMAN")[0]));

	}

	//std::cout << "OS memory 5.\n" << exec("free") << std::endl;
	// -- Run --
	if(setupReader.get("RUN_TYPE")[0] == "SimulatedTempering") {
		context.RunSimulatedTempering(context.getNofRounds(),
			std::stof(setupReader.get("TEMPERATURE_INI")[0]),
			std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
	}else if(setupReader.get("RUN_TYPE")[0] == "REX"){
		context.RunREX();
	}else{
		context.Run(context.getNofRounds(),
			std::stof(setupReader.get("TEMPERATURE_INI")[0]),
			std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
	}

	// Write final pdbs
	context.writeFinalPdb();

	//std::cout << "OS memory 6\n" << exec("free") << std::endl;
	//std::cout << "printStatus 1 " << std::endl;
	//context.printStatus();

	return 0;
} // END MAIN ////////////


///////////////////////////////////////////////////////////////////////////////
/////////////////////////////   MAIN   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//int main(int argc, char **argv)
//{
//
//	int retVal = run(argc, argv);
//	return retVal;
//
//}





