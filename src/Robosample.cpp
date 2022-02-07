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

//#ifndef ROBO_DEBUG_LEVEL01
//#define ROBO_DEBUG_LEVEL01
//#endif

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
		const int err = mkdir((outDir 
			+ "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (err == -1){
			std::cout << "Error creating " << outDir + "/pdbs" << std::endl;
			return false;
		}
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




///////////////////////////////////////////////////////////////////////////////
/////////////////////////////   MAIN   ////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char **argv)
{

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
	std::cout << "Molecule directory: " << molDir << std::endl;

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
	std::cout << "Requested " << requestedNofMols << " molecules\n";

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

	// Set Lennard-Jones mixing rule
	context.setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Add molecules based on the setup reader
	context.AddMolecules(requestedNofMols, setupReader);

	int finalNofMols = context.getNofMolecules();
	if(requestedNofMols != finalNofMols){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}
	std::cout << "Added " << finalNofMols << " molecules" << std::endl;

	// Loads parameters into DuMM, adopts compound by the CompoundSystem
	// and loads maps of indexes
	context.initializeWorlds(finalNofMols, setupReader);

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
		BaseSampler *p = context.addSampler(
			worldIx, setupReader.get("SAMPLER")[worldIx]);
	}

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

	// Check atom stations for debug purposes
	context.checkAtomStationsThroughDumm();
	
	// Load/store Mobilized bodies joint types in samplers
	context.loadMbxsToMobilities();

	// -- Setup REX --
	if(setupReader.get("RUN_TYPE")[0] == "REX"){
		// Look for trex.dat
		SetupReader rexReader;

		// Read input
		std::cout << "Reading " << setupReader.get("REX_FILE")[0] << "\n";
		rexReader.ReadSetup(setupReader.get("REX_FILE")[0]);
		rexReader.dump(true);

		// Use setup reader for REX file. TODO: adapt reader to format 
		size_t nofReplicas = rexReader.getNofKeys();
		context.setNofReplicas(nofReplicas);
		context.setNofThermodynamicStates(nofReplicas);

		// Add replicas
		for(size_t replica_k = 0; replica_k < nofReplicas; replica_k++){

			std::vector<std::string> vals_k = rexReader.get(std::to_string(replica_k));

			std::cout << "rexReader: "
				<< vals_k[0]
				<< std::endl;

			// Add a thermodynamic state
			SimTK::Real T = std::stod(rexReader.get(std::to_string(replica_k))[0]);
				
			context.addThermodynamicState(replica_k, T);

			// Add worlds to replica
			std::vector<int> worldIndexes_k;

			// Add the first world to all the replicas
			worldIndexes_k.push_back(0);

			// Add the rest of the worlds
			for(size_t i = 1; i < vals_k.size(); i++){
				std::cout << "Add worlds " << std::stoi(vals_k[i])
					<< " to replica " << replica_k
					<< std::endl;
			
				worldIndexes_k.push_back(std::stoi(vals_k[i]));
			}

			std::cout << "About to add replica\n" << std::flush;
			context.addReplica(replica_k, worldIndexes_k);

		}

		context.loadReplica2ThermoIxs();

		context.PrintReplicas();
		
	}
	//std::exit(0);

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

	//std::cout << "printStatus 1 " << std::endl;
	//context.printStatus();

	return 0;
} // END MAIN ////////////







