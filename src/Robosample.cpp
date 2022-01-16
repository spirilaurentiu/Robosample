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

int main(int argc, char **argv)
{

	// Read input into a SetupReader object
	SetupReader setupReader;
	if ( !LoadInputIntoSetupReader(argc, argv, setupReader) ){
		return 1;
	}
	//std::exit(0);

	// Declare global variables
	int currentWorldIx = 0;
	int round_mcsteps = 0;

	/////////// Create context ////////////

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
	
	// Adaptive Gibbs blocking
	context.setNofRoundsTillReblock(
		std::stoi((setupReader.get("ROUNDS_TILL_REBLOCK"))[0]));

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
	// TODO Move visualizer frequencies in Context in ValidateInput

	// Deal with visualizer. 
	SimTK::Real visualizerFrequency = -1;
	if ( setupReader.get("VISUAL")[0] == "FALSE" ){
		visualizerFrequency = -1;
	}else{
		visualizerFrequency = std::stod(setupReader.get("TIMESTEPS")[0]);
	}

	context.addEmptyWorlds(setupReader.get("WORLDS").size(),
		visualizerFrequency);

	int finalNofWorlds = context.getNofWorlds();
	if(requestedNofWorlds != finalNofWorlds){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}
	std::cout << "Added " << finalNofWorlds << " worlds" << std::endl;

	/////////////////////////////////////
	// 	STAGE EMPTY
	/////////////////////////////////////

	//////// FORCE FIELD //////////
	// Request threads 
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
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

	// Add filenames to Context filenames vectors
	// This has to be called before Worlds constructors so that
	// reserve will be called for molecules and topologies
	//int nofMols = static_cast<int>(setupReader.get("MOLECULES").size());

//	for(int molIx = 0; molIx < requestedNofMols; molIx++){
//		std::string topFN = 
//			setupReader.get("MOLECULES")[molIx] + std::string("/")
//			+ setupReader.get("PRMTOP")[molIx];
//
//		std::string crdFN = 
//			setupReader.get("MOLECULES")[molIx] + std::string("/")
//			+ setupReader.get("INPCRD")[molIx];
//
//		context.loadTopologyFile( topFN );
//
//		context.loadCoordinatesFile( crdFN );
//	}
	// Set top and crd FNs in Topology

	// Flexibility files 
//	for(int molIx = 0; molIx < requestedNofMols; molIx++){
//		for(int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
//			context.loadRigidBodiesSpecs( worldIx, molIx,
//				setupReader.get("MOLECULES")[molIx] + std::string("/")
//				+ setupReader.get("RBFILE")[(requestedNofMols * worldIx) + molIx]
//			);
//
//			context.loadFlexibleBondsSpecs( worldIx, molIx,
//				setupReader.get("MOLECULES")[molIx] + std::string("/")
//				+ setupReader.get("FLEXFILE")[(requestedNofMols * worldIx) + molIx]
//			);
//
//			context.setRegimen( worldIx, molIx,
//				setupReader.get("WORLDS")[worldIx] ); // TODO: delete from Topology
//		}
//	}

	// Add molecules to Worlds based on just read filenames
	context.AddMolecules(
		requestedNofMols,
		setupReader,
		setupReader.get("ROOTS"),
		setupReader.get("ROOT_MOBILITY"));

	int finalNofMols = context.getNofMolecules();
	if(requestedNofMols != finalNofMols){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}

	std::cout << "Added " << finalNofMols << " molecules" << std::endl;

	// Add membrane
	/*
	float memXWidth = std::stof(setupReader.get("MEMBRANE")[0]);
	float memYWidth = std::stof(setupReader.get("MEMBRANE")[1]);
	float memZWidth = std::stof(setupReader.get("MEMBRANE")[2]);
	int memResolution = std::stof(setupReader.get("MEMBRANE")[3]);
	bool haveMembrane = (memXWidth > SimTK::TinyReal) && (memYWidth > SimTK::TinyReal) && (memZWidth > SimTK::TinyReal);
	if(haveMembrane){
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
			(context.updWorld(worldIx))->addMembrane(memXWidth, memYWidth, memZWidth, memResolution);
		}
	}
	*/

	// Link the Compounds to Simbody System for all Worlds
	// This calls CompoundSystem::modelOneCompound function
	context.modelTopologies(setupReader.get("ROOT_MOBILITY"));


	// Add Fixman torque (Additional ForceSubsystem) if required
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			context.addFixmanTorque(worldIx);
		}
	}

	// Realize topology for all the Worlds
	context.realizeTopology();
	context.allocateReblockQsCacheQVectors();

	// Add samplers to the worlds
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		if(setupReader.get("SAMPLER")[worldIx] == "VV"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::HMC);
			pHMC(p)->setAlwaysAccept(true);
			pHMC(p)->setThermostat(ThermostatName::ANDERSEN);
		}else if(setupReader.get("SAMPLER")[worldIx] == "HMC"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::HMC);
			pHMC(p)->setAlwaysAccept(false);
			//pHMC(p)->setThermostat(ThermostatName::NONE);
			pHMC(p)->setThermostat(ThermostatName::ANDERSEN);
		}/*else if(setupReader.get("SAMPLER")[worldIx] == "LAHMC"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::LAHMC);
			pLAHMC(p)->setThermostat(ThermostatName::NONE);
		}else{
			BaseSampler *p = context.addSampler(worldIx, SamplerName::LAHMC);
			pLAHMC(p)->setThermostat(ThermostatName::NONE);
		}*/
	}

	// Set sampler parameters and initialize
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {

			// Set timesteps
			context.setTimestep(worldIx, samplerIx, std::stof(setupReader.get("TIMESTEPS")[worldIx]));

			// Activate Fixman potential if needed
			if(setupReader.get("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
				 context.useFixmanPotential(worldIx, samplerIx);
			}
		}
	}

	// This loop is just for check purposes (should be removed)
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			std::cout << "After setThermo world " << worldIx << " sampler " << samplerIx << "getThermostat: " ;
			std::cout << static_cast<int>(pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->getThermostat()) ;
			std::cout << std::endl;
		}
	}


	// Thermodynamics
	context.setTemperature(std::stof(setupReader.get("TEMPERATURE_INI")[0]));
	for(unsigned int worldIx = 0; worldIx < setupReader.get("WORLDS").size(); worldIx++){
		context.setNofSamplesPerRound(worldIx, std::stoi(setupReader.get("SAMPLES_PER_ROUND")[worldIx]));
		context.setNofMDStepsPerSample(worldIx, 0, std::stoi(setupReader.get("MDSTEPS")[worldIx]));
	}

	if(setupReader.find("MDSTEPS_STD")){
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
			for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {
				pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setMDStepsPerSampleStd(
						std::stoi(setupReader.get("MDSTEPS_STD")[worldIx]));
			}
		}
	}else{
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
			for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {
				pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setMDStepsPerSampleStd(0);
			}
		}
	}

	// NEW PARAMS
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {
			pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostTemperature(
				std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]));
			pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostMDSteps(
					std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]));
		}
	}

	// Set the seeds for reproducibility. Samplers have to be here already
	if( setupReader.find("SEED") ){
		if( !(setupReader.get("SEED").empty()) ){
			for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
				context.setSeed(worldIx, 0, std::stoi(setupReader.get("SEED")[worldIx] ));
			}
		}
	}

	// Initialize samplers
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			context.initializeSampler(worldIx, samplerIx);
		}
	}

	// Print thermodynamics
	for(unsigned int worldIx = 0; worldIx < setupReader.get("WORLDS").size(); worldIx++){
		std::cout << "MAIN World " << worldIx << " temperature = "
			<< context.getWorld(worldIx)->getTemperature()
			<< std::endl;
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			std::cout << "MAIN World " << worldIx
			<< " FixmanTorque temperature = "
			<< context.updWorld(worldIx)->updFixmanTorque()->getTemperature()
			<< std::endl;
		}
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			std::cout << "MAIN World " << worldIx << " Sampler " << samplerIx 
				<< " temperature = " << context.updWorld(worldIx)->updSampler(samplerIx)->getTemperature()
				<< " initial const state PE: " << std::setprecision(20)
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context.updAdvancedState(worldIx, samplerIx))
				<< " useFixmanPotential = "
				<< pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->isUsingFixmanPotential()
				<< std::endl;
		}

	}

	// Simulation parameters
	currentWorldIx = 0;
	round_mcsteps = 0;

	context.setNofRounds(std::stoi(setupReader.get("ROUNDS")[0]));

	if (setupReader.get("RUN_TYPE")[0] == "SimulatedTempering") {
		for (unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
			context.setNofBoostStairs(worldIx, std::stoi(setupReader.get("BOOST_STAIRS")[worldIx]));
		}
	}

	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){    
		round_mcsteps += context.getNofSamplesPerRound(worldIx);
	}

	// Set pdb writing frequency
	context.setPdbRestartFreq( std::stoi(setupReader.get("WRITEPDBS")[0]) );

	// Get atom indeces for geometry calculations
	if(setupReader.get("GEOMETRY")[0] == "TRUE"){
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
			std::vector<int> distanceIx;
			for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
				distanceIx.push_back(atoi(setupReader.get("DISTANCE")[i].c_str()));
			}
			// Get distances
			for(unsigned int ai = 0; ai < setupReader.get("DISTANCE").size() / 2; ai++){
				context.addDistance(worldIx, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
			}
	
			// Get dihedrals indeces
			std::vector<int> dihedralIx;
			for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
				dihedralIx.push_back(atoi(setupReader.get("DIHEDRAL")[i].c_str()));
			}
			// Get dihedrals
			for(unsigned int ai = 0; ai < setupReader.get("DIHEDRAL").size() / 4; ai++){
				context.addDihedral(worldIx, 0,
					dihedralIx[4*ai + 0], dihedralIx[4*ai + 1],
					dihedralIx[4*ai + 2], dihedralIx[4*ai + 3]);
			}
		}
	}

	// Write initial pdb for debug purposes: TODO remove 
	// - we need this to get compound atoms
	currentWorldIx = context.worldIndexes.front();
	SimTK::State& advancedState = (context.updWorld(currentWorldIx))->integ->updAdvancedState();

	// - Get molecules directory
	std::string path = setupReader.get("MOLECULES")[0];
	std::size_t lastSlashPos = path.find_last_of("/");
	std::string upDir = path.substr(0, lastSlashPos);
	std::string molDir = path.substr(lastSlashPos + 1, path.length());
	std::cout << "Molecule directory: " << molDir << std::endl;

	// - Get seed
	context.setOutputDir(setupReader.get("OUTPUT_DIR")[0] );
	context.setPdbPrefix(molDir
		+ std::to_string(context.updWorld(currentWorldIx)->updSampler(0)->getSeed()) 
		);

	// - Write
	constexpr int mc_step = -1;
	if(setupReader.get("WRITEPDBS")[0] == "TRUE"){
		(context.updWorld(currentWorldIx))->updateAtomListsFromCompound(advancedState);
		std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
		for(unsigned int mol_i = 0; mol_i < setupReader.get("MOLECULES").size(); mol_i++){
			((context.updWorld(currentWorldIx))->getTopology(mol_i)).writeAtomListPdb(context.getOutputDir(),
			"/pdbs/sb." + context.getPdbPrefix() + ".", ".pdb", 10, mc_step);
		}
	}

	// Alert user of CUDA environment variables
	if(SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT").empty()){
		std::cout << "CUDA_ROOT not set." << std::endl;
	}else{
		std::cout << "CUDA_ROOT set to " 
			<< SimTK::Pathname::getEnvironmentVariable("CUDA_ROOT") << std::endl;
	}

	// Get output printing frequency
	context.setPrintFreq( std::stoi(setupReader.get("PRINT_FREQ")[0]) );

	// Realize topology for all the Worlds
	context.realizeTopology();

	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		(context.updWorld(worldIx))->loadCompoundRelatedMaps();
		(context.updWorld(worldIx))->loadMobodsRelatedMaps();
	}

	// Membrane contacts. TODO: only one per world for now
	/*
	if(haveMembrane){
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
			(context.updWorld(worldIx))->addContacts( std::stoi(setupReader.get("CONTACTS")[worldIx]) );
		}
	}
	*/

	// Set ocnstraint. TODO: only one per world allowed for now
	/*
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		(context.updWorld(worldIx))->addConstraints( std::stoi(setupReader.get("CONSTRAINTS")[worldIx]) );
	}
	*/

	// U Scale Factors uses maps stored in Topology
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		(context.updWorld(worldIx))->setUScaleFactorsToMobods();
	}

	// Realize topology for all the Worlds
	context.realizeTopology();

	// Load/store Mobilized bodies joint types
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			SimTK::State& curAvancedState = (context.updWorld(worldIx))->integ->updAdvancedState();
			std::cout << "Loading mbx2mobility" << std::endl;
			(context.updWorld(worldIx)->updSampler(samplerIx))->loadMbx2mobility(curAvancedState);
		}
	}


	// -- Run --
	if(setupReader.get("RUN_TYPE")[0] == "SimulatedTempering") {
		context.RunSimulatedTempering(context.getNofRounds(),
					 std::stof(setupReader.get("TEMPERATURE_INI")[0]),
					 std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
	}else{
		context.Run(context.getNofRounds(),
					 std::stof(setupReader.get("TEMPERATURE_INI")[0]),
					 std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
	}

	// Write final pdbs
	for(unsigned int mol_i = 0; mol_i < setupReader.get("MOLECULES").size(); mol_i++){
		((context.updWorld(context.worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(
				context.getOutputDir(), "/pdbs/final." + context.getPdbPrefix() + ".", ".pdb", 10,
				context.getNofRounds());
	}

	//std::cout << "printStatus 1 " << std::endl;
	//context.printStatus();

	return 0;
} // END MAIN ////////////







