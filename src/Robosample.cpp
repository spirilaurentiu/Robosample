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

void PrintHelp() {
	std::cout <<
		"Usage: Robsample [options]\n" <<
		"Options:\n" <<
		"  -h, --help for help\n" <<
		"Usage: Robsample file\n";
}

int main(int argc, char **argv)
{
	// ios_base::sync_with_stdio(false);
    // std::cin.tie(nullptr);

	/////////////////////////////////////
	// 	NO SIMBODY OBJECTS YET
	/////////////////////////////////////

	// Check if there is any input
	if(argc < 2) {
		std::cout << "Error: not enough parameters to run. See help below.\n";
		PrintHelp();

		return 1;
	}
	else {
		auto arg = std::string(argv[1]);
		if("-h" == arg || "--help" == arg) {
			PrintHelp();
			return 1;
		}
	}

	// Initialize setup reader
	std::cout << "Got the following input:" << std::endl;
	SetupReader setupReader(argv[1]);
	setupReader.dump(true);

	// Declare global variables
	int currentWorldIx = 0;
	int round_mcsteps = 0;

	/////////// Create context ////////////

	// Create pdbs directory if necessary
	if( !SimTK::Pathname::fileExists(setupReader.get("OUTPUT_DIR")[0] + "/pdbs") ){
		const int err = mkdir((setupReader.get("OUTPUT_DIR")[0] 
			+ "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
		if (err == -1){
			std::cout << "Error creating " << setupReader.get("OUTPUT_DIR")[0] + "/pdbs" << std::endl;
			return 1;
		}
	}

	// Set the log filename
	std::string logFilename;
	if (setupReader.find("SEED")) {
		if (!(setupReader.get("SEED").empty())) {
			logFilename = setupReader.get("OUTPUT_DIR")[0] 
				+ std::string("/log.") + setupReader.get("SEED")[0];
		}
	} else {
		logFilename = "x";
	}

	Context context(setupReader, logFilename);

	int requestedNofWorlds = context.getNofWorlds();
	int requestedNofMols = context.getNofMolecules();
	//std::cout << "Requested " << requestedNofMols << " molecules\n";
	//std::exit(0);

	/////////// Add Worlds to context ////////////
	// Add Worlds to the Context. Every World instantiates a: 
	// CompoundSystem, SimbodyMatterSubsystem, GeneralForceSubsystem,
	// DuMMForceSubsystem, Integrator, TimeStepper and optionally:
	// DecorationSubsystem, Visualizer, VisuzlizerReporter,
	//  ParaMolecularDecorator
	// TODO Move visualizer frequencies in Context in ValidateInput
	for(unsigned int worldIx = 0;
		worldIx < setupReader.get("WORLDS").size(); 
		worldIx++){
		if(setupReader.get("VISUAL")[worldIx] == "TRUE"){
				context.AddWorld(true, std::stod(setupReader.get("TIMESTEPS")[worldIx]));
		}else{
			context.AddWorld(false);
		}
	}

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
//@    // Request threads 
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


//@
	// Add filenames to Context filenames vectors
	// This has to be called before Worlds constructors so that
	// reserve will be called for molecules and topologies
	//int nofMols = static_cast<int>(setupReader.get("MOLECULES").size());

	for(int molIx = 0; molIx < requestedNofMols; molIx++){
		context.loadTopologyFile(
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("PRMTOP")[molIx]
		);

		context.loadCoordinatesFile(
			setupReader.get("MOLECULES")[molIx] + std::string("/")
			+ setupReader.get("INPCRD")[molIx]
		);
	}

	for(int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){

		for(int molIx = 0; molIx < requestedNofMols; molIx++){
			context.loadRigidBodiesSpecs( worldIx, molIx,
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("RBFILE")[(requestedNofMols * worldIx) + molIx]
			);

			context.loadFlexibleBondsSpecs( worldIx, molIx,
				setupReader.get("MOLECULES")[molIx] + std::string("/")
				+ setupReader.get("FLEXFILE")[(requestedNofMols * worldIx) + molIx]
			);

			context.setRegimen( worldIx, molIx,
				setupReader.get("WORLDS")[worldIx] ); // TODO: delete from Topology
		}
	}

	// Add molecules to Worlds based on just read filenames
	context.AddMolecules(setupReader.get("ROOTS"),
		setupReader.get("ROOT_MOBILITY"));

	int finalNofMols = context.getNofMolecules();
	if(requestedNofMols != finalNofMols){
		std::cerr << "Something went wrong while adding the world\n";
		throw std::exception();
		std::exit(1);
	}

	std::cout << "Added " << finalNofMols << " molecules" << std::endl;

	// BEGIN MULMOL
	//for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
	//	(context.updWorld(worldIx))->realizeTopology();
	//}

/*
	// Membrane
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

	// Use OpenMM if possible
	//if(setupReader.get("OPENMM")[0] == "TRUE"){
	//	context.setUseOpenMMAcceleration(true);
	//}

	if(setupReader.get("OPENMM_CalcOnlyNonbonded")[0] == "TRUE"){
		context.setUseOpenMMCalcOnlyNonBonded(true);
	}
	else context.setUseOpenMMCalcOnlyNonBonded(false);

    for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
        // Only NoCutoff (0) and CutoffNonPeriodic(1) methods are supported. Additional 3 methods available in
        // OpenMM could be integrated as well.
        if(setupReader.get("NONBONDED_METHOD")[worldIx] == "1"){
            context.setNonbondedMethod(worldIx, 1);
            context.setNonbondedCutoff(worldIx, std::stod( setupReader.get("NONBONDED_CUTOFF")[worldIx] ) );
        }
    }



	// Set Lennard-Jones mixing rule
	//context.setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Link the Compounds to Simbody System for all Worlds
	context.modelTopologies(setupReader.get("ROOT_MOBILITY"));
	//context.printStatus();
	//std::exit(0);

	/////////////////////////////////////
	// 	STAGE EMPTY
	/////////////////////////////////////

	//std::cout << "printStatus 4 " << std::endl;
	//context.printStatus();
	//context.PrintMolmodelAndDuMMTypes();
	//context.PrintSimbodyMobods();
	//std::exit(0);

	// Add Fixman torque (Additional ForceSubsystem) if required
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			context.addFixmanTorque(worldIx);
		}
	}

/*
	// Set worlds force field scale factors
	for (unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		// Set force field scale factors.
		if(setupReader.get("FFSCALE")[worldIx] == "AMBER"){
			context.setAmberForceFieldScaleFactors(worldIx);
		}else{
			context.setGlobalForceFieldScaleFactor(worldIx,
					std::stod(setupReader.get("FFSCALE")[worldIx]));
		}
		// Set world GBSA scale factor
		context.setGbsaGlobalScaleFactor(worldIx,
				std::stod(setupReader.get("GBSA")[worldIx]));
	}
*/
/*
	// Request threads
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		context.setNumThreadsRequested(worldIx,
				std::stod(setupReader.get("THREADS")[worldIx]));
	}

	// Print the number of threads we got back
	context.PrintNumThreads();
*/

	// Realize topology for all the Worlds
	context.realizeTopology();

   // Add samplers to the worlds
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		if(setupReader.get("SAMPLER")[worldIx] == "VV"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::HMC);
			pHMC(p)->setThermostat(ThermostatName::ANDERSEN);
		}else if(setupReader.get("SAMPLER")[worldIx] == "HMC"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::HMC);
			pHMC(p)->setThermostat(ThermostatName::NONE);
		}/*else if(setupReader.get("SAMPLER")[worldIx] == "LAHMC"){
			BaseSampler *p = context.addSampler(worldIx, SamplerName::LAHMC);
			pLAHMC(p)->setThermostat(ThermostatName::NONE);
		}else{
			BaseSampler *p = context.addSampler(worldIx, SamplerName::LAHMC);
			pLAHMC(p)->setThermostat(ThermostatName::NONE);
		}*/
	}

	//std::cout << "printStatus 6 " << std::endl;
	//context.printStatus();
////////// DEBUG BEGIN
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {}}
////////// DEBUG END

	// Set sampler parameters and initialize
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++) {
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {
			// Set timesteps
			context.setTimestep(worldIx, samplerIx, std::stof(setupReader.get("TIMESTEPS")[worldIx]));

			// Set thermostats
			//pMC(context.updWorld(worldIx)->updSampler(samplerIx))->setThermostat(setupReader.get("THERMOSTAT")[worldIx]);
			pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostTemperature(
				std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]));
			pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostMDSteps(
					std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]));

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

	// Set the seeds for reproducibility. Samplers have to be here already
	if( setupReader.find("SEED") ){
		if( !(setupReader.get("SEED").empty()) ){
			for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
				context.setSeed(worldIx, 0, std::stoi(setupReader.get("SEED")[worldIx] ));
			}
		}
	}

	//std::cout << "printStatus 6.5 " << std::endl;
	//context.printStatus();
	// Initialize samplers
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			context.initializeSampler(worldIx, samplerIx);
		}
	}

	//std::cout << "printStatus 7 " << std::endl;
	//context.printStatus();
	// Print thermodynamics
	for(unsigned int worldIx = 0; worldIx < setupReader.get("WORLDS").size(); worldIx++){
		std::cout << "MAIN World " << worldIx << " temperature = " << context.getWorld(worldIx)->getTemperature() << std::endl;
		if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
			std::cout << "MAIN World " << worldIx << " FixmanTorque temperature = " << context.updWorld(worldIx)->updFixmanTorque()->getTemperature() << std::endl;
		}
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
			std::cout << "MAIN World " << worldIx << " Sampler " << samplerIx 
				<< " temperature = " << context.updWorld(worldIx)->updSampler(samplerIx)->getTemperature()
				<< " initial const state PE: " << std::setprecision(20)
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy((updWorld(worldIx))->integ->updAdvancedState())
				//<< (context.updWorld(worldIx))->forces->getMultibodySystem().calcPotentialEnergy(context.updAdvancedState(worldIx, samplerIx))
				<< " useFixmanPotential = " << pHMC(context.updWorld(worldIx)->updSampler(samplerIx))->isUsingFixmanPotential()
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

	// int total_mcsteps = round_mcsteps * context.getNofRounds();

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

	// // wth is this?
	// const SimTK::Compound * p_compounds[context.getNofWorlds()];
	// if(setupReader.get("GEOMETRY")[0] == "TRUE"){
	//     for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
	//         p_compounds[worldIx] = &((context.updWorld(worldIx))->getTopology(0));
	//     }
	// }
	// //

	// Update one round for the first regimen
	currentWorldIx = context.worldIndexes.front();
	SimTK::State& advancedState = (context.updWorld(currentWorldIx))->integ->updAdvancedState();

	// Find directory structure
	std::string path = setupReader.get("MOLECULES")[0];
	std::size_t lastSlashPos = path.find_last_of("/");
	std::string upDir = path.substr(0, lastSlashPos);
	std::string molDir = path.substr(lastSlashPos + 1, path.length());
	std::cout << "Molecule directory: " << molDir << std::endl;

	// Write pdb
	context.setOutputDir(setupReader.get("OUTPUT_DIR")[0] );
	//context.setPdbPrefix(setupReader.get("MOLECULES")[0]
	//	+ std::to_string(context.updWorld(currentWorldIx)->updSampler(0)->getSeed()) 
	//	);
	context.setPdbPrefix(molDir
		+ std::to_string(context.updWorld(currentWorldIx)->updSampler(0)->getSeed()) 
		);


	// Helper constant for step arithmetic
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

	//std::cout << "printStatus 8 " << std::endl;
	//context.printStatus();
	for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
		(context.updWorld(worldIx))->loadCompoundRelatedMaps();
		(context.updWorld(worldIx))->loadMobodsRelatedMaps();
	}

/*
	// Membrane contacts. TODO: only one per world for now
	if(haveMembrane){
		for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
			(context.updWorld(worldIx))->addContacts( std::stoi(setupReader.get("CONTACTS")[worldIx]) );
		}
	}

	// Set ocnstraint. TODO: only one per world allowed for now
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


	//std::cout << "printStatus 9 " << std::endl;
	//context.printStatus();

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







