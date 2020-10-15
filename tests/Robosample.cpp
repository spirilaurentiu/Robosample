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

    std::cout << "Setting Robosample up..." << std::endl;

    // Initialize setup reader
    SetupReader setupReader(argv[1]);
    setupReader.dump(true);

    const unsigned int nofWorlds = setupReader.get("WORLDS").size();
    std::cout << "Robosample nofWorlds " << nofWorlds << std::endl;

    // Create pdbs directory if necessary
    if( !SimTK::Pathname::fileExists(setupReader.get("OUTPUT_DIR")[0] + "/pdbs") ){
        const int err = mkdir((setupReader.get("OUTPUT_DIR")[0] + "/pdbs").c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        if (err == -1){
            std::cout << "Error creating " << setupReader.get("OUTPUT_DIR")[0] + "/pdbs" << std::endl;
            return 1;
        }
    }

    // Variables
    int currentWorldIx = 0;
    int round_mcsteps = 0;

    // Create context
    std::string logFilename;
    if (setupReader.find("SEED")) {
        if (!(setupReader.get("SEED").empty())) {
            logFilename = setupReader.get("OUTPUT_DIR")[0] + std::string("/log.") + setupReader.get("SEED")[0];
        }
    } else {
        logFilename = "x";
    }

    // Initialize a Context
    Context context(setupReader, logFilename);

    // Add empty Worlds to the Context
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.get("VISUAL")[worldIx] == "TRUE"){
            if(std::stod(setupReader.get("TIMESTEPS")[worldIx]) <= 0.0000001){
            	context.AddWorld(true, 0.001);
            }else{
            	context.AddWorld(true, std::stod(setupReader.get("TIMESTEPS")[worldIx]));
            }
        }else{
            context.AddWorld(false);
        }
    }

//@    // Request threads
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
        context.setNumThreadsRequested(worldIx,
                std::stod(setupReader.get("THREADS")[worldIx]));
    }

    // Print the number of threads we got back
    context.PrintNumThreads();

    // Set worlds force field scale factors
    for (unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
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
//@
    // Add filenames to Context filenames vectors
    unsigned int nofMols = setupReader.get("MOLECULES").size();
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for(unsigned int molIx = 0; molIx < nofMols; molIx++){

            context.loadTopologyFile( worldIx, molIx,
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                    + setupReader.get("PRMTOP")[molIx]
            );

            context.loadCoordinatesFile( worldIx, molIx,
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                    + setupReader.get("INPCRD")[molIx]
            );

            context.loadRigidBodiesSpecs( worldIx, molIx,
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                    + setupReader.get("RBFILE")[(nofMols * worldIx) + molIx]
            );

            context.loadFlexibleBondsSpecs( worldIx, molIx,
                    setupReader.get("MOLECULES")[molIx] + std::string("/")
                    + setupReader.get("FLEXFILE")[(nofMols * worldIx) + molIx]
            );

            context.setRegimen( worldIx, molIx,
                    setupReader.get("WORLDS")[worldIx] );


        }
    }

    //for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
            // @ DELETE FOR SINGLE MOLECULE
            // Add molecules to Worlds based on just read filenames
            //context.AddMolecules(setupReader.get("ROOTS")); // @
    //}
    // @ RESTORE FOR SINGLE MOLECULE
    // Add molecules to Worlds based on just read filenames
    context.AddMolecules(setupReader.get("ROOTS"));

    // BEGIN MULMOL
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        (context.updWorld(worldIx))->realizeTopology();
    }

/*    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        (context.updWorld(worldIx))->loadMobodsRelatedMaps();
    }*/
    // END MULMOL
    

    // Use OpenMM if possible
    if(setupReader.get("OPENMM")[0] == "TRUE"){
        context.setUseOpenMMAcceleration(true);
    }

    // Set Lennard-Jones mixing rule
    context.setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

    // Link the Compounds to Simbody System for all Worlds
    context.modelTopologies(setupReader.get("ROOT_MOBILITY"));

    // Add Fixman torque (Additional ForceSubsystem) if required
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        if(setupReader.get("FIXMAN_TORQUE")[worldIx] == "TRUE"){
            context.addFixmanTorque(worldIx);
        }
    }



/*
    // Set worlds force field scale factors
    for (unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
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
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
        context.setNumThreadsRequested(worldIx,
                std::stod(setupReader.get("THREADS")[worldIx]));
    }

    // Print the number of threads we got back
    context.PrintNumThreads();
*/

    // Realize topology for all the Worlds
    context.realizeTopology();

   // Add samplers to the worlds
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
	if(setupReader.get("SAMPLER")[worldIx] == "HMC"){
        	context.addSampler(worldIx, HMC);
	}else if(setupReader.get("SAMPLER")[worldIx] == "LAHMC"){
        	context.addSampler(worldIx, LAHMC);
	}else{
        	context.addSampler(worldIx, LAHMC);
	}
    }

////////// DEBUG BEGIN
	for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
		for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {}}
////////// DEBUG END

    // Set sampler parameters and initialize
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++) {
        for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++) {
            // Set timesteps
            context.setTimestep(worldIx, samplerIx, std::stod(setupReader.get("TIMESTEPS")[worldIx]));

            // Set thermostats
            pMC(context.updWorld(worldIx)->updSampler(samplerIx))->setThermostat(setupReader.get("THERMOSTAT")[worldIx]);
            pLAHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostTemperature(
                std::stof(setupReader.get("BOOST_TEMPERATURE")[worldIx]));
            pLAHMC(context.updWorld(worldIx)->updSampler(samplerIx))->setBoostMDSteps(
                    std::stoi(setupReader.get("BOOST_MDSTEPS")[worldIx]));

            // Activate Fixman potential if needed
            if(setupReader.get("FIXMAN_POTENTIAL")[worldIx] == "TRUE"){
                 context.useFixmanPotential(worldIx, samplerIx);
            }
        }
    }

    // This loop is just for check purposes (should be removed)
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
            std::cout << "After setThermo world " << worldIx << " sampler " << samplerIx << "getThermostat: " ;
            std::cout << pMC(context.updWorld(worldIx)->updSampler(samplerIx))->getThermostat() ;
            std::cout << std::endl;
        }
    }


    // Set thermodynamics
    for(unsigned int worldIx = 0; worldIx < setupReader.get("WORLDS").size(); worldIx++){
        context.setTemperature(worldIx, std::stof(setupReader.get("TEMPERATURE_INI")[worldIx]));
        context.setNofSamplesPerRound(worldIx, std::stoi(setupReader.get("SAMPLES_PER_ROUND")[worldIx]));
        context.setNofMDStepsPerSample(worldIx, 0, std::stoi(setupReader.get("MDSTEPS")[worldIx]));
    }

    // Set the seeds for reproducibility
    if( setupReader.find("SEED") ){
        if( !(setupReader.get("SEED").empty()) ){
            for(unsigned int worldIx = 0; worldIx < context.getNofWorlds(); worldIx++){
                context.setSeed(worldIx, 0, std::stoi(setupReader.get("SEED")[worldIx] ));
            }
        }
    }

    // Initialize samplers
    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        for (int samplerIx = 0; samplerIx < context.getWorld(worldIx)->getNofSamplers(); samplerIx++){
            context.initializeSampler(worldIx, samplerIx);
        }
    }

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
                << " useFixmanPotential = " << pMC(context.updWorld(worldIx)->updSampler(samplerIx))->isUsingFixmanPotential()
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
            int distanceIx[setupReader.get("DISTANCE").size()];
            for(unsigned int i = 0; i < setupReader.get("DISTANCE").size(); i++){
                distanceIx[i] = atoi(setupReader.get("DISTANCE")[i].c_str());
            }
            // Get distances
            for(unsigned int ai = 0; ai < setupReader.get("DISTANCE").size() / 2; ai++){
                context.addDistance(worldIx, 0, distanceIx[2*ai + 0], distanceIx[2*ai + 1]);
            }
    
            // Get dihedrals indeces
            int dihedralIx[setupReader.get("DIHEDRAL").size()];
            for(unsigned int i = 0; i < setupReader.get("DIHEDRAL").size(); i++){
                dihedralIx[i] = atoi(setupReader.get("DIHEDRAL")[i].c_str());
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

    // Write pdb
    context.setOutputDir(setupReader.get("OUTPUT_DIR")[0] );
    context.setPdbPrefix(setupReader.get("MOLECULES")[0]
        + std::to_string(context.updWorld(currentWorldIx)->updSampler(0)->getSeed()) 
        );


    // Helper constant for step arithmetic
    constexpr int mc_step = -1;
    if(setupReader.get("WRITEPDBS")[0] == "TRUE"){
        (context.updWorld(currentWorldIx))->updateAtomListsFromCompound(advancedState);
        std::cout << "Writing pdb  sb" << mc_step << ".pdb" << std::endl;
        for(unsigned int mol_i = 0; mol_i < setupReader.get("MOLECULES").size(); mol_i++){
            ((context.updWorld(currentWorldIx))->getTopology(mol_i)).writeAtomListPdb(context.getOutputDir(),
                                                                                       "/pdbs/sb." +
                                                                                       context.getPdbPrefix() + ".",
                                                                                       ".pdb", 10, mc_step);
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

    for(unsigned int worldIx = 0; worldIx < nofWorlds; worldIx++){
        (context.updWorld(worldIx))->loadMobodsRelatedMaps();
    }


    // -- Run --
    bool isWorldsOrderRandom = ((setupReader.get("RANDOM_WORLD_ORDER")[0] == "TRUE") ? true : false);
    if(setupReader.get("RUN_TYPE")[0] == "SimulatedTempering") {
        context.RunSimulatedTempering(context.getNofRounds(),
                     std::stof(setupReader.get("TEMPERATURE_INI")[0]),
                     std::stof(setupReader.get("TEMPERATURE_FIN")[0]));
    }else{
        context.Run(context.getNofRounds(),
                     std::stof(setupReader.get("TEMPERATURE_INI")[0]),
                     std::stof(setupReader.get("TEMPERATURE_FIN")[0]), isWorldsOrderRandom);
    }

    // Write final pdbs
    for(unsigned int mol_i = 0; mol_i < setupReader.get("MOLECULES").size(); mol_i++){
        ((context.updWorld(context.worldIndexes.front()))->getTopology(mol_i)).writeAtomListPdb(
                context.getOutputDir(), "/pdbs/final." + context.getPdbPrefix() + ".", ".pdb", 10,
                context.getNofRounds());
    }

    return 0;
} // END MAIN ////////////







