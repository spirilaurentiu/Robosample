#include "Robo.hpp"
#include "HMCSampler.hpp"
#include "World.hpp"

using namespace SimTK;

int main(int argc, char **argv)
{
	if(argc < 4){
		std::cout << "Error: not enough parameters to run.\n";
		return 1;
	}

	Real timestep = 0.002;
	bool wantVisual = std::stoi(argv[4]);
	Real visualFreq = timestep;
	int nofRounds = std::stoi(argv[1]);
	int noMDStepsPerSample = 20;
	World *world = new World(0, wantVisual, visualFreq);
	bool useFixmanPotential = false;
	bool useFixmanTorque = false;

	// SimTK::CompoundSystem *compoundSystem = world->compoundSystem;
	SimTK::SimbodyMatterSubsystem *matter = world->matter.get();
	SimTK::GeneralForceSubsystem *generalForceSubsystem = world->forces.get();
	SimTK::DuMMForceFieldSubsystem *dumm = world->forceField.get();
	SimTK::TimeStepper *timeStepper = world->ts.get();
	const SimTK::System *system = &(matter->getSystem());

	readAmberInput amberReader;
	amberReader.readAmberFiles("etanMob/ligand.inpcrd", "etanMob/ligand.prmtop");
	world->AddMolecule(&amberReader, "etanMob/ligand.rb", "etanMob/ligand.flex", "TD", "0");

	// Using OpenMM acceleration
	//world->updForceField()->setUseOpenMMAcceleration(true);

	// Set Lennard-Jones mixing rule
	world->updForceField()->setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Link the Compounds to Simbody System for all Worlds
    world->modelTopologies("Cartesian");

	// Add Fixman torque (Additional ForceSubsystem) if required
	if(std::string(argv[3]) == "FT"){
		useFixmanTorque = true;
		world->addFixmanTorque();
		std::cout << "Fixman torque added." << std::endl;
	}

	// Set worlds force field scale factors
	world->setGlobalForceFieldScaleFactor(0.0);

        // Set world GBSA scale factor
	world->setGbsaGlobalScaleFactor(0.0);


	// Request threads
    world->updForceField()->setUseMultithreadedComputation(true);
    world->updForceField()->setNumThreadsRequested(2);

	// Print the number of threads we got back
	std::cout << "World  requested "
        	<< world->updForceField()->getNumThreadsRequested()
        	<< " and got "
        	<< world->updForceField()->getNumThreadsInUse()
        	<< std::endl;

	// Realize topology for all the Worlds
	world->getCompoundSystem()->realizeTopology();

	// AddAdd samplers to the worlds
	world->addSampler(SamplerName::HMC);
	BaseSampler *pSampler = world->updSampler(0);

	// Set timestep
	(pHMC(pSampler))->setTimestep(timestep);

	// Set thermostats
	pMC(pSampler)->setThermostat(ThermostatName::ANDERSEN);
	pHMC(pSampler)->setBoostTemperature(1);
	pHMC(pSampler)->setBoostMDSteps(1);

	// Activate Fixman potential if needed
	if(std::string(argv[2]) == "FP"){
		useFixmanPotential = true;
		pMC(pSampler)->useFixmanPotential();
		std::cout << "Fixman potential added." << std::endl;
	}

	// Set thermodynamics
        world->setTemperature(3000);
	pHMC(pSampler)->setMDStepsPerSample(noMDStepsPerSample);

	// Set the seeds for reproducibility
	pSampler->setSeed(0);

	// Initialize samplers
	SimTK::State& advancedState = world->integ->updAdvancedState();
	pSampler->initialize( advancedState );

        // Reinitialize current sampler
        pSampler->reinitialize(advancedState);

	printf("NU: %d \n", advancedState.getNU());

	const Topology& topology = world->getTopology(0);

	// Geometry
	SimTK::Compound::AtomIndex AIx, BIx, CIx, DIx;
	SimTK::Vec3 Apos, Bpos, Cpos, Dpos;
	// Lin4
/*	AIx = SimTK::Compound::AtomIndex(1);
	BIx = SimTK::Compound::AtomIndex(0);
	CIx = SimTK::Compound::AtomIndex(2);
	DIx = SimTK::Compound::AtomIndex(3);
// */
	// Ethane
///*
	AIx = SimTK::Compound::AtomIndex(1);
	BIx = SimTK::Compound::AtomIndex(0);
	CIx = SimTK::Compound::AtomIndex(4);
	DIx = SimTK::Compound::AtomIndex(5);
// */

	// Thermodynamics
    Real RT = world->getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;
    // Real beta = 1.0 / RT;

    // // Random number generators - not sure if I need two
    std::mt19937 randomEngine = std::mt19937();
    std::uniform_real_distribution<double> uniformRealDistribution_0_2pi =
            std::uniform_real_distribution<double>(SimTK::Zero, 2*SimTK::Pi);
    std::uniform_real_distribution<double> uniformRealDistribution_mpi_pi =
            std::uniform_real_distribution<double>((-1)*SimTK::Pi, SimTK::Pi);
    std::uniform_real_distribution<double> uniformRealDistribution =
            std::uniform_real_distribution<double>(SimTK::Zero, SimTK::One);
    std::uniform_real_distribution<double> uniformRealDistribution_m1_1 =
            std::uniform_real_distribution<double>((-1)*SimTK::One, SimTK::One);
    std::normal_distribution<> gaurand = std::normal_distribution<>(0.0, 1.0);

	// Sampler variables
	Real pe_o, pe_set, pe_n, fix_o, fix_n, fix_set;
	Real ke_proposed, ke_lastAccepted, ke_n, etot_proposed, etot_o, etot_n, etot_set;
	bool acc = false;
        // Update
        for(int k = 0; k < nofRounds; k++){ // Iterate through samples
		// Propose
		system->realize(advancedState, SimTK::Stage::Position);

		//int t = 0;
		//for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		//    TVector[t] = SetTVector[t];
		//    t++;
		//}

		pe_o = pe_set; 
		fix_o = fix_set; 

		// Initialize velocities according to the Maxwell-Boltzmann distribution
		int nu = advancedState.getNU();
		double sqrtRT = std::sqrt(RT);
		Vector V(nu);
		Vector SqrtMInvV(nu);
		
		for (int i=0; i < nu; ++i){
		    V[i] = gaurand(randomEngine);
		}
		matter->multiplyBySqrtMInv(advancedState, V, SqrtMInvV);

		SqrtMInvV *= (sqrtRT); // Set stddev according to temperature

		// Raise the temperature
		advancedState.updU() = SqrtMInvV;
		
		//advancedState.updU() = 0.0;
		//advancedState.updU()[6] = 50.0; // Arbitrary

		// Realize velocity
		system->realize(advancedState, SimTK::Stage::Velocity);

		// Store the proposed energies
		// setProposedKE(matter->calcKineticEnergy(advancedState));
		ke_proposed = matter->calcKineticEnergy(advancedState);
		etot_proposed = pe_o + fix_o + ke_proposed;

		for(int MDstep = 0; MDstep < pHMC(pSampler)->getMDStepsPerSample(); MDstep++){
			timeStepper->stepTo(advancedState.getTime() + (timestep * MDstep));
			// Calculate geometric features
			//Apos = topology.calcAtomLocationInGroundFrame(advancedState, AIx);
			//Bpos = topology.calcAtomLocationInGroundFrame(advancedState, BIx);
			//Cpos = topology.calcAtomLocationInGroundFrame(advancedState, CIx);
			//Dpos = topology.calcAtomLocationInGroundFrame(advancedState, DIx);
			//std::cout << bDihedral(Apos, Bpos, Cpos, Dpos) << std::endl;
			//std::cout << advancedState.updU() << std::endl;
		}
		// END proposing move

		// Update
    bool acc;
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

    // Get new Fixman potential
    if(useFixmanPotential){
        fix_n = (pHMC(pSampler))->calcFixman(advancedState);
    }else{
        fix_n = 0.0;
    }

    // Get new kinetic energy
    system->realize(advancedState, SimTK::Stage::Velocity);
    ke_n = matter->calcKineticEnergy(advancedState);

    // Get new potential energy
//    if ( getThermostat() == ANDERSEN ){
        pe_n = generalForceSubsystem->getMultibodySystem().calcPotentialEnergy(advancedState);
//    }
//    else{
//        pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(advancedState); // ELIZA FULL
//    }

    // Calculate total energy
    if(useFixmanPotential){
        etot_n = pe_n + ke_n + fix_n;
        etot_proposed = pe_o + ke_proposed + fix_o;
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    // Decide and get a new sample
//    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
        acc = true;
        //std::cout << " acc" << std::endl;
        //setSetTVector(advancedState);
        pe_set = pe_n;
        fix_set = fix_n;
        ke_lastAccepted = ke_n;
        etot_set = pe_set + fix_set + ke_proposed;
//    }else { // Apply Metropolis correction
//        if ((!std::isnan(pe_n)) && ((etot_n < etot_proposed) ||
//             (rand_no < exp(-(etot_n - etot_proposed) * this->beta)))) { // Accept based on full energy
//             acc = true;
//             std::cout << " acc" << std::endl;
//             setSetTVector(advancedState);
//             pe_set = pe_n;
//             fix_set = fix_n;
//             ke_lastAccepted = ke_n;
//             etot_set = pe_set + fix_set + ke_proposed;
//
//        } else { // Reject
//             acc = false;
//             std::cout << " nacc" << std::endl;
//             assignConfFromSetTVector(advancedState);
//        }
//    }
/*
    std::cout << std::setprecision(5) << std::fixed;
    std::cout << "pe_o " << pe_o << " pe_n " << pe_n
        << " pe_nB " << generalForceSubsystem->getMultibodySystem().calcPotentialEnergy(advancedState)
        << " ke_prop " << ke_proposed << " ke_n " << ke_n
        << " fix_o " << fix_o << " fix_n " << fix_n << " "
        << " etot_n " << etot_n  << " etot_proposed " << etot_proposed
        << ' '
        ;
*/
		// Calculate geometric features
		Apos = topology.calcAtomLocationInGroundFrame(advancedState, AIx);
		Bpos = topology.calcAtomLocationInGroundFrame(advancedState, BIx);
		Cpos = topology.calcAtomLocationInGroundFrame(advancedState, CIx);
		Dpos = topology.calcAtomLocationInGroundFrame(advancedState, DIx);
		std::cout << bDihedral(Apos, Bpos, Cpos, Dpos) 
			//<< ' ' << advancedState.getQ() 
			//<< std::endl
		;

		//for(int i = 7; i < advancedState.getNQ(); i++){
		//	std::cout << ' ' << advancedState.getQ()[i];
		//}
		std::cout << std::endl;
		//pHMC(pSampler)->PrintDetailedEnergyInfo(advancedState);

		//world->matter->updMobilizedBody(MobilizedBodyIndex(0));
        }


	return 0;
}
