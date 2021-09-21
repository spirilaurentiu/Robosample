#include "Robo.hpp"
#include "HMCSampler.hpp"
#include "World.hpp"

using namespace SimTK;

int main(int argc, char **argv)
{
	if(argc < 7){
		std::cout << "Error: not enough parameters to run. Usage:\n";
		std::cout << argv[0] << " nofRounds (N)FP (N)FT (no)visual VV/HMC mdsteps\n"; 
		return 1;
	}

	Real timestep = 0.004;

	bool wantVisual = false;
	if( std::string(argv[4]) == "visual"){
		wantVisual = true;
	}
	Real visualFreq = timestep;
	int nofRounds = std::stoi(argv[1]);
	int noMDStepsPerSample = 0; 
	std::cout << "VISUAL " << wantVisual << std::endl;

	int worldIndex = 0;
	int requestedNofMols = 1;
	World *world = new World(worldIndex, requestedNofMols, wantVisual, visualFreq);

	bool useFixmanPotential = false;
	bool useFixmanTorque = false;

	// SimTK::CompoundSystem *compoundSystem = world->compoundSystem;
	SimTK::SimbodyMatterSubsystem *matter = world->matter.get();
	SimTK::GeneralForceSubsystem *generalForceSubsystem = world->forces.get();
	SimTK::DuMMForceFieldSubsystem *dumm = world->forceField.get();
	SimTK::TimeStepper *timeStepper = world->ts.get();
	const SimTK::System *system = &(matter->getSystem());

	readAmberInput amberReader;
	amberReader.readAmberFiles("lin4/ligand.rst7", "lin4/ligand.prmtop");
	//std::string externalJointType = "Ball";
	std::string externalJointType = "Weld";
	world->AddMolecule(&amberReader, "lin4/ligand.rb", "lin4/ligand.flex.middletd", "R0", "0", externalJointType);
	//amberReader.readAmberFiles("etanMob/ligand.rst7", "etanMob/ligand.prmtop");
	//world->AddMolecule(&amberReader, "etanMob/ligand.rb", "etanMob/ligand.flex", "R0", "0", "Free");

	// Using OpenMM acceleration
	//world->updForceField()->setUseOpenMMAcceleration(true);

	// Set Lennard-Jones mixing rule
	world->updForceField()->setVdwMixingRule(DuMMForceFieldSubsystem::LorentzBerthelot);

	// Link the Compounds to Simbody System for all Worlds
	world->modelTopologies("Free");

	// Add Fixman torque (Additional ForceSubsystem) if required
	if(std::string(argv[3]) == "FT"){
		useFixmanTorque = true;
		world->addFixmanTorque();
		std::cout << "Fixman torque added." << std::endl;
	}

	// Set worlds force field scale factors
	//world->setGlobalForceFieldScaleFactor(0.0);
	world->updForceField()->setBondStretchGlobalScaleFactor(0);
	world->updForceField()->setBondBendGlobalScaleFactor(0);
	world->updForceField()->setBondTorsionGlobalScaleFactor(0);
	world->updForceField()->setAmberImproperTorsionGlobalScaleFactor(0);
  
	world->updForceField()->setVdw12ScaleFactor(0);
	world->updForceField()->setVdw13ScaleFactor(0);
	world->updForceField()->setVdw14ScaleFactor(0);
	world->updForceField()->setVdw15ScaleFactor(0);
	world->updForceField()->setVdwGlobalScaleFactor(0);
  
	world->updForceField()->setCoulomb12ScaleFactor(0);
	world->updForceField()->setCoulomb13ScaleFactor(0);
	world->updForceField()->setCoulomb14ScaleFactor(0);
	world->updForceField()->setCoulomb15ScaleFactor(0);
	world->updForceField()->setCoulombGlobalScaleFactor(0);

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

	//if(std::string(argv[5]) == "VV"){
	(pHMC(pSampler))->setAlwaysAccept(true);
	//}else if(std::string(argv[5]) == "HMC"){
	//	(pHMC(pSampler))->setAlwaysAccept(false);
	//}else{
	//	std::cout << "Unknown sampler type: " << argv[5] << "\n" ;
	//	return 1;
	//}

	// Set thermostats
	pHMC(pSampler)->setThermostat(ThermostatName::ANDERSEN);
	pHMC(pSampler)->setBoostTemperature(1);
	pHMC(pSampler)->setBoostMDSteps(1);

	// Activate Fixman potential if needed
	if(std::string(argv[2]) == "FP"){
		useFixmanPotential = true;
		pHMC(pSampler)->useFixmanPotential();
		std::cout << "Fixman potential added." << std::endl;
	}

	// Set thermodynamics
	world->setTemperature(625);

	noMDStepsPerSample = std::stoi(argv[6]);
	std::cout << "Using MDSTEPS " << noMDStepsPerSample << std::endl;
	pHMC(pSampler)->setMDStepsPerSample(noMDStepsPerSample);

	// Set the seeds for reproducibility
	pSampler->setSeed(1);

	// Initialize samplers
	SimTK::State& advancedState = world->integ->updAdvancedState();
	pSampler->initialize( advancedState );

	// Reinitialize current sampler
	pSampler->reinitialize(advancedState);

	printf("NU: %d \n", advancedState.getNU());

	Topology& topology = world->updTopology(0);

	// Geometry
	SimTK::Compound::AtomIndex AIx, BIx, CIx, DIx;
	SimTK::Vec3 Apos, Bpos, Cpos, Dpos;

///*
	// Lin4
	AIx = SimTK::Compound::AtomIndex(1);
	BIx = SimTK::Compound::AtomIndex(0);
	CIx = SimTK::Compound::AtomIndex(2);
	DIx = SimTK::Compound::AtomIndex(3);
// */
	// Ethane
/*
	AIx = SimTK::Compound::AtomIndex(1);
	BIx = SimTK::Compound::AtomIndex(0);
	CIx = SimTK::Compound::AtomIndex(4);
	DIx = SimTK::Compound::AtomIndex(5);
// */

	// Thermodynamics
	Real RT = world->getTemperature() * SimTK_BOLTZMANN_CONSTANT_MD;
	std::cout << "RT= " << RT << std::endl;
	// Real beta = 1.0 / RT;

	std::mt19937 randomEngine = std::mt19937();

	// Random seed
	std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
	std::chrono::system_clock::duration dtn = tp.time_since_epoch();
	int64_t clockSeed = dtn.count();
	randomEngine.seed(clockSeed);

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
	Real pe_o, pe_set, pe_n, fix_o, fix_n, fix_set, extfix_o, extfix_n, extfix_set;
	Real ke_proposed, ke_lastAccepted, ke_n, etot_proposed, etot_o, etot_n, etot_set;
	bool acc = false;

	// Main loop
	for(int k = 0; k < nofRounds; k++){ // Iterate through samples

		// Realize position
		system->realize(advancedState, SimTK::Stage::Position);

		//int t = 0;
		//for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		//    TVector[t] = SetTVector[t];
		//    t++;
		//}

		// Set old quantities
		pe_o = pe_set; 
		fix_o = fix_set; 
		extfix_o = extfix_set; 

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

		advancedState.updU() = SqrtMInvV;
		
		//advancedState.updU() = 0.0;
		//advancedState.updU()[6] = 50.0; // Arbitrary

		// Realize velocity
		system->realize(advancedState, SimTK::Stage::Velocity);

		// Store the proposed energies
		// setProposedKE(matter->calcKineticEnergy(advancedState));
		ke_proposed = matter->calcKineticEnergy(advancedState);
		etot_proposed = pe_o + fix_o + extfix_o + ke_proposed;

		// Integrate trajectory
		for(int MDstep = 0; MDstep < pHMC(pSampler)->getMDStepsPerSample(); MDstep++){

			// Step
			timeStepper->stepTo(advancedState.getTime() + (timestep * MDstep));

			///////////////////////
			// Compute external Fixman MBAT term
			const SimTK::MobilizedBody& mobod1 = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(1));
			const SimTK::MobilizedBody* mobod1ptr = &mobod1;

			SimTK::Real pitch = 0.0;
			SimTK::Real FP = 0.0;
			SimTK::Real FT = 0.0;

			if((externalJointType == "Ball") || (externalJointType == "Free")){
				SimTK::Real w = 0;
				SimTK::Real x = 0;
				SimTK::Real y = 0;
				SimTK::Real z = 0;
				if(externalJointType == "Free"){
					const SimTK::Vec7& extQ = ((SimTK::MobilizedBody::Free *)mobod1ptr)->getQ(advancedState);
					w = extQ[0];
					x = extQ[1];
					y = extQ[2];
					z = extQ[3];
				}else if(externalJointType == "Ball"){
					const SimTK::Vec4& extQ = ((SimTK::MobilizedBody::Ball *)mobod1ptr)->getQ(advancedState);
					w = extQ[0];
					x = extQ[1];
					y = extQ[2];
					z = extQ[3];
				}
		
				// Normalize quaternion
				SimTK::Real quatnorm = std::sqrt((w*w) + (x*x) + (y*y) + (z*z));
				w /= quatnorm; x /= quatnorm; y /= quatnorm; z /= quatnorm;
				SimTK::Quaternion quat(w, x, y, z);
		
				// Compute pitch Euler angle 
				SimTK::Real sinPitch = 2*((w*y) - (z*x)); // SAME and FASTER
				SimTK::Real pitch = std::asin(sinPitch); 
				SimTK::Real cotPitch = 1 / std::tan(pitch);
		
				SimTK::Real FP = -0.5 * RT * std::log(sinPitch * sinPitch);
				SimTK::Real FT = RT * cotPitch;
			}else{
				pitch = 0.0;
				FP = 0.0;
				FT = 0.0;
			}
	
			// Print 
			std::cout << "test pitch FP FT dihedral " << pitch << " " << FP << " " << FT ;
			
			///////////////////////
			///////////////////////
			///////////////////////

			// Calculate geometric features
			Apos = topology.calcAtomLocationInGroundFrame(advancedState, AIx);
			Bpos = topology.calcAtomLocationInGroundFrame(advancedState, BIx);
			Cpos = topology.calcAtomLocationInGroundFrame(advancedState, CIx);
			Dpos = topology.calcAtomLocationInGroundFrame(advancedState, DIx);
			std::cout << " " << bDihedral(Apos, Bpos, Cpos, Dpos) << std::endl;
			//std::cout << advancedState.updU() << std::endl;

		}
		// END Integrate trajectory

		bool acc;
		SimTK::Real rand_no = uniformRealDistribution(randomEngine);

		// Get new Fixman potential
		if(useFixmanPotential){
			fix_n = (pHMC(pSampler))->calcFixman(advancedState);
			extfix_n = 0.5 * RT * topology.calcLogSineSqrGamma2(advancedState);
		}else{
			fix_n = 0;
			extfix_n = 0;
		}

		// Get new kinetic energy
		system->realize(advancedState, SimTK::Stage::Velocity);
		ke_n = matter->calcKineticEnergy(advancedState);

		// Get new potential energy
		pe_n = generalForceSubsystem->getMultibodySystem().calcPotentialEnergy(advancedState);
		//pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(advancedState); // ELIZA FULL

		// Calculate total energy
		if(useFixmanPotential){
			etot_n = pe_n + ke_n + fix_n + extfix_n;
			etot_proposed = pe_o + ke_proposed + fix_o + extfix_n;
		}else{
			etot_n = pe_n + ke_n;
			etot_proposed = pe_o + ke_proposed;
		}

		// Decide and get a new sample
		if(std::string(argv[5]) == "VV"){
			acc = true;
			//std::cout << " acc" << std::endl;
			//setSetTVector(advancedState);
			pe_set = pe_n;
			fix_set = fix_n;
			extfix_set = extfix_n;
			ke_lastAccepted = ke_n;
	
			etot_set = pe_set + fix_set + extfix_set + ke_proposed;
	
		}else if(std::string(argv[5]) == "HMC"){
			if ((!std::isnan(pe_n)) && ((etot_n < etot_proposed) ||
				(rand_no < exp(-(etot_n - etot_proposed) * (1 / RT))))) { // Accept based on full energy
				acc = true;
				std::cout << " acc" << std::endl;
				(pHMC(pSampler))->setSetTVector(advancedState);
				pe_set = pe_n;
				fix_set = fix_n;
				extfix_set = extfix_n;
				ke_lastAccepted = ke_n;
				etot_set = pe_set + fix_set + extfix_set + ke_proposed;
	
			} else { // Reject
				acc = false;
				std::cout << " nacc" << std::endl;
				(pHMC(pSampler))->assignConfFromSetTVector(advancedState);
			}
	    	}
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
		std::cout << "round dihedral " << bDihedral(Apos, Bpos, Cpos, Dpos) 
			//<< ' ' << advancedState.getQ() 
			//<< std::endl
		;

		//for(int i = 7; i < advancedState.getNQ(); i++){
		//	std::cout << ' ' << advancedState.getQ()[i];
		//}
		std::cout << std::endl;
		//pHMC(pSampler)->PrintDetailedEnergyInfo(advancedState);

		//world->matter->updMobilizedBody(MobilizedBodyIndex(0));
        } // END MAIN LOOP


	return 0;
}
