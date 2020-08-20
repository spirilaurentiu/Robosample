/**@file
Implementation of LAHMCSampler class. **/

#include "LAHMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
LAHMCSampler::LAHMCSampler(SimTK::CompoundSystem *argCompoundSystem
                                     ,SimTK::SimbodyMatterSubsystem *argMatter
                                     ,SimTK::Compound *argResidue
                                     ,SimTK::DuMMForceFieldSubsystem *argDumm
                                     ,SimTK::GeneralForceSubsystem *argForces
                                     ,SimTK::TimeStepper *argTimeStepper
				     ,unsigned int Kext
                                     )
    : Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    HMCSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
{
    this->useFixman = false;  
    this->fix_n = this->fix_o = 0.0;
    this->logSineSqrGamma2_n = this->logSineSqrGamma2_o = 0.0;
    this->residualEmbeddedPotential = 0.0;
    nofSamples = 0;
    this->alwaysAccept = false;
    this->timestep = 0.002; // ps
    this->temperature = 0.0;
    this->boostT = this->temperature;
    MDStepsPerSample = 0;
    proposeExceptionCaught = false;
    adaptTimestep = false;


    this->K = Kext;
    for(unsigned int i = 0; i < K + 1; i++){
    	pe_ns.push_back(SimTK::Infinity);
    	fix_ns.push_back(SimTK::Infinity);
    	logSineSqrGamma2_ns.push_back(SimTK::Infinity);
    	ke_ns.push_back(SimTK::Infinity);
    	etot_ns.push_back(SimTK::Infinity);
    }

    C.resize(this->K + 1, this->K + 1);
    Ctau.resize(this->K + 1, this->K + 1);
    P.resize(this->K + 1, this->K + 1); // Anti-diagonal identity

    //int k = -1;for(int i=0;i<someK;i++){for(int j=0;j<someK;j++){k++;C.set(i,j,k);}}
    int k = 1;
    for(int i = 0; i < this->K + 1; i++){
        for(int j = 0; j < this->K + 1; j++){
            k++;
            //C.set(i, j, 1.0 * k);
            //C.set(i, j, -1.0);

            C.set(i, j, SimTK::Infinity);
            Ctau.set(i, j, SimTK::Infinity);

            if(i == (this->K) - j){
                P.set(i, j, 1);
            }else{
                P.set(i, j, 0);
            }
        }
    }
    
    //Ctau = P * (C * P);
}

/** Destructor **/
LAHMCSampler::~LAHMCSampler()
{
}

/** Calculate sqrt(M) using Eigen. For debug purposes. **/
//void LAHMCSampler::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
//{
//    assert("!Not implemented");
//}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular matrix and it is computed 
multipling a set of orthonormal vectors with the sqrt(MInv). **/
//void LAHMCSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
//{
//    int nu = someState.getNU();
//    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");
//
//    SimTK::Vector V(nu);
//    SimTK::Vector SqrtMInvV(nu);
//
//    // Zero the matrix and the temporary vector
//    for (int i=0; i < nu; ++i){
//        V[i] = 0;
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = 0;
//        }
//    }
//
//    // Calculate the values inside the matrix
//    for (int i=0; i < nu; ++i){
//        V[i] = 1;
//        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = SqrtMInvV[j];
//        }
//        V[i] = 0;
//    }
//
//
//}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
This is lower triangular matrix and it is computed by multipling a set of
 orthonormal vectors with the sqrt(MInv) and transpose it. **/
//void LAHMCSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
//{
//    int nu = someState.getNU();
//    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");
//
//    SimTK::Vector V(nu);
//    SimTK::Vector SqrtMInvV(nu);
//
//    // Zero the matrix and the temporary vector
//    for (int i=0; i < nu; ++i){
//        V[i] = 0;
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(i, j) = 0;
//        }
//    }
//
//    // Calculate the values inside the matrix
//    for (int i=0; i < nu; ++i){
//        V[i] = 1;
//        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
//        for (int j=0; j < nu; ++j){
//            SqrtMInv(j, i) = SqrtMInvV[j];
//        }
//        V[i] = 0;
//    }
//
//}

/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void LAHMCSampler::setLastAcceptedKE(SimTK::Real inpKE)
{
    this->ke_lastAccepted = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void LAHMCSampler::setProposedKE(SimTK::Real inpKE)
{
    this->ke_proposed = inpKE;
}

/** Get/set the TimeStepper that manages the integrator **/
//const SimTK::TimeStepper * LAHMCSampler::getTimeStepper(void)
//{
//    return timeStepper;
//}

//SimTK::TimeStepper * LAHMCSampler::updTimeStepper(void)
//{
//    return timeStepper;
//}

//void LAHMCSampler::setTimeStepper(SimTK::TimeStepper * someTimeStepper)
//{
//    timeStepper = someTimeStepper;
//}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void LAHMCSampler::initialize(SimTK::State& someState )
{
    // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    int nu = someState.getNU();

    // Randomize configuration
    //if(randomizeConformation == true){
    //    system->realize(someState, SimTK::Stage::Position);
    //    int nq = someState.getNQ();
    //    SimTK::Vector QV(nq);
    //    for (int j=7; j < nq; ++j){
    //        QV[j] = uniformRealDistribution_mpi_pi(randomEngine);
    //    }
    //    someState.updQ() = QV;
    //}
    //

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
    setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
//r    this->useFixman = argUseFixman;
    if(useFixman){
        std::cout << "Look Ahead Hamiltonian Monte Carlo sampler: using Fixman potential." << std::endl;
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2( ((Topology *)residue)->calcLogSineSqrGamma2(someState));
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2(0.0);
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }

    // Initialize velocities to temperature
    // TODO Shouldn't be here
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
    this->etot_set = this->etot_proposed;

  
}

/** Same as initialize **/
//r void LAHMCSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
void LAHMCSampler::reinitialize(SimTK::State& someState)
{
     // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Set the simulation temperature
//r    setTemperature(argTemperature); // Needed for Fixman

    // Store the configuration
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

//std::cout << "reinitialize "
//<< dumm->CalcFullPotEnergyIncludingRigidBodies(someState) << std::endl;

    // Store potential energies
    //setOldPE(getPEFromEvaluator(someState));
    setSetPE(getOldPE());

    // Store Fixman potential
    if(useFixman){
        setOldFixman(calcFixman(someState));
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2( ((Topology *)residue)->calcLogSineSqrGamma2(someState));
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }else{
        setOldFixman(0.0);
        setSetFixman(getOldFixman());

        setOldLogSineSqrGamma2(0.0);
        setSetLogSineSqrGamma2(getOldLogSineSqrGamma2());
    }

    // Initialize velocities to temperature
    // TODO Shouldn't be here
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int j=0; j < nu; ++j){
        V[j] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
    SqrtMInvV *= sqrtRT; // Set stddev according to temperature
    someState.updU() = SqrtMInvV;
    system->realize(someState, SimTK::Stage::Velocity);

    // Store kinetic energies
    setProposedKE(matter->calcKineticEnergy(someState));
    setLastAcceptedKE(getProposedKE());

    // Store total energies
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
    this->etot_set = this->etot_proposed;

}

/** Store configuration **/
void LAHMCSampler::storeOldConfigurationAndPotentialEnergies(SimTK::State& someState){ // func
    system->realize(someState, SimTK::Stage::Position);

    int t = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        TVector[t] = SetTVector[t];
        t++;
    }

    // TODO: change the names from Old to Proposed and Set to lastAccepted
    pe_o = pe_set;
    fix_o = fix_set;
    logSineSqrGamma2_o = logSineSqrGamma2_set;
} // func
    
/** Initialize velocities according to the Maxwell-Boltzmann
distribution.  Coresponds to R operator in LAHMC **/
void LAHMCSampler::initializeVelocities(SimTK::State& someState){
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }
    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

    SqrtMInvV *= (sqrtRT); // Set stddev according to temperature

    // Raise the temperature
    someState.updU() = SqrtMInvV;

    // Realize velocity
    system->realize(someState, SimTK::Stage::Velocity);
}
    
/** Store the proposed energies **/
void LAHMCSampler::calcProposedKineticAndTotalEnergy(SimTK::State& someState){
    // setProposedKE(matter->calcKineticEnergy(someState));
    this->ke_proposed = matter->calcKineticEnergy(someState);
    
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
}
    
/** Apply the L operator **/
void LAHMCSampler::integrateTrajectory(SimTK::State& someState){    
    try {
        this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
    }catch(const std::exception&){
        proposeExceptionCaught = true;
        //std::cout << "HMC_stepTo_or_realizePosition_threw_exception";
        int i = 0;
        for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
            const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
            mobod.setQToFitTransform(someState, SetTVector[i]);
            i++;
        }
        system->realize(someState, SimTK::Stage::Position);
    }

}

/** Store new configuration and energy terms**/
void LAHMCSampler::calcNewConfigurationAndEnergies(SimTK::State& someState, int k)
{
    // Get new Fixman potential
    if(useFixman){
        fix_ns[k] = calcFixman(someState);
        logSineSqrGamma2_ns[k] = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }else{
        fix_ns[k] = 0.0;
        logSineSqrGamma2_ns[k] = 0.0;
    }

    // Get new kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    ke_ns[k] = matter->calcKineticEnergy(someState);

    // Get new potential energy
    pe_ns[k] = forces->getMultibodySystem().calcPotentialEnergy(someState);
    logSineSqrGamma2_ns[k] = ((Topology *)residue)->calcLogSineSqrGamma2(someState);

    // Calculate total energy
    if(useFixman){
        etot_ns[k] = pe_ns[k] + ke_ns[k] + fix_ns[k] - (0.5 * RT * logSineSqrGamma2_ns[k]);
        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
    }else{
        etot_ns[k] = pe_ns[k] + ke_ns[k];
        etot_proposed = pe_o + ke_proposed;
    }

    pe_n = pe_ns[k];
    fix_n = fix_ns[k];
    logSineSqrGamma2_n = logSineSqrGamma2_ns[k];
    ke_n = ke_ns[k];
    etot_n = etot_ns[k];

    PrintDetailedEnergyInfo(someState);
}

/*** Set C matrices entries to 0 ***/
void LAHMCSampler::resetCMatrices(void){
	for(int i = 0; i < this->K + 1; i++){
		for(int j = 0; j < this->K + 1; j++){
			C.set(i, j, SimTK::Infinity);
			Ctau.set(i, j, SimTK::Infinity);
		}
	}
}

/*** Set C matrices entry ***/
void LAHMCSampler::setCEntry(int i, int j, SimTK::Real entry){
	C.set(i, j, entry);
	Ctau.set(K - i, K - j, entry);
}

/*** Set C matrices entry ***/
void LAHMCSampler::setCtauEntry(int i, int j, SimTK::Real entry){
	Ctau.set(i, j, entry);
	C.set(K - i, K - j, entry);
}

/*** Metropolis-Hastings acception-rejection criterion  ***/
SimTK::Real LAHMCSampler::leap_prob(SimTK::State& someState, SimTK::Real E_o, SimTK::Real E_n)
{
	std::cout << "leap_prob " << E_o << " " << E_n ;

	SimTK::Real Ediff = E_o - E_n;
	SimTK::Real prob = exp(this->beta * Ediff);
	//if(prob > 1){
	//	return 1;
	//}else{
		return prob;
	//}
}

/*** Convert C indeces to Ctau indeces ***/
void LAHMCSampler::C_to_Ctau_Indeces(int C_i, int C_j, int &Ctau_i, int &Ctau_j, int currSize){
	Ctau_j = this->K - 1 - C_i;
	Ctau_i = this->K - 1 - C_j;
}

/*** Convert Ctau indeces to C indeces ***/
void LAHMCSampler::Ctau_to_C_Indeces(int Ctau_i, int Ctau_j, int &C_i, int &C_j, int currSize){
	int offset = (this->K+1) - currSize;
	C_j = offset + (currSize - 1) - Ctau_i;
	C_i = offset + (currSize - 1) - Ctau_j;
	
}

/*** Compute cumulative transition probabilities 
TODO Energy indeces ***/
SimTK::Real LAHMCSampler::leap_prob_recurse(SimTK::State& someState, int firstIx, int secondIx, bool dir_fwd)
{
	int currSize = secondIx - firstIx + 1;

	SimTK::Real Ediff;

	// Indeces between the two matrices
	int revFirstIx, revSecondIx;
	Ctau_to_C_Indeces(firstIx, secondIx, revFirstIx, revSecondIx, currSize);
	std::cout << "\nBEGIN Z_chain shape " << currSize << " " << dir_fwd 
		<< " " << firstIx << " " << secondIx 
		<< " " << revFirstIx << " " << revSecondIx << "\n";

	int upperCorner_i = firstIx;
	int upperCorner_j = secondIx;
	std::cout << upperCorner_i << " " << upperCorner_j << std::endl;

	// Check if already passed through this leaf
	SimTK::Real upperCorner;
	if(dir_fwd){
		upperCorner = C.get(upperCorner_i, upperCorner_j);
	}else{
		upperCorner = Ctau.get(upperCorner_i, upperCorner_j);
	}

	//if( upperCorner != SimTK::Infinity){
	if( upperCorner < 1){
		std::cout << "Already already visited this leaf" << std::endl;
		PrintBigMat(C, this->K + 1, this->K + 1, 6, std::string("C"));
		std::cout << "\nEND index gym size " << currSize << " " << dir_fwd 
			<< " " << firstIx << " " << secondIx 
			<< " " << revFirstIx << " " << revSecondIx << "\n";
		return upperCorner;
	}

	// Boltzmann probability
	if(currSize == 2){
		std::cout << "Do Metropolis Hastings.\n";
		
		SimTK::Real p_acc = 0.0;
    		if(!std::isnan(pe_ns[upperCorner_j])){
			if(dir_fwd){
				p_acc = leap_prob(someState, etot_ns[firstIx], etot_ns[secondIx]);
				setCEntry(upperCorner_i, upperCorner_j, p_acc);
				std::cout << " " << firstIx << " " << secondIx ;
			}else{
				p_acc = leap_prob(someState, etot_ns[revSecondIx], etot_ns[revFirstIx]);
				setCtauEntry(upperCorner_i, upperCorner_j, p_acc);
				std::cout << " " << revSecondIx << " " << revFirstIx ;
			}
			std::cout << " " << p_acc << std::endl;
		}
		PrintBigMat(C, this->K + 1, this->K + 1, 6, std::string("C"));
		std::cout << "\nEND index gym size " << currSize << " " << dir_fwd 
			<< " " << firstIx << " " << secondIx 
			<< " " << revFirstIx << " " << revSecondIx << "\n";
		return p_acc;
	}
	
	// Forward
	std::cout << "Reenter forward" << currSize << " " << dir_fwd;
	SimTK::Real cum_forward, cum_reverse;
	if(dir_fwd){
		std::cout << " ixs: " << firstIx << " " << secondIx << " to " << firstIx << " " << secondIx -1 << "\n";
		cum_forward = leap_prob_recurse(someState, firstIx, secondIx - 1, true);
	}else{
		std::cout << " ixs: " << firstIx << " " << secondIx << " to " << revFirstIx << " " << revSecondIx -1 << "\n";
		cum_forward = leap_prob_recurse(someState, revFirstIx, revSecondIx - 1, false);
	}
	std::cout << "cum_fwd = " << cum_forward << std::endl;
	// Backward
	std::cout << "Reenter reverse." << currSize << " " << dir_fwd; 
	if(dir_fwd){
		std::cout << " ixs: " << firstIx << " " << secondIx << " to " << revFirstIx << " " << revSecondIx -1 << "\n";
		cum_reverse = leap_prob_recurse(someState, revFirstIx, revSecondIx - 1, false);
	}else{
		std::cout << " ixs: " << firstIx << " " << secondIx << " to " << revFirstIx << " " << revSecondIx -1 << "\n";
		cum_reverse = leap_prob_recurse(someState, revFirstIx, revSecondIx - 1, true);
	}

	// Eq. 25
	std::cout << "Do LAHMC probability." << currSize << " " << dir_fwd 
		<< " " << firstIx << " " << secondIx 
		<< " " << revFirstIx << " " << revSecondIx << "\n";
	if(dir_fwd){
		Ediff = etot_ns[firstIx] - etot_ns[secondIx];
	}else{
		Ediff = etot_ns[revSecondIx] - etot_ns[revFirstIx];
	}
	SimTK::Real start_state_ratio = exp(this->beta * Ediff);

	SimTK::Real prob;
	prob = std::min(1 - cum_forward, start_state_ratio * (1 - cum_reverse));
	if(dir_fwd){
		std::cout << "E0 E1 start_state_ratio cum_fwd cum_rev prob " 
			<< etot_ns[firstIx] << " " << etot_ns[secondIx] << " " << start_state_ratio << " " << cum_forward << " " << cum_reverse << " " << prob << std::endl;
	}else{
		std::cout << "E0 E1 start_state_ratio cum_fwd cum_rev prob " 
			<< etot_ns[revSecondIx] << " " << etot_ns[revFirstIx] << " " << start_state_ratio << " " << cum_forward << " " << cum_reverse << " " << prob << std::endl;
	}

	SimTK::Real cumu = cum_forward + prob;

	if(dir_fwd){
		setCEntry(upperCorner_i, upperCorner_j, cumu);
	}else{
		setCtauEntry(upperCorner_i, upperCorner_j, cumu);
	}

	PrintBigMat(C, this->K + 1, this->K + 1, 6, std::string("C"));
	std::cout << "\nEND index gym size " << currSize << " " << dir_fwd 
		<< " " << firstIx << " " << secondIx 
		<< " " << revFirstIx << " " << revSecondIx << "\n";
	return cumu;
}


/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void LAHMCSampler::propose(SimTK::State& someState)
{

	PrintBigMat(C, this->K + 1, this->K + 1, 2, std::string("C"));
	std::cout << "===========\n";
	PrintBigMat(C.updBlock(0, 0, this->K + 1-1, this->K + 1-1), this->K + 1-1, this->K + 1-1, 2, "C[:-1,:-1]"); // 0 <= i < m, 0 <= j < n
	std::cout << "===========\n";
	PrintBigMat(Ctau.updBlock(0, 0, this->K + 1-1, this->K + 1-1), this->K + 1-1, this->K + 1-1, 2, "C[:0:-1,:0:-1]");
	std::cout << "===========\n";

    storeOldConfigurationAndPotentialEnergies(someState);

    initializeVelocities(someState);

    calcProposedKineticAndTotalEnergy(someState);

    // First set the starting point
    calcNewConfigurationAndEnergies(someState, 0);

    // Apply leap operator K times
    for(int k = 1; k <= K; k++){
        integrateTrajectory(someState);
        calcNewConfigurationAndEnergies(someState, k);
    }

	etot_ns[0] = 2;
	etot_ns[1] = 3;
	etot_ns[2] = 5;
	etot_ns[3] = 8;
	etot_ns[4] = 12;

	resetCMatrices();
	leap_prob_recurse(someState, 0, this->K, true);

}

/** Store new configuration and energy terms**/
void LAHMCSampler::setSetConfigurationAndEnergiesToNew(SimTK::State& someState)
{
    setSetTVector(someState);
    pe_set = pe_n;
    fix_set = fix_n;
    logSineSqrGamma2_set = logSineSqrGamma2_n;
    ke_lastAccepted = ke_n;
    etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;
}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
bool LAHMCSampler::update(SimTK::State& someState, SimTK::Real newBeta)
{

    // Declare variables
    bool acc;
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);


    // Decide and get a new sample
    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
        acc = true;
        std::cout << " acc" << std::endl;

	setSetConfigurationAndEnergiesToNew(someState);
	
        //setBeta(newBeta);

        ++acceptedSteps;
        acceptedStepsBuffer.push_back(1);
        acceptedStepsBuffer.pop_front();
    }else{ // Apply Metropolis correction
        if ( (proposeExceptionCaught == false) &&
                (!std::isnan(pe_n)) && ((etot_n < etot_proposed) ||
             (rand_no < exp(-(etot_n - etot_proposed) * this->beta)))) { // Accept based on full energy
             acc = true;
             std::cout << " acc" << std::endl;

	     setSetConfigurationAndEnergiesToNew(someState);

             ++acceptedSteps;
             acceptedStepsBuffer.push_back(1);
             acceptedStepsBuffer.pop_front();

        } else { // Reject
             acc = false;
             std::cout << " nacc" << std::endl;

             assignConfFromSetTVector(someState);

             proposeExceptionCaught = false;
             acceptedStepsBuffer.push_back(0);
             acceptedStepsBuffer.pop_front();
        }
    }

    ++nofSamples;
    return acc;

}

int LAHMCSampler::getMDStepsPerSample() const {
    return MDStepsPerSample;
}

void LAHMCSampler::setMDStepsPerSample(int mdStepsPerSample) {
    MDStepsPerSample = mdStepsPerSample;
}

/** Print detailed energy information **/
void LAHMCSampler::PrintDetailedEnergyInfo(SimTK::State& someState)
{

    std::cout << std::setprecision(5) << std::fixed;
    std::cout << "pe_o " << pe_o << " pe_n " << pe_n
        << " pe_nB " << getPEFromEvaluator(someState)
        << " ke_prop " << ke_proposed << " ke_n " << ke_n
        << " fix_o " << fix_o << " fix_n " << fix_n << " "
        << " logSineSqrGamma2_o " << logSineSqrGamma2_o << " logSineSqrGamma2_n " << logSineSqrGamma2_n << " "
        //<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
        << " ts " << timestep  << " exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
        << " etot_n " << etot_n  << " etot_proposed " << etot_proposed
        //<< std::endl
        ;
}


