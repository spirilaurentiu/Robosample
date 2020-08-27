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

    // // //
    CC.resize(this->K + 1, this->K + 1);

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

	SimTK::Real Ediff = E_o - E_n;
	SimTK::Real prob = exp(this->beta * Ediff);
	std::cout << "leap_prob " << E_o << " " << E_n << " " << prob << "\n" ;
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

/*** Return the antidiagonal transpose of a matrix ***/
SimTK::Matrix LAHMCSampler::reverseMatrix(SimTK::Matrix CC_loc)
{
	int M = CC_loc.nrow();
	int N = CC_loc.ncol();
	SimTK_ASSERT_ALWAYS(M == N, "Robosample: reverseMatrix function is designed for square matrices.");
	SimTK::Matrix CC_rev(M, N);

	for(int i = 0; i < M; i++){
		for(int j = 0; j < M; j++){
			CC_rev.set(M-i-1, M-j-1, CC_loc.get(i, j));
		}
	}

	return CC_rev;
}

/*** Return a submatrix of M with lesser rows and cols from the end ***/
SimTK::Matrix LAHMCSampler::extractFromTop(SimTK::Matrix CC_loc, int rowCut, int colCut)
{
	int M = CC_loc.nrow();
	int N = CC_loc.ncol();
	SimTK_ASSERT_ALWAYS(M == N, "Robosample: extractFromTop function is designed for square matrices.");
	SimTK_ASSERT_ALWAYS(rowCut <= M, "Robosample: extractFromTop: cannot cut more than the matrix size.");
	SimTK_ASSERT_ALWAYS(colCut <= N, "Robosample: extractFromTop: cannot cut more than the matrix size.");
	SimTK_ASSERT_ALWAYS((rowCut*colCut) >= 0, "Robosample: extractFromTop: indeces must both have the same sign.");

	int subM, subN;
	if(rowCut < 0){	
		subM = M + rowCut;
		subN = N + colCut;
	}else{
		subM = rowCut;
		subN = colCut;
	}

	SimTK::Matrix subCC(subM, subN);

	// Can't use operator= because Simbody resizes subCC
	for(int i = 0; i < subM; i++){
		for(int j = 0; j < subN; j++){
			subCC.set(i, j, CC_loc.get(i, j));
		}
	}

	return subCC;
}

/*** Copies different sizes matrices entries and avoids Simbody resize ***/
void LAHMCSampler::injectFromTop(const std::vector<SimTK::Real>& src, std::vector<SimTK::Real>& dest)
{
	int srcN = src.size();
	int destN = dest.size();

	// Only go up to the smallest dimension
	if(srcN < destN){ // Source is smaller
		for(int i = 0; i < srcN; i++){
			dest[i] = src[i];
		}
		
	}else{ // Destination is smaller
		for(int i = 0; i < destN; i++){
			dest[i] = src[i];
		}
	}

}

/*** Copies different sizes matrices entries and avoids Simbody resize ***/
void LAHMCSampler::injectFromTop(const SimTK::Matrix& src, SimTK::Matrix& dest)
{
	int srcM = src.nrow();
	int srcN = src.ncol();
	SimTK_ASSERT_ALWAYS(srcM == srcN, "Robosample: injectFromTop function is designed for square matrices.");
	int destM = dest.nrow();
	int destN = dest.ncol();
	SimTK_ASSERT_ALWAYS(destM == destN, "Robosample: injectFromTop function is designed for square matrices.");

	// Only go up to the smallest dimension
	if(srcM < destM){ // Source is smaller
		for(int i = 0; i < srcM; i++){
			for(int j = 0; j < srcN; j++){
				dest.set(i, j, src.get(i, j));
			}
		}
		
	}else{ // Destination is smaller
		for(int i = 0; i < destM; i++){
			for(int j = 0; j < destN; j++){
				dest.set(i, j, src.get(i, j));
			}
		}
	}

}

/*** Compute cumulative transition probabilities 
TODO Energy indeces ***/
SimTK::Matrix LAHMCSampler::leap_prob_recurse_hard(SimTK::State& someState, std::vector<SimTK::Real> Es, SimTK::Matrix CC)
{


	// Indeces between the two matrices

	int M = CC.nrow();
	int N = CC.ncol();	
	SimTK_ASSERT_ALWAYS(M == N, "Robosample: leap_prob_recurse is designed for square matrices.");
	SimTK_ASSERT_ALWAYS(M == (Es.size()), "Robosample: leap_prob_recurse: energy vector does not have the same size as C matrix.");
	std::cout << "\nBEGIN size " << M << "\n";

	// Check if already passed through this leaf
	SimTK::Real upperCorner = CC.get(M - 1, N - 1 );

	//if( upperCorner != SimTK::Infinity){
	if( upperCorner < 1){
		std::cout << "Already already visited this leaf" << std::endl;
		PrintBigMat(CC, M, N, 6, std::string(" C ") + std::to_string(M));
		std::cout << "\nEND size " << M << "\n";
		return CC;
	}

	// Boltzmann probability
	if(M == 2){
		std::cout << "Do Metropolis Hastings.\n";
		SimTK::Real p_acc = 0.0;
    		if(!std::isnan(Es[M - 1])){
				p_acc = leap_prob(someState, Es[0], Es[Es.size() - 1]);
				CC.set(0, M - 1, p_acc);
		}
		PrintBigMat(CC, M, N, 6, std::string(" C ") + std::to_string(M));
		std::cout << "\nEND index size " << M << "\n";
		return CC;
	}
	
	//std::cout << "\nDEBUG Es_1 after injectFromTop \n";
	//for(int i=0; i<Es_1.size(); i++){std::cout << Es_1[i] << " ";}
	//std::cout << "\nDEBUG ===============\n";
	SimTK::Real cum_forward, cum_reverse;
	// Forward ===============================
	std::cout << "Reenter forward";
	// Reduce size
	std::vector<SimTK::Real> Es_1(M - 1, SimTK::Infinity);
	injectFromTop(Es, Es_1);
	SimTK::Matrix CC_1 = extractFromTop(CC, -1, -1);

	// Reentry
	CC_1 = leap_prob_recurse_hard(someState, Es_1, CC_1);

	// Recover data
	//injectFromTop(Es_1, Es); // don't need recover
	injectFromTop(CC_1, CC);

	cum_forward = CC_1.get(0, CC_1.ncol() - 1);
	std::cout << " cum_fwd = " << cum_forward << std::endl;

	// Backward =============================
	std::cout << "Reenter reverse.";

	// Reverse
	SimTK::Matrix CC_rev(CC.nrow(), CC.ncol());
	CC_rev = reverseMatrix(CC);
	std::vector<SimTK::Real> Es_rev = Es; // deep copy
	std::reverse(Es_rev.begin(), Es_rev.end());

	// Reduce
	std::vector<SimTK::Real> Es_rev_1(M - 1, SimTK::Infinity);
	injectFromTop(Es_rev, Es_rev_1);
	SimTK::Matrix CC_rev_1 = extractFromTop(CC_rev, -1, -1);

	// Reentry
	CC_rev_1 = leap_prob_recurse_hard(someState, Es_rev_1, CC_rev_1);
	cum_reverse = CC_rev_1.get(0, M - 2);

        // Recover reversed
        //0000000000000000000000 Cl_rev[:-1,:-1,:] = Cl_rev_1 000000000000000000
	//injectFromTop(Es_rev_1, Es_rev); // don't need recover
	injectFromTop(CC_rev_1, CC_rev);

        // Recover = reverse the reversed
        //0000000000000000000000 C = Cl_rev[::-1,::-1,:] 000000000000000000000000
	//std::vector<SimTK::Real> Es_rev_rev = Es; // deep copy // don't need recover
	//std::reverse(Es_rev_rev.begin(), Es_rev_rev.end()); // don't need recover

	SimTK::Matrix CC_rev_rev(CC.nrow(), CC.ncol());
	CC_rev_rev = reverseMatrix(CC_rev);

	//injectFromTop(Es_rev_rev, Es); // don't need recover
	injectFromTop(CC_rev_rev, CC);

	PrintBigMat(CC, M, N, 6, std::string(" C ") + std::to_string(M));
	std::cout << " cum_rev = " << cum_reverse << std::endl;

	// Eq. 25
	std::cout << "Do LAHMC probability." << "\n";
	SimTK::Real Ediff;
	Ediff = Es[0] - Es[M - 1];
	SimTK::Real start_state_ratio = exp(this->beta * Ediff);

	SimTK::Real prob;
	prob = std::min(1 - cum_forward, start_state_ratio * (1 - cum_reverse));
	std::cout << "E0 E1 start_state_ratio cum_fwd cum_rev prob " 
		<< Es[0] << " " << Es[Es.size() - 1] 
		<< " " << start_state_ratio << " " << cum_forward << " " << cum_reverse << " " 
		<< prob << std::endl;

	SimTK::Real cumu = cum_forward + prob;
	CC.set(0, M - 1, cumu);

	PrintBigMat(CC, M, N, 6, std::string(" C ") + std::to_string(M));
	std::cout << "\nEND index size " << "\n";
	return CC;
}

// Good pieces of code
//Es_1.insert(Es_1.begin(), Es.begin(), Es.begin() + M - 1);

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
	//leap_prob_recurse(someState, 0, this->K, true);
	leap_prob_recurse_hard(someState, this->etot_ns, this->C);

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


