/**@file
Implementation of LAHMCSampler class. **/

#include "LAHMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
LAHMCSampler::LAHMCSampler(World *argWorld,
    SimTK::CompoundSystem *argCompoundSystem,
	SimTK::SimbodyMatterSubsystem *argMatter,
	//SimTK::Compound *argResidue,
	std::vector<Topology> &argTopologies,
	SimTK::DuMMForceFieldSubsystem *argDumm,
	SimTK::GeneralForceSubsystem *argForces,
	SimTK::TimeStepper *argTimeStepper,
    unsigned int Kext) :
        Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper),
        //MonteCarloSampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper),
        HMCSampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
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
    shouldAdaptTimestep = false;


    this->K = Kext;
    for(int i = 0; i < K + 1; i++){
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
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[mbx - 1] = TVector[mbx - 1] = mobod.getMobilizerTransform(someState);
    }

	// Initialize QsBuffer with zeros
	int nq = matter->getNQ(someState);
	int totSize = QsBufferSize * nq;
	for(int i = 0; i < totSize; i++){ 
		//QsBuffer.push_back(SimTK::Vector(nq, SimTK::Real(0)));
		QsBuffer.push_back(SimTK::Real(0));
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
        system->realize(someState, SimTK::Stage::Position);
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
SimTK::Real LAHMCSampler::MHAcceptProbability(SimTK::State&, SimTK::Real E_o, SimTK::Real E_n)
{
    // function args were SimTK::State& someState, SimTK::Real E_o, SimTK::Real E_n

	SimTK::Real Ediff = E_o - E_n;
	SimTK::Real prob = 1;
	if(Ediff < 0){
		prob = exp(this->beta * Ediff);
	}
	return prob;
}

/*** Convert C indeces to Ctau indeces ***/
void LAHMCSampler::C_to_Ctau_Indeces(int C_i, int C_j, int &Ctau_i, int &Ctau_j, int){
    // function args were int C_i, int C_j, int &Ctau_i, int &Ctau_j, int currSize

	Ctau_j = this->K - 1 - C_i;
	Ctau_i = this->K - 1 - C_j;
}

/*** Convert Ctau indeces to C indeces ***/
void LAHMCSampler::Ctau_to_C_Indeces(int Ctau_i, int Ctau_j, int &C_i, int &C_j, int currSize){
	int offset = (this->K+1) - currSize;
	C_j = offset + (currSize - 1) - Ctau_i;
	C_i = offset + (currSize - 1) - Ctau_j;
	
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
	std::size_t srcN = src.size();
	std::size_t destN = dest.size();

	// Only go up to the smallest dimension
	if(srcN < destN){ // Source is smaller
		for(std::size_t i = 0; i < srcN; i++){
			dest[i] = src[i];
		}
		
	}else{ // Destination is smaller
		for(std::size_t i = 0; i < destN; i++){
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
SimTK::Matrix& LAHMCSampler::leap_prob_recurse_hard(SimTK::State& someState, std::vector<SimTK::Real> Es, SimTK::Matrix& argCC)
{
	// Indices between the two matrices

	int M = argCC.nrow();
	int N = argCC.ncol();	
	SimTK_ASSERT_ALWAYS(M == N, "Robosample: leap_prob_recurse is designed for square matrices.");
	SimTK_ASSERT_ALWAYS(M == static_cast<int>(Es.size()), "Robosample: leap_prob_recurse: energy vector does not have the same size as C matrix."); // TODO 32 vs 64
	//std::cout << "\nBEGIN size " << M << "\n";

	// Check if already passed through this leaf
	SimTK::Real upperCorner = argCC.get(M - 1, N - 1 );

	//if( upperCorner != SimTK::Infinity){
	if( upperCorner < 1){
		//std::cout << "Already already visited this leaf" << std::endl;
		//PrintBigMat(argCC, M, N, 6, std::string(" C ") + std::to_string(M));
		return argCC;
	}

	// Boltzmann probability
	if(M == 2){
		//std::cout << "Do Metropolis Hastings.\n";
		SimTK::Real p_acc = 0.0;
    		if(!std::isnan(Es[M - 1])){
				p_acc = MHAcceptProbability(someState, Es[0], Es[Es.size() - 1]);
				argCC.set(0, M - 1, p_acc);
		}
		//PrintBigMat(argCC, M, N, 6, std::string(" C ") + std::to_string(M));
		return argCC;
	}
	
	//std::cout << "\nDEBUG Es_1 after injectFromTop \n";
	//for(int i=0; i<Es_1.size(); i++){std::cout << Es_1[i] << " ";}
	//std::cout << "\nDEBUG ===============\n";
	SimTK::Real cum_forward, cum_reverse;
	// Forward ===============================
	//std::cout << "Reenter forward";
	// Reduce size
	std::vector<SimTK::Real> Es_1(M - 1, SimTK::Infinity);
	injectFromTop(Es, Es_1);
	SimTK::Matrix CC_1 = extractFromTop(argCC, -1, -1);

	// Reentry
	CC_1 = leap_prob_recurse_hard(someState, Es_1, CC_1);

	// Recover data
	//injectFromTop(Es_1, Es); // don't need recover
	injectFromTop(CC_1, argCC);

	cum_forward = CC_1.get(0, CC_1.ncol() - 1);
	//std::cout << " cum_fwd = " << cum_forward << std::endl;

	// Backward =============================
	//std::cout << "Reenter reverse.";

	// Reverse
	SimTK::Matrix CC_rev(argCC.nrow(), argCC.ncol());
	CC_rev = reverseMatrix(argCC);
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

	SimTK::Matrix CC_rev_rev(argCC.nrow(), argCC.ncol());
	CC_rev_rev = reverseMatrix(CC_rev);

	//injectFromTop(Es_rev_rev, Es); // don't need recover
	injectFromTop(CC_rev_rev, argCC);

	//PrintBigMat(argCC, M, N, 6, std::string(" C ") + std::to_string(M));
	//std::cout << " cum_rev = " << cum_reverse << std::endl;

	// Eq. 25
	//std::cout << "Do LAHMC probability." << "\n";
	SimTK::Real Ediff;
	Ediff = Es[0] - Es[M - 1];
	SimTK::Real start_state_ratio = exp(this->beta * Ediff);

	SimTK::Real prob;
	prob = std::min(1 - cum_forward, start_state_ratio * (1 - cum_reverse));
	//std::cout << "E0 E1 start_state_ratio cum_fwd cum_rev prob " 
	//	<< Es[0] << " " << Es[Es.size() - 1] 
	//	<< " " << start_state_ratio << " " << cum_forward << " " << cum_reverse << " " 
	//	<< prob << std::endl;

	SimTK::Real cumu = cum_forward + prob;
	argCC.set(0, M - 1, cumu);

	//PrintBigMat(argCC, M, N, 6, std::string(" C ") + std::to_string(M));
	return argCC;
}

/** Acception rejection step **/
bool LAHMCSampler::accRejStep(SimTK::State& someState){
    //
    resetCMatrices(); // (K+1) dimension

    this->acc = false;

    SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    for(int kk = 0; kk < K; kk++){ // range(K)
	// Apply L
        integrateTrajectory(someState);

        calcNewConfigurationAndEnergies(someState, kk+1); // etot_ns[kk+1] set

	// Compute cumulative probability of doing this many leaps
	std::vector<SimTK::Real> subEs(&etot_ns[0], &etot_ns[kk+2]); // etot_ns[0 ... kk+1]
	//std::cout << "subEs[" << kk+1 << "] = " <<  etot_ns[kk+1] << std::endl;
	SimTK::Matrix CCk = extractFromTop(CC, kk+2, kk+2); // CC[:kk+2, :kk+2] 
	//PrintBigMat(CCk, kk+2, kk+2, 6, std::string(" CCk ") + std::to_string(kk+2));
	
	leap_prob_recurse_hard(someState, subEs, CCk);

	//std::cout << "rand_no " << rand_no << "; upper corner = " << CCk.get(0, kk+1) << std::endl;
	if(CCk.get(0, kk+1) >= rand_no){

		// Intermiediate steo - not necessary
		pe_n = pe_ns[kk+1];
    		fix_n = fix_ns[kk+1];
    		logSineSqrGamma2_n = logSineSqrGamma2_ns[kk+1];
    		ke_n = ke_ns[kk+1];
    		etot_n = etot_ns[kk+1];

		this->acc = true;
		std::cout << " acc" << std::endl;
		update(someState);

		injectFromTop(CCk, CC);
		break;
	}
    }

    // Anything left flip the momentum
    if(acc == false){
	std::cout << " nacc" << std::endl;
    	someState.updU() = -1.0 * someState.getU();
    }

    return this->acc;
}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
bool LAHMCSampler::propose(SimTK::State& someState)
{

    storeOldConfigurationAndPotentialEnergies(someState);

    initializeVelocities(someState);

    calcProposedKineticAndTotalEnergy(someState);

    // First set the starting point
    calcNewConfigurationAndEnergies(someState, 0);

    PrintDetailedEnergyInfo(someState);
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

/** Update. **/
void LAHMCSampler::update(SimTK::State& someState)
{
	setSetConfigurationAndEnergiesToNew(someState);
	++acceptedSteps;
	acceptedStepsBuffer.push_back(1);
	acceptedStepsBuffer.pop_front();
}

bool LAHMCSampler::sample_iteration(SimTK::State& someState)
{
	propose(someState);

	accRejStep(someState);
    
	++nofSamples;

	pushCoordinatesInR(someState);

	pushVelocitiesInRdot(someState);

	return this->acc;
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


