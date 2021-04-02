/**@file
Implementation of HMCSampler class. **/

#include "HMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"
#include "World.hpp"

//** Constructor **/
HMCSampler::HMCSampler(World *argWorld, SimTK::CompoundSystem *argCompoundSystem
                                     ,SimTK::SimbodyMatterSubsystem *argMatter

                                     //,SimTK::Compound *argResidue
				     ,std::vector<Topology *>& argTopologies

                                     ,SimTK::DuMMForceFieldSubsystem *argDumm
                                     ,SimTK::GeneralForceSubsystem *argForces
                                     ,SimTK::TimeStepper *argTimeStepper
                                     )
    //: Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    : Sampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper),
    //MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
    MonteCarloSampler(argWorld, argCompoundSystem, argMatter, argTopologies, argDumm, argForces, argTimeStepper)
{
	this->useFixman = false;  
	this->fix_n = this->fix_o = 0.0;
	this->logSineSqrGamma2_n = this->logSineSqrGamma2_o = 0.0;
	this->residualEmbeddedPotential = 0.0;
	nofSamples = 0;
	this->alwaysAccept = false;

	this->prevPrevAcceptance = SimTK::NaN;
	this->prevAcceptance = SimTK::NaN;
	this->acceptance = SimTK::NaN;

	
	this->prevPrevTimestep = SimTK::NaN; // ps
	this->prevTimestep = SimTK::NaN; // ps
	this->timestep = 0;

	this->temperature = 0.0;
	this->boostT = this->temperature;
	MDStepsPerSample = 0;
	proposeExceptionCaught = false;
	shouldAdaptTimestep = false;

	dR.resize(ndofs, 0);
	dRdot.resize(ndofs, 0);
}

/** Destructor **/
HMCSampler::~HMCSampler()
{
}

/** Calculate sqrt(M) using Eigen. For debug purposes. **/
void HMCSampler::calcNumSqrtMUpper(SimTK::State& someState, SimTK::Matrix& SqrtMUpper)
{
    assert("!Not implemented");
}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l = sqrt(D) * [I -JPsiK]. This is upper triangular matrix and it is computed 
multipling a set of orthonormal vectors with the sqrt(MInv). **/
void HMCSampler::calcSqrtMInvU(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvU: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = SqrtMInvV[j];
        }
        V[i] = 0;
    }


}

/** Calculate O(n2) the square root of the mass matrix inverse
denoted by Jain l* = [I -JPsiK]*sqrt(D) (adjoint of l).
This is lower triangular matrix and it is computed by multipling a set of
 orthonormal vectors with the sqrt(MInv) and transpose it. **/
void HMCSampler::calcSqrtMInvL(SimTK::State& someState, SimTK::Matrix& SqrtMInv)
{
    int nu = someState.getNU();
    assert((SqrtMInv.nrow() == nu) && (SqrtMInv.ncol() == nu) && "calcSqrtMInvL: passed matrix doesn't have nu x nu size.");

    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);

    // Zero the matrix and the temporary vector
    for (int i=0; i < nu; ++i){
        V[i] = 0;
        for (int j=0; j < nu; ++j){
            SqrtMInv(i, j) = 0;
        }
    }

    // Calculate the values inside the matrix
    for (int i=0; i < nu; ++i){
        V[i] = 1;
        matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);
        for (int j=0; j < nu; ++j){
            SqrtMInv(j, i) = SqrtMInvV[j];
        }
        V[i] = 0;
    }

}

/** Stores the accepted kinetic energy. This should be set right after a
move is accepted. It's a component of the total energy stored. **/
void HMCSampler::setLastAcceptedKE(SimTK::Real inpKE)
{
    this->ke_lastAccepted = inpKE;
}

/** Sets the proposed kinetic energy before the proposal. This should be
set right after the velocities are initialized. **/
void HMCSampler::setProposedKE(SimTK::Real inpKE)
{
    this->ke_proposed = inpKE;
}

/** Get/set the TimeStepper that manages the integrator **/
const SimTK::TimeStepper * HMCSampler::getTimeStepper(void)
{
    return timeStepper;
}

SimTK::TimeStepper * HMCSampler::updTimeStepper(void)
{
    return timeStepper;
}

void HMCSampler::setTimeStepper(SimTK::TimeStepper * someTimeStepper)
{
    timeStepper = someTimeStepper;
}

/** Seed the random number generator. Set simulation temperature,
velocities to desired temperature, variables that store the configuration
and variables that store the energies, both needed for the
acception-rejection step. Also realize velocities and initialize
the timestepper. **/
void HMCSampler::initialize(SimTK::State& someState )
{
    // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    timeStepper->initialize(compoundSystem->getDefaultState());

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
    if(useFixman){
        std::cout << "Hamiltonian Monte Carlo sampler: using Fixman potential." << std::endl;
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
//r void HMCSampler::reinitialize(SimTK::State& someState, SimTK::Real timestep, int nosteps, SimTK::Real argTemperature) 
void HMCSampler::reinitialize(SimTK::State& someState)
{
     // After an event handler has made a discontinuous change to the
    // Integrator's "advanced state", this method must be called to 
    // reinitialize the Integrator.
    //(this->timeStepper->updIntegrator()).reinitialize(SimTK::Stage::Topology, false);

    // Store the configuration
	// In our case, a MobilizedBody is an atom
    system->realize(someState, SimTK::Stage::Position);
    int i = 0;
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        SetTVector[i] = TVector[i] = mobod.getMobilizerTransform(someState);
        i++;
    }

    // Store potential energies
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

/** Get/Set the timestep for integration **/
float HMCSampler::getTimestep(void)
{
    return timestep;
}

void HMCSampler::setTimestep(float argTimestep)
{
    if(argTimestep <= 0){
        shouldAdaptTimestep = true;
    	timeStepper->updIntegrator().setFixedStepSize(-argTimestep);
        this->timestep = -argTimestep;
    }else{
        shouldAdaptTimestep = false;
    	timeStepper->updIntegrator().setFixedStepSize(argTimestep);
    	this->timestep = argTimestep;
    }
}

/** Get/Set boost temperature **/
SimTK::Real HMCSampler::getBoostTemperature(void)
{
    return this->boostT;
}

void HMCSampler::setBoostTemperature(SimTK::Real argT)
{
    this->boostT = argT;
    this->boostFactor = std::sqrt(this->boostT / this->temperature);
    this->unboostFactor = 1 / boostFactor;
    std::cout << "HMC: boost temperature: " << this->boostT << std::endl;
    std::cout << "HMC: boost velocity scale factor: " << this->boostFactor << std::endl;
}

void HMCSampler::setBoostMDSteps(int argMDSteps)
{
    this->boostMDSteps = argMDSteps;
    std::cout << "HMC: boost MD steps: " << this->boostMDSteps << std::endl;

}

/** Store configuration as a set of Transforms **/
void HMCSampler::storeOldConfigurationAndPotentialEnergies(SimTK::State& someState){ // func
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
void HMCSampler::initializeVelocities(SimTK::State& someState){
    int nu = someState.getNU();
    double sqrtRT = std::sqrt(RT);
    SimTK::Vector V(nu);
    SimTK::Vector SqrtMInvV(nu);
    for (int i=0; i < nu; ++i){
        V[i] = gaurand(randomEngine);
    }

    // RED ZONE
    ////////
    for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
        const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        const int mnu = mobod.getNumU(someState);
	const float scaleFactor = world->getMobodUScaleFactor(mbx);
        std::cout << "RED ZONE mbx scaleFactor uIxes " << int(mbx) << ' ' << scaleFactor;
	for(SimTK::UIndex uIx = mobod.getFirstUIndex(someState); uIx < mobod.getFirstUIndex(someState) + mobod.getNumU(someState); uIx++ ){
        	std::cout << ' ' << int(uIx) ;
		V[int(uIx)] *= scaleFactor;
		
        }
        std::cout << '\n';
    }

/*
    SimTK::Vector One(nu, 1.0);
    float dotProduct = 0.0;
    for (int i=0; i < nu; ++i){
        dotProduct += (One[i] * V[i]);
    }
    float phi = std::acos(dotProduct  / (V.norm() * One.norm()) );

    SimTK::Vector UScaleFactors(nu, 0.0);
    for (int i=0; i < nu; ++i){
        UScaleFactors[i] = i ;
    }

    dotProduct = 0.0;
    for (int i=0; i < nu; ++i){
        dotProduct += (UScaleFactors[i] * V[i]);
    }
    float psi = std::acos(dotProduct  / (V.norm() * UScaleFactors.norm()) );
    std::cout << "RED ZONE " << phi << ' ' << psi << '\n';
*/
    // RED ZONE END


    matter->multiplyBySqrtMInv(someState, V, SqrtMInvV);

    SqrtMInvV *= (sqrtRT); // Set stddev according to temperature

    // Raise the temperature
    someState.updU() = SqrtMInvV;

    // Realize velocity
    system->realize(someState, SimTK::Stage::Velocity);
}

/** Store the proposed energies **/
void HMCSampler::calcProposedKineticAndTotalEnergy(SimTK::State& someState){

    // Store proposed kinetic energy
    // setProposedKE(matter->calcKineticEnergy(someState));
    this->ke_proposed = matter->calcKineticEnergy(someState);

    // Store proposed total energy
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();
}

void HMCSampler::adaptTimestep(SimTK::State& someState)
{
		std::cout << std::endl;
		//std::cout << "Adapt: nofSamples= " << nofSamples << std::endl;
		if( (nofSamples % acceptedStepsBufferSize) == (acceptedStepsBufferSize-1) ){ // Do it only so often
			std::cout << "Adapt BEGIN: ts= " << timestep;

			// float stdAcceptance = 0.03;
			float idealAcceptance = 0.651;
			// float smallestAcceptance = 0.0001;
			// float timestepIncr = 0.001;

			// Compute acceptance in the buffer
			int sum = std::accumulate(acceptedStepsBuffer.begin(), acceptedStepsBuffer.end(), 0);

			SimTK::Real prevPrevPrevAcceptance = prevPrevAcceptance; // to reset
			prevPrevAcceptance = prevAcceptance;
			prevAcceptance = acceptance;
    		acceptance = float(sum) / float(acceptedStepsBufferSize);

			std::cout << " ppAcc= " << prevPrevAcceptance << " pAcc= " << prevAcceptance << " acc= " << acceptance << std::endl;

			if( !SimTK::isNaN(prevPrevAcceptance) ){ // Passed first two initial evaluations
				// Calculate gradients
				SimTK::Real a_n, a_n_1, a_n_2, t_n, t_n_1, t_n_2;
				a_n = acceptance; a_n_1 = prevAcceptance; a_n_2 = prevPrevAcceptance;
				t_n = timestep; t_n_1 = prevTimestep; t_n_2 = prevPrevTimestep;
	
				SimTK::Real F_n =   std::abs(a_n   - idealAcceptance);
				SimTK::Real F_n_1 = std::abs(a_n_1 - idealAcceptance);
				SimTK::Real F_n_2 = std::abs(a_n_2 - idealAcceptance);
	
				SimTK::Real dF_n   = F_n     - F_n_1;
				SimTK::Real dF_n_1 = F_n_1   - F_n_2;
				SimTK::Real dt_n = (t_n   - t_n_1);
				SimTK::Real dt_n_1 = (t_n_1   - t_n_2);
				SimTK::Real gradF_n   = dF_n / dt_n;
				SimTK::Real gradF_n_1 = dF_n_1 / dt_n_1;
	
				// Calculate gamma
				SimTK::Real gamma;
				SimTK::Real num, denom;
				num = std::abs(dt_n * (gradF_n - gradF_n_1));
				denom = (gradF_n - gradF_n_1) * (gradF_n - gradF_n_1);
				gamma = num / denom;
	
				// Increase or decrease timestep
				SimTK::Real newTimestep = t_n - (gamma * gradF_n);
				std::cout << "Adapt data: "
					<< " " << F_n  << " " << F_n_1  << " " << F_n_2 << " " << dF_n << " " << dF_n_1 
					<< " " << dt_n << " " << dt_n_1 << " " << gradF_n << " " << gradF_n_1
					<< " " << num << " " << denom << " " << gamma
					<< std::endl;
				if( (std::abs(prevTimestep - timestep) > 0.0001) || (newTimestep < 0.00089) || (SimTK::isNaN(newTimestep)) || ((acceptance - prevAcceptance) < 0.001)){
					std::cout << "newTimestep rejected" << std::endl; 
    					SimTK::Real r = uniformRealDistribution(randomEngine);
					if(acceptance < idealAcceptance){
						timestep = timestep - (0.3*r * timestep);
					}else{
						timestep = timestep + (0.3*r * timestep);
					}

					acceptance = prevAcceptance;
					prevAcceptance = prevPrevAcceptance;
					prevPrevAcceptance = prevPrevPrevAcceptance;
				}else{ // Reset acceptances
					std::cout << "newTimestep accepted" << std::endl; 
					prevPrevTimestep = prevTimestep;
					prevTimestep = timestep;
					timestep = newTimestep;
				}

    				timeStepper->updIntegrator().setFixedStepSize(timestep);
				std::cout << "Adapt END: "  << " ppTs= " << prevPrevTimestep << " pTs= " << prevTimestep << " ts= " << timestep << std::endl;
			}else{ // Alter the intial timesteps to get a valid dF_n next time

			//SimTK::ArticulatedInertia abi = matter->getArticulatedBodyInertia(someState, SimTK::MobilizedBodyIndex(2));
		    	//const SimTK::MobilizedBody& mobod2 = matter->getMobilizedBody(SimTK::MobilizedBodyIndex(2));
			//const SimTK::MassProperties mp2 = mobod2.getBodyMassProperties(someState);
			//const SimTK::Inertia i2 = mp2.calcInertia();
			//SimTK::Mat33 mi2 = i2.toMat33();
			//std::cout << std::endl;
			//std::cout << mi2(0,0) << " " << mi2(0,1) << " " << mi2(0,2) << std::endl;
			//std::cout << mi2(1,0) << " " << mi2(1,1) << " " << mi2(1,2) << std::endl;
			//std::cout << mi2(2,0) << " " << mi2(2,1) << " " << mi2(2,2) << std::endl;
			//std::cout << "det(2)^1/5 " << std::pow(SimTK::det(mi2), 0.2) << std::endl;

				prevPrevTimestep = prevTimestep;
				prevTimestep = timestep;
				if( SimTK::isNaN(prevPrevTimestep) ){
					timestep = 0.0009; // smaller than H vibration
    					timeStepper->updIntegrator().setFixedStepSize(timestep);
				}else{
					timestep = (prevTimestep + prevPrevTimestep) / 2;
    					timeStepper->updIntegrator().setFixedStepSize(timestep);
				}
				std::cout << "Adapt END: ppTs= " << prevPrevTimestep << " pTs= " << prevTimestep << " ts= " << timestep << std::endl;
			}

		} // is time to adapt 
}

/**  **/
void HMCSampler::adaptWorldBlocks(SimTK::State& someState){

    std::cout << std::setprecision(10) << std::fixed;
	int nq = someState.getNQ();
	// int totSize = QsBufferSize * nq;
	if( (nofSamples % QsBufferSize) == (QsBufferSize-1) ){
		std::cout << "Adapt blocks BEGIN: " << std::endl;
		//std::cout << "Print by row: " << std::endl;
		//for(int confi = 0; confi < QsBufferSize; confi++){
		//	for(int qi = 0; qi < nq; qi++){
		//		SimTK::Real x;
		//		int pos;
		//
		//		std::list<SimTK::Real>::iterator it = QsBuffer.begin();
		//		pos = (confi * nq) + qi;
		//		std::advance(it, pos);
		//		std::cout << *it << ' ';
		//	}
		//	std::cout << std::endl;
		//}

		//std::cout << "Print by column: " << std::endl;
		std::vector<std::vector<SimTK::Real>> QsBufferVec(nq, vector<SimTK::Real>(QsBufferSize));
		for(int qi = 0; qi < nq; qi++){

			std::vector<SimTK::Real> thisQ(QsBufferSize);

			//std::cout << "QsBuffer= " << std::endl;
			for(int confi = 0; confi < QsBufferSize; confi++){
				thisQ[confi] = QsBuffer[(confi * nq) + qi];
				//std::cout << QsBuffer[(confi * nq) + qi] << ' ';
			}
			//std::cout << std::endl;

			//std::cout << "thisQ= " << std::endl;
			//for(int confi = 0; confi < QsBufferSize; confi++){
			//	std::cout << thisQ[confi] << ' ';
			//}
			//std::cout << std::endl;

			QsBufferVec[qi] = thisQ;

		}

		//std::cout << "QsBufferVec= " << std::endl;
		//for(int qi = 0; qi < nq; qi++){
		//	for(int confi = 0; confi < QsBufferSize; confi++){
		//		std::cout << QsBufferVec[qi][confi] << ' ';
		//	}
		//	std::cout << std::endl;
		//}
		
		std::cout << "Compute correlation matrix: " << std::endl;
		SimTK::Real corr;
		for(int qi = 0; qi < nq; qi++){
			//for(int qj = 0; (qj < nq) && (qj < qi); qj++){
			for(int qj = 0; (qj < nq); qj++){
				if(qIndex2jointType[SimTK::QIndex(qi)] == ANGULAR360){
					if(qIndex2jointType[SimTK::QIndex(qj)] == ANGULAR360){
						corr = circCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else if(qIndex2jointType[SimTK::QIndex(qj)] == LINEAR){
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else{
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}
				}else if(qIndex2jointType[SimTK::QIndex(qi)] == LINEAR){
					if(qIndex2jointType[SimTK::QIndex(qj)] == LINEAR){
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else if(qIndex2jointType[SimTK::QIndex(qj)] == ANGULAR360){
						corr = circCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}else{
						corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
						std::cout << corr << ' ';
					}
				}else{
					corr = bCorr(QsBufferVec[qi], QsBufferVec[qj]);
					std::cout << corr << ' ';
				}
			}
			std::cout << std::endl;
		}

		std::cout << "Adapt blocks END: " << std::endl;
	}

}


/** Apply the L operator **/
void HMCSampler::integrateTrajectory(SimTK::State& someState){
    try {

        this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));

        system->realize(someState, SimTK::Stage::Position);

    }catch(const std::exception&){
        proposeExceptionCaught = true;
        int i = 0;
        for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
            const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
            mobod.setQToFitTransform(someState, SetTVector[i]);
            i++;
        }
        system->realize(someState, SimTK::Stage::Position);
    }

}

/**  **/
void HMCSampler::integrateTrajectoryOneStepAtATime(SimTK::State& someState
	//, std::function<void(void)> initialFunc
	//, std::function<void(void)> loopFunc
	//, std::function<void(void)> finalFunc
	){

	//// TODO: Intial function BEGIN ////
	//// Push initial R and Rdots
	//SimTK::Vec3 atomR;
	//SimTK::Vec3 atomRdot;
	//SimTK::Compound::AtomIndex aIx; 
	//for(int i = 0; i < topologies.size(); i++){
	//	for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
	//		aIx = ((topologies[i])->bAtomList[j]).getCompoundAtomIndex();
	//		atomR = (topologies[i])->calcAtomLocationInGroundFrame(someState, aIx);
	//		R.push_back(atomR[0]);
	//		R.push_back(atomR[1]);
	//		R.push_back(atomR[2]);
	//	}
	//}
	//// Initial function END ////

	//// Main loop ////
	for(int k = 0; k < MDStepsPerSample; k++){
		this->timeStepper->stepTo(someState.getTime() + (timestep));
                system->realize(someState, SimTK::Stage::Position);

		//// Calculate generalized quantities
    		//int nu = someState.getNU();
    		//SimTK::Vector MU(nu);
    		//for (int i=0; i < nu; ++i){
        	//	MU[i] = 1;
		//}
		//matter->multiplyByM(someState, someState.getU(), MU);
		//std::cout << "MU= " << MU << std::endl;
		
    		//int nu = someState.getNU();
    		//SimTK::Vector U(nu);
		//U = someState.getU();
		//std::cout << "U= " << U << std::endl;

		/////////////////////////


		//// TODO:Loop function BEGIN ////
		//	// Push current R and Rdot
		//	for(int i = 0; i < topologies.size(); i++){
		//		for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
		//			aIx = ((topologies[i])->bAtomList[j]).getCompoundAtomIndex();
		//			atomR = (topologies[i])->calcAtomLocationInGroundFrame(someState, aIx);
		//			R.push_back(atomR[0]);
		//			R.push_back(atomR[1]);
		//			R.push_back(atomR[2]);
		//			atomRdot = (topologies[i])->calcAtomVelocityInGroundFrame(someState, aIx);
		//			Rdot.push_back(atomRdot[0]);
		//			Rdot.push_back(atomRdot[1]);
		//			Rdot.push_back(atomRdot[2]);
		//		}
		//	}
		//
		//	// Calculate dR		
		//	for(unsigned int j = 0; j < (ndofs); j++){
		//		dR[j] = R[j + (ndofs)] - R[j];
		//	}
		//
		//	// Calculate MSD
		//	SimTK::Real MSD = 0;
		//	if(dR.size() >= (ndofs)){
		//		MSD += magSq(dR);
		//		MSD /= (ndofs);
		//	}
		//
		//	// Calculate RRdot	
		//	SimTK::Real RRdot = 0;
		//	if(dR.size() >= (ndofs)){
		//		std::vector<SimTK::Real> tempDR = dR;
		//		std::vector<SimTK::Real> tempRdot = Rdot;
		//
		//		normalize(tempDR);
		//		normalize(tempRdot);
		//
		//		for(unsigned int j = 0; j < ndofs; j++){
		//			RRdot += tempDR[j] * tempRdot[j];
		//		}
		//	}
		//	std::cout << std::setprecision(10) << std::fixed;
		//	std::cout << "k= " << k << " RRdot= " << RRdot
		//		<< " MSD= " << MSD << std::endl;
		//
		//	// Pop current R and Rdot
		//	for(int j = ((ndofs) - 1); j >= 0; --j){
		//		R.pop_back();
		//		Rdot.pop_back();
		//	}
		//// Loop function END ////

	} // END main loop

	//// TODO: Final function BEGIN ////
	//// Pop initial R
	//for(int j = ((ndofs) - 1); j >= 0; --j){
	//	R.pop_back();
	//}
	//// Final function END ////

}

void HMCSampler::geomDihedral(void){
    /* // INSTANT GEOMETRY
    SimTK::Vec3 a1pos, a2pos, a3pos, a4pos, a5pos;
    int a1, a2, a3, a4, a5;
    //DIHEDRAL 4 6 8 14 6 8 14 16
    //a1 = 16; a2 = 14; a3 = 0; a4 = 6; a5 = 8;
    a1 = 4; a2 = 6; a3 = 8; a4 = 14; a5 = 16;
    a1pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
    a2pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
    a3pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
    a4pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
    a5pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a5)));
    int distA1, distA2, distA3, distA4;
    SimTK::Vec3 distA1pos, distA2pos;
    SimTK::Vec3 distA3pos, distA4pos;
    distA1 = 2; distA2 = 17; distA3 = 6; distA4 = 17;
//		    distA1pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA1)));
//		    distA2pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA2)));
//		    distA3pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA3)));
//		    distA4pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(distA4)));
//            std::cout << " dihedral elements"
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a1))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a2))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a3))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a4))
//              << " " << residue->getAtomElement(SimTK::Compound::AtomIndex(a5))
//              << std::endl;
//            std::cout << " poss: " << a1pos << ' ' << a2pos << ' ' << a3pos << ' ' << a4pos << ' ';
    std::cout << "geom "  << bDihedral(a1pos, a2pos, a3pos, a4pos) ;
    std::cout << " "  << bDihedral(a2pos, a3pos, a4pos, a5pos) ;
    //std::cout << " "  << 10 * (distA1pos - distA2pos).norm() ; // Ang
    //std::cout << " "  << 10 * (distA3pos - distA4pos).norm() ; // Ang
    std::cout << std::endl;
    // */
    // INSTANT GEOMETRY END
}

/** Store new configuration and energy terms**/
void HMCSampler::calcNewConfigurationAndEnergies(SimTK::State& someState)
{
    // Get new Fixman potential
    if(useFixman){
        fix_n = calcFixman(someState);
        logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }else{
        fix_n = 0.0;
        logSineSqrGamma2_n = 0.0;
    }

    // Get new kinetic energy
    system->realize(someState, SimTK::Stage::Velocity);
    ke_n = matter->calcKineticEnergy(someState);

    // Get new potential energy
    pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
    // TODO: replace with the following after checking is the same thing
    //pe_n = compoundSystem->calcPotentialEnergy(someState);

    logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);

    // Calculate total energy
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }
}

/** Set energies and configuration to new state **/
void HMCSampler::setSetConfigurationAndEnergiesToNew(SimTK::State& someState)
{
    setSetTVector(someState);
    pe_set = pe_n;
    fix_set = fix_n;
    logSineSqrGamma2_set = logSineSqrGamma2_n;
    ke_lastAccepted = ke_n;
    etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;
}

void HMCSampler::setSetConfigurationAndEnergiesToOld(SimTK::State& someState)
{
	assignConfFromSetTVector(someState);
	proposeExceptionCaught = false;
	acceptedStepsBuffer.push_back(0);
	acceptedStepsBuffer.pop_front();
}

SimTK::Real HMCSampler::MHAcceptProbability(SimTK::State& someState, 
	SimTK::Real argEtot_proposed, SimTK::Real argEtot_n){

	if(argEtot_n < argEtot_proposed){
		return 1;
	}else{
		return exp(-(argEtot_n - argEtot_proposed) * this->beta);
	}
}

/** Acception rejection step **/
bool HMCSampler::accRejStep(SimTK::State& someState){

	// Decide and get a new sample
	if ( getThermostat() == ThermostatName::ANDERSEN ){ // MD with Andersen thermostat
		this->acc = true;
		std::cout << "\tsample accepted (always with andersen thermostat)" << std::endl;
		// update(someState); // smth is broken in here
	}else{ // Apply Metropolis-Hastings correction
		if (proposeExceptionCaught == false && !std::isnan(pe_n)) { 
			const SimTK::Real rand_no = uniformRealDistribution(randomEngine);
			if(rand_no < MHAcceptProbability(someState, etot_proposed, etot_n)){ 
				this->acc = true;
				std::cout << "\tsample accepted" << std::endl;
				update(someState);
			}else{ // Reject
				this->acc = false;
				std::cout << "\tsample rejected" << std::endl;
				setSetConfigurationAndEnergiesToOld(someState);
			}
		}
	}

	return this->acc;
}

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. **/
void HMCSampler::propose(SimTK::State& someState)
{

	// Store old configuration
	storeOldConfigurationAndPotentialEnergies(someState);

	// Initialize velocities according to the Maxwell-Boltzmann distribution
	initializeVelocities(someState);

	// Store the proposed energies 
	calcProposedKineticAndTotalEnergy(someState);

	// Adapt timestep
	if(shouldAdaptTimestep){
		adaptTimestep(someState);
	}
	
	// Adapt timestep
	bool shouldAdaptWorldBlocks = false;
	if(shouldAdaptWorldBlocks){
		adaptWorldBlocks(someState);
	}
	
	// Apply the L operator 
	integrateTrajectory(someState);
	//integrateTrajectoryOneStepAtATime(someState);

	calcNewConfigurationAndEnergies(someState);

	PrintDetailedEnergyInfo(someState);

}

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HMCSampler::update(SimTK::State& someState)
{
	setSetConfigurationAndEnergiesToNew(someState); 
    ++acceptedSteps;
    acceptedStepsBuffer.push_back(1);
    acceptedStepsBuffer.pop_front();
}

/** Push Cartesian coordinates into R vector stored in Sampler.
Return the size of R **/
int HMCSampler::pushCoordinatesInR(SimTK::State& someState)
{
	for(const auto& Topology : topologies){
		for(const auto& AtomList : Topology->bAtomList){
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto& atomR = Topology->calcAtomLocationInGroundFrame(someState, aIx);
			R.insert(R.end(), { atomR[0], atomR[1], atomR[2] });
		}
	}

	if(R.size() >= static_cast<std::size_t>(2 * ndofs)) {
		// for some weird reason, g++9.2 always evaluates the first if to true (45 >= 90)
		if(R.size() >= static_cast<std::size_t>(2 * ndofs)) {
			return -1;
		}

		for(unsigned int j = 0; j < ndofs; j++){
			dR[j] = R[j + ndofs] - R[j];
		}

		// // Transfer upper to lower half and cleanup the upper half of R
		// for(int j = ((ndofs) - 1); j >= 0; --j){
		// 	//std::cout << "R size= " << R.size() << " j= " << j << " j+ndofs= " << j+ndofs << std::endl;
		// 	R[j] = R[j + ndofs];
		// 	R.pop_back();
		// }

		// If stuff breaks, look up for the original code. This also compacts memory
		// See https://stackoverflow.com/questions/7351899/remove-first-n-elements-from-a-stdvector for more details
		std::vector<decltype(R)::value_type>(R.begin() + ndofs, R.end()).swap(R);
	}

}

/** Push velocities into Rdot vector stored in Sampler.
Return the size of Rdot **/
int HMCSampler::pushVelocitiesInRdot(SimTK::State& someState)
{
	for(const auto& Topology : topologies){
		for(const auto& AtomList : Topology->bAtomList){
			const auto aIx = AtomList.getCompoundAtomIndex();
			const auto& atomRdot = Topology->calcAtomVelocityInGroundFrame(someState, aIx);
			Rdot.insert(Rdot.end(), { atomRdot[0], atomRdot[1], atomRdot[2] });
		}
	}

	// Calculate dRdot		
	if(Rdot.size() >= static_cast<size_t>(2 * ndofs)) {
		// for some weird reason, g++9.2 always evaluates the first if to true (45 >= 90)
		if(Rdot.size() >= static_cast<size_t>(2 * ndofs)) {
			return -1;
		}

		for(int j = 0; j < ndofs; j++){
			dRdot[j] = Rdot[j + ndofs] - Rdot[j];
		}

		// // Transfer upper to lower half and cleanup the upper half of Rdot
		// for(int j = ((ndofs) - 1); j >= 0; --j){
		// 	std::cout << "Rdotdot size= " << Rdot.size() << " j= " << j << " j+ndofs= " << j+ndofs << std::endl;
		// 	Rdot[j] = Rdot[j + ndofs];
		// 	Rdot.pop_back();
		// }

		// If stuff breaks, look up for the original code. This also compacts memory
		// See https://stackoverflow.com/questions/7351899/remove-first-n-elements-from-a-stdvector for more details
		std::vector<decltype(Rdot)::value_type>(Rdot.begin() + ndofs, Rdot.end()).swap(Rdot);
	}
}

bool HMCSampler::sample_iteration(SimTK::State& someState)
{
    std::cout << std::setprecision(10) << std::fixed;

	propose(someState);

	accRejStep(someState);

	++nofSamples;

	// Add generalized coordinates to a buffer
	auto Q = someState.getQ(); // g++17 complains if this is auto& or const auto&
	QsBuffer.insert(QsBuffer.end(), Q.begin(), Q.end());
	QsBuffer.erase(QsBuffer.begin(), QsBuffer.begin() + Q.size());

	pushCoordinatesInR(someState);
	pushVelocitiesInRdot(someState);

	// Calculate MSD and RRdot to adapt the integration length
    std::cout << std::setprecision(10) << std::fixed;
	std::cout << "\tMSD= " << calculateMSD() << ", RRdot= " << calculateRRdot() << std::endl;
	
	return this->acc;
}


/** Calculate Mean Square Displacement based on stored R vectors **/
SimTK::Real HMCSampler::calculateMSD(void)
{
	SimTK::Real MSD = 0;
	if(dR.size() >= static_cast<size_t>(ndofs)){
		MSD += magSq(dR);
		MSD /= (ndofs);
	}
	return MSD;
}

/** Calculate RRdot based on stored R and Rdot vectors **/
SimTK::Real HMCSampler::calculateRRdot(void)
{
	SimTK::Real RRdot = 0;
	if(dR.size() >= static_cast<size_t>(ndofs)){
		std::vector<SimTK::Real> tempDR = dR;
		std::vector<SimTK::Real> tempRdot = Rdot;
		normalize(tempDR);
		normalize(tempRdot);

		for(int j = 0; j < ndofs; j++){
			RRdot += tempDR[j] * tempRdot[j];
		}
	}
	return RRdot;
}



int HMCSampler::getMDStepsPerSample() const {
    return MDStepsPerSample;
}

void HMCSampler::setMDStepsPerSample(int mdStepsPerSample) {
    MDStepsPerSample = mdStepsPerSample;
}

/** Print detailed energy information **/
void HMCSampler::PrintDetailedEnergyInfo(SimTK::State& someState)
{
    std::cout << std::setprecision(5) << std::fixed;
    std::cout
		<< "\tpe_o " << pe_o << ", pe_n " << pe_n << ", pe_nB " << getPEFromEvaluator(someState)
        << "\n\tke_prop " << ke_proposed << ", ke_n " << ke_n
        << "\n\tfix_o " << fix_o << ", fix_n " << fix_n
        << "\n\tlogSineSqrGamma2_o " << logSineSqrGamma2_o << ", logSineSqrGamma2_n " << logSineSqrGamma2_n
        //<< " detmbat_n " << detmbat_n //<< " detmbat_o " << detmbat_o << " "
        << "\n\tts " << timestep  << ", exp(bdE) " << exp(-(etot_n - etot_proposed) / RT)
        << "\n\tetot_n " << etot_n  << ", etot_proposed " << etot_proposed
        << std::endl;
}


/** Modifies Q randomly
 **/
void HMCSampler::perturbQ(SimTK::State& someState)
{
    // Perturb Q
    //SimTK::Real rand_no = uniformRealDistribution(randomEngine);
    //SimTK::Real rand_no = uniformRealDistribution_mpi_pi(randomEngine);
    int nq = someState.getNQ();
    //SimTK::Vector V(nq);
    for (int i=7; i < nq; ++i){
        //V[i] = uniformRealDistribution_mpi_pi(randomEngine);
        someState.updQ()[i] = uniformRealDistribution_mpi_pi(randomEngine);
    }
    std::cout << "perturbQ " << someState.getQ() << std::endl;
    system->realize(someState, SimTK::Stage::Position);

    // Get needed energies
    pe_o  = getOldPE();
    if(useFixman){
        fix_o = getOldFixman();
    }
    if(useFixman){
        fix_n = calcFixman(someState);
    }else{
        fix_n = 0.0;
    }

    pe_n = getPEFromEvaluator(someState); // OPENMM
    //std::cout << "Multibody PE " << getPEFromEvaluator(someState) << std::endl; // OPENMM
    //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL

    setSetTVector(someState);
    setSetPE(pe_n);
    setSetFixman(fix_n);
    ++acceptedSteps;
    assignConfFromSetTVector(someState);

    std::cout << someState.getNU() << ' ' << 1 << ' '
              //<< getSetPE() + getREP() << ' ' << getLastAcceptedKE()
              << getSetPE() << ' ' << 0
              << ' ' << getSetFixman() << ' ' << fix_o << ' ' << fix_n << ' ';

    // Keep track of how many MC trials have been done
    ++nofSamples;
}


/** Load the map of mobods to joint types **/
/*
void HMCSampler::loadMbx2mobility(SimTK::State& someState)
{
	// Lop through topologies
	for(int i = 0; i < topologies.size(); i++){
		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
	
			// Get body, parentBody
			SimTK::MobilizedBodyIndex mbx = (topologies[i])->getAtomMobilizedBodyIndex(aIx);
			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
			SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
			if(parentMbx != 0){
				// Get the neighbor atom in the parent mobilized body
				SimTK::Compound::AtomIndex chemParentAIx = (topologies[i])->getChemicalParent(matter, aIx);
			
				//std::cout << "mbx= " << mbx << " parentMbx= " << parentMbx
				//	<< " aIx= " << aIx << " chemParentAIx= " << chemParentAIx
				//	<< std::endl;
				// Get Top to parent frame
				SimTK::Compound::AtomIndex parentRootAIx = (topologies[i])->mbx2aIx[parentMbx];
			
				// Get mobility (joint type)
				bSpecificAtom *atom = (topologies[i])->updAtomByAtomIx(aIx);
				SimTK::BondMobility::Mobility mobility;
				bBond bond = (topologies[i])->getBond(
					(topologies[i])->getNumber(aIx), (topologies[i])->getNumber(chemParentAIx));
				mobility = bond.getBondMobility();
				mbx2mobility.insert(std::pair<SimTK::MobilizedBodyIndex, SimTK::BondMobility::Mobility>
                        		(mbx, mobility));
			
				(mobility == SimTK::BondMobility::Mobility::Torsion);
			
			} // Parent is not Ground
		
		} // END loop through atoms
	
	} // END loop through topologies
        for (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
                const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
                SimTK::QIndex qIx = mobod.getFirstQIndex(someState);
                int mobodNQ = mobod.getNumQ(someState);
                int mobodNU = mobod.getNumU(someState);
                //const SimTK::Transform T = mobod.getMobilizerTransform(someState);
                //const SimTK::MassProperties mp = mobod.getBodyMassProperties(someState);
                //const SimTK::UnitInertia unitInertia = mobod.getBodyUnitInertiaAboutBodyOrigin(someState);
                std::cout << "mbx= " << mbx << " mobility= " 
                //      << std::endl
                ;
                // Loop through body's Qs
                int internQIx = -1;
                for(int qi = qIx; qi < (mobodNQ + qIx); qi++){
                        internQIx++;
                        std::cout << " " << qi ;
                }
                std::cout << std::endl;
                // Loop through body's Us
                for(int ui = 0; ui < mobodNU; ui++){
                        //SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(ui));
                        //std::cout << "H_FMCol= " << H_FMCol << std::endl;
                }
        }
}
*/
