/**@file
Implementation of HMCSampler class. **/

#include "HMCSampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

//** Constructor **/
HMCSampler::HMCSampler(SimTK::CompoundSystem *argCompoundSystem
                                     ,SimTK::SimbodyMatterSubsystem *argMatter
                                     ,SimTK::Compound *argResidue
                                     ,SimTK::DuMMForceFieldSubsystem *argDumm
                                     ,SimTK::GeneralForceSubsystem *argForces
                                     ,SimTK::TimeStepper *argTimeStepper
                                     )
    : Sampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper),
    MonteCarloSampler(argCompoundSystem, argMatter, argResidue, argDumm, argForces, argTimeStepper)
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
    adaptTimestep = false;
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
        adaptTimestep = true;
    	timeStepper->updIntegrator().setFixedStepSize(-argTimestep);
        this->timestep = -argTimestep;
    }else{
        adaptTimestep = false;
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

/** It implements the proposal move in the Hamiltonian Monte Carlo
algorithm. It essentially propagates the trajectory after it stores
the configuration and energies. TODO: break in two functions:
initializeVelocities and propagate/integrate **/
void HMCSampler::propose(SimTK::State& someState)
{
    // Initialize configuration - not necessary unless we modify the
    // configuration in addition to velocities
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

    // Initialize velocities according to the Maxwell-Boltzmann distribution
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

    // Store the proposed energies
    // setProposedKE(matter->calcKineticEnergy(someState));
    this->ke_proposed = matter->calcKineticEnergy(someState);
    this->etot_proposed = getOldPE() + getProposedKE() + getOldFixman() + getOldLogSineSqrGamma2();


/*    // REC BUG
    std::cout << "HMC nbodies " << matter->getNumBodies() << std::endl;
    std::cout << "DuMM station_Bs before stepTo " << matter->getNumBodies() << std::endl;
    for (unsigned int i = 0; i < residue->getNumAtoms(); i++) {
        SimTK::Compound::AtomIndex aIx = (((Topology *)residue)->bAtomList[i]).getCompoundAtomIndex();
        SimTK::MobilizedBodyIndex mbx = residue->getAtomMobilizedBodyIndex(aIx);
        SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
        SimTK::Transform X_GB = mobod.getBodyTransform(someState);
        SimTK::DuMM::AtomIndex dAIx = residue->getDuMMAtomIndex(aIx);
        std::cout << "setAtomsLoc i aIx dAIx dumm.station_B gmol.locs " << i << " " << aIx
                  << " " << dAIx << " " << X_GB * dumm->getAtomStationOnBody(dAIx) << std::endl;
    } // REC BUG*/

    try {
        // IF NOT PRINT EVERY STEP // Integrate (propagate trajectory)
	if(adaptTimestep){
		std::cout << std::endl;
		//std::cout << "Adapt: nofSamples= " << nofSamples << std::endl;
		if( (nofSamples % acceptedStepsBufferSize) == (acceptedStepsBufferSize-1) ){ // Do it only so often
			std::cout << "Adapt BEGIN: ts= " << timestep;

			float stdAcceptance = 0.03;
			float idealAcceptance = 0.651;
			float smallestAcceptance = 0.0001;
			float timestepIncr = 0.001;

			// Compute acceptance in the buffer
			int sum = 0, sqSum = 0;
			for (std::list<int>::iterator p = acceptedStepsBuffer.begin(); p != acceptedStepsBuffer.end(); ++p){
        			sum += (int)*p;
				sqSum += (((int)*p) * ((int)*p));
			}
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
	
				//SimTK::Real F_n =   (a_n   - idealAcceptance) * (a_n   - idealAcceptance);
				//SimTK::Real F_n_1 = (a_n_1 - idealAcceptance) * (a_n_1 - idealAcceptance);
				//SimTK::Real F_n_2 = (a_n_2 - idealAcceptance) * (a_n_2 - idealAcceptance);
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

		}
        	this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));
	}else{
        	this->timeStepper->stepTo(someState.getTime() + (timestep*MDStepsPerSample));

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

	}
        system->realize(someState, SimTK::Stage::Position);

	/* // ELSE PRINT EVERY STEP
	for(int i = 0; i < MDStepsPerSample; i++){
		this->timeStepper->stepTo(someState.getTime() + (timestep));
                system->realize(someState, SimTK::Stage::Position);

                // / * // INSTANT GEOMETRY
                SimTK::Vec3 a1pos, a2pos, a3pos, a4pos, a5pos;
                int a1, a2, a3, a4, a5;
                a1 = 16; a2 = 14; a3 = 0; a4 = 6; a5 = 8;
                a1pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a1)));
                a2pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a2)));
                a3pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a3)));
                a4pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a4)));
                a5pos = ((Topology *)residue)->calcAtomLocationInGroundFrame(someState, SimTK::Compound::AtomIndex(SimTK::Compound::AtomIndex(a5)));
                int distA1, distA2, distA3, distA4;
                SimTK::Vec3 distA1pos, distA2pos;
                SimTK::Vec3 distA3pos, distA4pos;
                distA1 = 2; distA2 = 17; distA3 = 6; distA4 = 17;
                std::cout << "geom "  << bDihedral(a1pos, a2pos, a3pos, a4pos) ;
                std::cout << " "  << bDihedral(a2pos, a3pos, a4pos, a5pos) ;
                std::cout << std::endl;
                // * / // INSTANT GEOMETRY END

		/ *    if(useFixman){
		        fix_n = calcFixman(someState); logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
		    }else{
		        fix_n = 0.0; logSineSqrGamma2_n = 0.0;
		    }
		    // Get new kinetic energy
		    system->realize(someState, SimTK::Stage::Velocity); ke_n = matter->calcKineticEnergy(someState);
		    // Get new potential energy
		    if ( getThermostat() == ANDERSEN ){
		        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
			    //detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
			    logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
		    }
		    else{
		        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
		        //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL
			//detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
		        logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
		    }	
		    // Calculate total energy
		    if(useFixman){
		        etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
		        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
		    }else{
		        etot_n = pe_n + ke_n;
		        etot_proposed = pe_o + ke_proposed;
		    }
    		PrintDetailedEnergyInfo(someState); // * /
	} // */ 
	// END PRINT EVERY STEP

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

/** Main function that contains all the 3 steps of HMC.
Implements the acception-rejection step and sets the state of the
compound to the appropriate conformation wether it accepted or not. **/
void HMCSampler::update(SimTK::State& someState)
{
    SimTK::Real rand_no = uniformRealDistribution(randomEngine);

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
    //std::cout << " ke after timestepping " << this->ke_n << std::endl;

    // Get new potential energy
    if ( getThermostat() == ANDERSEN ){
        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
        //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL
	    //detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
	    logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }
    else{
        pe_n = forces->getMultibodySystem().calcPotentialEnergy(someState);
        //pe_n = dumm->CalcFullPotEnergyIncludingRigidBodies(someState); // ELIZA FULL
	    //detmbat_n = ((Topology *)residue)->calcLogDetMBAT(someState);
        logSineSqrGamma2_n = ((Topology *)residue)->calcLogSineSqrGamma2(someState);
    }

    // Calculate total energy
    if(useFixman){
        etot_n = pe_n + ke_n + fix_n - (0.5 * RT * logSineSqrGamma2_n);
        etot_proposed = pe_o + ke_proposed + fix_o - (0.5 * RT * logSineSqrGamma2_o);
    }else{
        etot_n = pe_n + ke_n;
        etot_proposed = pe_o + ke_proposed;
    }

    PrintDetailedEnergyInfo(someState);


    // Decide and get a new sample
    if ( getThermostat() == ANDERSEN ){ // MD with Andersen thermostat
        this->acc = true;
        std::cout << " acc" << std::endl;
        setSetTVector(someState);
        pe_set = pe_n;
        fix_set = fix_n;
        logSineSqrGamma2_set = logSineSqrGamma2_n;
        ke_lastAccepted = ke_n;
        etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;


        ++acceptedSteps;
        acceptedStepsBuffer.push_back(1);
        acceptedStepsBuffer.pop_front();
    }else{ // Apply Metropolis correction
        if ( (proposeExceptionCaught == false) &&
                (!std::isnan(pe_n)) && ((etot_n < etot_proposed) ||
             (rand_no < exp(-(etot_n - etot_proposed) * this->beta)))) { // Accept based on full energy
             //(rand_no < exp(-(pe_n - pe_o) * this->beta)))) { // Accept based on potential energy
             this->acc = true;
             std::cout << " acc" << std::endl;
             setSetTVector(someState);
             pe_set = pe_n;
             fix_set = fix_n;
             logSineSqrGamma2_set = logSineSqrGamma2_n;
             ke_lastAccepted = ke_n;
             etot_set = pe_set + fix_set + ke_proposed + logSineSqrGamma2_set;


             ++acceptedSteps;
             acceptedStepsBuffer.push_back(1);
             acceptedStepsBuffer.pop_front();

        } else { // Reject
             this->acc = false;
             std::cout << " nacc" << std::endl;
             assignConfFromSetTVector(someState);
             proposeExceptionCaught = false;
             acceptedStepsBuffer.push_back(0);
             acceptedStepsBuffer.pop_front();
        }
    }


	// MSD
	unsigned int natoms = ((Topology *)residue)->bAtomList.size();
	SimTK::Vec3 R;
	SimTK::Vec3 Rdot;
	SimTK::Compound::AtomIndex aIx; 
	SimTK::Compound c = *((Topology *)residue);

	// Load new coordinates
	for(unsigned int j = 0; j < natoms; j++){
		aIx = (((Topology *)residue)->bAtomList[j]).getCompoundAtomIndex();

		R = c.calcAtomLocationInGroundFrame(someState, aIx);
		Rs.push_back(R);

		Rdot = c.calcAtomVelocityInGroundFrame(someState, aIx);
		Rdots.push_back(Rdot);
	}
	//for(unsigned int j = 0; j < Rdots.size(); j++){std::cout << Rdots[j] << " " << (Rdots[j]).normalize() << std::endl;}

	// Calculate MSD
	SimTK::Real MSD = 0;
	SimTK::Real RRdot = 0;
	SimTK::Vec3 dR;

	if(Rs.size() == (2*natoms)){ // Don't delete initial coordinates
		for(unsigned int j = 0; j < natoms; j++){
			//std::cout << Rs[j + natoms] - Rs[j] << std::endl;
			//std::cout << "size= " << Rs.size() << " a= " << j + natoms << " b= " << j << std::endl;
			dR = Rs[j + natoms] - Rs[j];
			MSD += dR.normSqr();
		}
	
		MSD /= natoms;
		std::cout << "MSD= " << MSD << std::endl;

		// Compute RRdot
		for(unsigned int j = 0; j < natoms; j++){
			dR = Rs[j + natoms] - Rs[j];
			if(((Rdots[j]).norm() != 0) && (dR.norm() != 0)){
				//std::cout << "dR= " << dR << " " << dR.normalize() << std::endl;
				//std::cout << "Rdot" << Rdots[j] << " " << (Rdots[j]).normalize() << std::endl;
				RRdot += SimTK::dot(dR.normalize(), (Rdots[j]).normalize());
			}
		}
		std::cout << "RRdot= " << RRdot << std::endl;

		// Transfer upper to lower half and cleanup the upper half
		for(int j = (natoms - 1); j >= 0; --j){
			//std::cout << "size= " << Rs.size() << " a= " << j << " b= " << j  + natoms << std::endl;
			Rs[j] = Rs[j + natoms];
			Rs.pop_back();
		}
	}


	// Cleanup Rdots
	for(unsigned int j = 0; j < natoms; j++){
		Rdots.pop_back();
	}

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

    //++nofSamples;

}

bool HMCSampler::sample_iteration(SimTK::State& someState)
{
	propose(someState);
	update(someState);
	++nofSamples;
	return this->acc;
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



