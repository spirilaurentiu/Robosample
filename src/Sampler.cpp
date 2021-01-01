/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(SimTK::CompoundSystem *argCompoundSystem,
                 SimTK::SimbodyMatterSubsystem *argMatter,
                 //Topology *argResidue,

                 //SimTK::Compound *argResidue,
		 std::vector<Topology *>& argTopologies,

                 SimTK::DuMMForceFieldSubsystem *argDumm,
                 SimTK::GeneralForceSubsystem *argForces,
                 SimTK::TimeStepper *argTimeStepper)
{
	this->compoundSystem = argCompoundSystem;
	this->matter = argMatter;

	//this->residue = argResidue;
	this->topologies = argTopologies;
	assert(topologies.size());
	this->residue = topologies[0];

	// Set total number of atoms and dofs
	natoms = 0;
	for ( unsigned int i = 0; i < topologies.size(); i++){
		natoms += (topologies[i])->getNumAtoms();
	}
	int ThreeFrom3D = 3;
	ndofs = natoms * ThreeFrom3D;

	this->dumm = argDumm;
	this->forces = argForces;
	this->timeStepper = argTimeStepper;
	this->system = &(matter->getSystem());
}

// Destructor
Sampler::~Sampler(){
    ;
}

// Get set the seed
unsigned long long int Sampler::getSeed(void)
{
    return this->seed;
}

/** Store the value of the seed internally and also feed it to the random
number generator **/
void Sampler::setSeed(unsigned long long int argSeed)
{
    if(argSeed == 0){
        std::chrono::system_clock::time_point tp
            = std::chrono::system_clock::now();
        std::chrono::system_clock::duration dtn = tp.time_since_epoch();
        long clockSeed = dtn.count();
        randomEngine.seed( clockSeed );
        seed = clockSeed;
    }else{
        randomEngine.seed( argSeed );
        seed = argSeed;
    }

}

// Compute mass matrix determinant
SimTK::Real Sampler::calcMassDeterminant(SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

// Compute mass matrix determinant

SimTK::Real Sampler::calcMassDeterminant(const SimTK::State& state)
{
    int nu = state.getNU();
    SimTK::Vector V(nu);
    SimTK::Vector DetV(nu);
    SimTK::Real D0 = 1.0;
    matter->calcDetM(state, V, DetV, &D0);
    return D0;
}

/** Returns the number of MC trials done by this integrator. **/
int Sampler::getNofSamples(void){
    return nofSamples;
}


// Update - to be implemented by every specific sampler

//void Sampler::update(SimTK::State& somState){}

// TODO move
void Sampler::PrintSimbodyStateCache(SimTK::State& someState){
    std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
    for(int i = 0; i < someState.getNumSubsystems(); i++){
        std::cout
            << " Subsystem Name: "
            << someState.getSubsystemName(SimTK::SubsystemIndex(i))
            << " Stage: "
            << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
            << " Version: "
            << someState.getSubsystemVersion(SimTK::SubsystemIndex(i))
            << std::endl;
    }
}

void Sampler::initialize(SimTK::State& someState) {
    // Sampling
    nofSamples = 0;
}

void Sampler::reinitialize(SimTK::State& someState) {
    // Sampling
    nofSamples = 0;
}

// Getter for macroscopic temperature
SimTK::Real Sampler::getTemperature() const {
    return temperature;
}

/** Setter for macroscopic temperature. Also sets the RT and beta**/
void Sampler::setTemperature(SimTK::Real temperature) {
    Sampler::temperature = temperature;
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
    beta = 1.0 / RT;
}

/** Setter for macroscopic beta. Also sets the RT and temperature**/
void Sampler::setBeta(SimTK::Real argBeta) {
    this->beta = argBeta;
    Sampler::temperature = 1 / (beta * SimTK_BOLTZMANN_CONSTANT_MD);
    RT = this->temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

SimTK::Real Sampler::getRT() const {
    return RT;
}

SimTK::Real Sampler::getBeta() const {
    return beta;
}

SimTK::Real Sampler::generateRandomNumber(GmolRandDistributionType distributionType) {
    if(distributionType == GmolRandDistributionType::UNIFORM){
        return uniformRealDistribution(randomEngine);
    }else{
        return gaurand(randomEngine);
    }
}

/** Load the map of mobods to joint types **/
void Sampler::loadMbx2mobility(SimTK::State& someState)
{

	// Lop through topologies
	for(int i = 0; i < topologies.size(); i++){

		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
	
			if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin

				// Get body, parentBody
				SimTK::MobilizedBodyIndex mbx = (topologies[i])->getAtomMobilizedBodyIndex(aIx);
				const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
				const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
				SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
	
				if(parentMbx != 0){
					// Get the neighbor atom in the parent mobilized body
					SimTK::Compound::AtomIndex chemParentAIx = (topologies[i])->getChemicalParent(matter, aIx);
				
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
				
					//std::cout << "mbx= " << mbx << " parentMbx= " << parentMbx
					//	<< " aIx= " << aIx << " chemParentAIx= " << chemParentAIx
					//	<< " mobility " << mobility
					//	<< std::endl;
	
				} // Parent is not Ground

			} // if is a root atom
		
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
                //std::cout << "mbx= " << mbx << " mobility= " << mbx2mobility[mbx]
                //      << std::endl
                //;

		JointType jointType;
		SimTK::BondMobility::Mobility mobility = mbx2mobility[mbx];
		int qi, internQIx;
		switch(mobility) {
			case SimTK::BondMobility::Mobility::Free: ///< Unrestricted bond, permitting changes in stretch, bend, and torsion modes
                		//internQIx = -1;
                		//for(qi = qIx; qi < (mobodNQ + qIx); qi++){internQIx++;}
				jointType = QUATERNION_a;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_b;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_c;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_d;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Torsion: ///< Bond has fixed length and angles, but permits rotation about the bond axis
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Rigid: ///< Bond links both atoms to the same rigid unit
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::BallF: ///< Three rotational dofs. It allows  angle flexibility besides torsion.
				jointType = QUATERNION_a;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_b;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_c;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_d;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::BallM: ///< Three rotational dofs. It allows  angle flexibility besides torsion.
				jointType = QUATERNION_a;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_b;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_c;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = QUATERNION_d;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Cylinder: ///< Torsion plus translation along the bond
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Translation: ///< Three translational mobilities (Cartesian). // NEWMOB
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::FreeLine: ///< Three translational mobilities (Cartesian). // NEWMOB
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::LineOrientationF: ///< Two rotational mobilities // NEWMOB
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::LineOrientationM: ///< Two rotational mobilities // NEWMOB
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::UniversalM: ///< Cap de bara
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Spherical: ///< BAT coordinates
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = ANGULAR180;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::AnglePin:
				jointType = ANGULAR360;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			case SimTK::BondMobility::Mobility::Slider: ///< Translation along bond
				jointType = LINEAR;
				qIndex2jointType.insert( std::pair<SimTK::QIndex, JointType> (SimTK::QIndex(qi), jointType) );
				break;
			default:
				std::cout << "Warning: unknown joint type" << std::endl;
		}


                // Loop through body's Us
                for(int ui = 0; ui < mobodNU; ui++){
                        //SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(ui));
                        //std::cout << "H_FMCol= " << H_FMCol << std::endl;
                }

        }

}

