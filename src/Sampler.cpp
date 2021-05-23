/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(World *argWorld,
	SimTK::CompoundSystem *argCompoundSystem,
	SimTK::SimbodyMatterSubsystem *argMatter,
	//SimTK::Compound *argResidue,
	std::vector<Topology> &argTopologies,
	SimTK::DuMMForceFieldSubsystem *argDumm,
	SimTK::GeneralForceSubsystem *argForces,
	SimTK::TimeStepper *argTimeStepper) :
		world(argWorld),
		compoundSystem(argCompoundSystem),
		matter(argMatter),
		topologies(argTopologies),
		dumm(argDumm),
		forces(argForces),
		timeStepper(argTimeStepper),
		thermostat(ThermostatName::NONE),
		temperature(SimTK::Zero),
		RT(SimTK::Zero),
		beta(std::numeric_limits<SimTK::Real>::min()),
		seed(std::numeric_limits<uint32_t>::min()),
		nofSamples(std::numeric_limits<int>::min()),
		acc(false)
{
	assert(argCompoundSystem != nullptr);
	assert(argMatter != nullptr);
	assert(argDumm != nullptr);
	assert(argForces != nullptr);
	assert(argTimeStepper != nullptr);

	this->system = &argMatter->getSystem();

	//this->rootTopology = argResidue;
	assert(topologies.size() > 0);
	this->rootTopology = &topologies[0];

	// Set total number of atoms and dofs
	natoms = 0;
	for (const auto& topology: topologies){
		natoms += topology.getNumAtoms();
	}

	int ThreeFrom3D = 3;
	ndofs = natoms * ThreeFrom3D;
}

// Destructor
Sampler::~Sampler() {
}

// Get set the seed
uint32_t Sampler::getSeed() const
{
    return this->seed;
}

/** Store the value of the seed internally and also feed it to the random
number generator **/
void Sampler::setSeed(uint32_t argSeed)
{
    if(argSeed == 0){
		// TODO seed is 32 bit, but chrono returns 64 bit
        std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
        std::chrono::system_clock::duration dtn = tp.time_since_epoch();
        uint32_t clockSeed = static_cast<uint32_t>(dtn.count());
        randomEngine.seed(clockSeed);
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
int Sampler::getNofSamples(){
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

void Sampler::initialize(SimTK::State&) {
    // Sampling
	// function args were SimTK::State& someState
    nofSamples = 0;
}

void Sampler::reinitialize(SimTK::State&) {
    // Sampling
	// function args were SimTK::State& someState
    nofSamples = 0;
}

// Getter for macroscopic temperature
SimTK::Real Sampler::getTemperature() const {
	assert(!SimTK::isNumericallyEqual(temperature, SimTK::Zero));
    return temperature;
}

/** Setter for macroscopic temperature. Also sets the RT and beta**/
void Sampler::setTemperature(SimTK::Real inTemperature) {
    this->temperature = inTemperature;
    RT = this->temperature * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
    beta = 1.0 / RT;
}

/** Setter for macroscopic beta. Also sets the RT and temperature**/
void Sampler::setBeta(SimTK::Real argBeta) {
    this->beta = argBeta;
    this->temperature = 1 / (beta * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD));
    RT = this->temperature * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD);
}

SimTK::Real Sampler::getRT() const {
	assert(!SimTK::isNumericallyEqual(RT, SimTK::Zero));
    return RT;
}

SimTK::Real Sampler::getBeta() const {
	assert(!SimTK::isNumericallyEqual(beta, SimTK::Zero));
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
void Sampler::loadMbx2mobility(SimTK::State&)
{
	// function args were SimTK::State& someState

	// Lop through topologies
	for(auto& topology : topologies){

		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topology.getNumAtoms(); ++aIx){
	
			if(topology.getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin

				// Get body, parentBody
				const SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);
				const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
				const SimTK::MobilizedBody& parentMobod = mobod.getParentMobilizedBody();
				SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
	
				if(parentMbx != 0){
					// Get the neighbor atom in the parent mobilized body
					SimTK::Compound::AtomIndex chemParentAIx = topology.getChemicalParent(matter, aIx);
				
					// Get mobility (joint type)
					const auto& bond = topology.getBond(topology.getNumber(aIx), topology.getNumber(chemParentAIx));
					auto mobility = bond.getBondMobility();
	
					mbx2mobility.insert(std::make_pair(mbx, mobility));
				
					//std::cout << "mbx= " << mbx << " parentMbx= " << parentMbx
					//	<< " aIx= " << aIx << " chemParentAIx= " << chemParentAIx
					//	<< " mobility " << mobility
					//	<< std::endl;
	
				} // Parent is not Ground

			} // if is a root atom
		
		} // END loop through atoms
	
	} // END loop through topologies

    for (SimTK::MobilizedBodyIndex mbx(2); mbx < matter->getNumBodies(); ++mbx){
        // const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
        // SimTK::QIndex qIx = mobod.getFirstQIndex(someState);
        // int mobodNQ = mobod.getNumQ(someState);
        // int mobodNU = mobod.getNumU(someState);

        //const SimTK::Transform T = mobod.getMobilizerTransform(someState);
        //const SimTK::MassProperties mp = mobod.getBodyMassProperties(someState);
        //const SimTK::UnitInertia unitInertia = mobod.getBodyUnitInertiaAboutBodyOrigin(someState);
        //std::cout << "mbx= " << mbx << " mobility= " << mbx2mobility[mbx]
        //      << std::endl
        //;

		SimTK::BondMobility::Mobility mobility = mbx2mobility[mbx];
		int qi = 1; // TODO what value?
		switch(mobility) {
			///< Unrestricted bond, permitting changes in stretch, bend, and torsion modes
			case SimTK::BondMobility::Mobility::Free:
				//int internQIx = -1;
				//for(qi = qIx; qi < (mobodNQ + qIx); qi++){internQIx++;}

				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_a));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_b));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_c));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_d));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Bond has fixed length and angles, but permits rotation about the bond axis
			case SimTK::BondMobility::Mobility::Torsion:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				break;

			///< Bond links both atoms to the same rigid unit
			case SimTK::BondMobility::Mobility::Rigid:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Three rotational dofs. It allows  angle flexibility besides torsion.
			case SimTK::BondMobility::Mobility::BallF:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_a));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_b));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_c));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_d));
				break;

			///< Three rotational dofs. It allows  angle flexibility besides torsion.
			case SimTK::BondMobility::Mobility::BallM:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_a));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_b));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_c));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::QUATERNION_d));
				break;

			///< Torsion plus translation along the bond
			case SimTK::BondMobility::Mobility::Cylinder:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Three translational mobilities (Cartesian). // NEWMOB
			case SimTK::BondMobility::Mobility::Translation:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Three translational mobilities (Cartesian). // NEWMOB
			case SimTK::BondMobility::Mobility::FreeLine:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Two rotational mobilities // NEWMOB
			case SimTK::BondMobility::Mobility::LineOrientationF:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Two rotational mobilities // NEWMOB
			case SimTK::BondMobility::Mobility::LineOrientationM:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			///< Cap de bara
			case SimTK::BondMobility::Mobility::UniversalM:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				break;

			///< BAT coordinates
			case SimTK::BondMobility::Mobility::Spherical:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR180));
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			case SimTK::BondMobility::Mobility::AnglePin:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::ANGULAR360));
				break;

			///< Translation along bond
			case SimTK::BondMobility::Mobility::Slider:
				qIndex2jointType.insert(std::make_pair(SimTK::QIndex(qi), JointType::LINEAR));
				break;

			default:
				std::cout << "Warning: unknown joint type" << std::endl;
		}

		// // Loop through body's Us
		// for(int ui = 0; ui < mobodNU; ui++){
		// 	SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(ui));
		// 	std::cout << "H_FMCol= " << H_FMCol << std::endl;
		// }
    }

}

