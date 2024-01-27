/**@file
Implementation of Sampler class. **/

#include "Sampler.hpp"

// Includes to get the structure of additional classes

#include "Topology.hpp"

// Constructor

Sampler::Sampler(World &argWorld,
		SimTK::CompoundSystem &argCompoundSystem,
		SimTK::SimbodyMatterSubsystem &argMatter,
		std::vector<Topology> &argTopologies, 
		SimTK::DuMMForceFieldSubsystem &argDumm,
		SimTK::GeneralForceSubsystem &argForces,
		SimTK::TimeStepper &argTimeStepper) :
		world(&argWorld),
		compoundSystem(&argCompoundSystem),
		matter(&argMatter),
		topologies(argTopologies),
		dumm(&argDumm),
		forces(&argForces),
		timeStepper(&argTimeStepper),
		thermostat(ThermostatName::NONE),
		temperature(SimTK::Zero),
		RT(SimTK::Zero),
		beta(std::numeric_limits<SimTK::Real>::min()),
		seed(std::numeric_limits<uint32_t>::min()),
		nofSamples(std::numeric_limits<int>::min()),
		acc(false),
		system(&argMatter.getSystem()),
		alwaysAccept(false)
{
	//this->rootTopology = argResidue;
	assert(topologies.size() > 0);
	rootTopology = &topologies[0];

	// Set total number of atoms and dofs
	natoms = 0;
	for (const auto& topology: topologies){
		natoms += topology.getNumAtoms();
	}

	// Set total mass of the system to non-realistic valaue
	this->totalMass = 0;

	int ThreeFrom3D = 3;
	ndofs = natoms * ThreeFrom3D;

	int test_dofs = 3;
	gammarand = std::gamma_distribution<double>((test_dofs - 1.0) / 2.0, 1.0);

}

// Destructor
Sampler::~Sampler() {
}

// Random numbers

// Get set the seed
uint32_t Sampler::getSeed() const
{
    return this->seed;
}

/** Store the value of the seed internally and also feed it to the random
number generator **/
void Sampler::setSeed(uint32_t argSeed)
{
    if(argSeed == 0) {
		// We use chrono to get the time in ns.
        std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
        std::chrono::system_clock::duration dtn = tp.time_since_epoch();
        int64_t clockSeed = dtn.count();

		// Save the seed.
        this->seed = static_cast<RANDOM_ENGINE_INIT::result_type>(clockSeed);
    } else {
		// The seed is already provided, we just need to save it
        this->seed = argSeed;
    }

	// We need many bytes of random data to initialize a Mersenne Twister.
	// To ease reproducibility, we use one 32-bit seed to initialize a less powerful RNG.
	// Then, we use this RNG to initialize the state of MT.
	randomEngineInit.seed(seed);

	// Initial generator function.
	auto source = [this]() {
		// This LCG generates 64 bits.
		// According to this https://stackoverflow.com/a/18262495, the lower bits are of not-so-exceptional quality.
		// This is also stated in https://en.wikipedia.org/wiki/Linear_congruential_generator.
		// In order to mitigate this, only keep the highest 32 bits.
		// return static_cast<RANDOM_ENGINE_INIT::result_type>(randomEngineInit() >> sizeof(RANDOM_ENGINE_INIT::result_type)); // result is already 32 bit

		return static_cast<RANDOM_ENGINE_INIT::result_type>(randomEngineInit());
	};

	// Size for initial state bits
	constexpr std::size_t N = RANDOM_ENGINE::state_size * sizeof(RANDOM_ENGINE::result_type);
	constexpr auto Size = (N - 1) / sizeof(source()) + 1; // TODO is this correct?

	// Generate random data for seeding
	std::array<RANDOM_ENGINE_INIT::result_type, Size> RandomData;
	std::generate(RandomData.begin(), RandomData.end(), std::ref(source));

	// Seed the engine
	std::seed_seq SeedSeq(RandomData.begin(), RandomData.end());
	randomEngine.seed(SeedSeq);

	// According to Matsumoto [1][2], we should discard the first Size (or Size*2) numbers.
	// This isn't the case here, but we do it just in case.
	// [1] https://doi.org/10.1145/1276927.1276928
	// [2] https://stats.stackexchange.com/questions/436733/is-there-such-a-thing-as-a-good-bad-seed-in-pseudo-random-number-generation
	randomEngine.discard(Size * 2);
}

// Is the sampler always accepting the proposed moves
bool Sampler::getAlwaysAccept(void) const
{
    return alwaysAccept;
}

// Is the sampler always accepting the proposed moves
void Sampler::setAlwaysAccept(bool argAlwaysAccept)
{
    alwaysAccept = argAlwaysAccept;
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

	int test_dofs = 3;
	gammarand = std::gamma_distribution<double>((test_dofs - 1.0) / 2.0, 1.0);
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

SimTK::Real
Sampler::generateRandomNumber(GmolRandDistributionType distributionType)
{
    if(distributionType == GmolRandDistributionType::UNIFORM){
        return uniformRealDistribution(randomEngine);
    }else{
        return gaurand(randomEngine);
    }
}

SimTK::Quaternion Sampler::generateRandomQuaternion(void)
{
		double x,y,z, u,v,w, s;
		do {
			x = uniformRealDistribution_m1_1(randomEngine);
			y = uniformRealDistribution_m1_1(randomEngine);
			z = x*x + y*y;
		}while (z > 1);
		do {
			u = uniformRealDistribution_m1_1(randomEngine);
			v = uniformRealDistribution_m1_1(randomEngine);
			w = u*u + v*v;
		} while (w > 1);

		s = sqrt((1-z) / w);

		SimTK::Quaternion randQuat(x, y, s*u, s*v);

		return randQuat;
}


// Just for checking
void Sampler::checkAtomStationsThroughDumm(void)
{
	// Lop through topologies
	for(auto& topology : topologies){

		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topology.getNumAtoms(); ++aIx){

			//std::cout << "DEBUG getAtomLocationInMobilizedBodyFrame"
			//	<< " aIx " << aIx;

			//SimTK::Vec3 atomMobodStation =
			//	topology.getAtomLocationInMobilizedBodyFrame(aIx);
			SimTK::Vec3 atomMobodStationThroughDumm =
				topology.getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, *dumm);

			//std::cout << " through Compound " << atomMobodStation
			//	<< std::endl;
			//std::cout << " through DuMM " << atomMobodStationThroughDumm
			//	<< std::endl;
		}
	}

}

/** Load the map of mobods to joint types **/
// void Sampler::loadMbx2mobility(SimTK::State&) // SAFE
// The ideal way to do this would be to keep have MobilizedBody keep
// a Mobility name
void Sampler::loadMbx2mobility(int whichWorld) // DANGER
{

	// Lop through topologies
	for(auto& topology : topologies){

		// Loop through atoms
		for (SimTK::Compound::AtomIndex aIx(0); aIx < topology.getNumAtoms();
		++aIx)
		{

			// Get atom location in mobod
			SimTK::Vec3 atomLoc =
				topology.getAtomLocationInMobilizedBodyFrameThroughDumm(
					aIx, *dumm);

			// Mobod origin atom
			if(atomLoc == 0){

				// Get body, parentBody
				const SimTK::MobilizedBodyIndex mbx = 
					topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *dumm);
				const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
				const SimTK::MobilizedBody& parentMobod = mobod.getParentMobilizedBody();
				SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

				// Parent is not Ground
				if(parentMbx != 0){

					// Get the neighbor atom in the parent mobilized body
					SimTK::Compound::AtomIndex chemParentAIx =
						topology.getChemicalParent_IfIAmRoot(matter, aIx, *dumm);

					// Get mobility (joint type)
					const auto& bond = topology.getBond(
						topology.getNumber(aIx),
						topology.getNumber(chemParentAIx));

					auto mobility = bond.getBondMobility(whichWorld);

					// Insert mobility into the map
					mbx2mobility.insert(std::make_pair(mbx, mobility));

					//std::cout << "mbx= " << mbx << " parentMbx= " << parentMbx
					//	<< " aIx= " << aIx << " chemParentAIx= " << chemParentAIx
					//	<< " mobility " << mobility
					//	<< std::endl;

				} // parent is not Ground

			} // if is a root atom

		} // every atom

	} // every topology

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

			///< Torsion and translation
			case SimTK::BondMobility::Mobility::BendStretch:
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

		// Loop through body's Us
		// for(int ui = 0; ui < mobodNU; ui++){
		// 	SimTK::SpatialVec H_FMCol = mobod.getH_FMCol(someState, SimTK::MobilizerUIndex(ui));
		// 	std::cout << "H_FMCol= " << H_FMCol << std::endl;
		// }

    }

}


// Draws from von Mises distribution with concentration kK
SimTK::Real Sampler::generateVonMisesSample(SimTK::Real Mean, SimTK::Real K){

	// Declarations
	SimTK::Real tau = 0.0;
	SimTK::Real rho = 0.0;
	SimTK::Real r = 0.0;
	SimTK::Real z = 0.0;
	SimTK::Real f = 0.0;
	SimTK::Real c = 0.0;
	SimTK::Real step2Val = 0.0;
	SimTK::Real step3Val = 0.0;
	SimTK::Real step4Val = 0.0;

	SimTK::Real u1 = uniformRealDistribution(randomEngine);
	SimTK::Real u2 = uniformRealDistribution(randomEngine);
	SimTK::Real u3 = uniformRealDistribution(randomEngine);

	SimTK::Real theta = SimTK::NaN;

	// Prereq
	tau = 1 + std::sqrt(1 + 4*K*K);
	rho = (tau - std::sqrt(2*tau)) / (2*K);

	// Step 0 variables
	r = (1 + (rho*rho)) / (2 * rho);

	// Step 1 variables
	z = std::cos(SimTK_PI * u1);
	f = (1 + (r*z)) / (r + z);
	c = K * (r - f);

	// Step 2 variables
	step2Val = (c * (2 - c)) - u2;

	// Step 3 variables
	step3Val = std::log(c / u2) + 1 - c;

	step4Val = ((u3 - 0.5) < 0) ? -1 : 1;
	step4Val = step4Val * std::acos(f);

	for(int i = 0; i < 1000; i++){
		if(step2Val > 0.0){
			theta = step4Val;
			break;
		}else{
			if(step3Val < 0){
				;
			}else{
				theta = step4Val;
				break;
			}
		}
	}

	// Translate to the mean
	if( isnan(theta) ){
		std::cout << "Sampler::vonMises warning did not draw after 1000 tries.\n";
		return theta;
	}else{
		theta += Mean;
		if(theta > SimTK_PI){
			theta = theta - (2*SimTK_PI);
		}else if(theta < (-1.0 * SimTK_PI)){
			theta = (2*SimTK_PI) + theta;
		}

	}

	return theta;

}


// Draws from von Mises-Fisher distribution
std::vector<double>& Sampler::generateVonMisesFisherSample(std::vector<double>& X,
	double lambda)
{
	//int NDOFS = 3; // DELETE THIS
	std::cout << "Sampler::vonMisesFisher ndofs " << ndofs << "\n";

	double ndofs_1 = ndofs - 1;

	std::vector<double> V;
	V.resize(ndofs_1, 0.0);
	std::vector<double> X_lowdim;
	X_lowdim.resize(ndofs_1, 0.0);

	// STEP 0
	double b = -2.0 * lambda;
	b += std::sqrt((4 * (lambda*lambda)) + (ndofs_1*ndofs_1));
	b /= ndofs_1;

	double x0 = (1 - b) / (1 + b);

	double LOG = ndofs_1 * std::log(1 - (x0*x0));

	double c = (lambda*x0) + LOG;

	int nofTries = 10;
	for(int t = 0; t < nofTries; t++){
		// STEP 1
		double Z1 = gammarand(randomEngine);
		double Z2 = gammarand(randomEngine);
		double Z = Z1 / (Z1 + Z2);

		double U = uniformRealDistribution(randomEngine);

		double W_num = 1.0 - ((1.0 + b)*Z);
		double W_den = 1.0 - ((1.0 - b)*Z);
		double W = W_num / W_den;

		// STEP 2
		if( ((lambda*W) + LOG - c) >= std::log(U) ){
			// Generate uniform random vector on sphere ndofs - 1
			for(int j = 0; j < ndofs_1; j++){
				V[j] = gaurand(randomEngine);
			}
			bNormalizeInPlace(V);

			// Actual sample
			bMulByScalar(V, std::sqrt(1 - (W*W)), X_lowdim);
			for(int j = 1; j < ndofs; j++){
				X[j] = X_lowdim[j-1];
			}
			X[0] = W;

			break;
		}
	}

	return X;
}


// Draws from von Mises-Fisher distribution
double Sampler::generateChiSample(void)
{
	// TODO
	assert(!"Not implemented!");

	int NDOFS = 3;

	return 0;
}


// TODO revise param1 and param2
/**
 * Note that parameters here are not the parameters of the distribution
*/
SimTK::Real Sampler::convoluteVariable(SimTK::Real& passedVariable,
		std::string distrib,
		SimTK::Real param1, SimTK::Real param2, SimTK::Real param3)
{
	// TODO this string compare is slow
	// Bernoulli trial between the var and its inverse
	if(distrib == "BernoulliInverse"){
		SimTK::Real randomNumber_Unif;
		int randomSign;

		randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomSign = int(std::floor(randomNumber_Unif * 2.0) - 1.0);
		if(randomSign < 0){
			passedVariable = 1.0 / passedVariable;
		}
	}

	// Bernoulli trial between the var and its reciprocal
	else if(distrib == "BernoulliReciprocal"){
		SimTK::Real randomNumber_Unif;
		int randomSign;

		randomNumber_Unif = uniformRealDistribution(randomEngine);
		randomSign = int(std::floor(randomNumber_Unif * 2.0) - 1.0);
		if(randomSign < 0){
			passedVariable = -1.0 * passedVariable;
		}
	}

	// Run Gaussian distribution
	else if(distrib == "normal"){

		passedVariable = passedVariable + (gaurand(randomEngine) * param1);    
	}

	// Convolute a Gaussian distribution with the mean equal to the passed
	// variable
	else if(distrib == "truncNormal"){

		// The passed variable becomes the mean of a truncated normal
		SimTK::Real mean = passedVariable;
		bool flag = false;

		// Try to get a hit inside the limits
		for (int tz = 0; tz < 100; tz++){

			// Add Gaussian noise with std = param1
			passedVariable = passedVariable + (gaurand(randomEngine) * param1);

			// Check if the passed variable exceeds the limits
			if((passedVariable >= param2) && (passedVariable <= param3)){
				flag = true;
				break;
			}
		}

		// Return nan if no sample was generated
		if(!flag){
			passedVariable = SimTK::NaN;
		}

		/* if(var <= 0){
			var = mean;
		}
		if(var >= (2*mean)){
			var = mean ;
		} */

	}

	// Run bimodal Gaussian distribution
	else if(distrib == "bimodalNormal"){

		passedVariable = convoluteVariable(passedVariable, "BernoulliInverse");
		passedVariable = passedVariable + gaurand(randomEngine);

	}
	
	// Gamma distribution
	else if(distrib == "gamma"){
		passedVariable = gammarand(randomEngine);
	}else{
		std::cerr << "Sampler distribute variable: Unkown distribution\n";
	}

	return passedVariable;
}

SimTK::Real Sampler::convoluteVariable(
	std::vector<SimTK::Real>& vvar,
	std::string distrib,
	SimTK::Real param1, SimTK::Real param2)
{
	for(auto& var : vvar){
		convoluteVariable(var, distrib, param1, param2);
	}

	// TODO
    assert(!"What should we return here?");
    return std::numeric_limits<SimTK::Real>::min();
}

SimTK::Real Sampler::calcDeformationPotential(
		SimTK::Real& var,
		std::string distrib,
		SimTK::Real param1, SimTK::Real param2
)
{
	SimTK::Real retVal = 0.0;

	// Bernoulli trial between the var and its inverse
	if(distrib == "BernoulliInverse"){
		return 0;
	}
	// Bernoulli trial between the var and its reciprocal
	else if(distrib == "BernoulliReciprocal"){
		return 0;
	}
	// Run Gaussian distribution
	else if(distrib == "normal"){
		return 0;
	}
	// Run truncated Gaussian distribution
	else if(distrib == "truncNormal"){
		SimTK::Real var2 = var*var;

		SimTK::Real firstTerm = 0.0;
		SimTK::Real secondTerm = 0.0;

		firstTerm = (1.0 / (var2)) - var2;
		secondTerm = 2.0 * param1 * (var - (1.0/var));

		SimTK::Real numerator = firstTerm + secondTerm;
		SimTK::Real denominator = param2*param2;

		retVal = (1.0 / this->beta) * 0.5 * (numerator / denominator);

	}
	// Run bimodal Gaussian distribution
	else if(distrib == "bimodalNormal"){
		return 0;
	}
	// Gamma distribution
	else if(distrib == "gamma"){
		return 0;
	}

	return retVal;
	
}




