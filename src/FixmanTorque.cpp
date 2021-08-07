#include "Robo.hpp"
#include "FixmanTorque.hpp"


////////////////////////////
////// FIXMAN TORQUE //////////
////////////////////////////
FixmanTorque::FixmanTorque(SimTK::CompoundSystem *argCompoundSystem, SimTK::SimbodyMatterSubsystem& argMatter
					) : matter(argMatter){
	this->compoundSystem = argCompoundSystem;
	scaleFactor = 1.0;
	this->temperature = 0.0;
	this->RT = temperature * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD); // TODO double vs long double
	
}

void FixmanTorque::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>&,
						   SimTK::Vector_<SimTK::Vec3>&, SimTK::Vector& mobilityForces) const  
{
	// function args were
	// const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>& bodyForces, SimTK::Vector_<SimTK::Vec3>& particleForces, SimTK::Vector& mobilityForces

	// Compute Fixman torque
	int nu = state.getNU();
	SimTK::Vector V3(nu);
	SimTK::Vector V4(nu);
	SimTK::Real D0 = 1.0;
	matter.calcFixmanTorque(state, V3, V4, &D0);
	// end - Compute Fixman torque

	// Calculate geometric features fast
	//xstd::cout << " 1 " ;
	//xfor (SimTK::MobilizedBodyIndex mbx(2); mbx < matter.getNumBodies(); ++mbx){
	//x    const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);
	//x    const SimTK::MobilizedBody *p_mobod = &mobod;
	//x    std::cout << std::setprecision(10) << std::fixed <<  ((SimTK::MobilizedBody::Pin *)(p_mobod))->getAngle(state) << ' ' ;
	//x}

	//std::cout << " FT " ;
	int uslot = -1;
	for (SimTK::MobilizedBodyIndex mbx(0); mbx < matter.getNumBodies(); ++mbx){
		const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

		for(int k = 0; k < mobod.getNumU(state); k++){
			uslot++;

			//mobod.applyOneMobilityForce(state, k, scaleFactor * V4[uslot], mobilityForces);

			mobod.applyOneMobilityForce(state, k, (-1.0) * RT * V4[uslot], mobilityForces);
			//mobod.applyOneMobilityForce(state, k, (-1.0) * V4[uslot], mobilityForces);

			//xif(int(mbx) == 2){
			//    std::cout << std::setprecision(6) << std::fixed << (-1.0) * RT * V4[uslot] << ' '; 
			//x}
			//    << " to mbx " << std::setprecision(0) << int(mbx) << " slot " << uslot << " ";

			//std::cout << "Fixman torque scaleFactor " << scaleFactor << " " ;
			//std::cout << " temperature " << temperature << " -RT " << (-1.0) * RT << std::endl ;
		}
	}
	//std::cout << std::endl ;

	//const SimTK::Real q = knee.getOneQ(state, 0);
	//const SimTK::Real x = q < low ? q-low : (q > high ? q-high : 0);
	//knee.applyOneMobilityForce(state, 0, -k*x, mobilityForces);
}

FixmanTorque::~FixmanTorque(){}

// This should be carefully analyzed. Intended to be taken from somewhere else.
SimTK::Real FixmanTorque::calcPotentialEnergy(const SimTK::State&) const {
	// function args were const SimTK::State& state
	SimTK::Real energy = 0.0;
	return energy;
}

bool FixmanTorque::dependsOnlyOnPositions() const {
	return true;
}

SimTK::Real FixmanTorque::getTemperature(void)
{
	return this->temperature;
}

void FixmanTorque::setTemperature(SimTK::Real argTemperature)
{
	this->temperature = argTemperature;
	this->RT = temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

SimTK::Real FixmanTorque::getScaleFactor(void)
{
	return scaleFactor;
}

void FixmanTorque::setScaleFactor(SimTK::Real argScaleFactor)
{
	scaleFactor = argScaleFactor;
	std::cout << "FixmanTorque::setScaleFactor("<< this->scaleFactor << ")" << std::endl;
}



///////////////////////////////
////// END FIXMAN TORQUE //////
///////////////////////////////

////////////////////////////////////////
////// FIXMAN TORQUE EXTERNAL //////////
///////////////////////////////////////
FixmanTorqueExt::FixmanTorqueExt(SimTK::CompoundSystem *argCompoundSystem,
	SimTK::SimbodyMatterSubsystem& argMatter
	) : matter(argMatter){

	this->compoundSystem = argCompoundSystem;
	scaleFactor = 1.0;
	this->temperature = 0.0;
	this->RT = temperature * static_cast<SimTK::Real>(SimTK_BOLTZMANN_CONSTANT_MD); // TODO double vs long double
	
}

void FixmanTorqueExt::calcForce(const SimTK::State& state, SimTK::Vector_<SimTK::SpatialVec>&,
				SimTK::Vector_<SimTK::Vec3>&, SimTK::Vector& mobilityForces) const  
{
	// function args were

	// Compute Fixman torque
	int nu = state.getNU();

	// Get First mobod 
	const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(SimTK::MobilizedBodyIndex(1));
	int extnu = mobod.getNumU(state);

	if((extnu != 0) && (extnu != 6)){
		std::cerr << "External Fixman torque has been implemented only for 0 and 6 dofs.\n";
		throw std::exception();
		std::exit(1);
	}

	if(extnu == 6){ // Free mobilizer

		const SimTK::MobilizedBody* mobodptr = &mobod;
		
		// Get quaternion from root atom
		//bSpecificAtom *root = &(bAtomList[bSpecificAtomRootIndex]);
		//SimTK::Compound::AtomIndex aIx = root->getCompoundAtomIndex();
		//SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, aIx);
		//SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();
		////std::cout << "calcLogDetMBATGamma2Contribution quaternion " << quat << std::endl;
	
		//SimTK::Real w = quat[0];
		//SimTK::Real x = quat[1];
		//SimTK::Real y = quat[2];
		//SimTK::Real z = quat[3];
		//SimTK::Real sinPitch = 2 * (w * y - z * x);

		// Alternatively get quaternion from mobod 1
		const SimTK::Vec7& extQ = ((SimTK::MobilizedBody::Free *)mobodptr)->getQ(state);
		SimTK::Real w = extQ[0];
		SimTK::Real x = extQ[1];
		SimTK::Real y = extQ[2];
		SimTK::Real z = extQ[3];
		SimTK::Real sinPitch = 2 * (w * y - z * x);

		SimTK::Real cosPitch = 1 - (sinPitch * sinPitch);

		//std::cout << "External Q = ";
		//for(int j = 0; j < 7; j++){
		//	std::cout << extQ[j] << " ";
		//}
		//std::cout << std::endl;
	
		SimTK::Real torqueComponent = cosPitch / sinPitch;

		if(torqueComponent == SimTK::Infinity){
			torqueComponent = 10;
		}else if(torqueComponent == SimTK::Infinity){
			torqueComponent = -10;
		}

		//std::cout << "sinPitch, cosPitch torqueComponent " << sinPitch << " " << cosPitch << " " << torqueComponent << std::endl;

		torqueComponent *= (1.0) * RT; // internal Fixman sign reversed: Mbat is in the denominator
		
		int uslot = 1; // Pitch 
		//std::cout << "Applying external Fixman torque " << torqueComponent << " to slot " << uslot << ".\n";

		mobod.applyOneMobilityForce(state, uslot, torqueComponent, mobilityForces);
	}

}

FixmanTorqueExt::~FixmanTorqueExt(){}

// This should be carefully analyzed. Intended to be taken from somewhere else.
SimTK::Real FixmanTorqueExt::calcPotentialEnergy(const SimTK::State&) const {
	// function args were const SimTK::State& state
	SimTK::Real energy = 0.0;
	return energy;
}

bool FixmanTorqueExt::dependsOnlyOnPositions() const {
	return true;
}

SimTK::Real FixmanTorqueExt::getTemperature(void)
{
	return this->temperature;
}

void FixmanTorqueExt::setTemperature(SimTK::Real argTemperature)
{
	std::cout << "Setting T for external Fixman torque.\n";
	this->temperature = argTemperature;
	this->RT = temperature * SimTK_BOLTZMANN_CONSTANT_MD;
}

SimTK::Real FixmanTorqueExt::getScaleFactor(void)
{
	return scaleFactor;
}

void FixmanTorqueExt::setScaleFactor(SimTK::Real argScaleFactor)
{
	scaleFactor = argScaleFactor;
	std::cout << "FixmanTorqueExt::setScaleFactor("<< this->scaleFactor << ")" << std::endl;
}



////////////////////////////
////// END EXTERNAL FIXMAN TORQUE //////
////////////////////////////
