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

	if((extnu != 0) && (extnu != 3) && (extnu != 6)){
		std::cerr << "External Fixman torque has been implemented only for 0, 3 and 6 dofs.\n";
		std::cerr << "Others are just experimental.\n";
		//throw std::exception();
		//std::exit(1);
	}

	if((extnu == 6) || (extnu == 3)){ // Free or Ball mobilizer

		bool HEAVY_PRINT = false;

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
		SimTK::Real w = 0;
		SimTK::Real x = 0;
		SimTK::Real y = 0;
		SimTK::Real z = 0;

		if(extnu == 6){
			const SimTK::Vec7& extQ = ((SimTK::MobilizedBody::Free *)mobodptr)->getQ(state);
			w = extQ[0];
			x = extQ[1];
			y = extQ[2];
			z = extQ[3];
		}else if(extnu == 3){
			const SimTK::Vec4& extQ = ((SimTK::MobilizedBody::Ball *)mobodptr)->getQ(state);
			w = extQ[0];
			x = extQ[1];
			y = extQ[2];
			z = extQ[3];
		}
		//std::cout << "extQ= " << extQ << std::endl;

		// Normalize
		SimTK::Real quatnorm = std::sqrt((w*w) + (x*x) + (y*y) + (z*z));
		w /= quatnorm; x /= quatnorm; y /= quatnorm; z /= quatnorm;
		SimTK::Quaternion quat(w, x, y, z);

		// Euler angles Verification
		SimTK::Rotation R = SimTK::Rotation(quat);
		SimTK::Real ww = w*w;
		SimTK::Real xx = x*x;
		SimTK::Real yy = y*y;
		SimTK::Real zz = z*z;
		SimTK::Real phi = std::atan2(2*((w*x) + (y*z)), 1 - (2*(xx + yy)));
		SimTK::Real theta = std::asin(2*((w*y) - (z*x)));
		SimTK::Real psi = std::atan2(2*((w*z) + (x*y)), 1 - (2*(yy + zz)));

		// Rotation matrices
		//SimTK::Real s1 = std::sin(phi);   SimTK::Real c1 = std::cos(phi);
		//SimTK::Real s2 = std::sin(theta); SimTK::Real c2 = std::cos(theta);
		//SimTK::Real s3 = std::sin(psi);   SimTK::Real c3 = std::cos(psi);

		//std::cout << "Fixman torque Rotation Tait–Bryan Z1Y2X3 " << std::endl;
		//std::cout << c1*c2 << " " << (c1*s2*s3) - (c3*s1) << " " << (s1*s3) + (c1*c3*s2) << std::endl;
		//std::cout << c2*s1 << " " << (c1*c3) + (s1*s2*s3) << " " << (c3*s1*s2) - (c1*s3) << std::endl;
		//std::cout << -1*s2 << " " << c2*s3 << " " << c2*c3 << std::endl;

		//std::cout << "Fixman torque Rotation" << R << std::endl;
		//std::cout << "Fixman torque Rotation matrix " << R[2][2] << " " << R[1][2] << " " << R[0][2] << std::endl;
		//std::cout << "Fixman torque Rotation matrix " << R[2][1] << " " << R[1][1] << " " << R[0][1] << std::endl;
		//std::cout << "Fixman torque Rotation matrix " << R[2][0] << " " << R[1][0] << " " << R[0][0] << std::endl;

		// Test tan function
		//std::cout << "Test asin function\n";
		//for (float i = -1; i <= 1; i += 0.1){
		//	std::cout << "asin " << i << " " << std::asin(i) << std::endl;
		//}

		// Fixman torque
		//SimTK::Real sinPitch = -1.0 * R[2][0];
		SimTK::Real sinPitch = 2*((w*y) - (z*x)); // SAME and FASTER
		// std::cout << -1.0 * R[2][0] << " ?= " << *((w*y) - (z*x)); << std::endl;
		SimTK::Real pitch = std::asin(sinPitch); 
		SimTK::Real cosAsinPitch = std::cos(pitch);
		SimTK::Real cotPitch = 1 / std::tan(pitch);
		SimTK::Real torqueComponent = 0;


		// Fixman potential
		SimTK::Real extFixPot = std::log(sinPitch * sinPitch);
		// External Fixman potential cutoff
		if(extFixPot < -14.0){ // Around double precision log(0)
			extFixPot = -14.0;
		}
		if(HEAVY_PRINT){
			std::cout << "External Fixman potential= " << extFixPot << std::endl;
		}
		///////////////////////////////////////////


		// Angle is around 0
		SimTK::Real torUppLim =  1000;
		SimTK::Real torLowLim = -1000;
		if((sinPitch > 0) && (sinPitch <= SimTK::TinyReal)){
			torqueComponent = torUppLim;
		}else if((sinPitch >= (-1 * SimTK::TinyReal)) && (sinPitch < 0)){
			torqueComponent = torLowLim;
		}else{
			//SimTK::Real torqueComponent = cosAsinPitch / sinPitch; // wrong
			//SimTK::Real torqueComponent = 1 * cotPitch; // wrong
			torqueComponent = cosAsinPitch / sinPitch; // wrong
		}

		if(HEAVY_PRINT){
			std::cout << "before pitch sinPitch cotPitch cutcot torque " << pitch << " " << sinPitch << " " << cotPitch << " " ;
		}

		// Torque is to high/low log(14) e(14)
		if(torqueComponent > torUppLim){
			torqueComponent = torUppLim;
		}else if(torqueComponent < torLowLim){
			torqueComponent = torLowLim;
		}

		if(HEAVY_PRINT){
			std::cout << torqueComponent << " " ;
		}

		torqueComponent *= (1.0) * RT; // internal Fixman sign reversed: Mbat is in the denominator
		
		if(HEAVY_PRINT){
			std::cout << torqueComponent << std::endl;
		}

		//torqueComponent = 0.0; // TODO: DELETE TODO:

		int uslot = 1; // Pitch 
		//std::cout << "Applying external Fixman torque " << torqueComponent << " to slot " << uslot << ".\n";

		mobod.applyOneMobilityForce(state, uslot, torqueComponent, mobilityForces);

	}

//	else if(extnu == 3){ // Cartesian
//		const SimTK::MobilizedBody* mobodptr = &mobod;
//
//		const SimTK::Vec3& extQ = ((SimTK::MobilizedBody::Cartesian *)mobodptr)->getQ(state);
//		SimTK::Real x = extQ[0];
//		SimTK::Real y = extQ[1];
//		SimTK::Real z = extQ[2];
//		
//		SimTK::Real sinPitch = std::sin(y);
//		SimTK::Real pitch = std::asin(sinPitch); 
//		SimTK::Real cosAsinPitch = std::cos(pitch);
//		SimTK::Real cotPitch = 1 / std::tan(pitch);
//	
//		SimTK::Real torqueComponent = cotPitch; // wrong
//
//		std::cout << "before pitch sinPitch cotPitch cutcot torque " << pitch << " " << sinPitch << " " << cotPitch << " " ;
//		// Cutoff based on statistics (within 2 stds)
//		if(torqueComponent > 10){
//			torqueComponent = 10;
//		}else if(torqueComponent < -10){
//			torqueComponent = -10;
//		}
//
//		std::cout << torqueComponent << " " ;
//
//		torqueComponent = 0.0; // TODO: DELETE TODO:
//
//		torqueComponent *= (1.0) * RT; // internal Fixman sign reversed: Mbat is in the denominator
//		
//		std::cout << torqueComponent << std::endl;
//
//	}

	else if(extnu == 1){ // Pin
		
		SimTK::Real alpha = mobod.getOneQ(state, 0);
	}

	else if(extnu == 0){ // Weld
//		const SimTK::Transform & T = mobod.getMobilizerTransform(state);
//		SimTK::Rotation iniR = T.R();
//		std::cout << "initial Rotation " << iniR[2][2] << " " << iniR[1][2] << " " << iniR[0][2] << std::endl;
//		std::cout << "initial Rotation " << iniR[2][1] << " " << iniR[1][1] << " " << iniR[0][1] << std::endl;
//		std::cout << "initial Rotation " << iniR[2][0] << " " << iniR[1][0] << " " << iniR[0][0] << std::endl;
//
//		SimTK::Quaternion quat = T.R().convertRotationToQuaternion();
//
//		SimTK::Real w = quat[0];
//		SimTK::Real x = quat[1];
//		SimTK::Real y = quat[2];
//		SimTK::Real z = quat[3];
//
//		// Normalize
//		SimTK::Real quatnorm = std::sqrt((w*w) + (x*x) + (y*y) + (z*z));
//		w /= quatnorm; x /= quatnorm; y /= quatnorm; z /= quatnorm;
//		quat = SimTK::Quaternion(w, x, y, z);
//
//		// Euler angles Verification
//		SimTK::Rotation R = SimTK::Rotation(quat);
//		SimTK::Real ww = w*w;
//		SimTK::Real xx = x*x;
//		SimTK::Real yy = y*y;
//		SimTK::Real zz = z*z;
//		SimTK::Real phi = std::atan2(2*((w*x) + (y*z)), 1 - (2*(xx + yy)));
//		SimTK::Real theta = std::asin(2*((w*y) - (z*x)));
//		SimTK::Real psi = std::atan2(2*((w*z) + (x*y)), 1 - (2*(yy + zz)));
//
//		// Rotation matrices
//		SimTK::Real s1 = std::sin(phi);   SimTK::Real c1 = std::cos(phi);
//		SimTK::Real s2 = std::sin(theta); SimTK::Real c2 = std::cos(theta);
//		SimTK::Real s3 = std::sin(psi);   SimTK::Real c3 = std::cos(psi);
//
//		std::cout << "Fixman torque Rotation Tait–Bryan Z1Y2X3 " << std::endl;
//		std::cout << c1*c2 << " " << (c1*s2*s3) - (c3*s1) << " " << (s1*s3) + (c1*c3*s2) << std::endl;
//		std::cout << c2*s1 << " " << (c1*c3) + (s1*s2*s3) << " " << (c3*s1*s2) - (c1*s3) << std::endl;
//		std::cout << -1*s2 << " " << c2*s3 << " " << c2*c3 << std::endl;
//
//		//std::cout << "Fixman torque Rotation" << R << std::endl;
//		std::cout << "Fixman torque Rotation matrix " << R[2][2] << " " << R[1][2] << " " << R[0][2] << std::endl;
//		std::cout << "Fixman torque Rotation matrix " << R[2][1] << " " << R[1][1] << " " << R[0][1] << std::endl;
//		std::cout << "Fixman torque Rotation matrix " << R[2][0] << " " << R[1][0] << " " << R[0][0] << std::endl;
//
//		// Test tan function
//		//std::cout << "Test asin function\n";
//		//for (float i = -1; i <= 1; i += 0.1){
//		//	std::cout << "asin " << i << " " << std::asin(i) << std::endl;
//		//}
//
//		// Fixman torque
//		SimTK::Real sinPitch = -1.0 * R[2][0];
//		//SimTK::Real sinPitch = 2*((w*y) - (z*x)); // SAME and FASTER
//		// std::cout << -1.0 * R[2][0] << " ?= " << *((w*y) - (z*x)); << std::endl;
//		SimTK::Real pitch = std::asin(sinPitch); 
//		SimTK::Real cosAsinPitch = std::cos(pitch);
//		SimTK::Real cotPitch = 1 / std::tan(pitch);
//	
//		//SimTK::Real torqueComponent = cosAsinPitch / sinPitch; // wrong
//		SimTK::Real torqueComponent = cotPitch; // wrong
//
//		std::cout << "before pitch sinPitch cotPitch cutcot torque " << pitch << " " << sinPitch << " " << cotPitch << " " ;
//
//		// Cutoff based on statistics (within 2 stds)
//		if(torqueComponent > 10){
//			torqueComponent = 10;
//		}else if(torqueComponent < -10){
//			torqueComponent = -10;
//		}
//
//		std::cout << torqueComponent << " " ;
//
//		//torqueComponent = 0.0; // TODO: DELETE TODO:
//
//		torqueComponent *= (1.0) * RT; // internal Fixman sign reversed: Mbat is in the denominator
//		
//		std::cout << torqueComponent << std::endl;
//
//		int uslot = 1; // Pitch 
//		//std::cout << "Applying external Fixman torque " << torqueComponent << " to slot " << uslot << ".\n";
//
//		mobod.applyOneMobilityForce(state, uslot, torqueComponent, mobilityForces);
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
