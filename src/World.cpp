#include "World.hpp"



// TODO write a pdb writer for all the Compounds
// TODO move this in Topology since they work only for one Compound
void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix)
{
  double mult = 10000*advanced.getTime(); // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(const SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::Compound& c, SimTK::State& advanced,
		 const char *dirname, const char *prefix, int midlength, const char *sufix, double aTime)
{
  double mult = 10000*aTime; // pico to femto
  SimTK::PdbStructure  pdb(advanced, c);
  std::stringstream sstream;
  sstream<<dirname<<"/"<<prefix<<decimal_prefix(mult, std::pow(10, midlength))<<int(mult)<<sufix<<".pdb";
  string ofilename = sstream.str();
  std::filebuf fb;
  fb.open(ofilename.c_str(), std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); // automatically multiplies by ten (nm to A)
  fb.close();
}

void writePdb(SimTK::PdbStructure pdb, const char *FN)
{
  std::filebuf fb;
  fb.open(FN, std::ios::out);
  std::ostream os(&fb);
  pdb.write(os); //automatically multiplies by ten (nm to A)
  fb.close();
}

/** Print a Compound Cartesian coordinates as given by
	 * Compound::calcAtomLocationInGroundFrame **/
void World::printPoss(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vertex;
	std::cout<<"Positions:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
	}
}

void World::printVels(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vel;
	std::cout<<"Velocities:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
}

void World::printPossVels(const SimTK::Compound& c, SimTK::State& advanced)
{
	SimTK::Vec3 vertex, vel;
	std::cout<<"Positions:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vertex   = c.calcAtomLocationInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vertex[0]<<" "<<vertex[1]<<" "<<vertex[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
	std::cout<<"Velocities:"<<std::endl;
	for (SimTK::Compound::AtomIndex aIx(0); aIx < c.getNumAtoms(); ++aIx){
		vel      = c.calcAtomVelocityInGroundFrame(advanced, aIx);
		std::cout<<c.getAtomName(aIx)<<"="<<std::setprecision(8)<<std::fixed
		  <<"["<<vel[0]<<" "<<vel[1]<<" "<<vel[2]<<"]"<<std::endl;
	}
	std::cout<<std::endl;
}

/** Constructor. Initializes the following pretty much empty objects:
 *  - CompoundSystem,
 *	  - SimbodyMatterSubsystem, GeneralForceSubsystem, DecorationSubsystem,
 *		Visualizer, Visualizer::Reporter, DuMMForceFieldSubsystem,
 *  - Integrator with a TimeStepper on top **/
World::World(int worldIndex, bool isVisual, SimTK::Real visualizerFrequency)
{
	// Get an index from a higher caller
	ownWorldIndex = worldIndex;

	// Initialize Fixman torque flag
	_useFixmanTorque = false;
  
	// Molmodel System derived from Simbody System
	compoundSystem = new SimTK::CompoundSystem;

	// Simbody subsystems (Minimum requirements)
	matter = new SimTK::SimbodyMatterSubsystem(*compoundSystem);

	// It appears to be used only for getMultibodySystem.calcPotentialEnergy
	forces = new SimTK::GeneralForceSubsystem(*compoundSystem);

	// Initialize Molmodel default ForceSubsystem (DuMM)
	forceField = new SimTK::DuMMForceFieldSubsystem(*compoundSystem);

	// Contact system
	tracker = new ContactTrackerSubsystem(*compoundSystem);
	contactForces = new CompliantContactSubsystem(*compoundSystem, *tracker);
	//contactForces = new HuntCrossleyContact();
	contactForces->setTrackDissipatedEnergy(true);
	contactForces->setTransitionVelocity(1e-2);

	// Intialize an integrator and a TimeStepper to manage it
	integ = new SimTK::VerletIntegrator(*compoundSystem);
	ts = new SimTK::TimeStepper(*compoundSystem, *integ);

	// Set the visual flag and if true initialize a Decorations Subsystem,
	// a Visualizer and a Simbody EventReporter which interacts with the
	// Visualizer
	this->visual = isVisual;
	if(visual){
		decorations = new SimTK::DecorationSubsystem(*compoundSystem);
		visualizer = new SimTK::Visualizer(*compoundSystem);
		visualizerReporter = new SimTK::Visualizer::Reporter(*visualizer
				, visualizerFrequency);
		compoundSystem->addEventReporter( visualizerReporter );

		// Initialize a DecorationGenerator
		paraMolecularDecorator = new ParaMolecularDecorator(
				compoundSystem,
				matter,
				forceField,
				forces
		);

		visualizer->addDecorationGenerator(paraMolecularDecorator);

	}

	// Statistics
	moleculeCount = -1;
	nofSamples = 0;

	// Thermodynamics
	this->temperature = -1; // this leads to unusal behaviour hopefully

        clique1 = ContactSurface::createNewContactClique();
}

/** Creates Gmolmodel topologies objects and based on amberReader forcefield
 * adds parameters: defines Biotypes; - adds BAT parameters to DuMM. Also
 * creates decorations for visualizers **/
void World::AddMolecule(
		readAmberInput *amberReader,
		std::string rbFN,
		std::string flexFN,
		std::string regimenSpec,
		std::string argRoot)
{
	// Statistics
	moleculeCount++; // Used for unique names of molecules
 
	// Add a new molecule (Topology object which inherits Compound)
	// to the vector of molecules.
	// TODO: Why resName and moleculeName have to be the same?
	std::string moleculeName = regimenSpec + std::to_string(moleculeCount);
	Topology *top = new Topology(moleculeName);
	topologies.emplace_back(top);

	// Set atoms properties from a reader: number, name, element, initial
	// name, force field type, charge, coordinates, mass, LJ parameters
	(topologies.back())->SetGmolAtomPropertiesFromReader(amberReader);

	// Set bonds properties from reader: bond indeces, atom neighbours
	(topologies.back())->SetGmolBondingPropertiesFromReader(amberReader);

	// Set atoms Molmodel types (Compound::SingleAtom derived) based on
	// their valence
	//(topologies.back())->SetGmolAtomsMolmodelTypes();
	(topologies.back())->SetGmolAtomsMolmodelTypesTrial();

	// Add parameters from amberReader
	(topologies.back())->bAddAllParams(amberReader, *forceField);

	// Build the graph representing molecule's topology
	(topologies.back())->buildGraphAndMatchCoords(*forceField, std::stoi(argRoot));

	(topologies.back())->loadTriples();

	// Set flexibility according to the flexibility file
	(topologies.back())->setFlexibility(regimenSpec, flexFN);
	//(topologies.back())->PrintAtomList();

	// Set generalized velocities scale factors 
	(topologies.back())->setUScaleFactorsToBonds(flexFN);
	// Print Molmodel types
	//(topologies.back())->PrintMolmodelAndDuMMTypes(*forceField);

	// All ocate the vector of coordinates (DCD)
	Xs.resize(Xs.size() + topologies.back()->getNAtoms());
	Ys.resize(Ys.size() + topologies.back()->getNAtoms());
	Zs.resize(Zs.size() + topologies.back()->getNAtoms());
  
	// Add Topology to CompoundSystem and realize topology
	compoundSystem->adoptCompound( *(topologies.back()) );
	//std::cout<<"World::AddMolecule CompoundSystem adoptCompound "<< std::endl;
	(topologies.back())->setCompoundIndex(
			SimTK::CompoundSystem::CompoundIndex(
			 compoundSystem->getNumCompounds() - 1));

	// Add the new molecule to Decorators's vector of molecules
	if(visual){
		// We need copy here.
		paraMolecularDecorator->AddMolecule(topologies.back());
	}

}

/**  **/
void World::setUScaleFactorsToMobods(void)
{

	for(const auto& Topology : topologies){
		// Iterate bonds
		std::vector<bSpecificAtom>&  Atoms = Topology->bAtomList;
		//for(const auto& AtomList : Topology->bAtomList){
		for(const auto& Bond : Topology->bonds){
			SimTK::Compound::AtomIndex aIx1 = Atoms[Bond.i].getCompoundAtomIndex();
			SimTK::Compound::AtomIndex aIx2 = Atoms[Bond.j].getCompoundAtomIndex();

			SimTK::MobilizedBodyIndex mbx1 = Topology->getAtomMobilizedBodyIndex(aIx1);
			SimTK::MobilizedBodyIndex mbx2 = Topology->getAtomMobilizedBodyIndex(aIx2);

			const SimTK::MobilizedBody& mobod1 = matter->getMobilizedBody(mbx1);
			const SimTK::MobilizedBody& mobod2 = matter->getMobilizedBody(mbx2);

			int level1 = mobod1.getLevelInMultibodyTree();
			int level2 = mobod2.getLevelInMultibodyTree();

			if(level1 > level2){
				mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx1, Bond.getUScaleFactor()));
			}else if(level2 > level1){
				mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx2, Bond.getUScaleFactor()));
			}else{
				if(Bond.getUScaleFactor() != 0){
					std::cout << "World::setUScaleFactorsToMobods Warning: Trying to scale a bond inside a rigid body\n"; 
				}
			}


		}
	}
/*
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		Topology * molecule = topologies[i];
		bSpecificAtom *Atoms = molecule->bAtomList;
		bBond *Bonds = molecule->bonds;

		// Iterate bonds
		for(int j = 0; j < molecule->nbonds; j++){
			SimTK::Compound::AtomIndex aIx1 = Atoms[ Bonds[j].i ].getCompoundAtomIndex();
			SimTK::Compound::AtomIndex aIx2 = Atoms[ Bonds[j].j ].getCompoundAtomIndex();

			SimTK::MobilizedBodyIndex mbx = molecule->getAtomMobilizedBodyIndex(aIx1);
			SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
			Bonds[j].getUScaleFactor();
		}

		// Get parentship - use aIx2Mbx 
		// Add entry to mbx2uScale
	}
*/
}

/** Get U scale factor for the mobilized body **/
const float World::getMobodUScaleFactor(SimTK::MobilizedBodyIndex& mbx) const
{
	if(mbx2uScale.find(mbx) != mbx2uScale.end()){
		return mbx2uScale.at(mbx);
	}else{
		std::cout << "Warning: U scale factor for mobod " << int(mbx) << " not found.\n";
		return 1;
	}
}
//...............


/**  **/
void World::addMembrane(SimTK::Real xWidth, SimTK::Real yWidth, SimTK::Real zWidth, int resolution)
{

	SimTK::Real stiffness = 10000.0;
	SimTK::Real dissipation = 0.0;
	SimTK::Real staticFriction = 0.0;
	SimTK::Real dynamicFriction = 0.0;
	SimTK::Real viscousFriction = 0.0;

	//ContactGeometry::HalfSpace contactGeometry;
	//ContactGeometry::Sphere contactGeometry(1);


	PolygonalMesh mesh = PolygonalMesh::createBrickMesh(Vec3(xWidth, yWidth, zWidth), resolution);
	ContactGeometry::TriangleMesh contactGeometry(mesh);

	matter->Ground().updBody().addContactSurface(
		Transform(Rotation(-0.5 * SimTK::Pi, SimTK::ZAxis))
		//Transform()
		, ContactSurface(
			contactGeometry
			, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction)
			, 1.0)
	);

	DecorativeFrame contactGeometryDecoFrame;
	matter->Ground().updBody().addDecoration(Transform(), contactGeometryDecoFrame);

	//DecorativeSphere contactGeometryDeco(1);
	//matter->Ground().updBody().addDecoration(Transform(), contactGeometryDeco);

	//DecorativeMesh contactGeometryDeco(mesh);
	DecorativeMesh contactGeometryDeco(contactGeometry.createPolygonalMesh());
	matter->Ground().updBody().addDecoration(Transform(), contactGeometryDeco.setColor(Cyan).setOpacity(0.5));

}


// Get the number of molecules
int World::getNofMolecules(void)
{
	return (this->moleculeCount + 1);
}

/** Calls CompoundSystem.modelOneCompound which links the Compounds to the
 * Simbody subsystems and realizes Topology. To be called after setting all
 * Compounds properties. **/
void World::modelTopologies(std::string GroundToCompoundMobilizerType)
{
	// Model the Compounds one by one in case we want to attach different types
	// of Mobilizers to the Ground in the feature.
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		//SimTK::String GroundToCompoundMobilizerType = "Free";
		//SimTK::String GroundToCompoundMobilizerType = "Weld";
		//SimTK::String GroundToCompoundMobilizerType = "Cartesian";


		this->rootMobility = GroundToCompoundMobilizerType;

		//std::cout<<"World::ModelTopologies call to CompoundSystem::modelCompound " << i
		//         << " grounded with mobilizer " << GroundToCompoundMobilizerType << std::endl;

	if(i == 0){ // First compound
			compoundSystem->modelOneCompound(
					SimTK::CompoundSystem::CompoundIndex(i),
					GroundToCompoundMobilizerType);
	}else{
			compoundSystem->modelOneCompound(
					SimTK::CompoundSystem::CompoundIndex(i),
					"Free");
	}

		// Realize Topology
		//compoundSystem->realizeTopology(); // restore MULMOL
		//((this->topologies)[i])->loadMobodsRelatedMaps(); // restore MULMOL

	}

	// Realize Topology
	//compoundSystem->realizeTopology();
}


const SimTK::State& World::realizeTopology(void)
{
		const SimTK::State& returnState = compoundSystem->realizeTopology();

		//for ( unsigned int i = 0; i < this->topologies.size(); i++){
		//    ((this->topologies)[i])->loadMobodsRelatedMaps();
		//}

	return returnState;
}

const SimTK::State& World::addContacts(void)
{
//	for ( unsigned int i = 0; i < this->topologies.size(); i++){
//		SimTK::MobilizedBodyIndex mbx = ((this->topologies)[i])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(0));
//		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
//	}
	const Real stiffness = 10000.0; // stiffness in pascals
	const Real dissipation = 0.0;    // to turn off dissipation
	SimTK::Real staticFriction = 0.0;
	SimTK::Real dynamicFriction = 0.0;
	SimTK::Real viscousFriction = 0.0;

	int nofContactAtomIxs = 7;
	int contAIxs[nofContactAtomIxs] = {0, 405, 3099, 1545, 3927, 5047, 2282};
	//for ( int aIx = 0; aIx < nofContactAtomIxs; aIx++){
	//}
		Vec3 halfSize(0.3, 0.3, 0.3);

		SimTK::MobilizedBodyIndex 
		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[0]));
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco(mesh.createPolygonalMesh());
		mobod.updBody().addDecoration(Transform(), deco.setColor(Cyan).setOpacity(.6));
		mobod.updBody().addContactSurface(Transform(), ContactSurface(mesh, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[1])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(0));
		SimTK::MobilizedBody& mobod7 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh7(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco7(mesh7.createPolygonalMesh());
		mobod7.updBody().addDecoration(Transform(), deco7.setColor(Cyan).setOpacity(.6));
		mobod7.updBody().addContactSurface(Transform(), ContactSurface(mesh7, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));
/*
		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[1]));
		SimTK::MobilizedBody& mobod1 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh1(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco1(mesh1.createPolygonalMesh());
		mobod1.updBody().addDecoration(Transform(), deco1.setColor(Cyan).setOpacity(.6));
		mobod1.updBody().addContactSurface(Transform(), ContactSurface(mesh1, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[2]));
		SimTK::MobilizedBody& mobod2 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh2(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco2(mesh2.createPolygonalMesh());
		mobod2.updBody().addDecoration(Transform(), deco2.setColor(Cyan).setOpacity(.6));
		mobod2.updBody().addContactSurface(Transform(), ContactSurface(mesh2, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[3]));
		SimTK::MobilizedBody& mobod3 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh3(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco3(mesh3.createPolygonalMesh());
		mobod3.updBody().addDecoration(Transform(), deco3.setColor(Cyan).setOpacity(.6));
		mobod3.updBody().addContactSurface(Transform(), ContactSurface(mesh3, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[4]));
		SimTK::MobilizedBody& mobod4 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh4(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco4(mesh4.createPolygonalMesh());
		mobod4.updBody().addDecoration(Transform(), deco4.setColor(Cyan).setOpacity(.6));
		mobod4.updBody().addContactSurface(Transform(), ContactSurface(mesh4, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[5]));
		SimTK::MobilizedBody& mobod5 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh5(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco5(mesh5.createPolygonalMesh());
		mobod5.updBody().addDecoration(Transform(), deco5.setColor(Cyan).setOpacity(.6));
		mobod5.updBody().addContactSurface(Transform(), ContactSurface(mesh5, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));

		mbx = ((this->topologies)[0])->getAtomMobilizedBodyIndex(SimTK::Compound::AtomIndex(contAIxs[6]));
		SimTK::MobilizedBody& mobod6 = matter->updMobilizedBody(mbx);
		ContactGeometry::TriangleMesh mesh6(PolygonalMesh::createSphereMesh(0.3, 2));
		DecorativeMesh deco6(mesh6.createPolygonalMesh());
		mobod6.updBody().addDecoration(Transform(), deco6.setColor(Cyan).setOpacity(.6));
		mobod6.updBody().addContactSurface(Transform(), ContactSurface(mesh6, ContactMaterial(stiffness, dissipation, staticFriction, dynamicFriction, viscousFriction), 0.001));
*/


	const SimTK::State& returnState = compoundSystem->realizeTopology();
	return returnState;
}

void World::loadCompoundRelatedMaps(void)
{
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		((this->topologies)[i])->loadCompoundAtomIx2GmolAtomIx();
		//((this->topologies)[i])->printMaps();
	}
}


void World::loadMobodsRelatedMaps(void)
{
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		((this->topologies)[i])->loadMobodsRelatedMaps();
		//((this->topologies)[i])->printMaps();
	}
}

/** Set up Fixman torque **/
void World::addFixmanTorque()
{
	// Set flag
	_useFixmanTorque = true;

	// Alloc memory for FixmanTorque implementation and add to forces
	FixmanTorqueImpl = new FixmanTorque(compoundSystem, *matter);
	ExtForce = new Force::Custom(*forces, FixmanTorqueImpl);
}

/** Check if the Fixman torque flag is set **/
bool World::isUsingFixmanTorque(void)
{
	return _useFixmanTorque;
}

/** Get writble pointer to FixmanTorque implementation **/
FixmanTorque * World::updFixmanTorque(void)
{
	return FixmanTorqueImpl;
}

/** Get pointer to FixmanTorque implementation **/
FixmanTorque * World::getFixmanTorque(void) const
{
	return FixmanTorqueImpl;
}

// ----------------------
// --- Thermodynamics ---
// ----------------------

/** Get the World temperature **/
SimTK::Real World::getTemperature(void)
{
	return this->temperature;
}

/** Set this World temperature but also ths samplers and 
Fixman torque temperature. **/
void World::setTemperature(SimTK::Real argTemperature)
{
	// Set the temperature for the samplers
	for(unsigned int samplerIx = 0; samplerIx < samplers.size(); samplerIx++){
		samplers[samplerIx]->setTemperature(argTemperature);
	}

	// Set the temperature for this World
	this->temperature = argTemperature;

	// Set the temperature for the Fixman torque also
	if(_useFixmanTorque){
		FixmanTorqueImpl->setTemperature(this->temperature);
	}
}
//...............

//...................
// --- Simulation ---
//...................

/** Get/Set seed for reproducibility. **/
void World::setSeed(int whichSampler, unsigned long long int argSeed)
{
	samplers[whichSampler]->setSeed(argSeed);
}

unsigned long long int World::getSeed(int whichSampler)
{
	return samplers[whichSampler]->getSeed();
}


/** Amber like scale factors. **/
void World::setAmberForceFieldScaleFactors()
{
	forceField->setVdw12ScaleFactor(0.0);
	forceField->setVdw13ScaleFactor(0.0);
	forceField->setVdw14ScaleFactor(0.5); // RESTORE from OpenMM
	//forceField->setVdw14ScaleFactor(0.0); // for OpenMM
	forceField->setVdw15ScaleFactor(1.0);

//* RESTORE SAFETY  
	forceField->setCoulomb12ScaleFactor(0.0);
	forceField->setCoulomb13ScaleFactor(0.0);
	forceField->setCoulomb14ScaleFactor(0.8333333333); // RESTORE from OpenMM
	//forceField->setCoulomb14ScaleFactor(0.0); // for OpenMM
	forceField->setCoulomb15ScaleFactor(1.0);
	//forceField->setVdwMixingRule(
	//       SimTK::DuMMForceFieldSubsystem::LorentzBerthelot); */
 
   /* DANGER ! no electrostatics
	forceField->setCoulomb12ScaleFactor(0.0);
	forceField->setCoulomb13ScaleFactor(0.0);
	forceField->setCoulomb14ScaleFactor(0.0);
	forceField->setCoulomb15ScaleFactor(0.0);
	// */
}

/** Set a global scaling factor for all the terms in the forcefield **/
void World::setGlobalForceFieldScaleFactor(SimTK::Real scaleFactor)
{
  
	forceField->setBondStretchGlobalScaleFactor(scaleFactor);
	forceField->setBondBendGlobalScaleFactor(scaleFactor);
	forceField->setBondTorsionGlobalScaleFactor(scaleFactor);
	forceField->setAmberImproperTorsionGlobalScaleFactor(scaleFactor);
  
	forceField->setVdw12ScaleFactor(scaleFactor);
	forceField->setVdw13ScaleFactor(scaleFactor);
	forceField->setVdw14ScaleFactor(scaleFactor);
	forceField->setVdw15ScaleFactor(scaleFactor);
	forceField->setVdwGlobalScaleFactor(scaleFactor);
  
	forceField->setCoulomb12ScaleFactor(scaleFactor);
	forceField->setCoulomb13ScaleFactor(scaleFactor);
	forceField->setCoulomb14ScaleFactor(scaleFactor);
	forceField->setCoulomb15ScaleFactor(scaleFactor);
	forceField->setCoulombGlobalScaleFactor(scaleFactor);

}

/** Set GBSA implicit solvent scale factor. **/
void World::setGbsaGlobalScaleFactor(SimTK::Real scaleFactor)
{
	forceField->setGbsaGlobalScaleFactor(scaleFactor);
}

/** Get a writeble pointer to the DuMM force field **/
SimTK::DuMMForceFieldSubsystem * World::updForceField()
{
	return forceField;
}

//...................
// --- Statistics ---
//...................

/** How many samples do we have so far **/
int World::getNofSamples(void)
{
	// Zero it every time the user asks
	nofSamples = 0;
 
	// Gather samples from all the samplers
	for(size_t i = 0; i < samplers.size(); i++){
		nofSamples += (samplers[i])->getNofSamples();
	}

	return this->nofSamples;
}

/** How many Samplers does this World have. **/
int World::getNofSamplers(void)
{
	return samplers.size();
}

/** Add a sampler to this World using the specialized struct
for samplers names. **/
BaseSampler * World::addSampler(SamplerName samplerName)
{
	BaseSampler *p = NULL;
	if(samplerName == HMC){

        p = new HMCSampler(this,
                compoundSystem, matter, topologies,
                forceField, forces, ts
                );
        samplers.push_back(p);

	}else if(samplerName == LAHMC){

        p = new LAHMCSampler(this,
                compoundSystem, matter, topologies,
                forceField, forces, ts, 4
                );
        samplers.push_back(p);

	}

	//return samplers.size();
	return p;
}

// Get a sampler based on its position in the samplers vector
// TODO Use ampler polymorphism
const BaseSampler * World::getSampler(int which)
{
	return samplers[which];
}

// Get a writable sampler based on its position in the samplers vector
// TODO Use Sampler polymorphism
BaseSampler * World::updSampler(int which)
{
	return samplers[which];
}


/** Get a const reference to a molecule **/
const Topology& World::getTopology(int moleculeNumber) const{
	return *(topologies[moleculeNumber]);
}

/** Get a writble reference to the last molecule. **/
Topology& World::updTopology(int moleculeNumber){
	//return *(topologies.back());
	return *(topologies[moleculeNumber]);
}

/** Return a 2D vector representing all the coordinates of this World.
 * The first dimension represents the molecules (topologies) and the second
 * dimension (inner) represents the coordinates. The second inner dimension
 * type is pair of bSpecificAtom* and a Vec3. Thus, besides coordinates, it
 * contains all the information in bSpecificAtom as well. The bottleneck here
 * is the calcAtomLocationInGroundFrame from Compound.
 **/
std::vector < std::vector < std::pair <bSpecificAtom *, SimTK::Vec3> > >
		World::getAtomsLocationsInGround(const SimTK::State & state)
{
	std::vector< std::vector< std::pair<bSpecificAtom *, SimTK::Vec3> > >
			returnVector;

	// Iterate through topologies
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		// Iterate through atoms
		for(int j = 0; j < (topologies[i])->getNumAtoms(); j++){
			SimTK::Compound::AtomIndex compoundAtomIndex
				= topologies[i]->bAtomList[j].getCompoundAtomIndex();

			SimTK::Vec3 location = (topologies[i])
				    ->calcAtomLocationInGroundFrame(state, compoundAtomIndex);

			currentTopologyInfo.emplace_back(
				    &((topologies[i])->bAtomList[j]), location);
		}
		returnVector.emplace_back(currentTopologyInfo);
	}

	return returnVector;
}

/** Put coordinates into bAtomLists of Topologies.
 * When provided with a State, calcAtomLocationInGroundFrame
 * realizes Position and uses matter to calculate locations **/
void World::updateAtomListsFromCompound(const SimTK::State &state)
{
	// Iterate through topologies
	for ( unsigned int i = 0; i < this->topologies.size(); i++){
		// Iterate through atoms
		for(int j = 0; j < topologies[i]->getNumAtoms(); j++){
			SimTK::Compound::AtomIndex compoundAtomIndex = topologies[i]->bAtomList[j].getCompoundAtomIndex();
			SimTK::Vec3 location = topologies[i]->calcAtomLocationInGroundFrame(state, compoundAtomIndex);
			topologies[i]->bAtomList[j].setX(location[0]);
			topologies[i]->bAtomList[j].setY(location[1]);
			topologies[i]->bAtomList[j].setZ(location[2]);
		}
	}
}


/** Set Compound, MultibodySystem and DuMM configurations according to
some other World's atoms.
A body is composed of a root atom and other periferic atoms which have their own stations.
Unles is a fully flexible Cartesian world, the function has the following steps:
1. Set Compound
2. Set DuMM
3. Set Simbody bodies
	3.1. Transforms X_PF and X_BM
	3.2. Mass properties

**/ 
SimTK::State& World::setAtomsLocationsInGround(
		SimTK::State& someState,
		std::vector<std::vector<
		std::pair<bSpecificAtom *, SimTK::Vec3> > > otherWorldsAtomsLocations)
{

	// Get the total no of bodies in this world (each World has its own
	// SimbodyMatterSubsystem)
	int totalNofBodies = matter->getNumBodies();

	// Arrays of Transforms
	SimTK::Transform G_X_T;
	//SimTK::Transform T_X_root[totalNofBodies];

	// Loop through molecules/topologies
	for(unsigned int i = 0; i < otherWorldsAtomsLocations.size(); i++) {

		// Set the decorator
		if (visual == true) {
			paraMolecularDecorator->setAtomTargets(otherWorldsAtomsLocations[i]);
		}

		/////////////////
		// 1. COMPOUND
		/////////////////

		// Use Molmodel's Compound match functions to set the new 
		// Compound transforms
		std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
		for(unsigned int j = 0; j < otherWorldsAtomsLocations[i].size(); j++){
			SimTK::Compound::AtomIndex atomIndex = ((otherWorldsAtomsLocations[i][j]).first)->atomIndex;
			SimTK::Vec3 location = ((otherWorldsAtomsLocations[i][j]).second);
			atomTargets.insert(std::pair<SimTK::Compound::AtomIndex, SimTK::Vec3>(atomIndex, location));
		}

		topologies[i]->matchDefaultAtomChirality(atomTargets, 0.01, false);
		topologies[i]->matchDefaultBondLengths(atomTargets);
		topologies[i]->matchDefaultBondAngles(atomTargets);
		topologies[i]->matchDefaultDirections(atomTargets);
		topologies[i]->matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
		topologies[i]->matchDefaultTopLevelTransform(atomTargets);

		// Get the Ground to Top Transform
		G_X_T = topologies[i]->getTopLevelTransform();

		// Recalculate atom frames in top compound frames
		topologies[i]->calcTopTransforms();

		/////////////////
		// FULLY FLEXIBLE CARTESIAN WORLD
		// Parent body is always Ground
		/////////////////
		if(topologies[i]->getRegimen().at(0) == 'I'){

			/////////////////
			// 3. SIMBODY MATTER 
			//---------------
			// 3.1 MOBILIZED BODIES
			/////////////////

			// Parent (Ground) to mobile point M (atom) on body
			SimTK::Transform P_X_M[totalNofBodies]; // related to X_PFs

			// P_X_Ms are set to atom positions in Ground 
			for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
				const auto mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
				P_X_M[int(mbx)] = Transform(Rotation(), atomTargets[aIx]); 
			}

			// Set X_FM: backs up to stage Time
			for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				SimTK::Transform P_X_F = mobod.getInboardFrame(someState);
				SimTK::Transform F_X_P = ~P_X_F;
				SimTK::Transform F_X_M = F_X_P * P_X_M[int(mbx)];

				(mobod).setQToFitTransform(someState, F_X_M);
			}

		/////////////////
		// NON-CARTESIAN JOINT TYPE WORLDs
		// Bodies are linked in a tree with P-F-M-B links
		/////////////////
		}else{

			/////////////////
			// 2. DUMM 
			/////////////////
			// Get locations for DuMM
			SimTK::Vec3 locs[topologies[i]->getNumAtoms()];

			// Loop through atoms
			for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
				if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin
						locs[int(aIx)] = SimTK::Vec3(0); 
				}else{ // atom is not at body's origin
					SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
					SimTK::Transform G_X_root = G_X_T * topologies[i]->getTopTransform((topologies[i]->mbx2aIx)[mbx]);

					SimTK::Vec3 G_vchild = atomTargets[aIx];

					SimTK::Transform root_X_child = alignFlipAndTranslateFrameAlongXAxis(G_X_root, G_vchild);
					topologies[i]->bsetFrameInMobilizedBodyFrame(aIx, root_X_child);

					locs[int(aIx)] = root_X_child.p();
				}
			}

			// Set stations and AtomPLacements for atoms in DuMM
			// Loop through atoms
			for (SimTK::Compound::AtomIndex aIx(1); aIx < topologies[i]->getNumAtoms(); ++aIx){
				SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
					SimTK::DuMM::AtomIndex dAIx = topologies[i]->getDuMMAtomIndex(aIx);

					// Set station_B
					forceField->bsetAtomStationOnBody( dAIx, locs[int(aIx)] );
					forceField->bsetAllAtomStationOnBody( dAIx, locs[int(aIx)] ); // full

					// Set included atom
					forceField->updIncludedAtomStation(dAIx) = (locs[int(aIx)]);
					forceField->updAllAtomStation(dAIx) = (locs[int(aIx)]); // full

					// Atom placements in clusters
					forceField->bsetAtomPlacementStation(dAIx, mbx, locs[int(aIx)] );
			} /////////////////////////


			/////////////////
			// 3. SIMBODY MATTER 
			//---------------
			// 3.1 MOBILIZED BODIES TRANSFORMS
			/////////////////

			// Set X_PF and X_BM in Mobilized bodies
			// First get the atom 0 transform. Works for multiple mols because
			// every mol has Ground as parent
			SimTK::Transform P_X_F_1 = G_X_T * topologies[i]->getTopTransform(SimTK::Compound::AtomIndex(0)); // is this always true ??

			// Loop through the rest of the atoms and get P_X_F, B_X_M from the Compound transforms
			for (SimTK::Compound::AtomIndex aIx(0); aIx < topologies[i]->getNumAtoms(); ++aIx){
				if(topologies[i]->getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin

					// Get body, parentBody, parentAtom
					SimTK::MobilizedBodyIndex mbx = topologies[i]->getAtomMobilizedBodyIndex(aIx);
					SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
					const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
					SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

					if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground

						// Get mobod transforms
						std::vector<SimTK::Transform> mobodTs = topologies[i]->calcMobodTransforms(matter, aIx, someState);

						mobod.setDefaultInboardFrame(mobodTs[0]);
						mobod.setDefaultOutboardFrame(mobodTs[1]);

					} // END if parent not Ground
					else{ // parent is Ground
						mobod.setDefaultInboardFrame(P_X_F_1);
						mobod.setDefaultOutboardFrame(Transform());
					}
				} // END atom is at body's origin
			} // END loop through atoms


			/////////////////
			// 3. SIMBODY MATTER 
			//---------------
			// 3.2 MASS PROPERTIES
			/////////////////
			// Set mass properties for mobilized bodies
			// Loop through mobilized bodies
			for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
				DuMM::ClusterIndex clusterIx = forceField->bgetMobodClusterIndex(mbx);
				SimTK::MassProperties massProperties = forceField->calcClusterMassProperties(clusterIx);
				mobod.setDefaultMassProperties(massProperties);
			}

			this->compoundSystem->realizeTopology();
			someState = compoundSystem->updDefaultState();


		} // END TD regimen and all regimens

	} // END iterating through molecules/topologies

	this->compoundSystem->realize(someState, SimTK::Stage::Position);

	updateAtomListsFromCompound(someState);

	return someState;

}


/** Print information about Simbody systems. For debugging purpose. **/
void World::PrintSimbodyStateCache(SimTK::State& someState){
	std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
	for(int i = 0; i < someState.getNumSubsystems(); i++){
		std::cout << " Subsystem " << i << " Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
			<< " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
			<< " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i)) << std::endl;
	}
}

/** Fill Worlds Cartesian coordinates buffers. 
To be called before use of getXs, getYs or getZs **/
void World::updateCoordBuffers(void)
{
	int allAtIx = -1;
	for(unsigned int tIx = 0; tIx < topologies.size(); tIx++){
		for(int aIx = 0; aIx < topologies[tIx]->getNAtoms(); aIx++){
			allAtIx++;
			Xs[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getX();
			Ys[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getY();
			Zs[allAtIx] = (topologies[tIx]->bAtomList[aIx]).getZ();
		}
	}

}

/** Get the coordinates from buffers **/
std::vector<SimTK::Real> World::getXs(void)
{
	return Xs;
}

std::vector<SimTK::Real> World::getYs(void)
{
	return Ys;
}

std::vector<SimTK::Real> World::getZs(void)
{
	return Zs;
}

/** Destructor **/
World::~World(){
    if(this->visual == true){
        //delete paraMolecularDecorator;

        delete decorations;
		decorations = nullptr;

        delete visualizer;
		visualizer = nullptr;
        //delete vizReporter;
    }
    delete ts; ts = nullptr;
    delete integ; integ = nullptr;
    delete forceField; forceField = nullptr;
    if(_useFixmanTorque){
        delete ExtForce; ExtForce = nullptr;
    }
    delete matter; matter = nullptr;
    delete forces; forces = nullptr;
    delete compoundSystem; compoundSystem = nullptr;

    for(unsigned int i = 0; i < topologies.size(); i++){
        delete topologies[i]; topologies[i] = nullptr;
    }
    for(unsigned int i = 0; i < samplers.size(); i++){
        delete samplers[i]; samplers[i] = nullptr;
    }
}

/** Return own CompoundSystem **/
CompoundSystem *World::getCompoundSystem() const {
	return compoundSystem;
}

/** Set own Compound system **/
// TODO find a solution for the old one
void World::setCompoundSystem(CompoundSystem *compoundSystem) {
	World::compoundSystem = compoundSystem;
}



