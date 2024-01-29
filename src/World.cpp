#include "World.hpp"

#define CHECK_RECONSTRUCTION false

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
			//c.calcAtomLocationInGroundFrameThroughSimbody(
			//	aIx, *forceField, *matter, advanced);
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


//==============================================================================
//                   CLASS TaskSpace
//==============================================================================
/**
 *  Contains a Symbody task space and additional data
 **/
StationTaskLaurentiu::StationTaskLaurentiu(void)
{

}

//==============================================================================
//                             1. STRUCTURAL FUNCTIONS
//==============================================================================

/** Constructor. Initializes the following objects:
 *  - CompoundSystem,
 *	  - SimbodyMatterSubsystem, GeneralForceSubsystem, DecorationSubsystem,
 *		Visualizer, Visualizer::Reporter, DuMMForceFieldSubsystem,
 *  - Integrator with a TimeStepper on top **/
World::World(int worldIndex,
	int requestedNofMols,
	bool isVisual,
	SimTK::Real visualizerFrequency)
{
	// Get an index from a higher caller
	ownWorldIndex = worldIndex;

	// Molmodel System derived from Simbody System
	compoundSystem = std::make_unique<SimTK::CompoundSystem>();

	// Simbody subsystems (Minimum requirements)
	matter = std::make_unique<SimTK::SimbodyMatterSubsystem>(*compoundSystem);

	// It appears to be used only for getMultibodySystem.calcPotentialEnergy
	forces = std::make_unique<SimTK::GeneralForceSubsystem>(*compoundSystem);

	// Initialize Molmodel default ForceSubsystem (DuMM)
	forceField = std::make_unique<SimTK::DuMMForceFieldSubsystem>(*compoundSystem);

	// Contact system
	/*
	tracker = std::make_unique<ContactTrackerSubsystem>(*compoundSystem);
	contactForces = std::make_unique<CompliantContactSubsystem>(*compoundSystem, *tracker);
	contactForces->setTrackDissipatedEnergy(true);
	contactForces->setTransitionVelocity(1e-2);
    	clique1 = ContactSurface::createNewContactClique();
	*/

	// Intialize an integrator and a TimeStepper to manage it
	integ = std::make_unique<SimTK::VerletIntegrator>(*compoundSystem);
	ts = std::make_unique<SimTK::TimeStepper>(*compoundSystem, *integ);

	const SimTK::MultibodySystem& mbs = forces->getMultibodySystem();
	
	// Set the visual flag and if true initialize a Decorations Subsystem,
	// a Visualizer and a Simbody EventReporter which interacts with the
	// Visualizer
	this->visual = isVisual;
	if(visual){
		decorations = std::make_unique<SimTK::DecorationSubsystem>(*compoundSystem);
		visualizer = std::make_unique<SimTK::Visualizer>(*compoundSystem);
		visualizerReporter = std::make_unique<SimTK::Visualizer::Reporter>(
			*visualizer, std::abs(visualizerFrequency));

		compoundSystem->addEventReporter(visualizerReporter.get());

		if(contactForces){
			std::cout << "[WARNING] Victor check Teodor's contacts." << std::endl;
			visualizer->addDecorationGenerator(
				new ForceArrowGenerator(mbs, *contactForces));
		}else{
			std::cout << "[WARNING] Teodor's contacts." << std::endl;
		}
		

		// Initialize a DecorationGenerator
		paraMolecularDecorator = std::make_unique<ParaMolecularDecorator>(
			compoundSystem.get(),
			matter.get(),
			forceField.get(),
			forces.get()
		);

		visualizer->addDecorationGenerator(paraMolecularDecorator.get());
	}

	// Statistics
	moleculeCount = -1;

	// Thermodynamics
	this->temperature = -1; // this leads to unusal behaviour hopefully

	//topologies.reserve(requestedNofMols); // SAFE

	// currStage = this->getSimbodyMatterSubsystem()->getStage(
	// 	this->getSimbodyMatterSubsystem()->getSystem().getDefaultState()
	// );

}

/** Creates Gmolmodel topologies objects and based on amberReader forcefield
 * adds parameters: defines Biotypes; - adds BAT parameters to DuMM. Also
 * creates decorations for visualizers **/
void World::AddMolecule(
		readAmberInput *amberReader,
		//std::string rbFN,
		//std::string flexFN,
		//std::string regimenSpec,
		std::string argRoot
		//, std::string argRootMobility
		)
{
/*
	// Statistics
	moleculeCount++; // Used for unique names of molecules

	// Add a new molecule (Topology object which inherits Compound)
	// to the vector of molecules.
	// TODO: Why resName and moleculeName have to be the same?
	// TODO store molecule name in vector maybe
	//std::string moleculeName = regimenSpec + std::to_string(moleculeCount);
	std::string moleculeName = "MOL" + std::to_string(moleculeCount);

	roots.emplace_back(argRoot);

	//rootMobilities.emplace_back(argRootMobility); // TODO: move to setflexibilities
	rootMobilities.emplace_back("Pin"); // TODO: move to setflexibilities

	topologies.emplace_back(Topology{moleculeName}); // TODO is this ok?

	// Set atoms properties from a reader: number, name, element, initial
	// name, force field type, charge, coordinates, mass, LJ parameters
	topologies.back().SetGmolAtomPropertiesFromReader(amberReader);

	// Set bonds properties from reader: bond indeces, atom neighbours
	topologies.back().SetGmolBondingPropertiesFromReader(amberReader);

	// Set atoms Molmodel types (Compound::SingleAtom derived) based on
	// their valence
	//topologies.back().SetGmolAtomsMolmodelTypes();
	topologies.back().SetGmolAtomsMolmodelTypesTrial();
*/
}

// Add Biotypes
void World::AddBiotypes(int which, readAmberInput *amberReader)
{
/*
	//topologies.back().bAddBiotypes(amberReader); // SAFE
	topologies[which].bAddBiotypes(amberReader); // DANGER
*/
}

void World::generateDummParams(int which, readAmberInput *amberReader
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs)

{
	// Add DuMM parameters from amberReader
	((*topologies)[which]).generateDummParams(amberReader, *forceField,
	aClassParams2aClassId,
	allBondsACIxs, allAnglesACIxs, allDihedralsACIxs, allImpropersACIxs);
}

void World::transferDummParams(int which, readAmberInput *amberReader
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs
	)
{
	// Add DuMM parameters from amberReader
	((*topologies)[which]).transferDummParams(amberReader, *forceField,
		aClassParams2aClassId,
		allBondsACIxs, allAnglesACIxs, allDihedralsACIxs, allImpropersACIxs);
}

void World::BuildTopologyGraph(int which, std::string argRoot)
{
/*
	// Build the graph representing molecule's topology
	//topologies.back().buildGraphAndMatchCoords(std::stoi(argRoot)); // SAFE
	//topologies.back().loadTriples(); // SAFE
	topologies[which].buildGraphAndMatchCoords(std::stoi(argRoot)); // DANGER
	topologies[which].loadTriples(); // DANGER
*/
}

void World::AllocateCoordBuffers(int which)
{
	// All ocate the vector of coordinates (DCD)
	// TODO: where is this supposed to be?
	//Xs.resize(Xs.size() + topologies.back().getNAtoms()); // SAFE
	//Ys.resize(Ys.size() + topologies.back().getNAtoms()); // SAFE
	//Zs.resize(Zs.size() + topologies.back().getNAtoms()); // SAFE
	Xs.resize(Xs.size() + ((*topologies)[which]).getNAtoms()); // DANGER
	Ys.resize(Ys.size() + ((*topologies)[which]).getNAtoms()); // DANGER
	Zs.resize(Zs.size() + ((*topologies)[which]).getNAtoms()); // DANGER

}

/** It sets Compound BondFlexibilities . Also
 * creates decorations for visualizers
 **/
void World::SetBondFlexibilities(
		std::string flexFN,
		std::string regimenSpec,
		std::string argRootMobility,
		int which)
{
	// Set flexibility according to the flexibility file
	//((*topologies)[which]).PrintAtomList();
	((*topologies)[which]).setFlexibility(regimenSpec, flexFN, ownWorldIndex);

	rootMobilities[which] = argRootMobility; // NEW

	// Set generalized velocities scale factors
	((*topologies)[which]).setUScaleFactorsToBonds(flexFN);
	// Print Molmodel types
	//topologies.back().PrintMolmodelAndDuMMTypes(*forceField);

	// Allocate the vector of coordinates (DCD)
	// TODO: where is this supposed to be?
	//Xs.resize(Xs.size() + topologies[which].getNAtoms());
	//Ys.resize(Ys.size() + topologies[which].getNAtoms());
	//Zs.resize(Zs.size() + topologies[which].getNAtoms());
}

/** Adopts a topology **/
void World::adoptTopology(int which)
{
	// Add Topology to CompoundSystem and realize topology
	scout("World: adopting ") << which <<" " << "topology" << eol;
	compoundSystem->adoptCompound(((*topologies)[which]));

	// Sets the
	((*topologies)[which]).setCompoundIndex(
		SimTK::CompoundSystem::CompoundIndex(which));

	//topologies[which].setCompoundIndex(
	//		SimTK::CompoundSystem::CompoundIndex(
	//		 compoundSystem->getNumCompounds() - 1));

	// Add the Topology object to Decorators's vector of molecules
	if(visual){
		// We need copy here.
		//paraMolecularDecorator->AddMolecule(&(topologies[which])); // SAFE
		paraMolecularDecorator->AddMolecule( &((*topologies)[which]) ); // DANGER
	}

}

/** Calls CompoundSystem.modelOneCompound which links the Compounds to the
 * Simbody subsystems and realizes Topology. To be called after setting all
 * Compounds properties. **/

// UNUSABLE

void World::modelTopologies(std::string GroundToCompoundMobilizerType)
{
	// Model the Compounds one by one in case we want to attach different types
	// of Mobilizers to the Ground in the feature.
	for ( std::size_t i = 0; i < this->topologies->size(); i++){

			Topology& topology = (*topologies)[i];

			compoundSystem->modelOneCompound(
				SimTK::CompoundSystem::CompoundIndex(i),
				topology.atomFrameCache,
				rootMobilities[i]);

		 std::cout<<"World::ModelTopologies call to CompoundSystem::modelCompound " << i
		         << " grounded with mobilizer " << rootMobilities[i] << std::endl;

		//std::cout << "World::ModelTopologies " <<
		for(std::size_t k = 0; k < (*topologies)[i].getNumAtoms(); k++){
			SimTK::Compound::AtomIndex aIx = (((*topologies)[i]).bAtomList[k]).getCompoundAtomIndex();
			SimTK::MobilizedBodyIndex mbx = ((*topologies)[i]).getAtomMobilizedBodyIndex(aIx);
			std::cout << "k aIx " << k << " " << aIx << " " << mbx << std::endl << std::flush;
		}

		// // Realize Topology
		//compoundSystem->realizeTopology(); // restore MULMOL
		// ((this->topologies)[i])->loadMobodsRelatedMaps(); // restore MULMOL

	}

	// // Realize Topology
	// compoundSystem->realizeTopology();
}

// Print recommended timesteps. We need and advanced State here
void World::PrintInitialRecommendedTimesteps(void)
{
	SimTK::State& someState = integ->updAdvancedState();
	int nu = matter->getNU(someState);

	/* SimTK::Real CartesianTimestepInv = 1.0 / 0.001; // ps
	SimTK::Vector V(nu, CartesianTimestepInv);
	SimTK::Vector timesteps(nu);

	matter->multiplyBySqrtMInv(someState, V, timesteps);

	for (int j=0; j < nu; ++j){
		timesteps[j] = 1.0 / timesteps[j];
	}

	for (int j=0; j < nu; ++j){
		std::cout << timesteps[j] << " ";
	}
	std::cout << std::endl; */

	/* SimTK::Matrix M;
	matter->calcM(someState, M);
	std::cout << M << " "; */

	SimTK::Real CartesianTimestep = 0.001;
	SimTK::Vector V(matter->getNumBodies() - 1, CartesianTimestep);

	// Get Mobods mass properties
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		SimTK::SpatialInertia sp = mobod.getBodySpatialInertiaInGround(someState);

		V[int(mbx) - 1] *= sp.getMass();

	}

	std::cout << SimTK::min(V) << " ";

	/* for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		SimTK::SpatialInertia sp = mobod.getBodySpatialInertiaInGround(someState);

		std::cout << sp.getMass() << " ";

	} */

}



//==============================================================================
//                   TaskSpace Functions
//==============================================================================

/**
 * Allocate memory for a task space consisting of a set of body indeces
 * station on the bodies expresed in both guest (target) and host 
 * and the difference between them
*/
void World::addTaskSpaceLS(void)
{
	//StationTaskLaurentiu stationTask;

	int guestTopology = 1;
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	int topi = -1;
	for(auto& topology : (*topologies)){
		topi++;

		if(topi == guestTopology){

			// Guest atoms iteration
			for (int bAtomIx : bAtomIxs_guest) {
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);

				onBodyB.emplace_back(mbx);
				taskStationPInGuest.emplace_back(SimTK::Vec3());
				taskStationPInHost.emplace_back(SimTK::Vec3());
				taskDeltaStationP.emplace_back(SimTK::Vec3());

				if(this->visual == true){
					paraMolecularDecorator->loadArrow(SimTK::Vec3(0), SimTK::Vec3(0));
				}

			}
		}
	}
}

/**
 * Update target task space
*/
void World::updateTaskSpace(const State& someState)
{

	int hostTopology = 0;
	int guestTopology = 1;

	std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	// Get stations in host
	int topi = -1;
	for(auto& topology : (*topologies)){
		topi++;
		if(topi == hostTopology){

			// Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_host) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				taskStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);
			}
		}
	}

	// Get stations in guest
	topi = -1;
	for(auto& topology : (*topologies)){
		topi++;

		if(topi == guestTopology){

			// Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_guest) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).getCompoundAtomIndex();
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				taskStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);

				taskDeltaStationP[tz] = taskStationPInHost[tz] - taskStationPInGuest[tz];

				if(this->visual == true){
					paraMolecularDecorator->updateArrow(tz, taskStationPInGuest[tz], taskStationPInGuest[tz] + taskDeltaStationP[tz]);
				}
			}
		}
	}	

}

/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceStationPInGuest(void)
{
	return taskStationPInGuest;
}
/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceStationPInHost(void)
{
	return taskStationPInHost;
}

/**
 * Get the difference between the station task and the target
*/
SimTK::Array_<SimTK::Vec3>& 
World::getTaskSpaceDeltaStationP(void)
{
	return taskDeltaStationP;
}

/**
 * Calc Station Jacobian JS
*/
void World::calcStationJacobian(
	const State&                           someState,
	SimTK::Matrix_<SimTK::Vec3>&                      JS) const
{
		matter->calcStationJacobian(someState, onBodyB, taskStationPInGuest, JS);

		std::cout << "Task Bodies ";
		std::cout << onBodyB << std::endl;
		std::cout << "Task Stations ";
		std::cout << taskStationPInGuest << std::endl;
		std::cout << "Station Jacobian ";
		std::cout << JS << std::endl;

		//matter->calcBiasForStationJacobian(someState, onBodyB, stationPInB, JSDotu);
}


//=============================================================================
//                   CONSTRAINTS
//=============================================================================

/**
 * Add contact constraints to specific bodies.
 **/
void World::addRodConstraint(State& someState)
{

	int hostTopology = 0;
	int guestTopology = 1;

	std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology
	rodBodies.emplace_back(std::make_pair(SimTK::MobilizedBody(), SimTK::MobilizedBody()));
	conStationPInGuest.emplace_back(SimTK::Vec3());
	conStationPInHost.emplace_back(SimTK::Vec3());
	conDeltaStationP.emplace_back(SimTK::Vec3());

	// Get stations in host
	int topi = -1;
	for(auto& topology : (*topologies)){
		topi++;
		if(topi == hostTopology){

			/* // Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_host) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);

				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				rodBodies[0].first = mbx;

				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				conStationPInHost[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);
			} */
		}
	}

	// Get stations in guest
	topi = -1;
	for(auto& topology : (*topologies)){
		topi++;

		if(topi == guestTopology){

			/* // Atoms
			int tz = -1;
			for (int bAtomIx : bAtomIxs_guest) {
				tz++;
				SimTK::Compound::AtomIndex aIx = (topology.bAtomList[bAtomIx]).compoundAtomIndex;
				SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
				
				SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);

				rodBodies[0].second = mbx;
				
				SimTK::Transform X_GB = mobod.getBodyTransform(someState);
				SimTK::Vec3 B_aLoc = topology.getAtomLocationInMobilizedBodyFrame(aIx);
				conStationPInGuest[tz] = X_GB.p() + ((X_GB.R()) * B_aLoc);

			} */
		}
	}

/* 	std::cout << "Rod bodies: " << rodBodies[0].first
		<< " " << rodBodies[0].second << std::endl << std::flush;

	rodConstraints.emplace_back( SimTK::Constraint::Rod(
		matter->updMobilizedBody(rodBodies[0].first),  SimTK::Vec3(),
		matter->updMobilizedBody(rodBodies[0].second), SimTK::Vec3(), 0.1) ); */

	//rodConstraints.back().enable(someState);

}

/** Add speed constraints to specific bodies.
TODO:use number of mobilities. TODO: Solve if **/
const SimTK::State& World::addSpeedConstraint(int prmtopIndex)
{

	int hostTopology = 0;
	int guestTopology = 1;

	std::vector<int> bAtomIxs_host = {4}; // atoms on host topology
	std::vector<int> bAtomIxs_guest = {29}; // atoms on target topology

	if(prmtopIndex >= 0){
		std::cout << "Adding constraint to atom with prmtop index "
			<< prmtopIndex << "\n" ;
		SimTK::MobilizedBodyIndex mbx =
			((*topologies)[0]).getAtomMobilizedBodyIndexThroughDumm(
				SimTK::Compound::AtomIndex(prmtopIndex), *forceField
		);
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		SimTK::Constraint::ConstantSpeed B3291ConstraintU1(
			mobod, SimTK::MobilizerUIndex(0), 0);
		if(matter->getNumBodies() > 5000){
			SimTK::Constraint::ConstantSpeed B3291ConstraintU2(
				mobod, SimTK::MobilizerUIndex(1), 0);
			SimTK::Constraint::ConstantSpeed B3291ConstraintU3(
				mobod, SimTK::MobilizerUIndex(2), 0);
		}
	}

	const SimTK::State& returnState = compoundSystem->realizeTopology();
	return returnState;
}

//=============================================================================
//                   CONTACTS Functions
//=============================================================================

/** Add contact surfaces to bodies **/
void World::addContacts(const std::vector<int>& prmtopIndex, const int topologyIx, 
						const SimTK::ContactCliqueId cliqueId)
{
	// Iterate through the prmtopIndex and add the appropriate spheres

	for(int prmtopIx : (prmtopIndex)){
		if (prmtopIx == -1) {
			std::cout << "World::addContacts: prmtop index -1 found. Skipping topologyIx " << topologyIx << std::endl;
			break;
		}

		std::cout << "Adding contacts with membrane to atom with prmtop index " <<
			prmtopIx << " to topology " << topologyIx << " in clique " << cliqueId << std::endl;



		const Real stiffness = 10000; // stiffness in pascals
		const Real dissipation = 0.0; 
		SimTK::Real staticFriction = 0.0;
		SimTK::Real dynamicFriction = 0.0;
		SimTK::Real viscousFriction = 0.0;

		// Map from AmberAtomIndex to CompoundAtomIndex 
		SimTK::Compound::AtomIndex cAIx = ((*topologies)[topologyIx]).
										bAtomList[prmtopIx].getCompoundAtomIndex();

		SimTK::MobilizedBodyIndex
		mbx = ((*topologies)[topologyIx]).getAtomMobilizedBodyIndexThroughDumm(
			cAIx, *forceField);
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		
		// Need to get vdw radius for said atom
		const SimTK::DuMM::AtomIndex dAIx = ((*topologies)[topologyIx]).getDuMMAtomIndex(cAIx);
		Real vdwRadius = forceField->getAtomRadius(dAIx);

		std::cout << "Atom element: " << forceField->getAtomElement(dAIx)
				  << " Atom radius (nm): " << forceField->getAtomRadius(dAIx) << std::endl;

		PolygonalMesh sphereMesh;
    	sphereMesh = PolygonalMesh::createSphereMesh(vdwRadius,1);
			
		// We need to account for the position of the actual atom, so we
		// must also apply a translation from the Origin of the MoBod to the
		// location of the station coressponding to said atom.

		SimTK::Vec3 atomPos = ((*topologies)[topologyIx]).getAtomLocationInMobilizedBodyFrameThroughDumm
								(cAIx,getForceField());


		ContactGeometry::TriangleMesh sphere(sphereMesh);
		mobod.updBody().addContactSurface(Transform(atomPos), 
			ContactSurface(sphere,
				ContactMaterial(stiffness, dissipation,
			staticFriction, dynamicFriction, viscousFriction),
				0.1).joinClique(cliqueId));

		if (visual == true) {
			DecorativeMesh showSphere(sphere.createPolygonalMesh());
			mobod.updBody().addDecoration(Transform(atomPos), 
				showSphere.setColor(Cyan).setOpacity(0.2));
			//TODO: Remove this in production
			mobod.updBody().addDecoration(Transform(atomPos), 
				showSphere.setColor(Gray).setRepresentation(DecorativeGeometry::DrawWireframe));
		}
	}
}

// CONTACT DEBUG
/*
int numForces = updWorld(currentWorldIx)->contactForces->getNumContactForces(
	currentAdvancedState);
SimTK::Real dissEnergy = updWorld(currentWorldIx)->contactForces->
	getDissipatedEnergy(currentAdvancedState);
bool hasDefaultForceGenerator =
	updWorld(currentWorldIx)->contactForces->hasDefaultForceGenerator();

const MultibodySystem & mbs =
	updWorld(currentWorldIx)->contactForces->getMultibodySystem();
int nofMobods = mbs.getMatterSubsystem().getNumBodies();

const ContactTrackerSubsystem & cts =
	updWorld(currentWorldIx)->contactForces->getContactTrackerSubsystem();
int ctsNofSurfaces = cts.getNumSurfaces();

std::cout << "CONTACT INFO:"
	<< " #forces= " << numForces	const SimTK::State& addContacts(int prmtopIndex);

	const SimTK::State& realizeTopology();

	<< " dissEnergy= " << dissEnergy
	<< " hasDefaultForceGenerator= " << hasDefaultForceGenerator
	<< " #mobods= " << nofMobods
	<< " ctsNofSurfaces= " << ctsNofSurfaces
<< std::endl;
*/
// CONTACT DEBUG enD
/** Assign a scale factor for generalized velocities to every mobilized
body **/

//=============================================================================
//                   MEMBRANE Functions
//=============================================================================

/** Add a membrane represented by a contact surface **/
void World::addMembrane(const SimTK::Real halfThickness)
{

	SimTK::Real stiffness = 10000;
	SimTK::Real dissipation = 0;
	SimTK::Real staticFriction = 0.0;
	SimTK::Real dynamicFriction = 0.0;
	SimTK::Real viscousFriction = 0.0;

 	matter->Ground().updBody().addContactSurface(
		Transform(Rotation(-0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, halfThickness)),
		ContactSurface(
		ContactGeometry::HalfSpace(),
		ContactMaterial(stiffness, dissipation,
			staticFriction, dynamicFriction, viscousFriction))
			.joinClique(SimTK::ContactCliqueId(1))
			.joinClique(SimTK::ContactCliqueId(2))
			.joinClique(SimTK::ContactCliqueId(3))
			);

	matter->Ground().updBody().addContactSurface(
		Transform(Rotation(-0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, -halfThickness)),
		ContactSurface(
		ContactGeometry::HalfSpace(),
		ContactMaterial(stiffness, dissipation,
			staticFriction, dynamicFriction, viscousFriction))
			.joinClique(SimTK::ContactCliqueId(0))
			.joinClique(SimTK::ContactCliqueId(2))
			.joinClique(SimTK::ContactCliqueId(3)));
	
	matter->Ground().updBody().addContactSurface(
		Transform(Rotation(0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, halfThickness)),
		ContactSurface(
		ContactGeometry::HalfSpace(),
		ContactMaterial(stiffness, dissipation,
			staticFriction, dynamicFriction, viscousFriction))
			.joinClique(SimTK::ContactCliqueId(0))
			.joinClique(SimTK::ContactCliqueId(1))
			.joinClique(SimTK::ContactCliqueId(3))
			);

	matter->Ground().updBody().addContactSurface(
		Transform(Rotation(0.5 * SimTK::Pi, SimTK::YAxis), Vec3(0, 0, -halfThickness)),
		ContactSurface(
		ContactGeometry::HalfSpace(),
		ContactMaterial(stiffness, dissipation,
			staticFriction, dynamicFriction, viscousFriction))
			.joinClique(SimTK::ContactCliqueId(0))
			.joinClique(SimTK::ContactCliqueId(1))
			.joinClique(SimTK::ContactCliqueId(2))
			);


	if (visual == true) {
		DecorativeFrame contactGeometryDecoFrame;
		matter->Ground().updBody().addDecoration(
		Transform(),
        DecorativeBrick(Vec3(10,10,halfThickness)).setColor(Orange).setOpacity(0.25));
		
	}

}





/** Realize Topology for this World **/
const SimTK::State& World::realizeTopology()
{
	const SimTK::State& returnState = compoundSystem->realizeTopology();

	// for ( unsigned int i = 0; i < this->topologies.size(); i++){
	//    ((this->topologies)[i])->loadMobodsRelatedMaps();
	// }

	return returnState;
}

/** 
 * Assign a scale factor for generalized velocities to every mobilized
 * body
 */
void World::setUScaleFactorsToMobods(void)
{

	//for(auto& topology : topologies){ // SAFE
	for(auto& topology : (*topologies)){ // DANGER
		// Iterate bonds

		//for(const auto& AtomList : topology.bAtomList){
		for(auto& Bond : topology.bonds){
			SimTK::Compound::AtomIndex aIx1 = topology.bAtomList[Bond.i].getCompoundAtomIndex();
			SimTK::Compound::AtomIndex aIx2 = topology.bAtomList[Bond.j].getCompoundAtomIndex();

			// SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndexFromMap(aIx1, ownWorldIndex);
			// SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndexFromMap(aIx2, ownWorldIndex);

			SimTK::MobilizedBodyIndex mbx1 = topology.getAtomMobilizedBodyIndexThroughDumm(aIx1, *forceField);
			SimTK::MobilizedBodyIndex mbx2 = topology.getAtomMobilizedBodyIndexThroughDumm(aIx2, *forceField);

			//std::cout
			//<< "World::setUScaleFactorsToMobods aIx1 aIx2 mbx1 mbx2 "
			//<< aIx1 << " " << aIx2 << " "
			//<< mbx1 << " " << mbx2 << std::endl;

			const SimTK::MobilizedBody& mobod1 = matter->getMobilizedBody(mbx1);
			const SimTK::MobilizedBody& mobod2 = matter->getMobilizedBody(mbx2);

			int level1 = mobod1.getLevelInMultibodyTree();
			int level2 = mobod2.getLevelInMultibodyTree();

			if(level1 > level2){
				mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx1, Bond.getUScaleFactor(ownWorldIndex)));
			}else if(level2 > level1){
				mbx2uScale.insert( std::pair< SimTK::MobilizedBodyIndex, float > (mbx2, Bond.getUScaleFactor(ownWorldIndex)));
			}else{
				if(Bond.getUScaleFactor(ownWorldIndex) != 0){
					//std::cout << "World::setUScaleFactorsToMobods Warning: Trying to scale a bond inside a rigid body\n";
				}
			}


		}
	}
}

/** Load CompoundAtomIndex to Gmolmodel atom index map **/
void World::loadCompoundRelatedMaps()
{
	//for (auto& topology : topologies){ // SAFE
	for (auto& topology : (*topologies)){ // DANGER
		topology.loadCompoundAtomIx2GmolAtomIx();
		// topology.printMaps();
	}
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void World::loadMbx2AIxMap(){

	for (auto& topology : (*topologies)){

		// Iterate through atoms and get their MobilizedBodyIndeces
		for (unsigned int i = 0; i < topology.getNumAtoms(); ++i) {

			// Get atomIndex from atomList
			SimTK::Compound::AtomIndex aIx = (topology.bAtomList[i]).getCompoundAtomIndex();

			// Get MobilizedBodyIndex from CompoundAtom
			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);

			if(!aIx.isValid() || !mbx.isValid())
			{
				std::cout << "[WARNING]: Tried adding invalid index in loadMbx2AIxMap() aIx: " << aIx << ", mbx: " << mbx << std::endl;
				continue;
			}

			// Map mbx2aIx contains only atoms at the origin of mobods
			if (topology.getAtomLocationInMobilizedBodyFrame(aIx) == 0) {
				mbx2aIx.insert(std::make_pair(mbx, aIx));
			}

		} // atoms
	} // topologies

}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void World::loadMbx2AIxMap_SP_NEW(){

	for (auto& topology : (*topologies)){

		// Iterate through atoms and get their MobilizedBodyIndeces
		for (unsigned int aCnt = 0; aCnt < topology.getNumAtoms(); ++aCnt) {

			// Get atomIndex from atomList
			SimTK::Compound::AtomIndex aIx = (topology.subAtomList[aCnt]).getCompoundAtomIndex();

			// Get MobilizedBodyIndex from CompoundAtom
			SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(aIx);

			if(!aIx.isValid() || !mbx.isValid())
			{
				std::cout << "[WARNING]: Tried adding invalid index in loadMbx2AIxMap() aIx: " << aIx << ", mbx: " << mbx << std::endl;
				continue;
			}

			// Map mbx2aIx contains only atoms at the origin of mobods
			if (topology.getAtomLocationInMobilizedBodyFrame(aIx) == 0) {
				mbx2aIx.insert(std::make_pair(mbx, aIx));
			}

		} // atoms
	} // topologies

}

void World::loadMobodsRelatedMaps()
{
	//for (auto& topology : topologies){ // SAFE
	for (auto& topology : (*topologies)){ // DANGER
		topology.loadAIx2MbxMap();
		loadMbx2AIxMap();
		// topology.printMaps();
	}
}

// Allocate space for containers that keep statistics if we're doing any
void World::allocateStatsContainers(void)
{
	// Arccos of the X_PF first entry which should contain 
	// the angle of rotation on X for the BendStretch joint
	acosX_PF00.resize(matter->getNumBodies() - 1);
	acosX_PF00_means.resize(matter->getNumBodies() - 1);

	// The norm of the translation vector of X_BM which should
	// contain the bond length for the BendStretch joint
	normX_BMp.resize(matter->getNumBodies() - 1);
	normX_BMp_means.resize(matter->getNumBodies() - 1);
}

// Get the (potential) energy transfer
// If any of the Q, U or tau is actively modifyied by the sampler
// the Jacobian of that transformation will be included too
SimTK::Real World::getWorkOrHeat(void)
{
	// Accumulate in this variable
	SimTK::Real retValue = 0.0;

	// Get the energy transfer from all the samplers
	for(auto& sampler : this->samplers){

		// Get the potential energy difference
		retValue += 
			( getSampler(0)->getNewPE() - getSampler(0)->getOldPE() );

		// Get Fixman potential difference
		retValue +=
			( getSampler(0)->getNewFixman() - getSampler(0)->getOldFixman());

		/* // Get the Q modifying samplers Jacobians
		if(sampler->getDistortOpt() < 0){
			retValue -= 
				sampler->getDistortJacobianDetLog();
		} */
		
	}
	
	return retValue;

}

// Get the (potential) energy transfer in the form of work
// If any of the Q, U or tau is actively modifyied by the sampler
// the Jacobian of that transformation will be included too
SimTK::Real World::getWork(void)
{
	// Accumulate in this variable
	SimTK::Real retValue = 0.0;

	// Get the energy transfer from all the samplers
	for(auto& sampler : this->samplers){

		if(sampler->getDistortOpt() < 0){

			// Get the potential energy difference
			retValue += 
				( getSampler(0)->getNewPE() - getSampler(0)->getOldPE() );

			// Get Fixman potential difference
			retValue +=
				( getSampler(0)->getNewFixman() - getSampler(0)->getOldFixman());

			/* // Get the Jacobians
			retValue -= 
				sampler->getDistortJacobianDetLog(); */
		}
		
	}
	
	return retValue;

}

/*
 * Shift all the generalized coordinates
 */
void World::getTransformsStatistics(SimTK::State& someState)
{
	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const Transform& X_PF = mobod.getInboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const Transform& X_FM = mobod.getMobilizerTransform(someState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const Transform& X_BM = mobod.getOutboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		// Get BAT coordinate "angle"
		/*
		/ cos t
		| 
		|
		\
		*/
		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		/* std::cout << "World " << ownWorldIndex << " " 
			//<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			//<< "angleMean " << int(mbx) - 1 << " "
			//<< angleMean * (180 / SimTK::Pi) << " "
			//<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */

	}

}

// Print bond lengths and angle bends
void World::traceBendStretch(SimTK::State& someState){
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){
		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		
		// Get mobod inboard frame X_PF
		const Transform& X_PF = mobod.getInboardFrame(someState);
		//PrintTransform(X_PF, 4, "X_PF " + std::to_string(int(mbx)));

		// Get mobod inboard frame X_FM measured and expressed in P
		const Transform& X_FM = mobod.getMobilizerTransform(someState);

		// Get mobod inboard frame X_BM
		const Transform& X_BM = mobod.getOutboardFrame(someState);
		//PrintTransform(X_BM, 4, "X_BM " + std::to_string(int(mbx)));

		// Get BAT coordinates B and A
		SimTK::Vec3 bondVector = X_BM.p();
		trace( bondVector.norm(),  std::acos(X_PF.R()(0)(0)));
		//trace("X_FM");
		//PrintTransform(X_FM, 10);
	
	}
}

// Print X_PF means
void World::PrintAcosX_PFs(void)
{
	int i = -1;
	for(const auto &xpf : acosX_PF00 ){
		i += 1;
		//std::cout << "X_PF " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
		std::cout << "acosX_PF " << i << " " << xpf << std::endl;

	}
}

// Print X_PF means
void World::PrintNormX_BMs(void)
{
	int i = -1;
	for(const auto &xbm : normX_BMp ){
		i += 1;
		std::cout << "normX_BM " << i << " " << xbm << std::endl;
	}
}

// Print X_PF means
void World::PrintAcosX_PFMeans(void)
{
	int i = -1;
	for(const auto &xpf : acosX_PF00_means ){
		i += 1;
		//std::cout << "X_PFMean " << i << " " << xpf * (180 / SimTK::Pi) << std::endl;
		std::cout << "acosX_PFMean " << i << " " << xpf << std::endl;
	}
}

// Print X_PF means
void World::PrintNormX_BMMeans(void)
{
	int i = -1;
	for(const auto &xbm : normX_BMp_means ){
		i += 1;
		std::cout << "normX_BMMean " << i << " " << xbm << std::endl;
	}
}

// REORIENT

SimTK::Transform& World::getReorientTransformInAnotherBody(
	const State &someState,
	const MobilizedBody &inBodyA, const MobilizedBody &ofBodyB,
	const SimTK::Transform &reorientAB,
	SimTK::Transform& X_FMprim)
{

	SimTK::Transform X_MB = ~(ofBodyB.getOutboardFrame(someState));
	SimTK::Transform X_FM = ofBodyB.getMobilizerTransform(someState);
	SimTK::Transform X_AB = 
		ofBodyB.findBodyTransformInAnotherBody(someState, inBodyA);

	SimTK::Transform X_BBprim = (~X_AB) * reorientAB;
	X_FMprim = X_FM * X_MB * X_BBprim * X_MB;

	return X_FMprim;
}

//...............

/**
 * Set X_PF, X_BM means
*/
void World::setTransformsMeans(const std::vector<SimTK::Real>& givenX_PF,
const std::vector<SimTK::Real>& givenX_BM)
{
	// Update acosX_PF00 means
	int i = -1;
	for(auto &xpf : acosX_PF00_means ){
		i += 1;
		xpf = givenX_PF[i]; 
	}
	
	// Update normX_BMp means
	i = -1;
	for(auto &xbm : normX_BMp_means ){
		i += 1;
		xbm = givenX_BM[i];
	}

}

/**
 * Set X_PF, X_BM means to initial values
*/
void World::setTransformsMeansToIni(void)
{
	const SimTK::State& defaultState = matter->getSystem().getDefaultState();

	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const Transform& X_PF = mobod.getInboardFrame(defaultState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const Transform& X_FM = mobod.getMobilizerTransform(defaultState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const Transform& X_BM = mobod.getOutboardFrame(defaultState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		/* SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		std::cout 
			<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			<< "angleMean " << int(mbx) - 1 << " "
			<< angleMean * (180 / SimTK::Pi) << " "
			<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */

	}

}

/*
 * Shift all the generalized coordinates
 */
void World::setTransformsMeansToCurrent(SimTK::State& someState)
{
	// Get generalized coordinates Q template values. These are the values that
	// Q is extending. In the case of an AnglePin mobilizer, Q is extending an 
	// <(P_x, F_x) angle.

	// Get bonds and angles values
	for (SimTK::MobilizedBodyIndex mbx(1);
		mbx < matter->getNumBodies();
		++mbx){

		// Get mobod
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);

		// Get mobod inboard frame X_PF
		const Transform& X_PF = mobod.getInboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_PF\n" << X_PF << std::endl;

		// Get mobod inboard frame X_FM measured and expressed in P
		const Transform& X_FM = mobod.getMobilizerTransform(someState);
		//std::cout << "mobod " << mbx << " X_FM\n" << X_FM << std::endl;

		// Get mobod inboard frame X_BM
		const Transform& X_BM = mobod.getOutboardFrame(someState);
		//std::cout << "mobod " << mbx << " X_BM\n" << X_BM << std::endl;
		//std::cout << "mobod " << mbx << " X_PM\n" << X_PF * X_FM * (~X_BM) << std::endl;

		SimTK::Vec3 bondVector = X_BM.p();
		acosX_PF00_means[int(mbx) - 1] = std::acos(X_PF.R()(0)(0));
		normX_BMp_means[int(mbx) - 1] = bondVector.norm();

		// Print something for now
		/* SimTK::Real bond = normX_BMp[int(mbx) - 1];
		SimTK::Real bondMean = normX_BMp_means[int(mbx) - 1];
		SimTK::Real angle = acosX_PF00[int(mbx) - 1];
		SimTK::Real angleMean = acosX_PF00_means[int(mbx) - 1];

		std::cout 
			<< "bondMean " << int(mbx) - 1 << " " << bondMean << " "
			<< "bond " << int(mbx) - 1 << " " << bond << " "
			<< "angleMean " << int(mbx) - 1 << " "
			<< angleMean * (180 / SimTK::Pi) << " "
			<< "angle " << int(mbx) - 1 << " " << angle * (180 / SimTK::Pi) << " "
			<< std::endl; */

	}

}

/**
 * Set bonds values
*/
void World::setTransformsMeansToMin(readAmberInput &amberReader)
{
	const SimTK::State& defaultState = matter->getSystem().getDefaultState();

	// Set bonds and angles values
	for(int bondIndex = 0; bondIndex < amberReader.getNumberBonds(); bondIndex++){

		int prm_a_1 = amberReader.getBondsAtomsIndex1(bondIndex);
		int prm_a_2 = amberReader.getBondsAtomsIndex2(bondIndex);
		
		//std::cout << "setTransformsStatisticsToMin atomIxs " << a_1 << " " << a_2 << " ";

		for (auto& topology : (*topologies)){

			//bSpecificAtom * gAtom = topology.bAtomList[prm_a_1];
			bool rinClosing = topology.bonds[bondIndex].isRingClosing();

			SimTK::Compound::AtomIndex aIx_1 = topology.bAtomList[prm_a_1].getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex dAIx_1 = topology.getDuMMAtomIndex(aIx_1);
			const SimTK::MobilizedBodyIndex mbx_1 = forceField->getAtomBody(dAIx_1);
			SimTK::Compound::AtomIndex aIx_2 = topology.bAtomList[prm_a_2].getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex dAIx_2 = topology.getDuMMAtomIndex(aIx_2);
			const SimTK::MobilizedBodyIndex mbx_2 = forceField->getAtomBody(dAIx_2);

			const SimTK::MobilizedBodyIndex mbx = mbx_2;

			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			const Transform& X_PF = mobod.getInboardFrame(defaultState);
			const Transform& X_BM = mobod.getOutboardFrame(defaultState);

			std::cout << "bond redundancy " << prm_a_1 << " " << prm_a_2 << " " 
				<< int(mbx_1) << " " <<  int(mbx_2)  << std::endl;

			if(rinClosing){
				std::cout << "ring closing\n";
			}
			else if((int(mbx_2) == 1)){
				normX_BMp_means[int(mbx) - 1] = X_BM.p().norm();
				
			}else{
				normX_BMp_means[int(mbx) - 1] =
					amberReader.getBondsEqval(bondIndex) / 10.0; // Ang to nano
			}

			//if(mobod.getNumQ(defaultState) == 2){}else{}
			//std::cout << " mbx " << int(mbx) << " normX_BMMean[" << int(mbx) - 1 << "]= " << normX_BMp_means[int(mbx) - 1] << std::endl;

		}
	}	

	// Set angles values
	for(int angleIndex = 0; angleIndex < amberReader.getNumberAngles(); angleIndex++){

		int prm_a_1 = amberReader.getAnglesAtomsIndex1(angleIndex);
		int prm_a_2 = amberReader.getAnglesAtomsIndex2(angleIndex);
		int prm_a_3 = amberReader.getAnglesAtomsIndex3(angleIndex);

		for (auto& topology : (*topologies)){

			SimTK::Compound::AtomIndex aIx_1 = topology.bAtomList[prm_a_1].getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex dAIx_1 = topology.getDuMMAtomIndex(aIx_1);
			const SimTK::MobilizedBodyIndex mbx_1 = forceField->getAtomBody(dAIx_1);
			SimTK::Compound::AtomIndex aIx_2 = topology.bAtomList[prm_a_2].getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex dAIx_2 = topology.getDuMMAtomIndex(aIx_2);
			const SimTK::MobilizedBodyIndex mbx_2 = forceField->getAtomBody(dAIx_2);
			SimTK::Compound::AtomIndex aIx_3 = topology.bAtomList[prm_a_3].getCompoundAtomIndex();
			SimTK::DuMM::AtomIndex dAIx_3 = topology.getDuMMAtomIndex(aIx_3);
			const SimTK::MobilizedBodyIndex mbx_3 = forceField->getAtomBody(dAIx_3);

			const SimTK::MobilizedBodyIndex mbx = mbx_3;

			const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
			const Transform& X_PF = mobod.getInboardFrame(defaultState);
			const Transform& X_BM = mobod.getOutboardFrame(defaultState);

			std::cout << "angle redundancy " << prm_a_1 << " " << prm_a_2 << " " << prm_a_3 << " "
								<< int(mbx_1) << " " <<  int(mbx_2) << " " <<  int(mbx_3) << std::endl;

			if( (int(mbx_3) == 1) ){
				acosX_PF00_means[int(mbx) - 1] =
					std::acos(X_PF.R()(0)(0));
			}else{
				acosX_PF00_means[int(mbx) - 1] =
					amberReader.getAnglesEqval(angleIndex);
			}
				
			//if(mobod.getNumQ(defaultState) == 2){}else{}

		}
	}	


}

/**
 * Update X_PF, X_BM means
*/
void World::updateTransformsMeans(SimTK::State& someState)
{
	int nofSamples = getNofSamples() + 1;
	//std::cout << "Nof samples " << nofSamples << std::endl;

	// Useful vars
	SimTK::Real N_1overN = 9999, NInv = 9999;

	if(nofSamples == 1){
		for(unsigned int k = 0; k < acosX_PF00_means.size(); k++){
			acosX_PF00_means[k] = acosX_PF00[k];
		}
		for(unsigned int k = 0; k < normX_BMp_means.size(); k++){
			normX_BMp_means[k] = normX_BMp[k];
			//std::cout << "World " << ownWorldIndex << " bondUpdMean " << k << " " << normX_BMp_means[k] << std::endl;
		}
	}else{
		if(nofSamples == 2){
			N_1overN = NInv = 0.5;
		}else{
			// Update useful vars
			SimTK::Real N_1 = nofSamples - 1.0;
			N_1overN = N_1 / nofSamples;
			NInv = 1.0 / nofSamples;
		}
		//std::cout << "updateX_PFMeans check " << " "
		//	<<  N_1overN << " " <<  NInv  << " " << std::flush;

		// Update acosX_PF00 means
		int i = -1;
		for(auto &xpf : acosX_PF00_means ){
			i += 1;
			xpf = (N_1overN * xpf) + (NInv * acosX_PF00.at(i)); 
		}
		
		// Update normX_BMp means
		i = -1;
		for(auto &xbm : normX_BMp_means ){
			i += 1;
			xbm = (N_1overN * xbm) + (NInv * normX_BMp.at(i));
			//std::cout << "World " << ownWorldIndex << " bondUpdMean " << i << " " << xbm << std::endl;
		}
	}

}

// Get X_PF means
std::vector<SimTK::Real>& World::getX_PFMeans(void)
{
	return acosX_PF00_means;
}

// Get X_BM means
std::vector<SimTK::Real>& World::getX_BMMeans(void)
{
	return normX_BMp_means;
}

/**
 * Calculate bond length and angle deviations from their means
*/ 
void World::calcBendStretchDeviations(SimTK::State& someState,
	std::vector<SimTK::Real>& X_PFdiffs,
	std::vector<SimTK::Real>& X_BMdiffs
)
{
	// Make sure it has 
	X_PFdiffs.resize(this->acosX_PF00_means.size(), 0.0);
	X_BMdiffs.resize(this->normX_BMp_means.size(), 0.0);

	// 
	for(unsigned int k = 0; k < X_PFdiffs.size(); k++){
		X_PFdiffs[k] = this->acosX_PF00[k] - this->acosX_PF00_means[k];
	}
	for(unsigned int k = 0; k < X_BMdiffs.size(); k++){
		X_BMdiffs[k] = this->normX_BMp[k] - this->normX_BMp_means[k];
		/* std::cout << "World " << ownWorldIndex << " bondDiff " << k << " "
		//<< this->normX_BMp[k] << " " << this->normX_BMp_means[k] << " "
		<< X_BMdiffs[k] << std::endl; */
	}

}

// Get the number of molecules
int World::getNofMolecules() const
{
	return (this->moleculeCount + 1);
}

/** Get MobilizedBody to AtomIndex map **/
std::map< SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex >
World::getMbx2aIx(){
	return mbx2aIx;
}

std::size_t World::getNofMobilizedBodies() const{
	return mbx2aIx.size();
}

/** Get U scale factor for the mobilized body **/
SimTK::Real World::getMobodUScaleFactor(SimTK::MobilizedBodyIndex& mbx) const
{
	if(!mbx2uScale.empty()){
		if(mbx2uScale.find(mbx) != mbx2uScale.end()){
			return mbx2uScale.at(mbx);
		}else{
			//std::cout << "Warning: U scale factor for mobod " << int(mbx) << " not found.\n";
			return 1;
		}
	}else{
		return 1;
	}
}

/** Print atom to MobilizedBodyIndex and bond to Compound::Bond index maps **/
void World::printMaps(void)
{
	for (auto& topology : (*topologies)){
		topology.printMaps();
	}
}

SimTK::Vec3 World::calcAtomLocationInGroundFrameThroughOMM(const SimTK::DuMM::AtomIndex& )
{
	assert("NOT IMPLEMENTED!");
	return SimTK::Vec3(-1, -1, -1);
}

//==============================================================================
//                             2. Inter-world functions.
//==============================================================================
// Pass configurations between Worlds


// RANDOM_WALK functions
void World::setTopologyIXs(std::vector<int> argTopologyIXs){
	topologyIXs = argTopologyIXs;
}

void World::setAmberAtomIXs(std::vector<std::vector<int>> argAmberAtomIXs){
	amberAtomIXs = argAmberAtomIXs;
}

SimTK::Vec3 
World::getGeometricCenterOfSelection(const SimTK::State & state 
									 //const std::vector<int>& topologyIx, 
									 //const std::vector<std::vector<int>>& amberAtomList
									 )
{

	// return Vec3
	SimTK::Vec3 geometricCenter={0,0,0};
	// We could just divide by the size of amberAtomList
	// but this works even if the user *mistakenly* repeats
	// indices
	int nOfPoints=0;


	// Just a quick check, to skip unnecessary computation in case of 
	// user error.
	if (amberAtomIXs.size() == 0) {
		std::cerr << "Warning: getGeometricCenterOfSelection called with amberAtomList of size 0" << std::endl;
		return geometricCenter;
	}
	
	std::cout << "topologies atoms size " << topologyIXs.size() << " " << topologyIXs.size() << std::endl;

	for (int i = 0; i < topologyIXs.size(); i++) {
		const auto& topology = (*topologies)[topologyIXs[i]];
		const auto& atoms = amberAtomIXs[i];

		int amberIx=0;

		// Iterate through atoms in said topology and check 	
		// if they are in the list
		for (auto& atom : topology.bAtomList) {
			if (std::find(atoms.begin(), atoms.end(), amberIx) != atoms.end()){
				// found
				// Get Compound atom index
				auto compoundAtomIndex = atom.getCompoundAtomIndex();
				// Get DuMM atom index
				const SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(compoundAtomIndex);
				// Get Mobilized Body index
				const MobilizedBodyIndex mobilizedBodyIndex = forceField->getAtomBody(dAIx);
				// Get DuMM Atom Station on its body.
				const Vec3 dAS_B = forceField->getAtomStationOnBody(dAIx);
				// Re-Express in G
				const SimTK::MobilizedBody& mobod_A = matter->getMobilizedBody(mobilizedBodyIndex);
				const SimTK::Vec3 dAS_G = mobod_A.findStationLocationInGround(state, dAS_B);
				/* const Transform& X_GP = mobod_A.getBodyTransform(state);
				const SimTK::Vec3 dAS_G = X_GP*dAS_B; */
				geometricCenter += dAS_G;
				nOfPoints += 1;

/* 				std::cout << "amberIx: " << amberIx << " dAIx: " << dAIx
				<< " MobilizedBodyIndex: " << mobilizedBodyIndex 
				<< " dAS_G: " << dAS_G << " nOfPoints: " << nOfPoints
				<< std::endl; */
			}
			
			amberIx += 1;
		}
	}

	// This can probably be done better, but is it clearer?
	for(int i=0;i<3;++i)
		geometricCenter[i] = geometricCenter[i] / nOfPoints;
	std::cout << "geometricCenter : " << geometricCenter << std::endl;

	return geometricCenter;

}

/** Get the current Compound Cartesian coords.
* Return a 2D vector representing all the coordinates of this World.
 * The first dimension represents the molecules (topologies) and the second
 * dimension (inner) represents the coordinates. The second inner dimension
 * type is pair of bSpecificAtom* and a Vec3. Thus, besides coordinates, it
 * contains all the information in bSpecificAtom as well. The bottleneck here
 * is the calcAtomLocationInGroundFrame from Compound.
 **/


std::vector< std::vector<
std::pair <bSpecificAtom *, SimTK::Vec3 > > >
World::getAtomsLocationsInGround(SimTK::State & state)
{
	// Return vector
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
		returnVector;

	// Iterate through topologies
	for (auto& topology : (*topologies)){
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve(topology.bAtomList.size());

		// Iterate through atoms
		for (auto& atom : topology.bAtomList) {

			// Get Compound atom index
			auto compoundAtomIndex = atom.getCompoundAtomIndex();

			// Get DuMM atom index too for OpenMM
			const SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(compoundAtomIndex);

			SimTK::Vec3 location;

			if(samplers[0]->getIntegratorName() == IntegratorName::OMMVV){
				// ELIZA
				// location = calcAtomLocationInGroundFrameThroughOMM(dAIx);
				location = forceField->calcAtomLocationInGroundFrameThroughOMM(dAIx);
			}else{
				location = 
				topology.calcAtomLocationInGroundFrameThroughSimbody(
					compoundAtomIndex, *forceField, *matter, state);
			}

			// std::cout << "getAtomsLocationsInGround " << dAIx << " " << location << std::endl;

			currentTopologyInfo.emplace_back(&atom, location);
		}

		returnVector.emplace_back(currentTopologyInfo);
	}

	return returnVector;
}

/** Order: bAtomList
 * Almost the same thing as above but it takes the current integrator state
 * and returns a reference to atomsLocations
 **/
std::vector< std::vector< 
std::pair <bSpecificAtom *, SimTK::Vec3 > > >
World::getCurrentAtomsLocationsInGround(void)
{
	// Get the state from the integrator
	SimTK::State& state = integ->updAdvancedState();

	// Return vector
	std::vector<std::vector<std::pair <bSpecificAtom *, SimTK::Vec3>>>
		returnVector;

	// Iterate through topologies
	for (auto& topology : (*topologies)){ 
		std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>> currentTopologyInfo;
		currentTopologyInfo.reserve(topology.bAtomList.size());

		// Iterate through atoms
		for (auto& atom : topology.bAtomList) {

			// Get Compound atom index
			auto compoundAtomIndex = atom.getCompoundAtomIndex();

			// Get DuMM atom index too for OpenMM
			const SimTK::DuMM::AtomIndex dAIx = topology.getDuMMAtomIndex(compoundAtomIndex);

			// Calculate atom location in Ground
			SimTK::Vec3 location;

			// if(samplers[0]->getIntegratorName() == IntegratorName::OMMVV){
			// 	// ELIZA
			// 	location = calcAtomLocationInGroundFrameThroughOMM(dAIx);
			// }else{
				location = 
				topology.calcAtomLocationInGroundFrameThroughSimbody(
					compoundAtomIndex, *forceField, *matter, state);
			//}

			currentTopologyInfo.emplace_back(&atom, location);
		}

		returnVector.emplace_back(currentTopologyInfo);
	}

	return returnVector;
}


/** Nice print helper for get/setAtomsLocations */
void World::PrintAtomsLocations(const std::vector<std::vector<
	std::pair<bSpecificAtom *, SimTK::Vec3> > >& someAtomsLocations)
{	
	std::cout << "myAtomsLocations[0]" << std::endl;
	for(std::size_t j = 0; j < someAtomsLocations[0].size(); j++){
		int compoundAtomIndex = someAtomsLocations[0][j].first->getCompoundAtomIndex();
		auto loc = someAtomsLocations[0][j].second;

		printf("%d %.10f %.10f %.10f\n", compoundAtomIndex, loc[0], loc[1], loc[2]);
	}
}

/**
 * Write coordinates to a rst7 file
 * Very costly
*/
void World::WriteRst7FromTopology(std::string FN)

{
	updateAtomListsFromCompound(integ->updAdvancedState());

	FILE *File = fopen(FN.c_str(), "w+");

	int Natoms = 0;
	for(auto& topology : (*topologies)){
		Natoms += topology.getNAtoms();
	}

	fprintf(File, "TITLE: Created by Robosample with %d atoms\n", Natoms);
	fprintf(File, "%6d\n", Natoms);

	int atomCnt = -1;
	for(auto& topology : (*topologies)){
		for (auto& atom : topology.bAtomList) {			++atomCnt;
			fprintf(File, "%12.7f%12.7f%12.7f", 
				atom.getX() * 10.0, atom.getY() * 10.0, atom.getZ() * 10.0);

			if(atomCnt % 2 == 1){
				fprintf(File, "\n");
			}
		}
	}

	if(atomCnt % 2 == 0){
		fprintf(File, "\n");
	}

	fflush(File);
	fclose(File);
}

/** Print transformation geometries */
void World::PrintFullTransformationGeometry(const SimTK::State& someState,
		bool x_pf_r, bool x_fm_r, bool x_bm_r,
		bool x_pf_p, bool x_fm_p, bool x_bm_p)
{
	std::cout << std::fixed << std::setprecision(3);
	std::cout << "X_PF X_FM X_MB\n";
	for (SimTK::MobilizedBodyIndex mbx(1); mbx < matter->getNumBodies(); ++mbx){
		SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
		const Transform& X_PF = mobod.getDefaultInboardFrame();
		const Transform& X_FM = mobod.getMobilizerTransform(someState);
		const Transform& X_BM = mobod.getOutboardFrame(someState);
		const Transform& X_MB = ~X_BM;

		if( x_pf_r && x_fm_r && x_bm_r && x_pf_p && x_fm_p && x_bm_p){
			std::cout << "mobod " << int(mbx) << std::endl
				<< std::fixed << std::setprecision(3);

			printf("%9.3f %9.3f %9.3f %9.3f ",  X_PF.toMat44()[0][0], X_PF.toMat44()[0][1], X_PF.toMat44()[0][2], X_PF.toMat44()[0][3]);
			printf("%9.3f %9.3f %9.3f %9.3f ",  X_FM.toMat44()[0][0], X_FM.toMat44()[0][1], X_FM.toMat44()[0][2], X_FM.toMat44()[0][3]);
			printf("%9.3f %9.3f %9.3f %9.3f\n", X_MB.toMat44()[0][0], X_MB.toMat44()[0][1], X_MB.toMat44()[0][2], X_MB.toMat44()[0][3]);

			printf("%9.3f %9.3f %9.3f %9.3f ",  X_PF.toMat44()[1][0], X_PF.toMat44()[1][1], X_PF.toMat44()[1][2], X_PF.toMat44()[1][3]);
			printf("%9.3f %9.3f %9.3f %9.3f ",  X_FM.toMat44()[1][0], X_FM.toMat44()[1][1], X_FM.toMat44()[1][2], X_FM.toMat44()[1][3]);
			printf("%9.3f %9.3f %9.3f %9.3f\n", X_MB.toMat44()[1][0], X_MB.toMat44()[1][1], X_MB.toMat44()[1][2], X_MB.toMat44()[1][3]);

			printf("%9.3f %9.3f %9.3f %9.3f ",  X_PF.toMat44()[2][0], X_PF.toMat44()[2][1], X_PF.toMat44()[2][2], X_PF.toMat44()[2][3]);
			printf("%9.3f %9.3f %9.3f %9.3f ",  X_FM.toMat44()[2][0], X_FM.toMat44()[2][1], X_FM.toMat44()[2][2], X_FM.toMat44()[2][3]);
			printf("%9.3f %9.3f %9.3f %9.3f\n", X_MB.toMat44()[2][0], X_MB.toMat44()[2][1], X_MB.toMat44()[2][2], X_MB.toMat44()[2][3]);

			printf("%9.3f %9.3f %9.3f %9.3f ",  X_PF.toMat44()[3][0], X_PF.toMat44()[3][1], X_PF.toMat44()[3][2], X_PF.toMat44()[3][3]);
			printf("%9.3f %9.3f %9.3f %9.3f ",  X_FM.toMat44()[3][0], X_FM.toMat44()[3][1], X_FM.toMat44()[3][2], X_FM.toMat44()[3][3]);
			printf("%9.3f %9.3f %9.3f %9.3f\n", X_MB.toMat44()[3][0], X_MB.toMat44()[3][1], X_MB.toMat44()[3][2], X_MB.toMat44()[3][3]);
			
		}else if(x_pf_p && x_fm_p && x_bm_p){
			std::cout << X_PF.p()[0] << " " << X_PF.p()[1] << " " << X_PF.p()[2] << " ";
			std::cout << X_FM.p()[0] << " " << X_FM.p()[1] << " " << X_FM.p()[2] << " ";
			std::cout << X_BM.p()[0] << " " << X_BM.p()[1] << " " << X_BM.p()[2] << " ";
	
		}else if(x_pf_r){
			PrintMat33(X_PF.R().toMat33(), 3, "X_PF.R");
		}
		else if(x_fm_r){
			PrintMat33(X_FM.R().toMat33(), 3, "X_FM.R");
		}
		else if(x_bm_r){
			PrintMat33(X_BM.R().toMat33(), 3, "X_BM.R");
		}

		else if(x_pf_p){
			std::cout << X_PF.p()[0] << " " << X_PF.p()[1] << " " << X_PF.p()[2] << " ";
		}
		else if(x_fm_p){
			std::cout << X_FM.p()[0] << " " << X_FM.p()[1] << " " << X_FM.p()[2] << " ";
		}
		else if(x_bm_p){
			std::cout << X_BM.p()[0] << " " << X_BM.p()[1] << " " << X_BM.p()[2] << " ";
		}

		std::cout << std::endl;

	}
}



/** Put coordinates into bAtomLists of Topologies.
 * When provided with a State, calcAtomLocationInGroundFrame
 * realizes Position and uses matter to calculate locations **/
void World::updateAtomListsFromCompound(const SimTK::State &state)
{
	// Iterate through topologies
	for (auto& topology : (*topologies)){

		// Iterate through atoms
		for (auto& atom : topology.bAtomList) {

			const auto compoundAtomIndex = atom.getCompoundAtomIndex();
			SimTK::Vec3 location =
				topology.calcAtomLocationInGroundFrameThroughSimbody(
					compoundAtomIndex, *forceField, *matter, state);

			atom.setX(location[0]);
			atom.setY(location[1]);
			atom.setZ(location[2]);
			atom.setCartesians(location);

			// std::cout << "updateAtomListsFromCompound (after f_x_m, ix= " << compoundAtomIndex << ") " << atom.getX() << ", " << atom.getY() << ", " << atom.getZ() << std::endl;
		}
	}
}


/**
 * RMSD function
*/
SimTK::Real World::RMSD(
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations
		) const
{

	assert( srcWorldsAtomsLocations.size() == destWorldsAtomsLocations.size() && 
		!("RMSD different size") );

	SimTK::Real rmsdVal = 0;
	int localNofAtoms = 0;

	for(std::size_t i = 0; i < srcWorldsAtomsLocations.size(); i++){
		for(std::size_t j = 0; j < srcWorldsAtomsLocations[i].size(); j++){

			auto srcLoc = srcWorldsAtomsLocations[i][j].second;
			auto destLoc = destWorldsAtomsLocations[i][j].second;

			SimTK::Real sqX = (srcLoc[0] - destLoc[0]);
			sqX = sqX * sqX;
			SimTK::Real sqY = (srcLoc[1] - destLoc[1]);
			sqY = sqY * sqY;
			SimTK::Real sqZ = (srcLoc[2] - destLoc[2]);
			sqZ = sqZ * sqZ;

			rmsdVal += (sqX + sqY + sqZ);

			localNofAtoms++;

		}
	}

	rmsdVal /= SimTK::Real(localNofAtoms);
	rmsdVal = std::sqrt(rmsdVal);

	return rmsdVal;

}

/**
 * Maximum distance between two corresponding atoms
*/
std::pair<int, SimTK::Real> World::maxAtomDeviation(
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		 srcWorldsAtomsLocations,
	const std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3> > >&
		destWorldsAtomsLocations
		) const
{

	assert( srcWorldsAtomsLocations.size() == destWorldsAtomsLocations.size() && 
		!("RMSD different size") );

	SimTK::Real distVal = 0;
	SimTK::Real maxVal = 0;
	int atomCnt = -1;
	int pairIndex = -1;

	for(std::size_t i = 0; i < srcWorldsAtomsLocations.size(); i++){
		for(std::size_t j = 0; j < srcWorldsAtomsLocations[i].size(); j++){

			atomCnt++;

			auto srcLoc = srcWorldsAtomsLocations[i][j].second;
			auto destLoc = destWorldsAtomsLocations[i][j].second;

			SimTK::Real sqX = (srcLoc[0] - destLoc[0]);
			sqX = sqX * sqX;
			SimTK::Real sqY = (srcLoc[1] - destLoc[1]);
			sqY = sqY * sqY;
			SimTK::Real sqZ = (srcLoc[2] - destLoc[2]);
			sqZ = sqZ * sqZ;

			SimTK::Real dist = std::sqrt(sqX + sqY + sqZ);
			
			if(dist > distVal){
				distVal = dist;
				pairIndex = atomCnt;
			}

		}
	}

	return std::pair<int, SimTK::Real>(pairIndex, distVal);

}


/** 
 * Set Compound, MultibodySystem and DuMM configurations according to
 * some other World's atoms.
 * A body is composed of a root atom and other periferic atoms which have their own stations.
 * Unles is a fully flexible Cartesian world, the function has the following steps:
 * 1. Set Compound
 * 2. Set DuMM
 * 3. Set Simbody bodies
 *  3.1. Transforms X_PF and X_BM
 *  3.2. Mass properties
 * 
**/
SimTK::State& World::setAtomsLocationsInGround(
		SimTK::State& someState,
		const std::vector<std::vector<
		std::pair<bSpecificAtom *, SimTK::Vec3> > >& otherWorldsAtomsLocations)
{

	// ========================================================================
	// Check reconstruction
	// ========================================================================
	SimTK::Real pe_init = SimTK::NaN;
	bool reconstructFail = false;

	// Try to reconstruct multiple times
	int nofReconstructAttempts = 1;
	for(int reCnt = 0; reCnt < nofReconstructAttempts; reCnt++){

	pe_init = SimTK::NaN;
	reconstructFail = false;

	if(CHECK_RECONSTRUCTION){
		pe_init = forces->getMultibodySystem().calcPotentialEnergy(someState);
	}
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	// Get the total no of bodies in this world (each World has its own
	// SimbodyMatterSubsystem)
	int totalNofBodies = matter->getNumBodies();

	// Arrays of Transforms
	SimTK::Transform G_X_T;
	//SimTK::Transform T_X_root[totalNofBodies];

	// Loop through molecules/topologies
	for(std::size_t topoIx = 0; topoIx < otherWorldsAtomsLocations.size(); topoIx++) {

					// Convenient vars
					Topology& currTopology = (*topologies)[topoIx];
					int currNAtoms = currTopology.getNumAtoms();

					// Set the decorator
					if (visual == true) {
						paraMolecularDecorator->setAtomTargets(otherWorldsAtomsLocations[topoIx]);
					}

					///////////////////////////////////////////////////////////
					// 1. COMPOUND
					///////////////////////////////////////////////////////////

					// Use Molmodel's Compound match functions to set the new
					// Compound transforms
					std::map<SimTK::Compound::AtomIndex, SimTK::Vec3> atomTargets;
					for(std::size_t j = 0; j < otherWorldsAtomsLocations[topoIx].size(); j++){
						auto atomIndex = otherWorldsAtomsLocations[topoIx][j].first->getCompoundAtomIndex();
						auto location = otherWorldsAtomsLocations[topoIx][j].second;
						atomTargets.insert(std::make_pair(atomIndex, location));
					}

					//std::cout << "Match start." << "\n" << std::flush;
					bool flipAllChirality = false;
					currTopology.matchDefaultAtomChirality(atomTargets, 0.01, flipAllChirality);
					//std::cout << "matchDefaultAtomChirality done. " << "\n" << std::flush;
					currTopology.matchDefaultBondLengths(atomTargets);
					//std::cout << "matchDefaultBondLengths done. " << "\n" << std::flush;
					currTopology.matchDefaultBondAngles(atomTargets);
					//std::cout << "matchDefaultBondAngles done. " << "\n" << std::flush;
					currTopology.matchDefaultDirections(atomTargets);
					//std::cout << "matchDefaultDirections done. " << "\n" << std::flush;
					currTopology.matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
					//std::cout << "matchDefaultDefaultDihedralAngles done. " << "\n" << std::flush;
					currTopology.matchDefaultTopLevelTransform(atomTargets);
					//std::cout << "matchDefaultDefaultTopLevelTransform done. " << "\n" << std::flush;

					// Get the Ground to Top Transform
					G_X_T = currTopology.getTopLevelTransform();

					// Recalculate atom frames in top compound frame
					//std::cout << "Calculate defaultAtomFrames start ...." << "\n" << std::flush;
					currTopology.calcTopTransforms();
					//std::cout << "defaultAtomFrames done" << "\n" << std::flush;

					///////////////////////////////////////////////////////////
					// 2.1 MORE COMPOUND FOR DUMM
					///////////////////////////////////////////////////////////
					// Get locations for DuMM
					SimTK::Vec3 locs[currTopology.getNumAtoms()];

					// Loop through atoms
					for (SimTK::Compound::AtomIndex aIx(1); aIx < currTopology.getNumAtoms(); ++aIx){
						if(currTopology.getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, *forceField) == 0){ // atom is at body's origin
								locs[int(aIx)] = SimTK::Vec3(0);
						}else{ // atom is not at body's origin
							SimTK::MobilizedBodyIndex mbx = currTopology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
							SimTK::Transform G_X_root = G_X_T * currTopology.getTopTransform(mbx2aIx[mbx]); // DANGER

							SimTK::Vec3 G_vchild = atomTargets[aIx];

							SimTK::Transform root_X_child = alignFlipAndTranslateFrameAlongXAxis(G_X_root, G_vchild);
							currTopology.bsetFrameInMobilizedBodyFrame(aIx, root_X_child);

							locs[int(aIx)] = root_X_child.p();
						}
					}

					//std::cout << "MORE COMPOUND FOR DUMM done" << "\n" << std::flush;

					///////////////////////////////////////////////////////////
					// 2.2 DUMM
					///////////////////////////////////////////////////////////
					// Set stations and AtomPLacements for atoms in DuMM
					// Loop through atoms
					for (SimTK::Compound::AtomIndex aIx(1); aIx < currTopology.getNumAtoms(); ++aIx){
						SimTK::MobilizedBodyIndex mbx =
							currTopology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);

						//std::cout << "mbx " << mbx << "\n" << std::flush;

						SimTK::DuMM::AtomIndex dAIx = currTopology.getDuMMAtomIndex(aIx);

						//std::cout << "dAIx " << dAIx << "\n" << std::flush;

						// Set station_B
						forceField->bsetAtomStationOnBody( dAIx, locs[int(aIx)] );
						forceField->bsetAllAtomStationOnBody( dAIx, locs[int(aIx)] ); // full

						//std::cout << "bset " << locs[int(aIx)] << "\n" << std::flush;
						//std::cout << "dumm nofatoms " << forceField->getNumAtoms() << "\n" << std::flush;

						// Set included atom
						forceField->updIncludedAtomStation(dAIx) = (locs[int(aIx)]);
						//std::cout << "updIncludedAtomStation " << "\n" << std::flush;
						forceField->updAllAtomStation(dAIx) = (locs[int(aIx)]); // full
						//std::cout << "updAllAtomStation " << "\n" << std::flush;


						// Atom placements in clusters
						forceField->bsetAtomPlacementStation(dAIx, mbx, locs[int(aIx)] );

						//std::cout << "bsetAtomPl " << "\n" << std::flush;

					} /////////////////////////

					//std::cout << "DUMM done" << "\n" << std::flush;

					///////////////////////////////////////////////////////////
					// 3. SIMBODY MATTER
					//---------------------------------------------------------
					// 3.1 MOBILIZED BODIES TRANSFORMS
					///////////////////////////////////////////////////////////

					// Set X_PF and X_BM in Mobilized bodies. First get the atom 0
					// transform. Works for multiple mols because every mol has Ground
					// as parent
					SimTK::Transform P_X_F_1 = G_X_T * 
						currTopology.getTopTransform(SimTK::Compound::AtomIndex(0)); // is this always true ??

					// Loop through the rest of the atoms and get P_X_F, B_X_M from
					// the Compound transforms
					for (SimTK::Compound::AtomIndex aIx(0);
					aIx < currTopology.getNumAtoms(); ++aIx)
					{
						// Is atom at the body's origin
						if(currTopology.getAtomLocationInMobilizedBodyFrameThroughDumm(
						aIx, *forceField) == 0)
						{

							// Get body, parentBody, parentAtom
							SimTK::MobilizedBodyIndex mbx =
								currTopology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
							SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
							const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
							SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

							if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground

								// Get mobod transforms
								std::vector<SimTK::Transform> mobodTs = 
									calcMobodToMobodTransforms( 
										currTopology,
										aIx,
										someState);

								mobod.setDefaultInboardFrame(mobodTs[0]);
								mobod.setDefaultOutboardFrame(mobodTs[1]);

								// std::cout << "P_X_F " << mobodTs[0] << std::endl;
								// std::cout << "B_X_M " << mobodTs[1] << std::endl;

							} // END if parent not Ground
							else{ // parent is Ground
								mobod.setDefaultInboardFrame(P_X_F_1);
								mobod.setDefaultOutboardFrame(Transform());
							}
						} // END atom is at body's origin
					} // END loop through atoms

					///////////////////////////////////////////////////////////
					// 3. SIMBODY MATTER
					//---------------------------------------------------------
					// 3.2 MASS PROPERTIES
					///////////////////////////////////////////////////////////
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

		} // END all topologies (topoIx)

		// Realize Position 
		//this->compoundSystem->realize(someState, SimTK::Stage::Position);

	///////////////////////////////////////////////////////////
	// 4. SIMBODY MATTER
	//---------------------------------------------------------
	// 4.1 Qs and X_FMs
	///////////////////////////////////////////////////////////
	// After realizePosition try to set the Qs (X_FMs) too

	// Loop through molecules/topologies
	for(std::size_t topoIx = 0; topoIx < otherWorldsAtomsLocations.size(); topoIx++) {
	
					// Convenient vars
					Topology& currTopology = (*topologies)[topoIx];
					int currNAtoms = currTopology.getNumAtoms();
	
					// Top transform: is this always true ??
					SimTK::Transform P_X_F_1 = G_X_T * 
						currTopology.getTopTransform(SimTK::Compound::AtomIndex(0));
	
					// Loop through the rest of the atoms and get F_X_M from
					// the Compound transforms
					for (SimTK::Compound::AtomIndex aIx(0); aIx < currNAtoms; ++aIx)
					{
						// Is atom at the body's origin
						if(currTopology.getAtomLocationInMobilizedBodyFrameThroughDumm(
						aIx, *forceField) == 0)
						{
	
							// Get body, parentBody, parentAtom
							SimTK::MobilizedBodyIndex mbx =
								currTopology.getAtomMobilizedBodyIndexThroughDumm(aIx, *forceField);
							SimTK::MobilizedBody& mobod = matter->updMobilizedBody(mbx);
							const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
							SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();
	
							if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
	

								SimTK::Transform X_FM = calcX_FMTransforms( 
										currTopology, aIx, someState);

								// Get mobod transforms	
								mobod.setQToFitTransform(someState, X_FM);
								//std::cout << " F_X_M " << int(aIx) << " "  <<  int(mbx) << " " << X_FM << std::endl;
	
							}else{ // parent is Ground
	
								mobod.setQToFitTransform(someState, Transform());
	
							}
						} // END atom is at body's origin
					} // END loop through atoms
	} // END loop through topologies

	// Realize Position 
	this->compoundSystem->realize(someState, SimTK::Stage::Position);

	// PRINT
	//if(ownWorldIndex == 1){
		// PrintFullTransformationGeometry(someState);
		// std::cout << "Sleeping... " << std::flush;
		// std::this_thread::sleep_for(std::chrono::milliseconds(5000));
		// std::cout << "done.\n" << std::flush;
	//}

	// Copy everything to OpenMM too
	if(getSampler(0)->getIntegratorName() == IntegratorName::OMMVV){
		updSampler(0)->Simbody_To_OMM_setAtomsLocationsCartesian(someState); // COMPLETE
	}

	// ========================================================================
	// Check reconstruction
	// ========================================================================
	if(CHECK_RECONSTRUCTION){
		SimTK::Real deltaPE =
			forces->getMultibodySystem().calcPotentialEnergy(someState) -
			pe_init;

			if(!this->updSampler(0)->checkDistortionBasedOnE(deltaPE)){
				std::cout << "deltaPE " <<  deltaPE
					<< " from setAtomsLocations check distortion: ";

				std::vector<std::vector<std::pair<bSpecificAtom *, SimTK::Vec3>>>
					cal = getCurrentAtomsLocationsInGround();

				SimTK::Real rmsd =
					RMSD(otherWorldsAtomsLocations,
					getCurrentAtomsLocationsInGround());

				std::pair<int, SimTK::Real> mad =
					maxAtomDeviation(otherWorldsAtomsLocations,
					getCurrentAtomsLocationsInGround());

				std::cout << "[WARNING] RMSD (nm) " << rmsd << std::endl;
				std::cout << "[WARNING] MAD (nm) " << mad.first << " " << mad.second << std::endl;
				PrintFullTransformationGeometry(someState, false, false, false, true, true, true);
				/* for(std::size_t i = 0; i < otherWorldsAtomsLocations.size(); i++){
					std::cout << "otherWorldsAtomsLocations[" << i << "]" << std::endl;
					for(std::size_t j = 0; j < otherWorldsAtomsLocations[i].size(); j++){
						auto compoundAtomIndex = otherWorldsAtomsLocations[i][j].first->getCompoundAtomIndex();
						auto loc = otherWorldsAtomsLocations[i][j].second;
						auto calloc = cal[i][j].second;
						printf("%d %.10f %.10f %.10f %.10f %.10f %.10f\n", int(compoundAtomIndex),
							loc[0], loc[1], loc[2], calloc[0], calloc[1], calloc[2]);
					}
				}
				std::cout << "=================" << std::endl << std::flush; */

				reconstructFail = true;

			}else{
				reconstructFail = false;
				break; // goto return
			}
	}
	// ------------------------------------------------------------------------
	// ------------------------------------------------------------------------

	} // multiple reconstruction tries

	return someState;

}

/**
 * This function is only intended for root atoms
*/
std::vector<SimTK::Transform>
World::calcMobodToMobodTransforms(
	Topology& topology,
	SimTK::Compound::AtomIndex rootAIx,
	const SimTK::State& someState)
{

	// Get body and parentBody
	// don't know why it works
	//SimTK::MobilizedBodyIndex mbx = topology.getAtomMobilizedBodyIndex(rootAIx);
 	// this should be the correct version
	SimTK::MobilizedBodyIndex mbx =
		topology.getAtomMobilizedBodyIndexThroughDumm(rootAIx, *forceField);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Get the neighbor atom in the parent mobilized body
	SimTK::Compound::AtomIndex chemParentAIx =
		topology.getChemicalParent_IfIAmRoot(matter.get(), rootAIx, *forceField);


	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC =
	  topology.getDefaultBondCenterFrameInOtherBondCenterFrame(
		rootAIx, chemParentAIx);  

	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform(rootAIx);

	// Get Top to parent frame
	SimTK::Compound::AtomIndex parentRootAIx = mbx2aIx[parentMbx];
	// Origin of the parent mobod
	SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;
	
	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle =
		topology.bgetDefaultInboardDihedralAngle(rootAIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get inboard bond length
	SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(rootAIx);
	SimTK::Transform InboardLength_mZAxis
		= SimTK::Transform(Rotation(), Vec3(0, 0, -inboardBondlength));

	// Samuel Flores' terminology
	SimTK::Transform M_X_pin =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

	// Get the old PxFxMxB transform
	SimTK::Transform oldX_PB =
		Proot_X_T * T_X_root
		* InboardDihedral_XAxis * X_to_Z
		* InboardDihedral_ZAxis * Z_to_X;

	// B_X_Ms
	SimTK::Transform B_X_M = X_to_Z; // aka M_X_pin
	SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
	SimTK::Transform B_X_M_pin 		= X_parentBC_childBC * X_to_Z;
	SimTK::Transform B_X_M_univ 	= X_parentBC_childBC * Y_to_Z;

	// P_X_Fs = old P_X_B * B_X_M
	SimTK::Transform P_X_F 			= oldX_PB * B_X_M;
	SimTK::Transform P_X_F_anglePin = oldX_PB * B_X_M_anglePin;
	SimTK::Transform P_X_F_pin 		= oldX_PB * B_X_M_pin;
	SimTK::Transform P_X_F_univ 	= oldX_PB * B_X_M;

	//SimTK::Transform B_X_M_spheric  = X_parentBC_childBC * X_to_Z * InboardLength_mZAxis;
	//SimTK::Transform B_X_M_spheric  = Transform();
	//SimTK::Transform B_X_M_spheric = B_X_M_pin;
	SimTK::Transform B_X_M_spheric  = X_parentBC_childBC * X_to_Z;
	
	//SimTK::Transform P_X_F_spheric  = oldX_PB * B_X_M_pin;
	//SimTK::Transform P_X_F_spheric = Transform();
	//SimTK::Transform P_X_F_spheric = P_X_F_pin;
	SimTK::Transform P_X_F_spheric = oldX_PB * B_X_M_pin;
	
	// Get mobility (joint type)
	bSpecificAtom *atom = topology.updAtomByAtomIx(rootAIx);
	SimTK::BondMobility::Mobility mobility;
	bBond bond = 	topology.getBond(topology.getNumber(rootAIx),
					topology.getNumber(chemParentAIx));
	mobility = bond.getBondMobility(ownWorldIndex);

	bool anglePin_OR = (
		   (mobility == SimTK::BondMobility::Mobility::AnglePin)
		|| (mobility == SimTK::BondMobility::Mobility::Slider)
		|| (mobility == SimTK::BondMobility::Mobility::BendStretch)
	);

	if( (anglePin_OR) && ((atom->neighbors).size() == 1)){
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};

	}else if( (anglePin_OR) && ((atom->neighbors).size() != 1)){
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};

	}else if(	(mobility == SimTK::BondMobility::Mobility::Torsion)
			||	(mobility == SimTK::BondMobility::Mobility::Cylinder)
	){
		return std::vector<SimTK::Transform> {P_X_F_pin, B_X_M_pin};

	}else if((mobility == SimTK::BondMobility::Mobility::BallM)
	|| (mobility == SimTK::BondMobility::Mobility::Rigid)
	|| (mobility == SimTK::BondMobility::Mobility::Translation) // Cartesian
	){
		return std::vector<SimTK::Transform> {P_X_F, B_X_M};

	// Spherical
	}else if((mobility == SimTK::BondMobility::Mobility::Spherical)
	){
		return std::vector<SimTK::Transform> {
			P_X_F_spheric,
			B_X_M_spheric
		};

	}else{
		std::cout << "Warning: unknown mobility\n";
		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};
	}

	//// The Molmodel notation
	// M0 and Mr are actually used for F
	//SimTK::Transform root_X_M0 = InboardDihedral_XAxis; // name used in Molmodel
	//SimTK::Transform T_X_M0 = T_X_root[int(mbx)] * root_X_M0;
	//SimTK::Transform Proot_X_M0 = Proot_X_T * T_X_M0;
	//Transform oldX_PF = Proot_X_M0 * M_X_pin;
	//Transform oldX_BM = M_X_pin;
	//Transform oldX_MB = ~oldX_BM;
	//Transform oldX_FM = InboardDihedral_ZAxis;
	//Transform oldX_PB = oldX_PF * oldX_FM * oldX_MB;
}

/**
 * Calc X_FM transforms for reconstruction for root atoms
*/
SimTK::Transform World::calcX_FMTransforms(
	Topology& topology,
	SimTK::Compound::AtomIndex rootAIx,
	const SimTK::State& someState)
{

	// X_FM return value
	SimTK::Transform X_FM;


	// Get body and parentBody
	SimTK::MobilizedBodyIndex mbx =
		topology.getAtomMobilizedBodyIndexThroughDumm(rootAIx, *forceField);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Get the neighbor atom in the parent mobilized body
	SimTK::Compound::AtomIndex chemParentAIx =
		topology.getChemicalParent_IfIAmRoot(matter.get(), rootAIx, *forceField);

	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC =
	  topology.getDefaultBondCenterFrameInOtherBondCenterFrame(
		rootAIx, chemParentAIx);  

	// Get Top frame
	SimTK::Transform T_X_root = topology.getTopTransform(rootAIx);

	// Get Top to parent frame
	SimTK::Compound::AtomIndex parentRootAIx = mbx2aIx[parentMbx];
	SimTK::Transform T_X_Proot = topology.getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;



	// chemical parent atom
	SimTK::Transform T_X_chemProot = topology.getTopTransform(chemParentAIx);
	

	// BEGIN GET ANGLE
	SimTK::Compound::AtomIndex chemGrandParentIx;
	SimTK::Transform T_X_grand;

	if(chemParentAIx > 0){

		chemGrandParentIx = topology.getInboardAtomIndex(chemParentAIx);

		T_X_grand = topology.getTopTransform(chemGrandParentIx);

		SimTK::Vec3 v1 = (~(T_X_root.R())) * T_X_grand.p();
		SimTK::Vec3 v2 = (~(T_X_root.R())) * T_X_chemProot.p();
		SimTK::Vec3 v3 = (~(T_X_root.R()))  * T_X_root.p();

		SimTK::Real bondAngle = bAngle(v2, v1, v3);

		//if(ownWorldIndex == 1){
			std::cout << "World::calcMobodToMobodTransforms chemGrandParentIx chemParentAIx rootAIx angle "
				<< chemGrandParentIx << " " << chemParentAIx << " " << rootAIx << " " 
				//<< std::endl << T_X_grand << T_X_chemProot << T_X_root << std::endl
				//<< std::endl << topology.getTopTransform(chemGrandParentIx) << topology.getTopTransform(chemParentAIx) << topology.getTopTransform(rootAIx) << std::endl
				//<< "============================="
				//<< std::endl << topology.getTopTransform(Compound::AtomIndex(2)) << topology.getTopTransform(Compound::AtomIndex(4)) << topology.getTopTransform(Compound::AtomIndex(7)) << std::endl
				//<< v1 << " " << v2 << " " << v3 << " "
				<< bondAngle << std::endl;
		//}

	}
	// END GET ANGLE	





	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle =
		topology.bgetDefaultInboardDihedralAngle(rootAIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get inboard bond length
	SimTK::Real inboardBondlength = topology.bgetDefaultInboardBondLength(rootAIx);
	SimTK::Transform InboardLength_mZAxis
		= SimTK::Transform(Rotation(), Vec3(0, 0, -inboardBondlength));

	// Samuel Flores' terminology
	SimTK::Transform M_X_pin =
		SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis);

	// Get the old PxFxMxB transform
	SimTK::Transform oldX_PB =
		Proot_X_T * T_X_root
		* InboardDihedral_XAxis * X_to_Z
		* InboardDihedral_ZAxis * Z_to_X;

	// B_X_Ms
	SimTK::Transform B_X_M = X_to_Z; // aka M_X_pin
	SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
	SimTK::Transform B_X_M_pin 		= X_parentBC_childBC * X_to_Z;
	SimTK::Transform B_X_M_univ 	= X_parentBC_childBC * Y_to_Z;

	// P_X_Fs = old P_X_B * B_X_M
	SimTK::Transform P_X_F 			= oldX_PB * B_X_M;
	SimTK::Transform P_X_F_anglePin = oldX_PB * B_X_M_anglePin;
	SimTK::Transform P_X_F_pin 		= oldX_PB * B_X_M_pin;
	SimTK::Transform P_X_F_univ 	= oldX_PB * B_X_M;

	// Get mobility (joint type)
	bSpecificAtom *atom = topology.updAtomByAtomIx(rootAIx);
	SimTK::BondMobility::Mobility mobility;
	bBond bond = 	topology.getBond(topology.getNumber(rootAIx),
					topology.getNumber(chemParentAIx));
	mobility = bond.getBondMobility(ownWorldIndex);

	// Convenient bool
	bool anglePin_OR = (
		   (mobility == SimTK::BondMobility::Mobility::AnglePin)
		|| (mobility == SimTK::BondMobility::Mobility::Slider)
		|| (mobility == SimTK::BondMobility::Mobility::BendStretch)
	);

	// Set X_FM value
	if( (anglePin_OR) && ((atom->neighbors).size() == 1)){

		X_FM = Transform();
		return X_FM;

	}else if( (anglePin_OR) && ((atom->neighbors).size() != 1)){

		X_FM = Transform();
		return X_FM;

	}else if(	(mobility == SimTK::BondMobility::Mobility::Torsion)
			||	(mobility == SimTK::BondMobility::Mobility::Cylinder)
	){

		X_FM = Transform();
		return X_FM;

	}else if((mobility == SimTK::BondMobility::Mobility::BallM)
	|| (mobility == SimTK::BondMobility::Mobility::Rigid)
	|| (mobility == SimTK::BondMobility::Mobility::Translation)
	){
		
		X_FM = Transform();
		return X_FM;

	}else if((mobility == SimTK::BondMobility::Mobility::Spherical)
	){

		X_FM = Transform();
		//X_FM = InboardLength_mZAxis;
		return X_FM;

	}else{

		X_FM = Transform();
		std::cout << "Warning: unknown mobility\n";
		return X_FM;
	}



}










/** Set up Fixman torque **/
void World::addFixmanTorque()
{
	// Set flag
	assert(!isUsingFixmanTorque());
	useFixmanTorque = true;

	// Alloc memory for FixmanTorque implementation and add to forces
	FixmanTorqueImpl = new FixmanTorque(matter.get());
	FixmanTorqueForce = std::make_unique<Force::Custom>(*forces, FixmanTorqueImpl);

	FixmanTorqueExtImpl = new FixmanTorqueExt(matter.get());
	FixmanTorqueExtForce = std::make_unique<Force::Custom>(*forces, FixmanTorqueExtImpl);

	// FixmanTorqueImpl = new FixmanTorque(matter.get());			//	
	// FixmanTorqueForce = new Force::Custom(*forces, FixmanTorqueImpl);			//
	// FixmanTorqueExtImpl = new FixmanTorqueExt(matter.get());	//
	// FixmanTorqueExtForce = new Force::Custom(*forces, FixmanTorqueExtImpl);		//

	// for (int i = 0; i < 10; i++) {
	// 	controller.push_back(std::make_unique<SimTK::ConformationalController>(*forces, *matter, SimTK::MobilizedBodyIndex(i), SimTK::Vec3(0,0,0)));
	// 	controlForce.push_back(std::make_unique<SimTK::Force::Custom>(*forces, controller.back().get()));
	// }
}

/** Check if the Fixman torque flag is set **/
bool World::isUsingFixmanTorque() const
{
	return useFixmanTorque;
}

/** Get writble pointer to FixmanTorque implementation **/
FixmanTorque * World::updFixmanTorque()
{
	assert(isUsingFixmanTorque());
	// return FixmanTorqueImpl.get();
	return FixmanTorqueImpl;
}

/** Get pointer to FixmanTorque implementation **/
FixmanTorque * World::getFixmanTorque() const
{
	assert(isUsingFixmanTorque());
	// return FixmanTorqueImpl.get();
	return FixmanTorqueImpl;
}

// ----------------------
// --- Thermodynamics ---
// ----------------------

/** Get the World temperature **/
SimTK::Real World::getTemperature()
{
	return this->temperature;
}

/** Set this World temperature but also ths samplers and
Fixman torque temperature. **/
void World::setTemperature(SimTK::Real argTemperature)
{
	// Set the temperature for the samplers
	// TODO unused
	for (auto& sampler: samplers) {
		sampler->setTemperature(argTemperature);
	}

	// Set the temperature for this World
	this->temperature = argTemperature;

	// Set the temperature for the Fixman torque also
	if(useFixmanTorque){
		FixmanTorqueImpl->setTemperature(this->temperature);
		FixmanTorqueExtImpl->setTemperature(this->temperature);
	}
}
//...............

/** Set this World temperature but also ths samplers and
Fixman torque temperature. **/
void World::setBoostTemperature(SimTK::Real argTemperature)
{
	// Set the temperature for the samplers
	for (auto& sampler: samplers) {
		sampler->setBoostTemperature(argTemperature);
	}
}
//...............

//...................
// --- Simulation ---
//...................

/** Get/Set seed for reproducibility. **/
void World::setSeed(uint32_t argSeed)
{
	randomEngine = buildRandom32(argSeed);
	forceField->setOpenMMseed(randomEngine());
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


//	std::cout << "GLOBAL SCALE FACTORS SET TO 0.0\n";
//	forceField->setBondTorsionGlobalScaleFactor(0);
//	forceField->setAmberImproperTorsionGlobalScaleFactor(0);
//
//	forceField->setVdw12ScaleFactor(0);
//	forceField->setVdw13ScaleFactor(0);
//	forceField->setVdw14ScaleFactor(0);
//	forceField->setVdw15ScaleFactor(0);
//	forceField->setVdwGlobalScaleFactor(0);
//
//	forceField->setCoulomb12ScaleFactor(0);
//	forceField->setCoulomb13ScaleFactor(0);
//	forceField->setCoulomb14ScaleFactor(0);
//	forceField->setCoulomb15ScaleFactor(0);
//	forceField->setCoulombGlobalScaleFactor(0);

}

/** Set GBSA implicit solvent scale factor. **/
void World::setGbsaGlobalScaleFactor(SimTK::Real scaleFactor)
{
	forceField->setGbsaGlobalScaleFactor(scaleFactor);
}

/** Get a writeble pointer to the DuMM force field **/
SimTK::DuMMForceFieldSubsystem& World::getForceField()
{
	return *forceField;
}

/** Get a writeble pointer to the DuMM force field **/
SimTK::DuMMForceFieldSubsystem * World::updForceField()
{
	return forceField.get();
}

//...................
// --- Statistics ---
//...................

/** How many samples do we have so far **/
std::size_t World::getNofSamples() const
{
	// Zero it every time the user asks
	std::size_t nofSamples = 0;

	// Gather samples from all the samplers
	for(size_t i = 0; i < samplers.size(); i++){
		nofSamples += (samplers[i])->getNofSamples();
	}

	return nofSamples;
}

/** How many Samplers does this World have. **/
std::size_t World::getNofSamplers() const
{
	return samplers.size();
}

/*!
 * <!-- Add a sampler to this World using the specialized struct
for samplers names. -->
*/
bool World::addSampler(SamplerName samplerName,
	const std::string& generatorName,
	const std::string& integratorName,
	const std::string& thermostatName,
	SimTK::Real timestep,
	int mdStepsPerSample,
	int mdStepsPerSampleStd, 
	bool useFixmanPotential)
{
	// We only use HMCSampler for now
	if(samplerName == SamplerName::HMC) {

		// Construct a new sampler
		samplers.emplace_back(std::make_unique<HMCSampler>(*this, *compoundSystem, *matter, *topologies, *forceField, *forces, *ts));
		
		// Set sampler parameters
		samplers.back()->setSampleGenerator(generatorName);
		samplers.back()->setIntegratorName(integratorName);
		samplers.back()->setThermostat(thermostatName);
		samplers.back()->setTemperature(this->temperature); // TODO where???
		samplers.back()->setTimestep(timestep); // TODO should error when negative
		samplers.back()->setMDStepsPerSample(mdStepsPerSample);
		samplers.back()->setMDStepsPerSampleStd(mdStepsPerSampleStd);
		samplers.back()->setSeed(randomEngine);

		// TODO should this be inherited from parent world?
		if (useFixmanPotential) {
			samplers.back()->useFixmanPotential();
		}
	} else {
		std::cerr << "Unknown sampler name" << std::endl;
		return false;
	}

	// Initialize the sampler
	// This does not care about passed parameters
	SimTK::State& worldAdvancedState = integ->updAdvancedState();
	samplers.back()->initialize(worldAdvancedState);

	return true;
}

void World::useOpenMM(bool ommvv, SimTK::Real boostTemp, SimTK::Real timestep) {
	forceField->setUseOpenMMAcceleration(true);

	if (ommvv) {
		forceField->setUseOpenMMIntegration(true);
		forceField->setUseOpenMMCalcOnlyNonBonded(false);
		forceField->setDuMMTemperature(boostTemp);
		forceField->setOpenMMstepsize(timestep);
	} else {
		forceField->setUseOpenMMCalcOnlyNonBonded(true);
	}
}

// Get a sampler based on its position in the samplers vector
// TODO Use ampler polymorphism (was const BaseSampler *)
BaseSampler * World::getSampler(std::size_t which) const
{
	return samplers[which].get();
}

// Get a writable sampler based on its position in the samplers vector
// TODO Use Sampler polymorphism
BaseSampler * World::updSampler(std::size_t which)
{
	return samplers[which].get();
}

/** Get a const reference to a molecule **/
const Topology& World::getTopology(std::size_t moleculeNumber) const{
	return (*topologies)[moleculeNumber];
}

/** Get a writble reference to the last molecule. **/
Topology& World::updTopology(std::size_t moleculeNumber){
	//return *topologies.back();
	return (*topologies)[moleculeNumber];
}

// DOESN'T WORK WITH OPENMM
SimTK::Real World::CalcFullPotentialEnergyIncludingRigidBodies(void)
{
	SimTK::State& currentAdvancedState = integ->updAdvancedState();
	updateAtomListsFromCompound(currentAdvancedState);

	// Set old potential energy of the new world via DuMM !!!
	return forceField->CalcFullPotEnergyIncludingRigidBodies(currentAdvancedState);// DOESN'T WORK WITH OPENMM
}

// 
SimTK::Real World::CalcPotentialEnergy(void)
{
	SimTK::State& currentAdvancedState = integ->updAdvancedState();
	updateAtomListsFromCompound(currentAdvancedState);

	// Set old potential energy of the new world via DuMM !!!
	return forces->getMultibodySystem().calcPotentialEnergy(currentAdvancedState);
}

// Calculate Fixman potential
SimTK::Real World::calcFixman(void)
{
    SimTK::State& currentAdvancedState = integ->updAdvancedState();
    updateAtomListsFromCompound(currentAdvancedState); // for det(MBAT)
	SimTK::Real Fixman = updSampler(0)->calcFixman(currentAdvancedState);
	return Fixman;
}

/**
 *  Generate a proposal
 **/
bool World::generateProposal(void)
{
	// Update Robosample bAtomList
	SimTK::State& currentAdvancedState = integ->updAdvancedState();
	updateAtomListsFromCompound(currentAdvancedState);

	// Print message to identify this World
	std::cout << "World " << ownWorldIndex 
		<< ", NU " << currentAdvancedState.getNU() << ":\n";

	// GENERATE a proposal
	bool validated = updSampler(0)->reinitialize(currentAdvancedState);	
	validated = updSampler(0)->propose(currentAdvancedState) && validated;

	return validated;
}

/**
 *  Generate a number of samples
 * */
bool World::generateSamples(int howMany)
{
	// Update Robosample bAtomList
	SimTK::State& currentAdvancedState = integ->updAdvancedState();
	updateAtomListsFromCompound(currentAdvancedState);

	// Print message to identify this World
	std::cout << "World " << ownWorldIndex 
		<< ", NU " << currentAdvancedState.getNU() << ":\n";

	// GENERATE the requested number of samples
	bool validated = updSampler(0)->reinitialize(currentAdvancedState);

	for(int k = 0; k < howMany; k++) {
		validated = updSampler(0)->sample_iteration(currentAdvancedState) && validated;
	}

	// Return the number of accepted samples
	return validated;
}

/** Print information about Simbody systems. For debugging purpose. **/
void World::PrintSimbodyStateCache(SimTK::State& someState){
	std::cout << " System Stage: " << someState.getSystemStage() << std::endl;
	for(int i = 0; i < someState.getNumSubsystems(); i++){
		std::cout << " Subsystem " << i
			<< " Name: " << someState.getSubsystemName(SimTK::SubsystemIndex(i))
			<< " Stage: " << someState.getSubsystemStage(SimTK::SubsystemIndex(i))
			<< " Version: " << someState.getSubsystemVersion(SimTK::SubsystemIndex(i))
			<< std::endl;
	}
}

/** Fill Worlds Cartesian coordinates buffers.
To be called before use of getXs, getYs or getZs **/
void World::updateCoordBuffers()
{
	int allAtIx = -1;
	for(std::size_t tIx = 0; tIx < topologies->size(); tIx++){
		for(int aIx = 0; aIx < ((*topologies)[tIx]).getNAtoms(); aIx++){
			allAtIx++;
			Xs[allAtIx] = (((*topologies)[tIx]).bAtomList[aIx]).getX();
			Ys[allAtIx] = (((*topologies)[tIx]).bAtomList[aIx]).getY();
			Zs[allAtIx] = (((*topologies)[tIx]).bAtomList[aIx]).getZ();
		}
	}

}

/** Get the coordinates from buffers **/
std::vector<SimTK::Real> World::getXs()
{
	return Xs;
}

std::vector<SimTK::Real> World::getYs()
{
	return Ys;
}

std::vector<SimTK::Real> World::getZs()
{
	return Zs;
}

/** Return own CompoundSystem **/
CompoundSystem *World::getCompoundSystem() const {
	return compoundSystem.get();
}

/** Set own Compound system **/
// TODO find a solution for the old one
void World::setCompoundSystem(CompoundSystem *inCompoundSystem) {
	compoundSystem.reset(inCompoundSystem);
}

// void World::initializeTaskSpace(SimTK::CompoundSystem &compoundSystem, SimTK::GeneralForceSubsystem& force, SimTK::SimbodyMatterSubsystem& matter) {
	

// 	// compoundSystem.realizeTopology();
// 	// compoundSystem.realize()
// }

// void World::getLocationsForTaskSpace() {

// }


void World::setSamplesPerRound(int samples) {
	samplesPerRound = samples;
}

int World::getSamplesPerRound() const {
	return samplesPerRound;
}

// void World::setDistortOption(int distort) {
// 	distortOption = distort;
// }

// int World::getDistortOption() const {
// 	return distortOption;
// }

