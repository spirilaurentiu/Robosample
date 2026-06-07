// Key Folders

Molmodel/ Simbody01/ openmm/ and ./

// ### File: tests_inputs/simulate.py
import robosample

run_type = getattr(robosample.RunType, args.runType)
context = robosample.Context(args.name, args.seed, 0, 1, run_type, 1, 0)

context.loadAmberSystem(args.top, args.rst7)

# World means a Simbody robot that represents a molecule

context.addWorld(False, 1, robosample.RootMobility.CARTESIAN, flexes_Cart, True, False, 0)
context.addWorld(FIXMAN_TORQ, 1, robosample.RootMobility.WELD, flexibilities, True, False, 0)

context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
context.getWorld(worldIx).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

context.addReplicasAndLoadCoordinates(args.top, args.rstDir, nof_replicas)

for replIx in range(nof_replicas):
    context.addThermodynamicState(replIx,
                temperatures[replIx],
                accept_reject_modes,
                distort_options,
                distort_args,
                flow,
                work,
                integrators,
                worldIndexes,
                timesteps,
                mdsteps)

context.Initialize()

context.RunREX(args.equilSteps, args.prodSteps)

// ### File: src/PyBind11.cpp

Contains all the bindings from Python. For Example

py::class_<Context>(m, "Context").def(py::init<const std::string&, uint32_t, uint32_t, uint32_t, RUN_TYPE, uint32_t, uint32_t>())


// ### File: include/Topology.hpp

class Topology : public SimTK::Compound{
public:

	Topology(); [...]}

// ### File: src/Context.cpp

enum class RUN_TYPE : int {[...]};
Context::Context(const std::string& baseName_arg, uint32_t seed, uint32_t threads, uint32_t nofRoundsTillReblock, RUN_TYPE runType, uint32_t swapFreq, uint32_t swapFixmanFreq)
{[...]}

void Context::loadAmberSystem(const std::string& prmtop, const std::string& inpcrd) {
	
	AmberReader reader;
	reader.readAmberFiles(inpcrd, prmtop);

	loadAtoms(reader);
	loadBonds(reader);
	loadAngles(reader);
	loadTorsions(reader);

	// Calculate InternalCoordinates BONDS / BAT graphs
	calc_Gmolmodel_Graph(); 

	// Construct a Compound for every atom
	setAtomsCompounds();

	// Biotype will be used to look up molecular
	// force field specific parameters for an atom type
	addBiotypes();

	// // Calculate InternalCoordinates BONDS / BAT graphs
	// calc_Gmolmodel_Graph();

	// bBonds to BAT / InternalCoordinates bonds map
	load_BONDS_to_bonds( internCoords.getBonds() );

	// InternalCoordinates -> Molmodel Compound graphs with bondAtom
	build_Molmodel_AcyclicGraphs();

	// Close rings
	addRingClosingBonds_All(); 

	// Generate Topologies sub_array views (also sort bonds - not BONDS)
	generateTopologiesSubarrays();

	// Match Compounds configurations to atoms Cartesian coords
	matchDefaultConfigurations();
}


// Set Gmolmodel atoms properties from a reader: number, name, element, etc
void Context::loadAtoms(const AmberReader& reader) {


	// Iterate through atoms and set as much as possible from amberReader
	for(int aCnt = 0; aCnt < natoms; aCnt++) {

		atoms[aCnt].setDummAtomClassIndex(SimTK::DuMM::AtomClassIndex(aCnt));
		atoms[aCnt].setChargedAtomTypeIndex(SimTK::DuMM::ChargedAtomTypeIndex(aCnt));

		constexpr SimTK::Real chargeMultiplier = 18.2223;
		atoms[aCnt].setCharge(reader.getAtomsCharge(aCnt) / chargeMultiplier);

		// Set coordinates in nm (AMBER uses Angstroms)
		atoms[aCnt].setCartesians(
			reader.getAtomsXcoord(aCnt) / 10.0, );

	}
}

void Context::loadBonds(const AmberReader& reader) {
	
	// Allocate memory for bonds list
	nbonds = reader.getNumberBonds();
	bonds.reserve(nbonds);

	// Iterate bonds and get atom indeces
	// This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	for(int bCnt = 0; bCnt < nbonds; bCnt++) {

		bond.setForceK( reader.getBondsForceK(bCnt) );
		bonds.push_back(bond);
		atoms[bond.i].addNeighborIndex(bonds[bCnt].j);
		atoms[bond.j].addNeighborIndex(bonds[bCnt].i);
	}

}

// Set dumm angle list to be used in World generateDummParams.
void Context::loadAngles(const AmberReader& reader) {

	for (int i = 0; i < reader.getNumberAngles(); i++) {
		DUMM_ANGLE angle;
		angle.first = reader.getAnglesAtomsIndex1(i);

		angle.k = reader.getAnglesForceK(i);

		angle.equil = static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * angle.equil));

		dummAngles.push_back(angle);
	}
}

// Get BAT graphs
void Context::calc_Gmolmodel_Graph(){

	// Find a root in the unvisited atoms and build BAT graphs
	nofMols = 0;
	while( internCoords.computeRoot( getAtoms() )){ // find a root

		// Compute the new molecule's BAT coordinates
		internCoords.computeBAT( getAtoms() );

	}

	internCoords.computeLevelsAndOffsets( getAtoms() );

}


// Create a Comopund for each atom
void Context::setAtomsCompounds() {
	for(auto& atom : atoms) {
		const std::string& currAtomName = atom.getName();
		const int atomicNumber = atom.getAtomicNumber();

		// Create Compound
		atom.setAtomCompound(
			elementCache.getElement(atomicNumber, mass));

	}
}


// Define biotypes for each atom
void Context::addBiotypes() {

	// Set a residue name
	std::string resName = "MOL0"; // TODO: delete

	// Define atoms' Biotypes and BiotypeIndexes with their indeces and names
	int aCnt = -1;
	for(auto& atom : atoms) {

		SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
			elementCache.getElement(atom.getAtomicNumber(), atom.getMass()),
			atom.getNBonds(),
			atom.getResidueName().c_str(),
			atom.getName().c_str(),
			SimTK::Ordinality::Any
		);

	}
}


// Build Molmodel graphs with bondAtom
void Context::build_Molmodel_AcyclicGraphs(void){
	for(unsigned int molIx = 0; molIx < nofMols; molIx++){

		Topology topology("MOL" + std::to_string(++moleculeCount));

		const int rootAmberIx = internCoords.getRoot( molIx ).first;

		setRootAtom( topology, rootAmberIx );

		buildAcyclicGraph(topology, rootAmberIx, molIx);

		topologies.push_back(topology);

	}
}


// Match Compounds configurations to atoms Cartesian coords
void Context::matchDefaultConfigurations(void){

	for(unsigned int molIx = 0; molIx < nofMols; molIx++){
		matchDefaultConfigurationFromAtomsCoords(topology, molIx);
	}

}


// Assign Compound coordinates by matching bAtomList coordinates
void Context::matchDefaultConfigurationFromAtomsCoords(Topology& topology, int molIx){
	std::map<Compound::AtomIndex, SimTK::Vec3> atomTargets;
	array_view<std::vector<bSpecificAtom>::iterator>& topoSubAtomList =
		topology.subAtomList;

	for(int ix = 0; ix < topology.getNumAtoms(); ++ix){
		
		SimTK::Vec3 atomCoords(
			topoSubAtomList[ix].getX(),
);

		atomTargets.insert(std::pair<Compound::AtomIndex, SimTK::Vec3> (
			topoSubAtomList[ix].getCompoundAtomIndex(),
			atomCoords));

	}

	// Match coordinates

	topology.setTopLevelTransform(Transform(Rotation(), topLevelShift));

	topology.matchDefaultBondLengths(atomTargets);
	topology.matchDefaultAtomChirality(atomTargets, 0.01, flipAllChirality);
	topology.matchDefaultBondAngles(atomTargets);
	topology.matchDefaultDirections(atomTargets);
	topology.matchDefaultDihedralAngles(atomTargets, SimTK::Compound::DistortPlanarBonds);
	topology.matchDefaultTopLevelTransform(atomTargets);

}


