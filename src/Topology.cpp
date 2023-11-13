/**@file
Implementation of Topology class. **/

#include "Topology.hpp"

using namespace SimTK;

void Topology::BAT() {
	// std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>BAT BEGIN" << std::endl;

	// InternalCoordinates bat;
	// bat.compute(bAtomList);

	// for (const auto& b : bat.getBonds()) {
	// 	std::cout << "bond " << b.first << " " << b.second << std::endl;
	// }

	// for (const auto& a : bat.getAngles()) {
	// 	std::cout << "angle " << a.first << " " << a.second << " " << a.third << " " << std::endl;
	// }

	// for (const auto& t : bat.getTorsions()) {
	// 	std::cout << "torsion " << t.first << " " << t.second << " " << t.third << " " << t.fourth << std::endl;
	// }

	// for (int i = 0; i < bAtomList.size(); i++) {
	// 	std::cout << i << " " << bat.amber2BAT(i) << " " << bat.BAT2amber(bat.amber2BAT(i)) << std::endl;
	// }

	// std::cout << ">>>>>>>>>>>>>>>>>>>>>>>>BAT END" << std::endl;
}

void Topology::load(const readAmberInput& reader) {
	// natoms = reader.getNumberAtoms();
	// bAtomList.resize(natoms);

	// std::cout << "natoms " << natoms << std::endl;

	// // iterate through atoms as read from the prmtop file
	// for(int amberId = 0; amberId < natoms; amberId++){
	// 	bSpecificAtom atom;

	// 	// Assign a "unique" name. The generator is however limited.
	// 	// AAAA, AAAB, AAAC, AAAD, AAAF etc - why do we need this?
	// 	const std::string atomName = GetUniqueName(amberId);
	// 	atom.setName(atomName);

	// 	const int atomicNumber = reader.getAtomicNumber(amberId);
	// 	const int atomicMass = reader.getAtomsMass(amberId);

	// 	std::string elementName, elementSymbol;
	// 	GetElement(atomicNumber, elementName, elementSymbol);

	// 	// assign amber id from prmtop
	// 	atom.setNumber(amberId);

	// 	// N, H1, OX etc
	// 	auto atomType = reader.getAtomsName(amberId);
	// 	ROBOSAMPLE::trim(atomType);
	// 	atom.setInName(atomType);
	// 	atom.setFfType(atomType);

	// 	// Set charge as it is used in Amber 						????????????????????????
	// 	constexpr SimTK::Real chargeMultiplier = 18.2223;
	// 	atom.setCharge(reader.getAtomsCharge(amberId) / chargeMultiplier);

	// 	// Set coordinates in nm
	// 	atom.setX(reader.getAtomsXcoord(amberId) / 10.0);
	// 	atom.setY(reader.getAtomsYcoord(amberId) / 10.0);
	// 	atom.setZ(reader.getAtomsZcoord(amberId) / 10.0);

	// 	// Set mass
	// 	atom.setMass(reader.getAtomsMass(amberId));
	// 	atom.setAtomicNumber(atomicNumber);

	// 	// Set Lennard-Jones parameters
	// 	atom.setVdwRadius(reader.getAtomsRVdW(amberId));
	// 	atom.setLJWellDepth(reader.getAtomsEpsilon(amberId));

	// 	// Set residue name and index
	// 	atom.residueName = std::string("UNK");
	// 	atom.residueIndex = 1;

	// 	atom.setAtomCompoundType(atomName, atomicNumber, elementName, elementSymbol, atomicMass);
		
	// 	bAtomList.push_back(atom);
	// } // END atom properties

	// // SetGmolBondingPropertiesFromReader();
	// assert( (!bAtomList.empty()) &&
	// "Topology::loadAtomAndBondInfoFromReader: atom list empty.");

	// // Alloc memory for bonds list
	// nbonds = reader.getNumberBonds();
	// bonds.resize(nbonds);

	// // Iterate bonds and get atom indeces
	// // This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	// for(int i=0; i<nbonds; i++){
	// 	bBond& b = bonds[i];
	// 	b.setIndex(i);
	// 	b.i = reader.getBondsAtomsIndex1(i);
	// 	b.j = reader.getBondsAtomsIndex2(i);
	// 	// bonds.push_back(b);

	// 	bAtomList[b.i].addNeighbor(&bAtomList[b.j]);
	// 	bAtomList[b.i].addBond(&bonds[i]);

	// 	bAtomList[b.j].addNeighbor(&bAtomList[b.i]);
	// 	bAtomList[b.j].addBond(&bonds[i]);

	// 	std::cout << bonds[i].i << " " << bonds[i].j << std::endl;
	// }

	// // Assign the number of bonds an atom has and set the number of freebonds
	// // equal to the number of bonds for now
	// for (int i = 0; i < natoms; i++) {
	// 	bAtomList[i].nbonds = bAtomList[i].bondsInvolved.size();
	// 	bAtomList[i].freebonds = bAtomList[i].bondsInvolved.size();
	// }

	// // SetGmolAtomsCompoundTypes();
	// const Angle TetrahedralAngle = 109.47 * Deg2Rad;

	// // Set Gmolmodel name and element and inboard length
	// for(auto& atom : bAtomList) {
	// 	// Add BondCenters
	// 	if(atom.nbonds > 0){
	// 		if (atom.nbonds == 1){
	// 			atom.addFirstBondCenter("bond1", atom.name);
	// 		} else {
	// 			atom.addFirstTwoBondCenters("bond1", "bond2", atom.name, UnitVec3(1, 0, 0), UnitVec3(-0.5, 0.866025, 0.0));
	// 			if (atom.nbonds > 2) {
	// 				atom.addLeftHandedBondCenter("bond3", atom.name, TetrahedralAngle, TetrahedralAngle);
	// 			}
	// 			if (atom.nbonds > 3){
	// 				atom.addRightHandedBondCenter("bond4", atom.name, TetrahedralAngle, TetrahedralAngle);
	// 			}
	// 		}

	// 		// Set the inboard BondCenter
	// 		atom.setInboardBondCenter("bond1");
	// 		atom.setDefaultInboardBondLength(0.19);
	// 	}

	// 	// shouldn't they be non-redundant?
	// 	SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
	// 		// atom.compoundSingleAtom->getAtomElement(atom.getName()),
	// 		SimTK::Element(
	// 			  atom.getAtomicNumber(),
	// 			  atom.getElem().c_str(),
	// 			  atom.getElem().c_str(),
	// 			  atom.getMass()),
	// 		atom.getNBonds(),
	// 		atom.residueName.c_str(),
	// 		atom.getName().c_str(),
	// 		SimTK::Ordinality::Any
	// 	);

	// 	atom.setBiotypeIndex(biotypeIndex);

	// 	std::cout << atom.getAtomicNumber() << " " << atom.getMass() << std::endl;

	// 	// Assign atom's biotype as a composed name: name + force field type
	// 	//std::string temp(bAtomList[i].name); // @ Restore MULMOL
	// 	std::string temp(name + atom.name); // DEL MULMOL
	// 	temp += atom.fftype;
	// 	atom.setBiotype(temp);

	// 	std::cout << "Added Biotype " << temp << " with BiotypeIndex " << biotypeIndex << std::endl;
	// }

	// std::cout << "FINISHED" << std::endl;

	// // bat
}

/** Default constructor.Sets the name of this molecule to 'no_name '.
The name has no particular function and is not guaranteed to be unique **/
Topology::Topology(){
	this->name = std::string("no_name");
	this->setCompoundName((this->name));
}

/** Constructor that sets the name of the molecule. The name has no particular
function and is not guaranteed to be unique **/
Topology::Topology(std::string nameOfThisMolecule){
	this->name = nameOfThisMolecule;
	this->setCompoundName((this->name));
}

/** Default destructor. It deallocates bAtomType of every atom in the bAtomList
because we want to allow the valence to change during the simulation
e.g. semi-grand canonical ensemble. **/
Topology::~Topology() {
	for (auto& atom : bAtomList) {
		atom.destroy();
	}
}

/** Set Gmolmodel atoms properties from a reader: number, name, element,
 * initial name, force field type, charge, coordinates, mass, LJ parameters.
 * 1-to-1 correspondence between prmtop and Gmolmodel.
 * This does not set anything in Compund or DuMM.  **/
void Topology::SetGmolAtomPropertiesFromReader(readAmberInput *amberReader)
{
	// Alloc memory for atoms and bonds list
	natoms = amberReader->getNumberAtoms();
	bAtomList.resize(natoms);

	// Iterate through atoms and set as much as possible from amberReader
	for(int i = 0; i < natoms; i++){

		// Assign an index like in prmtop
		bAtomList[i].setNumber(i);

		// This is the name of the atom in the .prmtop file
		// Examples: "O1", "C1", "C2", "H1", "H10"
		const std::string initialName = amberReader->getAtomsName(i);
		
		// Set element
		const int atomicNumber = amberReader->getAtomicNumber(i);
 		bAtomList[i].setAtomicNumber(atomicNumber);
		bAtomList[i].setElem(elementCache.getSymbolByAtomicNumber(atomicNumber));

		// Assign a "unique" name. The generator is however limited.
		bAtomList[i].generateName(i);

		// Store the initial name from prmtop
		bAtomList[i].setInName(initialName);

		// Set force field atom type
		bAtomList[i].setFfType(initialName);

		// Set charge as it is used in Amber
		SimTK::Real chargeMultiplier = 18.2223;
		bAtomList[i].setCharge(amberReader->getAtomsCharge(i) / chargeMultiplier);

		// Set coordinates in nm
		bAtomList[i].setX(amberReader->getAtomsXcoord(i) / 10.0);
		bAtomList[i].setY(amberReader->getAtomsYcoord(i) / 10.0);
		bAtomList[i].setZ(amberReader->getAtomsZcoord(i) / 10.0);
		bAtomList[i].setCartesians(
			amberReader->getAtomsXcoord(i) / 10.0,
			amberReader->getAtomsYcoord(i) / 10.0,
			amberReader->getAtomsZcoord(i) / 10.0 );

		// Set mass
		const int mass = amberReader->getAtomsMass(i);
		bAtomList[i].setMass(mass);

		// Set Lennard-Jones parameters
		bAtomList[i].setVdwRadius(amberReader->getAtomsRVdW(i));
		bAtomList[i].setLJWellDepth(amberReader->getAtomsEpsilon(i));

		// Set residue name and index
		bAtomList[i].setResidueName(std::string("UNK"));
		bAtomList[i].residueIndex = 1;

		// Add element to cache
		elementCache.addElement(atomicNumber, mass);
	}
}

/** Set bonds properties from reader: bond indeces, atom neighbours.
 *  1-to-1 correspondence between prmtop and Gmolmodel.
 **/
void Topology::SetGmolBondingPropertiesFromReader(readAmberInput *amberReader)
{
	assert( (!bAtomList.empty()) && "Topology::loadAtomAndBondInfoFromReader: atom list empty.");

	// Alloc memory for bonds list
	nbonds = amberReader->getNumberBonds();
	bonds.resize(nbonds);

	// Iterate bonds and get atom indeces
	// This establishes a 1-to-1 correspondence between prmtop and Gmolmodel
	for(int i=0; i<nbonds; i++){
		bonds[i].setIndex(i);
		bonds[i].i = amberReader->getBondsAtomsIndex1(i);
		bonds[i].j = amberReader->getBondsAtomsIndex2(i);

		// BAD! this will break if we invalidate the bonds vector
		bAtomList[bonds[i].i].addNeighbor(&bAtomList[bonds[i].j]);
		bAtomList[bonds[i].i].addBond(&bonds[i]);

		bAtomList[bonds[i].j].addNeighbor(&bAtomList[bonds[i].i]);
		bAtomList[bonds[i].j].addBond(&bonds[i]);
	}

	// Assign the number of bonds an atom has and set the number of freebonds
	// equal to the number of bonds for now
	for(int i = 0; i < natoms ; i++) {
		bAtomList[i].setNbonds(bAtomList[i].bondsInvolved.size());
		bAtomList[i].setFreebonds(bAtomList[i].bondsInvolved.size());
	}
}

/** Set atoms Molmodel types (Compound::SingleAtom derived) based on
 * their valence **/
void Topology::SetGmolAtomsCompoundTypes(){

	Angle TetrahedralAngle = 109.47 * Deg2Rad;

	// Set Gmolmodel name and element and inboard length
	for(auto& atom : bAtomList) {
		
		// bond centers must have been already set
		const std::string& currAtomName = atom.getName();
		const int atomicNumber = atom.getAtomicNumber();
		const int mass = atom.getMass();
		atom.setAtomCompoundType(currAtomName, elementCache.getElement(atomicNumber, mass));
		
		// Add BondCenters
		const int currAtomNBonds = atom.getNBonds();
		if(currAtomNBonds > 0){
			if (currAtomNBonds == 1){
				atom.addFirstBondCenter("bond1", currAtomName);
			} else {
				atom.addFirstTwoBondCenters("bond1", "bond2", currAtomName, UnitVec3(1, 0, 0), UnitVec3(-0.5, 0.866025, 0.0));
				if (currAtomNBonds > 2) {
					atom.addLeftHandedBondCenter("bond3", currAtomName, TetrahedralAngle, TetrahedralAngle);
				}
				if (currAtomNBonds > 3){
					atom.addRightHandedBondCenter("bond4", currAtomName, TetrahedralAngle, TetrahedralAngle);
				}
			}

			// Set the inboard BondCenter
			atom.setInboardBondCenter("bond1");
			atom.setDefaultInboardBondLength(0.19);
		}
	}
}

/** Reads information from a readAmberInput object and put it in
 * bAtomList and bonds lists **/
void Topology::loadAtomAndBondInfoFromReader(readAmberInput *amberReader)
{
	// Set atoms properties from a reader: number, name, element, initial
	// name, force field type, charge, coordinates, mass, LJ parameters
	SetGmolAtomPropertiesFromReader(amberReader);

	// Set bonds properties from reader: bond indeces, atom neighbours
	SetGmolBondingPropertiesFromReader(amberReader);

	// Set atoms Molmodel types (Compound::SingleAtom derived) based on
	// their valence
	SetGmolAtomsCompoundTypes();
}

/** Print atom and bonds list with details**/
void Topology::PrintAtomList(int whichWorld)
{
	// // Atoms
	// std::cout<<"Topology::PrintAtomList\n";
	// for(unsigned int i = 0; i < bAtomList.size(); i++){
	// 	bAtomList[i].Print(whichWorld);
	// }

/* 	// Bonds
	for(unsigned int i = 0; i < bonds.size(); i++){
		bonds[i].Print(whichWorld);
	} */
}

/** Biotype is a Molmodel hook that is usually used to look up molecular
force field specific parameters for an atom type. Gmolmodel defines a
new Biotype for each atom. The only thing that is specified is the element
with info about name, atomic number, valence and mass. **/
void Topology::bAddBiotypes(
	readAmberInput *amberReader
)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Iterate atoms and define Biotypes with their indeces and names
	for(auto& atom : bAtomList){
		SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
			elementCache.getElement(atom.getAtomicNumber(), atom.getMass()),
			atom.getNBonds(),
			resName.c_str(),
			(atom.getName()).c_str(),
			SimTK::Ordinality::Any
		);

		atom.setBiotypeIndex(biotypeIndex);

		// Assign atom's biotype as a composed name: name + force field type
		std::string biotype = resName + atom.getName() + atom.getFftype();
		atom.setBiotype(biotype);
	}
}

/** The following functions are used to build the molecular graph using bonding
information from bonds list and bondsInvolved list of each atom in bAtomList.
**/

/** The actual recursive function that builds the graph **/
void Topology::buildAcyclicGraph(bSpecificAtom *node, bSpecificAtom *previousNode)
{
	// The base atom has to be set once Molmodel
	baseSetFlag = 0;

	// Only process unvisited nodes
	if( node->wasVisited() ){
		return;
	}

	// Mark the depth of the recursivity
	++nofProcesses;

	// Mark Gmolmodel bond and create bond in Molmodel
	for(std::vector<bBond *>::iterator bondsInvolvedIter = (node->bondsInvolved).begin();
		bondsInvolvedIter != (node->bondsInvolved).end(); ++bondsInvolvedIter)
	{
		// Check if there is a bond between prevnode and node based on bonds
		// read from amberReader
		if ((*bondsInvolvedIter)->isThisMe(node->getNumber(), previousNode->getNumber()) ) {
			(*bondsInvolvedIter)->setVisited(1);

			// Skip the first step as we don't have yet two atoms
			if (nofProcesses != 1) {

				// The first bond is special in Molmodel and has to be
				// treated differently. Set a base atom first
				if (nofProcesses == 2) {
					if (baseSetFlag == 0) {
						this->setBaseAtom(previousNode->getSingleAtom());
						this->setAtomBiotype(previousNode->getName(), (this->name), previousNode->getName());
						this->convertInboardBondCenterToOutboard();
						baseSetFlag = 1;
					}
				}

				// Bond current node by the previous (Compound function)
				std::stringstream parentBondCenterPathName;
				if (previousNode->getNumber() == baseAtomNumber) {
					parentBondCenterPathName << previousNode->getName()
						<< "/bond" << previousNode->getFreebonds();
				} else {
					parentBondCenterPathName << previousNode->getName()
						<< "/bond" << (previousNode->getNBonds() - previousNode->getFreebonds() + 1);
				}

				// THIS IS WHERE WE PERFORM THE ACTUAL BONDING
				// (Compound::SingleAtom&, BondCenterPathName, Length, Angle
				std::string debugString = parentBondCenterPathName.str();
				this->bondAtom(node->getSingleAtom(),
						(parentBondCenterPathName.str()).c_str(), 0.149, 0);

				// Set the final Biotype
				this->setAtomBiotype(node->getName(), (this->name).c_str(), node->getName());

				// Set bSpecificAtom atomIndex to the last atom added to bond
				node->setCompoundAtomIndex(getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1));


				// The only time we have to set atomIndex to the previous node
				if (nofProcesses == 2) {
					previousNode->setCompoundAtomIndex(getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0));
				}

				// Set bBond Molmodel Compound::BondIndex
				(*bondsInvolvedIter)->setBondIndex(Compound::BondIndex(getNumBonds() - 1));
				std::pair<SimTK::Compound::BondIndex, int> pairToBeInserted(
						Compound::BondIndex(getNumBonds() - 1),
						(*bondsInvolvedIter)->getIndex()
				);

				bondIx2GmolBond.insert(pairToBeInserted);

				GmolBond2bondIx.insert( std::pair<int, SimTK::Compound::BondIndex>(
						(*bondsInvolvedIter)->getIndex(),
						Compound::BondIndex(getNumBonds() - 1)
				) );

				//std::cout << "DEBUG inserted into GmolBond2bondIx bondIx2GmolBond  " 
				//<< GmolBond2bondIx.size() << " " << bondIx2GmolBond.size() 
				//<< std::endl << std::flush;

				// Drop the number of available bonds
				/* --(previousNode->freebonds);
				--(node->freebonds); */

				previousNode->decrFreebonds();
				node->decrFreebonds();

				// Bond was inserted in Molmodel Compound. Get out and search
				// the next bond
				break;

			}
		}
	}

	// Mark the node as visited
	node->setVisited(1);

	// Set the previous node to this node
	previousNode = node;

	// Go to the next node. Choose it from his neighbours.
	for(unsigned int i = 0; i < (node->neighbors).size(); i++) {
		buildAcyclicGraph((node->neighbors)[i], previousNode);
	}

}

/** After building the acyclic molecular tree close the remaining bonds **/
void Topology::addRingClosingBonds() {
	// Consider all remaining bonds ring closing bonds and close them
	for(int i=0; i<nbonds; i++){
		if(bonds[i].isVisited() == 0){

			bSpecificAtom *leftNode  =  &(bAtomList[bonds[i].i]);
			bSpecificAtom *rightNode =  &(bAtomList[bonds[i].j]);

			std::stringstream sbuff;
			if(leftNode->getNumber() == baseAtomNumber){
				sbuff << leftNode->getName() << "/bond" << leftNode->getFreebonds();
			}else{
				sbuff << leftNode->getName() << "/bond"
					<< (leftNode->getNBonds() - leftNode->getFreebonds() + 1);
			}

			std::stringstream otsbuff;
			if(rightNode->getNumber() == baseAtomNumber){
				otsbuff << rightNode->getName() << "/bond" << rightNode->getFreebonds();
			}else{
				otsbuff << rightNode->getName() << "/bond"
					<< (rightNode->getNBonds() - rightNode->getFreebonds() + 1);
			}

			this->addRingClosingBond(
					(sbuff.str()).c_str(),
					(otsbuff.str()).c_str(),
					0.14,
					109*Deg2Rad,
					BondMobility::Rigid // CHANGE
			);

			// Set bBond Molmodel Compound::BondIndex // TODO is this necessary ?
			bonds[i].setBondIndex(Compound::BondIndex(getNumBonds() - 1));
			bonds[i].setAsRingClosing();
			std::pair<SimTK::Compound::BondIndex, int> pairToBeInserted(
					Compound::BondIndex(getNumBonds() - 1),
					bonds[i].getIndex()
			);

			bondIx2GmolBond.insert(pairToBeInserted);

			GmolBond2bondIx.insert( std::pair<int, SimTK::Compound::BondIndex>(
					bonds[i].getIndex(),
					Compound::BondIndex(getNumBonds() - 1)
			) );
			////////////////////////////////////////////

			// Compound::setAtomBiotype(Compound::AtomPathName,
			// biotypeResidueName, biotypeAtomName
			// SimTK::Ordinality::Residue = SimTK::Ordinality::Any)
			this->setAtomBiotype(leftNode->getName(), (this->name), leftNode->getName());
			this->setAtomBiotype(rightNode->getName(), (this->name), rightNode->getName());

			/* --leftNode->freebonds;
			--rightNode->freebonds; */

			leftNode->decrFreebonds();
			rightNode->decrFreebonds();


		}
	}
}

/** Match Default configuration with the coordinates loaded from
 * the input reader **/
void Topology::matchDefaultConfigurationWithAtomList(
		SimTK::Compound::MatchStratagem matchStratagem)
{
	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].getCompoundAtomIndex(), vec));
	}

	matchDefaultConfiguration(atomTargets, matchStratagem, true, 150.0);
	std::cout << "Gmolmodel Topology match done" << std::endl;
}

/** Builds the molecular tree, closes the rings, matches the configuration
on the graph using using Molmodels matchDefaultConfiguration and sets the
general flexibility of the molecule. **/
// TODO: break in two functions
void Topology::buildGraphAndMatchCoords(int argRoot)
{
	// Initialize all atoms and bonds to unvisited
	for (int i = 0; i < natoms; i++) {
		bAtomList[i].setVisited(0);
	}
	for (int i = 0; i < nbonds; i++) {
		bonds[i].setVisited(0);
	}

	// Find an atom to be the root. It has to have more than one bond
	bSpecificAtom *root = nullptr;
	if ((static_cast<size_t>(argRoot) > bAtomList.size()) || (bAtomList[argRoot].getNBonds() > 1)) {
		baseAtomNumber = argRoot;
		root = &(bAtomList[argRoot]);
		bSpecificAtomRootIndex = argRoot;
	}else {
		std::cout << "Root atom will be chosen by Gmolmodel...  ";
		int baseAtomListIndex = 0;
		for (int i = 0; i < natoms; i++) {
			if (bAtomList[i].getNBonds() > 1) {
				baseAtomListIndex = i;
				std::cout << "done. Root chosen " << i << std::endl;
				break;
			}
		}

		root = &(bAtomList[baseAtomListIndex]);
		bSpecificAtomRootIndex = baseAtomListIndex;
		baseAtomNumber = root->getNumber();
	}

	// Build the graph
	if(bAtomList.size() == 1){
		this->setBaseAtom(bAtomList[0].getSingleAtom());
		(bAtomList[0]).setCompoundAtomIndex(SimTK::Compound::AtomIndex(0));
		std::cout << "Topology::buildGraphAndMatcoords single atom done\n" << std::flush;
	}else{
		nofProcesses = 0;
		buildAcyclicGraph(root, root);
		std::cout << "Topology::buildGraphAndMatcoords buildAcyclicGraph done\n" << std::flush;
	}

	// Close the remaining bonds
	addRingClosingBonds();
    std::cout << "Topology::buildGraphAndMatcoords  addRingClosingBonds done\n" << std::flush;

	// Build the conformation
	matchDefaultConfigurationWithAtomList(SimTK::Compound::Match_Exact);

	// Now that everything is built, initialize aIx2TopTransforms map
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		aIx2TopTransform.insert(std::make_pair((bAtomList[i]).getCompoundAtomIndex(), SimTK::Transform()));
	}

}


/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
setBiotypeChargedAtomType for every atom. These Molmodel functions contain
information regarding the force field parameters. **/
void Topology::generateDummAtomClasses(
		std::string resName
	, readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
)
{
	// std::cout << "Topology::generateDummAtomClasses START\n";

	// Declarations
	std::vector<bool> founditInDuMM(
		amberReader->getNumberAtoms(), false);

	// Iterate Gmolmodel atoms and define AtomClasses based on
	for(int i = 0; i < amberReader->getNumberAtoms(); i++){

		// Get AtomClass parameters from bAtomList info
		AtomClassParams atomClassParams(
			bAtomList[i].getAtomicNumber(),
			bAtomList[i].getNBonds(),
			bAtomList[i].getVdwRadius() / 10.0, // nm
			bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
		);

		SimTK::DuMM::AtomClassIndex dummAtomClassIndex;
		std::string atomClassName;

		// Check if we already have this Atom Class in this Topology
		bool founditInThisTopology = false;
		std::map<AtomClassParams, AtomClassId>::const_iterator aCP2aCI_iter;
		for (aCP2aCI_iter  = aClassParams2aClassId.begin();
			 aCP2aCI_iter != aClassParams2aClassId.end();
		   ++aCP2aCI_iter)
		{
			//if(aCP2aCI_iter->first == atomClassParams)
			if((aCP2aCI_iter->second).name == bAtomList[i].getFftype())
			{
				founditInThisTopology = true;
				break;
			}
		}

		// Further search in DuMM
		if(!founditInThisTopology){
			std::string str(bAtomList[i].getFftype());
			const SimTK::String simtkstr(str);
			founditInDuMM[i] = dumm.hasAtomClass(simtkstr);
		}

		/* std::cout << "Topology::generateDummAtomClasses try AtomClassName " 
			<< " " << bAtomList[i].getFftype() << " "
			<< founditInThisTopology << " " << founditInDuMM[i];
		atomClassParams.dump(); */

		// we don't have this AtomClass
		if ( (!founditInThisTopology) && (!founditInDuMM[i]) ) {

			// Get an AtomClass index from DuMM
			dummAtomClassIndex = dumm.getNextUnusedAtomClassIndex();
		
			// Define an AtomClass name
			atomClassName = bAtomList[i].getFftype();

			//std::cout << "!foundit aCIx atomClassName " << aCIx << " "
			//	<< atomClassName  << std::endl;

			// Define an AtomClass
			dumm.defineAtomClass(dummAtomClassIndex,
				atomClassName.c_str(),
				atomClassParams.atomicNumber,
				atomClassParams.valence,
				atomClassParams.vdwRadius,
				atomClassParams.LJWellDepth
			);

			/* std::cout << "Topology::generateDummAtomClasses insert AtomClassIndex AtomClassName " 
				<< dummAtomClassIndex << " " << atomClassName ;
			atomClassParams.dump(); */

			// Insert an entry in our map too
			AtomClassId atomClassId(dummAtomClassIndex, atomClassName);
			aClassParams2aClassId.insert( std::make_pair(atomClassParams, atomClassId) );

		}else{ // we already have this AtomClass

			if(!founditInDuMM[i]){
				//aClassParams2aClassId.at(atomParams);
				//AtomClassId& atomClassId = aClassParams2aClassId.at(atomParams);
				const AtomClassId& atomClassId = aCP2aCI_iter->second;

				dummAtomClassIndex = atomClassId.dummAtomClassIndex;

				//atomClassName = atomClassId.name;
		
				//std::cout << "foundit  aCIx atomClassName " << aCIx << " " << atomClassName << std::endl;
			}else{
				dummAtomClassIndex = dumm.getAtomClassIndex(bAtomList[i].getFftype());
				//atomClassName = bAtomList[i].getFftype();
			}

		}

		// Insert AtomClass index in Gmolmodel atom list too
		bAtomList[i].setDummAtomClassIndex(dummAtomClassIndex);

	} // --------------------------- DONE AmberReader atoms



	// Define ChargedAtomTypeIndeces
	SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
	std::string chargedAtomTypeName;

	// Iterate through atoms and define DuMM charged atom types
	for(int k = 0; k < amberReader->getNumberAtoms(); k++){

		// Get a ChargedAtomType index
		chargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
		bAtomList[k].setChargedAtomTypeIndex(chargedAtomTypeIndex);

		// Define a chargedAtomType name
		chargedAtomTypeName =  resName;
		chargedAtomTypeName += bAtomList[k].getBiotype();

		// Define a ChargedAtomType (AtomClass with a charge)
		dumm.defineChargedAtomType(
		chargedAtomTypeIndex,
		chargedAtomTypeName.c_str(),
		bAtomList[k].getDummAtomClassIndex(),
		bAtomList[k].charge
		);
		/*std::cout << "Defined chargedAtomType " << chargedAtomTypeName 
			<< " with chargedAtomTypeIndex " << chargedAtomTypeIndex
			<< std::endl; */

		// Associate a ChargedAtomTypeIndex with a Biotype index
		dumm.setBiotypeChargedAtomType(
		bAtomList[k].getChargedAtomTypeIndex(),
		bAtomList[k].getBiotypeIndex()
		);

	}

}

//TODO: DELETE // NO LONGER NEEDED
/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
setBiotypeChargedAtomType for every atom. These Molmodel functions contain
information regarding the force field parameters. **/
void Topology::transferDummAtomClasses(
		std::string resName
	, readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId				
)
{

	SimTK::DuMM::AtomClassIndex aCIx;
	std::string atomClassName;

	// Iterate through AtomClasses map and put AtomClasses in Dumm
	std::map<AtomClassParams, AtomClassId>::const_iterator it;
	for (	it = aClassParams2aClassId.begin();
		it != aClassParams2aClassId.end(); ++it){

		const AtomClassParams& atomParams = it->first;
		const AtomClassId& atomClassId = it->second;

		aCIx = atomClassId.dummAtomClassIndex;
		atomClassName = atomClassId.name;

		/* std::cout << "Topology::transferAtomClasses " 
			<< aCIx << " " << atomClassName ;
		atomParams.dump(); */

		// Define an AtomClass
		dumm.defineAtomClass(aCIx, atomClassName.c_str(),
			atomParams.atomicNumber,
			atomParams.valence,
			atomParams.vdwRadius,
			atomParams.LJWellDepth
		);

	}
	// --------------------------- DONE
}

void Topology::transferDummChargedAtomClasses(
			std::string resName
	, readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId				
)
{

	// Define ChargedAtomTypeIndeces
	SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex;
	std::string chargedAtomTypeName;

	// Iterate through atoms and define DuMM charged atom types
	for(int k = 0; k < amberReader->getNumberAtoms(); k++){
		// Get a ChargedAtomType index
		chargedAtomTypeIndex = dumm.getNextUnusedChargedAtomTypeIndex();
		bAtomList[k].setChargedAtomTypeIndex(chargedAtomTypeIndex);

		// Define a chargedAtomType name
		chargedAtomTypeName =  resName;
		chargedAtomTypeName += bAtomList[k].getBiotype();

		// Define a ChargedAtomType (AtomClass with a charge)
		dumm.defineChargedAtomType(
		  chargedAtomTypeIndex,
		  chargedAtomTypeName.c_str(),
		  bAtomList[k].getDummAtomClassIndex(),
		  bAtomList[k].charge
		);
		//std::cout << "Defined chargedAtomType " << chargedAtomTypeName << " with chargedAtomTypeIndex " << chargedAtomTypeIndex << std::endl;

		// Associate a ChargedAtomTypeIndex with a Biotype index
		dumm.setBiotypeChargedAtomType(
		  bAtomList[k].getChargedAtomTypeIndex(),
		  bAtomList[k].getBiotypeIndex()
		);

	}

}

/** Print Molmodel specific types as introduced in Gmolmodel **/
const void Topology::PrintMolmodelAndDuMMTypes(
	SimTK::DuMMForceFieldSubsystem& dumm) const
{
	std::cout << "Print Molmodel And DuMM Types:" << std::endl 
	<< std::flush;
	for(size_t i = 0; i < bAtomList.size(); i++){
		std::cout << " list ix " << i
			<< " CompoundAtomIndex " << bAtomList[i].getCompoundAtomIndex()
			//<< " DuMMAtomIndex " << getDuMMAtomIndex(bAtomList[i].getCompoundAtomIndex())
			<< " biotypename " << bAtomList[i].getBiotype()
			<< " name " << bAtomList[i].getName()
			<< " BiotypeIndex " << bAtomList[i].getBiotypeIndex()
			<< " ChargedAtomTypeIndex "<< bAtomList[i].getChargedAtomTypeIndex()
			<< " AtomClassIx " << bAtomList[i].getDummAtomClassIndex()
			<< " partialChargeInE " << bAtomList[i].charge
			<< " chargedAtomTypeIndex "
			<< bAtomList[i].getChargedAtomTypeIndex()
			<< " DuMM VdW Radius "
			<< dumm.getVdwRadius(bAtomList[i].getDummAtomClassIndex())
			<< " DuMM VdW Well Depth "
			<< dumm.getVdwWellDepth(bAtomList[i].getDummAtomClassIndex())
			<< std::endl << std::flush;
	}
}

/** Calls DuMM defineBondStretch to define bonds parameters. **/
void Topology::bAddDummBondParams(std::string, 
	readAmberInput *amberReader, 
	SimTK::DuMMForceFieldSubsystem& dumm,
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs)
{

	// Keep track of inserted AtomClass pairs
	//std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;

	// Iterate through bonds and define their parameters
	// Suppose or try to have the same order as the reader
	for(int t = 0; t < amberReader->getNumberBonds(); t++){
		
		// Generate a pair of atom classes for this bond
		std::vector<SimTK::DuMM::AtomClassIndex> thisBondACIxs;
		thisBondACIxs.push_back( SimTK::DuMM::AtomClassIndex(
			(bAtomList[bonds[t].i]).getDummAtomClassIndex()) );
		thisBondACIxs.push_back( SimTK::DuMM::AtomClassIndex(
			(bAtomList[bonds[t].j]).getDummAtomClassIndex()) );

		// Check if we already have this bond
		bool foundit = false;
		for(auto& row:allBondsACIxs){
			if ( IsTheSameBond (thisBondACIxs, row) ){
				foundit = true;	break;
			}
		}

		if (  !foundit ){ // bond was not found

			/* std::cout << "Topology::bAddDummBondParams " 
				<< thisBondACIxs[0] << " "
				<< thisBondACIxs[1] << std::endl; */
		
			dumm.defineBondStretch_KA(
				(bAtomList[bonds[t].i]).getDummAtomClassIndex(),
				(bAtomList[bonds[t].j]).getDummAtomClassIndex(),
				amberReader->getBondsForceK(t),  //k1
				amberReader->getBondsEqval(t)   //equil1
			);

			// Put the entry in our map too
			allBondsACIxs.push_back(thisBondACIxs);

		}
	}
}

/** Calls DuMM defineBondStretch to define bonds parameters. **/
void Topology::bAddDummBondParams(std::string, 
	readAmberInput *amberReader, 
	SimTK::DuMMForceFieldSubsystem& dumm)
{
/* 	// Keep track of inserted AtomClass pairs
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allBondsACIxs;

	bAddDummBondParams(this->name, amberReader, dumm, allBondsACIxs); */

}


/** Calls DuMM defineBondBend to define angle parameters. **/
void Topology::bAddDummAngleParams(std::string,
	readAmberInput *amberReader,
	SimTK::DuMMForceFieldSubsystem& dumm,
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs)
{

	// Keep track of inserted AtomClass pairs
	//std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allAnglesACIxs;

	// Iterate angles and define their parameters
	for(int t = 0; t < amberReader->getNumberAngles(); t++){

		// Generate a triple of atom class indexes for this angle
		std::vector<SimTK::DuMM::AtomClassIndex> thisAngleACIxs;
		thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
			bAtomList[amberReader->getAnglesAtomsIndex1(t)].getDummAtomClassIndex()) );
		thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
			bAtomList[amberReader->getAnglesAtomsIndex2(t)].getDummAtomClassIndex()) );
		thisAngleACIxs.push_back( SimTK::DuMM::AtomClassIndex(
			bAtomList[amberReader->getAnglesAtomsIndex3(t)].getDummAtomClassIndex()) );

		// Check if we already have this angle
		bool foundit = false;
		for(auto& row:allAnglesACIxs){
			if ( IsTheSameAngle (thisAngleACIxs, row) ){
				foundit = true;	break;
			}
		}

		if (  !foundit ){ // angle was not found

			/* std::cout << "Topology::bAddDummAngleParams " 
				<< thisAngleACIxs[0] << " "
				<< thisAngleACIxs[1] << " "
				<< thisAngleACIxs[2] << std::endl; */
		
			dumm.defineBondBend_KA(
				bAtomList[amberReader->getAnglesAtomsIndex1(t)].getDummAtomClassIndex(),
				bAtomList[amberReader->getAnglesAtomsIndex2(t)].getDummAtomClassIndex(),
				bAtomList[amberReader->getAnglesAtomsIndex3(t)].getDummAtomClassIndex(),
				amberReader->getAnglesForceK(t),
				static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE 
					* amberReader->getAnglesEqval(t))) // TODO 32 vs 64 bit
			);

			// Put the entry in our map too
			allAnglesACIxs.push_back(thisAngleACIxs);
		}
	}
}

/** Calls DuMM defineBondTorsion for 1, 2, 3 and 4 periodicities **/
void Topology::bAddDummTorsionParams(
	  std::string resName
	, readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>&
		allDihedralsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& 
		allImpropersACIxs

)
{
	// Keep track of inserted AtomClass pairs
	//std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> allDihedralsACIxs;

	std::cout << "Topology::bAddDummTorsionParams allDihedralsACIxs\n";
	//PrintCppVector(allDihedralsACIxs);

	// Amber reader dihedrals vector
	std::vector<std::pair<int, int>> pairStartAndLens =
		amberReader->getPairStartAndLen();

	for(unsigned int index=0; index<pairStartAndLens.size(); index++){

		// Get start and len of this dihedral
		int first    = pairStartAndLens[index].first;
		int t        = pairStartAndLens[index].first;
		int numberOf = pairStartAndLens[index].second;

		// Iterate through this dihedral definitions
		//for(int t = first; t < (first + numberOf); t++){

			// Get AtomClass indeces first
			SimTK::DuMM::AtomClassIndex aCIx1 =
				bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx2 =
				bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx3 =
				bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getDummAtomClassIndex();
			SimTK::DuMM::AtomClassIndex aCIx4 =
				bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getDummAtomClassIndex();

			/* std::cout << "Topology::bAddDummDihedralParams checking " 
				<< amberReader->getDihedralsAtomsIndex1(t) << " " << amberReader->getDihedralsAtomsIndex2(t) << " "
				<< amberReader->getDihedralsAtomsIndex3(t) << " " << amberReader->getDihedralsAtomsIndex4(t) << " :| "
				<< first << " " << t << " " << numberOf << " | " << " forceK: " 
				<< SimTK_KCAL_TO_KJOULE * amberReader->getDihedralsForceK(t) << " period " 
				<< amberReader->getDihedralsPeriod(t)
				<< " phase: " << amberReader->getDihedralsPhase(t) << " | "
				<< aCIx1 << " " << aCIx2 << " " << aCIx3 << " " << aCIx4
				<< std::endl; */

			// Check if a quad of atom indices is a normal dihedral
			// or an improper dihedral, by checking if consecutive
			// atoms are bonded 

			int amber_aIx_1 = amberReader->getDihedralsAtomsIndex1(t);
			int amber_aIx_2 = amberReader->getDihedralsAtomsIndex2(t);
			int amber_aIx_3 = amberReader->getDihedralsAtomsIndex3(t);
			int amber_aIx_4 = amberReader->getDihedralsAtomsIndex4(t);

			bool dihedral=false;
			bool improper=true;

			if (checkBond(amber_aIx_1, amber_aIx_2) &&
				checkBond(amber_aIx_2, amber_aIx_3) &&
				checkBond(amber_aIx_3, amber_aIx_4))
			{
				dihedral = true;
				improper = false;
			}
			
			// Generate a quad of atom class indexes for this dihedral, 
			// regardless of whether it's a torsion or an improper
			std::vector<SimTK::DuMM::AtomClassIndex> thisDihedralACIxs;
			thisDihedralACIxs.push_back(aCIx1);
			thisDihedralACIxs.push_back(aCIx2);
			thisDihedralACIxs.push_back(aCIx3);
			thisDihedralACIxs.push_back(aCIx4);


/* 			// Check if this is a regular dihedral angle or an improper,
			// by checking if the atoms are bonded in the order they appear.
			int amber_aIx_1 = amberReader->getDihedralsAtomsIndex1(t);
			int amber_aIx_2 = amberReader->getDihedralsAtomsIndex2(t);
			int amber_aIx_3 = amberReader->getDihedralsAtomsIndex3(t);
			int amber_aIx_4 = amberReader->getDihedralsAtomsIndex4(t);

			if (checkBond(amber_aIx_1, amber_aIx_2) &&
				checkBond(amber_aIx_2, amber_aIx_3) &&
				checkBond(amber_aIx_3, amber_aIx_4))
			{
				
				printf ("%d %d %d %d is a normal angle\n", amber_aIx_1, amber_aIx_2, amber_aIx_3, amber_aIx_4);
			}
			else
			{
				// call defineAmberImproperTorsion_KA if this is the case
				printf("%d %d %d %d is an improper angle\n", amber_aIx_1, amber_aIx_2, amber_aIx_3, amber_aIx_4);
			} */

			if (dihedral){
				/* for (int tempVar=0; tempVar<numberOf; tempVar++)
					std::cout << "found dihedral angle!\n"; */
				// If it is a normal dihedral, check if we have it in our
				// dihedral list
				bool foundit = false;
				//printf("allDihedralsACIxs.size: %d\n", allDihedralsACIxs.size());

				for(auto& row:allDihedralsACIxs)
				{
					if ( IsTheSameTorsion (thisDihedralACIxs, row))
					{
						foundit = true;	break;
					}
				}

				if (!foundit){ // dihedral was not found

					/* std::cout << "Topology::bAddDummDihedralParams insert " 
						<< thisDihedralACIxs[0] << " " << thisDihedralACIxs[1] << " "
						<< thisDihedralACIxs[2] << " " << thisDihedralACIxs[3] << " := " 
						<< numberOf << " " << amberReader->getDihedralsPeriod(t) 
						<< " ;| " << amberReader->getDihedralsForceK(t) 
						<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))); */
					/* if(numberOf == 2){std::cout
						<< " | " << amberReader->getDihedralsForceK(t+1)
						<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)));}
					if(numberOf == 3){std::cout
						<< " | " << amberReader->getDihedralsForceK(t+2)
						<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2)));}
					if(numberOf ==4){std::cout
						<< " | " << amberReader->getDihedralsForceK(t+3)
						<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+3)));}
					if(numberOf ==5){std::cout
						<< " | " << amberReader->getDihedralsForceK(t+4)
						<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+4)));}
					std::cout << std::endl << std::flush; */

					// Define the dihedrals
					if(numberOf == 1){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)))
						);
					}
					else if(numberOf == 2){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
							static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)))
						);
					}
					else if(numberOf == 3){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
							static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))),
							static_cast<int>(amberReader->getDihedralsPeriod(t + 2)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t + 2),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2)))
						);
					}else if (numberOf == 4){
						dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
							static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
							static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t + 1),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))),
							static_cast<int>(amberReader->getDihedralsPeriod(t + 2)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t + 2),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2))),
							static_cast<int>(amberReader->getDihedralsPeriod(t+ 3)), // TODO wants int, returns double
							amberReader->getDihedralsForceK(t+3),
							static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+3)))
						);
					}else if (numberOf == 5){
						// dumm.defineBondTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
						// 	static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
						// 	amberReader->getDihedralsForceK(t),
						// 	static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
						// 	static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
						// 	amberReader->getDihedralsForceK(t + 1),
						// 	static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))),
						// 	static_cast<int>(amberReader->getDihedralsPeriod(t + 2)), // TODO wants int, returns double
						// 	amberReader->getDihedralsForceK(t + 2),
						// 	static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2))),
						// 	static_cast<int>(amberReader->getDihedralsPeriod(t+ 3)), // TODO wants int, returns double
						// 	amberReader->getDihedralsForceK(t+3),
						// 	static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+3))),
						// 	static_cast<int>(amberReader->getDihedralsPeriod(t+ 4)), // TODO wants int, returns double
						// 	amberReader->getDihedralsForceK(t+4),
						// 	static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+4)))
						// );
					}
					
					// Add the dihedral to the list of impropers.
					allDihedralsACIxs.push_back(thisDihedralACIxs);

				} // END of not foundit
		
			}
			
			if (improper){
			// If it is an improper dihedral, we check if it exitsts, without
			// checking for the reverse (order matters for impropers)

			bool foundit = false;
			for(auto& row:allImpropersACIxs){
				if ((thisDihedralACIxs == row))
					{
					foundit = true;	
					break;
					}
				}
			
			if (!foundit){ // improper was not found

				/* std::cout << "*** improper inserted: " 
					<< thisDihedralACIxs[0] << " " << thisDihedralACIxs[1] << " "
					<< thisDihedralACIxs[2] << " " << thisDihedralACIxs[3] << " := " 
					<< numberOf << " " << amberReader->getDihedralsPeriod(t) 
					<< " ;| " << amberReader->getDihedralsForceK(t) 
					<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)));
				if(numberOf == 2){std::cout
					<< " | " << amberReader->getDihedralsForceK(t+1)
					<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)));}
				if(numberOf == 3){std::cout
					<< " | " << amberReader->getDihedralsForceK(t+2)
					<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2)));}
				if(numberOf > 3){std::cout
					<< " | " << amberReader->getDihedralsForceK(t+3)
					<< " " << (ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+3)));}
				std::cout << std::endl << std::flush; */

				// Define the dihedrals
				if(numberOf == 1){
					dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
						static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)))
					);
				}
				else if(numberOf == 2){
					dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
						static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
						static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t + 1),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)))
					);
				}
				else if(numberOf == 3){
					dumm.defineAmberImproperTorsion_KA(aCIx1, aCIx2, aCIx3, aCIx4,
						static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
						static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t + 1),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))),
						static_cast<int>(amberReader->getDihedralsPeriod(t + 2)), // TODO wants int, returns double
						amberReader->getDihedralsForceK(t + 2),
						static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2)))
					);
				}
				// Add the improper to the list of impropers.
				allImpropersACIxs.push_back(thisDihedralACIxs);
			}
			}	


	}
}

/** Adds force field parameters read by the inputReader **/
void Topology::generateDummParams(
	readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs

)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Add types
	generateDummAtomClasses(resName, amberReader, dumm, aClassParams2aClassId);

	// Add parameters
	bAddDummBondParams(resName, amberReader, dumm, allBondsACIxs);
	bAddDummAngleParams(resName, amberReader, dumm, allAnglesACIxs);
	bAddDummTorsionParams(resName, amberReader, dumm, allDihedralsACIxs, allImpropersACIxs);
}

/** Transfer already generated force field parameters to DuMM **/
void Topology::transferDummParams(
	readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
	, std::map<AtomClassParams, AtomClassId>& aClassParams2aClassId
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allBondsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allAnglesACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allDihedralsACIxs
	, std::vector<std::vector<SimTK::DuMM::AtomClassIndex>>& allImpropersACIxs

)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Add types
	//transferDummAtomClasses(resName, amberReader, dumm, aClassParams2aClassId);
	transferDummChargedAtomClasses(resName, amberReader, dumm, aClassParams2aClassId);

	// Add parameters
	bAddDummBondParams(resName, amberReader, dumm, allBondsACIxs);
	bAddDummAngleParams(resName, amberReader, dumm, allAnglesACIxs);
	bAddDummTorsionParams(resName, amberReader, dumm, allDihedralsACIxs, allImpropersACIxs);
}

bool Topology::checkIfTripleUnorderedAreEqual(
		std::vector<SimTK::Compound::AtomIndex> &first,
		std::vector<SimTK::Compound::AtomIndex> &second)
{
	if(first == second){
		return true;
	}
	else if(
			(first[0] == second[2]) &&
			(first[1] == second[1]) &&
			(first[2] == second[0])
			){
		return true;
	}
	else{
		return false;
	}

}

// Helper function for calcLogDetMBATAnglesContribution
// Finds all triple runs - TODO VERY INEFFICIENT
void Topology::loadTriples()
{
	// Assign Compound coordinates by matching bAtomList coordinates
	//std::cout << "Topology triples: " << std::endl ;
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].getCompoundAtomIndex(), vec));
	}
	std::vector< std::vector<Compound::AtomIndex> > bondedAtomRuns =
	getBondedAtomRuns(3, atomTargets);

	// Find root bAtomList index
	// int rootIx;
	int ix = -1;
	for(const auto& atom: bAtomList){
		ix++;
		// if(atom.getCompoundAtomIndex() == 0){
		// 	rootIx = ix;
		// 	break;
		// }
	}

	// Find neighbour with maximum atomIndex
	int maxAIx = -1;
	Compound::AtomIndex aIx;
	for(auto atom: bAtomList[ix].neighbors){
		aIx = atom->getCompoundAtomIndex();
		if(aIx > maxAIx){
			maxAIx = aIx;
		}
	}

	//std::cout << "==========================" << std::endl;

	int flag;
	int bIx = -1;
	for(auto bAR: bondedAtomRuns){ // Iterate bondedAtomRuns
		//std::cout << "checking " ;
		//for(auto aIx: bAR){std::cout << " " << aIx;}; std::cout << std::endl;

		bIx++;
		flag = 0;
		for(auto tripleEntry: triples){ // Iterate triples gathered so far
			if(checkIfTripleUnorderedAreEqual(bAR, tripleEntry)){
				flag = 1;
				break;
			}
		} // END Iterate triples gathered so far

		if(!flag){ // Not found in gathered triples
			if((bAR[0] < bAR[1]) || (bAR[2] < bAR[1]) // Only level changing branches
			|| ((bAR[1] == 0) && (bAR[2] == maxAIx)) // except for the root atom
			){
				triples.push_back(bAR);
				//for(auto aIx: triples.back()){std::cout << " " << aIx;}; std::cout << std::endl;
			}
		}

	} // END Iterate bondedAtomRuns
}

// Numerically unstable around -pi, 0 and pi due to the log(0)
SimTK::Real Topology::calcLogSineSqrGamma2(const SimTK::State &quatState)
{
	bSpecificAtom *root = &(bAtomList[bSpecificAtomRootIndex]);
	SimTK::Compound::AtomIndex aIx = root->getCompoundAtomIndex();
	SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, aIx);
	SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();

	SimTK::Real w = quat[0];
	SimTK::Real x = quat[1];
	SimTK::Real y = quat[2];
	SimTK::Real z = quat[3];
	SimTK::Real sinPitch = 2 * ((w * y) - (z * x));

	//std::cout << std::setprecision(20) << std::fixed;
	SimTK::Real pitch = std::asin(sinPitch);
	std::cout << "Topology pitch " << pitch << std::endl;

	SimTK::Real result = std::log(sinPitch * sinPitch);

	// Quick and dirty
	if(result < -14.0){ // Around double precision log(0)
		result = -14.0;
	}

	// More elaborate fix
	//TODO: Given a difference compute the limit of this log

	return result;
}


SimTK::Real Topology::calcLogDetMBATGamma2Contribution(const SimTK::State& quatState){
	//State& eulerState;
	//matter.convertToEulerAngles(quatState, eulerState);
	//std::cout << "calcLogDetMBATGamma2Contribution quaternionState " << quatState << std::endl;
	//std::cout << "calcLogDetMBATGamma2Contribution	  eulerState " << eulerState << std::endl;

	bSpecificAtom *root = &(bAtomList[bSpecificAtomRootIndex]);
	SimTK::Compound::AtomIndex aIx = root->getCompoundAtomIndex();
	SimTK::Transform X = calcAtomFrameInGroundFrame(quatState, aIx);
	SimTK::Quaternion quat = (X.R()).convertRotationToQuaternion();
	//std::cout << "calcLogDetMBATGamma2Contribution quaternion " << quat << std::endl;

	SimTK::Real w = quat[0];
	SimTK::Real x = quat[1];
	SimTK::Real y = quat[2];
	SimTK::Real z = quat[3];
	SimTK::Real sinPitch = 2 * (w * y - z * x);

	std::cout << std::setprecision(20) << std::fixed;
	std::cout << "sinpitch " << sinPitch << std::endl;
	std::cout << "sinpitchsq" << sinPitch * sinPitch << std::endl;

	SimTK::Real pitch = std::asin(sinPitch);
	std::cout << "pitch " << pitch << std::endl;
	//if(pitch < 0){
	//	pitch = pitch + (2*SimTK_PI);
	//	std::cout << "sin converted pitch " << std::sin(pitch) << std::endl;
	//}

	if(sinPitch < SimTK::Eps){ // consider using SimTK::Eps
		return -SimTK::Infinity;
	}
	SimTK::Real result = std::log(sinPitch * sinPitch);
	return result;
}


SimTK::Real Topology::calcLogDetMBATDistsContribution(const SimTK::State&){
	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].getCompoundAtomIndex(), vec));
	}

	//std::cout << "Topology::calcLogDetMBATDistsContribution dists: " ;
	SimTK::Real result = 0.0;

	for(auto bond: bonds){
		//SimTK::Vec3 vec0 = calcAtomLocationInGroundFrame(someState, triple[0]);
		//SimTK::Vec3 vec1 = calcAtomLocationInGroundFrame(someState, triple[1]);
		//SimTK::Vec3 vec2 = calcAtomLocationInGroundFrame(someState, triple[2]);

		SimTK::Vec3 atom1pos = SimTK::Vec3(bAtomList[bond.i].getX(), bAtomList[bond.i].getY(), bAtomList[bond.i].getZ());
		SimTK::Vec3 atom2pos = SimTK::Vec3(bAtomList[bond.j].getX(), bAtomList[bond.j].getY(), bAtomList[bond.j].getZ());
		//SimTK::Real distSqr = (atom2pos - atom1pos).normSqr(); // funny results ?
		SimTK::Real dist = std::sqrt(
				std::pow(atom2pos[0] - atom1pos[0], 2) +
				std::pow(atom2pos[1] - atom1pos[1], 2) +
				std::pow(atom2pos[2] - atom1pos[2], 2));

		//std::cout << "atom " << bond.j << " 2*logDistSqr " << std::log(distSqr) + std::log(distSqr) << " " ;
		//std::cout << "dist " << bond.j << " = " << dist << " " ;

		result = result + (4.0 * std::log(dist));

	}
	//std::cout << std::endl;

	return result;
}


SimTK::Real Topology::calcLogDetMBATDistsMassesContribution(const SimTK::State&)
{

	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
			Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
			atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].getCompoundAtomIndex(), vec));
	}

	//std::cout << "Topology::calcLogDetMBATDistsMassesContribution dists squared: " ;
	SimTK::Real result = 3.0*std::log(bAtomList[bonds[0].i].mass);
	//std::cout << "atom " << bonds[0].i << " 1stlogmass " << std::log(bAtomList[bonds[0].i].mass) << " " ;
	for(const auto& bond: bonds){
		//SimTK::Vec3 vec0 = calcAtomLocationInGroundFrame(someState, triple[0]);
		//SimTK::Vec3 vec1 = calcAtomLocationInGroundFrame(someState, triple[1]);
		//SimTK::Vec3 vec2 = calcAtomLocationInGroundFrame(someState, triple[2]);

		SimTK::Vec3 atom1pos = SimTK::Vec3(bAtomList[bond.i].getX(), bAtomList[bond.i].getY(), bAtomList[bond.i].getZ());
		SimTK::Vec3 atom2pos = SimTK::Vec3(bAtomList[bond.j].getX(), bAtomList[bond.j].getY(), bAtomList[bond.j].getZ());
		SimTK::Real distSqr = (atom2pos - atom1pos).normSqr();

		//std::cout << "atom " << bond.j << " 2*logDistSqr+mass " << std::log(distSqr) + std::log(distSqr) + std::log(bAtomList[bond.j].mass) << " " ;

		result = result + std::log(distSqr) + std::log(distSqr) + (3.0*std::log(bAtomList[bond.j].mass));

	}
	//std::cout << std::endl;

	return result;
}

/* Calculate angle contribution at the MBAT determinant:
TODO: calculate atomTargets only once: getAtomLocationsInGround
TODO: realize position
*/
SimTK::Real Topology::calcLogDetMBATAnglesContribution(const SimTK::State&){
	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
			Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
			atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].getCompoundAtomIndex(), vec));
	}

	//std::cout << "Topology::calcLogDetMBATAnglesContribution angles: " ;
	SimTK::Real result = 0.0;
	for(auto triple: triples){
		//SimTK::Vec3 vec0 = calcAtomLocationInGroundFrame(someState, triple[0]);
		//SimTK::Vec3 vec1 = calcAtomLocationInGroundFrame(someState, triple[1]);
		//SimTK::Vec3 vec2 = calcAtomLocationInGroundFrame(someState, triple[2]);

		SimTK::UnitVec3 v1(atomTargets.find(triple[0])->second - atomTargets.find(triple[1])->second);
		SimTK::UnitVec3 v2(atomTargets.find(triple[2])->second - atomTargets.find(triple[1])->second);

		SimTK::Real dotProduct = dot(v1, v2);
		assert(dotProduct < 1.1);
		assert(dotProduct > -1.1);
		if (dotProduct > 1.0) dotProduct = 1.0;
		if (dotProduct < -1.0) dotProduct = -1.0;
		// SimTK::Real angle = std::acos(dotProduct);
		//std::cout << SimTK_RADIAN_TO_DEGREE * angle << " " ;

		SimTK::Real sinSqAngle = 1 - (dotProduct * dotProduct);

		result = result + std::log(sinSqAngle);

	}
	//std::cout << std::endl;

	return result;
}


SimTK::Real Topology::calcLogDetMBATMassesContribution(const SimTK::State&)
{
	// function args were const SimTK::State& someState

	//std::cout << "Topology::calcLogDetMBATMassesContribution masses: " ;
	SimTK::Real result = 0.0;
	for(const auto& atom: bAtomList){
		//std::cout << 3.0 * std::log(atom.mass) << " " ;
		//std::cout << atom.mass << " " ;
		result += 3.0 * std::log(atom.mass);
	}
	//std::cout << std::endl;

	return result;
}

SimTK::Real Topology::calcLogDetMBAT(const SimTK::State& someState)
{
	SimTK::Real gamma2Contribution = calcLogDetMBATGamma2Contribution(someState);
	SimTK::Real distsContribution = calcLogDetMBATDistsContribution(someState);
	//SimTK::Real distsMassesContribution = calcLogDetMBATDistsMassesContribution(someState);
	SimTK::Real anglesContribution = calcLogDetMBATAnglesContribution(someState);
	SimTK::Real massesContribution = calcLogDetMBATMassesContribution(someState);

	//std::cout << std::setprecision(20) << std::fixed;
	//std::cout << " gamma dists masses angles contributions: "
	//	<< gamma2Contribution << " "
	//	<< distsContribution << " "
	//	<< massesContribution << " "
	//	//<< distsMassesContribution << " "
	//	<< anglesContribution << std::endl;

	return gamma2Contribution + distsContribution + anglesContribution + massesContribution;
}

SimTK::Real Topology::calcLogDetMBATInternal(const SimTK::State& someState)
{
	SimTK::Real distsContribution = calcLogDetMBATDistsContribution(someState);
	SimTK::Real anglesContribution = calcLogDetMBATAnglesContribution(someState);
	SimTK::Real massesContribution = calcLogDetMBATMassesContribution(someState);

	//std::cout << std::setprecision(20) << std::fixed;
	//std::cout << "dists masses angles contributions: "
	//          << distsContribution << " "
	//          << massesContribution << " "
	//          << anglesContribution << std::endl;

	return distsContribution + anglesContribution + massesContribution;
}


/**
 ** Interface **
 **/

/** Get the number of atoms in the molecule **/
int Topology::getNAtoms() const{
	return getNumAtoms();
}

/** Get the number of bonds in the molecule **/
int Topology::getNBonds() const{
	assert(!"Not implemented."); throw std::exception();

	return std::numeric_limits<int>::max();
}

/** Get a pointer to an atom object in the atom list inquiring
by number **/
bSpecificAtom * Topology::getAtomByNumber(int) const {
	// function args were int number
	assert(!"Not implemented."); throw std::exception();

	return nullptr;
}

/** Get a pointer to an atom object in the atom list inquiring
by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
// TODO: Optimize use CompoundAtomIx2GmolAtomIx instead
bSpecificAtom * Topology::updAtomByAtomIx(int aIx) {
	for (int i = 0; i < natoms; i++){
		if(bAtomList[i].getCompoundAtomIndex() == aIx){
			return &bAtomList[i];
		}
	}

	return nullptr;
}

/** Get a pointer to an atom object in the atom list inquiring
by atom name. **/
bSpecificAtom * Topology::getAtomByName(std::string) const {
	// function args were std::string name
	assert(!"Not implemented."); throw std::exception();

	return nullptr;
}

/** Get the neighbours in the graph. **/
std::vector<bSpecificAtom *> Topology::getNeighbours(int) const {
	assert(!"Not implemented."); throw std::exception();

	return {};
}

/* Check if a1 and a2 are bonded */
bool Topology::checkBond(int a1, int a2)
{
	for(int i = 0; i < nbonds; i++){
		if( (bonds[i]).isThisMe(a1, a2) )
		{
			return true;
		}
	}
	return false;
}

/** **/
const bBond& Topology::getBond(int a1, int a2) const
{
	// TODO is this fast enough?
	const auto bond = std::find_if(bonds.begin(), bonds.end(), [a1, a2](bBond b) { return b.isThisMe(a1, a2); });

#ifndef NDEBUG
    std::string assert_string("No bond with these atom indeces found: " + to_string(a1) + " " + to_string(a2));
	assert(bond != bonds.end() && assert_string.size());
#endif

	return *bond;
}

/** Get bond order. **/
int Topology::getBondOrder(int, int) const
{
	assert(!"Not implemented.");

	return 0;
}



/** Get regimen **/
std::string Topology::getRegimen(){
	return this->regimen;
}

/** Set scale factors for U entries according to flexibility file **/
void Topology::setUScaleFactorsToBonds(std::string flexFN)
{
	//
	std::string line;
	std::ifstream F(flexFN);

	while(F.good()){
		std::getline(F, line);
		//std::cout << "setUScaleFactors line " << line << std::endl;
		if(!line.empty()){
			if(line.at(0) == '#'){
				continue;
			}

			std::istringstream iss(line);
			std::string word;
			std::vector<std::string> lineWords;

			while(iss >> word){
				lineWords.push_back(std::move(word));
			}
			if(lineWords.size() >= 4 ){
				// if ( (lineWords[3]).find_first_not_of (".+-0123456789") == ..... ) // slow variant
				for(int i = 0; i < nbonds; i++){
					if(bonds[i].isThisMe(std::stoi(lineWords[0]), std::stoi(lineWords[1])) ){
						if(lineWords[3][0] != '#'){
							bonds[i].addUScaleFactor(std::stof(lineWords[3]));
						}else{
							bonds[i].addUScaleFactor(1.0);
						}
						break;
					}
				}
			}
		}
	}

}


/** Set regimen according to input file **/
// Compound doesn't care about Worlds. We keep multiple Bond::Mobilities in bBond
void Topology::setFlexibility(std::string argRegimen, std::string flexFN, int whichWorld)
{

	if(argRegimen == "IC"){
		for (unsigned int r=0 ; r<getNumBonds(); r++){
			setBondMobility(BondMobility::Free, Compound::BondIndex(r));
			bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
					BondMobility::Free); //OLDMOB
			/*setBondMobility(BondMobility::Translation, Compound::BondIndex(r));
			bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
					BondMobility::Translation); //NEWMOB*/
		}
	}else if(argRegimen == "TD") {
		for (unsigned int r = 0; r < getNumBonds(); r++) {
			setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
			bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
					BondMobility::Torsion);
		}
	}else if(argRegimen == "BA"){
		for (unsigned int r=0 ; r<getNumBonds(); r++){

			int  firstSpecificAtomIndex = bonds[bondIx2GmolBond[Compound::BondIndex(r)]].i;
			int secondSpecificAtomIndex = bonds[bondIx2GmolBond[Compound::BondIndex(r)]].j;

			if((std::string(bAtomList[ firstSpecificAtomIndex].getFftype()) == "CA") ||
			   (std::string(bAtomList[secondSpecificAtomIndex].getFftype() )== "CA")){
				setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
				bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
						BondMobility::Torsion);
			} else if( (bAtomList[ firstSpecificAtomIndex].getNBonds() < 3) ||
				(bAtomList[secondSpecificAtomIndex].getNBonds() < 3)){
				setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
				bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
						BondMobility::Torsion);
			}else{
				setBondMobility(BondMobility::BallF, Compound::BondIndex(r));
				bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
						BondMobility::BallF);
			}

			if(bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].isRingClosing()){
				bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].addBondMobility(
						BondMobility::Rigid);
			}

		}

	//}else if((argRegimen.at(0) == 'R') || (argRegimen.at(0) == 'S')){
	}else{

		// Set all Compound and Topology bonds to rigid
		for (unsigned int r=0 ; r<getNumBonds(); r++){
			setBondMobility(BondMobility::Rigid, SimTK::Compound::BondIndex(r));
		}
		for (unsigned int r=0 ; r<getNumBonds(); r++){
			bonds[r].addBondMobility(BondMobility::Rigid); // TODO: Change to rigid
		}

		// Get flexible bonds from file. Numbering starts at 0 in prmtop
		std::string line;
		std::ifstream F(flexFN);

		//printMaps();
		//std::cout << "DEBUG GmolBond2bondIx " << std::endl << std::flush;
		//std::cout << "GmolBond2bondIx size " << GmolBond2bondIx.size() << std::endl;
		//std::cout << "GmolBond2bondIx:" << std::endl;
		//for(unsigned int i = 0; i < nbonds; i++){
		//	std::cout << i << ' ' << GmolBond2bondIx.at(i) << std::endl;
		//}

		while(F.good()){
			std::getline(F, line);
			if(!line.empty()){
				if(line.at(0) == '#'){
					continue;
				}

				std::istringstream iss(line);
				std::string word;
				std::vector<std::string> lineWords;

				while(iss >> word){
					lineWords.push_back(std::move(word));
				}
				if(lineWords.size() >= 2 ){ // TODO: Check if it should be 3
					for(int i = 0; i < nbonds; i++){

						if(bonds[i].isThisMe(
						    std::stoi(lineWords[0]), std::stoi(lineWords[1])) ){

						    //std::cout << "DEBUG bond about to be set " << i << " ";

						    if(lineWords.size() == 2) {
						        bonds[i].setBondMobility(BondMobility::Torsion, whichWorld);
						        setBondMobility(BondMobility::Torsion,
						                        GmolBond2bondIx.at(i));
						        break;
						    }else{

								//std::cout << " lineWords[2] " << lineWords[2] << " ";

						        if(lineWords[2] == "Free"){
						            bonds[i].setBondMobility(BondMobility::Free, whichWorld);
						            setBondMobility(BondMobility::Free,
						                            GmolBond2bondIx.at(i));
						        }else if(lineWords[2] == "FreeLine") {
						            bonds[i].setBondMobility(BondMobility::FreeLine, whichWorld);
						            setBondMobility(BondMobility::FreeLine,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "LineOrientationF") {
						            bonds[i].setBondMobility(BondMobility::LineOrientationF, whichWorld);
						            setBondMobility(BondMobility::LineOrientationF,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "LineOrientationM") {
						            bonds[i].setBondMobility(BondMobility::LineOrientationM, whichWorld);
						            setBondMobility(BondMobility::LineOrientationM,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if((lineWords[2] == "Translation") || (lineWords[2] == "Cartesian")) {
						            bonds[i].setBondMobility(BondMobility::Translation, whichWorld);
						            setBondMobility(BondMobility::Translation,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if((lineWords[2] == "Pin") || (lineWords[2] == "Torsion")) {
						            bonds[i].setBondMobility(BondMobility::Torsion, whichWorld);
						            setBondMobility(BondMobility::Torsion,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "AnglePin") {
						            bonds[i].setBondMobility(BondMobility::AnglePin, whichWorld);
						            setBondMobility(BondMobility::AnglePin,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "Slider") {
						            bonds[i].setBondMobility(BondMobility::Slider, whichWorld);
						            setBondMobility(BondMobility::Slider,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "Cylinder") {
						            bonds[i].setBondMobility(BondMobility::Cylinder, whichWorld);
						            setBondMobility(BondMobility::Cylinder,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "BendStretch") {
						            bonds[i].setBondMobility(BondMobility::BendStretch, whichWorld);
						            setBondMobility(BondMobility::BendStretch,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "UniversalM") {
						            bonds[i].setBondMobility(BondMobility::UniversalM, whichWorld);
						            setBondMobility(BondMobility::UniversalM,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "Spherical") {
						            bonds[i].setBondMobility(BondMobility::Spherical, whichWorld);
						            setBondMobility(BondMobility::Spherical,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "BallF") {
						            bonds[i].setBondMobility(BondMobility::BallF, whichWorld);
						            setBondMobility(BondMobility::BallF,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if(lineWords[2] == "BallM") {
						            bonds[i].setBondMobility(BondMobility::BallM, whichWorld);
						            setBondMobility(BondMobility::BallM,
						                            GmolBond2bondIx.at(i));
						            break;
						        }else if((lineWords[2] == "Rigid") || (lineWords[2] == "Weld")){
						            bonds[i].setBondMobility(BondMobility::Rigid, whichWorld);
						            setBondMobility(BondMobility::Rigid,
						                            GmolBond2bondIx.at(i));
						        }else{
						            bonds[i].setBondMobility(BondMobility::Torsion, whichWorld);
						            setBondMobility(BondMobility::Torsion,
						                            GmolBond2bondIx.at(i));
						            break;
						        }


						    } // if nof words is different than 2

						} // if bond found

					} // end iterate bonds

				} // if more than 2 words

			} // if line not empty

		} // end read file

/*		std::cout << "Assigned mobilities:" << std::endl;
		for(unsigned int i = 0; i < nbonds; i++){
			std::cout << i << ' ' << GmolBond2bondIx.at(i) << " " << bonds[i].getBondMobility() << std::endl;
		}*/


	} // RB
	//PrintAtomList(0);
	this->regimen = argRegimen;
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void Topology::loadAIx2MbxMap()
{

	// If the map is empty fill with empty vectors first
	if(aIx2mbx.empty()){
		// Iterate through atoms and get their MobilizedBodyIndeces
		for (unsigned int i = 0; i < getNumAtoms(); ++i) {

			// Get atomIndex from atomList
			SimTK::Compound::AtomIndex aIx = (bAtomList[i]).getCompoundAtomIndex();

			// Insert
			aIx2mbx.insert(
				std::pair< SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex> >
					(aIx, std::vector<SimTK::MobilizedBodyIndex>())
			);
		}

	}

	// Iterate through atoms and get their MobilizedBodyIndeces
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {

		// Get atomIndex from atomList
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).getCompoundAtomIndex();

		// Get MobilizedBodyIndex from CompoundAtom
		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);

		// Insert
		//aIx2mbx.insert(
		//		std::pair<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>
		//		(aIx, mbx));
		aIx2mbx[aIx].emplace_back(mbx);
	}


}


/** Compound AtomIndex to bAtomList number **/
void Topology::loadCompoundAtomIx2GmolAtomIx()
{
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).getCompoundAtomIndex();
		int gmolIx = (bAtomList[i]).getNumber();

		CompoundAtomIx2GmolAtomIx.insert(
			std::pair<SimTK::Compound::AtomIndex, int>
			(aIx, gmolIx));
	}
}

/**  **/
int Topology::getNumber(SimTK::Compound::AtomIndex aIx)
{
	return CompoundAtomIx2GmolAtomIx[aIx];
}


/** Calculate all atom frames in top frame. It avoids calling
calcDefaultAtomFrameInCompoundFrame multiple times. This has
to be called every time the coordinates change though. **/
void Topology::calcTopTransforms()
{
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).getCompoundAtomIndex();
		aIx2TopTransform[aIx] = calcDefaultAtomFrameInCompoundFrame(aIx);
	}
}

/**  **/
void Topology::printTopTransforms()
{
	std::cout << "Topology TopTransforms " << std::endl;
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).getCompoundAtomIndex();
		std::cout << aIx << " " << aIx2TopTransform[aIx] << std::endl;
	}
}

SimTK::Transform Topology::getTopTransform(SimTK::Compound::AtomIndex aIx)
{
	return aIx2TopTransform[aIx];
}

SimTK::Transform Topology::calcDefaultAtomFrameInCompoundFrameThroughDuMM(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm,
	SimTK::SimbodyMatterSubsystem& matter,
	const SimTK::State& someState)
{
	// TODO
	assert(!"Not implemented!");
	return {};

	// SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	// const  SimTK::MobilizedBodyIndex mbx = dumm.getAtomBody(dAIx);
	// const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

	// // Body transforms
	// const Transform&    X_GB = mobod.getBodyTransform(someState);
	// const Rotation&     R_GB = X_GB.R();
	// const Vec3&         p_GB = X_GB.p();

	// // Atom
	// SimTK::Vec3 station = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);

	// // return
	// //const Vec3 p_BS_G = R_GB * station;
	// //return X;

}

// Return mbx by calling DuMM functions
SimTK::MobilizedBodyIndex Topology::getAtomMobilizedBodyIndexThroughDumm(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{
	SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	return dumm.getAtomBody(dAIx);
}

// Get atom location on mobod through DuMM functions
SimTK::Vec3 Topology::getAtomLocationInMobilizedBodyFrameThroughDumm(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{
	SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	return dumm.getAtomStationOnBody(dAIx);
}

SimTK::Vec3 Topology::calcAtomLocationInGroundFrameThroughSimbody(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm,
	SimTK::SimbodyMatterSubsystem& matter,
	const SimTK::State& someState)
{
	const SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm);
	const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

	const Transform&    X_GB = mobod.getBodyTransform(someState);
	const Rotation&     R_GB = X_GB.R();
	const Vec3&         p_GB = X_GB.p();

	SimTK::Vec3 station = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);

	const Vec3 p_BS_G = R_GB * station;
	return p_GB + p_BS_G;

}

// Retunr mbx from an olresdy saved map inside Topology
SimTK::MobilizedBodyIndex Topology::getAtomMobilizedBodyIndexFromMap(
	SimTK::Compound::AtomIndex aIx, int whichWorld)
{
	if(!aIx2mbx.empty()){
		if(!((aIx2mbx[aIx]).empty())){
			return (aIx2mbx[aIx])[whichWorld];
		}else{
			std::cerr << "Topology::getAtomMobilizedBodyIndexFromMap: aIx2mbx for "
				<< aIx <<  "atom not yet loaded.\n";
			throw std::exception();
			std::exit(1);
		}
	}else{
		std::cerr << "Topology::getAtomMobilizedBodyIndexFromMap: aIx2mbx not yet loaded.\n";
		throw std::exception();
		std::exit(1);
	}
}

/** Print maps **/
void Topology::printMaps()
{
/*
	std::cout << "Topology " << name << " maps " << std::endl;
	std::cout << "mbx2aIx:" << std::endl;
	map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>::const_iterator mbx2aIxIt;
	for(mbx2aIxIt = mbx2aIx.begin();
	   mbx2aIxIt != mbx2aIx.end(); ++mbx2aIxIt)
	{
		std::cout << "mbx " << mbx2aIxIt->first
			<< " atomIndex " << mbx2aIxIt->second << std::endl;
	}
*/
	std::cout << "Topology map aIx2mbx:" << std::endl;
	map<SimTK::Compound::AtomIndex, std::vector<SimTK::MobilizedBodyIndex>>::const_iterator aIx2mbxIt;
	for(aIx2mbxIt = aIx2mbx.begin();
	   aIx2mbxIt != aIx2mbx.end(); ++aIx2mbxIt)
	{
		std::cout << "atomIndex " << aIx2mbxIt->first << " mbxs:";
		for (auto val : (aIx2mbxIt->second)){
			std::cout << " " << val ;
		}
		std::cout << std::endl << std::flush;
	}

	std::cout << "Topology map CompoundAtomIx2GmolAtomIx:" << std::endl;
	std::map< SimTK::Compound::AtomIndex, int >::const_iterator aIx2gmolaIxIt;
	for(aIx2gmolaIxIt = CompoundAtomIx2GmolAtomIx.begin();
	   aIx2gmolaIxIt != CompoundAtomIx2GmolAtomIx.end(); ++aIx2gmolaIxIt)
	{
		std::cout << "atomIndex " << aIx2gmolaIxIt->first
			<< " gmolaIx " << aIx2gmolaIxIt->second
			<< std::endl << std::flush;
	}
/*
	map<SimTK::Compound::BondIndex, int>::const_iterator bondIx2GmolBondIt;
	for(bondIx2GmolBondIt = bondIx2GmolBond.begin();
	   bondIx2GmolBondIt != bondIx2GmolBond.end(); ++bondIx2GmolBondIt)
	{
		std::cout << "Compound bondIndex " << bondIx2GmolBondIt->first
			<< " bBond index " << bondIx2GmolBondIt->second << std::endl;
	}

	map<int, SimTK::Compound::BondIndex>::const_iterator GmolBond2bondIxIt;
	for(GmolBond2bondIxIt = GmolBond2bondIx.begin();
		GmolBond2bondIxIt != GmolBond2bondIx.end(); ++GmolBond2bondIxIt)
	{
		std::cout << "bBond index " << GmolBond2bondIxIt->first
			<< " Compound index " << GmolBond2bondIxIt->second << std::endl;
	}
*/
}

/** Write a pdb with bAtomList coordinates and inNames **/
void Topology::writeAtomListPdb(std::string dirname,
	std::string prefix,
	std::string sufix,
	int maxNofDigits,
	int index) const
{
	// Using floor here is no buneo because the index can be zero
	std::string zeros("");
	int nofDigits = static_cast<int>(std::to_string(index).size());
	if(maxNofDigits > nofDigits){
		zeros = std::string(maxNofDigits - nofDigits, '0');
	}

	std::stringstream sstream;
	sstream << dirname << "/"
		<< prefix << zeros << std::to_string(index) << sufix;
	string ofilename = sstream.str();

	FILE *oF = fopen (ofilename.c_str(),"w");
	if (oF) {
		// Pdb lines
		for(int i = 0; i < getNumAtoms(); i++){
			fprintf(oF, "%-6s%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  %4.2f%6.2f          %2s\n"
				, "ATOM"                 // record
				, i                      // index
				, bAtomList[i].getInName().c_str()  // name
				, "UNK"                  // residue name
				, 'A'                    // chain
				, 1                      // residue index
				, 10.0*bAtomList[i].getX()    // x in A
				, 10.0*bAtomList[i].getY()    // y in A
				, 10.0*bAtomList[i].getZ()    // z in A
				, 1.0                    // occupancy
				, 0.0                    // beta factor
				, "  ");                 // element
		}

		fclose(oF);
		std::cout << "\tTopology written to '" << ofilename << "'\n";
	} else {
		std::cout << "FAILED TO OPEN '" << ofilename << "' TO WRIE!\n";
	}
}

/** Get bAtomList coordinates coordinates **/
void Topology::getCoordinates(
		std::vector<SimTK::Real> Xs,
		std::vector<SimTK::Real> Ys,
		std::vector<SimTK::Real> Zs)
{
	assert(Xs.size() == static_cast<size_t>(getNumAtoms()));
	assert(Ys.size() == static_cast<size_t>(getNumAtoms()));
	assert(Zs.size() == static_cast<size_t>(getNumAtoms()));
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Xs[ix] = bAtomList[ix].getX();
		Ys[ix] = bAtomList[ix].getY();
		Zs[ix] = bAtomList[ix].getZ();
	}
}

/** Get own CompoundIndex in CompoundSystem **/
const CompoundSystem::CompoundIndex &Topology::getCompoundIndex() const
{
	return compoundIndex;
}

/** Set the compoundIndex which is the position in the vector of Compounds
 * of the CompoundSystem **/
void Topology::setCompoundIndex(
	const CompoundSystem::CompoundIndex &compoundIndex)
{
	//Topology::compoundIndex = compoundIndex;
	this->compoundIndex = compoundIndex;
}


/** Get the neighbor atom in the parent mobilized body.
TODO: No chemical parent for satelite atoms or first atom. **/
SimTK::Compound::AtomIndex
Topology::getChemicalParent(
	SimTK::SimbodyMatterSubsystem *matter,
	//std::unique_ptr<SimTK::SimbodyMatterSubsystem> matter,
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm)
{

	SimTK::Compound::AtomIndex chemParentAIx;
	int gmolAtomIndex = -111111;

	//if(getAtomLocationInMobilizedBodyFrame(aIx) == 0){ // atom is at body's origin // SAFE
	if(getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm) == 0){ // atom is at body's origin // DANGER

		// Get body, parentBody, parentAtom
		//SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx); // SAFE
		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndexThroughDumm(aIx, dumm); // DANGER
		const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
		const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
		SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

		if(parentMobod.getMobilizedBodyIndex() != 0){ // parent not Ground
			// Find the true bSpecificAtom (CHEMICAL) parent
			bSpecificAtom *originSpecAtom = nullptr;
			originSpecAtom = updAtomByAtomIx(aIx); //TODO: optimize

			// TODO: Check is neighbors and bondsInvolved are redundant
			// Loop through neighbor atoms (bSpecificAtom)
			for(unsigned int k = 0; k < (originSpecAtom->neighbors).size(); k++) {

				// Loop through bonds that this atom is involved in (bBond);
				for(std::vector<bBond *>::iterator bondsInvolvedIter = (originSpecAtom->bondsInvolved).begin();
					bondsInvolvedIter != (originSpecAtom->bondsInvolved).end(); ++bondsInvolvedIter){

					// Check if this neighbor is involved in this bond
					if( (*bondsInvolvedIter)->isThisMe(
						originSpecAtom->getNumber(), originSpecAtom->neighbors[k]->getNumber()
						) ){

						Compound::AtomIndex candidateChemParentAIx = originSpecAtom->neighbors[k]->getCompoundAtomIndex();

						// Check if neighbor atom's mobod is a parent mobod
						//if(getAtomMobilizedBodyIndex(candidateChemParentAIx) == parentMbx){ // SAFE
						if(getAtomMobilizedBodyIndexThroughDumm(candidateChemParentAIx, dumm) == parentMbx){ // DANGER

							if(!(*bondsInvolvedIter)->isRingClosing()){ // No ring atoms are allowed
								chemParentAIx = candidateChemParentAIx;
								gmolAtomIndex = originSpecAtom->neighbors[k]->getNumber();
								break;
							}
						}
					}
				}
			}
		}
		
	}else{
		std::cout << "Warning: requiring chemical parent for non-root atom\n";
		bSpecificAtom *originSpecAtom = nullptr;
		originSpecAtom = updAtomByAtomIx(aIx); //TODO: optimize
		for(unsigned int k = 0; k < (originSpecAtom->neighbors).size(); k++) {
			Compound::AtomIndex candidateChemParentAIx = originSpecAtom->neighbors[k]->getCompoundAtomIndex();
			//if(getAtomLocationInMobilizedBodyFrame(candidateChemParentAIx) == 0){ // atom is at body's origin // SAFE
			if(getAtomLocationInMobilizedBodyFrameThroughDumm(candidateChemParentAIx, dumm) == 0){ // atom is at body's origin // DANGER
				chemParentAIx = candidateChemParentAIx;
				gmolAtomIndex = originSpecAtom->neighbors[k]->getNumber();
				break;
			}
		}
	}

	return chemParentAIx;
}



