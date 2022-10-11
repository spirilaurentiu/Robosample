/**@file
Implementation of Topology class. **/

#include "Topology.hpp"

using namespace SimTK;

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
Topology::~Topology(){
	for(size_t i = 0; i < bAtomList.size(); i++){
		delete bAtomList[i].bAtomType;
	}
}

/** Set atoms properties from a reader: number, name, element, initial
 * name, force field type, charge, coordinates, mass, LJ parameters **/
void Topology::SetGmolAtomPropertiesFromReader(readAmberInput *amberReader)
{
	// Alloc memory for atoms and bonds list
	natoms = amberReader->getNumberAtoms();
	bAtomList.resize(natoms);

	// Initialize atom variables
	for(int i = 0; i < natoms; i++){
		bAtomList[i].Zero();
	}

	// Declare handy variables
	std::string str_buf;

	// Iterate through atoms and set as much as possible from amberReader
	for(int i = 0; i < natoms; i++){

		// Assign an index like in prmtop
		bAtomList[i].setNumber(i);

		// Assign element from the first letter of the name
		str_buf = amberReader->getAtomsName(i);
		unsigned int strix;
		for (strix = 0; strix < str_buf.length(); strix++){
			if(str_buf.at(strix) != ' '){
				break;
			}
		}
		bAtomList[i].setElem(str_buf.at(strix));

		// Assign a "unique" name. The generator is however limited.
		bAtomList[i].setName(GetUniqueName(i));

		// Store the initial name from prmtop
		bAtomList[i].setInName(str_buf);

		// Set atom type
		str_buf = amberReader->getAtomsNameAlias(i);
		ROBOSAMPLE::trim(str_buf);
		bAtomList[i].setFftype(str_buf);

		// Set charge as it is used in Amber
		SimTK::Real chargeMultiplier = 18.2223;
		bAtomList[i].setCharge(amberReader->getAtomsCharge(i) / chargeMultiplier);

		// Set coordinates in nm
		bAtomList[i].setX(amberReader->getAtomsXcoord(i) / 10.0);
		bAtomList[i].setY(amberReader->getAtomsYcoord(i) / 10.0);
		bAtomList[i].setZ(amberReader->getAtomsZcoord(i) / 10.0);

		// Set mass
		bAtomList[i].setMass(amberReader->getAtomsMass(i));

		// Set Lennard-Jones parameters
		bAtomList[i].setVdwRadius(amberReader->getAtomsRVdW(i));
		bAtomList[i].setLJWellDepth(amberReader->getAtomsEpsilon(i));

		// Set residue name and index
		bAtomList[i].residueName = std::string("UNK");
		bAtomList[i].residueIndex = 1;

	} // END atom properties
}

/** Set bonds properties from reader: bond indeces, atom neighbours **/
void Topology::SetGmolBondingPropertiesFromReader(readAmberInput *amberReader)
{
	assert( (!bAtomList.empty()) &&
	"Topology::loadAtomAndBondInfoFromReader: atom list empty.");

	// Alloc memory for bonds list
	nbonds = amberReader->getNumberBonds();
	bonds.resize(nbonds);

	// Iterate through bonds and get atom indeces
	for(int i=0; i<nbonds; i++){
		bonds[i].setIndex(i);
		bonds[i].i = amberReader->getBondsAtomsIndex1(i);
		bonds[i].j = amberReader->getBondsAtomsIndex2(i);
	}

	// Assign the number of bonds an atom has and set the number of freebonds
	// equal to the number of bonds for now
	for(int i = 0; i < natoms ; i++){
		bAtomList[i].nbonds = 0;
		for(int j = 0; j < nbonds; j++){
			if((bAtomList[i].number == bonds[j].i) || \
			   (bAtomList[i].number == bonds[j].j)){
				++bAtomList[i].nbonds;
				++bAtomList[i].freebonds;
			}
		}
	}

	// Assign neighbors and bonds involved for each atom
	// which translates into pushing bSpecificAtom * and bBond *
	// into their apropriate vectors
	for(int i=0; i<nbonds; i++){
		(bAtomList[ bonds[i].i  ]).addNeighbor( &(bAtomList[ bonds[i].j  ]) );
		(bAtomList[ bonds[i].i  ]).addBond( &(bonds[i]) );

		(bAtomList[ bonds[i].j  ]).addNeighbor( &(bAtomList[ bonds[i].i  ]) );
		(bAtomList[ bonds[i].j  ]).addBond( &(bonds[i]) );
	}
}

/** Set atoms Molmodel types (Compound::SingleAtom derived) based on
 * their valence **/
void Topology::SetGmolAtomsMolmodelTypesTrial(){

	// Set Gmolmodel name and element and inboard length
	for(int i = 0; i < (natoms); i++) {
		if(toupper(bAtomList[i].elem) == 'H'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom( bAtomList[i].name, Element( 1, "Hydrogen", "H", bAtomList[i].getMass()) );
			bAtomList[i].setAtomicNumber(1);
		}else if(toupper(bAtomList[i].elem) == 'O'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name, Element(8, "Oxygen", "O", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(8);
		}else if(toupper(bAtomList[i].elem) == 'F'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name, Element(9, "Fluorine", "F", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(9);
		}else if(toupper(bAtomList[i].elem) == 'I'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name, Element(53, "Iodine", "I", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(53);
		}else if(toupper(bAtomList[i].elem) == 'N'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name, Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(7);
		}else if(toupper(bAtomList[i].elem) == 'C'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name,  Element(6, "Carbon", "C", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(6);
		}else if(toupper(bAtomList[i].elem) == 'S'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(16);
		}else if(toupper(bAtomList[i].elem) == 'P'){
			bAtomList[i].bAtomType = new
					Compound::SingleAtom(bAtomList[i].name,  Element(15, "Phosphorus", "P", bAtomList[i].getMass()));
			bAtomList[i].setAtomicNumber(15);
		}

	}

	Angle TetrahedralAngle = 109.47 * Deg2Rad;
	for(int i = 0; i < (natoms); i++) {
		(bAtomList[i].bAtomType)->setCompoundName("SingleAtom");
		// Add BondCenters
		if(bAtomList[i].nbonds == 1){
			(bAtomList[i].bAtomType)->addFirstBondCenter("bond1", bAtomList[i].name);
		}else{
			(bAtomList[i].bAtomType)->addFirstTwoBondCenters("bond1", "bond2", bAtomList[i].name, UnitVec3(1, 0, 0), UnitVec3(-0.5, 0.866025, 0.0)); // FEBREM
			//(bAtomList[i].bAtomType)->addFirstTwoBondCenters("bond1", "bond2", bAtomList[i].name, UnitVec3(0.866025, 0.0, -0.5), UnitVec3(0.0, 0.0, 1.0)); // REFEBREM
			if(bAtomList[i].nbonds > 2){
				(bAtomList[i].bAtomType)->addLeftHandedBondCenter("bond3", bAtomList[i].name, TetrahedralAngle, TetrahedralAngle);
			}
			if(bAtomList[i].nbonds > 3){
				(bAtomList[i].bAtomType)->addRightHandedBondCenter("bond4", bAtomList[i].name, TetrahedralAngle, TetrahedralAngle);
			}
		}

		/*if(bAtomList[i].nbonds > 1){
			//(bAtomList[i].bAtomType)->addLeftHandedBondCenter("bond2", bAtomList[i].name, TetrahedralAngle, 0);
			(bAtomList[i].bAtomType)->addSecondBondCenter("bond2", bAtomList[i].name, TetrahedralAngle);
		}*/

		// Set the inboard BondCenter
		(bAtomList[i].bAtomType)->setInboardBondCenter("bond1");
		(bAtomList[i].bAtomType)->setDefaultInboardBondLength(0.19);
	}
}

void Topology::SetGmolAtomsMolmodelTypes()
{
	// ---------------------------------------------
	// Set every atom's (SimTK::Compound::SingleAtom *) to it's
	// appropriate element and assign it's Compound::AtomName to unique name
	// Every atom is derived from SingleAtom in turn derived from
	// Compound with one atom (AtomIndex 0)
	// Also set atom forecfield type
	// TODO: Bromine and Clorine and others
	// ---------------------------------------------
	for(int i = 0; i < (natoms); i++){
		// Atoms with one bond
		if(bAtomList[i].nbonds == 1){
			if(toupper(bAtomList[i].elem) == 'H'){
				bAtomList[i].bAtomType = new UnivalentAtom(bAtomList[i].name,
						                                   SimTK::Element( 1, "Hydrogen", "H", bAtomList[i].getMass() ));
				bAtomList[i].setAtomicNumber(1);
			}
				/*else if((toupper(bAtomList[i].name[0]) == 'C') && (toupper(bAtomList[i].name[0]) == 'L')){
					bAtomList[i].bAtomType = new
					UnivalentAtom(bAtomList[i].name, Element(17, "Chlorine", "Cl", bAtomList[i].getMass()));
					bAtomList[i].setAtomicNumber(17);
				}*/
			else if(toupper(bAtomList[i].elem) == 'O'){
				bAtomList[i].bAtomType = new UnivalentAtom(bAtomList[i].name,
						                                   Element(8, "Oxygen", "O", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(8);
			}
			else if(toupper(bAtomList[i].elem) == 'F'){
				bAtomList[i].bAtomType = new
						UnivalentAtom(bAtomList[i].name, Element(9, "Fluorine", "F", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(9);
			}
				/*
				else if((toupper(bAtomList[i].name[0]) == 'B') && (toupper(bAtomList[i].name[0]) == 'R')){
				  bAtomList[i].bAtomType = new
					UnivalentAtom(bAtomList[i].name, Element(35, "Bromine", "Br", bAtomList[i].getMass()));
				  bAtomList[i].setAtomicNumber(35);
				}
				*/
			else if(toupper(bAtomList[i].elem) == 'I'){
				bAtomList[i].bAtomType = new
						UnivalentAtom(bAtomList[i].name, Element(53, "Iodine", "I", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(53);
			}
			else if(toupper(bAtomList[i].elem) == 'N'){
				bAtomList[i].bAtomType = new
						UnivalentAtom(bAtomList[i].name, Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(7);
			}
			bAtomList[i].bAtomType->setDefaultInboardBondLength(0.1112); // Just for initial construction
		}
			// Atoms with two bonds
		else if (bAtomList[i].nbonds == 2){
			if(toupper(bAtomList[i].elem) == 'H'){
				bAtomList[i].bAtomType = new
						BivalentAtom(bAtomList[i].name, Element(1, "Hydrogen", "H", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(1);
			}
			else if(toupper(bAtomList[i].elem) == 'C'){
				bAtomList[i].bAtomType = new
						BivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(6);
			}
			else if(toupper(bAtomList[i].elem) == 'O'){
				bAtomList[i].bAtomType = new
						BivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()),
						             109.47*Deg2Rad);
				bAtomList[i].setAtomicNumber(8);
			}
			else if(toupper(bAtomList[i].elem) == 'N'){
				bAtomList[i].bAtomType = new
						BivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(7);
			}
			else if(toupper(bAtomList[i].elem) == 'S'){
				bAtomList[i].bAtomType = new
						BivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()),
						             109.47*Deg2Rad);
				bAtomList[i].setAtomicNumber(16);
			}
			bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
		}
			// Atoms with three bonds
		else if (bAtomList[i].nbonds == 3){
			if(toupper(bAtomList[i].elem) == 'C'){
				bAtomList[i].bAtomType = new
						TrivalentAtom(bAtomList[i].name, Element(6, "Carbon", "C", bAtomList[i].getMass()),
						              120*Deg2Rad, 120*Deg2Rad
				);
				bAtomList[i].setAtomicNumber(6);
			}
			else if(toupper(bAtomList[i].elem) == 'O'){
				bAtomList[i].bAtomType = new
						TrivalentAtomTetra(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(8);
			}
			else if(toupper(bAtomList[i].elem) == 'N'){
				bAtomList[i].bAtomType = new
						TrivalentAtomTetra(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(7);
			}
			else if(toupper(bAtomList[i].elem) == 'S'){
				bAtomList[i].bAtomType = new
						TrivalentAtomTetra(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(16);
			}
			else if(toupper(bAtomList[i].elem) == 'P'){
				bAtomList[i].bAtomType = new
						TrivalentAtomTetra(bAtomList[i].name,  Element(15, "Phosphorus", "P", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(15);
			}
			bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
		}
			// Atoms with four bonds
		else if (bAtomList[i].nbonds == 4){
			if(toupper(bAtomList[i].elem) == 'C'){
				bAtomList[i].bAtomType = new
						QuadrivalentAtom(bAtomList[i].name,  Element(6, "Carbon", "C", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(6);
			}
			else if(toupper(bAtomList[i].elem) == 'O'){
				bAtomList[i].bAtomType = new
						QuadrivalentAtom(bAtomList[i].name,  Element(8, "Oxygen", "O", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(8);
			}
			else if(toupper(bAtomList[i].elem) == 'N'){
				bAtomList[i].bAtomType = new
						QuadrivalentAtom(bAtomList[i].name,  Element(7, "Nitrogen", "N", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(7);
			}
			else if(toupper(bAtomList[i].elem) == 'S'){
				bAtomList[i].bAtomType = new
						QuadrivalentAtom(bAtomList[i].name,  Element(16, "Sulfur", "S", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(16);
			}
			else if(toupper(bAtomList[i].elem) == 'P'){
				bAtomList[i].bAtomType = new
						QuadrivalentAtom(bAtomList[i].name,  Element(15, "Phosphorus", "P", bAtomList[i].getMass()));
				bAtomList[i].setAtomicNumber(15);
			}
			bAtomList[i].bAtomType->setDefaultInboardBondLength(0.19);
		}


	} // Finish assigning Compound::SingleAtoms
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
	SetGmolAtomsMolmodelTypes();
}

/** Print atom and bonds list with details**/
void Topology::PrintAtomList(int whichWorld)
{
	// Atoms
	std::cout<<"Topology::PrintAtomList\n";
	for(unsigned int i = 0; i < bAtomList.size(); i++){
		bAtomList[i].Print(whichWorld);
	}

	// Bonds
	for(unsigned int i = 0; i < bonds.size(); i++){
		bonds[i].Print(whichWorld);
	}
}

/** Biotype is a Molmodel hook that is usually used to look up molecular
force field specific parameters for an atom type. Gmolmodel defines a
new Biotype for each atom. The only thing that is specified is the element
with info about name, atomic number, valence and mass. **/
void Topology::bAddBiotypes(
	  //std::string resName,
	readAmberInput *amberReader
	//, SimTK::DuMMForceFieldSubsystem& dumm
)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Iterate through atoms and define Biotypes based on resname
	for(int i = 0; i < amberReader->getNumberAtoms(); i++){
		SimTK::BiotypeIndex biotypeIndex = SimTK::Biotype::defineBiotype(
			  SimTK::Element(
				  bAtomList[i].getAtomicNumber(),
				  (std::to_string(bAtomList[i].getElem())).c_str(),
				  (std::to_string(bAtomList[i].getElem())).c_str(),
				  bAtomList[i].getMass()),
			bAtomList[i].getNBonds(),
			resName.c_str(),
			(bAtomList[i].getName()).c_str(),
			SimTK::Ordinality::Any
		);

		bAtomList[i].setBiotypeIndex(biotypeIndex);

		// Assign atom's biotype as a composed name: name + force field type
		//std::string temp(bAtomList[i].name); // @ Restore MULMOL
		std::string temp(name + bAtomList[i].name); // DEL MULMOL
		temp += bAtomList[i].fftype;
		bAtomList[i].setBiotype(temp);
		//std::cout << "Added Biotype " << temp << " with BiotypeIndex " << biotypeIndex << std::endl;
	}
}

/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
setBiotypeChargedAtomType for every atom. These Molmodel functions contain
information regarding the force field parameters. **/
void Topology::generateDummAtomClasses(
				  std::string resName
				, readAmberInput *amberReader
				, SimTK::DuMMForceFieldSubsystem& dumm
)
{
/*
	// Define AtomClasses
	SimTK::DuMM::AtomClassIndex aCIx;

	// Iterate through amberReader atoms and define AtomClasses
	for(int i = 0; i < amberReader->getNumberAtoms(); i++){
		// Get an AtomClass index
		aCIx = dumm.getNextUnusedAtomClassIndex();
		bAtomList[i].setAtomClassIndex(aCIx);

		// Define an AtomClass name
		std::string atomClassName =
				//std::string("top") // restore MULMOL
				name // del MULMOL
				+ resName
				+ bAtomList[i].getFftype()
				+ std::string("_")
				+ std::to_string(bAtomList[i].getNumber());

		// Define an AtomClass (has info about van der Waals)
		dumm.defineAtomClass(
			aCIx,
			atomClassName.c_str(),
			bAtomList[i].getAtomicNumber(), // int atomicNumber
			bAtomList[i].getNBonds(), // expected valence
			bAtomList[i].getVdwRadius() / 10.0, // nm
			bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
		);

		//std::cout << "Defined AtomClass " << atomClassName << " with atomClassIndex " << aCIx << std::endl;
	}
*/

	// Iterate through Gmolmodel atoms and define AtomClasses
	for(int i = 0; i < amberReader->getNumberAtoms(); i++){

		// Get AtomClass parameters from bSpecificAtom info
		AtomClassParams atomParams(
			bAtomList[i].getAtomicNumber(),
			bAtomList[i].getNBonds(),
			bAtomList[i].getVdwRadius() / 10.0, // nm
			bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
		);

		SimTK::DuMM::AtomClassIndex aCIx;
		std::string atomClassName;

		// Check if we already have this Atom Class
		bool foundit = false;
		std::map<AtomClassParams, AtomClassId>::const_iterator it;
		for (	it = aClassParams2aClassId.begin();
			it != aClassParams2aClassId.end(); ++it){
			if(it->first == atomParams){
				foundit = true;
				break;
			}
		}

		if ( !foundit ) { // we don't have this AtomClass
			// Get an AtomClass index from DuMM
			aCIx = dumm.getNextUnusedAtomClassIndex();
		
			// Define an AtomClass name
			atomClassName = bAtomList[i].getFftype();

			//std::cout << "!foundit aCIx atomClassName " << aCIx << " " << atomClassName  << std::endl;
			// Define an AtomClass
			dumm.defineAtomClass(aCIx, atomClassName.c_str(),
				atomParams.atomicNumber,
				atomParams.valence,
				atomParams.vdwRadius,
				atomParams.LJWellDepth
			);

			std::cout << "Topology::bAddDummAtomClasses "; atomParams.dump();

			// Insert an entry in our map too
			AtomClassId atomClassId(aCIx, atomClassName);
			aClassParams2aClassId.insert( std::make_pair(atomParams, atomClassId) );

		}else{ // we already have this AtomClass

			//aClassParams2aClassId.at(atomParams);
			//AtomClassId& atomClassId = aClassParams2aClassId.at(atomParams);
			const AtomClassId& atomClassId = it->second;

			aCIx = atomClassId.index;
			atomClassName = atomClassId.name;
	
			//std::cout << "foundit  aCIx atomClassName " << aCIx << " " << atomClassName << std::endl;
		}

		// Insert AtomClass index in Gmolmodel atom list too
		bAtomList[i].setAtomClassIndex(aCIx);

	} // --------------------------- DONE



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
		chargedAtomTypeName += bAtomList[k].biotype;

		// Define a ChargedAtomType (AtomClass with a charge)
		dumm.defineChargedAtomType(
		  chargedAtomTypeIndex,
		  chargedAtomTypeName.c_str(),
		  bAtomList[k].getAtomClassIndex(),
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

/** It calls DuMMs defineAtomClass, defineChargedAtomTye and
setBiotypeChargedAtomType for every atom. These Molmodel functions contain
information regarding the force field parameters. **/
void Topology::transferDummAtomClasses(
				  std::string resName
				, readAmberInput *amberReader
				, SimTK::DuMMForceFieldSubsystem& dumm
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

		aCIx = atomClassId.index;
		atomClassName = atomClassId.name;

		std::cout << "Topology::transferAtomClasses " 
			<< aCIx << " " << atomClassName ;
		atomParams.dump();

		// Define an AtomClass
		dumm.defineAtomClass(aCIx, atomClassName.c_str(),
			atomParams.atomicNumber,
			atomParams.valence,
			atomParams.vdwRadius,
			atomParams.LJWellDepth
		);

	}
	// --------------------------- DONE


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
		chargedAtomTypeName += bAtomList[k].biotype;

		// Define a ChargedAtomType (AtomClass with a charge)
		dumm.defineChargedAtomType(
		  chargedAtomTypeIndex,
		  chargedAtomTypeName.c_str(),
		  bAtomList[k].getAtomClassIndex(),
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
const void Topology::PrintMolmodelAndDuMMTypes(SimTK::DuMMForceFieldSubsystem& dumm) const
{
	std::cout << "Print Molmodel And DuMM Types:" << std::endl;
	for(size_t i = 0; i < bAtomList.size(); i++){
		std::cout << " list ix " << i
			<< " CompoundAtomIndex " << bAtomList[i].getCompoundAtomIndex()
			<< " biotypename " << bAtomList[i].biotype
			<< " name " << bAtomList[i].name
			<< " BiotypeIndex " << bAtomList[i].getBiotypeIndex()
			<< " ChargedAtomTypeIndex "<< bAtomList[i].getChargedAtomTypeIndex()
			<< " AtomClassIx " << bAtomList[i].getAtomClassIndex()
			<< " partialChargeInE " << bAtomList[i].charge
			<< " chargedAtomTypeIndex "
			<< bAtomList[i].getChargedAtomTypeIndex()
			<< " DuMM VdW Radius "
			<< dumm.getVdwRadius(bAtomList[i].getAtomClassIndex())
			<< " DuMM VdW Well Depth "
			<< dumm.getVdwWellDepth(bAtomList[i].getAtomClassIndex())
			<< std::endl;
	}
}

/** Calls DuMM defineBondStretch to define bonds parameters. **/
void Topology::bAddDummBondParams(std::string, readAmberInput *amberReader, SimTK::DuMMForceFieldSubsystem& dumm)
{

	// Keep track of inserted AtomClass pairs
	std::vector<std::vector<SimTK::DuMM::AtomClassIndex>> ixs;

	// Iterate through bonds and define their parameters
	// Suppose or try to have the same order as the reader
	for(int t = 0; t < amberReader->getNumberBonds(); t++){
		
		// Generate a pair of atom classes for this bond
		std::vector<SimTK::DuMM::AtomClassIndex> bondStretchIxs;
		bondStretchIxs.push_back( SimTK::DuMM::AtomClassIndex(
			(bAtomList[bonds[t].i]).getAtomClassIndex()) );
		bondStretchIxs.push_back( SimTK::DuMM::AtomClassIndex(
			(bAtomList[bonds[t].j]).getAtomClassIndex()) );

		// Check if we already have this bond
		bool foundit = false;
		for(auto& row:ixs){
			if ((bondStretchIxs[0] == row[0])
			 && (bondStretchIxs[1] == row[1])){
				foundit = true;	break;
			}
		}

		if (  !foundit ){ // bond was not found
		
			dumm.defineBondStretch_KA(
				(bAtomList[bonds[t].i]).getAtomClassIndex(),
				(bAtomList[bonds[t].j]).getAtomClassIndex(),
				amberReader->getBondsForceK(t),  //k1
				amberReader->getBondsEqval(t)   //equil1
			);

		}
	}
}

/** Calls DuMM defineBondBend to define angle parameters. **/
void Topology::bAddDummAngleParams(std::string, readAmberInput *amberReader, SimTK::DuMMForceFieldSubsystem& dumm)
{
	// function args were std::string resName

	// Iterate through angles and define their parameters
	for(int t = 0; t < amberReader->getNumberAngles(); t++){
		dumm.defineBondBend_KA(
			bAtomList[amberReader->getAnglesAtomsIndex1(t)].getAtomClassIndex(),
			bAtomList[amberReader->getAnglesAtomsIndex2(t)].getAtomClassIndex(),
			bAtomList[amberReader->getAnglesAtomsIndex3(t)].getAtomClassIndex(),
			amberReader->getAnglesForceK(t),
			static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE 
				* amberReader->getAnglesEqval(t))) // TODO 32 vs 64 bit
		);
	}
}

/** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
void Topology::bAddDummTorsionParams(
	  std::string resName
	, readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
)
{
	std::vector<std::pair<int, int>> pairStartAndLens = amberReader->getPairStartAndLen();

	for(unsigned int index=0; index<pairStartAndLens.size(); index++){

		int first    = pairStartAndLens[index].first;
		int numberOf = pairStartAndLens[index].second;

		for(int t = first; t < (first + numberOf); t++){
			if(numberOf == 1){
				dumm.defineBondTorsion_KA(
					bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
					static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
					amberReader->getDihedralsForceK(t),
					static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)))
				);
			}
			else if(numberOf == 2){
				dumm.defineBondTorsion_KA(
					bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
					static_cast<int>(amberReader->getDihedralsPeriod(t)), // TODO wants int, returns double
					amberReader->getDihedralsForceK(t),
					static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))),
					static_cast<int>(amberReader->getDihedralsPeriod(t + 1)), // TODO wants int, returns double
					amberReader->getDihedralsForceK(t + 1),
					static_cast<SimTK::Real>(ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)))
				);
			}
			else if(numberOf == 3){
				dumm.defineBondTorsion_KA(
					bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
					bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
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
		}
	}
}

/** Adds force field parameters read by the inputReader **/
void Topology::generateDummParams(
	readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Add types
	generateDummAtomClasses(resName, amberReader, dumm);

	// Add parameters
	bAddDummBondParams(resName, amberReader, dumm);
	bAddDummAngleParams(resName, amberReader, dumm);
	bAddDummTorsionParams(resName, amberReader, dumm);
}

/** Transfer already generated force field parameters to DuMM **/
void Topology::transferDummParams(
	readAmberInput *amberReader
	, SimTK::DuMMForceFieldSubsystem& dumm
)
{
	// We don't have any residues. The whole molecule is one residue
	std::string resName = this->name;

	// Add types
	transferDummAtomClasses(resName, amberReader, dumm);

	// Add parameters
	bAddDummBondParams(resName, amberReader, dumm);
	bAddDummAngleParams(resName, amberReader, dumm);
	bAddDummTorsionParams(resName, amberReader, dumm);
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

void Topology::loadTriples()
{
	// Assign Compound coordinates by matching bAtomList coordinates
	//std::cout << "Topology triples: " << std::endl ;
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
	}
	std::vector< std::vector<Compound::AtomIndex> > bondedAtomRuns =
	getBondedAtomRuns(3, atomTargets);

	// Find root bAtomList index
	// int rootIx;
	int ix = -1;
	for(auto atom: bAtomList){
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
SimTK::Real Topology::calcLogSineSqrGamma2(const SimTK::State &quatState) {
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
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
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


SimTK::Real Topology::calcLogDetMBATDistsMassesContribution(const SimTK::State&){
	// function args were const SimTK::State& someState

	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
			Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
			atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
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
			atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
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


SimTK::Real Topology::calcLogDetMBATMassesContribution(const SimTK::State&){
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

/** **/
const bBond& Topology::getBond(int a1, int a2) const
{
	for(int i = 0; i < nbonds; i++){
		if( (bonds[i]).isThisMe(a1, a2) ){
			return bonds[i];
		}
	}

	std::string assert_string("No bond with these atom indeces found: " + to_string(a1) + " " + to_string(a2) + ". Exiting.");
	std::cout << assert_string << std::endl;
	assert(0);

	// TODO what should we return? std::optional? this std::exception looks sketchy
	// anyways, return {}; returns a local and it's not good
	throw std::exception();
	return {};
}

/** Get bond order. **/
int Topology::getBondOrder(int, int) const
{
	assert(!"Not implemented."); throw std::exception();

	return 0;
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
	if( node->visited ){
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
		if ((*bondsInvolvedIter)->isThisMe(node->number, previousNode->number) ) {
			(*bondsInvolvedIter)->setVisited(1);

			// Skip the first step as we don't have yet two atoms
			if (nofProcesses != 1) {

				// The first bond is special in Molmodel and has to be
				// treated differently. Set a base atom first
				if (nofProcesses == 2) {
					if (baseSetFlag == 0) {
						this->setBaseAtom(*(previousNode->bAtomType));
						this->setAtomBiotype(previousNode->name, (this->name), previousNode->getName());
						this->convertInboardBondCenterToOutboard();
						baseSetFlag = 1;
					}
				}

				// Bond current node by the previous (Compound function)
				std::stringstream parentBondCenterPathName;
				if (previousNode->number == baseAtomNumber) {
					parentBondCenterPathName << previousNode->name
						<< "/bond" << previousNode->freebonds;
				} else {
					parentBondCenterPathName << previousNode->name
						<< "/bond" << (previousNode->nbonds - previousNode->freebonds + 1);
				}

				// THIS IS WHERE WE PERFORM THE ACTUAL BONDING
				// (Compound::SingleAtom&, BondCenterPathName, Length, Angle
				std::string debugString = parentBondCenterPathName.str();
				this->bondAtom(*node->bAtomType,
						(parentBondCenterPathName.str()).c_str(), 0.149, 0);

				// Set the final Biotype
				this->setAtomBiotype(node->name, (this->name).c_str(), node->getName());

				// Set bSpecificAtom atomIndex to the last atom added to bond
				node->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 1);


				// The only time we have to set atomIndex to the previous node
				if (nofProcesses == 2) {
					previousNode->atomIndex = getBondAtomIndex(Compound::BondIndex(getNumBonds() - 1), 0);
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

				//std::cout << "DEBUG inserted into GmolBond2bondIx bondIx2GmolBond  " << GmolBond2bondIx.size() << " " << bondIx2GmolBond.size() << std::endl << std::flush;

				// Drop the number of available bonds
				--previousNode->freebonds;
				--node->freebonds;

				// Bond was inserted in Molmodel Compound. Get out and search
				// the next bond
				break;

			}
		}
	}

	// Mark the node as visited
	node->visited = 1;

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
			if(leftNode->number == baseAtomNumber){
				sbuff << leftNode->name << "/bond" << leftNode->freebonds;
			}else{
				sbuff << leftNode->name << "/bond"
					<< (leftNode->nbonds - leftNode->freebonds + 1);
			}

			std::stringstream otsbuff;
			if(rightNode->number == baseAtomNumber){
				otsbuff << rightNode->name << "/bond" << rightNode->freebonds;
			}else{
				otsbuff << rightNode->name << "/bond"
					<< (rightNode->nbonds - rightNode->freebonds + 1);
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
			this->setAtomBiotype(leftNode->name, (this->name), leftNode->getName());
			this->setAtomBiotype(rightNode->name, (this->name), rightNode->getName());

			--leftNode->freebonds;
			--rightNode->freebonds;

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
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
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
		baseAtomNumber = root->number;
	}

	// Build the graph
	nofProcesses = 0;
	buildAcyclicGraph(root, root);
    std::cout << "Topology::buildGraphAndMatcoords buildAcyclicGraph done\n" << std::flush;

	// Close the remaining bonds
	addRingClosingBonds();
    std::cout << "Topology::buildGraphAndMatcoords  addRingClosingBonds done\n" << std::flush;

	// Build the conformation
	matchDefaultConfigurationWithAtomList(SimTK::Compound::Match_Exact);

	// Now that everything is built, initialize aIx2TopTransforms map
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		aIx2TopTransform.insert(std::make_pair((bAtomList[i]).atomIndex, SimTK::Transform()));
	}

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
void Topology::setFlexibility(std::string argRegimen, std::string flexFN, int whichWorld){

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
	//PrintAtomList();
	this->regimen = argRegimen;
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void Topology::loadAIx2MbxMap(){

	// If the map is empty fill with empty vectors first
	if(aIx2mbx.empty()){
		// Iterate through atoms and get their MobilizedBodyIndeces
		for (unsigned int i = 0; i < getNumAtoms(); ++i) {

			// Get atomIndex from atomList
			SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;

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
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;

		// Get MobilizedBodyIndex from CompoundAtom
		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);

		// Insert
		//aIx2mbx.insert(
		//		std::pair<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>
		//		(aIx, mbx));
		aIx2mbx[aIx].emplace_back(mbx);
	}


}

//void Topology::loadMbx2AIxMap(){
//	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
//
//		// Get atomIndex from atomList
//		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;
//
//		// Get MobilizedBodyIndex from CompoundAtom
//		SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);
//
//		// Map mbx2aIx contains only atoms at the origin of mobods
//		//std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex >
//		//        pairToBeInserted(mbx, aIx);
//		if (getAtomLocationInMobilizedBodyFrame(aIx) == 0) {
//			mbx2aIx.insert(
//				std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>
//				(mbx, aIx));
//		}
//
//	}
//
//}


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
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;
		aIx2TopTransform[aIx] = calcDefaultAtomFrameInCompoundFrame(aIx);
	}
}

/**  **/
void Topology::printTopTransforms()
{
	std::cout << "Topology TopTransforms " << std::endl;
	for (unsigned int i = 0; i < getNumAtoms(); ++i) {
		SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;
		std::cout << aIx << " " << aIx2TopTransform[aIx] << std::endl;
	}
}

SimTK::Transform Topology::getTopTransform(SimTK::Compound::AtomIndex aIx)
{
	return aIx2TopTransform[aIx];
}

SimTK::Transform& Topology::calcDefaultAtomFrameInCompoundFrameThroughDuMM(
	SimTK::Compound::AtomIndex aIx,
	SimTK::DuMMForceFieldSubsystem& dumm,
	SimTK::SimbodyMatterSubsystem& matter,
	const SimTK::State& someState)
{
	SimTK::DuMM::AtomIndex dAIx = getDuMMAtomIndex(aIx);
	const  SimTK::MobilizedBodyIndex mbx = dumm.getAtomBody(dAIx);
	const SimTK::MobilizedBody& mobod = matter.getMobilizedBody(mbx);

	// Body transforms
	const Transform&    X_GB = mobod.getBodyTransform(someState);
	const Rotation&     R_GB = X_GB.R();
	const Vec3&         p_GB = X_GB.p();

	// Atom
	SimTK::Vec3 station = getAtomLocationInMobilizedBodyFrameThroughDumm(aIx, dumm);

	// return
	//const Vec3 p_BS_G = R_GB * station;
	//return X;

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
	if(int nofDigits = static_cast<int>(std::to_string(index).size());
	maxNofDigits > nofDigits){
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
				, bAtomList[i].inName  // name
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
const CompoundSystem::CompoundIndex &Topology::getCompoundIndex() const {
	return compoundIndex;
}

/** Set the compoundIndex which is the position in the vector of Compounds
 * of the CompoundSystem **/
void Topology::setCompoundIndex(const CompoundSystem::CompoundIndex &compoundIndex) {
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
					if((*bondsInvolvedIter)->isThisMe(originSpecAtom->getNumber(), originSpecAtom->neighbors[k]->getNumber())){

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

// This function is only intended for root atoms!!
/*
std::vector<SimTK::Transform>
Topology::calcMobodToMobodTransforms(
	SimTK::SimbodyMatterSubsystem *matter,
	SimTK::Compound::AtomIndex aIx,
	const SimTK::State& someState,
	SimTK::DuMMForceFieldSubsystem& dumm,
	int whichWorld)
{
	// There is no P_X_F and B_X_M inside a body.
	assert(getAtomLocationInMobilizedBodyFrame(aIx) == 0);

	// Get body, parentBody
	SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);
	const SimTK::MobilizedBody& mobod = matter->getMobilizedBody(mbx);
	const SimTK::MobilizedBody& parentMobod =  mobod.getParentMobilizedBody();
	SimTK::MobilizedBodyIndex parentMbx = parentMobod.getMobilizedBodyIndex();

	// Get the neighbor atom in the parent mobilized body
	//SimTK::Compound::AtomIndex chemParentAIx = getChemicalParent(matter, aIx); // SAFE
	SimTK::Compound::AtomIndex chemParentAIx = getChemicalParent(matter, aIx, dumm); // DANGER

	// Get parent-child BondCenters relationship
	SimTK::Transform X_parentBC_childBC =
			getDefaultBondCenterFrameInOtherBondCenterFrame(aIx, chemParentAIx);

	// Get Top to parent frame
	SimTK::Compound::AtomIndex parentRootAIx = mbx2aIx[parentMbx];
	SimTK::Transform T_X_Proot = getTopTransform(parentRootAIx);
	SimTK::Transform Proot_X_T = ~T_X_Proot;

	// Get inboard dihedral angle
	SimTK::Angle inboardBondDihedralAngle = bgetDefaultInboardDihedralAngle(aIx);
	SimTK::Transform InboardDihedral_XAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::XAxis);
	SimTK::Transform InboardDihedral_ZAxis
		= SimTK::Rotation(inboardBondDihedralAngle, SimTK::ZAxis);

	// Get the old PxFxMxB transform
	SimTK::Transform X_to_Z = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis); // aka M_X_pin
	SimTK::Transform Z_to_X = ~X_to_Z;
	SimTK::Transform oldX_PB =
		(~T_X_Proot) * getTopTransform(aIx)
		* InboardDihedral_XAxis * X_to_Z
		* InboardDihedral_ZAxis * Z_to_X;

	// B_X_Ms
	SimTK::Transform B_X_M = SimTK::Rotation(-90*SimTK::Deg2Rad, SimTK::YAxis); // aka M_X_pin
	SimTK::Transform B_X_M_anglePin = X_parentBC_childBC;
	SimTK::Transform B_X_M_univ = X_parentBC_childBC
		* SimTK::Transform(Rotation(-90*Deg2Rad, XAxis)); // Move rotation axis Y to Z

	// P_X_Fs = old P_X_B * B_X_M
	SimTK::Transform P_X_F = oldX_PB * B_X_M;
	SimTK::Transform P_X_F_anglePin = oldX_PB * B_X_M_anglePin;
	SimTK::Transform P_X_F_univ = oldX_PB * B_X_M;

	// Get mobility (joint type)
	bSpecificAtom *atom = updAtomByAtomIx(aIx);
	SimTK::BondMobility::Mobility mobility;
	bBond bond = getBond(getNumber(aIx), getNumber(chemParentAIx));
	mobility = bond.getBondMobility(whichWorld);

	bool pinORslider =
		(mobility == SimTK::BondMobility::Mobility::Torsion)
		|| (mobility == SimTK::BondMobility::Mobility::AnglePin)
		|| (mobility == SimTK::BondMobility::Mobility::Slider);

	if( (pinORslider) && ((atom->neighbors).size() == 1)){

		return std::vector<SimTK::Transform> {P_X_F_anglePin, B_X_M_anglePin};

	}else if( (pinORslider) && ((atom->neighbors).size() != 1)){

		return std::vector<SimTK::Transform> {P_X_F, B_X_M};

	}else if((mobility == SimTK::BondMobility::Mobility::BallM)
	|| (mobility == SimTK::BondMobility::Mobility::Rigid)
	|| (mobility == SimTK::BondMobility::Mobility::Translation) // Cartesian
	){

		return std::vector<SimTK::Transform> {P_X_F, B_X_M};

	}else{

		std::cout << "Warning: unknown mobility." << std::endl;
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
*/


