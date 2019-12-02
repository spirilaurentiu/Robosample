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
        boost::trim(str_buf);
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
void Topology::PrintAtomList()
{
    // Atoms
    std::cout<<"Topology::PrintAtomList\n";
    for(unsigned int i = 0; i < bAtomList.size(); i++){
        bAtomList[i].Print();
    }

    // Bonds
    for(unsigned int i = 0; i < bAtomList.size(); i++){
        bonds[i].Print();
    }
}

/** Biotype is a Molmodel hook that is usually used to look up molecular
force field specific parameters for an atom type. Gmolmodel defines a
new Biotype for each atom. The only thing that is specified is the element
with info about name, atomic number, valence and mass. **/
void Topology::bAddBiotypes(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
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
void Topology::bAddAtomClasses(
                  std::string resName
                , readAmberInput *amberReader
                , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Define AtomClasses
    SimTK::DuMM::AtomClassIndex aCIx;

    // Iterate through amberReader atoms and define AtomClasses
    for(int i = 0; i < amberReader->getNumberAtoms(); i++){
        // Get an AtomClass index
        aCIx = dumm.getNextUnusedAtomClassIndex();
        bAtomList[i].setAtomClassIndex(aCIx);

        // Define an AtomClass name
        const char* atomClassName = (
                //std::string("top") // restore MULMOL
                name // del MULMOL 
                + resName
                + bAtomList[i].getFftype()
                + std::string("_")
                + std::to_string(bAtomList[i].getNumber()) ).c_str();

        // Define an AtomClass (has info about van der Waals)
        dumm.defineAtomClass(
            aCIx,
            atomClassName,
            bAtomList[i].getAtomicNumber(), // int atomicNumber
            bAtomList[i].getNBonds(), // expected valence
            bAtomList[i].getVdwRadius() / 10.0, // nm
            bAtomList[i].getLJWellDepth() * 4.184 // kcal to kJ
        );

        //std::cout << "Defined AtomClass " << atomClassName << " with atomClassIndex " << aCIx << std::endl;
    }

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
void Topology::PrintMolmodelAndDuMMTypes(SimTK::DuMMForceFieldSubsystem& dumm)
{
    std::cout << "Print Molmodel And DuMM Types:" << std::endl;
    for(size_t i = 0; i < bAtomList.size(); i++){
        std::cout << " list ix " << i
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
void Topology::bAddBondParams(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Iterate through bonds and define their parameters
    // Suppose or try to have the same order as the reader
    for(int t = 0; t < amberReader->getNumberBonds(); t++){
        dumm.defineBondStretch_KA(
            (bAtomList[bonds[t].i]).getAtomClassIndex(), 
            (bAtomList[bonds[t].j]).getAtomClassIndex(),
            amberReader->getBondsForceK(t),  //k1
            amberReader->getBondsEqval(t)   //equil1
        );

	//std::cout << "Defined bond stretch between aCIx1 aCIx2 k b0 " 
	//	<< (bAtomList[bonds[t].i]).getAtomClassIndex() << " "
	//	<< (bAtomList[bonds[t].j]).getAtomClassIndex() << " "
	//	<< amberReader->getBondsForceK(t) << " "
	//	<< amberReader->getBondsEqval(t) << std::endl;

    }
}

/** Calls DuMM defineBondBend to define angle parameters. **/
void Topology::bAddAngleParams(
      std::string resName
    , readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // Iterate through angles and define their parameters
    for(int t = 0; t < amberReader->getNumberAngles(); t++){
        dumm.defineBondBend_KA(
            bAtomList[amberReader->getAnglesAtomsIndex1(t)].getAtomClassIndex(),
            bAtomList[amberReader->getAnglesAtomsIndex2(t)].getAtomClassIndex(),
            bAtomList[amberReader->getAnglesAtomsIndex3(t)].getAtomClassIndex(),
            amberReader->getAnglesForceK(t),
            ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getAnglesEqval(t))
        );
    }
}

/** Calls DuMM defineBondTorsion for 1, 2 and 3 periodicities **/
void Topology::bAddTorsionParams(
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
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),   
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t))
                );
            }
            else if(numberOf == 2){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1), 
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1))
                );
            }
            else if(numberOf == 3){
                dumm.defineBondTorsion_KA(
                    bAtomList[amberReader->getDihedralsAtomsIndex1(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex2(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex3(t)].getAtomClassIndex(),
                    bAtomList[amberReader->getDihedralsAtomsIndex4(t)].getAtomClassIndex(),
                    amberReader->getDihedralsPeriod(t),   amberReader->getDihedralsForceK(t),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t)),
                    amberReader->getDihedralsPeriod(t+1), amberReader->getDihedralsForceK(t+1),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+1)),
                    amberReader->getDihedralsPeriod(t+2), amberReader->getDihedralsForceK(t+2),
                    ANG_360_TO_180(SimTK_RADIAN_TO_DEGREE * amberReader->getDihedralsPhase(t+2))
                );
            }
        }
    }
}

/** Adds force field parameters read by the inputReader **/
void Topology::bAddAllParams(
    readAmberInput *amberReader
    , SimTK::DuMMForceFieldSubsystem& dumm
)
{
    // We don't have any residues. The whole molecule is one residue
    std::string resName = this->name;

    // Add types
    bAddBiotypes(resName, amberReader, dumm); 
    bAddAtomClasses(resName, amberReader, dumm);

    // Add parameters
    bAddBondParams(resName, amberReader, dumm); 
    bAddAngleParams(resName, amberReader, dumm); 
    bAddTorsionParams(resName, amberReader, dumm); 
}

void Topology::loadTriples(void)
{
	// Assign Compound coordinates by matching bAtomList coordinates
	std::map<AtomIndex, Vec3> atomTargets;
	for(int ix = 0; ix < getNumAtoms(); ++ix){
		Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
		atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
	}
	std::vector< std::vector<Compound::AtomIndex> > bondedAtomRuns =
	getBondedAtomRuns(3, atomTargets);	

	//std::cout << "Topology triples: " << std::endl ;
	int bIx = -1;
	for(auto bAR: bondedAtomRuns){
		bIx++;

		if(bAR[0] < bAR[2]){
			triples.push_back(bAR);
			for(auto aIx: triples.back()){
				//std::cout << " " << aIx;
			}
			//std::cout << std::endl;
		}
	}
}

SimTK::Real Topology::calcLogDetMBATGamma2Contribution(const SimTK::State& quatState){
	//State& eulerState;
	//matter.convertToEulerAngles(quatState, eulerState);
	//std::cout << "calcLogDetMBATGamma2Contribution quaternionState " << quatState << std::endl;
	//std::cout << "calcLogDetMBATGamma2Contribution      eulerState " << eulerState << std::endl;

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
	//std::cout << "sin pitch " << sinPitch << std::endl;

	//SimTK::Real pitch = std::asin(sinPitch);
	//std::cout << "pitch " << pitch << std::endl;
	//if(pitch < 0){
	//	pitch = pitch + (2*SimTK_PI);
	//	std::cout << "sin converted pitch " << std::sin(pitch) << std::endl;
	//}

	SimTK::Real result = std::log(sinPitch * sinPitch);

	return result;
}

SimTK::Real Topology::calcLogDetMBATDistsMassesContribution(const SimTK::State& someState){
        // Assign Compound coordinates by matching bAtomList coordinates
        std::map<AtomIndex, Vec3> atomTargets;
        for(int ix = 0; ix < getNumAtoms(); ++ix){
                Vec3 vec(bAtomList[ix].getX(), bAtomList[ix].getY(), bAtomList[ix].getZ());
                atomTargets.insert(pair<AtomIndex, Vec3> (bAtomList[ix].atomIndex, vec));
        }

	//std::cout << "Topology::calcLogDetMBATDistsContribution dists squared: " ;
	SimTK::Real result = std::log(bAtomList[bonds[0].i].mass);
	for(auto bond: bonds){
		//SimTK::Vec3 vec0 = calcAtomLocationInGroundFrame(someState, triple[0]);
		//SimTK::Vec3 vec1 = calcAtomLocationInGroundFrame(someState, triple[1]);
		//SimTK::Vec3 vec2 = calcAtomLocationInGroundFrame(someState, triple[2]);

		SimTK::Vec3 atom1pos = SimTK::Vec3(bAtomList[bond.i].getX(), bAtomList[bond.i].getY(), bAtomList[bond.i].getZ());
		SimTK::Vec3 atom2pos = SimTK::Vec3(bAtomList[bond.j].getX(), bAtomList[bond.j].getY(), bAtomList[bond.j].getZ());
		SimTK::Real distSqr = (atom2pos - atom1pos).normSqr();

		//std::cout << "atom " << bond.j << " c " << std::log(distSqr) + std::log(distSqr) + std::log(bAtomList[bond.j].mass) << " " ;

		result = result + std::log(distSqr) + std::log(distSqr) + std::log(bAtomList[bond.j].mass);

	}
	//std::cout << std::endl;

	return result;
}

/* Calculate angle contribution at the MBAT determinant: 
TODO: calculate atomTargets only once: getAtomLocationsInGround
*/
SimTK::Real Topology::calcLogDetMBATAnglesContribution(const SimTK::State& someState){
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
		SimTK::Real angle = std::acos(dotProduct);
		//std::cout << angle << " " ;

		SimTK::Real sinSqAngle = 1 - (dotProduct * dotProduct);

		result = result + std::log(sinSqAngle);

	}
	//std::cout << std::endl;

	return result;
}

SimTK::Real Topology::calcLogDetMBATMassesContribution(const SimTK::State& someState)
{
	std::cout << "Topology::calcLogDetMBATMassesContribution masses: " ;
	SimTK::Real result = 1.0;
	for(auto atom: bAtomList){
		std::cout << atom.mass << " " ;
		result = result * (atom.mass * atom.mass * atom.mass);
	}
	std::cout << std::endl;
}

SimTK::Real Topology::calcLogDetMBAT(const SimTK::State& someState)
{
	SimTK::Real gamma2Contribution = calcLogDetMBATGamma2Contribution(someState);
	SimTK::Real distsMassesContribution = calcLogDetMBATDistsMassesContribution(someState);
	SimTK::Real anglesContribution = calcLogDetMBATAnglesContribution(someState);
	//SimTK::Real massesContribution = calcLogDetMBATMassesContribution(someState);

	//std::cout << std::setprecision(20) << std::fixed;
	//std::cout << " gamma distsMasses angles contributions: " 
	//	<< gamma2Contribution << " " 
	//	<< distsMassesContribution << " " 
	//	<< anglesContribution << std::endl;

	return gamma2Contribution + distsMassesContribution + anglesContribution;
}




/**
 ** Interface **
 **/

/** Get the number of atoms in the molecule **/
int Topology::getNAtoms(void) const{
    return getNumAtoms();
}

/** Get the number of bonds in the molecule **/
int Topology::getNBonds(void) const{
    assert(!"Not implemented.");
}

/** Get a pointer to an atom object in the atom list inquiring
by number **/
bSpecificAtom * Topology::getAtomByNumber(int number) const{
    assert(!"Not implemented.");
}

/** Get a pointer to an atom object in the atom list inquiring
by its Molmodel assigned atom index (SimTK::Compound::AtomIndex) .**/
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
bSpecificAtom * Topology::getAtomByName(std::string name) const{assert(!"Not implemented.");}

/** Get the neighbours in the graph. **/
std::vector<bSpecificAtom *> Topology::getNeighbours(int) const{assert(!"Not implemented.");}

/** **/
bBond * Topology::getBond(int, int) const{assert(!"Not implemented.");}

/** Get bond order. **/
int Topology::getBondOrder(int, int) const{assert(!"Not implemented.");}


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

                // Perform the actual bonding
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
void Topology::buildGraphAndMatchCoords(
        SimTK::DuMMForceFieldSubsystem &dumm,
        int argRoot
) {

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
        std::cout << "Root atom will be chosen by Gmolmodel." << std::endl;
        int baseAtomListIndex = 0;
        for (int i = 0; i < natoms; i++) {
            if (bAtomList[i].getNBonds() > 1) {
                baseAtomListIndex = i;
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

    // Close the remaining bonds
    addRingClosingBonds();

    // Build the conformation
    matchDefaultConfigurationWithAtomList(SimTK::Compound::Match_Exact);

    // Implement flexibility/rigidity specifications
    //setRegimen(regimenSpec, flexFN);

}

/** Get regimen **/
std::string Topology::getRegimen(){
    return this->regimen;
}

/** Set regimen according to input file **/
void Topology::setFlexibility(std::string argRegimen, std::string flexFN){
    
    if(argRegimen.at(0) == 'I'){
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Free, Compound::BondIndex(r));
            bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                    BondMobility::Free);
        }
    }else if(argRegimen == "TD") {
        for (unsigned int r = 0; r < getNumBonds(); r++) {
            setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
            bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                    BondMobility::Torsion);
        }
    }else if(argRegimen == "BA"){
        for (unsigned int r=0 ; r<getNumBonds(); r++){

            int  firstSpecificAtomIndex = bonds[bondIx2GmolBond[Compound::BondIndex(r)]].i;
            int secondSpecificAtomIndex = bonds[bondIx2GmolBond[Compound::BondIndex(r)]].j;

            if((std::string(bAtomList[ firstSpecificAtomIndex].getFftype()) == "CA") ||
               (std::string(bAtomList[secondSpecificAtomIndex].getFftype() )== "CA")){
                setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
                bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                        BondMobility::Torsion);
            } else if( (bAtomList[ firstSpecificAtomIndex].getNBonds() < 3) ||
                (bAtomList[secondSpecificAtomIndex].getNBonds() < 3)){
                setBondMobility(BondMobility::Torsion, Compound::BondIndex(r));
                bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                        BondMobility::Torsion);
            }else{
                setBondMobility(BondMobility::Ball, Compound::BondIndex(r));
                bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                        BondMobility::Ball);
            }

            if(bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].isRingClosing()){
                bonds[bondIx2GmolBond.at(Compound::BondIndex(r))].setBondMobility(
                        BondMobility::Rigid);
            }

        }

    }else if(argRegimen.at(0) == 'R'){

        // Set all Compound and Topology bonds to rigid
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            setBondMobility(BondMobility::Rigid, SimTK::Compound::BondIndex(r));
        }
        for (unsigned int r=0 ; r<getNumBonds(); r++){
            bonds[r].setBondMobility(BondMobility::Rigid); // TODO: Change to rigid
        }

        // Get flexible bonds from file. Numbering starts at 0 in prmtop
        std::string line;
        std::ifstream F(flexFN);

        //printMaps();
/*        std::cout << "GmolBond2bondIx " << GmolBond2bondIx.size() << std::endl;
        std::cout << "GmolBond2bondIx:" << std::endl;
        for(unsigned int i = 0; i < nbonds; i++){
            std::cout << i << ' ' << GmolBond2bondIx.at(i) << std::endl;
        }*/

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
                if(lineWords.size() >= 2 ){
                    for(int i = 0; i < nbonds; i++){
                        if(bonds[i].isThisMe(
                            std::stoi(lineWords[0]), std::stoi(lineWords[1])) ){
                            if(lineWords.size() == 2) {
                                bonds[i].setBondMobility(BondMobility::Torsion);
                                setBondMobility(BondMobility::Torsion,
                                                GmolBond2bondIx.at(i));
                                break;
                            }else{
                                if(lineWords[2] == "Free"){
                                    bonds[i].setBondMobility(BondMobility::Free);
                                    setBondMobility(BondMobility::Free,
                                                    GmolBond2bondIx.at(i));
                                }else if((lineWords[2] == "Pin") || (lineWords[2] == "Torsion")) {
                                    bonds[i].setBondMobility(BondMobility::Torsion);
                                    setBondMobility(BondMobility::Torsion,
                                                    GmolBond2bondIx.at(i));
                                    break;
                                }else if(lineWords[2] == "Cylinder") {
                                    bonds[i].setBondMobility(BondMobility::Cylinder);
                                    setBondMobility(BondMobility::Cylinder,
                                                    GmolBond2bondIx.at(i));
                                    break;
                                }else if(lineWords[2] == "Ball") {
                                    bonds[i].setBondMobility(BondMobility::Ball);
                                    setBondMobility(BondMobility::Ball,
                                                    GmolBond2bondIx.at(i));
                                    break;
                                }else if(lineWords[2] == "Rigid"){
                                    bonds[i].setBondMobility(BondMobility::Rigid);
                                    setBondMobility(BondMobility::Rigid,
                                                    GmolBond2bondIx.at(i));
                                }else{
                                    bonds[i].setBondMobility(BondMobility::Torsion);
                                    setBondMobility(BondMobility::Torsion,
                                                    GmolBond2bondIx.at(i));
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

/*        std::cout << "Assigned mobilities:" << std::endl;
        for(unsigned int i = 0; i < nbonds; i++){
            std::cout << i << ' ' << GmolBond2bondIx.at(i) << " " << bonds[i].getBondMobility() << std::endl;
        }*/


    } // RB

    this->regimen = argRegimen;
}

/** Create MobilizedBodyIndex vs Compound::AtomIndex maps **/
void Topology::loadMobodsRelatedMaps(void){

    // Iterate through atoms and get their MobilizedBodyIndeces
    //for (SimTK::Compound::AtomIndex aIx(bAtomList[0].atomIdex); aIx < getNumAtoms(); ++aIx){
    for (unsigned int i = 0; i < getNumAtoms(); ++i){
	SimTK::Compound::AtomIndex aIx = (bAtomList[i]).atomIndex;
        std::cout << "Topology::loadMobodsRelatedMaps atomIndex for atom " << i << " = " << (bAtomList[i]).atomIndex << std::endl;
	
        // Map mbx2aIx contains only atoms at the origin of mobods
        SimTK::MobilizedBodyIndex mbx = getAtomMobilizedBodyIndex(aIx);
        //std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex >
        //        pairToBeInserted(mbx, aIx);
        mbx2aIx.insert(
                std::pair<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>
                (mbx, aIx));

        // Map aIx is redundant in MobilizedBodyIndeces
        aIx2mbx.insert(
                std::pair<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>
                (aIx, mbx));
    }

}


/** Print maps **/
void Topology::printMaps(void)
{
    std::cout << "Topology " << name << " maps " << std::endl;
    std::cout << "mbx2aIx:" << std::endl;
    map<SimTK::MobilizedBodyIndex, SimTK::Compound::AtomIndex>::const_iterator mbx2aIxIt;
    for(mbx2aIxIt = mbx2aIx.begin();
       mbx2aIxIt != mbx2aIx.end(); ++mbx2aIxIt)
    {
        std::cout << "mbx " << mbx2aIxIt->first
            << " atomIndex " << mbx2aIxIt->second << std::endl;
    }
    std::cout << "aIx2mbx:" << std::endl;
    map<SimTK::Compound::AtomIndex, SimTK::MobilizedBodyIndex>::const_iterator aIx2mbxIt;
    for(aIx2mbxIt = aIx2mbx.begin();
       aIx2mbxIt != aIx2mbx.end(); ++aIx2mbxIt)
    {
        std::cout << "atomIndex " << aIx2mbxIt->first
            << " mbx " << aIx2mbxIt->second << std::endl;
    }

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

}

/** Write a pdb with bAtomList coordinates and inNames **/
void Topology::writeAtomListPdb(std::string dirname, std::string prefix,
                                std::string sufix, int maxNofDigits, int index) const
{
    int nofDigits = (int) floor(log10(index));
    std::string zeros("");
    if(maxNofDigits > nofDigits){
        for(int i = 0; i < (maxNofDigits - nofDigits); i++){
            zeros += std::string("0");
        }
    }
    std::stringstream sstream;
    sstream << dirname << "/" << prefix << zeros << std::to_string(index) << sufix;
    string ofilename = sstream.str();
    //std::cout << "Topology writePdb to " << ofilename << std::endl;

    FILE *oF = fopen (ofilename.c_str(),"w");
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
    Topology::compoundIndex = compoundIndex;
}


