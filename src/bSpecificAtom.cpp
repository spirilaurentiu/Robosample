#include "bSpecificAtom.hpp"

using namespace SimTK;

/****
 * bSpecificAtom
 ****/
bSpecificAtom::bSpecificAtom(){
    nbonds = 0;
    freebonds = 0;
    number = 0;
    bZeroCharArray(name, 5);
    bZeroCharArray(inName, 5);
    bZeroCharArray(fftype, 20);
    //bZeroCharArray(biotype, 20);
    x = -999;
    y = -999;
    z = -999;
    visited = 0;
    charge = 0.0;
}

/** We do not deallocate bAtomType here. We leave this task to the Topology
class that owns this atom in order to allow the number of atoms connected to
this to change (e.g. semi-grand canonical ensemble). **/
bSpecificAtom::~bSpecificAtom(){
}

// Init function
void bSpecificAtom::Zero(void){
    nbonds = 0;
    freebonds = 0;
    number = 0;
    bZeroCharArray(name, 5);
    bZeroCharArray(inName, 5);
    bZeroCharArray(fftype, 20);
    //bZeroCharArray(biotype, 20);
    x = -999;
    y = -999;
    z = -999;
    visited = 0;
    bAtomType = NULL;
}

void bSpecificAtom::Print(void)
{
    std::cout<<"bSpecificAtom Print: nbonds "<<nbonds<<" freebonds "<<freebonds<<" name "<<name<<" inName "<<inName
        <<" number "<<number<<" atomIndex  "<<atomIndex<<" elem "<<elem<<" atomicNumber "<<atomicNumber<<" x "<<x<<" y "<< y<<" z "<<z
        <<" mass "<<mass<<" vdwRadius  "<<vdwRadius<<" LJWellDepth  "<<LJWellDepth<<" fftype "<<fftype
        <<" atomClassIndex  "<<atomClassIndex<<" biotype (useless) "<< biotype << " biotypeIndex " << biotypeIndex 
        //<< " bAtomType "<< bAtomType 
        <<" charge "<<charge<<" mobile "<<mobile<<" visited "<<visited<<std::endl;

    std::cout << "Neighbors:";
    for(std::vector<bSpecificAtom *>::iterator it = neighbors.begin();
    it != neighbors.end(); ++it){
        std::cout << ' ' << (*it)->number;
    }
    std::cout << std::endl;

    std::cout << "Bonds Involved:" << std::endl;
    for(std::vector<bBond *>::iterator it = bondsInvolved.begin();
    it != bondsInvolved.end(); ++it){
        (*it)->Print();
    }

    // Residue info
    //std::string residueName;
    //long int residueIndex;
    //std::string chain;
    //int moleculeIndex;


}


// bSpecificAtom Interface

// Returns the # of bonds this atom is involved in
int bSpecificAtom::getNBonds(void){
    return this->nbonds;
}

// Returns the # of available bonds left for this atom
int bSpecificAtom::getFreebonds(void)
{
    return this->freebonds;
}

// Returns atom name
std::string bSpecificAtom::getName(void)
{
    return this->name;
}

// Returns atom initial name as received from the reader
std::string bSpecificAtom::getInName(void)
{
    return this->inName;
}

// Returns atom number
int bSpecificAtom::getNumber(void)
{
    return this->number;
}

// Returns atom element
char bSpecificAtom::getElem(void)
{
    return this->elem;
}

// Returns X Cartesian coordinate
SimTK::Real bSpecificAtom::getX(void) const
{
    return this->x;
}

// Returns Y Cartesian coordinate
SimTK::Real bSpecificAtom::getY(void) const
{
    return this->y;
}

// Returns Z Cartesian coordinate
SimTK::Real bSpecificAtom::getZ(void) const
{
    return this->z;
}

// Returns force field atom type
std::string bSpecificAtom::getFftype(void)
{
    return this->fftype;
}

// Returns atom Molmodel biotype name
std::string bSpecificAtom::getBiotype(void)
{
    return this->biotype;
}

// Returns Molmodel atom type as SimTK::Compound::SingleAtom *
SimTK::Compound::SingleAtom * bSpecificAtom::getBAtomType(void)
{
    return this->bAtomType;
}

// Get atom's index as held in Compound
SimTK::Compound::AtomIndex bSpecificAtom::getCompoundAtomIndex(void)
{
    return this->atomIndex;
}

//
SimTK::Real bSpecificAtom::getCharge(void){
    assert(!"Not implemented");
}

//
int bSpecificAtom::getIsMobile(void)
{
    assert(!"Not implemented");
}

//
int bSpecificAtom::getIsVisited(void)
{
    assert(!"Not implemented");
}

void bSpecificAtom::setNbonds(int)
{
    assert(!"Not implemented");
}

void bSpecificAtom::setFreebonds(int){assert(!"Not implemented");}

// Set atom unique name
void bSpecificAtom::setName(std::string inpName){
    strncpy(this->name, inpName.c_str(), 5);
}

// Set initial name

void bSpecificAtom::setInName(std::string inpInName){
    strncpy(this->inName, inpInName.c_str(), 4);
}

// Set number
void bSpecificAtom::setNumber(int inpNumber){
    this->number = inpNumber;
}

// Set element
void bSpecificAtom::setElem(char inpElem){
    this->elem = inpElem;
}

// Set the X coordinate
void bSpecificAtom::setX(SimTK::Real inpX){
    this->x = inpX;
}

// Set the Y coordinate
void bSpecificAtom::setY(SimTK::Real inpY){
    this->y = inpY;
}

// Set the Z coordinate
void bSpecificAtom::setZ(SimTK::Real inpZ){
    this->z = inpZ;
}

// Get atomic mass
SimTK::mdunits::Mass bSpecificAtom::getMass(void)
{
    return this->mass;
}

// Set atomic mass
void bSpecificAtom::setMass(SimTK::Real inpMass)
{
    this->mass = inpMass;
}

// Set force field atom type
void bSpecificAtom::setFftype(std::string inpFftype){
    //strncpy(this->fftype, inpFftype.c_str(), 20);
    sprintf(this->fftype, "%s", inpFftype.c_str());
}

// Set atom Biotype name - dangerous
void bSpecificAtom::setBiotype(std::string inpBiotype)
{
    biotype = inpBiotype;
}

// Set atom Biotype name - dangerous
void bSpecificAtom::setBiotype(const char * inpBiotype)
{
    biotype = inpBiotype;
}

// Set 
void bSpecificAtom::setBAtomType(SimTK::Compound::SingleAtom *){assert(!"Not implemented");}

void bSpecificAtom::setCompoundAtomIndex(SimTK::Compound::AtomIndex inpAtomIndex)
{
    this->atomIndex = inpAtomIndex;
}

// Set charge
void bSpecificAtom::setCharge(SimTK::Real inpCharge){
    this->charge = inpCharge;
}


void bSpecificAtom::setIsMobile(int){assert(!"Not implemented");}

/** Set the number of times this atom was visited during the construction of
 * the graph **/
void bSpecificAtom::setVisited(int argVisited)
{
    this->visited = argVisited;
}

// Get the atom class index
DuMM::AtomClassIndex bSpecificAtom::getAtomClassIndex(void)
{
    return this->atomClassIndex;
}

// Set the atom class index
void bSpecificAtom::setAtomClassIndex(DuMM::AtomClassIndex inpAtomClassIndex)
{
    this->atomClassIndex = inpAtomClassIndex;
}

// Get the atomic number
int bSpecificAtom::getAtomicNumber(void)
{
    return this->atomicNumber;
}

// Set the atomic number
void bSpecificAtom::setAtomicNumber(int inpAtomicNumber)
{
    this->atomicNumber = inpAtomicNumber;
}

//
void bSpecificAtom::setVdwRadius(SimTK::Real inpVdwRadius)
{
    this->vdwRadius = inpVdwRadius;
}

//
SimTK::Real bSpecificAtom::getVdwRadius(void)
{
    return this->vdwRadius;
}

//
void bSpecificAtom::setLJWellDepth(SimTK::Real inpLJWellDepth)
{
    this->LJWellDepth = inpLJWellDepth;
}

//
SimTK::Real bSpecificAtom::getLJWellDepth(void)
{
    return this->LJWellDepth;
}

// Get BiotypeIndex
SimTK::BiotypeIndex bSpecificAtom::getBiotypeIndex(void)
{
    return biotypeIndex;
}

// Set BiotypeIndex
void bSpecificAtom::setBiotypeIndex(SimTK::BiotypeIndex argBiotypeIndex)
{
    this->biotypeIndex = argBiotypeIndex;
}

// Add an atom pointer to the vector of atom neighbors
void bSpecificAtom::addNeighbor(bSpecificAtom *someNeighbor)
{
    neighbors.push_back(someNeighbor);
}

// Add a bond this atom is involved in to the vector of bonds
void bSpecificAtom::addBond(bBond *someBond)
{
    bondsInvolved.push_back(someBond);
}

const DuMM::ChargedAtomTypeIndex bSpecificAtom::getChargedAtomTypeIndex() const {
    return chargedAtomTypeIndex;
}

void bSpecificAtom::setChargedAtomTypeIndex(const SimTK::DuMM::ChargedAtomTypeIndex cAIx) {
    bSpecificAtom::chargedAtomTypeIndex = cAIx;
}

/********************
 *     FUNCTIONS
 * ******************/
int bAtomAssign(MolAtom *dest, const bSpecificAtom *src)
{
  dest->name = new char[5];
  strncpy(dest->name, src->name, 4);

  dest->num = src->number;

  switch(src->elem){
    case 'C': dest->type = MOL_ATOM_ELEMENT_CARBON; break;
    case 'H': dest->type = MOL_ATOM_ELEMENT_HYDROGEN; break;
    case 'N': dest->type = MOL_ATOM_ELEMENT_NITROGEN; break;
    case 'O': dest->type = MOL_ATOM_ELEMENT_OXYGEN; break;
    case 'S': dest->type = MOL_ATOM_ELEMENT_SULFUR; break;
    case 'P': dest->type = MOL_ATOM_ELEMENT_PHOSPHORUS; break;
    default: dest->type = MOL_ATOM_ELEMENT_UNKNOWN; break;
  }

  dest->pos[0] = src->x;
  dest->pos[1] = src->y;
  dest->pos[2] = src->z;

  return 0;
}




