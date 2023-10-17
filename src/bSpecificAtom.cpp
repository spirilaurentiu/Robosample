#include "bSpecificAtom.hpp"

using namespace SimTK;

// Init function
void bSpecificAtom::Zero(){
    biotype.clear();
    residueName.clear();
    chain.clear();
    neighbors.clear();
    bondsInvolved.clear();

    std::fill(begin(fftype), end(fftype), '\0');

    compoundSingleAtom = nullptr;

    charge = 0.0;
    residueIndex = std::numeric_limits<long int>::min();

    mass = std::numeric_limits<SimTK::Real>::min();
    vdwRadius = std::numeric_limits<SimTK::Real>::min();
    LJWellDepth = std::numeric_limits<SimTK::Real>::min();

    std::fill(begin(name), end(name), '\0');
    std::fill(begin(inName), end(inName), '\0');

    x = std::numeric_limits<SimTK::Real>::min();
    y = std::numeric_limits<SimTK::Real>::min();
    z = std::numeric_limits<SimTK::Real>::min();

    nbonds = 0;
    freebonds = 0;
    number = 0;
    atomicNumber = 0;
    mobile = 0;
    visited = 0;
    moleculeIndex = std::numeric_limits<int>::min();

    elem = "";
}

void bSpecificAtom::Print(int whichWorld)
{
    // std::cout<<"bSpecificAtom Print: nbonds "<<nbonds<<" freebonds "<<freebonds<<" name "<< name <<" inName "<< inName
    //     <<" number "<<number<<" atomIndex  "<<compoundAtomIndex<<" elem "<<elem<<" atomicNumber "<<atomicNumber<<" x "<<x<<" y "<< y<<" z "<<z
    //     <<" mass "<<mass<<" vdwRadius  "<<vdwRadius<<" LJWellDepth  "<<LJWellDepth<<" fftype "<< fftype
    //     <<" atomClassIndex  "<<dummAtomClassIndex<<" biotype "<< biotype << " biotypeIndex " << biotypeIndex 
    //     //<< " bAtomType "<< bAtomType 
    //     <<" charge "<<charge<<" mobile "<<mobile<<" visited "<<visited<<std::endl;


/*     std::cout << "Neighbors:";
    for(std::vector<bSpecificAtom *>::iterator it = neighbors.begin();
    it != neighbors.end(); ++it){
        std::cout << ' ' << (*it)->number;
    }
    std::cout << std::endl;

    std::cout << "Bonds Involved:" << std::endl;
    for(std::vector<bBond *>::iterator it = bondsInvolved.begin();
    it != bondsInvolved.end(); ++it){
        (*it)->Print(whichWorld);
    } */

    // Residue info
    //std::string residueName;
    //long int residueIndex;
    //std::string chain;
    //int moleculeIndex;
}


// bSpecificAtom Interface

// Returns the # of bonds this atom is involved in
int bSpecificAtom::getNBonds() const
{
    return this->nbonds;
}

// Returns the # of available bonds left for this atom
int bSpecificAtom::getFreebonds() const
{
    return this->freebonds;
}

// Returns atom name
std::string bSpecificAtom::getName() const
{
    return name;
}

// Returns atom initial name as received from the reader
std::string bSpecificAtom::getInName() const
{
    return inName;
}

// Returns atom number
int bSpecificAtom::getNumber() const
{
    return this->number;
}

// Returns atom element
std::string bSpecificAtom::getElem() const
{
    return this->elem;
}

// Returns X Cartesian coordinate
SimTK::Real bSpecificAtom::getX() const
{
    return this->x;
}

// Returns Y Cartesian coordinate
SimTK::Real bSpecificAtom::getY() const
{
    return this->y;
}

// Returns Z Cartesian coordinate
SimTK::Real bSpecificAtom::getZ() const
{
    return this->z;
}

// Returns force field atom type
std::string bSpecificAtom::getFftype() const
{
    return fftype;
}

// Returns atom Molmodel biotype name
std::string bSpecificAtom::getBiotype() const
{
    return this->biotype;
}

// Returns Molmodel atom type as SimTK::Compound::SingleAtom *
SimTK::Compound::SingleAtom * bSpecificAtom::getBAtomType() const
{
    return this->compoundSingleAtom;
}

// Get atom's index as held in Compound
SimTK::Compound::AtomIndex bSpecificAtom:: getCompoundAtomIndex() const
{
    return this->compoundAtomIndex;
}

//
SimTK::Real bSpecificAtom::getCharge() const {
    assert(!"Not implemented");
    throw std::exception();

    return SimTK::NaN;
}

//
int bSpecificAtom::getIsMobile() const
{
    assert(!"Not implemented");
    throw std::exception();
    
    return std::numeric_limits<int>::min();
}

//
int bSpecificAtom::getIsVisited() const
{
    assert(!"Not implemented");
    throw std::exception();
    
    return std::numeric_limits<int>::min();
}

void bSpecificAtom::setNbonds(int)
{
    assert(!"Not implemented");
    throw std::exception();
}

void bSpecificAtom::setFreebonds(int)
{
    assert(!"Not implemented");
    throw std::exception();
}

// Set atom unique name
void bSpecificAtom::setName(std::string inpName){
    std::copy(inpName.begin(), inpName.end(), begin(name));
}

// Set initial name

void bSpecificAtom::setInName(std::string inpInName){
    std::copy(inpInName.begin(), inpInName.end(), begin(inName));
}

// Set number
void bSpecificAtom::setNumber(int inpNumber){
    this->number = inpNumber;
}

// Set element
void bSpecificAtom::setElem(std::string inpElem){
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
SimTK::mdunits::Mass bSpecificAtom::getMass() const
{
    return this->mass;
}

// Set atomic mass
void bSpecificAtom::setMass(SimTK::Real inpMass)
{
    this->mass = inpMass;
}

// Set force field atom type
void bSpecificAtom::setFfType(std::string inpFftype)
{
    std::copy(inpFftype.begin(), inpFftype.end(), begin(fftype));
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
void bSpecificAtom::setBAtomType(SimTK::Compound::SingleAtom *)
{
    assert(!"Not implemented");
    throw std::exception();
}

void bSpecificAtom::setCompoundAtomIndex(SimTK::Compound::AtomIndex inpAtomIndex)
{
    this->compoundAtomIndex = inpAtomIndex;
}

// Set charge
void bSpecificAtom::setCharge(SimTK::Real inpCharge){
    this->charge = inpCharge;
}


void bSpecificAtom::setIsMobile(int)
{
    assert(!"Not implemented");
    throw std::exception();
}

/** Set the number of times this atom was visited during the construction of
 * the graph **/
void bSpecificAtom::setVisited(int argVisited)
{
    this->visited = argVisited;
}

// Get the atom class index
DuMM::AtomClassIndex bSpecificAtom::getDummAtomClassIndex() const
{
    return this->dummAtomClassIndex;
}

// Set the atom class index
void bSpecificAtom::setDummAtomClassIndex(DuMM::AtomClassIndex inpAtomClassIndex)
{
    this->dummAtomClassIndex = inpAtomClassIndex;
}

// Get the atomic number
int bSpecificAtom::getAtomicNumber() const
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
SimTK::Real bSpecificAtom::getVdwRadius() const
{
    return this->vdwRadius;
}

//
void bSpecificAtom::setLJWellDepth(SimTK::Real inpLJWellDepth)
{
    this->LJWellDepth = inpLJWellDepth;
}

//
SimTK::Real bSpecificAtom::getLJWellDepth() const
{
    return this->LJWellDepth;
}

// Get BiotypeIndex
SimTK::BiotypeIndex bSpecificAtom::getBiotypeIndex() const
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

DuMM::ChargedAtomTypeIndex bSpecificAtom::getChargedAtomTypeIndex() const {
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
/*     dest->name = new char[5]; // TODO who deletes this?

    std::copy(dest->name, dest->name + 5, src->getName().begin());

    dest->num = src->getNumber();

    switch(src->getElem()){
        case 'C': dest->type = MOL_ATOM_ELEMENT_CARBON; break;
        case 'H': dest->type = MOL_ATOM_ELEMENT_HYDROGEN; break;
        case 'N': dest->type = MOL_ATOM_ELEMENT_NITROGEN; break;
        case 'O': dest->type = MOL_ATOM_ELEMENT_OXYGEN; break;
        case 'S': dest->type = MOL_ATOM_ELEMENT_SULFUR; break;
        case 'P': dest->type = MOL_ATOM_ELEMENT_PHOSPHORUS; break;
        default: dest->type = MOL_ATOM_ELEMENT_UNKNOWN; break;
    }

    // TODO they use float, we use double
    dest->pos[0] = static_cast<float>(src->getX());
    dest->pos[1] = static_cast<float>(src->getY());
    dest->pos[2] = static_cast<float>(src->getZ()); */

    return 0;
}




