#include "bSpecificAtom.hpp"

using namespace SimTK;

void bSpecificAtom::setAtomCompoundType(const SimTK::Element &element) {
    const std::string& atomName = getName();
    const int numAtomBonds = getNBonds();
    
    switch (numAtomBonds) {
        case 0:
            compoundSingleAtom = new SimTK::Compound::SingleAtom(atomName, element);
            break;
        case 1:
            compoundSingleAtom = new SimTK::UnivalentAtom(atomName, element);
            break;
        case 2:
            compoundSingleAtom = new SimTK::BivalentAtom(atomName, element);
            break;
        case 3:
            compoundSingleAtom = new SimTK::TrivalentAtom(atomName, element); // 120 rad
            break;
        case 4:
            compoundSingleAtom = new SimTK::QuadrivalentAtom(atomName, element); // 109.47 rad
            break;
        default:
            std::cerr << "[ERROR] Atom " << atomName << " has " << numAtomBonds << " bonds, which is not supported." << std::endl;
            throw std::exception();
    }

    // q: how do i test that this works?
    // a: use the following code
    // SimTK::CompoundSystem system;
    // SimTK::DuMMForceFieldSubsystem dumm(system);
    // SimTK::DuMM::ChargedAtomTypeIndex chargedAtomTypeIndex = dumm.defineChargedAtomType(compoundSingleAtom->getAtomName(), element, 0.0, 0.0);
    // cout << "Charged atom type index: " << chargedAtomTypeIndex << endl;
    // cout << "Charged atom type name: " << dumm.getChargedAtomTypeName(chargedAtomTypeIndex) << endl;
    // cout << "Charged atom type element: " << dumm.getChargedAtomTypeElement(chargedAtomTypeIndex) << endl;
    // cout << "Charged atom type charge: " << dumm.getChargedAtomTypeCharge(chargedAtomTypeIndex) << endl;
    // cout << "Charged atom type radius: " << dumm.getChargedAtomTypeRadius(chargedAtomTypeIndex) << endl;
    // cout << "Charged atom type well depth: " << dumm.getChargedAtomTypeWellDepth(chargedAtomTypeIndex) << endl;
    // cout << "Charged atom type index: " << dumm.getChargedAtomTypeIndex(compoundSingleAtom->getAtomName(), element) << endl;

    // compoundSingleAtom->setAtomClassIndex(dumm.getAtomClassIndex(compoundSingleAtom->getAtomName(), element));
    // compoundSingleAtom->setChargedAtomTypeIndex(dumm.getChargedAtomTypeIndex(compoundSingleAtom->getAtomName(), element));
    // compoundSingleAtom->setCompoundName(compoundSingleAtom->getAtomName());
    // compoundSingleAtom->setDefaultInboardBondLength(0.19);

    compoundSingleAtom->setCompoundName("SingleAtom");

    if (numAtomBonds > 0) {
        compoundSingleAtom->setDefaultInboardBondLength(0.19);
    }
}

void bSpecificAtom::destroy() {
    if (compoundSingleAtom != nullptr) {
        delete compoundSingleAtom;
        compoundSingleAtom = nullptr;
    }
}

void bSpecificAtom::Print(int whichWorld)
{
    std::cout<<"bSpecificAtom Print: nbonds "<<nbonds<<" freebonds "<<freebonds<<" name "<< name <<" inName "<< inName
        <<" number "<<number<<" atomIndex  "<<compoundAtomIndex<<" elem "<<elem<<" atomicNumber "<<atomicNumber
        << " x " << x << " y " << y << " z " << z
        <<" Cartesians " << Cartesians[0] << " " << Cartesians[1] << " " << Cartesians[2] 
        <<" mass "<<mass<<" vdwRadius  "<<vdwRadius<<" LJWellDepth  "<<LJWellDepth<<" fftype "<< fftype
        <<" atomClassIndex  "<<dummAtomClassIndex<<" biotype "<< biotype << " biotypeIndex " << biotypeIndex 
        //<< " bAtomType "<< bAtomType 
        <<" charge "<<charge
        //<< " mobile "<<mobile
        <<" visited "<<visited
        <<" moleculeIndex "<<moleculeIndex
        <<std::endl;


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

// Returns the Cartesian coordinates
const double* bSpecificAtom::getCartesians() const
{
    return Cartesians;
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

const SimTK::Compound::SingleAtom& bSpecificAtom::getSingleAtom() const {
    return *compoundSingleAtom;
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
// int bSpecificAtom::getIsMobile() const
// {
//     assert(!"Not implemented");
//     throw std::exception();
    
//     return std::numeric_limits<int>::min();
// }

//
bool bSpecificAtom::wasVisited() const
{
    if( this->visited > 0){
        return true;
    }else{
        return false;
    }
}

int bSpecificAtom::getVisited() const
{
    return this->visited;
}

/** Set the number of times this atom was visited during the construction of
 * the graph **/
void bSpecificAtom::setVisited(const int& argVisited)
{
    this->visited = argVisited;
}

void bSpecificAtom::setNbonds(const int& nbondsArg)
{
     this->nbonds = nbondsArg;
}

void bSpecificAtom::setFreebonds(const int& freebondsArg)
{
     this->freebonds = freebondsArg;
}

void bSpecificAtom::incrFreebonds(void)
{
     (this->freebonds)++;
}

void bSpecificAtom::decrFreebonds(void)
{
     (this->freebonds)--;
}

// Set atom unique name
void bSpecificAtom::generateName(int nameCounter) {
    std::string string_name;
    std::string aStr, bStr, cStr, dStr;
    int a=65, b=65, c=65, d=65;
    int aRest=0, bRest=0, cRest=0;

    a = int(nameCounter / std::pow(25, 3));
    aStr = (char)(a + 65);
    aRest = nameCounter % int(std::pow(25, 3));

    b = int(aRest / std::pow(25, 2));
    bStr = (char)(b + 65);
    bRest = aRest % int(std::pow(25, 2));

    c = int(bRest / std::pow(25, 1));
    cStr = (char)(c + 65);
    cRest = bRest % int(std::pow(25, 1));

    d = int(cRest / std::pow(25, 0));
    dStr = (char)(d + 65);

    name = aStr + bStr + cStr + dStr;
}

// Set initial name
void bSpecificAtom::setInName(const std::string& inName){
    this->inName = inName;
}

// Set number
void bSpecificAtom::setNumber(int inpNumber){
    this->number = inpNumber;
}

// Set element
void bSpecificAtom::setElem(const std::string& elem){
    this->elem = elem;
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

// Set the Cartesian coordinates
void bSpecificAtom::setCartesians(
    SimTK::Real inpX, SimTK::Real inpY, SimTK::Real inpZ)
{
    Cartesians[0] = inpX;
    Cartesians[1] = inpY;
    Cartesians[2] = inpZ;
}

// // Set the Cartesian coordinates
// void bSpecificAtom::setCartesians(
//     double* inpXYZ)
// {
//     assert( sizeof(inpXYZ) / sizeof(inpXYZ[0]) == 3 );
//     Cartesians[0] = inpXYZ[0];
//     Cartesians[1] = inpXYZ[1];
//     Cartesians[2] = inpXYZ[2];

// }

// Set the Cartesian coordinates
void bSpecificAtom::setCartesians(
    SimTK::Vec3 inpXYZ)
{
    Cartesians[0] = inpXYZ[0];
    Cartesians[1] = inpXYZ[1];
    Cartesians[2] = inpXYZ[2];
}

// Set the Cartesian coordinates
void bSpecificAtom::setCartesians(
    OpenMM::Vec3 inpXYZ)
{
    Cartesians[0] = inpXYZ[0];
    Cartesians[1] = inpXYZ[1];
    Cartesians[2] = inpXYZ[2];
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
void bSpecificAtom::setFfType(const std::string& fftype)
{
    this->fftype = fftype;
}

// Set atom Biotype name - dangerous
void bSpecificAtom::setBiotype(const std::string& biotype)
{
    this->biotype = biotype;
}

void bSpecificAtom::setCompoundAtomIndex(SimTK::Compound::AtomIndex inpAtomIndex)
{
    this->compoundAtomIndex = inpAtomIndex;
}

// Set charge
void bSpecificAtom::setCharge(SimTK::Real inpCharge){
    this->charge = inpCharge;
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

// /********************
//  *     FUNCTIONS
//  * ******************/
// int bAtomAssign(MolAtom *dest, const bSpecificAtom *src)
// {
// /*     dest->name = new char[5]; // TODO who deletes this?

//     std::copy(dest->name, dest->name + 5, src->getName().begin());

//     dest->num = src->getNumber();

//     switch(src->getElem()){
//         case 'C': dest->type = MOL_ATOM_ELEMENT_CARBON; break;
//         case 'H': dest->type = MOL_ATOM_ELEMENT_HYDROGEN; break;
//         case 'N': dest->type = MOL_ATOM_ELEMENT_NITROGEN; break;
//         case 'O': dest->type = MOL_ATOM_ELEMENT_OXYGEN; break;
//         case 'S': dest->type = MOL_ATOM_ELEMENT_SULFUR; break;
//         case 'P': dest->type = MOL_ATOM_ELEMENT_PHOSPHORUS; break;
//         default: dest->type = MOL_ATOM_ELEMENT_UNKNOWN; break;
//     }

//     // TODO they use float, we use double
//     dest->pos[0] = static_cast<float>(src->getX());
//     dest->pos[1] = static_cast<float>(src->getY());
//     dest->pos[2] = static_cast<float>(src->getZ()); */

//     return 0;
// }

void bSpecificAtom::setCompoundName(const SimTK::Compound::Name& name) {
    compoundSingleAtom->setCompoundName(name);
}

// Getter and setter for the residueName property
std::string bSpecificAtom::getResidueName() const {
    return residueName;
}

void bSpecificAtom::setResidueName(const std::string& value) {
    residueName = value; // Directly assign the new value
}


const int bSpecificAtom::getMoleculeIndex() const
{
    return moleculeIndex;
}

void bSpecificAtom::setMoleculeIndex(int molIxArg)
{
    this->moleculeIndex = molIxArg;
}

void bSpecificAtom::setParentNumber(int n) {
    parentNumber = n;
}

int bSpecificAtom::getParentNumber() const {
    return parentNumber;
};

void bSpecificAtom::setName(const std::string& name) {
    this->name = name;
}
