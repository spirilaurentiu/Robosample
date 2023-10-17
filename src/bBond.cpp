#include "bBond.hpp"

/*
#ifndef DEBUG_LEVEL01
#define DEBUG_LEVEL01
#endif

#ifndef DEBUG_LEVEL02
#define DEBUG_LEVEL02
#endif
*/

/****
 * intpair
 ****/
intpair::intpair(int inI, int inJ) : i(inI), j(inJ) {
}

bool intpair::operator==(const intpair& rhs){
  return (
    ((this->i == rhs.i) && (this->j == rhs.j)) ||
    ((this->i == rhs.j) && (this->j == rhs.i))
  );
}

bool intpair::operator!=(const intpair& rhs){
  return (
    ((this->i != rhs.i) || (this->j != rhs.j)) &&
    ((this->i != rhs.j) || (this->j != rhs.i))
  );
}

bool intpair::isTheSameAs(const intpair& rhs){
  return (
    ((this->i == rhs.i) && (this->j == rhs.j)) ||
    ((this->i == rhs.j) && (this->j == rhs.i))
  );
}

void intpair::swap(){
  std::swap(i, j);
}

void intpair::dump(){
  std::cout<<i<<' '<<j<<std::endl;
}

std::string intpair::getString(){
  std::stringstream ret;
  ret<<i<<' '<<j;
  return ret.str();
}

/********************************/


/****
 * bBond
 ****/
bBond::bBond(int a, int b) : intpair(a, b){
}

bool bBond::isInRing() const {
  return this->inring;
}

bool bBond::isRingClosing() const {
  return this->ring_closing;
}

/*
bool bBond::isRigid() const {
  if(this->rigid == 1)
    return true;
  else
    return false;
}
*/

SimTK::BondMobility::Mobility bBond::getBondMobility(int which) const {
  return mobilities[which];
}

void bBond::addBondMobility(SimTK::BondMobility::Mobility argmobility)
{
	mobilities.emplace_back(argmobility);
}

void bBond::setBondMobility(SimTK::BondMobility::Mobility argmobility, int which)
{
    mobilities[which] = argmobility;
}

void bBond::updBondMobility(SimTK::BondMobility::Mobility argmobility, int which)
{
    mobilities[which] = argmobility;
}

float bBond::getUScaleFactor(int which) const
{
  
  if(which >= uScaleFactors.size()){
    std::cout << "[WARNING] bBond::getUScaleFactor() " << which << " not found. Returning 1" << std::endl;
    return 1.0;
  }else{
    return uScaleFactors[which];
  }

}

void bBond::addUScaleFactor(float argUScaleFactor)
{
	uScaleFactors.emplace_back(argUScaleFactor);
}

void bBond::setUScaleFactor(int which, float argUScaleFactor)
{
	uScaleFactors[which] = argUScaleFactor;
}

void bBond::updUScaleFactor(int which, float argUScaleFactor)
{
	uScaleFactors[which] = argUScaleFactor;
}


int bBond::ringNo() const{
  return this->ring_no;
}

void bBond::setInRing(){
  this->inring = true;
}

void bBond::setAsRingClosing(){
  this->ring_closing = true;
}
/*
void bBond::setAsRigid(){
  this->rigid = 1;
}
*/
void bBond::setRingNo(int rn){
  this->ring_no = rn;
}

SimTK::Compound::BondIndex bBond::getBondIndex() const {
  return bondIndex;
}
  
void bBond::setBondIndex(SimTK::Compound::BondIndex otherIx){
  this->bondIndex = otherIx;
}

// Print bond variables
void bBond::Print(int whichWorld)
{
    std::cout << "bBond Print: i " << i << " j " << j << " mobility " << mobilities[whichWorld]
        << " inring " << inring << " ring_no " << ring_no
        << " ring_closing " << ring_closing 
        << " visited " << visited 
        << " uScaleFactor[0] " << uScaleFactors[0]
	<< std::endl;
}

// Return true if this is set as the first bond in Compound
bool bBond::isFirst() const {
    return _isFirst;
}

// Set this bond as the first one in a Compound
void bBond::setAsFirst()
{
    _isFirst = true;
}

// Check if this bond contains the arguments. If they are swapped return -1
int bBond::isThisMe(int argFirst, int argSecond) const
{
    if( (argFirst == this->i) && (argSecond == this->j) ){
        return 1;
    }
    else if( (argSecond == this->i) && (argFirst == this->j) ){
        return -1;
    }
    else{
        return 0;
    }
    return 0;
}

// Set the number of times this bond was visited
void bBond::setVisited(int argVisited)
{
    this->visited = argVisited;
}

// Return the number of times this bond was visited
int bBond::isVisited() const
{
	return this->visited;
}

// Gmolmodel indices (prmtop)
void bBond::setIndex(int someOtherIndex)
{
    myindex = someOtherIndex;
}

int bBond::getIndex() const
{
    return myindex;
}






