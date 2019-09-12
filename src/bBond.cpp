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
intpair::intpair(){
    i = 0; j = 0;
  }
intpair::intpair(int inI, int inJ){
    this->i = inI;
    this->j = inJ;
}
intpair::~intpair(){}

bool intpair::operator==(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

bool intpair::operator!=(const intpair *other){
  return (
    ((this->i != other->i) || (this->j != other->j)) &&
    ((this->i != other->j) || (this->j != other->i))
  );
}

bool intpair::isTheSameAs(const intpair *other){
  return (
    ((this->i == other->i) && (this->j == other->j)) ||
    ((this->i == other->j) && (this->j == other->i))
  );
}

void intpair::swap(void){
  int inter;
  inter = this->i;
  this->i = this->j;
  this->j = inter;
}

void intpair::dump(void){
  std::cout<<i<<' '<<j<<std::endl;
}

std::string intpair::getString(void){
  std::stringstream ret;
  ret<<i<<' '<<j;
  return ret.str();
}

/********************************/


/****
 * bBond
 ****/
bBond::bBond(void) : intpair(){
  inring = 0;
  //rigid = 0;
  mobility = SimTK::BondMobility::Mobility::Free;
  ring_closing = 0; // later
  ring_no = 0; // later
  _isFirst = false;
  visited = 0;
  bondIndex = SimTK::Compound::BondIndex(99999999);
}

bBond::bBond(int a, int b) : intpair(a, b){
  inring = 0;
  //rigid = 0;
  mobility = SimTK::BondMobility::Mobility::Free;
  ring_closing = 0;
  ring_no = 0; // later
  _isFirst = false;
  visited = 0;
  bondIndex = SimTK::Compound::BondIndex(99999999);
}

bBond::~bBond(void){;}

bool bBond::isInRing(void){
  if(this->inring == 1)
    return true;
  else
    return false;
}

bool bBond::isRingClosing(void){
  if(this->ring_closing == 1)
    return true;
  else
    return false;
}

/*
bool bBond::isRigid(void){
  if(this->rigid == 1)
    return true;
  else
    return false;
}
*/

SimTK::BondMobility::Mobility bBond::getBondMobility(void)
{
    return mobility;
}

void bBond::setBondMobility(SimTK::BondMobility::Mobility argmobility)
{
    mobility = argmobility;
}


int bBond::ringNo(void){
  return this->ring_no;
}

void bBond::setInRing(void){
  this->inring = 1;
}

void bBond::setAsRingClosing(void){
  this->ring_closing = 1;
}
/*
void bBond::setAsRigid(void){
  this->rigid = 1;
}
*/
void bBond::setRingNo(int rn){
  this->ring_no = rn;
}

SimTK::Compound::BondIndex bBond::getBondIndex(void){
  return bondIndex;
}
  
void bBond::setBondIndex(SimTK::Compound::BondIndex otherIx){
  this->bondIndex = otherIx;
}

// Print bond variables
void bBond::Print(void)
{
    std::cout << "i " << i << " j " << j << " mobility " << SimTK::BondMobility::Rigid
        << " inring " << inring << " ring_no " << ring_no
        << " ring_closing " << ring_closing 
        << " visited " << visited << std::endl;
}

// Return true if this is set as the first bond in Compound
bool bBond::isFirst(void)
{
    return _isFirst;
}

// Set this bond as the first one in a Compound
void bBond::setAsFirst(void)
{
    _isFirst = true;
}

// Check if this bond contains the arguments. If they are swapped return -1
int bBond::isThisMe(int argFirst, int argSecond)
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
int bBond::isVisited(void)
{
    return this->visited;
}

// Gmolmodel indices (prmtop)
void bBond::setIndex(int someOtherIndex)
{
    myindex = someOtherIndex;
}

int bBond::getIndex(void)
{
    return myindex;
}

