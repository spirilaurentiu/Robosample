#ifndef __BBOND__
#define __BBOND__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the bMoleculeReader class and additional heloer classes
 **/

#include "bgeneral.hpp"
#include "Robo.hpp"
#include "Simbody.h"
#include "Molmodel.h"

//==============================================================================
//                           CLASS intpair
//==============================================================================
/** 
 * Intpair Class is a two int vector used for connectivity definition in MoleculeReader.
**/
class intpair{
 public:
  int i; int j; // These will correspond to bSpecificAtom.number
  intpair();
  intpair(int inI, int inJ);
  ~intpair();

  bool operator==(const intpair *other);
  bool operator!=(const intpair *other);
  bool isTheSameAs(const intpair *other);
  void swap(void);
  void dump(void);
  std::string getString(void);
};

//==============================================================================
//                           CLASS Bond
//==============================================================================
/** 
 * Bond Class used for connectivity definition in MoleculeReader.
**/
class bBond : public intpair{
 private:
  int visited;
  int inring;
  int ring_closing;
  //int rigid;
  SimTK::BondMobility::Mobility mobility;
  int ring_no;
  SimTK::Compound::BondIndex bondIndex;
  bool _isFirst;
  int myindex;

 public:
  bBond(void);
  bBond(int a, int b);
  ~bBond(void);

  bool isInRing(void);
  bool isRingClosing(void);
  //bool isRigid(void);
  SimTK::BondMobility::Mobility getBondMobility(void);
  int ringNo(void);

  void setInRing(void);
  void setAsRingClosing(void);
  //void setAsRigid(void);
  void setBondMobility(SimTK::BondMobility::Mobility someMobility);
  void setRingNo(int rn);

  SimTK::Compound::BondIndex getBondIndex(void);
  void setBondIndex(SimTK::Compound::BondIndex otherIx);

  // Gmolmodel indices (prmtop)
  void setIndex(int);
  int getIndex(void);

  void Print(void);

  bool isFirst(void);
  void setAsFirst(void);

  int isThisMe(int argFirst, int argSecond);
  void setVisited(int);
  int isVisited(void);
};



#endif  //__BBOND__


