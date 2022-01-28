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
	void swap();
	void dump();
	std::string getString();
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
	std::vector<SimTK::BondMobility::Mobility> mobilities;
	int ring_no;
	SimTK::Compound::BondIndex bondIndex;
	bool _isFirst;
	int myindex;
	std::vector<float> uScaleFactors;

 public:
	bBond();
	bBond(int a, int b);
	~bBond();

	bool isInRing();
	bool isRingClosing();
	//bool isRigid();
	SimTK::BondMobility::Mobility getBondMobility(int whichWorld) const;
	int ringNo();

	void setInRing();
	void setAsRingClosing();
	//void setAsRigid();

	void addBondMobility(SimTK::BondMobility::Mobility someMobility);	
	void setBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
	void updBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
	void setRingNo(int rn);

	SimTK::Compound::BondIndex getBondIndex();
	void setBondIndex(SimTK::Compound::BondIndex otherIx);

	// Gmolmodel indices (prmtop)
	void setIndex(int);
	int getIndex();

	void Print(int whichWorld);

	bool isFirst();
	void setAsFirst();

	int isThisMe(int argFirst, int argSecond) const ;

	void setVisited(int);
	int isVisited();

	float getUScaleFactor(int) const;
	void addUScaleFactor(float);
	void setUScaleFactor(int, float);
	void updUScaleFactor(int, float);

};



#endif  //__BBOND__


