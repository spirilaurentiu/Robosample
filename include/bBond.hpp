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
// #include "Simbody.h"
#include "Molmodel.h"

//==============================================================================
//                           CLASS intpair
//==============================================================================
/** 
 * Intpair Class is a two int vector used for connectivity definition in MoleculeReader.
**/
class intpair{
	public:
		intpair() = default;
		intpair(int inI, int inJ);

		bool operator==(const intpair& rhs);
		bool operator!=(const intpair& rhs);
		bool isTheSameAs(const intpair& rhs);
		void swap();
		void dump();
		std::string getString();

		// These will correspond to bSpecificAtom.number
		int i = 0;
		int j = 0;
};

//==============================================================================
//                           CLASS Bond
//==============================================================================
/** 
 * Bond Class used for connectivity definition in MoleculeReader.
**/
class bBond : public intpair{
	private:
		std::vector<SimTK::BondMobility::Mobility> mobilities;
		// std::vector<float> uScaleFactors = std::vector<float>(1, 1.0);
		std::vector<float> uScaleFactors = { 1.0f };
		SimTK::Compound::BondIndex bondIndex = SimTK::Compound::BondIndex(99999999);
		int ring_no = 0;
		int myindex = -1;
		//int rigid;
		bool visited = false;
		bool inring = false;
		bool ring_closing = false;
		bool _isFirst = false;

	public:
		bBond() = default;
		bBond(int a, int b);

		bool isInRing() const;
		bool isRingClosing() const;
		//bool isRigid() const;
		SimTK::BondMobility::Mobility getBondMobility(int whichWorld) const;
		int ringNo() const;

		void setInRing();
		void setAsRingClosing();
		//void setAsRigid();

		void addBondMobility(SimTK::BondMobility::Mobility someMobility);	
		void setBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
		void updBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
		void setRingNo(int rn);

		SimTK::Compound::BondIndex getBondIndex() const;
		void setBondIndex(SimTK::Compound::BondIndex otherIx);

		// Gmolmodel indices (prmtop)
		void setIndex(int);
		int getIndex() const;

		void Print(int whichWorld);

		bool isFirst() const;
		void setAsFirst();

		int isThisMe(int argFirst, int argSecond) const;

		void setVisited(int);
		int isVisited() const;

		float getUScaleFactor(int) const;
		void addUScaleFactor(float);
		void setUScaleFactor(int, float);
		void updUScaleFactor(int, float);

};



#endif  //__BBOND__


