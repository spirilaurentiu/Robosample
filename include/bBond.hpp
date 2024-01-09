#pragma once

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file
 * This defines the bMoleculeReader class and additional heloer classes
 **/

// #include "bgeneral.hpp"
// #include "Robo.hpp"
// #include "Simbody.h"
#include "Molmodel.h"

// //==============================================================================
// //                           CLASS intpair
// //==============================================================================
// /** 
//  * Intpair Class is a two int vector used for connectivity definition in MoleculeReader.
// **/
// class intpair{
// 	public:
// 		intpair() = default;
// 		intpair(int inI, int inJ);

// 		bool operator==(const intpair& rhs);
// 		bool operator!=(const intpair& rhs);
// 		bool isTheSameAs(const intpair& rhs);
// 		void swap();
// 		void dump();
// 		std::string getString();

// 		// These will correspond to bSpecificAtom.number
// 		int i = 0;
// 		int j = 0;
// };

//==============================================================================
//                           CLASS Bond
//==============================================================================
/** 
 * Bond Class used for connectivity definition in MoleculeReader.
**/
class bBond /* : public intpair */ {
	private:
		std::vector<SimTK::BondMobility::Mobility> mobilities;
		// std::vector<float> uScaleFactors = std::vector<float>(1, 1.0);
		std::vector<float> uScaleFactors = { 1.0f };
		SimTK::Compound::BondIndex bondIndex = std::numeric_limits<SimTK::Compound::BondIndex>::max();
		// int ring_no = 0;
		int myindex = std::numeric_limits<int>::min();

		SimTK::Real forceK = std::numeric_limits<SimTK::Real>::min();
		SimTK::Real forceEquil = std::numeric_limits<SimTK::Real>::min();

	public:
		// These will correspond to bSpecificAtom.number
		int i = std::numeric_limits<int>::min();
		int j = std::numeric_limits<int>::min();

	private:
		//int rigid;
		bool visited = false;
		// bool inring = false;
		bool ring_closing = false;
		// bool _isFirst = false;

	public:
		// bBond() = default;
		// bBond(int a, int b);

		void setForceK(SimTK::Real forceK);
		SimTK::Real getForceK() const;

		void setForceEquil(SimTK::Real forceEquil);
		SimTK::Real getForceEquil() const;

		// bool isInRing() const;
		bool isRingClosing() const;
		//bool isRigid() const;
		SimTK::BondMobility::Mobility getBondMobility(int whichWorld) const;
		// int ringNo() const;

		// void setInRing();
		void setAsRingClosing();
		//void setAsRigid();

		void addBondMobility(SimTK::BondMobility::Mobility someMobility);	
		void setBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
		void updBondMobility(SimTK::BondMobility::Mobility someMobility, int whichWorld);
		// void setRingNo(int rn);

		SimTK::Compound::BondIndex getBondIndex() const;
		void setBondIndex(SimTK::Compound::BondIndex otherIx);

		// Gmolmodel indices (prmtop)
		void setIndex(int);
		int getIndex() const;

		void Print(int whichWorld);

		// bool isFirst() const;
		// void setAsFirst();

		int isThisMe(int argFirst, int argSecond) const;

		void setVisited(int);
		int isVisited() const;

		float getUScaleFactor(int) const;
		void addUScaleFactor(float);
		void setUScaleFactor(int, float);
		// void updUScaleFactor(int, float);
};

struct DUMM_ANGLE {
	int first = std::numeric_limits<int>::min();
	int second = std::numeric_limits<int>::min();
	int third = std::numeric_limits<int>::min();

	SimTK::Real k = std::numeric_limits<SimTK::Real>::min();
	SimTK::Real equil = std::numeric_limits<SimTK::Real>::min();
};

struct DUMM_TORSION {
	int first = std::numeric_limits<int>::min();
	int second = std::numeric_limits<int>::min();
	int third = std::numeric_limits<int>::min();
	int fourth = std::numeric_limits<int>::min();
	bool improper = false;

	// How many impropers with these four indices are present here
	int num = 0; 

	// These values are filled according to num (see above)
	int period[4] {};
	SimTK::Real k[4] {};
	SimTK::Real phase[4] {};
};
