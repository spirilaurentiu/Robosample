#ifndef __INTERNALCOORDINATES__
#define __INTERNALCOORDINATES__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file This defines the Internal Coordinates class
 **/

#include "bSpecificAtom.hpp"

using namespace SimTK;

//==============================================================================
//                               struct BOND
//==============================================================================
struct BOND {
	int first = -1,
		second = -1;

	bool operator==(const BOND& rhs) const;
	bool operator!=(const BOND& rhs) const;

	void Print(void) const {
		std::cout << first << " " << second;
	}
};
//------------------------------------------------------------------------------


//==============================================================================
//                               struct ANGLE
//==============================================================================
struct ANGLE {
	int first = -1,
		second = -1,
		third = -1;

	bool operator==(const ANGLE& rhs) const;
	bool operator!=(const ANGLE& rhs) const;

	void Print(void) const {
		std::cout << first << " " << second << " " << third;
	}	
};


//==============================================================================
//                               struct TORSION
//==============================================================================
struct TORSION {
	int first = -1,
		second = -1,
		third = -1,
		fourth = -1;

	bool operator==(const TORSION& rhs) const;
	bool operator!=(const TORSION& rhs) const;

	void Print(void) const {
		std::cout << first << " " << second << " " << third << " " << fourth;
	}
};
//------------------------------------------------------------------------------


//==============================================================================
//                               Amber ATOM
//==============================================================================
struct AmberAtom {
	SimTK::Real mass = -1.0;
	int bonds = -1,
		amberId = -1,
		parent = -1;

	AmberAtom(const bSpecificAtom& b);
	AmberAtom(const bSpecificAtom* b);

	bool operator==(const AmberAtom& rhs);
	bool operator!=(const AmberAtom& rhs);
};
//------------------------------------------------------------------------------

//==============================================================================
//                               BAT ATOM
//==============================================================================
struct BAT_ATOM {
	int amberIndex = -1;
	int parent = -1;
};
//------------------------------------------------------------------------------


//==============================================================================
//                           CLASS InternalCoordinates
//==============================================================================
/** 
 * This defines the Internal Coordinates class. Code inspired from MDAnalysis
 * package, analysis.bat algorithm written by Soohaeng Yoo Willow and David Minh
 * and derived from:
 * Chia-En Chang, Michael J. Potter, and Michael K. Gilson. * Calculation of 
 * molecular configuration integrals.  * The Journal of Physical Chemistry B,
 * 107(4):1048–1055, 2003. doi:10.1021/jp027149c.
 * Simon Hikiri, Takashi Yoshidome, and Mitsunori Ikeguchi. Computational methods
 * for configurational entropy using internal and cartesian coordinates. Journal
 * of Chemical Theory and Computation, 12(12):5990–6000, 2016. PMID: 27951672.
 * David D. L. Minh. Alchemical grid dock (algdock): binding free energy calculations
 * between flexible ligands and rigid receptors. Journal of Computational Chemistry,
 * 41(7):715–730, 2020. doi:https://doi.org/10.1002/jcc.26036.
**/

// ecternal atom i alay lighter. if equal, petite id -> big id


class InternalCoordinates {
public:

    /** @name Main functions. This function should be called before accesing 
	 any interface functions **/

    /**@{**/

	/**
	* @brief This function computes the BAT graph in two formats
	* @param bAtomList Atom list read from Amber.
	*/
	void compute(const std::vector<bSpecificAtom>& bAtomList);

	/**
	* @brief This function computes the BAT graph
	* @param bAtomList Atom list read from Amber.
	*/
	void computeBATOld(const std::vector<bSpecificAtom>& bAtomList);
	void computeBAT(const std::vector<bSpecificAtom>& bAtomList);

	/**
	* @brief This function prints the BAT graph
	*/
	void PrintBAT(void){
		std::vector<BOND>::iterator bIt;
		std::vector<ANGLE>::iterator aIt;
		std::vector<TORSION>::iterator tIt;

		std::cout << "BAT\n";
		for(bIt = perMolBonds.back().begin(); bIt < perMolBonds.back().end(); bIt++){
			std::cout << "bond " ;
			bIt->Print();
			std::cout << std::endl;
		}
		for(aIt = perMolAngles.back().begin(); aIt < perMolAngles.back().end(); aIt++){
			std::cout << "angle " ;
			aIt->Print();
			std::cout << std::endl;
		}
		for(tIt = perMolTorsions.back().begin(); tIt < perMolTorsions.back().end(); tIt++){
			std::cout << "torsion " ;
			tIt->Print();
			std::cout << std::endl;
		}
		//std::cout << std::endl;

	}

	void updateVisited(std::vector<bSpecificAtom>& bAtomList);

	/**
	* @brief This function prints the BAT graph
	*/
	void PrintRoot(void){

		std::vector<TORSION>::iterator tIt;

		std::cout << "BAT root\n";
		for(tIt = roots.begin(); tIt < roots.end(); tIt++){
			std::cout << "torsion " ;
			tIt->Print();
			std::cout << std::endl;
		}
		//std::cout << std::endl;

	}

	/**	
	* @brief Puts BAT graph in level / offset format
	* @param bAtomList Atom list read from Amber.
	*/
	void computeLevelsAndOffsets(const std::vector<bSpecificAtom>& bAtomList);
    
	/**@}**/

    /** @name Interface functions
	*/

    /**@{**/


    const TORSION& getLastRoot() const;
	const std::vector<BOND>& getBonds() const;
	const std::vector<BOND>& getRingClosingBonds() const;
	const std::vector<ANGLE>& getAngles() const;
	const std::vector<TORSION>& getTorsions() const;

	/**	
	* @brief Return the molecular graph in a level / offset format which
	* contains a table of : level, offset, atom and parent
	* @return 2D vector of BAT_ATOM (Amber index)
	*/
	const std::vector<std::vector<BAT_ATOM>>& getLevelGraph() const;

	/**	
	* @brief
	* @return
	*/
	const std::vector<std::vector<int>>& getGraphedAtoms() const;

	/**	
	* @brief
	* @return
	*/
	const std::vector<std::vector<BOND>>& getGraphedBonds() const;

	/**@}**/

    /** @name Cartesian to BAT and inverse onversion functions
	*/

    /**@{**/
	
	/**	
	* @brief
	* @param Amber index
	* @return BAT index
	*/
	int amber2BAT(int amberIx) const;

	/**	
	* @brief
	* @param BAT index
	* @return Amber index
	*/
	int BAT2amber(int bat) const;
	/**@}**/

    bool computeRoot(const std::vector<bSpecificAtom>& bAtomList);

private:

	void sortByMass(std::vector<AmberAtom>& v, bool reverse);

	bool isSelected(const AmberAtom& a) const;
	void selectAtom(const AmberAtom& a);
	void selectAtom(int a);

	bool isTerminal(const AmberAtom& a) const;

	std::vector<BOND> allbonds, ringClosingBonds;
	std::vector<ANGLE> allangles;
	std::vector<TORSION> alltorsions;

	std::vector<std::vector<BOND>> perMolBonds; // Laurentiu
	std::vector<std::vector<ANGLE>> perMolAngles; // Laurentiu
	std::vector<std::vector<TORSION>> perMolTorsions; // Laurentiu


	std::vector<std::vector<BAT_ATOM>> levelGraph;
	std::vector<std::vector<int>> graphedAtoms;
	std::vector<std::vector<BOND>> graphedBonds;

	std::vector<AmberAtom> a0_list, a1_list, a2_list, a3_list;

	std::vector<int> selectedAtoms;
	std::vector<int> indexMap;

    TORSION root;
	std::vector<TORSION> roots;
};

#endif //__INTERNALCOORDINATES__