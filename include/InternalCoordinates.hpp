#ifndef __INTERNALCOORDINATES__
#define __INTERNALCOORDINATES__

/* -------------------------------------------------------------------------- *
 *                                gMolModel                                   *
 * -------------------------------------------------------------------------- *
 *  This is an extension of SimTK Molmodel package.                           *
 * -------------------------------------------------------------------------- */
/** @file This defines the Internal Coordinates class
 **/

#include "Topology.hpp"

using namespace SimTK;

struct BOND {
	int first = -1,
		second = -1;
};

struct ANGLE {
	int first = -1,
		second = -1,
		third = -1;
};

struct TORSION {
	int first = -1,
		second = -1,
		third = -1,
		fourth = -1;
};

struct AmberAtom {
	SimTK::Real mass = -1.0;
	int bonds = -1,
		amberId = -1;

	AmberAtom(const bSpecificAtom& b);
	AmberAtom(const bSpecificAtom* b);

	bool operator==(const AmberAtom& rhs);
	bool operator!=(const AmberAtom& rhs);
};

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

class InternalCoordinates {
public:
	void compute(const std::vector<bSpecificAtom>& bAtomList);

    const TORSION& getRoot() const;
	const std::vector<BOND>& getBonds() const;
	const std::vector<ANGLE>& getAngles() const;
	const std::vector<TORSION>& getTorsions() const;

	int amber2BAT(int amberIx) const;
	int BAT2amber(int bat) const;

private:
    void computeRoot(const std::vector<bSpecificAtom>& bAtomList);

	void sortByMass(std::vector<AmberAtom>& v, bool reverse);

	bool isSelected(const AmberAtom& a) const;
	void selectAtom(const AmberAtom& a);
	void selectAtom(int a);

	bool isTerminal(const AmberAtom& a) const;

	std::vector<BOND> bonds;
	std::vector<ANGLE> angles;
	std::vector<TORSION> torsions;

	std::vector<AmberAtom> a0_list, a1_list, a2_list, a3_list;

	std::vector<int> selectedAtoms;
	std::vector<int> indexMap;

    TORSION root;
};

#endif //__INTERNALCOORDINATES__