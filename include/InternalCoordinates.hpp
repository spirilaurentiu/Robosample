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