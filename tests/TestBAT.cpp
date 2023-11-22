#include "InternalCoordinates.hpp"
#include "bSpecificAtom.hpp"
#include "readAmberInput.hpp"
#include "Context.hpp"

#include <vector>
#include <string>

class TEST_CONTEXT : public Context {
public:
    std::vector<bSpecificAtom> getAtoms() {
        return atoms;
    }
};

struct LEVEL_OFFSET_ATOM {
    int level = -1,
        offset = -1,
        atom = -1,
        parent = -1;
};

template<typename T>
bool equal(const std::vector<T>& a, const std::vector<T>& b) {
    if (a.size() != b.size()) {
        return false;
    }

    if (!std::equal(a.begin(), a.end(), b.begin())) {
        return false;
    }

    return true;
}

bool test2But() {
    TEST_CONTEXT c;
    c.loadAmberSystem("2but/ligand.prmtop", "2but/ligand.rst7");

    InternalCoordinates ic;
    ic.compute(c.getAtoms());

    std::vector<BOND> theoretical_bonds = {
        { 14, 0 },
        { 0, 1 },
        { 5, 1 },
        { 2, 1 },
        { 6, 2 },
        { 7, 2 },
        { 4, 2 },
        { 11, 4 },
        { 12, 4 },
        { 13, 4 },
        { 3, 1 },
        { 8, 3 },
        { 9, 3 },
        { 10, 3 }
    };

    if (!equal(ic.getBonds(), theoretical_bonds)) {
        std::cout << "Incorrect bonds in 2but." << std::endl;
        return false;
    }

    std::vector<BOND> theoreticalRingClosingBonds;
    if (!equal(ic.getRingClosingBonds(), theoreticalRingClosingBonds)) {
        std::cout << "Incorrect ring closing bonds in 2but." << std::endl;
        return false;
    }

    std::vector<LEVEL_OFFSET_ATOM> levelOffset = {
        { 0, 0, 14, -1 },
        { 1, 0, 0, 14 },
        { 2, 0, 1, 0 },
        { 3, 0, 5, 1 },
        { 3, 1, 2, 1 },
        { 3, 2, 3, 1 },
        { 4, 0, 6, 2 },
        { 4, 1, 7, 2 },
        { 4, 2, 8, 3 },
        { 4, 3, 9, 3 },
        { 4, 4, 10, 3 },
        { 4, 5, 4, 2 },
        { 5, 0, 11, 4 },
        { 5, 1, 12, 4 },
        { 5, 2, 13, 4 },
    };


    std::vector<ANGLE> theoretical_angles = {
        { 5, 1, 0 },
        { 2, 1, 0 },
        { 6, 2, 1 },
        { 7, 2, 1 },
        { 4, 2, 1 },
        { 11, 4, 2 },
        { 12, 4, 2 },
        { 13, 4, 2 },
        { 3, 1, 2 },
        { 8, 3, 1 },
        { 9, 3, 1 },
        { 10, 3, 1 },
    };

    if (!equal(ic.getAngles(), theoretical_angles)) {
        std::cout << "Incorrect angles in 2but." << std::endl;
        return false;
    }

    std::vector<TORSION> theoretical_torsions = {
        { 5, 1, 0, 14 },
        { 2, 1, 0, 14 },
        { 6, 2, 1, 5 },
        { 7, 2, 1, 5 },
        { 4, 2, 1, 5 },
        { 11, 4, 2, 6 },
        { 12, 4, 2, 6 },
        { 13, 4, 2, 6 },
        { 3, 1, 2, 6 },
        { 8, 3, 1, 5 },
        { 9, 3, 1, 5 },
        { 10, 3, 1, 5 }
    };

    if (!equal(ic.getTorsions(), theoretical_torsions)) {
        std::cout << "Incorrect torsions in 2but." << std::endl;
        return false;
    }
    
    return true;
}

bool test3cycles() {
    TEST_CONTEXT c;
    c.loadAmberSystem("3cycles/3cycles.prmtop", "3cycles/3cycles.rst7");

    InternalCoordinates ic;
    ic.compute(c.getAtoms());

    std::vector<BOND> theoretical_bonds = {
        { 57, 55 },
        { 55, 38 },
        { 39, 38 },
        { 40, 38 },
        { 41, 40 },
        { 42, 40 },
        { 43, 40 },
        { 44, 43 },
        { 53, 43 },
        { 45, 44 },
        { 46, 44 },
        { 54, 53 },
        { 51, 53 },
        { 47, 46 },
        { 48, 46 },
        { 52, 51 },
        { 49, 48 },
        { 50, 49 },
        { 56, 55 },
        { 36, 38 },
        { 37, 36 },
        { 34, 36 },
        { 18, 34 },
        { 19, 18 },
        { 20, 18 },
        { 21, 20 },
        { 22, 20 },
        { 23, 20 },
        { 24, 23 },
        { 32, 23 },
        { 25, 24 },
        { 26, 24 },
        { 33, 32 },
        { 30, 32 },
        { 27, 26 },
        { 28, 26 },
        { 31, 30 },
        { 29, 28 },
        { 35, 34 },
        { 16, 18 },
        { 17, 16 },
        { 14, 16 },
        { 12, 14 },
        { 13, 12 },
        { 9, 12 },
        { 10, 9 },
        { 11, 9 },
        { 6, 9 },
        { 7, 6 },
        { 8, 6 },
        { 3, 6 },
        { 4, 3 },
        { 5, 3 },
        { 0, 3 },
        { 1, 0 },
        { 2, 0 },
        { 15, 14 }
    };

    if (!equal(ic.getBonds(), theoretical_bonds)) {
        std::cout << "Incorrect bonds in 3cycles." << std::endl;
        return false;
    }

    std::vector<BOND> theoreticalRingClosingBonds = {
        { 12, 0 },
        { 30, 28 },
        { 51, 48 }
    };

    if (!equal(ic.getRingClosingBonds(), theoreticalRingClosingBonds)) {
        std::cout << "Incorrect ring closing bonds in 3cycles." << std::endl;
        return false;
    }

    std::vector<ANGLE> theoretical_angles = {
        { 39, 38, 55 },
        { 40, 38, 55 },
        { 41, 40, 38 },
        { 42, 40, 38 },
        { 43, 40, 38 },
        { 44, 43, 40 },
        { 53, 43, 40 },
        { 45, 44, 43 },
        { 46, 44, 43 },
        { 54, 53, 43 },
        { 51, 53, 43 },
        { 47, 46, 44 },
        { 48, 46, 44 },
        { 52, 51, 48 },
        { 49, 48, 46 },
        { 50, 49, 48 },
        { 56, 55, 38 },
        { 36, 38, 40 },
        { 37, 36, 38 },
        { 34, 36, 38 },
        { 18, 34, 36 },
        { 19, 18, 34 },
        { 20, 18, 34 },
        { 21, 20, 18 },
        { 22, 20, 18 },
        { 23, 20, 18 },
        { 24, 23, 20 },
        { 32, 23, 20 },
        { 25, 24, 23 },
        { 26, 24, 23 },
        { 33, 32, 23 },
        { 30, 32, 23 },
        { 27, 26, 24 },
        { 28, 26, 24 },
        { 31, 30, 28 },
        { 29, 28, 26 },
        { 35, 34, 18 },
        { 16, 18, 20 },
        { 17, 16, 18 },
        { 14, 16, 18 },
        { 12, 14, 16 },
        { 13, 12, 14 },
        { 9, 12, 14 },
        { 10, 9, 12 },
        { 11, 9, 12 },
        { 6, 9, 12 },
        { 7, 6, 9 },
        { 8, 6, 9 },
        { 3, 6, 9 },
        { 4, 3, 6 },
        { 5, 3, 6 },
        { 0, 3, 6 },
        { 1, 0, 3 },
        { 2, 0, 3 },
        { 15, 14, 12 }
    };

    if (!equal(ic.getAngles(), theoretical_angles)) {
        std::cout << "Incorrect angles in 3cycles." << std::endl;
        return false;
    }

    std::vector<TORSION> theoretical_torsions = {
        { 39, 38, 55, 57 },
        { 40, 38, 55, 57 },
        { 41, 40, 38, 39 },
        { 42, 40, 38, 39 },
        { 43, 40, 38, 39 },
        { 44, 43, 40, 41 },
        { 53, 43, 40, 41 },
        { 45, 44, 43, 40 },
        { 46, 44, 43, 40 },
        { 54, 53, 43, 40 },
        { 51, 53, 43, 40 },
        { 47, 46, 44, 45 },
        { 48, 46, 44, 45 },
        { 52, 51, 48, 46 },
        { 49, 48, 46, 47 },
        { 50, 49, 48, 46 },
        { 56, 55, 38, 39 },
        { 36, 38, 40, 41 },
        { 37, 36, 38, 39 },
        { 34, 36, 38, 39 },
        { 18, 34, 36, 37 },
        { 19, 18, 34, 36 },
        { 20, 18, 34, 36 },
        { 21, 20, 18, 19 },
        { 22, 20, 18, 19 },
        { 23, 20, 18, 19 },
        { 24, 23, 20, 21 },
        { 32, 23, 20, 21 },
        { 25, 24, 23, 20 },
        { 26, 24, 23, 20 },
        { 33, 32, 23, 20 },
        { 30, 32, 23, 20 },
        { 27, 26, 24, 25 },
        { 28, 26, 24, 25 },
        { 31, 30, 28, 26 },
        { 29, 28, 26, 27 },
        { 35, 34, 18, 19 },
        { 16, 18, 20, 21 },
        { 17, 16, 18, 19 },
        { 14, 16, 18, 19 },
        { 12, 14, 16, 17 },
        { 13, 12, 14, 16 },
        { 9, 12, 14, 16 },
        { 10, 9, 12, 13 },
        { 11, 9, 12, 13 },
        { 6, 9, 12, 13 },
        { 7, 6, 9, 10 },
        { 8, 6, 9, 10 },
        { 3, 6, 9, 10 },
        { 4, 3, 6, 7 },
        { 5, 3, 6, 7 },
        { 0, 3, 6, 7 },
        { 1, 0, 3, 4 },
        { 2, 0, 3, 4 },
        { 15, 14, 12, 13 }
    };

    if (!equal(ic.getTorsions(), theoretical_torsions)) {
        std::cout << "Incorrect torsions in 3cycles." << std::endl;
        return false;
    }
    
    return true;
}

int main() {
    if (!test2But()) {
        return 1;
    }

    // if (!test3cycles()) {
    //     return 1;
    // }

    // endConstruction()

    return 0;
}