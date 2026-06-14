#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "RobotModel.hpp"

enum NonbondedMethod : std::uint8_t {
    // No cutoff is applied to nonbonded interactions. The full set of N^2 interactions is computed exactly.
    // This necessarily means that periodic boundary conditions cannot be used. This is the default.
    NoCutoff = 0,

    // Interactions beyond the cutoff distance are ignored. Coulomb interactions closer than the cutoff
    // distance are modified using the reaction field method.
    CutoffNonPeriodic,

    // Periodic boundary conditions are used, so that each particle interacts only with the nearest periodic
    // copy of each other particle. Interactions beyond the cutoff distance are ignored. Coulomb interactions
    // closer than the cutoff distance are modified using the reaction field method.
    CutoffPeriodic,

    // Periodic boundary conditions are used, and Ewald summation is used to compute the interaction of each
    // particle with all periodic copies of every other particle.
    Ewald,

    // Periodic boundary conditions are used, and Particle-Mesh Ewald (PME) summation is used to compute the
    // interaction of each particle with all periodic copies of every other particle.
    PME,
};

struct SystemTopology {
    // -------------------------------------------------------------------------
    // Molecule ranges  -- [begin, end) into each array above; length == numMolecules
    // -------------------------------------------------------------------------

    int numMolecules{0};
    std::vector<int> atomsBegin;
    std::vector<int> atomsEnd;
    std::vector<int> bondsBegin;
    std::vector<int> bondsEnd;
    std::vector<int> anglesBegin;
    std::vector<int> anglesEnd;
    std::vector<int> periodicTorsionsBegin;
    std::vector<int> periodicTorsionsEnd;
    std::vector<int> harmonicTorsionsBegin;
    std::vector<int> harmonicTorsionsEnd;
    std::vector<int> zMatrixBegin;
    std::vector<int> zMatrixEnd;
    std::vector<int> ureyBradleyBegin;
    std::vector<int> ureyBradleyEnd;
    std::vector<int> scaling14Begin;
    std::vector<int> scaling14End;
    std::vector<int> exclusionBegin;
    std::vector<int> exclusionEnd;
    std::vector<int> atomsRootIndex;
    std::vector<RootMobility> rootMobilities;

    // -------------------------------------------------------------------------
    // Atom arrays  (BFS order, unit-converted)
    // -------------------------------------------------------------------------

    int numAtoms{0};
    std::vector<std::string> atomsUniqueName;
    std::vector<int> atomsNonbondedIndex;
    std::vector<int> atomsPrmtopIndex;

    std::vector<std::string> atomsElementName;
    std::vector<std::string> atomsElementSymbol;
    std::vector<double> atomsMass;    ///< Mass.                               [daltons]
    std::vector<double> atomsCharge;  ///< Partial charge.                     [e]
    std::vector<double> atomsSigma;   ///< Lennard-Jones sigma (vdW radius).   [nm]
    std::vector<double> atomsEpsilon; ///< Lennard-Jones epsilon (well depth). [kJ/mol]
    std::vector<double> atomsRadius;  ///< GBSA solvent radius.                [nm]
    std::vector<double> atomsScreen;  ///< OBC screening factor.               [dimensionless]
    std::vector<double> atomsX;       ///< x coordinate (reference structure). [nm]
    std::vector<double> atomsY;       ///< y coordinate (reference structure). [nm]
    std::vector<double> atomsZ;       ///< z coordinate (reference structure). [nm]
    std::vector<int> atomsAtomicNumber;
    std::vector<int> atomsNumBondsInvolved;

    // -------------------------------------------------------------------------
    // Bond arrays  (BFS order)
    // -------------------------------------------------------------------------

    int numBonds{0};                      ///< Total bonds (tree + ring-closing).
    std::vector<int> bondsI;              ///< BFS index of bond endpoint atom 1.
    std::vector<int> bondsJ;              ///< BFS index of bond endpoint atom 2.
    std::vector<int> bondsMoleculeIndex;  ///< Molecule index this bond belongs to.
    std::vector<bool> bondsRingClosing;   ///< Whether this bond is a ring-closing bond.
    std::vector<double> bondsStiffness;   ///< Harmonic force constant.         [kJ/mol/nm^2]
    std::vector<double> bondsEquilibrium; ///< Equilibrium bond length.         [nm]

    // -------------------------------------------------------------------------
    // Angle arrays
    // -------------------------------------------------------------------------

    int numAngles{0};
    std::vector<int> anglesI;
    std::vector<int> anglesJ;
    std::vector<int> anglesK;
    std::vector<double> anglesEquilibrium; ///< Equilibrium angle.               [rad]
    std::vector<double> anglesStiffness;   ///< Harmonic force constant.         [kJ/mol/rad^2]

    // -------------------------------------------------------------------------
    // Periodic torsion arrays
    // -------------------------------------------------------------------------

    int numPeriodicTorsions{0};
    std::vector<bool> periodicTorsionsImproper;
    std::vector<int> periodicTorsionsI;
    std::vector<int> periodicTorsionsJ;
    std::vector<int> periodicTorsionsK;
    std::vector<int> periodicTorsionsL;
    std::vector<int> periodicTorsionsN;            ///< Periodicity.
    std::vector<double> periodicTorsionsPhase;     ///< Phase offset.             [rad]
    std::vector<double> periodicTorsionsStiffness; ///< Force constant.           [kJ/mol]

    // -------------------------------------------------------------------------
    // Harmonic torsion arrays
    // -------------------------------------------------------------------------

    int numHarmonicTorsions{0};
    std::vector<int> harmonicTorsionsI;
    std::vector<int> harmonicTorsionsJ;
    std::vector<int> harmonicTorsionsK;
    std::vector<int> harmonicTorsionsL;
    std::vector<double> harmonicTorsionsStiffness; ///< Force constant.           [kJ/mol]
    std::vector<double> harmonicTorsionsPhase;     ///< Equilibrium angle.        [rad]

    // -------------------------------------------------------------------------
    // Z-matrix arrays  (tree-traversal order, rooted at atomsRootIndex)
    // Length n; sentinel value -1 fills inapplicable leading rows (see below).
    // -------------------------------------------------------------------------

    int numZMatrixRows{0};
    std::vector<int> zMatrixI; ///< Global atom index of the atom placed at row r.
    std::vector<int> zMatrixJ; ///< Bond-length reference atom.  Row 0   : -1 (root).
    std::vector<int> zMatrixK; ///< Bond-angle reference atom.   Rows 0-1: -1.
    std::vector<int> zMatrixL; ///< Dihedral reference atom.     Rows 0-2: -1.

    // -------------------------------------------------------------------------
    // Urey-Bradley 1-3 interactions
    // -------------------------------------------------------------------------

    int numUreyBradley{0};
    std::vector<int> ureyBradleyI;              ///< Global index of atom 1 (outer atom of angle i-j-k).
    std::vector<int> ureyBradleyK;              ///< Global index of atom 3 (outer atom of angle i-j-k).
    std::vector<double> ureyBradleyStiffness;   ///< Harmonic force constant. [kJ/mol/nm^2]
    std::vector<double> ureyBradleyEquilibrium; ///< Nominal 1-3 distance.                           [nm]

    // -------------------------------------------------------------------------
    // 1-4 pair scaling
    // -------------------------------------------------------------------------

    int numScaling14{0};
    std::vector<int> scaling14I;                ///< Global index of atom 1 (first atom of dihedral i-j-k-l).
    std::vector<int> scaling14L;                ///< Global index of atom 4 (last  atom of dihedral i-j-k-l).
    std::vector<double> scaling14ChargeProduct; ///< q1*q4 pre-scaled by the 1-4 electrostatic factor.  [e^2]
    std::vector<double> scaling14Epsilon;       ///< Combined LJ well depth.                          [kJ/mol]
    std::vector<double> scaling14Sigma;         ///< Combined LJ radius.                              [nm]

    // -------------------------------------------------------------------------
    // Exclusions
    // -------------------------------------------------------------------------

    int numExclusions{0};
    std::vector<int> exclusionI; ///< Global index of atom 1.
    std::vector<int> exclusionJ; ///< Global index of atom 2.

    // -------------------------------------------------------------------------
    // CMAP torsion arrays
    // -------------------------------------------------------------------------

    int cmapGridSize{0};
    std::vector<double> cmapGridEnergy;

    std::vector<int> cmapTorsionMapIndex;
    std::vector<int> cmapTorsionA1;
    std::vector<int> cmapTorsionA2;
    std::vector<int> cmapTorsionA3;
    std::vector<int> cmapTorsionA4;
    std::vector<int> cmapTorsionB1;
    std::vector<int> cmapTorsionB2;
    std::vector<int> cmapTorsionB3;
    std::vector<int> cmapTorsionB4;

    bool hasNBfix = false;
    int numNBTypes = 0;
    std::vector<double> aCoef;
    std::vector<double> bCoef;

    bool useGBSAOBC2 = false;
    double gbsaSolventDielectric = 78.5;
    double gbsaSoluteDielectric = 1.0;

    NonbondedMethod nonbondedMethod = NonbondedMethod::NoCutoff;
    double nonbondedCutoff = 1.2;

    double thermostatTemperature = 300.0;
    double collisionFrequency = 1.0;
    int seed = 0;
};