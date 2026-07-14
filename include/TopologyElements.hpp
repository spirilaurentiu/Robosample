#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "RobotModel.hpp"

/**
 * @brief Nonbonded interaction treatment requested for the OpenMM system,
 *        mirroring OpenMM's NonbondedForce method enum.
 * @note The periodic methods (CutoffPeriodic, Ewald, PME) require
 *       SystemTopology::boxVectors to be set; NoCutoff forbids periodic boundaries.
 */
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

/**
 * @brief Structure-of-arrays payload describing one molecular system at the
 *        Python<->engine boundary: atoms, bonded terms, nonbonded parameters,
 *        z-matrix, periodic box, and virtual sites.
 *
 * @par Producer / consumer and lifetime
 * Owned by @c Context as a public member for the whole run. Filled in place from
 * Python (the producer) after parsing the prmtop/inpcrd, then read by
 * World::buildModel / ModelBuilder and the OpenMM system builder (the consumers).
 * It is an input snapshot: the engine treats it as read-only after the build.
 *
 * @par Index and array conventions
 * Atom arrays are indexed by the global BFS atom index (the array position is the
 * index, matching the OpenMM particle order). All parallel arrays in one section
 * share that section's index (e.g. every @c atoms* array is indexed by atom;
 * every @c bonds* array by bond). The @c *Begin / @c *End vectors are CSR-style
 * half-open ranges [begin,end) per molecule (length @c numMolecules), slicing the
 * per-interaction arrays; @c atomsBegin[mol] doubles as molecule @p mol's root atom.
 *
 * @par Units (INV-3, consistent MD system)
 * Lengths nm, angles radians, energies kJ/mol, charge in elementary units,
 * mass in daltons, temperature in kelvin. Per-field units are annotated inline.
 *
 * @note The sub-struct split of this god-struct is deferred behind a decision
 *       record; it is documented here as it stands.
 */
struct SystemTopology {
    // -------------------------------------------------------------------------
    // Molecule ranges  -- [begin, end) into each array above; length == numMolecules
    // -------------------------------------------------------------------------

    /** @brief Number of molecules; the length of every @c *Begin / @c *End range
     *  array and of @c rootMobilities. */
    int numMolecules{0};
    std::vector<int> atomsBegin; ///< Per-molecule first atom index; also the molecule's root atom.
    std::vector<int> atomsEnd;   ///< Per-molecule one-past-last atom index.
    std::vector<int> bondsBegin;             ///< Per-molecule [begin,end) into the @c bonds* arrays.
    std::vector<int> bondsEnd;
    std::vector<int> anglesBegin;            ///< Per-molecule [begin,end) into the @c angles* arrays.
    std::vector<int> anglesEnd;
    std::vector<int> periodicTorsionsBegin;  ///< Per-molecule [begin,end) into the @c periodicTorsions* arrays.
    std::vector<int> periodicTorsionsEnd;
    std::vector<int> harmonicTorsionsBegin;  ///< Per-molecule [begin,end) into the @c harmonicTorsions* arrays.
    std::vector<int> harmonicTorsionsEnd;
    std::vector<int> zMatrixBegin;           ///< Per-molecule [begin,end) into the @c zMatrix* arrays.
    std::vector<int> zMatrixEnd;
    std::vector<int> ureyBradleyBegin;       ///< Per-molecule [begin,end) into the @c ureyBradley* arrays.
    std::vector<int> ureyBradleyEnd;
    std::vector<int> scaling14Begin;         ///< Per-molecule [begin,end) into the @c scaling14* arrays.
    std::vector<int> scaling14End;
    std::vector<int> exclusionBegin;         ///< Per-molecule [begin,end) into the @c exclusion* arrays.
    std::vector<int> exclusionEnd;
    std::vector<int> atomsRootIndex;         ///< Per-molecule root atom index (BFS root of each molecule's tree).
    std::vector<JointType> rootMobilities;   ///< Per-molecule joint attaching that molecule's root to Ground.

    // -------------------------------------------------------------------------
    // Atom arrays  (BFS order, unit-converted)
    // -------------------------------------------------------------------------

    int numAtoms{0};                              ///< Number of atoms; length of every @c atoms* array.
    std::vector<std::string> atomsUniqueName;     ///< Per-atom unique name (analysis/IO label).
    std::vector<int> atomsNonbondedIndex;         ///< Per-atom index into the nonbonded parameter tables.
    std::vector<int> atomsPrmtopIndex;            ///< Per-atom original prmtop index (pre-BFS-reorder).

    std::vector<std::string> atomsElementName;    ///< Per-atom element name.
    std::vector<std::string> atomsElementSymbol;  ///< Per-atom element symbol.
    std::vector<double> atomsMass;    ///< Mass.                               [daltons]
    std::vector<double> atomsCharge;  ///< Partial charge.                     [e]
    std::vector<double> atomsSigma;   ///< Lennard-Jones sigma (vdW radius).   [nm]
    std::vector<double> atomsEpsilon; ///< Lennard-Jones epsilon (well depth). [kJ/mol]
    std::vector<double> atomsRadius;  ///< GBSA solvent radius.                [nm]
    std::vector<double> atomsScreen;  ///< OBC screening factor.               [dimensionless]
    std::vector<double> atomsX;       ///< x coordinate (reference structure). [nm]
    std::vector<double> atomsY;       ///< y coordinate (reference structure). [nm]
    std::vector<double> atomsZ;       ///< z coordinate (reference structure). [nm]
    std::vector<int> atomsAtomicNumber;    ///< Per-atom atomic number Z.
    std::vector<int> atomsNumBondsInvolved; ///< Per-atom count of bonds incident on the atom.

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

    int cmapGridSize{0};                     ///< Side length of each square CMAP energy grid.
    std::vector<double> cmapGridEnergy;      ///< Flattened CMAP correction grids (kJ/mol).

    std::vector<int> cmapTorsionMapIndex;    ///< Per-CMAP-term index of the grid it uses.
    std::vector<int> cmapTorsionA1;          ///< First torsion, atom 1 (global index).
    std::vector<int> cmapTorsionA2;          ///< First torsion, atom 2.
    std::vector<int> cmapTorsionA3;          ///< First torsion, atom 3.
    std::vector<int> cmapTorsionA4;          ///< First torsion, atom 4.
    std::vector<int> cmapTorsionB1;          ///< Second torsion, atom 1.
    std::vector<int> cmapTorsionB2;          ///< Second torsion, atom 2.
    std::vector<int> cmapTorsionB3;          ///< Second torsion, atom 3.
    std::vector<int> cmapTorsionB4;          ///< Second torsion, atom 4.

    bool hasNBfix = false;                   ///< Whether explicit NBFIX off-diagonal LJ tables are present.
    int numNBTypes = 0;                      ///< Number of Lennard-Jones atom types (side of the aCoef/bCoef tables).
    std::vector<double> aCoef;               ///< NBFIX A-coefficient table, numNBTypes^2 (LJ r^-12 term).
    std::vector<double> bCoef;               ///< NBFIX B-coefficient table, numNBTypes^2 (LJ r^-6 term).

    bool useGBSAOBC2 = false;                ///< Enable GBSA-OBC2 implicit solvent.
    double gbsaSolventDielectric = 78.5;     ///< GBSA solvent dielectric constant.
    double gbsaSoluteDielectric = 1.0;       ///< GBSA solute dielectric constant.

    NonbondedMethod nonbondedMethod = NonbondedMethod::NoCutoff; ///< Nonbonded treatment for the OpenMM system.
    double nonbondedCutoff = 1.2;            ///< Nonbonded cutoff distance (nm); used by the cutoff/periodic methods.

    // -------------------------------------------------------------------------
    // Periodic box (explicit solvent)
    // -------------------------------------------------------------------------
    // Three REDUCED lattice vectors in OpenMM's lower-triangular convention
    //   a = (ax, 0,  0 ),  b = (bx, by, 0 ),  c = (cx, cy, cz),
    // stored row-major as 9 doubles [a.x a.y a.z  b.x b.y b.z  c.x c.y c.z], in
    // nm. Filled directly from ParmEd's parm.box_vectors (already reduced), so no
    // a/b/c/alpha/beta/gamma reduction happens on the C++ side. EMPTY for any
    // non-periodic method; REQUIRED (length 9) for CutoffPeriodic / Ewald / PME.
    // The box is set on the OpenMM System (setDefaultPeriodicBoxVectors) BEFORE
    // the Context is created, which is mandatory for PME.
    std::vector<double> boxVectors;

    // Reciprocal-space accuracy for Ewald/PME (OpenMM setEwaldErrorTolerance).
    // Ignored by the non-Ewald methods. 5e-4 is OpenMM's usual default.
    double ewaldErrorTolerance = 5.0e-4;

    // -------------------------------------------------------------------------
    // Virtual sites (extra points: massless particles placed by real atoms)
    // -------------------------------------------------------------------------
    // 4-point water (OPC, TIP4P family) carries a massless EP that holds the
    // negative charge; its position is the affine combination
    //   r_site = w1*r_a1 + w2*r_a2 + w3*r_a3,   w1 + w2 + w3 = 1
    // (a 3-particle AVERAGE site). These MUST be declared to OpenMM via
    // setVirtualSite(ThreeParticleAverageSite) so the integrator skips them and
    // redistributes their force onto the parents; otherwise OpenMM freezes the
    // massless particle in place and the EP detaches as the molecule moves.
    // Indices are global/BFS atom indices (the OpenMM particle order). Only the
    // 3-particle average type is represented here (covers OPC/TIP4P/-Ew/-2005);
    // out-of-plane sites (e.g. TIP5P) would need an additional type.
    int numVirtualSites{0};        ///< Number of 3-particle-average virtual sites; length of every @c vs* array.
    std::vector<int> vsSite;       ///< global index of the massless EP particle.
    std::vector<int> vsAtom1;      ///< parent atom 1 (global index).
    std::vector<int> vsAtom2;      ///< parent atom 2 (global index).
    std::vector<int> vsAtom3;      ///< parent atom 3 (global index).
    std::vector<double> vsWeight1; ///< weight on parent 1.
    std::vector<double> vsWeight2; ///< weight on parent 2.
    std::vector<double> vsWeight3; ///< weight on parent 3.

    double thermostatTemperature = 300.0; ///< Target temperature (K) for the OpenMM thermostat.
    double collisionFrequency = 1.0;      ///< Thermostat collision/friction frequency (1/ps).
    int seed = 0;                         ///< RNG seed handed to the OpenMM integrator/thermostat.
};