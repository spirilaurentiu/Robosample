#include <Python.h>

#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>

#include "molmodel/internal/Compound.h"

#include "CompoundSystem.h"
#include "Context.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc() = "Robosample bindings";

    py::enum_<NonbondedMethod>(m, "NonbondedMethod", R"doc(
Methods for computing nonbonded interactions.
)doc")
        .value("NO_CUTOFF", NonbondedMethod::NoCutoff, R"doc(
No cutoff is applied to nonbonded interactions.

The full set of N^2 interactions is computed exactly.
Periodic boundary conditions cannot be used.
This is the default.
)doc")
        .value("CUTOFF_NON_PERIODIC", NonbondedMethod::CutoffNonPeriodic, R"doc(
Interactions beyond the cutoff distance are ignored.

Coulomb interactions closer than the cutoff distance
are modified using the reaction field method.
)doc")
        .value("CUTOFF_PERIODIC", NonbondedMethod::CutoffPeriodic, R"doc(
Periodic boundary conditions are used.

Each particle interacts only with the nearest periodic
copy of each other particle.

Interactions beyond the cutoff distance are ignored.

Coulomb interactions closer than the cutoff distance
are modified using the reaction field method.
)doc")
        .value("EWALD", NonbondedMethod::Ewald, R"doc(
Periodic boundary conditions are used.
                    dst.extend(x if x == ZM_SENTINEL else x + atom_off for x in src)

Ewald summation is used to compute the interaction
of each particle with all periodic copies of every
other particle.
)doc")
        .value("PME", NonbondedMethod::PME, R"doc(
Periodic boundary conditions are used.

Particle-Mesh Ewald (PME) summation is used to compute
the interaction of each particle with all periodic
copies of every other particle.
)doc")
        .export_values();

    py::enum_<SimTK::RootMobility>(m, "RootMobility")
        .value("FREE", SimTK::RootMobility::Free)
        .value("CARTESIAN", SimTK::RootMobility::Cartesian)
        .value("WELD", SimTK::RootMobility::Weld)
        .value("FREE_LINE", SimTK::RootMobility::FreeLine)
        .value("BALL", SimTK::RootMobility::Ball)
        .value("PIN", SimTK::RootMobility::Pin);

    py::enum_<SimTK::BondMobility::Mobility>(m, "BondMobility")
        .value("FREE", SimTK::BondMobility::Mobility::Free)
        .value("TORSION", SimTK::BondMobility::Mobility::Torsion)
        .value("RIGID", SimTK::BondMobility::Mobility::Rigid)
        .value("BALL_F", SimTK::BondMobility::Mobility::BallF)
        .value("BALL_M", SimTK::BondMobility::Mobility::BallM)
        .value("CYLINDER", SimTK::BondMobility::Mobility::Cylinder)
        .value("TRANSLATION", SimTK::BondMobility::Mobility::Translation)
        .value("FREE_LINE", SimTK::BondMobility::Mobility::FreeLine)
        .value("LINE_ORIENTATION_F", SimTK::BondMobility::Mobility::LineOrientationF)
        .value("LINE_ORIENTATION_M", SimTK::BondMobility::Mobility::LineOrientationM)
        .value("UNIVERSAL_M", SimTK::BondMobility::Mobility::UniversalM)
        .value("SPHERICAL", SimTK::BondMobility::Mobility::Spherical)
        .value("ANGLE_PIN", SimTK::BondMobility::Mobility::AnglePin)
        .value("BEND_STRETCH", SimTK::BondMobility::Mobility::BendStretch)
        .value("SLIDER", SimTK::BondMobility::Mobility::Slider)
        .value("ORTHO_SPHERICAL", SimTK::BondMobility::Mobility::OrthoSpherical);

    py::enum_<RunType>(m, "RunType")
        .value("DEFAULT", RunType::Default, "Standard simulation run.")
        .value("REMC", RunType::REMC, "Replica Exchange Monte Carlo.")
        .value("RENEMC", RunType::RENEMC, "Replica Exchange Non-Equilibrium Monte Carlo.")
        .value("RENE", RunType::RENE, "Replica Exchange Non-Equilibrium.")
        .value("REBASONTOP", RunType::REBASONTOP, "Replica exchange on top of a base simulation.");

    py::enum_<SamplerName>(m, "SamplerName")
        .value("EMPTY", SamplerName::Empty)
        .value("HMC", SamplerName::HMC);

    py::enum_<AcceptRejectMode>(m, "AcceptRejectMode")
        .value("ALWAYS_ACCEPT", AcceptRejectMode::AlwaysAccept)
        .value("METROPOLIS_HASTINGS", AcceptRejectMode::MetropolisHastings);

    py::enum_<IntegratorType>(m, "IntegratorType")
        .value("EMPTY", IntegratorType::Empty)
        .value("VERLET", IntegratorType::Verlet)
        .value("EULER", IntegratorType::Euler)
        .value("EULER2", IntegratorType::Euler2)
        .value("CPODES", IntegratorType::CPodes)
        .value("RUNGE_KUTTA", IntegratorType::RungeKutta)
        .value("RUNGE_KUTTA2", IntegratorType::RungeKutta2)
        .value("RUNGE_KUTTA3", IntegratorType::RungeKutta3)
        .value("RUNGE_KUTTA_FELDBERG", IntegratorType::RungeKuttaFeldberg)
        .value("BEND_STRETCH", IntegratorType::BendStretch)
        .value("OMMVV", IntegratorType::OpenMMVelocityVerlet)
        .value("BOUND_WALK", IntegratorType::BoundWalk)
        .value("BOUND_HMC", IntegratorType::BoundHMC)
        .value("STATIONS_TASK", IntegratorType::StationsTask)
        .value("NOF_INTEGRATORS", IntegratorType::NofIntegrators);

    py::enum_<ThermostatName>(m, "ThermostatName")
        .value("NONE", ThermostatName::None)
        .value("ANDERSEN", ThermostatName::Andersen)
        .value("BERENDSEN", ThermostatName::Berendsen)
        .value("LANGEVIN", ThermostatName::Langevin)
        .value("NOSE_HOOVER", ThermostatName::NoseHoover);

    py::class_<BondFlexibility>(m, "BondFlexibility")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &BondFlexibility::globalIndex1)
        .def_readwrite("globalIndex2", &BondFlexibility::globalIndex2)
        .def_readwrite("uniqueAtomName1", &BondFlexibility::uniqueAtomName1)
        .def_readwrite("uniqueAtomName2", &BondFlexibility::uniqueAtomName2)
        .def_readwrite("mobility", &BondFlexibility::mobility);

    py::class_<Context>(m, "Context")
        .def(py::init<const std::string&, std::int32_t>(), py::arg("base_name"), py::arg("seed"))
        .def("addReplica", &Context::addReplica, "Add an empty replica to the context.")
        .def("addThermodynamicState",
             &Context::addThermodynamicState,
             "Add an empty themodynamic state to the context.")
        .def("run_rex",
             &Context::RunREX,
             py::arg("run_type"),
             py::arg("num_equilibration_rounds"),
             py::arg("num_production_rounds"),
             py::arg("write_frequency"),
             py::arg("write_to_stdio"),
             "Run replica exchange.")
        .def("set_world_temperatures",
             &Context::setWorldTemperatures,
             py::arg("world_temperatures_in_k"),
             "Set the temperatures of the worlds in the context.")
        .def("setVerbose", &Context::setVerbose, "Control if you want extraneous output to cout.")
        .def("setPdbRestartFreq", &Context::setPdbRestartFreq, "Set the PDB restart frequency.")
        .def("setNonbonded", &Context::setNonbonded, "Set nonbonded method and cutoff.")
        .def("setGBSAOptions", &Context::setGBSAOptions, "Set GBSA-OBC2 options.")
        .def("loadAmberSystem", &Context::loadAmberSystem, "Load an AMBER system.")
        .def("initialize_openmm", &Context::initializeOpenMM, "Load an OpenMM system from components.")
        .def("calculate_openmm_energy",
             &Context::calculatePotentialEnergy,
             py::arg("worldIndex"),
             "Calculate the OpenMM energy of the current state for a specific world index.")
        .def("add_world",
             &Context::addWorld,
             py::arg("fixman_torque"),
             py::arg("samples_per_round"),
             py::arg("roll_flexibilities"),
             py::arg("want_spatial_force_history"),
             R"doc(
            Add an empty world.

            Args:
                roll_flexibilities: A list of lists of BondFlexibility objects, 
                                    e.g., [[rb.BondFlexibility(), ...], [...]]
                want_spatial_force_history: A boolean indicating whether to track spatial force history.
                )doc")
        .def("getWorld",
             py::overload_cast<std::size_t>(&Context::getWorld, py::const_),
             py::return_value_policy::reference)
        .def("getWorlds",
             py::overload_cast<>(&Context::getWorlds, py::const_),
             py::return_value_policy::reference);

    py::class_<World>(m, "World")
        .def("add_sampler",
             &World::addSampler,
             py::arg("sampler_name"),
             py::arg("integrator_type"),
             py::arg("thermostat_name"),
             py::arg("use_fixman_potential"),
             py::arg("use_nuts"),
             "Add a sampler to the world.");


    py::class_<SystemTopology>(m, "SystemTopology")
        .def(py::init<>())

        // -------------------------------------------------------------------------
        // Molecule ranges
        // -------------------------------------------------------------------------

        .def_readwrite("num_molecules", &SystemTopology::numMolecules, "Total number of molecules.")
        .def_readwrite("atoms_begin",
                       &SystemTopology::atomsBegin,
                       "Begin index into atoms arrays for each molecule. Length equals num_molecules.")
        .def_readwrite("atoms_end",
                       &SystemTopology::atomsEnd,
                       "End index (exclusive) into atoms arrays for each molecule.")
        .def_readwrite("bonds_begin",
                       &SystemTopology::bondsBegin,
                       "Begin index into bonds arrays for each molecule.")
        .def_readwrite("bonds_end",
                       &SystemTopology::bondsEnd,
                       "End index (exclusive) into bonds arrays for each molecule.")
        .def_readwrite("angles_begin",
                       &SystemTopology::anglesBegin,
                       "Begin index into angles arrays for each molecule.")
        .def_readwrite("angles_end",
                       &SystemTopology::anglesEnd,
                       "End index (exclusive) into angles arrays for each molecule.")
        .def_readwrite("periodic_torsions_begin",
                       &SystemTopology::periodicTorsionsBegin,
                       "Begin index into periodic torsions arrays for each molecule.")
        .def_readwrite("periodic_torsions_end",
                       &SystemTopology::periodicTorsionsEnd,
                       "End index (exclusive) into periodic torsions arrays for each molecule.")
        .def_readwrite("harmonic_torsions_begin",
                       &SystemTopology::harmonicTorsionsBegin,
                       "Begin index into harmonic torsions arrays for each molecule.")
        .def_readwrite("harmonic_torsions_end",
                       &SystemTopology::harmonicTorsionsEnd,
                       "End index (exclusive) into harmonic torsions arrays for each molecule.")
        .def_readwrite("z_matrix_begin",
                       &SystemTopology::zMatrixBegin,
                       "Begin index into Z-matrix arrays for each molecule.")
        .def_readwrite("z_matrix_end",
                       &SystemTopology::zMatrixEnd,
                       "End index (exclusive) into Z-matrix arrays for each molecule.")
        .def_readwrite("urey_bradley_begin",
                       &SystemTopology::ureyBradleyBegin,
                       "Begin index into Urey-Bradley arrays for each molecule.")
        .def_readwrite("urey_bradley_end",
                       &SystemTopology::ureyBradleyEnd,
                       "End index (exclusive) into Urey-Bradley arrays for each molecule.")
        .def_readwrite("scaling14_begin",
                       &SystemTopology::scaling14Begin,
                       "Begin index into 1-4 scaling arrays for each molecule.")
        .def_readwrite("scaling14_end",
                       &SystemTopology::scaling14End,
                       "End index (exclusive) into 1-4 scaling arrays for each molecule.")
        .def_readwrite("exclusion_begin",
                       &SystemTopology::exclusionBegin,
                       "Begin index into exclusions arrays for each molecule.")
        .def_readwrite("exclusion_end",
                       &SystemTopology::exclusionEnd,
                       "End index (exclusive) into exclusions arrays for each molecule.")
        .def_readwrite("atoms_root_index",
                       &SystemTopology::atomsRootIndex,
                       "Per-molecule index of the root atom for Z-matrix tree traversal.")

        // -------------------------------------------------------------------------
        // Atom arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_atoms", &SystemTopology::numAtoms, "Total number of atoms.")
        .def_readwrite("atoms_unique_name",
                       &SystemTopology::atomsUniqueName,
                       "Unique atom name string (e.g. \"ALA1_CA_3\").")
        .def_readwrite("atoms_class_names",
                       &SystemTopology::atomsClassNames,
                       "Force-field atom class name for each atom.")
        .def_readwrite("atoms_charged_type_names",
                       &SystemTopology::atomsChargedTypeNames,
                       "Charged atom type name for each atom.")

        .def_readwrite("atoms_compound_atom_index",
                       &SystemTopology::atomsCompoundAtomIndex,
                       "Index of the atom within its compound (residue) for each atom.")
        .def_readwrite("atoms_class_index",
                       &SystemTopology::atomsAtomClassIndex,
                       "Index into the nonbonded parameter table for this atom's class.")
        .def_readwrite("atoms_charged_atom_type_index",
                       &SystemTopology::atomsChargedAtomTypeIndex,
                       "Index into the nonbonded parameter table for this atom's charged type.")

        .def_readwrite("atoms_nonbonded_index",
                       &SystemTopology::atomsNonbondedIndex,
                       "Index into the nonbonded parameter table for this atom.")
        .def_readwrite("atoms_element_name",
                       &SystemTopology::atomsElementName,
                       "Full element name string (e.g. \"Carbon\").")
        .def_readwrite("atoms_element_symbol",
                       &SystemTopology::atomsElementSymbol,
                       "Element symbol string (e.g. \"C\").")
        .def_readwrite("atoms_mass", &SystemTopology::atomsMass, "Mass in daltons [Da].")
        .def_readwrite("atoms_charge",
                       &SystemTopology::atomsCharge,
                       "Partial charge in units of the proton charge [e].")
        .def_readwrite("atoms_sigma",
                       &SystemTopology::atomsSigma,
                       "Lennard-Jones sigma (van der Waals radius) in nanometers [nm].")
        .def_readwrite("atoms_epsilon",
                       &SystemTopology::atomsEpsilon,
                       "Lennard-Jones epsilon (well depth) in kilojoules per mole [kJ/mol].")
        .def_readwrite("atoms_radius",
                       &SystemTopology::atomsRadius,
                       "GBSA implicit-solvent radius in nanometers [nm].")
        .def_readwrite("atoms_screen", &SystemTopology::atomsScreen, "OBC screening factor (dimensionless).")
        .def_readwrite("atoms_x",
                       &SystemTopology::atomsX,
                       "x coordinate of the reference structure in nanometers [nm].")
        .def_readwrite("atoms_y",
                       &SystemTopology::atomsY,
                       "y coordinate of the reference structure in nanometers [nm].")
        .def_readwrite("atoms_z",
                       &SystemTopology::atomsZ,
                       "z coordinate of the reference structure in nanometers [nm].")
        .def_readwrite("atoms_atomic_number",
                       &SystemTopology::atomsAtomicNumber,
                       "Atomic number (number of protons in the nucleus) for each atom.")
        .def_readwrite(
            "atoms_num_bonds_involved",
            &SystemTopology::atomsNumBondsInvolved,
            "Number of bonds this atom is involved in (used for determining terminal vs. internal atoms).")
        .def_readwrite("root_mobilities",
                       &SystemTopology::rootMobilities,
                       "Root mobility for each molecule, used to determine the mobility of the root atom in "
                       "each molecule.")

        // -------------------------------------------------------------------------
        // Bond arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_bonds",
                       &SystemTopology::numBonds,
                       "Total number of bonds (tree bonds plus ring-closing bonds).")
        .def_readwrite("bonds_i", &SystemTopology::bondsI, "BFS index of bond endpoint atom 1.")
        .def_readwrite("bonds_j", &SystemTopology::bondsJ, "BFS index of bond endpoint atom 2.")
        .def_readwrite("bonds_molecule_index",
                       &SystemTopology::bondsMoleculeIndex,
                       "Index of the molecule this bond belongs to.")
        .def_readwrite("bonds_ring_closing",
                       &SystemTopology::bondsRingClosing,
                       "True if this bond closes a ring (i.e. is not a tree bond).")
        .def_readwrite(
            "bonds_stiffness",
            &SystemTopology::bondsStiffness,
            "Harmonic bond force constant in kilojoules per mole per nanometer squared [kJ/mol/nm^2].")
        .def_readwrite("bonds_equilibrium",
                       &SystemTopology::bondsEquilibrium,
                       "Equilibrium bond length in nanometers [nm].")

        // -------------------------------------------------------------------------
        // Angle arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_angles", &SystemTopology::numAngles, "Total number of angles.")
        .def_readwrite("angles_i", &SystemTopology::anglesI, "Global index of angle atom 1 (outer).")
        .def_readwrite("angles_j", &SystemTopology::anglesJ, "Global index of angle atom 2 (central).")
        .def_readwrite("angles_k", &SystemTopology::anglesK, "Global index of angle atom 3 (outer).")
        .def_readwrite("angles_equilibrium",
                       &SystemTopology::anglesEquilibrium,
                       "Equilibrium bond angle in radians [rad].")
        .def_readwrite(
            "angles_stiffness",
            &SystemTopology::anglesStiffness,
            "Harmonic angle force constant in kilojoules per mole per radian squared [kJ/mol/rad^2].")

        // -------------------------------------------------------------------------
        // Periodic torsion arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_periodic_torsions",
                       &SystemTopology::numPeriodicTorsions,
                       "Total number of periodic torsion terms.")
        .def_readwrite("periodic_torsions_improper",
                       &SystemTopology::periodicTorsionsImproper,
                       "True if the torsion is an improper (out-of-plane) term.")
        .def_readwrite("periodic_torsions_i",
                       &SystemTopology::periodicTorsionsI,
                       "Global index of torsion atom 1.")
        .def_readwrite("periodic_torsions_j",
                       &SystemTopology::periodicTorsionsJ,
                       "Global index of torsion atom 2.")
        .def_readwrite("periodic_torsions_k",
                       &SystemTopology::periodicTorsionsK,
                       "Global index of torsion atom 3.")
        .def_readwrite("periodic_torsions_l",
                       &SystemTopology::periodicTorsionsL,
                       "Global index of torsion atom 4.")
        .def_readwrite("periodic_torsions_n",
                       &SystemTopology::periodicTorsionsN,
                       "Torsion periodicity (integer, dimensionless).")
        .def_readwrite("periodic_torsions_phase",
                       &SystemTopology::periodicTorsionsPhase,
                       "Phase offset in radians [rad].")
        .def_readwrite("periodic_torsions_stiffness",
                       &SystemTopology::periodicTorsionsStiffness,
                       "Force constant in kilojoules per mole [kJ/mol].")

        // -------------------------------------------------------------------------
        // Harmonic torsion arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_harmonic_torsions",
                       &SystemTopology::numHarmonicTorsions,
                       "Total number of harmonic torsion terms.")
        .def_readwrite("harmonic_torsions_i",
                       &SystemTopology::harmonicTorsionsI,
                       "Global index of torsion atom 1.")
        .def_readwrite("harmonic_torsions_j",
                       &SystemTopology::harmonicTorsionsJ,
                       "Global index of torsion atom 2.")
        .def_readwrite("harmonic_torsions_k",
                       &SystemTopology::harmonicTorsionsK,
                       "Global index of torsion atom 3.")
        .def_readwrite("harmonic_torsions_l",
                       &SystemTopology::harmonicTorsionsL,
                       "Global index of torsion atom 4.")
        .def_readwrite("harmonic_torsions_stiffness",
                       &SystemTopology::harmonicTorsionsStiffness,
                       "Force constant in kilojoules per mole [kJ/mol].")
        .def_readwrite("harmonic_torsions_phase",
                       &SystemTopology::harmonicTorsionsPhase,
                       "Equilibrium dihedral angle in radians [rad].")

        // -------------------------------------------------------------------------
        // Z-matrix arrays
        // -------------------------------------------------------------------------

        .def_readwrite("num_z_matrix_rows",
                       &SystemTopology::numZMatrixRows,
                       "Number of Z-matrix rows (equals number of atoms).")
        .def_readwrite("z_matrix_i",
                       &SystemTopology::zMatrixI,
                       "Global atom index of the atom placed at row r.")
        .def_readwrite(
            "z_matrix_j",
            &SystemTopology::zMatrixJ,
            "Bond-length reference atom at row r. Row 0 holds sentinel -1 (root has no bond reference).")
        .def_readwrite("z_matrix_k",
                       &SystemTopology::zMatrixK,
                       "Bond-angle reference atom at row r. Rows 0-1 hold sentinel -1.")
        .def_readwrite("z_matrix_l",
                       &SystemTopology::zMatrixL,
                       "Dihedral reference atom at row r. Rows 0-2 hold sentinel -1.")

        // -------------------------------------------------------------------------
        // Urey-Bradley 1-3 interactions
        // -------------------------------------------------------------------------

        .def_readwrite("num_urey_bradley",
                       &SystemTopology::numUreyBradley,
                       "Total number of Urey-Bradley 1-3 terms.")
        .def_readwrite("urey_bradley_i",
                       &SystemTopology::ureyBradleyI,
                       "Global index of atom 1 (outer atom of angle i-j-k).")
        .def_readwrite("urey_bradley_k",
                       &SystemTopology::ureyBradleyK,
                       "Global index of atom 3 (outer atom of angle i-j-k).")
        .def_readwrite("urey_bradley_stiffness",
                       &SystemTopology::ureyBradleyStiffness,
                       "Harmonic force constant in kilojoules per mole per nanometer squared [kJ/mol/nm^2].")
        .def_readwrite("urey_bradley_equilibrium",
                       &SystemTopology::ureyBradleyEquilibrium,
                       "Nominal 1-3 distance in nanometers [nm].")

        // -------------------------------------------------------------------------
        // 1-4 pair scaling
        // -------------------------------------------------------------------------

        .def_readwrite("num_scaling14",
                       &SystemTopology::numScaling14,
                       "Total number of 1-4 pair scaling terms.")
        .def_readwrite("scaling14_i",
                       &SystemTopology::scaling14I,
                       "Global index of atom 1 (first atom of dihedral i-j-k-l).")
        .def_readwrite("scaling14_l",
                       &SystemTopology::scaling14L,
                       "Global index of atom 4 (last atom of dihedral i-j-k-l).")
        .def_readwrite("scaling14_charge_product",
                       &SystemTopology::scaling14ChargeProduct,
                       "q1 times q4 pre-scaled by the 1-4 electrostatic factor, in units of proton charge "
                       "squared [e^2].")
        .def_readwrite("scaling14_epsilon",
                       &SystemTopology::scaling14Epsilon,
                       "Combined Lennard-Jones well depth in kilojoules per mole [kJ/mol].")
        .def_readwrite("scaling14_sigma",
                       &SystemTopology::scaling14Sigma,
                       "Combined Lennard-Jones radius in nanometers [nm].")

        // -------------------------------------------------------------------------
        // Exclusions
        // -------------------------------------------------------------------------

        .def_readwrite("num_exclusions",
                       &SystemTopology::numExclusions,
                       "Total number of non-bonded exclusion pairs.")
        .def_readwrite("exclusion_i",
                       &SystemTopology::exclusionI,
                       "Global index of atom 1 of the excluded pair.")
        .def_readwrite("exclusion_j",
                       &SystemTopology::exclusionJ,
                       "Global index of atom 2 of the excluded pair.")

        // -------------------------------------------------------------------------
        // CMAP torsion arrays
        // -------------------------------------------------------------------------

        .def_readwrite("cmap_grid_size",
                       &SystemTopology::cmapGridSize,
                       "Dimension of each CMAP grid (grid has cmap_grid_size x cmap_grid_size points).")
        .def_readwrite("cmap_grid_energy",
                       &SystemTopology::cmapGridEnergy,
                       "Flattened CMAP energy grid in row-major order in kilojoules per mole [kJ/mol].")
        .def_readwrite("cmap_torsion_map_index",
                       &SystemTopology::cmapTorsionMapIndex,
                       "Index into the CMAP grid table for each torsion pair.")
        .def_readwrite("cmap_torsion_a1", &SystemTopology::cmapTorsionA1, "Global index of torsion A atom 1.")
        .def_readwrite("cmap_torsion_a2", &SystemTopology::cmapTorsionA2, "Global index of torsion A atom 2.")
        .def_readwrite("cmap_torsion_a3", &SystemTopology::cmapTorsionA3, "Global index of torsion A atom 3.")
        .def_readwrite("cmap_torsion_a4", &SystemTopology::cmapTorsionA4, "Global index of torsion A atom 4.")
        .def_readwrite("cmap_torsion_b1", &SystemTopology::cmapTorsionB1, "Global index of torsion B atom 1.")
        .def_readwrite("cmap_torsion_b2", &SystemTopology::cmapTorsionB2, "Global index of torsion B atom 2.")
        .def_readwrite("cmap_torsion_b3", &SystemTopology::cmapTorsionB3, "Global index of torsion B atom 3.")
        .def_readwrite("cmap_torsion_b4", &SystemTopology::cmapTorsionB4, "Global index of torsion B atom 4.")

        // -------------------------------------------------------------------------
        // NBfix pairwise corrections
        // -------------------------------------------------------------------------

        .def_readwrite("has_nbfix",
                       &SystemTopology::hasNBfix,
                       "True if NBfix pairwise corrections are present.")
        .def_readwrite("num_nb_types",
                       &SystemTopology::numNBTypes,
                       "Number of distinct nonbonded atom types.")
        .def_readwrite("a_coef",
                       &SystemTopology::aCoef,
                       "Lennard-Jones A coefficients (repulsive term) for each type pair, stored as a flat "
                       "upper-triangular matrix in kJ/mol*nm^12.")
        .def_readwrite("b_coef",
                       &SystemTopology::bCoef,
                       "Lennard-Jones B coefficients (attractive term) for each type pair, stored as a flat "
                       "upper-triangular matrix in kJ/mol*nm^6.")

        // -------------------------------------------------------------------------
        // GBSA implicit solvent
        // -------------------------------------------------------------------------

        .def_readwrite("use_gbsa_obc2",
                       &SystemTopology::useGBSAOBC2,
                       "True if the GBSA OBC2 implicit solvent model is active.")
        .def_readwrite(
            "gbsa_solvent_dielectric",
            &SystemTopology::gbsaSolventDielectric,
            "Relative dielectric constant of the solvent (dimensionless). Default 78.5 (water at 298 K).")
        .def_readwrite("gbsa_solute_dielectric",
                       &SystemTopology::gbsaSoluteDielectric,
                       "Relative dielectric constant of the solute interior (dimensionless). Default 1.0.")

        // -------------------------------------------------------------------------
        // Nonbonded settings
        // -------------------------------------------------------------------------

        .def_readwrite("nonbonded_method",
                       &SystemTopology::nonbondedMethod,
                       "Algorithm used to evaluate nonbonded interactions.")
        .def_readwrite("nonbonded_cutoff",
                       &SystemTopology::nonbondedCutoff,
                       "Distance cutoff for nonbonded interactions in nanometers [nm].")

        // -------------------------------------------------------------------------
        // Thermostat
        // -------------------------------------------------------------------------

        .def_readwrite("thermostat_temperature",
                       &SystemTopology::thermostatTemperature,
                       "Target temperature for the Langevin thermostat in Kelvin [K].")
        .def_readwrite("collision_frequency",
                       &SystemTopology::collisionFrequency,
                       "Langevin collision frequency (friction coefficient) in inverse picoseconds [ps^-1].")
        .def_readwrite("seed",
                       &SystemTopology::seed,
                       "Random number seed for the integrator. 0 means system-generated.");
}
