#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Context.hpp"
#include "OpenMMContext.hpp"
#include "TopologyElements.hpp"
#include "World.hpp"

namespace py = pybind11;

PYBIND11_MODULE(robo_bindings, m) {
    m.doc() = "Robosample C++ bindings (robo_bindings)";

    // -- Enums ---------------------------------------------------------------
    py::enum_<RootMobility>(m, "RootMobility")
        .value("FREE", RootMobility::Free)
        .value("CARTESIAN", RootMobility::Cartesian)
        .value("WELD", RootMobility::Weld)
        .value("FREE_LINE", RootMobility::FreeLine)
        .value("BALL", RootMobility::Ball)
        .value("PIN", RootMobility::Pin);

    py::enum_<NonbondedMethod>(m, "NonbondedMethod")
        .value("NoCutoff", NonbondedMethod::NoCutoff)
        .value("CutoffNonPeriodic", NonbondedMethod::CutoffNonPeriodic)
        .value("CutoffPeriodic", NonbondedMethod::CutoffPeriodic)
        .value("Ewald", NonbondedMethod::Ewald)
        .value("PME", NonbondedMethod::PME);

    py::enum_<AcceptRejectMode>(m, "AcceptRejectMode")
        .value("AlwaysAccept", AcceptRejectMode::AlwaysAccept)
        .value("MetropolisHastings", AcceptRejectMode::MetropolisHastings);

    // Proposal type of a world's sampler (a docking world auto-selects RigidKick).
    py::enum_<MoveType>(m, "MoveType")
        .value("MdHmc", MoveType::MdHmc)
        .value("RigidKick", MoveType::RigidKick);

    py::enum_<BondMobility>(m, "BondMobility")
        .value("Rigid", BondMobility::Rigid)
        .value("Torsion", BondMobility::Torsion)
        .value("Free", BondMobility::Free)
        .value("Ball", BondMobility::Ball)
        .value("Pin", BondMobility::Pin)
        .value("Slider", BondMobility::Slider)
        .value("Cylinder", BondMobility::Cylinder)
        .value("BendStretch", BondMobility::BendStretch);

    py::class_<Selection>(m, "Selection").def(py::init<>());

    // A flexibility world. add_sampler returns *this so Python can chain. Note:
    // boostMDSteps was removed (it was never read). New optional kwargs:
    //   sphere_radius : docking binding-sphere radius [nm] (RigidKick worlds).
    //   use_fixman    : include the Fixman potential + logSineSqr in the
    //                   acceptance H of torsional (internal-coordinate) worlds.
    py::class_<World>(m, "World")
        .def("add_sampler",
             &World::add_sampler,
             py::arg("timeStep"),
             py::arg("mdSteps"),
             py::arg("acceptRejectMode"),
             py::arg("use_nuts"),
             py::arg("sphere_factor") = 1.0,     // auto R = R_receptor + factor*R_ligand
             py::arg("use_fixman") = py::none(), // None => auto: on for non-Cartesian
             py::arg("always_kick") = false,     // perturb every round (else: only when COM left sphere)
             py::arg("clash_threshold") = 10.0,  // reject if |peNew|>factor*|pePre| (relative, dimensionless)
             py::return_value_policy::reference,
             "Configure this world's sampler; returns the world for chaining. "
             "sphere_factor scales the auto-sized per-ligand binding sphere "
             "(R = R_receptor + sphere_factor*R_ligand). The docking kick relocates a "
             "ligand only when its COM leaves the sphere (always_kick=True perturbs every "
             "round). A proposed pose is rejected -- in ALL modes, including AlwaysAccept "
             "-- if its potential energy is non-finite or |PE| exceeds clash_threshold, so "
             "overlapping geometry never passes. use_fixman=None auto-enables "
             "Fixman+logSineSqr on non-Cartesian worlds.")
        .def_property_readonly("index", &World::index)
        .def_property_readonly("is_cartesian", &World::isCartesian)
        .def_property_readonly("is_docking", &World::isDocking);

    py::class_<OpenMMContext::ForceGroupEnergy>(m, "ForceGroupEnergy")
        .def_readonly("group", &OpenMMContext::ForceGroupEnergy::group)
        .def_readonly("name", &OpenMMContext::ForceGroupEnergy::name)
        .def_readonly("energy", &OpenMMContext::ForceGroupEnergy::energy);

    // -- SystemTopology (SoA) ------------------------------------------------
    py::class_<SystemTopology>(m, "SystemTopology")
        .def(py::init<>())
        .def_readwrite("num_molecules", &SystemTopology::numMolecules)
        .def_readwrite("atoms_begin", &SystemTopology::atomsBegin)
        .def_readwrite("atoms_end", &SystemTopology::atomsEnd)
        .def_readwrite("bonds_begin", &SystemTopology::bondsBegin)
        .def_readwrite("bonds_end", &SystemTopology::bondsEnd)
        .def_readwrite("angles_begin", &SystemTopology::anglesBegin)
        .def_readwrite("angles_end", &SystemTopology::anglesEnd)
        .def_readwrite("periodic_torsions_begin", &SystemTopology::periodicTorsionsBegin)
        .def_readwrite("periodic_torsions_end", &SystemTopology::periodicTorsionsEnd)
        .def_readwrite("harmonic_torsions_begin", &SystemTopology::harmonicTorsionsBegin)
        .def_readwrite("harmonic_torsions_end", &SystemTopology::harmonicTorsionsEnd)
        .def_readwrite("z_matrix_begin", &SystemTopology::zMatrixBegin)
        .def_readwrite("z_matrix_end", &SystemTopology::zMatrixEnd)
        .def_readwrite("urey_bradley_begin", &SystemTopology::ureyBradleyBegin)
        .def_readwrite("urey_bradley_end", &SystemTopology::ureyBradleyEnd)
        .def_readwrite("scaling14_begin", &SystemTopology::scaling14Begin)
        .def_readwrite("scaling14_end", &SystemTopology::scaling14End)
        .def_readwrite("exclusion_begin", &SystemTopology::exclusionBegin)
        .def_readwrite("exclusion_end", &SystemTopology::exclusionEnd)
        .def_readwrite("atoms_root_index", &SystemTopology::atomsRootIndex)
        .def_readwrite("root_mobilities", &SystemTopology::rootMobilities)
        .def_readwrite("num_atoms", &SystemTopology::numAtoms)
        .def_readwrite("atoms_unique_name", &SystemTopology::atomsUniqueName)
        .def_readwrite("atoms_nonbonded_index", &SystemTopology::atomsNonbondedIndex)
        .def_readwrite("atoms_prmtop_index", &SystemTopology::atomsPrmtopIndex)
        .def_readwrite("atoms_element_name", &SystemTopology::atomsElementName)
        .def_readwrite("atoms_element_symbol", &SystemTopology::atomsElementSymbol)
        .def_readwrite("atoms_mass", &SystemTopology::atomsMass)
        .def_readwrite("atoms_charge", &SystemTopology::atomsCharge)
        .def_readwrite("atoms_sigma", &SystemTopology::atomsSigma)
        .def_readwrite("atoms_epsilon", &SystemTopology::atomsEpsilon)
        .def_readwrite("atoms_radius", &SystemTopology::atomsRadius)
        .def_readwrite("atoms_screen", &SystemTopology::atomsScreen)
        .def_readwrite("atoms_x", &SystemTopology::atomsX)
        .def_readwrite("atoms_y", &SystemTopology::atomsY)
        .def_readwrite("atoms_z", &SystemTopology::atomsZ)
        .def_readwrite("atoms_atomic_number", &SystemTopology::atomsAtomicNumber)
        .def_readwrite("atoms_num_bonds_involved", &SystemTopology::atomsNumBondsInvolved)
        .def_readwrite("num_bonds", &SystemTopology::numBonds)
        .def_readwrite("bonds_i", &SystemTopology::bondsI)
        .def_readwrite("bonds_j", &SystemTopology::bondsJ)
        .def_readwrite("bonds_molecule_index", &SystemTopology::bondsMoleculeIndex)
        .def_readwrite("bonds_ring_closing", &SystemTopology::bondsRingClosing)
        .def_readwrite("bonds_stiffness", &SystemTopology::bondsStiffness)
        .def_readwrite("bonds_equilibrium", &SystemTopology::bondsEquilibrium)
        .def_readwrite("num_angles", &SystemTopology::numAngles)
        .def_readwrite("angles_i", &SystemTopology::anglesI)
        .def_readwrite("angles_j", &SystemTopology::anglesJ)
        .def_readwrite("angles_k", &SystemTopology::anglesK)
        .def_readwrite("angles_equilibrium", &SystemTopology::anglesEquilibrium)
        .def_readwrite("angles_stiffness", &SystemTopology::anglesStiffness)
        .def_readwrite("num_periodic_torsions", &SystemTopology::numPeriodicTorsions)
        .def_readwrite("periodic_torsions_improper", &SystemTopology::periodicTorsionsImproper)
        .def_readwrite("periodic_torsions_i", &SystemTopology::periodicTorsionsI)
        .def_readwrite("periodic_torsions_j", &SystemTopology::periodicTorsionsJ)
        .def_readwrite("periodic_torsions_k", &SystemTopology::periodicTorsionsK)
        .def_readwrite("periodic_torsions_l", &SystemTopology::periodicTorsionsL)
        .def_readwrite("periodic_torsions_n", &SystemTopology::periodicTorsionsN)
        .def_readwrite("periodic_torsions_phase", &SystemTopology::periodicTorsionsPhase)
        .def_readwrite("periodic_torsions_stiffness", &SystemTopology::periodicTorsionsStiffness)
        .def_readwrite("num_harmonic_torsions", &SystemTopology::numHarmonicTorsions)
        .def_readwrite("harmonic_torsions_i", &SystemTopology::harmonicTorsionsI)
        .def_readwrite("harmonic_torsions_j", &SystemTopology::harmonicTorsionsJ)
        .def_readwrite("harmonic_torsions_k", &SystemTopology::harmonicTorsionsK)
        .def_readwrite("harmonic_torsions_l", &SystemTopology::harmonicTorsionsL)
        .def_readwrite("harmonic_torsions_stiffness", &SystemTopology::harmonicTorsionsStiffness)
        .def_readwrite("harmonic_torsions_phase", &SystemTopology::harmonicTorsionsPhase)
        .def_readwrite("num_z_matrix_rows", &SystemTopology::numZMatrixRows)
        .def_readwrite("z_matrix_i", &SystemTopology::zMatrixI)
        .def_readwrite("z_matrix_j", &SystemTopology::zMatrixJ)
        .def_readwrite("z_matrix_k", &SystemTopology::zMatrixK)
        .def_readwrite("z_matrix_l", &SystemTopology::zMatrixL)
        .def_readwrite("num_urey_bradley", &SystemTopology::numUreyBradley)
        .def_readwrite("urey_bradley_i", &SystemTopology::ureyBradleyI)
        .def_readwrite("urey_bradley_k", &SystemTopology::ureyBradleyK)
        .def_readwrite("urey_bradley_stiffness", &SystemTopology::ureyBradleyStiffness)
        .def_readwrite("urey_bradley_equilibrium", &SystemTopology::ureyBradleyEquilibrium)
        .def_readwrite("num_scaling14", &SystemTopology::numScaling14)
        .def_readwrite("scaling14_i", &SystemTopology::scaling14I)
        .def_readwrite("scaling14_l", &SystemTopology::scaling14L)
        .def_readwrite("scaling14_charge_product", &SystemTopology::scaling14ChargeProduct)
        .def_readwrite("scaling14_epsilon", &SystemTopology::scaling14Epsilon)
        .def_readwrite("scaling14_sigma", &SystemTopology::scaling14Sigma)
        .def_readwrite("num_exclusions", &SystemTopology::numExclusions)
        .def_readwrite("exclusion_i", &SystemTopology::exclusionI)
        .def_readwrite("exclusion_j", &SystemTopology::exclusionJ)
        .def_readwrite("cmap_grid_size", &SystemTopology::cmapGridSize)
        .def_readwrite("cmap_grid_energy", &SystemTopology::cmapGridEnergy)
        .def_readwrite("cmap_torsion_map_index", &SystemTopology::cmapTorsionMapIndex)
        .def_readwrite("cmap_torsion_a1", &SystemTopology::cmapTorsionA1)
        .def_readwrite("cmap_torsion_a2", &SystemTopology::cmapTorsionA2)
        .def_readwrite("cmap_torsion_a3", &SystemTopology::cmapTorsionA3)
        .def_readwrite("cmap_torsion_a4", &SystemTopology::cmapTorsionA4)
        .def_readwrite("cmap_torsion_b1", &SystemTopology::cmapTorsionB1)
        .def_readwrite("cmap_torsion_b2", &SystemTopology::cmapTorsionB2)
        .def_readwrite("cmap_torsion_b3", &SystemTopology::cmapTorsionB3)
        .def_readwrite("cmap_torsion_b4", &SystemTopology::cmapTorsionB4)
        .def_readwrite("has_nb_fix", &SystemTopology::hasNBfix)
        .def_readwrite("num_nb_types", &SystemTopology::numNBTypes)
        .def_readwrite("a_coef", &SystemTopology::aCoef)
        .def_readwrite("b_coef", &SystemTopology::bCoef)
        .def_readwrite("use_gbsa_obc2", &SystemTopology::useGBSAOBC2)
        .def_readwrite("gbsa_solvent_dielectric", &SystemTopology::gbsaSolventDielectric)
        .def_readwrite("gbsa_solute_dielectric", &SystemTopology::gbsaSoluteDielectric)
        .def_readwrite("nonbonded_method", &SystemTopology::nonbondedMethod)
        .def_readwrite("nonbonded_cutoff", &SystemTopology::nonbondedCutoff)
        .def_readwrite("thermostat_temperature", &SystemTopology::thermostatTemperature)
        .def_readwrite("collision_frequency", &SystemTopology::collisionFrequency)
        .def_readwrite("seed", &SystemTopology::seed);

    // -- Context -------------------------------------------------------------
    py::class_<Context>(m, "Context")
        .def(py::init<std::string, std::uint32_t>(), py::arg("base_name"), py::arg("seed"))
        .def_readwrite("system_topology", &Context::systemTopology)
        .def("initialize",
             &Context::initialize,
             py::arg("temperatures") = std::vector<double>{},
             "Build OpenMM, set the replica temperature ladder, seed coordinates.")
        .def("add_cartesian_world",
             &Context::addCartesianWorld,
             py::return_value_policy::reference_internal,
             "Add a Cartesian (OpenMM-MD) world; returns it for .add_sampler(...).")
        .def("add_robotic_world",
             &Context::addRoboticWorld,
             py::arg("selection"),
             py::return_value_policy::reference_internal,
             "Add an internal-coordinate (torsional) world for the given selection.")
        .def("add_torsional_world",
             &Context::addRoboticWorld,
             py::arg("selection"),
             py::return_value_policy::reference_internal,
             "Alias of add_robotic_world.")
        .def("add_docking_world",
             &Context::addDockingWorld,
             py::arg("ligand_molecule_indices"),
             py::return_value_policy::reference_internal,
             "Add a rigid-body docking world: the listed ligand molecules get Free roots, "
             "every other molecule is welded/rigid. Chain .add_sampler(sphere_radius=...).")
        .def("build_flexibilities",
             &Context::buildFlexibilities,
             py::arg("bonds"),
             py::arg("mobility"),
             py::arg("flag"),
             "Build a per-bond mobility selection (bonds=None => all eligible).")
        .def("run_rex",
             &Context::runREX,
             py::arg("equil_rounds"),
             py::arg("prod_rounds"),
             py::arg("write_freq"),
             py::arg("verbose"),
             "Run replica exchange: Gibbs sweep over worlds + adjacent swaps.")
        .def("set_root_mobility",
             &Context::setRootMobility,
             py::arg("molecule_index"),
             py::arg("mobility"),
             "Override a molecule's root attachment to Ground.")
        .def("set_mts",
             &Context::setMTS,
             py::arg("enabled"),
             py::arg("inner_substeps") = 4,
             "Enable r-RESPA multiple-timestep OpenMM MD (Cartesian world): slow forces once per "
             "outer step, fast bonded forces inner_substeps times. Call before initialize().")
        .def("initialize_openmm",
             &Context::initializeOpenMM,
             "Build the single OpenMM System/Context from system_topology.")
        .def("calc_openmm_potential_energy",
             &Context::calcOpenMMPotentialEnergy,
             "Set the reference coordinates and return the OpenMM potential energy [kJ/mol].")
        .def(
            "set_separate_force_groups",
            [](Context& /*ctx*/, bool enabled) {
                OpenMMContext::get().setSeparateForceGroups(enabled);
            },
            py::arg("enabled"),
            "Enable/disable separate OpenMM force groups for each Force.")
        .def("calc_openmm_potential_energy_by_group",
             &Context::computePotentialEnergyByGroup,
             "Compute potential energy by OpenMM force group, returning (group, name, energy) tuples.");
}