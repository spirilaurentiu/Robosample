#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Context.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc() = "Robosample bindings";



    py::enum_<ROOT_MOBILITY>(m, "RootMobility")
        .value("FREE", ROOT_MOBILITY::FREE)
        .value("CARTESIAN", ROOT_MOBILITY::CARTESIAN)
        .value("WELD", ROOT_MOBILITY::WELD)
        .value("FREE_LINE", ROOT_MOBILITY::FREE_LINE)
        .value("BALL", ROOT_MOBILITY::BALL)
        .value("PIN", ROOT_MOBILITY::PIN);

    py::enum_<BondMobility::Mobility>(m, "BondMobility")
        .value("Free", BondMobility::Mobility::Free)
        .value("Torsion", BondMobility::Mobility::Torsion)
        .value("Rigid", BondMobility::Mobility::Rigid)
        .value("BallF", BondMobility::Mobility::BallF)
        .value("BallM", BondMobility::Mobility::BallM)
        .value("Cylinder", BondMobility::Mobility::Cylinder)
        .value("Translation", BondMobility::Mobility::Translation)
        .value("FreeLine", BondMobility::Mobility::FreeLine)
        .value("LineOrientationF", BondMobility::Mobility::LineOrientationF)
        .value("LineOrientationM", BondMobility::Mobility::LineOrientationM)
        .value("UniversalM", BondMobility::Mobility::UniversalM)
        .value("Spherical", BondMobility::Mobility::Spherical)
        .value("AnglePin", BondMobility::Mobility::AnglePin)
        .value("BendStretch", BondMobility::Mobility::BendStretch)
        .value("Slider", BondMobility::Mobility::Slider)
        .value("OrthoSpherical", BondMobility::Mobility::OrthoSpherical);

    py::enum_<RUN_TYPE>(m, "RunType")
        .value("DEFAULT", RUN_TYPE::DEFAULT)
        .value("REMC", RUN_TYPE::REMC)
        .value("RENEMC", RUN_TYPE::RENEMC)
        .value("RENE", RUN_TYPE::RENE);

    py::enum_<SamplerName>(m, "SamplerName")
        .value("EMPTY", SamplerName::EMPTY)
        .value("MC", SamplerName::MC)
        .value("HMC", SamplerName::HMC)
        .value("LAHMC", SamplerName::LAHMC);

    py::enum_<AcceptRejectMode>(m, "AcceptRejectMode")
        .value("AlwaysAccept", AcceptRejectMode::AlwaysAccept)
        .value("MetropolisHastings", AcceptRejectMode::MetropolisHastings);

    py::enum_<IntegratorType>(m, "IntegratorType")
        .value("EMPTY", IntegratorType::EMPTY)
        .value("VERLET", IntegratorType::VERLET)
        .value("EULER", IntegratorType::EULER)
        .value("EULER2", IntegratorType::EULER2)
        .value("CPODES", IntegratorType::CPODES)
        .value("RUNGEKUTTA", IntegratorType::RUNGEKUTTA)
        .value("RUNGEKUTTA2", IntegratorType::RUNGEKUTTA2)
        .value("RUNGEKUTTA3", IntegratorType::RUNGEKUTTA3)
        .value("RUNGEKUTTAFELDBERG", IntegratorType::RUNGEKUTTAFELDBERG)
        .value("BENDSTRETCH", IntegratorType::BENDSTRETCH)
        .value("OMMVV", IntegratorType::OMMVV)
        .value("BOUND_WALK", IntegratorType::BOUND_WALK)
        .value("BOUND_HMC", IntegratorType::BOUND_HMC)
        .value("STATIONS_TASK", IntegratorType::STATIONS_TASK)
        .value("NOF_INTEGRATORS", IntegratorType::NOF_INTEGRATORS);

    py::enum_<ThermostatName>(m, "ThermostatName")
        .value("NONE", ThermostatName::NONE)
        .value("ANDERSEN", ThermostatName::ANDERSEN)
        .value("BERENDSEN", ThermostatName::BERENDSEN)
        .value("LANGEVIN", ThermostatName::LANGEVIN)
        .value("NOSE_HOOVER", ThermostatName::NOSE_HOOVER);

    py::class_<BOND_FLEXIBILITY>(m, "BondFlexibility")
        .def(py::init<>())
        .def(py::init<int, int, BondMobility::Mobility>())
        .def_readwrite("i", &BOND_FLEXIBILITY::i)
        .def_readwrite("j", &BOND_FLEXIBILITY::j)
        .def_readwrite("mobility", &BOND_FLEXIBILITY::mobility);

    py::class_<AtomDefinition>(m, "AtomDefinition")
        .def(py::init<>())
        .def_readwrite("global_index", &AtomDefinition::globalIndex)
        .def_readwrite("prmtop_index", &AtomDefinition::prmtopIndex)
        .def_readwrite("molecule_index", &AtomDefinition::moleculeIndex)
        .def_readwrite("residue_index", &AtomDefinition::residueIndex)
        .def_readwrite("atom_class_name", &AtomDefinition::atomClassName)
        .def_readwrite("atom_class_index", &AtomDefinition::atomClassIndex)
        .def_readwrite("charged_atom_type_name", &AtomDefinition::chargedAtomTypeName)
        .def_readwrite("charged_atom_type_index", &AtomDefinition::chargedAtomTypeIndex)
        .def_readwrite("residue_name", &AtomDefinition::residueName)
        .def_readwrite("unique_atom_name", &AtomDefinition::uniqueAtomName)
        .def_readwrite("neighbors_global_indices", &AtomDefinition::neighborsGlobalIndices)
        .def_readwrite("root", &AtomDefinition::root)
        .def_readwrite("atomic_number", &AtomDefinition::atomicNumber)
        .def_readwrite("charge_in_e", &AtomDefinition::chargeInE)
        .def_readwrite("mass_in_daltons", &AtomDefinition::massInDaltons)
        .def_readwrite("vdw_radius_nm", &AtomDefinition::vdwRadiusInNm)
        .def_readwrite("sigma_nm", &AtomDefinition::sigmaInNm)
        .def_readwrite("vdw_well_depth_kj", &AtomDefinition::vdwWellDepthInKJ)
        .def_readwrite("x_nm", &AtomDefinition::x_nm)
        .def_readwrite("y_nm", &AtomDefinition::y_nm)
        .def_readwrite("z_nm", &AtomDefinition::z_nm);

    py::class_<BondStretchDefinition>(m, "BondStretchDefinition")
        .def(py::init<>())
        .def_readwrite("parentAtomGlobalIndex", &BondStretchDefinition::parentAtomGlobalIndex)
        .def_readwrite("childAtomGlobalIndex", &BondStretchDefinition::childAtomGlobalIndex)
        .def_readwrite("bondGlobalIndex", &BondStretchDefinition::bondGlobalIndex)
        .def_readwrite("moleculeIndex", &BondStretchDefinition::moleculeIndex)
        .def_readwrite("ringClosing", &BondStretchDefinition::ringClosing)
        .def_readwrite("stiffnessInKJPerNmSq", &BondStretchDefinition::stiffnessInKJPerNmSq)
        .def_readwrite("nominalLengthInNm", &BondStretchDefinition::nominalLengthInNm);

    py::class_<BondBendDefinition>(m, "BondBendDefinition")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &BondBendDefinition::globalIndex1)
        .def_readwrite("globalIndex2", &BondBendDefinition::globalIndex2)
        .def_readwrite("globalIndex3", &BondBendDefinition::globalIndex3)
        .def_readwrite("stiffnessInKJPerRadSq", &BondBendDefinition::stiffnessInKJPerRadSq)
        .def_readwrite("nominalAngleInDeg", &BondBendDefinition::nominalAngleInDeg);

    py::class_<BondTorsionDefinition>(m, "BondTorsionDefinition")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &BondTorsionDefinition::globalIndex1)
        .def_readwrite("globalIndex2", &BondTorsionDefinition::globalIndex2)
        .def_readwrite("globalIndex3", &BondTorsionDefinition::globalIndex3)
        .def_readwrite("globalIndex4", &BondTorsionDefinition::globalIndex4)
        .def_readwrite("improper", &BondTorsionDefinition::improper)
        .def_readwrite("ampInKJ_1", &BondTorsionDefinition::ampInKJ_1)
        .def_readwrite("phaseInDegrees_1", &BondTorsionDefinition::phaseInDegrees_1)
        .def_readwrite("periodicity_1", &BondTorsionDefinition::periodicity_1)
        .def_readwrite("ampInKJ_2", &BondTorsionDefinition::ampInKJ_2)
        .def_readwrite("phaseInDegrees_2", &BondTorsionDefinition::phaseInDegrees_2)
        .def_readwrite("periodicity_2", &BondTorsionDefinition::periodicity_2)
        .def_readwrite("ampInKJ_3", &BondTorsionDefinition::ampInKJ_3)
        .def_readwrite("phaseInDegrees_3", &BondTorsionDefinition::phaseInDegrees_3)
        .def_readwrite("periodicity_3", &BondTorsionDefinition::periodicity_3)
        .def_readwrite("ampInKJ_4", &BondTorsionDefinition::ampInKJ_4)
        .def_readwrite("phaseInDegrees_4", &BondTorsionDefinition::phaseInDegrees_4)
        .def_readwrite("periodicity_4", &BondTorsionDefinition::periodicity_4)
        .def_readwrite("ampInKJ_5", &BondTorsionDefinition::ampInKJ_5)
        .def_readwrite("phaseInDegrees_5", &BondTorsionDefinition::phaseInDegrees_5)
        .def_readwrite("periodicity_5", &BondTorsionDefinition::periodicity_5);

    py::class_<Atom>(m, "Atom")
        .def(py::init<const AtomDefinition&>(), py::arg("spec"));

    py::class_<BondStretch>(m, "BondStretch")
        .def(py::init<const BondStretchDefinition&>(), py::arg("spec"));
        
    py::class_<BondBend>(m, "BondBend")
        .def(py::init<const BondBendDefinition&>(), py::arg("spec"));

    py::class_<BondTorsion>(m, "BondTorsion")
        .def(py::init<const BondTorsionDefinition&>(), py::arg("spec"));

    py::class_<Context>(m, "Context")
        .def(py::init<const std::string&, uint32_t, uint32_t, uint32_t, RUN_TYPE, uint32_t, uint32_t>())
        .def("addReplica", &Context::addReplica, "Add an empty replica to the context.")
        .def("addThermodynamicState", &Context::addThermodynamicState, "Add an empty themodynamic state to the context.")
        .def("Initialize", py::overload_cast<>(&Context::Initialize), "Initializes the context after all worlds and replicas have been set.")
        .def("RunREX", &Context::RunREX, "Run replica exchange.")
        .def("setVerbose", &Context::setVerbose, "Control if you want extraneous output to cout.")
        .def("setPdbRestartFreq", &Context::setPdbRestartFreq, "Set the PDB restart frequency.")
        .def("setPrintFreq", &Context::setPrintFreq, "Set the print frequency.")
        .def("setNonbonded", &Context::setNonbonded, "Set nonbonded method and cutoff.")
        .def("setGBSA", &Context::setGBSA, "Set GBSA.")
        .def("loadAmberSystem", &Context::loadAmberSystem, "Load an Amber system.")
        .def("addWorld", &Context::addWorld, "Add an empty world.")
        .def("getWorld", (World& (Context::*)(std::size_t which)) &Context::getWorld, py::return_value_policy::reference, "Gets the world at the specified index.");

    py::class_<World>(m, "World")
        .def("setFlexibilities", &World::setFlexibilities, "Set the flexibilities of the bonds.")
        .def("addSampler", &World::addSampler, "Add a sampler to the world.")
        .def("setRollFlexibilities", &World::setRollFlexibilities, "Set rolling sequence.")
        .def("getRollFlexibilities", &World::getRollFlexibilities, "Get rolling sequence.");
}
