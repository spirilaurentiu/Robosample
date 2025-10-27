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

    py::class_<AtomSpec>(m, "AtomSpec")
        .def(py::init<>())
        .def_readwrite("globalIndex", &AtomSpec::globalIndex)
        .def_readwrite("moleculeIndex", &AtomSpec::moleculeIndex)
        .def_readwrite("residueIndex", &AtomSpec::residueIndex)
        .def_readwrite("atomClassIndex", &AtomSpec::atomClassIndex)
        .def_readwrite("chargedAtomTypeIndex", &AtomSpec::chargedAtomTypeIndex)
        .def_readwrite("atomName", &AtomSpec::atomName)
        .def_readwrite("residueName", &AtomSpec::residueName)
        .def_readwrite("atomClassName", &AtomSpec::atomClassName)
        .def_readwrite("chargedAtomName", &AtomSpec::chargedAtomName)
        .def_readwrite("neighborsGlobalIndices", &AtomSpec::neighborsGlobalIndices)
        .def_readwrite("availableBonds", &AtomSpec::availableBonds)
        .def_readwrite("root", &AtomSpec::root)
        .def_readwrite("atomicNumber", &AtomSpec::atomicNumber)
        .def_readwrite("chargeInE", &AtomSpec::chargeInE)
        .def_readwrite("massInDaltons", &AtomSpec::massInDaltons)
        .def_readwrite("vdwRadiusInNm", &AtomSpec::vdwRadiusInNm)
        .def_readwrite("vdwWellDepthInKJ", &AtomSpec::vdwWellDepthInKJ)
        .def_readwrite("x", &AtomSpec::x)
        .def_readwrite("y", &AtomSpec::y)
        .def_readwrite("z", &AtomSpec::z);

    py::class_<BondLinkSpec>(m, "BondLinkSpec")
        .def(py::init<>())
        .def_readwrite("parentAtomGlobalIndex", &BondLinkSpec::parentAtomGlobalIndex)
        .def_readwrite("childAtomGlobalIndex", &BondLinkSpec::childAtomGlobalIndex)
        .def_readwrite("bondGlobalIndex", &BondLinkSpec::bondGlobalIndex)
        .def_readwrite("moleculeIndex", &BondLinkSpec::moleculeIndex)
        .def_readwrite("ringClosing", &BondLinkSpec::ringClosing)
        .def_readwrite("forceK", &BondLinkSpec::forceK)
        .def_readwrite("forceEquil", &BondLinkSpec::forceEquil);

    py::class_<Atom>(m, "Atom")
        .def(py::init<const AtomSpec&>(), py::arg("spec"));

    py::class_<BondLink>(m, "BondLink")
        .def(py::init<const BondLinkSpec&>(), py::arg("spec"));

    py::class_<BondAngle>(m, "BondAngle")
        .def(py::init<>())
        .def(py::init<int, int, int, SimTK::Real, SimTK::Real>(), 
             py::arg("firstGlobalIndex"), py::arg("secondGlobalIndex"), py::arg("thirdGlobalIndex"), py::arg("k"), py::arg("equil"));

    py::class_<BondTorsion>(m, "BondTorsion")
        .def(py::init<>())
        .def(py::init<int, int, int, int, bool, const std::array<SimTK::Real, 4>&, const std::array<SimTK::Real, 4>&, const std::array<int, 4>&>(), 
             py::arg("firstGlobalIndex"), py::arg("secondGlobalIndex"), py::arg("thirdGlobalIndex"), py::arg("fourthGlobalIndex"), py::arg("improper"), py::arg("k"), py::arg("phase"), py::arg("period"));

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
