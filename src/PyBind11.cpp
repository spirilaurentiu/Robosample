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
        .value("Slider", BondMobility::Mobility::Slider);

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

    py::enum_<IntegratorName>(m, "IntegratorName")
        .value("Empty", IntegratorName::None)
        .value("Verlet", IntegratorName::Verlet)
        .value("Euler", IntegratorName::Euler)
        .value("Euler2", IntegratorName::Euler2)
        .value("Cpodes", IntegratorName::Cpodes)
        .value("RungeKutta", IntegratorName::RungeKutta)
        .value("RungeKutta2", IntegratorName::RungeKutta2)
        .value("RungeKutta3", IntegratorName::RungeKutta3)
        .value("RungeKuttaFeldberg", IntegratorName::RungeKuttaFeldberg)
        .value("BendStretch", IntegratorName::BendStretch)
        .value("OMMVV", IntegratorName::OMMVV)
        .value("BoundWalk", IntegratorName::BoundWalk)
        .value("BoundHMC", IntegratorName::BoundHMC)
        .value("StationsTask", IntegratorName::StationsTask)
        .value("NofIntegrators", IntegratorName::NofIntegrators);

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

    py::class_<Context>(m, "Context")
        .def(py::init<std::string, uint32_t, uint32_t, uint32_t, RUN_TYPE, uint32_t, uint32_t>())
        .def("addWorld", &Context::addWorld, "Add an empty world")
        .def("getWorld", (World& (Context::*)(std::size_t which)) &Context::getWorld, py::return_value_policy::reference, "Run the simulation")
        .def("loadAmberSystem", &Context::loadAmberSystem, "Load an Amber system")
        .def("RunREXNew", &Context::RunREXNew, "Run the simulation")
        .def("setPdbPrefix", &Context::setPdbPrefix, "Set the prefix for the PDB files")
        .def("setOutput", &Context::setOutput, "Set the output directory")
        .def("setRequiredNofRounds", &Context::setRequiredNofRounds, "Set the required number of rounds")
        .def("setPdbRestartFreq", &Context::setPdbRestartFreq, "Set the PDB restart frequency")
        .def("setPrintFreq", &Context::setPrintFreq, "Set the print frequency")
        .def("setGBSA", &Context::setGBSA, "Set the GBSA global scale factor (0 for vaccum, 1 for water - implicit solvent)")
        .def("setRunType", &Context::setRunType, "Set the run type")
        .def("initialize", &Context::initialize, "Initialize the simulation")
        .def("addReplica", &Context::addReplica, "Add a replica")
        .def("addThermodynamicState", &Context::addThermodynamicState, "Add a thermodynamic state");

    py::class_<World>(m, "World")
        .def("setFlexibilities", &World::setFlexibilities, "Set the flexibilities of the bonds")
        .def("addSampler", &World::addSampler, "Add a sampler to the world");
}
