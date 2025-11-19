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

    // py::class_<AtomClassDefinition>(m, "AtomClassDefinition")
    //     .def(py::init<>())
    //     .def_readwrite("atomTypeName", &AtomClassDefinition::atomTypeName)
    //     .def_readwrite("vdwRadiusInNm", &AtomClassDefinition::vdwRadiusInNm)
    //     .def_readwrite("vdwWellDepthInKJ", &AtomClassDefinition::vdwWellDepthInKJ)
    //     .def_readwrite("atomClassIndex", &AtomClassDefinition::atomClassIndex)
    //     .def_readwrite("atomicNumber", &AtomClassDefinition::atomicNumber)
    //     .def_readwrite("expectedValence", &AtomClassDefinition::expectedValence);

    // py::class_<ChargedAtomTypeDefinition>(m, "ChargedAtomTypeDefinition")
    //     .def(py::init<>())
    //     .def_readwrite("biotypeAtomName", &ChargedAtomTypeDefinition::biotypeAtomName)
    //     .def_readwrite("biotypeResidueName", &ChargedAtomTypeDefinition::biotypeResidueName)
    //     .def_readwrite("partialChargeInE", &ChargedAtomTypeDefinition::partialChargeInE)
    //     .def_readwrite("chargedAtomTypeIndex", &ChargedAtomTypeDefinition::chargedAtomTypeIndex)
    //     .def_readwrite("atomClassIndex", &ChargedAtomTypeDefinition::atomClassIndex);

    // py::class_<BondStretchDefinition>(m, "BondStretchDefinition")
    //     .def(py::init<>())
    //     .def_readwrite("atomClassIndex1", &BondStretchDefinition::atomClassIndex1)
    //     .def_readwrite("atomClassIndex2", &BondStretchDefinition::atomClassIndex2)
    //     .def_readwrite("stiffnessInKJperNmSq", &BondStretchDefinition::stiffnessInKJperNmSq)
    //     .def_readwrite("nominalLengthInNm", &BondStretchDefinition::nominalLengthInNm);

    // py::class_<BondBendDefinition>(m, "BondBendDefinition")
    //     .def(py::init<>())
    //     .def_readwrite("atomClassIndex1", &BondBendDefinition::atomClassIndex1)
    //     .def_readwrite("atomClassIndex1", &BondBendDefinition::atomClassIndex2)
    //     .def_readwrite("atomClassIndex1", &BondBendDefinition::atomClassIndex3)
    //     .def_readwrite("atomClassIndex1", &BondBendDefinition::stiffnessInKJPerRadSq)
    //     .def_readwrite("atomClassIndex1", &BondBendDefinition::nominalAngleInDeg);

    // py::class_<BondTorsionDefinition>(m, "BondTorsionDefinition")
    //     .def(py::init<>())
    //     .def_readwrite("atomClassIndex1", &BondTorsionDefinition::atomClassIndex1)
    //     .def_readwrite("atomClassIndex2", &BondTorsionDefinition::atomClassIndex2)
    //     .def_readwrite("atomClassIndex3", &BondTorsionDefinition::atomClassIndex3)
    //     .def_readwrite("atomClassIndex4", &BondTorsionDefinition::atomClassIndex4)
    //     .def_readwrite("ampInKJ", &BondTorsionDefinition::ampInKJ)
    //     .def_readwrite("phaseInDegrees", &BondTorsionDefinition::phaseInDegrees)
    //     .def_readwrite("periodicity", &BondTorsionDefinition::periodicity)
    //     .def_readwrite("improper", &BondTorsionDefinition::improper);

    py::class_<AtomDefinition>(m, "AtomDefinition")
        .def(py::init<>())
        .def_readwrite("globalIndex", &AtomDefinition::globalIndex)
        .def_readwrite("moleculeIndex", &AtomDefinition::moleculeIndex)
        .def_readwrite("residueIndex", &AtomDefinition::residueIndex)
        .def_readwrite("atomClassIndex", &AtomDefinition::atomClassIndex)
        .def_readwrite("chargedAtomTypeIndex", &AtomDefinition::chargedAtomTypeIndex)
        .def_readwrite("biotypeAtomName", &AtomDefinition::biotypeAtomName)
        .def_readwrite("atomClassName", &AtomDefinition::atomClassName)
        .def_readwrite("chargedAtomName", &AtomDefinition::chargedAtomName)
        .def_readwrite("residueName", &AtomDefinition::residueName)
        .def_readwrite("uniqueAtomName", &AtomDefinition::uniqueAtomName)
        .def_readwrite("neighborsGlobalIndices", &AtomDefinition::neighborsGlobalIndices)
        .def_readwrite("root", &AtomDefinition::root)
        .def_readwrite("atomicNumber", &AtomDefinition::atomicNumber)
        .def_readwrite("chargeInE", &AtomDefinition::chargeInE)
        .def_readwrite("massInDaltons", &AtomDefinition::massInDaltons)
        .def_readwrite("vdwRadiusInNm", &AtomDefinition::vdwRadiusInNm)
        .def_readwrite("sigmaInNm", &AtomDefinition::sigmaInNm)
        .def_readwrite("vdwWellDepthInKJ", &AtomDefinition::vdwWellDepthInKJ)
        .def_readwrite("x", &AtomDefinition::x)
        .def_readwrite("y", &AtomDefinition::y)
        .def_readwrite("z", &AtomDefinition::z);

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
        .def_readwrite("ampInKJ", &BondTorsionDefinition::ampInKJ)
        .def_readwrite("phaseInDegrees", &BondTorsionDefinition::phaseInDegrees)
        .def_readwrite("periodicity", &BondTorsionDefinition::periodicity)
        .def_readwrite("improper", &BondTorsionDefinition::improper);

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
