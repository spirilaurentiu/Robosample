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

    py::enum_<SimTK::BondMobility::Mobility>(m, "BondMobility")
        .value("Free", SimTK::BondMobility::Mobility::Free)
        .value("Torsion", SimTK::BondMobility::Mobility::Torsion)
        .value("Rigid", SimTK::BondMobility::Mobility::Rigid)
        .value("BallF", SimTK::BondMobility::Mobility::BallF)
        .value("BallM", SimTK::BondMobility::Mobility::BallM)
        .value("Cylinder", SimTK::BondMobility::Mobility::Cylinder)
        .value("Translation", SimTK::BondMobility::Mobility::Translation)
        .value("FreeLine", SimTK::BondMobility::Mobility::FreeLine)
        .value("LineOrientationF", SimTK::BondMobility::Mobility::LineOrientationF)
        .value("LineOrientationM", SimTK::BondMobility::Mobility::LineOrientationM)
        .value("UniversalM", SimTK::BondMobility::Mobility::UniversalM)
        .value("Spherical", SimTK::BondMobility::Mobility::Spherical)
        .value("AnglePin", SimTK::BondMobility::Mobility::AnglePin)
        .value("BendStretch", SimTK::BondMobility::Mobility::BendStretch)
        .value("Slider", SimTK::BondMobility::Mobility::Slider)
        .value("OrthoSpherical", SimTK::BondMobility::Mobility::OrthoSpherical);

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
        .def(py::init<int, int, SimTK::BondMobility::Mobility>())
        .def_readwrite("i", &BOND_FLEXIBILITY::i)
        .def_readwrite("j", &BOND_FLEXIBILITY::j)
        .def_readwrite("mobility", &BOND_FLEXIBILITY::mobility);

    py::class_<RoboAtomDefinition>(m, "RoboAtomDefinition")
        .def(py::init<>())
        .def_readwrite("global_index", &RoboAtomDefinition::globalIndex)
        .def_readwrite("prmtop_index", &RoboAtomDefinition::prmtopIndex)
        .def_readwrite("compound_atom_index", &RoboAtomDefinition::compoundAtomIndex)
        .def_readwrite("molecule_index", &RoboAtomDefinition::moleculeIndex)
        .def_readwrite("residue_index", &RoboAtomDefinition::residueIndex)
        .def_readwrite("atom_class_name", &RoboAtomDefinition::atomClassName)
        .def_readwrite("atom_class_index", &RoboAtomDefinition::atomClassIndex)
        .def_readwrite("charged_atom_type_name", &RoboAtomDefinition::chargedAtomTypeName)
        .def_readwrite("charged_atom_type_index", &RoboAtomDefinition::chargedAtomTypeIndex)
        .def_readwrite("residue_name", &RoboAtomDefinition::residueName)
        .def_readwrite("unique_atom_name", &RoboAtomDefinition::uniqueAtomName)
        .def_readwrite("neighbors_global_indices", &RoboAtomDefinition::neighborsGlobalIndices)
        .def_readwrite("root", &RoboAtomDefinition::root)
        .def_readwrite("atomic_number", &RoboAtomDefinition::atomicNumber)
        .def_readwrite("charge_in_e", &RoboAtomDefinition::chargeInE)
        .def_readwrite("mass_in_daltons", &RoboAtomDefinition::massInDaltons)
        .def_readwrite("vdw_radius_nm", &RoboAtomDefinition::vdwRadiusInNm)
        .def_readwrite("sigma_nm", &RoboAtomDefinition::sigmaInNm)
        .def_readwrite("vdw_well_depth_kj", &RoboAtomDefinition::vdwWellDepthInKJ)
        .def_readwrite("solventRadiusInNm", &RoboAtomDefinition::solventRadiusInNm)
        .def_readwrite("screen", &RoboAtomDefinition::screen)
        .def_readwrite("x_nm", &RoboAtomDefinition::x_nm)
        .def_readwrite("y_nm", &RoboAtomDefinition::y_nm)
        .def_readwrite("z_nm", &RoboAtomDefinition::z_nm);

    py::class_<RoboBondStretchDefinition>(m, "RoboBondStretchDefinition")
        .def(py::init<>())
        .def_readwrite("parentAtomGlobalIndex", &RoboBondStretchDefinition::parentAtomGlobalIndex)
        .def_readwrite("childAtomGlobalIndex", &RoboBondStretchDefinition::childAtomGlobalIndex)
        .def_readwrite("parentCompoundAtomIndex", &RoboBondStretchDefinition::parentCompoundAtomIndex)
        .def_readwrite("childCompoundAtomIndex", &RoboBondStretchDefinition::childCompoundAtomIndex)
        .def_readwrite("bondGlobalIndex", &RoboBondStretchDefinition::bondGlobalIndex)
        .def_readwrite("moleculeIndex", &RoboBondStretchDefinition::moleculeIndex)
        .def_readwrite("ringClosing", &RoboBondStretchDefinition::ringClosing)
        .def_readwrite("stiffnessInKJPerNmSq", &RoboBondStretchDefinition::stiffnessInKJPerNmSq)
        .def_readwrite("nominalLengthInNm", &RoboBondStretchDefinition::nominalLengthInNm);

    py::class_<RoboBondBendDefinition>(m, "RoboBondBendDefinition")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &RoboBondBendDefinition::globalIndex1)
        .def_readwrite("globalIndex2", &RoboBondBendDefinition::globalIndex2)
        .def_readwrite("globalIndex3", &RoboBondBendDefinition::globalIndex3)
        .def_readwrite("compoundAtomIndex1", &RoboBondBendDefinition::compoundAtomIndex1)
        .def_readwrite("compoundAtomIndex2", &RoboBondBendDefinition::compoundAtomIndex2)
        .def_readwrite("compoundAtomIndex3", &RoboBondBendDefinition::compoundAtomIndex3)
        .def_readwrite("moleculeIndex", &RoboBondBendDefinition::moleculeIndex)
        .def_readwrite("stiffnessInKJPerRadSq", &RoboBondBendDefinition::stiffnessInKJPerRadSq)
        .def_readwrite("nominalAngleInDeg", &RoboBondBendDefinition::nominalAngleInDeg);

    py::class_<RoboBondTorsionDefinition>(m, "RoboBondTorsionDefinition")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &RoboBondTorsionDefinition::globalIndex1)
        .def_readwrite("globalIndex2", &RoboBondTorsionDefinition::globalIndex2)
        .def_readwrite("globalIndex3", &RoboBondTorsionDefinition::globalIndex3)
        .def_readwrite("globalIndex4", &RoboBondTorsionDefinition::globalIndex4)
        .def_readwrite("compoundAtomIndex1", &RoboBondTorsionDefinition::compoundAtomIndex1)
        .def_readwrite("compoundAtomIndex2", &RoboBondTorsionDefinition::compoundAtomIndex2)
        .def_readwrite("compoundAtomIndex3", &RoboBondTorsionDefinition::compoundAtomIndex3)
        .def_readwrite("compoundAtomIndex4", &RoboBondTorsionDefinition::compoundAtomIndex4)
        .def_readwrite("moleculeIndex", &RoboBondTorsionDefinition::moleculeIndex)
        .def_readwrite("improper", &RoboBondTorsionDefinition::improper)
        .def_readwrite("ampInKJ_1", &RoboBondTorsionDefinition::ampInKJ_1)
        .def_readwrite("phaseInDegrees_1", &RoboBondTorsionDefinition::phaseInDegrees_1)
        .def_readwrite("periodicity_1", &RoboBondTorsionDefinition::periodicity_1)
        .def_readwrite("ampInKJ_2", &RoboBondTorsionDefinition::ampInKJ_2)
        .def_readwrite("phaseInDegrees_2", &RoboBondTorsionDefinition::phaseInDegrees_2)
        .def_readwrite("periodicity_2", &RoboBondTorsionDefinition::periodicity_2)
        .def_readwrite("ampInKJ_3", &RoboBondTorsionDefinition::ampInKJ_3)
        .def_readwrite("phaseInDegrees_3", &RoboBondTorsionDefinition::phaseInDegrees_3)
        .def_readwrite("periodicity_3", &RoboBondTorsionDefinition::periodicity_3)
        .def_readwrite("ampInKJ_4", &RoboBondTorsionDefinition::ampInKJ_4)
        .def_readwrite("phaseInDegrees_4", &RoboBondTorsionDefinition::phaseInDegrees_4)
        .def_readwrite("periodicity_4", &RoboBondTorsionDefinition::periodicity_4)
        .def_readwrite("ampInKJ_5", &RoboBondTorsionDefinition::ampInKJ_5)
        .def_readwrite("phaseInDegrees_5", &RoboBondTorsionDefinition::phaseInDegrees_5)
        .def_readwrite("periodicity_5", &RoboBondTorsionDefinition::periodicity_5);

    py::class_<OpenMMEnergyComponents>(m, "OpenMMEnergyComponents")
        .def(py::init<>())
        .def_readwrite("totalEnergy", &OpenMMEnergyComponents::totalEnergy)
        .def_readwrite("harmonicBondForce", &OpenMMEnergyComponents::harmonicBondForce)
        .def_readwrite("harmonicAngleForce", &OpenMMEnergyComponents::harmonicAngleForce)
        .def_readwrite("periodicTorsionForce", &OpenMMEnergyComponents::periodicTorsionForce)
        .def_readwrite("nonbondedForce", &OpenMMEnergyComponents::nonbondedForce)
        .def_readwrite("andersenThermostat", &OpenMMEnergyComponents::andersenThermostat)
        .def_readwrite("gbsaObcForce", &OpenMMEnergyComponents::gbsaObcForce);

    py::class_<RoboAtom>(m, "RoboAtom")
        .def(py::init<const RoboAtomDefinition&>(), py::arg("spec"));

    py::class_<RoboBondStretch>(m, "RoboBondStretch")
        .def(py::init<const RoboBondStretchDefinition&>(), py::arg("spec"));
        
    py::class_<RoboBondBend>(m, "RoboBondBend")
        .def(py::init<const RoboBondBendDefinition&>(), py::arg("spec"));

    py::class_<RoboBondTorsion>(m, "RoboBondTorsion")
        .def(py::init<const RoboBondTorsionDefinition&>(), py::arg("spec"));

    py::class_<IteratorPair>(m, "IteratorPair")
        .def(py::init<>())
        .def_readwrite("begin", &IteratorPair::begin)
        .def_readwrite("last", &IteratorPair::last);

    py::class_<Context>(m, "Context")
        .def(py::init<const std::string&, uint32_t, uint32_t, uint32_t, RUN_TYPE, uint32_t, uint32_t, bool>())
        .def("addReplica", &Context::addReplica, "Add an empty replica to the context.")
        .def("addThermodynamicState", &Context::addThermodynamicState, "Add an empty themodynamic state to the context.")
        .def("Initialize", py::overload_cast<>(&Context::Initialize), "Initializes the context after all worlds and replicas have been set.")
        .def("RunREX", &Context::RunREX, "Run replica exchange.")
        .def("setVerbose", &Context::setVerbose, "Control if you want extraneous output to cout.")
        .def("setPdbRestartFreq", &Context::setPdbRestartFreq, "Set the PDB restart frequency.")
        .def("setPrintFreq", &Context::setPrintFreq, "Set the print frequency.")
        .def("setNonbonded", &Context::setNonbonded, "Set nonbonded method and cutoff.")
        .def("setGBSA", &Context::setGBSA, "Set GBSA.")
        .def("loadAmberSystem", &Context::loadAmberSystem, "Load an AMBER system.")
        .def("addWorld", &Context::addWorld, "Add an empty world.")

        // Binds the const version of getWorld
        // For non-const version: def("getWorld", py::overload_cast<std::size_t>(&Context::getWorld), py::return_value_policy::reference)
        .def("getWorld", py::overload_cast<std::size_t>(&Context::getWorld, py::const_), py::return_value_policy::reference)
        .def("getWorlds", py::overload_cast<>(&Context::getWorlds, py::const_), py::return_value_policy::reference)
        .def("getInitialOpenMMEnergyComponents", &Context::getInitialOpenMMEnergyComponents, "Get energy components from OpenMM.");

    py::class_<World>(m, "World")
        .def("setFlexibilities", &World::setFlexibilities, "Set the flexibilities of the bonds.")
        .def("addSampler", &World::addSampler, "Add a sampler to the world.")
        .def("setRollFlexibilities", &World::setRollFlexibilities, "Set rolling sequence.")
        .def("getRollFlexibilities", &World::getRollFlexibilities, "Get rolling sequence.")
        .def("getMatchAtomTargetLocationsResiduals", &World::getMatchAtomTargetLocationsResiduals, "Get residuals from matchAtomTargetLocations during testing.")
        .def("getCumulativeCartesianDisplacements", &World::getCumulativeCartesianDisplacements, "Get cumulative Cartesian displacements.")
        .def("getCumulativeBondDisplacements", &World::getCumulativeBondDisplacements, "Get cumulative bond displacements.")
        .def("getCumulativeAngleDisplacements", &World::getCumulativeAngleDisplacements, "Get cumulative angle displacements.")
        .def("getCumulativeTorsionDisplacements", &World::getCumulativeTorsionDisplacements, "Get cumulative torsion displacements.")
        .def("getAcceptanceRMSD", &World::getAcceptanceRMSD, "Get RMSD values for accepted/rejected moves during testing.")
        .def("getRigidBodyBondRMSDInNm", &World::getRigidBodyBondRMSDInNm, "Get rigid body bond RMSD in nm during testing.")
        .def("getRigidBodyAngleDriftInRad", &World::getRigidBodyAngleDriftInRad, "Get rigid body angle drift in rad during testing.")
        .def("getRigidBodyProperTorsionDriftInRad", &World::getRigidBodyProperTorsionDriftInRad, "Get rigid body proper torsion drift in rad during testing.")
        .def("getRigidBodyImproperTorsionDriftInRad", &World::getRigidBodyImproperTorsionDriftInRad, "Get rigid body improper torsion drift in rad during testing.");
}
