#include <Python.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "Context.hpp"

namespace py = pybind11;

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc() = "Robosample bindings";

    py::enum_<NonbondedMethod>(m, "NonbondedMethod")
        .value("NoCutoff", NonbondedMethod::NoCutoff)
        .value("CutoffNonPeriodic", NonbondedMethod::CutoffNonPeriodic);

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

    py::class_<BondFlexibility>(m, "BondFlexibility")
        .def(py::init<>())
        .def_readwrite("globalIndex1", &BondFlexibility::globalIndex1)
        .def_readwrite("globalIndex2", &BondFlexibility::globalIndex2)
        .def_readwrite("uniqueAtomName1", &BondFlexibility::uniqueAtomName1)
        .def_readwrite("uniqueAtomName2", &BondFlexibility::uniqueAtomName2)
        .def_readwrite("mobility", &BondFlexibility::mobility);

    py::class_<RoboAtomPhysics>(m, "RoboAtomPhysics")
        .def(py::init<SimTK::Real, SimTK::Real, SimTK::Real, SimTK::Real, SimTK::Real, SimTK::Real, SimTK::Real>(),
            py::arg("charge_e"),
            py::arg("mass_daltons"),
            py::arg("vdw_radius_nm"),
            py::arg("vdw_well_depth_kj"),
            py::arg("sigma_nm"),
            py::arg("solvent_radius_nm"),
            py::arg("screen")
        )
        .def_readwrite("charge_e", &RoboAtomPhysics::chargeInE)
        .def_readwrite("mass_daltons", &RoboAtomPhysics::massInDaltons)
        .def_readwrite("vdw_radius_nm", &RoboAtomPhysics::vdwRadiusInNm)
        .def_readwrite("vdw_well_depth_kj", &RoboAtomPhysics::vdwWellDepthInKJ)
        .def_readwrite("sigma_nm", &RoboAtomPhysics::sigmaInNm)
        .def_readwrite("solvent_radius_nm", &RoboAtomPhysics::solventRadiusInNm)
        .def_readwrite("screen", &RoboAtomPhysics::screen);

    py::class_<RoboAtomIdentity>(m, "RoboAtomIdentity")
        .def(py::init<std::string, std::string, std::string, std::string, int, int, int, int, int, int, int, int>(),
            py::arg("unique_name"),
            py::arg("residue_name"),
            py::arg("atom_class_name"),
            py::arg("charged_atom_type_name"),
            py::arg("global_index"),
            py::arg("prmtop_index"),
            py::arg("molecule_index"),
            py::arg("residue_index"),
            py::arg("nonbonded_index"),
            py::arg("compound_atom_index"),
            py::arg("atom_class_index"),
            py::arg("charged_atom_type_index")
        )
        .def_readwrite("unique_name", &RoboAtomIdentity::uniqueAtomName)
        .def_readwrite("residue_name", &RoboAtomIdentity::residueName)
        .def_readwrite("atom_class_name", &RoboAtomIdentity::atomClassName)
        .def_readwrite("charged_atom_type_name", &RoboAtomIdentity::chargedAtomTypeName)
        .def_readwrite("global_index", &RoboAtomIdentity::globalIndex)
        .def_readwrite("prmtop_index", &RoboAtomIdentity::prmtopIndex)
        .def_readwrite("molecule_index", &RoboAtomIdentity::moleculeIndex)
        .def_readwrite("residue_index", &RoboAtomIdentity::residueIndex)
        .def_readwrite("nonbonded_index", &RoboAtomIdentity::nonbondedIndex)
        .def_readwrite("compound_atom_index", &RoboAtomIdentity::compoundAtomIndex)
        .def_readwrite("atom_class_index", &RoboAtomIdentity::atomClassIndex)
        .def_readwrite("charged_atom_type_index", &RoboAtomIdentity::chargedAtomTypeIndex);

    py::class_<RoboAtomElement>(m, "RoboAtomElement")
        .def(py::init<std::string, std::string, int>(),
            py::arg("element_name"),
            py::arg("element_symbol"),
            py::arg("atomic_number")
        )
        .def_readwrite("elementName", &RoboAtomElement::elementName)
        .def_readwrite("elementSymbol", &RoboAtomElement::elementSymbol)
        .def_readwrite("atomicNumber", &RoboAtomElement::atomicNumber);

    py::class_<RoboAtomConnectivity>(m, "RoboAtomConnectivity")
        .def(py::init<std::vector<int>, bool>(),
            py::arg("neighbors_global_indices"),
            py::arg("root")
        )
        .def_readwrite("neighbors_global_indices", &RoboAtomConnectivity::neighborsGlobalIndices)
        .def_readwrite("root", &RoboAtomConnectivity::root);

    py::class_<RoboAtom>(m, "RoboAtom")
        .def(py::init<RoboAtomIdentity, RoboAtomElement, RoboAtomPhysics, RoboAtomConnectivity, std::array<SimTK::Real, 3>>(),
            py::arg("identity"),
            py::arg("element_info"),
            py::arg("physics"),
            py::arg("connectivity"),
            py::arg("position")
        )
        .def_readwrite("identity", &RoboAtom::identity)
        .def_readwrite("element_info", &RoboAtom::elementInfo)
        .def_readwrite("physics", &RoboAtom::physics)
        .def_readwrite("connectivity", &RoboAtom::connectivity)
        .def_readwrite("position", &RoboAtom::position);

    py::class_<RoboBond>(m, "RoboBond")
        .def(py::init<std::array<int, 2>, std::array<int, 2>, SimTK::Real, SimTK::Real, int, bool>(),
            py::arg("global_indices"),
            py::arg("compound_atom_indices"),
            py::arg("stiffness_in_kj_per_nm_sq"),
            py::arg("nominal_length_in_nm"),
            py::arg("molecule_index"),
            py::arg("ring_closing")
        )
        .def_readwrite("global_indices", &RoboBond::globalIndices)
        .def_readwrite("compound_atom_indices", &RoboBond::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboBond::moleculeIndex)
        .def_readwrite("ring_closing", &RoboBond::ringClosing)
        .def_readwrite("stiffness_in_kj_per_nm_sq", &RoboBond::stiffnessInKJPerNmSq)
        .def_readwrite("nominal_length_in_nm", &RoboBond::nominalLengthInNm)
        ;

    py::class_<RoboAngle>(m, "RoboAngle")
        .def(py::init<std::array<int, 3>, std::array<int, 3>, int, SimTK::Real, SimTK::Real>(),
            py::arg("global_indices"),
            py::arg("compound_atom_indices"),
            py::arg("molecule_index"),
            py::arg("stiffness_in_kj_per_rad_sq"),
            py::arg("nominal_angle_in_deg")
        )
        .def_readwrite("global_indices", &RoboAngle::globalIndices)
        .def_readwrite("compound_atom_indices", &RoboAngle::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboAngle::moleculeIndex)
        .def_readwrite("stiffness_in_kj_per_rad_sq", &RoboAngle::stiffnessInKJPerRadSq)
        .def_readwrite("nominal_angle_in_deg", &RoboAngle::nominalAngleInDeg);

    py::class_<RoboPeriodicTorsionTerm>(m, "RoboPeriodicTorsionTerm")
        .def(py::init<SimTK::Real, SimTK::Real, int>(),
            py::arg("amplitude_kj"),
            py::arg("phase_deg"),
            py::arg("periodicity")
        )
        .def_readwrite("amplitude_kj", &RoboPeriodicTorsionTerm::amplitudeKJ)
        .def_readwrite("phase_deg", &RoboPeriodicTorsionTerm::phaseDeg)
        .def_readwrite("periodicity", &RoboPeriodicTorsionTerm::periodicity);

    py::class_<RoboPeriodicTorsion>(m, "RoboPeriodicTorsion")
        .def(
            py::init< std::array<int, 4>, std::array<int, 4>, int, bool, std::vector<RoboPeriodicTorsionTerm>>(),
            py::arg("global_indices"),
            py::arg("compound_atom_indices"),
            py::arg("molecule_index"),
            py::arg("improper"),
            py::arg("terms")    
        );

    py::class_<RoboHarmonicImproperTorsion>(m, "RoboHarmonicImproperTorsion")
        .def(
            py::init< std::array<int, 4>, std::array<int, 4>, int, SimTK::Real, SimTK::Real>(),
            py::arg("global_indices"),
            py::arg("compound_atom_indices"),
            py::arg("molecule_index"),
            py::arg("stiffness_in_kj_per_rad_sq"),
            py::arg("nominal_angle_in_rad")    
        );

    py::class_<CMAPGrid>(m, "CMAPGrid")
        .def(py::init<>())
        .def_readwrite("size", &CMAPGrid::size)
        .def_readwrite("energy", &CMAPGrid::energy);

    py::class_<CMAPTorsion>(m, "CMAPTorsion")
        .def(py::init<>())
        .def_readwrite("mapIndex", &CMAPTorsion::mapIndex)
        .def_readwrite("a1", &CMAPTorsion::a1)
        .def_readwrite("a2", &CMAPTorsion::a2)
        .def_readwrite("a3", &CMAPTorsion::a3)
        .def_readwrite("a4", &CMAPTorsion::a4)
        .def_readwrite("b1", &CMAPTorsion::b1)
        .def_readwrite("b2", &CMAPTorsion::b2)
        .def_readwrite("b3", &CMAPTorsion::b3)
        .def_readwrite("b4", &CMAPTorsion::b4);

    py::class_<UreyBradley>(m, "UreyBradley")
        .def(py::init<>())
        .def_readwrite("a1", &UreyBradley::a1)
        .def_readwrite("a3", &UreyBradley::a3)
        .def_readwrite("stiffness_in_kj_per_nm_sq", &UreyBradley::stiffnessInKJPerNmSq)
        .def_readwrite("nominal_length_in_nm", &UreyBradley::nominalLengthInNm);

    py::class_<Scaling14>(m, "Scaling14")
        .def(py::init<>())
        .def(py::init<int, int, SimTK::Real, SimTK::Real, SimTK::Real>(),
            py::arg("a1"),
            py::arg("a4"),
            py::arg("charge_product"),
            py::arg("epsilon"),
            py::arg("sigma")
        )
        .def_readwrite("a1", &Scaling14::a1)
        .def_readwrite("a4", &Scaling14::a4)
        .def_readwrite("charge_product", &Scaling14::chargeProduct)
        .def_readwrite("epsilon", &Scaling14::epsilon)
        .def_readwrite("sigma", &Scaling14::sigma);

    py::class_<Exclusion>(m, "Exclusion")
        .def(py::init<>())
        .def(py::init<int, int>(),
            py::arg("a1"),
            py::arg("a2")
        )
        .def_readwrite("a1", &Exclusion::a1)
        .def_readwrite("a2", &Exclusion::a2);

    py::enum_<TopologyRangeType>(m, "TopologyRangeType")
        .value("EMPTY", TopologyRangeType::ATOM)
        .value("BOND", TopologyRangeType::BOND)
        .value("ANGLE", TopologyRangeType::ANGLE)
        .value("PERIODIC_TORSION", TopologyRangeType::PERIODIC_TORSION)
        .value("IMPROPER_HARMONIC_TORSION", TopologyRangeType::IMPROPER_HARMONIC_TORSION);

    py::class_<TopologyRange>(m, "TopologyRange")
        .def(py::init<std::vector<int>>(), py::arg("startCounts"))
        .def("close", &TopologyRange::close, py::arg("endCounts"));

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
        .def("setGBSAOptions", &Context::setGBSAOptions, "Set GBSA-OBC2 options.")
        .def("loadAmberSystem", &Context::loadAmberSystem, "Load an AMBER system.")
        .def("initialize_openmm", &Context::initializeOpenMM, "Load an OpenMM system from components.")
        .def("addWorld", &Context::addWorld, "Add an empty world.")

        // Binds the const version of getWorld
        // For non-const version: def("getWorld", py::overload_cast<std::size_t>(&Context::getWorld), py::return_value_policy::reference)
        .def("getWorld", py::overload_cast<std::size_t>(&Context::getWorld, py::const_), py::return_value_policy::reference)
        .def("getWorlds", py::overload_cast<>(&Context::getWorlds, py::const_), py::return_value_policy::reference);

    py::class_<World>(m, "World")
        .def("addSampler", &World::addSampler, "Add a sampler to the world.")
        .def("getMatchAtomTargetLocationsResiduals", &World::getMatchAtomTargetLocationsResiduals, "Get residuals from matchAtomTargetLocations during testing.")
        .def("getCumulativeCartesianDisplacements", &World::getCumulativeCartesianDisplacements, "Get cumulative Cartesian displacements.")
        .def("getCumulativeBondDisplacements", &World::getCumulativeBondDisplacements, "Get cumulative bond displacements.")
        .def("getCumulativeAngleDisplacements", &World::getCumulativeAngleDisplacements, "Get cumulative angle displacements.")
        .def("getCumulativeTorsionDisplacements", &World::getCumulativeTorsionDisplacements, "Get cumulative torsion displacements.")
        .def("getOpenMMCumulativeCartesianDisplacements", &World::getOpenMMCumulativeCartesianDisplacements, "Get cumulative Cartesian displacements between OpenMM and SimTK during testing.")
        .def("getOpenMMCumulativeBondDisplacements", &World::getOpenMMCumulativeBondDisplacements, "Get cumulative bond displacements between OpenMM and SimTK during testing.")
        .def("getOpenMMCumulativeAngleDisplacements", &World::getOpenMMCumulativeAngleDisplacements, "Get cumulative angle displacements between OpenMM and SimTK during testing.")
        .def("getOpenMMCumulativeTorsionDisplacements", &World::getOpenMMCumulativeTorsionDisplacements, "Get cumulative torsion displacements between OpenMM and SimTK during testing.")
        .def("getAcceptanceRMSD", &World::getAcceptanceRMSD, "Get RMSD values for accepted/rejected moves during testing.")
        .def("getRigidBodyBondRMSDInNm", &World::getRigidBodyBondRMSDInNm, "Get rigid body bond RMSD in nm during testing.")
        .def("getRigidBodyAngleDriftInRad", &World::getRigidBodyAngleDriftInRad, "Get rigid body angle drift in rad during testing.")
        .def("getRigidBodyProperTorsionDriftInRad", &World::getRigidBodyProperTorsionDriftInRad, "Get rigid body proper torsion drift in rad during testing.")
        .def("getRigidBodyImproperTorsionDriftInRad", &World::getRigidBodyImproperTorsionDriftInRad, "Get rigid body improper torsion drift in rad during testing.");
}
