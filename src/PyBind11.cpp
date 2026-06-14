#include <Python.h>

#include <optional>
#include <pybind11/cast.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <vector>

#include "molmodel/internal/Compound.h"

#include "BondCenter.hpp"
#include "CompoundSystem.h"
#include "Context.hpp"
#include "Rotation.h"
#include "SmallMatrix.h"
#include "TopologyElements.hpp"
#include "World.hpp"

namespace py = pybind11;

using VectorInt = std::vector<int>;
PYBIND11_MAKE_OPAQUE(VectorInt);
PYBIND11_MAKE_OPAQUE(std::vector<RoboAtom>);
PYBIND11_MAKE_OPAQUE(std::vector<RoboBond>);
PYBIND11_MAKE_OPAQUE(std::vector<RoboAngle>);
PYBIND11_MAKE_OPAQUE(std::vector<RoboPeriodicTorsion>);
PYBIND11_MAKE_OPAQUE(std::vector<RoboHarmonicImproperTorsion>);
PYBIND11_MAKE_OPAQUE(std::vector<CMAPGrid>);
PYBIND11_MAKE_OPAQUE(std::vector<CMAPTorsion>);
PYBIND11_MAKE_OPAQUE(std::vector<UreyBradley>);
PYBIND11_MAKE_OPAQUE(std::vector<Scaling14>);
PYBIND11_MAKE_OPAQUE(std::vector<Exclusion>);
PYBIND11_MAKE_OPAQUE(std::vector<TopologyRange>);
PYBIND11_MAKE_OPAQUE(std::vector<SimTK::RootMobility>);

auto transform_to_numpy(const SimTK::Transform& transform) -> py::array_t<SimTK::Real> {
    py::array_t<SimTK::Real> array({3, 4});
    auto buffer = array.mutable_unchecked<2>();

    const auto& rotation = transform.R();
    const auto& position = transform.p();

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            buffer(i, j) = rotation[i][j];
        }
        buffer(i, 3) = position[i];
    }

    return array;
}

auto numpy_to_transform(const py::array_t<SimTK::Real>& array) -> SimTK::Transform {
    if (array.ndim() != 2 || array.shape(0) != 3 || array.shape(1) != 4) {
        throw std::runtime_error("Expected a (3,4) array for Transform");
    }

    auto buffer = array.unchecked<2>();

    SimTK::Mat33 rotation;
    SimTK::Vec3 position;

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            rotation[i][j] = buffer(i, j);
        }
        position[i] = buffer(i, 3);
    }

    return {SimTK::Rotation(rotation), position};
}

namespace pybind11::detail {

template <>
struct type_caster<SimTK::UnitVec3> {
    public:
    /**
     * Modernized Name:
     * By using a list-like type name, stub-gen knows what to expect.
     * If you want it to show up as 'UnitVec3' in stubs, you MUST
     * also export a dummy class or type alias in your main module.
     */
    PYBIND11_TYPE_CASTER(SimTK::UnitVec3, _("Annotated[list[float], FixedSize(3)]"));

    /**
     * Python -> C++ (Standardized)
     */
    bool load(handle src, bool convert) {
        if (!src) {
            return false;
        }

        // Try to cast to a sequence (works for list, tuple, or numpy array)
        if (!py::isinstance<py::sequence>(src)) {
            return false;
        }

        auto seq = py::reinterpret_borrow<py::sequence>(src);
        if (seq.size() != 3) {
            return false;
        }

        try {
            SimTK::Vec3 raw_vec;
            for (size_t i = 0; i < 3; ++i) {
                raw_vec[i] = seq[i].cast<SimTK::Real>();
            }

            // SimTK::UnitVec3 handles normalization and validation
            value = SimTK::UnitVec3(raw_vec);
            return true;
        } catch (...) {
            return false;
        }
    }

    /**
     * C++ -> Python
     * Returns a tuple (more "modern" for fixed-size mathematical vectors
     * as they are immutable, matching the spirit of UnitVec3).
     */
    static auto cast(const SimTK::UnitVec3& src, return_value_policy /* policy */, handle /* parent */)
        -> handle {
        py::tuple vector(3);
        vector[0] = py::cast(src[0]);
        vector[1] = py::cast(src[1]);
        vector[2] = py::cast(src[2]);
        return vector.release();
    }
};

} // namespace pybind11::detail

PYBIND11_MODULE(MODULE_NAME, m) {
    m.doc() = "Robosample bindings";

    py::bind_vector<VectorInt>(m, "VectorInt").def("__repr__", [](const VectorInt&) -> std::string {
        return "VectorInt()";
    });
    py::implicitly_convertible<py::iterable, VectorInt>();

    py::class_<SimTK::ReferenceIndices>(m, "ReferenceIndices")
        .def(py::init<int, int, int>(), py::arg("zero"), py::arg("one"), py::arg("two"))
        .def_readwrite("zero", &SimTK::ReferenceIndices::zero)
        .def_readwrite("one", &SimTK::ReferenceIndices::one)
        .def_readwrite("two", &SimTK::ReferenceIndices::two);

    py::enum_<SimTK::BondCenter::Chirality>(m, "BondCenterChirality")
        .value("RightHanded", SimTK::BondCenter::Chirality::RightHanded)
        .value("LeftHanded", SimTK::BondCenter::Chirality::LeftHanded)
        .value("Planar", SimTK::BondCenter::Chirality::Planar);

    // m.def("triple_product", &SimTK::tripleProduct, "");
    // m.def("plane_normal", &SimTK::planeNormal, "");
    // m.def("is_chirality_mismatch", &SimTK::isChiralityMismatch, "");
    // m.def("signed_plane_deviation", &SimTK::signedPlaneDeviation, "");
    // m.def("exceeds_planarity_threshold", &SimTK::exceedsPlanarityThreshold, "");
    // m.def("chirality_from_plane_deviation", &SimTK::chiralityFromPlaneDeviation, "");
    // m.def("flipped_chirality", &SimTK::flippedChirality, "");
    // m.def("resolve_reference_indices", &SimTK::resolveReferenceIndices, "");
    // m.def("is_bond_chirality_mismatch", &SimTK::isBondChiralityMismatch, "");

    m.def("calculate_log_sum_exp2", &calculateLogSumExp2, "");
    m.def("calculate_angle_in_rad",
          [](const py::array_t<SimTK::Real>& pos0,
             const py::array_t<SimTK::Real>& pos1,
             const py::array_t<SimTK::Real>& pos2) -> SimTK::Real {
              SimTK::Vec3 vec0(pos0.at(0), pos0.at(1), pos0.at(2));
              SimTK::Vec3 vec1(pos1.at(0), pos1.at(1), pos1.at(2));
              SimTK::Vec3 vec2(pos2.at(0), pos2.at(1), pos2.at(2));
              return calculateAngleInRad(vec0, vec1, vec2);
          });
    m.def("calculate_dihedral_in_rad",
          [](const py::array_t<SimTK::Real>& pos0,
             const py::array_t<SimTK::Real>& pos1,
             const py::array_t<SimTK::Real>& pos2,
             const py::array_t<SimTK::Real>& pos3) -> SimTK::Real {
              SimTK::Vec3 vec0(pos0.at(0), pos0.at(1), pos0.at(2));
              SimTK::Vec3 vec1(pos1.at(0), pos1.at(1), pos1.at(2));
              SimTK::Vec3 vec2(pos2.at(0), pos2.at(1), pos2.at(2));
              SimTK::Vec3 vec3(pos3.at(0), pos3.at(1), pos3.at(2));
              return calculateDihedralInRad(vec0, vec1, vec2, vec3);
          });
    m.def("calculate_mag_sq", &calculateMagSq, "");
    m.def("normalize_in_place", &normalizeInPlace, "");
    m.def("multiply_by_scalar", &multiplyByScalar, "");
    m.def("safe_log_sine_sqr", &safeLogSineSqr, "");

    m.def("align_flip_and_translate_frame_along_x_axis",
          [](const py::array_t<SimTK::Real>& gTransform_F1,
             const py::array_t<SimTK::Real>& gPoint_v1) -> py::array_t<SimTK::Real> {
              auto TX = numpy_to_transform(gTransform_F1);
              SimTK::Vec3 vv(gPoint_v1.at(0), gPoint_v1.at(1), gPoint_v1.at(2));

              auto out = alignFlipAndTranslateFrameAlongXAxis(TX, vv);
              return transform_to_numpy(out);
          });

    py::class_<SimTK::DuMM::AtomClassIndex>(m, "AtomClassIndex")
        .def(py::init<int>())
        .def("__int__",
             [](const SimTK::DuMM::AtomClassIndex& i) {
                 return int(i);
             })
        .def("__repr__", [](const SimTK::DuMM::AtomClassIndex& i) {
            return "AtomClassIndex(" + std::to_string(int(i)) + ")";
        });

    py::class_<SimTK::DuMM::ChargedAtomTypeIndex>(m, "ChargedAtomTypeIndex")
        .def(py::init<int>())
        .def("__int__",
             [](const SimTK::DuMM::ChargedAtomTypeIndex& i) {
                 return int(i);
             })
        .def("__repr__", [](const SimTK::DuMM::ChargedAtomTypeIndex& i) {
            return "ChargedAtomTypeIndex(" + std::to_string(int(i)) + ")";
        });

    py::class_<SimTK::Compound::AtomIndex>(m, "CompoundAtomIndex")
        .def(py::init<int>())
        .def("__int__",
             [](const SimTK::Compound::AtomIndex& i) {
                 return int(i);
             })
        .def("__repr__", [](const SimTK::Compound::AtomIndex& i) {
            return "CompoundAtomIndex(" + std::to_string(int(i)) + ")";
        });

    py::class_<SimTK::Vec3>(m, "Vec3")
        .def(py::init<>())
        .def(py::init<SimTK::Real, SimTK::Real, SimTK::Real>())
        .def(py::init([](std::vector<SimTK::Real> v) {
            if (v.size() != 3) {
                throw py::value_error("Vec3 must have 3 elements");
            }
            return new SimTK::Vec3(v[0], v[1], v[2]);
        }))
        .def("__getitem__",
             [](const SimTK::Vec3& v, int i) {
                 if (i < 0 || i >= 3) {
                     throw py::index_error();
                 }
                 return v[i];
             })
        .def("__setitem__",
             [](SimTK::Vec3& v, int i, SimTK::Real val) {
                 if (i < 0 || i >= 3) {
                     throw py::index_error();
                 }
                 v[i] = val;
             })
        .def("__repr__", [](const SimTK::Vec3& v) {
            return "Vec3(" + std::to_string(v[0]) + ", " + std::to_string(v[1]) + ", " + std::to_string(v[2])
                   + ")";
        });

    py::enum_<NonbondedMethod>(m, "NonbondedMethod")
        .value("NoCutoff", NonbondedMethod::NoCutoff)
        .value("CutoffNonPeriodic", NonbondedMethod::CutoffNonPeriodic)
        .export_values();

    py::enum_<SimTK::RootMobility>(m, "RootMobility")
        .value("Free", SimTK::RootMobility::Free)
        .value("Cartesian", SimTK::RootMobility::Cartesian)
        .value("Weld", SimTK::RootMobility::Weld)
        .value("FreeLine", SimTK::RootMobility::FreeLine)
        .value("Ball", SimTK::RootMobility::Ball)
        .value("Pin", SimTK::RootMobility::Pin);

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
        .value("DEFAULT", RUN_TYPE::Default)
        .value("REMC", RUN_TYPE::REMC)
        .value("RENEMC", RUN_TYPE::RENEMC)
        .value("RENE", RUN_TYPE::RENE);

    py::enum_<SamplerName>(m, "SamplerName")
        .value("EMPTY", SamplerName::Empty)
        .value("MC", SamplerName::MC)
        .value("HMC", SamplerName::HMC)
        .value("LAHMC", SamplerName::LAHMC);

    py::enum_<AcceptRejectMode>(m, "AcceptRejectMode")
        .value("AlwaysAccept", AcceptRejectMode::AlwaysAccept)
        .value("MetropolisHastings", AcceptRejectMode::MetropolisHastings);

    py::enum_<IntegratorType>(m, "IntegratorType")
        .value("EMPTY", IntegratorType::Empty)
        .value("VERLET", IntegratorType::Verlet)
        .value("EULER", IntegratorType::Euler)
        .value("EULER2", IntegratorType::Euler2)
        .value("CPODES", IntegratorType::CPodes)
        .value("RUNGEKUTTA", IntegratorType::RungeKutta)
        .value("RUNGEKUTTA2", IntegratorType::RungeKutta2)
        .value("RUNGEKUTTA3", IntegratorType::RungeKutta3)
        .value("RUNGEKUTTAFELDBERG", IntegratorType::RungeKuttaFeldberg)
        .value("BENDSTRETCH", IntegratorType::BendStretch)
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

    py::class_<RoboAtomPhysics>(m, "RoboAtomPhysics")
        .def(py::init<SimTK::Real,
                      SimTK::Real,
                      SimTK::Real,
                      SimTK::Real,
                      SimTK::Real,
                      SimTK::Real,
                      SimTK::Real>(),
             py::arg("charge_e"),
             py::arg("mass_daltons"),
             py::arg("vdw_radius_nm"),
             py::arg("vdw_well_depth_kj"),
             py::arg("sigma_nm"),
             py::arg("solvent_radius_nm"),
             py::arg("screen"))
        .def_readwrite("charge_e", &RoboAtomPhysics::chargeInE)
        .def_readwrite("mass_daltons", &RoboAtomPhysics::massInDaltons)
        .def_readwrite("vdw_radius_nm", &RoboAtomPhysics::vdwRadiusInNm)
        .def_readwrite("vdw_well_depth_kj", &RoboAtomPhysics::vdwWellDepthInKJ)
        .def_readwrite("sigma_nm", &RoboAtomPhysics::sigmaInNm)
        .def_readwrite("solvent_radius_nm", &RoboAtomPhysics::solventRadiusInNm)
        .def_readwrite("screen", &RoboAtomPhysics::screen);

    py::class_<RoboAtomIdentity>(m, "RoboAtomIdentity")
        .def(py::init<std::string,
                      std::string,
                      std::string,
                      std::string,
                      int,
                      int,
                      int,
                      int,
                      int,
                      int,
                      int,
                      int>(),
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
             py::arg("charged_atom_type_index"))
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
             py::arg("atomic_number"))
        .def_readwrite("elementName", &RoboAtomElement::elementName)
        .def_readwrite("elementSymbol", &RoboAtomElement::elementSymbol)
        .def_readwrite("atomicNumber", &RoboAtomElement::atomicNumber);

    py::class_<RoboAtomConnectivity>(m, "RoboAtomConnectivity")
        .def(py::init<VectorInt, bool>(), py::arg("neighbors_global_indices"), py::arg("root"))
        .def_readwrite("neighbors_global_indices", &RoboAtomConnectivity::neighborsGlobalIndices)
        .def_readwrite("root", &RoboAtomConnectivity::root);

    py::class_<RoboAtom>(m, "RoboAtom")
        .def(py::init<RoboAtomIdentity,
                      RoboAtomElement,
                      RoboAtomPhysics,
                      RoboAtomConnectivity,
                      std::array<SimTK::Real, 3>>(),
             py::arg("identity"),
             py::arg("element_info"),
             py::arg("physics"),
             py::arg("connectivity"),
             py::arg("position"))
        .def_readwrite("identity", &RoboAtom::identity)
        .def_readwrite("element_info", &RoboAtom::elementInfo)
        .def_readwrite("physics", &RoboAtom::physics)
        .def_readwrite("connectivity", &RoboAtom::connectivity)
        .def_readwrite("position", &RoboAtom::position);

    py::class_<RoboBond>(m, "RoboBond")
        .def(py::init<std::array<int, 2>,
                      std::array<int, 2>,
                      std::array<int, 2>,
                      SimTK::Real,
                      SimTK::Real,
                      int,
                      bool,
                      const std::string&>(),
             py::arg("global_indices"),
             py::arg("prmtop_indices"),
             py::arg("compound_atom_indices"),
             py::arg("stiffness_in_kj_per_nm_sq"),
             py::arg("nominal_length_in_nm"),
             py::arg("molecule_index"),
             py::arg("ring_closing"),
             py::arg("dihedral_type"))
        .def_readwrite("global_indices", &RoboBond::globalIndices)
        .def_readwrite("prmtop_indices", &RoboBond::prmtopIndices)
        .def_readwrite("compound_atom_indices", &RoboBond::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboBond::moleculeIndex)
        .def_readwrite("ring_closing", &RoboBond::ringClosing)
        .def_readwrite("stiffness_in_kj_per_nm_sq", &RoboBond::stiffnessInKJPerNmSq)
        .def_readwrite("nominal_length_in_nm", &RoboBond::nominalLengthInNm)
        .def_readwrite("dihedral_type", &RoboBond::dihedralType);

    py::class_<RoboAngle>(m, "RoboAngle")
        .def(py::init<std::array<int, 3>,
                      std::array<int, 3>,
                      std::array<int, 3>,
                      int,
                      SimTK::Real,
                      SimTK::Real>(),
             py::arg("global_indices"),
             py::arg("prmtop_indices"),
             py::arg("compound_atom_indices"),
             py::arg("molecule_index"),
             py::arg("stiffness_in_kj_per_rad_sq"),
             py::arg("nominal_angle_in_deg"))
        .def_readwrite("global_indices", &RoboAngle::globalIndices)
        .def_readwrite("prmtop_indices", &RoboAngle::prmtopIndices)
        .def_readwrite("compound_atom_indices", &RoboAngle::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboAngle::moleculeIndex)
        .def_readwrite("stiffness_in_kj_per_rad_sq", &RoboAngle::stiffnessInKJPerRadSq)
        .def_readwrite("nominal_angle_in_deg", &RoboAngle::nominalAngleInDeg);

    py::class_<RoboPeriodicTorsionTerm>(m, "RoboPeriodicTorsionTerm")
        .def(py::init<SimTK::Real, SimTK::Real, int>(),
             py::arg("amplitude_kj"),
             py::arg("phase_deg"),
             py::arg("periodicity"))
        .def_readwrite("amplitude_kj", &RoboPeriodicTorsionTerm::amplitudeKJ)
        .def_readwrite("phase_deg", &RoboPeriodicTorsionTerm::phaseDeg)
        .def_readwrite("periodicity", &RoboPeriodicTorsionTerm::periodicity);

    py::class_<RoboPeriodicTorsion>(m, "RoboPeriodicTorsion")
        .def(py::init<std::array<int, 4>,
                      std::array<int, 4>,
                      std::array<int, 4>,
                      int,
                      bool,
                      std::vector<RoboPeriodicTorsionTerm>>(),
             py::arg("global_indices"),
             py::arg("prmtop_indices"),
             py::arg("compound_atom_indices"),
             py::arg("molecule_index"),
             py::arg("improper"),
             py::arg("terms"))
        .def_readwrite("terms", &RoboPeriodicTorsion::terms)
        .def_readwrite("global_indices", &RoboPeriodicTorsion::globalIndices)
        .def_readwrite("prmtop_indices", &RoboPeriodicTorsion::prmtopIndices)
        .def_readwrite("compound_atom_indices", &RoboPeriodicTorsion::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboPeriodicTorsion::moleculeIndex)
        .def_readwrite("improper", &RoboPeriodicTorsion::improper);

    py::class_<RoboHarmonicImproperTorsion>(m, "RoboHarmonicImproperTorsion")
        .def(py::init<std::array<int, 4>,
                      std::array<int, 4>,
                      std::array<int, 4>,
                      int,
                      SimTK::Real,
                      SimTK::Real>(),
             py::arg("global_indices"),
             py::arg("prmtop_indices"),
             py::arg("compound_atom_indices"),
             py::arg("molecule_index"),
             py::arg("stiffness_in_kj_per_rad_sq"),
             py::arg("nominal_angle_in_rad"))
        .def_readwrite("global_indices", &RoboHarmonicImproperTorsion::globalIndices)
        .def_readwrite("prmtop_indices", &RoboHarmonicImproperTorsion::prmtopIndices)
        .def_readwrite("compound_atom_indices", &RoboHarmonicImproperTorsion::compoundAtomIndices)
        .def_readwrite("molecule_index", &RoboHarmonicImproperTorsion::moleculeIndex)
        .def_readwrite("stiffness_in_kj_per_rad_sq", &RoboHarmonicImproperTorsion::stiffnessInKJPerRadSq)
        .def_readwrite("nominal_angle_in_rad", &RoboHarmonicImproperTorsion::nominalAngleInRad);

    py::class_<ZMatrixRow>(m, "ZMatrixRow")
        .def(py::init<std::array<int, 4>, std::array<int, 4>, int>(),
             py::arg("global_indices"),
             py::arg("compound_atom_indices"),
             py::arg("molecule_index"))
        .def_readwrite("global_indices", &ZMatrixRow::globalIndices)
        .def_readwrite("compound_atom_indices", &ZMatrixRow::compoundAtomIndices)
        .def_readwrite("molecule_index", &ZMatrixRow::moleculeIndex);

    py::class_<CMAPGrid>(m, "CMAPGrid")
        .def(py::init<>())
        .def_readwrite("size", &CMAPGrid::size)
        .def_readwrite("energy", &CMAPGrid::energy);

    py::class_<CMAPTorsion>(m, "CMAPTorsion")
        .def(py::init<>())
        .def(py::init([](int mapIndex, int a1, int a2, int a3, int a4, int b1, int b2, int b3, int b4)
                          -> CMAPTorsion {
                 return CMAPTorsion{mapIndex, a1, a2, a3, a4, b1, b2, b3, b4};
             }),
             py::arg("mapIndex"),
             py::arg("torsion_a_atom_1_global_index"),
             py::arg("torsion_a_atom_2_global_index"),
             py::arg("torsion_a_atom_3_global_index"),
             py::arg("torsion_a_atom_4_global_index"),
             py::arg("torsion_b_atom_1_global_index"),
             py::arg("torsion_b_atom_2_global_index"),
             py::arg("torsion_b_atom_3_global_index"),
             py::arg("torsion_b_atom_4_global_index"))
        .def_readwrite("mapIndex", &CMAPTorsion::mapIndex)
        .def_readwrite("torsion_a_atom_1_global_index", &CMAPTorsion::torsionAAtom1GlobalIndex)
        .def_readwrite("torsion_a_atom_2_global_index", &CMAPTorsion::torsionAAtom2GlobalIndex)
        .def_readwrite("torsion_a_atom_3_global_index", &CMAPTorsion::torsionAAtom3GlobalIndex)
        .def_readwrite("torsion_a_atom_4_global_index", &CMAPTorsion::torsionAAtom4GlobalIndex)
        .def_readwrite("torsion_b_atom_1_global_index", &CMAPTorsion::torsionBAtom1GlobalIndex)
        .def_readwrite("torsion_b_atom_2_global_index", &CMAPTorsion::torsionBAtom2GlobalIndex)
        .def_readwrite("torsion_b_atom_3_global_index", &CMAPTorsion::torsionBAtom3GlobalIndex)
        .def_readwrite("torsion_b_atom_4_global_index", &CMAPTorsion::torsionBAtom4GlobalIndex);

    py::class_<UreyBradley>(m, "UreyBradley")
        .def(py::init<>())
        .def(py::init([](int a1, int a3, SimTK::Real k, SimTK::Real r0) -> UreyBradley {
                 return UreyBradley{a1, a3, k, r0};
             }),
             py::arg("atom_1_global_index"),
             py::arg("atom_3_global_index"),
             py::arg("stiffness_in_kj_per_nm_sq"),
             py::arg("nominal_length_in_nm"))
        .def_readwrite("atom_1_global_index", &UreyBradley::atom1GlobalIndex)
        .def_readwrite("atom_3_global_index", &UreyBradley::atom3GlobalIndex)
        .def_readwrite("stiffness_in_kj_per_nm_sq", &UreyBradley::stiffnessInKJPerNmSq)
        .def_readwrite("nominal_length_in_nm", &UreyBradley::nominalLengthInNm);

    py::class_<Scaling14>(m, "Scaling14")
        .def(py::init<>())
        .def(py::init([](int a1, int a4, SimTK::Real chargeProduct, SimTK::Real epsilon, SimTK::Real sigma)
                          -> Scaling14 {
                 return Scaling14{a1, a4, chargeProduct, epsilon, sigma};
             }),
             py::arg("atom_1_global_index"),
             py::arg("atom_4_global_index"),
             py::arg("charge_product"),
             py::arg("epsilon"),
             py::arg("sigma"))
        .def_readwrite("atom_1_global_index", &Scaling14::atom1GlobalIndex)
        .def_readwrite("atom_4_global_index", &Scaling14::atom4GlobalIndex)
        .def_readwrite("charge_product", &Scaling14::chargeProduct)
        .def_readwrite("epsilon", &Scaling14::epsilon)
        .def_readwrite("sigma", &Scaling14::sigma);

    py::class_<Exclusion>(m, "Exclusion")
        .def(py::init<>())
        .def(py::init([](int a1, int a2) -> Exclusion {
                 return Exclusion{a1, a2};
             }),
             py::arg("atom_1_global_index"),
             py::arg("atom_2_global_index"))
        .def_readwrite("atom_1_global_index", &Exclusion::atom1GlobalIndex)
        .def_readwrite("atom_2_global_index", &Exclusion::atom2GlobalIndex);

    // System topology contains many vectors, so we bind them all as vector types to get list-like behavior
    py::bind_vector<std::vector<TopologyRange>>(m, "VectorTopologyRange")
        .def("__repr__", [](const std::vector<TopologyRange>&) -> std::string {
            return "VectorTopologyRange()";
        });

    py::bind_vector<std::vector<RoboAtom>>(m, "VectorRoboAtom")
        .def("__repr__", [](const std::vector<RoboAtom>&) -> std::string {
            return "VectorRoboAtom()";
        });

    py::bind_vector<std::vector<RoboBond>>(m, "VectorRoboBond")
        .def("__repr__", [](const std::vector<RoboBond>&) -> std::string {
            return "VectorRoboBond()";
        });

    py::bind_vector<std::vector<RoboAngle>>(m, "VectorRoboAngle")
        .def("__repr__", [](const std::vector<RoboAngle>&) -> std::string {
            return "VectorRoboAngle()";
        });

    py::bind_vector<std::vector<RoboPeriodicTorsion>>(m, "VectorRoboPeriodicTorsion")
        .def("__repr__", [](const std::vector<RoboPeriodicTorsion>&) -> std::string {
            return "VectorRoboPeriodicTorsion()";
        });

    py::bind_vector<std::vector<RoboHarmonicImproperTorsion>>(m, "VectorRoboHarmonicImproperTorsion")
        .def("__repr__", [](const std::vector<RoboHarmonicImproperTorsion>&) -> std::string {
            return "VectorRoboHarmonicImproperTorsion()";
        });

    py::bind_vector<std::vector<CMAPGrid>>(m, "VectorCMAPGrid")
        .def("__repr__", [](const std::vector<CMAPGrid>&) -> std::string {
            return "VectorCMAPGrid()";
        });

    py::bind_vector<std::vector<CMAPTorsion>>(m, "VectorCMAPTorsion")
        .def("__repr__", [](const std::vector<CMAPTorsion>&) -> std::string {
            return "VectorCMAPTorsion()";
        });

    py::bind_vector<std::vector<UreyBradley>>(m, "VectorUreyBradley")
        .def("__repr__", [](const std::vector<UreyBradley>&) -> std::string {
            return "VectorUreyBradley()";
        });

    py::bind_vector<std::vector<Scaling14>>(m, "VectorScaling14")
        .def("__repr__", [](const std::vector<Scaling14>&) -> std::string {
            return "VectorScaling14()";
        });

    py::bind_vector<std::vector<Exclusion>>(m, "VectorExclusion")
        .def("__repr__", [](const std::vector<Exclusion>&) -> std::string {
            return "VectorExclusion()";
        });

    py::bind_vector<std::vector<SimTK::RootMobility>>(m, "VectorRootMobility")
        .def("__repr__", [](const std::vector<SimTK::RootMobility>&) -> std::string {
            return "VectorRootMobility()";
        });

    // Allow construction from plain Python lists
    py::implicitly_convertible<py::iterable, std::vector<TopologyRange>>();
    py::implicitly_convertible<py::iterable, std::vector<RoboAtom>>();
    py::implicitly_convertible<py::iterable, std::vector<RoboBond>>();
    py::implicitly_convertible<py::iterable, std::vector<RoboAngle>>();
    py::implicitly_convertible<py::iterable, std::vector<RoboPeriodicTorsion>>();
    py::implicitly_convertible<py::iterable, std::vector<RoboHarmonicImproperTorsion>>();
    py::implicitly_convertible<py::iterable, std::vector<CMAPGrid>>();
    py::implicitly_convertible<py::iterable, std::vector<CMAPTorsion>>();
    py::implicitly_convertible<py::iterable, std::vector<UreyBradley>>();
    py::implicitly_convertible<py::iterable, std::vector<Scaling14>>();
    py::implicitly_convertible<py::iterable, std::vector<Exclusion>>();
    py::implicitly_convertible<py::iterable, std::vector<SimTK::RootMobility>>();

    py::class_<SystemTopology>(m, "SystemTopology")
        .def(py::init<>())
        .def(py::init(
                 [](const std::optional<VectorInt>& root_indices,
                    const std::optional<std::vector<TopologyRange>>& ranges,
                    const std::optional<std::vector<RoboAtom>>& atoms,
                    const std::optional<std::vector<RoboBond>>& bonds,
                    const std::optional<std::vector<RoboAngle>>& angles,
                    const std::optional<std::vector<RoboPeriodicTorsion>>& propers,
                    const std::optional<std::vector<RoboHarmonicImproperTorsion>>& impropers,
                    const std::optional<std::vector<CMAPGrid>>& grids,
                    const std::optional<std::vector<CMAPTorsion>>& cmap_torsions,
                    const std::optional<std::vector<UreyBradley>>& urey_bradleys,
                    const std::optional<std::vector<Scaling14>>& scaling14s,
                    const std::optional<std::vector<Exclusion>>& exclusions,
                    const std::optional<std::vector<SimTK::RootMobility>>& rootMobilities) -> SystemTopology {
                     return SystemTopology{root_indices.value_or(VectorInt()),
                                           ranges.value_or(std::vector<TopologyRange>()),
                                           atoms.value_or(std::vector<RoboAtom>()),
                                           bonds.value_or(std::vector<RoboBond>()),
                                           angles.value_or(std::vector<RoboAngle>()),
                                           propers.value_or(std::vector<RoboPeriodicTorsion>()),
                                           impropers.value_or(std::vector<RoboHarmonicImproperTorsion>()),
                                           grids.value_or(std::vector<CMAPGrid>()),
                                           cmap_torsions.value_or(std::vector<CMAPTorsion>()),
                                           urey_bradleys.value_or(std::vector<UreyBradley>()),
                                           scaling14s.value_or(std::vector<Scaling14>()),
                                           exclusions.value_or(std::vector<Exclusion>()),
                                           rootMobilities.value_or(std::vector<SimTK::RootMobility>())};
                 }),
             py::arg("root_atom_global_indices") = py::none(),
             py::arg("topology_ranges") = py::none(),
             py::arg("atoms") = py::none(),
             py::arg("bonds") = py::none(),
             py::arg("angles") = py::none(),
             py::arg("periodic_torsions") = py::none(),
             py::arg("harmonic_improper_torsions") = py::none(),
             py::arg("cmap_grids") = py::none(),
             py::arg("cmap_torsions") = py::none(),
             py::arg("urey_bradleys") = py::none(),
             py::arg("scaling14s") = py::none(),
             py::arg("exclusions") = py::none(),
             py::arg("root_mobilities") = py::none())
        .def_readwrite("root_atom_global_indices", &SystemTopology::rootAtomGlobalIndices)
        .def_readwrite("topology_ranges", &SystemTopology::topologyRanges)
        .def_readwrite("atoms", &SystemTopology::atoms)
        .def_readwrite("bonds", &SystemTopology::bonds)
        .def_readwrite("angles", &SystemTopology::angles)
        .def_readwrite("periodic_torsions", &SystemTopology::periodicTorsions)
        .def_readwrite("harmonic_improper_torsions", &SystemTopology::harmonicImproperTorsions)
        .def_readwrite("cmap_grids", &SystemTopology::cmapGrids)
        .def_readwrite("cmap_torsions", &SystemTopology::cmapTorsions)
        .def_readwrite("urey_bradleys", &SystemTopology::ureyBradleys)
        .def_readwrite("scaling14s", &SystemTopology::scaling14s)
        .def_readwrite("exclusions", &SystemTopology::exclusions)
        .def_readwrite("root_mobilities", &SystemTopology::rootMobilities);

    py::class_<ForceFieldParams>(m, "ForceFieldParams")
        .def(py::init<>())
        .def(py::init([](bool has_nbfix,
                         int num_types,
                         std::optional<std::vector<SimTK::Real>> a_coef,
                         std::optional<std::vector<SimTK::Real>> b_coef,
                         bool use_gbsaobc2,
                         SimTK::Real gbsa_solvent_dielectric,
                         SimTK::Real gbsa_solute_dielectric,
                         std::optional<NonbondedMethod> nonbonded_method,
                         SimTK::Real nonbonded_cutoff_in_nm) -> ForceFieldParams {
                 return ForceFieldParams{has_nbfix,
                                         num_types,
                                         a_coef.value_or(std::vector<SimTK::Real>()),
                                         b_coef.value_or(std::vector<SimTK::Real>()),
                                         use_gbsaobc2,
                                         gbsa_solvent_dielectric,
                                         gbsa_solute_dielectric,
                                         nonbonded_method.value_or(NonbondedMethod::NoCutoff),
                                         nonbonded_cutoff_in_nm};
             }),
             py::arg("has_nbfix") = false,
             py::arg("num_types") = 0,
             py::arg("a_coef") = py::none(),
             py::arg("b_coef") = py::none(),
             py::arg("use_gbsaobc2") = true,
             py::arg("gbsa_solvent_dielectric") = 78.5,
             py::arg("gbsa_solute_dielectric") = 1.0,
             py::arg("nonbonded_method") = py::none(),
             py::arg("nonbonded_cutoff_in_nm") = 1.2)
        .def_readwrite("has_nbfix", &ForceFieldParams::hasNBfix)
        .def_readwrite("num_types", &ForceFieldParams::numTypes)
        .def_readwrite("a_coef", &ForceFieldParams::aCoef)
        .def_readwrite("b_coef", &ForceFieldParams::bCoef)
        .def_readwrite("use_gbsaobc2", &ForceFieldParams::useGBSAOBC2)
        .def_readwrite("gbsa_solvent_dielectric", &ForceFieldParams::gbsaSolventDielectric)
        .def_readwrite("gbsa_solute_dielectric", &ForceFieldParams::gbsaSoluteDielectric)
        .def_readwrite("nonbonded_method", &ForceFieldParams::nonbondedMethod)
        .def_readwrite("nonbonded_cutoff_in_nm", &ForceFieldParams::nonbondedCutoffInNm);

    py::class_<SimulationSettings>(m, "SimulationSettings")
        .def(py::init<>())
        .def(py::init([](SimTK::Real temp, SimTK::Real freq, int seed) -> SimulationSettings {
                 return SimulationSettings{temp, freq, seed};
             }),
             py::arg("thermostat_temperature_in_k") = 300.0,
             py::arg("collision_frequency") = 1.0,
             py::arg("seed") = 0)
        .def_readwrite("seed", &SimulationSettings::seed)
        .def_readwrite("thermostat_temperature_in_k", &SimulationSettings::thermostatTemperatureInK)
        .def_readwrite("collision_frequency", &SimulationSettings::collisionFrequency);

    py::enum_<TopologyRangeType>(m, "TopologyRangeType")
        .value("EMPTY", TopologyRangeType::Atom)
        .value("BOND", TopologyRangeType::Bond)
        .value("ANGLE", TopologyRangeType::Angle)
        .value("PERIODIC_TORSION", TopologyRangeType::PeriodicTorsion)
        .value("IMPROPER_HARMONIC_TORSION", TopologyRangeType::ImproperHarmonicTorsion);

    py::class_<TopologyRange>(m, "TopologyRange")
        .def(py::init<VectorInt>(), py::arg("startCounts"))
        .def("close", &TopologyRange::close, py::arg("endCounts"));

    py::class_<CoordinateTransferError>(m, "CoordinateTransferError")
        .def_readwrite("matchResiduals", &CoordinateTransferError::matchResiduals)
        .def_readwrite("cartesian", &CoordinateTransferError::cartesian)
        .def_readwrite("cartesianMax", &CoordinateTransferError::cartesianMax)
        .def_readwrite("bonds", &CoordinateTransferError::bonds)
        .def_readwrite("bondsMax", &CoordinateTransferError::bondsMax)
        .def_readwrite("angles", &CoordinateTransferError::angles)
        .def_readwrite("anglesMax", &CoordinateTransferError::anglesMax)
        .def_readwrite("properDihedrals", &CoordinateTransferError::properDihedrals)
        .def_readwrite("properDihedralsMax", &CoordinateTransferError::properDihedralsMax)
        .def_readwrite("improperDihedrals", &CoordinateTransferError::improperDihedrals)
        .def_readwrite("improperDihedralsMax", &CoordinateTransferError::improperDihedralsMax);

    py::class_<Context>(m, "Context")
        .def(py::init<const std::string&, uint32_t, uint32_t, RUN_TYPE, uint32_t, uint32_t, bool>())
        .def("getAtomNameByPrmtopIndex",
             &Context::getAtomNameByPrmtopIndex,
             py::arg("prmtopIndex"),
             "Get the unique atom name for a given prmtop index.")
        .def("addReplica", &Context::addReplica, "Add an empty replica to the context.")
        .def("addThermodynamicState",
             &Context::addThermodynamicState,
             "Add an empty themodynamic state to the context.")
        .def("validate_context",
             &Context::validateContext,
             "Validates all worlds and replicas in the context.")
        .def("run_rex",
             &Context::RunREX,
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
                sphere_radius_in_nm: A real number indicating the sphere radius in nanometers.
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
             py::arg("sphere_radius_in_nm"),
             py::arg("use_fixman_potential"),
             py::arg("use_nuts"),
             "Add a sampler to the world.")
        .def("get_coordinate_transfer_errors",
             &World::getCoordinateTransferErrors,
             "Get the coordinate transfer errors for all samplers in the world.")
        .def("has_rigid_body_violations",
             &World::hasRigidBodyViolations,
             py::arg("timeStep"),
             py::arg("numSteps"),
             "Checks for rigid body violations.");
}
