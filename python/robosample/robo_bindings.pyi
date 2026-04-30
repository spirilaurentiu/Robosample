"""
Robosample bindings
"""
from __future__ import annotations
import numpy
import pybind11_stubgen.typing_ext
import typing
__all__: list[str] = ['AcceptRejectMode', 'AtomClassIndex', 'BondCenterChirality', 'BondFlexibility', 'BondMobility', 'CMAPGrid', 'CMAPTorsion', 'ChargedAtomTypeIndex', 'CompoundAtomIndex', 'Context', 'CoordinateTransferError', 'Exclusion', 'IntegratorType', 'NonbondedMethod', 'ReferenceIndices', 'RoboAngle', 'RoboAtom', 'RoboAtomConnectivity', 'RoboAtomElement', 'RoboAtomIdentity', 'RoboAtomPhysics', 'RoboBond', 'RoboHarmonicImproperTorsion', 'RoboPeriodicTorsion', 'RoboPeriodicTorsionTerm', 'RootMobility', 'RunType', 'SamplerName', 'Scaling14', 'ThermostatName', 'TopologyRange', 'TopologyRangeType', 'UreyBradley', 'Vec3', 'World', 'ZMatrixRow', 'align_flip_and_translate_frame_along_x_axis', 'calculate_angle_in_rad', 'calculate_dihedral_in_rad', 'calculate_log_sum_exp2', 'calculate_mag_sq', 'chirality_from_plane_deviation', 'exceeds_planarity_threshold', 'flipped_chirality', 'is_bond_chirality_mismatch', 'is_chirality_mismatch', 'multiply_by_scalar', 'normalize_in_place', 'plane_normal', 'resolve_reference_indices', 'safe_log_sine_sqr', 'signed_plane_deviation', 'triple_product']
class AcceptRejectMode:
    """
    Members:
    
      AlwaysAccept
    
      MetropolisHastings
    """
    AlwaysAccept: typing.ClassVar[AcceptRejectMode]  # value = <AcceptRejectMode.AlwaysAccept: 0>
    MetropolisHastings: typing.ClassVar[AcceptRejectMode]  # value = <AcceptRejectMode.MetropolisHastings: 1>
    __members__: typing.ClassVar[dict[str, AcceptRejectMode]]  # value = {'AlwaysAccept': <AcceptRejectMode.AlwaysAccept: 0>, 'MetropolisHastings': <AcceptRejectMode.MetropolisHastings: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class AtomClassIndex:
    def __init__(self, arg0: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
class BondCenterChirality:
    """
    Members:
    
      RightHanded
    
      LeftHanded
    
      Planar
    """
    LeftHanded: typing.ClassVar[BondCenterChirality]  # value = <BondCenterChirality.LeftHanded: 1>
    Planar: typing.ClassVar[BondCenterChirality]  # value = <BondCenterChirality.Planar: 2>
    RightHanded: typing.ClassVar[BondCenterChirality]  # value = <BondCenterChirality.RightHanded: 0>
    __members__: typing.ClassVar[dict[str, BondCenterChirality]]  # value = {'RightHanded': <BondCenterChirality.RightHanded: 0>, 'LeftHanded': <BondCenterChirality.LeftHanded: 1>, 'Planar': <BondCenterChirality.Planar: 2>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class BondFlexibility:
    globalIndex1: int
    globalIndex2: int
    mobility: BondMobility
    uniqueAtomName1: str
    uniqueAtomName2: str
    def __init__(self) -> None:
        ...
class BondMobility:
    """
    Members:
    
      Free
    
      Torsion
    
      Rigid
    
      BallF
    
      BallM
    
      Cylinder
    
      Translation
    
      FreeLine
    
      LineOrientationF
    
      LineOrientationM
    
      UniversalM
    
      Spherical
    
      AnglePin
    
      BendStretch
    
      Slider
    
      OrthoSpherical
    """
    AnglePin: typing.ClassVar[BondMobility]  # value = <BondMobility.AnglePin: 13>
    BallF: typing.ClassVar[BondMobility]  # value = <BondMobility.BallF: 4>
    BallM: typing.ClassVar[BondMobility]  # value = <BondMobility.BallM: 5>
    BendStretch: typing.ClassVar[BondMobility]  # value = <BondMobility.BendStretch: 14>
    Cylinder: typing.ClassVar[BondMobility]  # value = <BondMobility.Cylinder: 6>
    Free: typing.ClassVar[BondMobility]  # value = <BondMobility.Free: 1>
    FreeLine: typing.ClassVar[BondMobility]  # value = <BondMobility.FreeLine: 8>
    LineOrientationF: typing.ClassVar[BondMobility]  # value = <BondMobility.LineOrientationF: 9>
    LineOrientationM: typing.ClassVar[BondMobility]  # value = <BondMobility.LineOrientationM: 10>
    OrthoSpherical: typing.ClassVar[BondMobility]  # value = <BondMobility.OrthoSpherical: 16>
    Rigid: typing.ClassVar[BondMobility]  # value = <BondMobility.Rigid: 3>
    Slider: typing.ClassVar[BondMobility]  # value = <BondMobility.Slider: 15>
    Spherical: typing.ClassVar[BondMobility]  # value = <BondMobility.Spherical: 12>
    Torsion: typing.ClassVar[BondMobility]  # value = <BondMobility.Torsion: 2>
    Translation: typing.ClassVar[BondMobility]  # value = <BondMobility.Translation: 7>
    UniversalM: typing.ClassVar[BondMobility]  # value = <BondMobility.UniversalM: 11>
    __members__: typing.ClassVar[dict[str, BondMobility]]  # value = {'Free': <BondMobility.Free: 1>, 'Torsion': <BondMobility.Torsion: 2>, 'Rigid': <BondMobility.Rigid: 3>, 'BallF': <BondMobility.BallF: 4>, 'BallM': <BondMobility.BallM: 5>, 'Cylinder': <BondMobility.Cylinder: 6>, 'Translation': <BondMobility.Translation: 7>, 'FreeLine': <BondMobility.FreeLine: 8>, 'LineOrientationF': <BondMobility.LineOrientationF: 9>, 'LineOrientationM': <BondMobility.LineOrientationM: 10>, 'UniversalM': <BondMobility.UniversalM: 11>, 'Spherical': <BondMobility.Spherical: 12>, 'AnglePin': <BondMobility.AnglePin: 13>, 'BendStretch': <BondMobility.BendStretch: 14>, 'Slider': <BondMobility.Slider: 15>, 'OrthoSpherical': <BondMobility.OrthoSpherical: 16>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class CMAPGrid:
    energy: list[float]
    size: int
    def __init__(self) -> None:
        ...
class CMAPTorsion:
    a1: int
    a2: int
    a3: int
    a4: int
    b1: int
    b2: int
    b3: int
    b4: int
    mapIndex: int
    def __init__(self) -> None:
        ...
class ChargedAtomTypeIndex:
    def __init__(self, arg0: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
class CompoundAtomIndex:
    def __init__(self, arg0: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
class Context:
    def __init__(self, arg0: str, arg1: int, arg2: int, arg3: RunType, arg4: int, arg5: int, arg6: bool) -> None:
        ...
    def addReplica(self) -> None:
        """
        Add an empty replica to the context.
        """
    def addThermodynamicState(self, arg0: float, arg1: list[AcceptRejectMode], arg2: list[int], arg3: list[str], arg4: list[int], arg5: list[int], arg6: list[IntegratorType], arg7: list[int], arg8: list[float], arg9: list[int]) -> None:
        """
        Add an empty themodynamic state to the context.
        """
    def add_world(self, arg0: bool, arg1: int, arg2: RootMobility, arg3: list[list[BondFlexibility]]) -> None:
        """
        Add an empty world.
        """
    def calculate_openmm_energy(self, worldIndex: int) -> float:
        """
        Calculate the OpenMM energy of the current state for a specific world index.
        """
    def getAtomNameByPrmtopIndex(self, prmtopIndex: int) -> str:
        """
        Get the unique atom name for a given prmtop index.
        """
    def getWorld(self, arg0: int) -> World:
        ...
    def getWorlds(self) -> list[World]:
        ...
    def initialize_openmm(self, arg0: list[RoboAtom], arg1: list[RoboBond], arg2: list[RoboAngle], arg3: list[RoboPeriodicTorsion], arg4: list[RoboHarmonicImproperTorsion], arg5: list[CMAPGrid], arg6: list[CMAPTorsion], arg7: list[UreyBradley], arg8: bool, arg9: int, arg10: list[float], arg11: list[float], arg12: list[Exclusion], arg13: list[Scaling14]) -> bool:
        """
        Load an OpenMM system from components.
        """
    def loadAmberSystem(self, arg0: list[int], arg1: list[RoboAtom], arg2: list[RoboBond], arg3: list[RoboAngle], arg4: list[RoboPeriodicTorsion], arg5: list[RoboHarmonicImproperTorsion], arg6: list[TopologyRange], arg7: list[ZMatrixRow]) -> None:
        """
        Load an AMBER system.
        """
    def run_rex(self, num_equilibration_rounds: int, num_production_rounds: int, write_frequency: int, write_to_stdio: bool) -> None:
        """
        Run replica exchange.
        """
    def setGBSAOptions(self, arg0: bool, arg1: float, arg2: float) -> None:
        """
        Set GBSA-OBC2 options.
        """
    def setNonbonded(self, arg0: NonbondedMethod, arg1: float) -> None:
        """
        Set nonbonded method and cutoff.
        """
    def setPdbRestartFreq(self, arg0: int) -> None:
        """
        Set the PDB restart frequency.
        """
    def setVerbose(self, arg0: bool) -> None:
        """
        Control if you want extraneous output to cout.
        """
    def validate_context(self) -> bool:
        """
        Validates all worlds and replicas in the context.
        """
class CoordinateTransferError:
    angles: float
    anglesMax: float
    bonds: float
    bondsMax: float
    cartesian: float
    cartesianMax: float
    improperDihedrals: float
    improperDihedralsMax: float
    matchResiduals: list[float]
    properDihedrals: float
    properDihedralsMax: float
class Exclusion:
    a1: int
    a2: int
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, a1: int, a2: int) -> None:
        ...
class IntegratorType:
    """
    Members:
    
      EMPTY
    
      VERLET
    
      EULER
    
      EULER2
    
      CPODES
    
      RUNGEKUTTA
    
      RUNGEKUTTA2
    
      RUNGEKUTTA3
    
      RUNGEKUTTAFELDBERG
    
      BENDSTRETCH
    
      OMMVV
    
      BOUND_WALK
    
      BOUND_HMC
    
      STATIONS_TASK
    
      NOF_INTEGRATORS
    """
    BENDSTRETCH: typing.ClassVar[IntegratorType]  # value = <IntegratorType.BENDSTRETCH: 9>
    BOUND_HMC: typing.ClassVar[IntegratorType]  # value = <IntegratorType.BOUND_HMC: 12>
    BOUND_WALK: typing.ClassVar[IntegratorType]  # value = <IntegratorType.BOUND_WALK: 11>
    CPODES: typing.ClassVar[IntegratorType]  # value = <IntegratorType.CPODES: 4>
    EMPTY: typing.ClassVar[IntegratorType]  # value = <IntegratorType.EMPTY: 0>
    EULER: typing.ClassVar[IntegratorType]  # value = <IntegratorType.EULER: 2>
    EULER2: typing.ClassVar[IntegratorType]  # value = <IntegratorType.EULER2: 3>
    NOF_INTEGRATORS: typing.ClassVar[IntegratorType]  # value = <IntegratorType.NOF_INTEGRATORS: 14>
    OMMVV: typing.ClassVar[IntegratorType]  # value = <IntegratorType.OMMVV: 10>
    RUNGEKUTTA: typing.ClassVar[IntegratorType]  # value = <IntegratorType.RUNGEKUTTA: 5>
    RUNGEKUTTA2: typing.ClassVar[IntegratorType]  # value = <IntegratorType.RUNGEKUTTA2: 6>
    RUNGEKUTTA3: typing.ClassVar[IntegratorType]  # value = <IntegratorType.RUNGEKUTTA3: 7>
    RUNGEKUTTAFELDBERG: typing.ClassVar[IntegratorType]  # value = <IntegratorType.RUNGEKUTTAFELDBERG: 8>
    STATIONS_TASK: typing.ClassVar[IntegratorType]  # value = <IntegratorType.STATIONS_TASK: 13>
    VERLET: typing.ClassVar[IntegratorType]  # value = <IntegratorType.VERLET: 1>
    __members__: typing.ClassVar[dict[str, IntegratorType]]  # value = {'EMPTY': <IntegratorType.EMPTY: 0>, 'VERLET': <IntegratorType.VERLET: 1>, 'EULER': <IntegratorType.EULER: 2>, 'EULER2': <IntegratorType.EULER2: 3>, 'CPODES': <IntegratorType.CPODES: 4>, 'RUNGEKUTTA': <IntegratorType.RUNGEKUTTA: 5>, 'RUNGEKUTTA2': <IntegratorType.RUNGEKUTTA2: 6>, 'RUNGEKUTTA3': <IntegratorType.RUNGEKUTTA3: 7>, 'RUNGEKUTTAFELDBERG': <IntegratorType.RUNGEKUTTAFELDBERG: 8>, 'BENDSTRETCH': <IntegratorType.BENDSTRETCH: 9>, 'OMMVV': <IntegratorType.OMMVV: 10>, 'BOUND_WALK': <IntegratorType.BOUND_WALK: 11>, 'BOUND_HMC': <IntegratorType.BOUND_HMC: 12>, 'STATIONS_TASK': <IntegratorType.STATIONS_TASK: 13>, 'NOF_INTEGRATORS': <IntegratorType.NOF_INTEGRATORS: 14>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class NonbondedMethod:
    """
    Members:
    
      NoCutoff
    
      CutoffNonPeriodic
    """
    CutoffNonPeriodic: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.CutoffNonPeriodic: 1>
    NoCutoff: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.NoCutoff: 0>
    __members__: typing.ClassVar[dict[str, NonbondedMethod]]  # value = {'NoCutoff': <NonbondedMethod.NoCutoff: 0>, 'CutoffNonPeriodic': <NonbondedMethod.CutoffNonPeriodic: 1>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class ReferenceIndices:
    one: int
    two: int
    zero: int
    def __init__(self, zero: int, one: int, two: int) -> None:
        ...
class RoboAngle:
    compound_atom_indices: typing.Annotated[list[CompoundAtomIndex], pybind11_stubgen.typing_ext.FixedSize(3)]
    global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]
    molecule_index: int
    nominal_angle_in_deg: float
    prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)]
    stiffness_in_kj_per_rad_sq: float
    def __init__(self, global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)], prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)], compound_atom_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(3)], molecule_index: int, stiffness_in_kj_per_rad_sq: float, nominal_angle_in_deg: float) -> None:
        ...
class RoboAtom:
    connectivity: RoboAtomConnectivity
    element_info: RoboAtomElement
    identity: RoboAtomIdentity
    physics: RoboAtomPhysics
    position: Vec3
    def __init__(self, identity: RoboAtomIdentity, element_info: RoboAtomElement, physics: RoboAtomPhysics, connectivity: RoboAtomConnectivity, position: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> None:
        ...
class RoboAtomConnectivity:
    neighbors_global_indices: list[int]
    root: bool
    def __init__(self, neighbors_global_indices: list[int], root: bool) -> None:
        ...
class RoboAtomElement:
    atomicNumber: int
    elementName: str
    elementSymbol: str
    def __init__(self, element_name: str, element_symbol: str, atomic_number: int) -> None:
        ...
class RoboAtomIdentity:
    atom_class_index: AtomClassIndex
    atom_class_name: str
    charged_atom_type_index: ChargedAtomTypeIndex
    charged_atom_type_name: str
    compound_atom_index: CompoundAtomIndex
    global_index: int
    molecule_index: int
    nonbonded_index: int
    prmtop_index: int
    residue_index: int
    residue_name: str
    unique_name: str
    def __init__(self, unique_name: str, residue_name: str, atom_class_name: str, charged_atom_type_name: str, global_index: int, prmtop_index: int, molecule_index: int, residue_index: int, nonbonded_index: int, compound_atom_index: int, atom_class_index: int, charged_atom_type_index: int) -> None:
        ...
class RoboAtomPhysics:
    charge_e: float
    mass_daltons: float
    screen: float
    sigma_nm: float
    solvent_radius_nm: float
    vdw_radius_nm: float
    vdw_well_depth_kj: float
    def __init__(self, charge_e: float, mass_daltons: float, vdw_radius_nm: float, vdw_well_depth_kj: float, sigma_nm: float, solvent_radius_nm: float, screen: float) -> None:
        ...
class RoboBond:
    compound_atom_indices: typing.Annotated[list[CompoundAtomIndex], pybind11_stubgen.typing_ext.FixedSize(2)]
    dihedral_type: str
    global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(2)]
    molecule_index: int
    nominal_length_in_nm: float
    prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(2)]
    ring_closing: bool
    stiffness_in_kj_per_nm_sq: float
    def __init__(self, global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(2)], prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(2)], compound_atom_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(2)], stiffness_in_kj_per_nm_sq: float, nominal_length_in_nm: float, molecule_index: int, ring_closing: bool, dihedral_type: str) -> None:
        ...
class RoboHarmonicImproperTorsion:
    compound_atom_indices: typing.Annotated[list[CompoundAtomIndex], pybind11_stubgen.typing_ext.FixedSize(4)]
    global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]
    molecule_index: int
    nominal_angle_in_rad: float
    prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]
    stiffness_in_kj_per_rad_sq: float
    def __init__(self, global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], compound_atom_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], molecule_index: int, stiffness_in_kj_per_rad_sq: float, nominal_angle_in_rad: float) -> None:
        ...
class RoboPeriodicTorsion:
    compound_atom_indices: typing.Annotated[list[CompoundAtomIndex], pybind11_stubgen.typing_ext.FixedSize(4)]
    global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]
    improper: bool
    molecule_index: int
    prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]
    terms: typing.Annotated[list[RoboPeriodicTorsionTerm], pybind11_stubgen.typing_ext.FixedSize(5)]
    def __init__(self, global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], prmtop_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], compound_atom_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], molecule_index: int, improper: bool, terms: list[RoboPeriodicTorsionTerm]) -> None:
        ...
class RoboPeriodicTorsionTerm:
    amplitude_kj: float
    periodicity: int
    phase_deg: float
    def __init__(self, amplitude_kj: float, phase_deg: float, periodicity: int) -> None:
        ...
class RootMobility:
    """
    Members:
    
      FREE
    
      CARTESIAN
    
      WELD
    
      FREE_LINE
    
      BALL
    
      PIN
    """
    BALL: typing.ClassVar[RootMobility]  # value = <RootMobility.BALL: 4>
    CARTESIAN: typing.ClassVar[RootMobility]  # value = <RootMobility.CARTESIAN: 1>
    FREE: typing.ClassVar[RootMobility]  # value = <RootMobility.FREE: 0>
    FREE_LINE: typing.ClassVar[RootMobility]  # value = <RootMobility.FREE_LINE: 3>
    PIN: typing.ClassVar[RootMobility]  # value = <RootMobility.PIN: 5>
    WELD: typing.ClassVar[RootMobility]  # value = <RootMobility.WELD: 2>
    __members__: typing.ClassVar[dict[str, RootMobility]]  # value = {'FREE': <RootMobility.FREE: 0>, 'CARTESIAN': <RootMobility.CARTESIAN: 1>, 'WELD': <RootMobility.WELD: 2>, 'FREE_LINE': <RootMobility.FREE_LINE: 3>, 'BALL': <RootMobility.BALL: 4>, 'PIN': <RootMobility.PIN: 5>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class RunType:
    """
    Members:
    
      DEFAULT
    
      REMC
    
      RENEMC
    
      RENE
    """
    DEFAULT: typing.ClassVar[RunType]  # value = <RunType.DEFAULT: 0>
    REMC: typing.ClassVar[RunType]  # value = <RunType.REMC: 1>
    RENE: typing.ClassVar[RunType]  # value = <RunType.RENE: 3>
    RENEMC: typing.ClassVar[RunType]  # value = <RunType.RENEMC: 2>
    __members__: typing.ClassVar[dict[str, RunType]]  # value = {'DEFAULT': <RunType.DEFAULT: 0>, 'REMC': <RunType.REMC: 1>, 'RENEMC': <RunType.RENEMC: 2>, 'RENE': <RunType.RENE: 3>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class SamplerName:
    """
    Members:
    
      EMPTY
    
      MC
    
      HMC
    
      LAHMC
    """
    EMPTY: typing.ClassVar[SamplerName]  # value = <SamplerName.EMPTY: 0>
    HMC: typing.ClassVar[SamplerName]  # value = <SamplerName.HMC: 2>
    LAHMC: typing.ClassVar[SamplerName]  # value = <SamplerName.LAHMC: 3>
    MC: typing.ClassVar[SamplerName]  # value = <SamplerName.MC: 1>
    __members__: typing.ClassVar[dict[str, SamplerName]]  # value = {'EMPTY': <SamplerName.EMPTY: 0>, 'MC': <SamplerName.MC: 1>, 'HMC': <SamplerName.HMC: 2>, 'LAHMC': <SamplerName.LAHMC: 3>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class Scaling14:
    a1: int
    a4: int
    charge_product: float
    epsilon: float
    sigma: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, a1: int, a4: int, charge_product: float, epsilon: float, sigma: float) -> None:
        ...
class ThermostatName:
    """
    Members:
    
      NONE
    
      ANDERSEN
    
      BERENDSEN
    
      LANGEVIN
    
      NOSE_HOOVER
    """
    ANDERSEN: typing.ClassVar[ThermostatName]  # value = <ThermostatName.ANDERSEN: 1>
    BERENDSEN: typing.ClassVar[ThermostatName]  # value = <ThermostatName.BERENDSEN: 2>
    LANGEVIN: typing.ClassVar[ThermostatName]  # value = <ThermostatName.LANGEVIN: 3>
    NONE: typing.ClassVar[ThermostatName]  # value = <ThermostatName.NONE: 0>
    NOSE_HOOVER: typing.ClassVar[ThermostatName]  # value = <ThermostatName.NOSE_HOOVER: 4>
    __members__: typing.ClassVar[dict[str, ThermostatName]]  # value = {'NONE': <ThermostatName.NONE: 0>, 'ANDERSEN': <ThermostatName.ANDERSEN: 1>, 'BERENDSEN': <ThermostatName.BERENDSEN: 2>, 'LANGEVIN': <ThermostatName.LANGEVIN: 3>, 'NOSE_HOOVER': <ThermostatName.NOSE_HOOVER: 4>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class TopologyRange:
    def __init__(self, startCounts: list[int]) -> None:
        ...
    def close(self, endCounts: list[int]) -> None:
        ...
class TopologyRangeType:
    """
    Members:
    
      EMPTY
    
      BOND
    
      ANGLE
    
      PERIODIC_TORSION
    
      IMPROPER_HARMONIC_TORSION
    """
    ANGLE: typing.ClassVar[TopologyRangeType]  # value = <TopologyRangeType.ANGLE: 2>
    BOND: typing.ClassVar[TopologyRangeType]  # value = <TopologyRangeType.BOND: 1>
    EMPTY: typing.ClassVar[TopologyRangeType]  # value = <TopologyRangeType.EMPTY: 0>
    IMPROPER_HARMONIC_TORSION: typing.ClassVar[TopologyRangeType]  # value = <TopologyRangeType.IMPROPER_HARMONIC_TORSION: 4>
    PERIODIC_TORSION: typing.ClassVar[TopologyRangeType]  # value = <TopologyRangeType.PERIODIC_TORSION: 3>
    __members__: typing.ClassVar[dict[str, TopologyRangeType]]  # value = {'EMPTY': <TopologyRangeType.EMPTY: 0>, 'BOND': <TopologyRangeType.BOND: 1>, 'ANGLE': <TopologyRangeType.ANGLE: 2>, 'PERIODIC_TORSION': <TopologyRangeType.PERIODIC_TORSION: 3>, 'IMPROPER_HARMONIC_TORSION': <TopologyRangeType.IMPROPER_HARMONIC_TORSION: 4>}
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: int) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: int) -> None:
        ...
    def __str__(self) -> str:
        ...
    @property
    def name(self) -> str:
        ...
    @property
    def value(self) -> int:
        ...
class UreyBradley:
    a1: int
    a3: int
    nominal_length_in_nm: float
    stiffness_in_kj_per_nm_sq: float
    def __init__(self) -> None:
        ...
class Vec3:
    def __getitem__(self, arg0: int) -> float:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: float, arg1: float, arg2: float) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: list[float]) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def __setitem__(self, arg0: int, arg1: float) -> None:
        ...
class World:
    def add_sampler(self, sampler_name: SamplerName, integrator_type: IntegratorType, thermostat_name: ThermostatName, use_fixman_potential: bool, use_nuts: bool) -> bool:
        """
        Add a sampler to the world.
        """
    def get_coordinate_transfer_errors(self) -> list[CoordinateTransferError]:
        """
        Get the coordinate transfer errors for all samplers in the world.
        """
    def has_rigid_body_violations(self, timeStep: float, numSteps: int) -> bool:
        """
        Checks for rigid body violations.
        """
class ZMatrixRow:
    compound_atom_indices: typing.Annotated[list[CompoundAtomIndex], pybind11_stubgen.typing_ext.FixedSize(4)]
    global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)]
    molecule_index: int
    def __init__(self, global_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], compound_atom_indices: typing.Annotated[list[int], pybind11_stubgen.typing_ext.FixedSize(4)], molecule_index: int) -> None:
        ...
def align_flip_and_translate_frame_along_x_axis(arg0: numpy.ndarray[numpy.float64], arg1: numpy.ndarray[numpy.float64]) -> numpy.ndarray[numpy.float64]:
    ...
def calculate_angle_in_rad(arg0: numpy.ndarray[numpy.float64], arg1: numpy.ndarray[numpy.float64], arg2: numpy.ndarray[numpy.float64]) -> float:
    ...
def calculate_dihedral_in_rad(arg0: numpy.ndarray[numpy.float64], arg1: numpy.ndarray[numpy.float64], arg2: numpy.ndarray[numpy.float64], arg3: numpy.ndarray[numpy.float64]) -> float:
    ...
def calculate_log_sum_exp2(arg0: float, arg1: float) -> float:
    ...
def calculate_mag_sq(arg0: list[float]) -> float:
    ...
def chirality_from_plane_deviation(arg0: float) -> BondCenterChirality:
    ...
def exceeds_planarity_threshold(arg0: float, arg1: float) -> bool:
    ...
def flipped_chirality(arg0: BondCenterChirality) -> BondCenterChirality:
    ...
def is_bond_chirality_mismatch(arg0: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg1: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg2: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg3: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg4: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg5: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> bool:
    ...
def is_chirality_mismatch(arg0: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg1: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg2: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg3: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg4: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg5: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> bool:
    ...
def multiply_by_scalar(arg0: list[float], arg1: float, arg2: list[float]) -> None:
    ...
def normalize_in_place(arg0: list[float]) -> None:
    ...
def plane_normal(arg0: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg1: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]:
    ...
def resolve_reference_indices(arg0: list[int]) -> ReferenceIndices:
    ...
def safe_log_sine_sqr(arg0: float) -> float:
    ...
def signed_plane_deviation(arg0: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg1: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> float:
    ...
def triple_product(arg0: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg1: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)], arg2: typing.Annotated[list[float], pybind11_stubgen.typing_ext.FixedSize(3)]) -> float:
    ...
