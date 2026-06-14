"""
Robosample bindings
"""
from __future__ import annotations
import numpy
import pybind11_stubgen.typing_ext
import typing
__all__: list[str] = ['AcceptRejectMode', 'AtomClassIndex', 'BondCenterChirality', 'BondFlexibility', 'BondMobility', 'CMAPGrid', 'CMAPTorsion', 'ChargedAtomTypeIndex', 'CompoundAtomIndex', 'Context', 'CoordinateTransferError', 'CutoffNonPeriodic', 'Exclusion', 'ForceFieldParams', 'IntegratorType', 'NoCutoff', 'NonbondedMethod', 'ReferenceIndices', 'RoboAngle', 'RoboAtom', 'RoboAtomConnectivity', 'RoboAtomElement', 'RoboAtomIdentity', 'RoboAtomPhysics', 'RoboBond', 'RoboHarmonicImproperTorsion', 'RoboPeriodicTorsion', 'RoboPeriodicTorsionTerm', 'RootMobility', 'RunType', 'SamplerName', 'Scaling14', 'SimulationSettings', 'SystemTopology', 'ThermostatName', 'TopologyRange', 'TopologyRangeType', 'UreyBradley', 'Vec3', 'VectorCMAPGrid', 'VectorCMAPTorsion', 'VectorExclusion', 'VectorInt', 'VectorRoboAngle', 'VectorRoboAtom', 'VectorRoboBond', 'VectorRoboHarmonicImproperTorsion', 'VectorRoboPeriodicTorsion', 'VectorRootMobility', 'VectorScaling14', 'VectorTopologyRange', 'VectorUreyBradley', 'World', 'ZMatrixRow', 'align_flip_and_translate_frame_along_x_axis', 'calculate_angle_in_rad', 'calculate_dihedral_in_rad', 'calculate_log_sum_exp2', 'calculate_mag_sq', 'multiply_by_scalar', 'normalize_in_place', 'safe_log_sine_sqr']
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
    mapIndex: int
    torsion_a_atom_1_global_index: int
    torsion_a_atom_2_global_index: int
    torsion_a_atom_3_global_index: int
    torsion_a_atom_4_global_index: int
    torsion_b_atom_1_global_index: int
    torsion_b_atom_2_global_index: int
    torsion_b_atom_3_global_index: int
    torsion_b_atom_4_global_index: int
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, mapIndex: int, torsion_a_atom_1_global_index: int, torsion_a_atom_2_global_index: int, torsion_a_atom_3_global_index: int, torsion_a_atom_4_global_index: int, torsion_b_atom_1_global_index: int, torsion_b_atom_2_global_index: int, torsion_b_atom_3_global_index: int, torsion_b_atom_4_global_index: int) -> None:
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
    def addThermodynamicState(self, arg0: float, arg1: list[AcceptRejectMode], arg2: VectorInt, arg3: list[str], arg4: VectorInt, arg5: VectorInt, arg6: list[IntegratorType], arg7: VectorInt, arg8: list[float], arg9: VectorInt) -> None:
        """
        Add an empty themodynamic state to the context.
        """
    def add_world(self, fixman_torque: bool, samples_per_round: int, roll_flexibilities: list[list[BondFlexibility]], want_spatial_force_history: bool) -> None:
        """
                    Add an empty world.
        
                    Args:
                        roll_flexibilities: A list of lists of BondFlexibility objects, 
                                            e.g., [[rb.BondFlexibility(), ...], [...]]
                        want_spatial_force_history: A boolean indicating whether to track spatial force history.
                        sphere_radius_in_nm: A real number indicating the sphere radius in nanometers.
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
    def initialize_openmm(self) -> bool:
        """
        Load an OpenMM system from components.
        """
    def loadAmberSystem(self, arg0: SystemTopology, arg1: ForceFieldParams, arg2: SimulationSettings, arg3: list[ZMatrixRow]) -> None:
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
    def set_world_temperatures(self, world_temperatures_in_k: list[float]) -> None:
        """
        Set the temperatures of the worlds in the context.
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
    atom_1_global_index: int
    atom_2_global_index: int
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: int, atom_2_global_index: int) -> None:
        ...
class ForceFieldParams:
    a_coef: list[float]
    b_coef: list[float]
    gbsa_solute_dielectric: float
    gbsa_solvent_dielectric: float
    has_nbfix: bool
    nonbonded_cutoff_in_nm: float
    nonbonded_method: NonbondedMethod
    num_types: int
    use_gbsaobc2: bool
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, has_nbfix: bool = False, num_types: int = 0, a_coef: list[float] | None = None, b_coef: list[float] | None = None, use_gbsaobc2: bool = True, gbsa_solvent_dielectric: float = 78.5, gbsa_solute_dielectric: float = 1.0, nonbonded_method: NonbondedMethod | None = None, nonbonded_cutoff_in_nm: float = 1.2) -> None:
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
    neighbors_global_indices: VectorInt
    root: bool
    def __init__(self, neighbors_global_indices: VectorInt, root: bool) -> None:
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
    
      Free
    
      Cartesian
    
      Weld
    
      FreeLine
    
      Ball
    
      Pin
    """
    Ball: typing.ClassVar[RootMobility]  # value = <RootMobility.Ball: 4>
    Cartesian: typing.ClassVar[RootMobility]  # value = <RootMobility.Cartesian: 1>
    Free: typing.ClassVar[RootMobility]  # value = <RootMobility.Free: 0>
    FreeLine: typing.ClassVar[RootMobility]  # value = <RootMobility.FreeLine: 3>
    Pin: typing.ClassVar[RootMobility]  # value = <RootMobility.Pin: 5>
    Weld: typing.ClassVar[RootMobility]  # value = <RootMobility.Weld: 2>
    __members__: typing.ClassVar[dict[str, RootMobility]]  # value = {'Free': <RootMobility.Free: 0>, 'Cartesian': <RootMobility.Cartesian: 1>, 'Weld': <RootMobility.Weld: 2>, 'FreeLine': <RootMobility.FreeLine: 3>, 'Ball': <RootMobility.Ball: 4>, 'Pin': <RootMobility.Pin: 5>}
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
    atom_1_global_index: int
    atom_4_global_index: int
    charge_product: float
    epsilon: float
    sigma: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: int, atom_4_global_index: int, charge_product: float, epsilon: float, sigma: float) -> None:
        ...
class SimulationSettings:
    collision_frequency: float
    seed: int
    thermostat_temperature_in_k: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, thermostat_temperature_in_k: float = 300.0, collision_frequency: float = 1.0, seed: int = 0) -> None:
        ...
class SystemTopology:
    angles: VectorRoboAngle
    atoms: VectorRoboAtom
    bonds: VectorRoboBond
    cmap_grids: VectorCMAPGrid
    cmap_torsions: VectorCMAPTorsion
    exclusions: VectorExclusion
    harmonic_improper_torsions: VectorRoboHarmonicImproperTorsion
    periodic_torsions: VectorRoboPeriodicTorsion
    root_atom_global_indices: VectorInt
    root_mobilities: VectorRootMobility
    scaling14s: VectorScaling14
    topology_ranges: VectorTopologyRange
    urey_bradleys: VectorUreyBradley
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, root_atom_global_indices: VectorInt | None = None, topology_ranges: VectorTopologyRange | None = None, atoms: VectorRoboAtom | None = None, bonds: VectorRoboBond | None = None, angles: VectorRoboAngle | None = None, periodic_torsions: VectorRoboPeriodicTorsion | None = None, harmonic_improper_torsions: VectorRoboHarmonicImproperTorsion | None = None, cmap_grids: VectorCMAPGrid | None = None, cmap_torsions: VectorCMAPTorsion | None = None, urey_bradleys: VectorUreyBradley | None = None, scaling14s: VectorScaling14 | None = None, exclusions: VectorExclusion | None = None, root_mobilities: VectorRootMobility | None = None) -> None:
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
    def __init__(self, startCounts: VectorInt) -> None:
        ...
    def close(self, endCounts: VectorInt) -> None:
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
    atom_1_global_index: int
    atom_3_global_index: int
    nominal_length_in_nm: float
    stiffness_in_kj_per_nm_sq: float
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: int, atom_3_global_index: int, stiffness_in_kj_per_nm_sq: float, nominal_length_in_nm: float) -> None:
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
class VectorCMAPGrid:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorCMAPGrid:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> CMAPGrid:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorCMAPGrid) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[CMAPGrid]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: CMAPGrid) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorCMAPGrid) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: CMAPGrid) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorCMAPGrid) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: CMAPGrid) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> CMAPGrid:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> CMAPGrid:
        """
        Remove and return the item at index ``i``
        """
class VectorCMAPTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorCMAPTorsion:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> CMAPTorsion:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorCMAPTorsion) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[CMAPTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: CMAPTorsion) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorCMAPTorsion) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: CMAPTorsion) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorCMAPTorsion) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: CMAPTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> CMAPTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> CMAPTorsion:
        """
        Remove and return the item at index ``i``
        """
class VectorExclusion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorExclusion:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> Exclusion:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorExclusion) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[Exclusion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: Exclusion) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorExclusion) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: Exclusion) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorExclusion) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: Exclusion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> Exclusion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> Exclusion:
        """
        Remove and return the item at index ``i``
        """
class VectorInt:
    __hash__: typing.ClassVar[None] = None
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: int) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    def __eq__(self, arg0: VectorInt) -> bool:
        ...
    @typing.overload
    def __getitem__(self, s: slice) -> VectorInt:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> int:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorInt) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[int]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: VectorInt) -> bool:
        ...
    @typing.overload
    def __repr__(self) -> str:
        """
        Return the canonical string representation of this list.
        """
    @typing.overload
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: int) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorInt) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: int) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: int) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: VectorInt) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: int) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> int:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> int:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: int) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
class VectorRoboAngle:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRoboAngle:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RoboAngle:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRoboAngle) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RoboAngle]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RoboAngle) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRoboAngle) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RoboAngle) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorRoboAngle) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RoboAngle) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboAngle:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RoboAngle:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboAtom:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRoboAtom:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RoboAtom:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRoboAtom) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RoboAtom]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RoboAtom) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRoboAtom) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RoboAtom) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorRoboAtom) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RoboAtom) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboAtom:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RoboAtom:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboBond:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRoboBond:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RoboBond:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRoboBond) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RoboBond]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RoboBond) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRoboBond) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RoboBond) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorRoboBond) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RoboBond) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboBond:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RoboBond:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboHarmonicImproperTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRoboHarmonicImproperTorsion:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RoboHarmonicImproperTorsion:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRoboHarmonicImproperTorsion) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RoboHarmonicImproperTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RoboHarmonicImproperTorsion) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRoboHarmonicImproperTorsion) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RoboHarmonicImproperTorsion) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorRoboHarmonicImproperTorsion) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RoboHarmonicImproperTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboHarmonicImproperTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RoboHarmonicImproperTorsion:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboPeriodicTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRoboPeriodicTorsion:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RoboPeriodicTorsion:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRoboPeriodicTorsion) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RoboPeriodicTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RoboPeriodicTorsion) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRoboPeriodicTorsion) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RoboPeriodicTorsion) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorRoboPeriodicTorsion) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RoboPeriodicTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboPeriodicTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RoboPeriodicTorsion:
        """
        Remove and return the item at index ``i``
        """
class VectorRootMobility:
    __hash__: typing.ClassVar[None] = None
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: RootMobility) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    def __eq__(self, arg0: VectorRootMobility) -> bool:
        ...
    @typing.overload
    def __getitem__(self, s: slice) -> VectorRootMobility:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> RootMobility:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorRootMobility) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[RootMobility]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: VectorRootMobility) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: RootMobility) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorRootMobility) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: RootMobility) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: RootMobility) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: VectorRootMobility) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: RootMobility) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RootMobility:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> RootMobility:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: RootMobility) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
class VectorScaling14:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorScaling14:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> Scaling14:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorScaling14) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[Scaling14]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: Scaling14) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorScaling14) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: Scaling14) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorScaling14) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: Scaling14) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> Scaling14:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> Scaling14:
        """
        Remove and return the item at index ``i``
        """
class VectorTopologyRange:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorTopologyRange:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> TopologyRange:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorTopologyRange) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[TopologyRange]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: TopologyRange) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorTopologyRange) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: TopologyRange) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorTopologyRange) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: TopologyRange) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> TopologyRange:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> TopologyRange:
        """
        Remove and return the item at index ``i``
        """
class VectorUreyBradley:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: int) -> None:
        """
        Delete the list elements at index ``i``
        """
    @typing.overload
    def __delitem__(self, arg0: slice) -> None:
        """
        Delete list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, s: slice) -> VectorUreyBradley:
        """
        Retrieve list elements using a slice object
        """
    @typing.overload
    def __getitem__(self, arg0: int) -> UreyBradley:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: VectorUreyBradley) -> None:
        """
        Copy constructor
        """
    @typing.overload
    def __init__(self, arg0: typing.Iterable) -> None:
        ...
    def __iter__(self) -> typing.Iterator[UreyBradley]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: int, arg1: UreyBradley) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorUreyBradley) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: UreyBradley) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    @typing.overload
    def extend(self, L: VectorUreyBradley) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: typing.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: int, x: UreyBradley) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> UreyBradley:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: int) -> UreyBradley:
        """
        Remove and return the item at index ``i``
        """
class World:
    def add_sampler(self, sampler_name: SamplerName, integrator_type: IntegratorType, thermostat_name: ThermostatName, sphere_radius_in_nm: typing.SupportsFloat | typing.SupportsIndex, use_fixman_potential: bool, use_nuts: bool) -> bool:
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
def multiply_by_scalar(arg0: list[float], arg1: float, arg2: list[float]) -> None:
    ...
def normalize_in_place(arg0: list[float]) -> None:
    ...
def safe_log_sine_sqr(arg0: float) -> float:
    ...
CutoffNonPeriodic: NonbondedMethod  # value = <NonbondedMethod.CutoffNonPeriodic: 1>
NoCutoff: NonbondedMethod  # value = <NonbondedMethod.NoCutoff: 0>
