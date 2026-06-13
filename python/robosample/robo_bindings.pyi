"""
Robosample bindings
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['AcceptRejectMode', 'AtomClassIndex', 'BondCenterChirality', 'BondFlexibility', 'BondMobility', 'CMAPGrid', 'CMAPTorsion', 'ChargedAtomTypeIndex', 'CompoundAtomIndex', 'Context', 'CutoffNonPeriodic', 'Exclusion', 'ForceFieldParams', 'IntegratorType', 'NoCutoff', 'NonbondedMethod', 'ReferenceIndices', 'RoboAngle', 'RoboAtom', 'RoboAtomConnectivity', 'RoboAtomElement', 'RoboAtomIdentity', 'RoboAtomPhysics', 'RoboBond', 'RoboHarmonicImproperTorsion', 'RoboPeriodicTorsion', 'RoboPeriodicTorsionTerm', 'RootMobility', 'RunType', 'SamplerName', 'Scaling14', 'SimulationSettings', 'SystemTopology', 'ThermostatName', 'TopologyRange', 'TopologyRangeType', 'UreyBradley', 'Vec3', 'VectorCMAPGrid', 'VectorCMAPTorsion', 'VectorExclusion', 'VectorInt', 'VectorRoboAngle', 'VectorRoboAtom', 'VectorRoboBond', 'VectorRoboHarmonicImproperTorsion', 'VectorRoboPeriodicTorsion', 'VectorRootMobility', 'VectorScaling14', 'VectorTopologyRange', 'VectorUreyBradley', 'World', 'ZMatrixRow', 'align_flip_and_translate_frame_along_x_axis', 'calculate_angle_in_rad', 'calculate_dihedral_in_rad', 'calculate_log_sum_exp2', 'calculate_mag_sq', 'multiply_by_scalar', 'normalize_in_place', 'safe_log_sine_sqr']
class AcceptRejectMode:
    """
    Members:
    
      AlwaysAccept
    
      MetropolisHastings
    """
    AlwaysAccept: typing.ClassVar[AcceptRejectMode]  # value = <AcceptRejectMode.AlwaysAccept: 0>
    MetropolisHastings: typing.ClassVar[AcceptRejectMode]  # value = <AcceptRejectMode.MetropolisHastings: 1>
    __members__: typing.ClassVar[dict[str, AcceptRejectMode]]  # value = {'AlwaysAccept': <AcceptRejectMode.AlwaysAccept: 0>, 'MetropolisHastings': <AcceptRejectMode.MetropolisHastings: 1>}
    @typing.overload
    def __eq__(self, other: AcceptRejectMode) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: AcceptRejectMode) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __init__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: BondCenterChirality) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: BondCenterChirality) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    mobility: BondMobility
    uniqueAtomName1: str
    uniqueAtomName2: str
    def __init__(self) -> None:
        ...
    @property
    def globalIndex1(self) -> int:
        ...
    @globalIndex1.setter
    def globalIndex1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def globalIndex2(self) -> int:
        ...
    @globalIndex2.setter
    def globalIndex2(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: BondMobility) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: BondMobility) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __init__(self) -> None:
        ...
    @property
    def energy(self) -> list[float]:
        ...
    @energy.setter
    def energy(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def size(self) -> int:
        ...
    @size.setter
    def size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class CMAPTorsion:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, mapIndex: typing.SupportsInt | typing.SupportsIndex, torsion_a_atom_1_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_a_atom_2_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_a_atom_3_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_a_atom_4_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_b_atom_1_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_b_atom_2_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_b_atom_3_global_index: typing.SupportsInt | typing.SupportsIndex, torsion_b_atom_4_global_index: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def mapIndex(self) -> int:
        ...
    @mapIndex.setter
    def mapIndex(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_a_atom_1_global_index(self) -> int:
        ...
    @torsion_a_atom_1_global_index.setter
    def torsion_a_atom_1_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_a_atom_2_global_index(self) -> int:
        ...
    @torsion_a_atom_2_global_index.setter
    def torsion_a_atom_2_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_a_atom_3_global_index(self) -> int:
        ...
    @torsion_a_atom_3_global_index.setter
    def torsion_a_atom_3_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_a_atom_4_global_index(self) -> int:
        ...
    @torsion_a_atom_4_global_index.setter
    def torsion_a_atom_4_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_b_atom_1_global_index(self) -> int:
        ...
    @torsion_b_atom_1_global_index.setter
    def torsion_b_atom_1_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_b_atom_2_global_index(self) -> int:
        ...
    @torsion_b_atom_2_global_index.setter
    def torsion_b_atom_2_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_b_atom_3_global_index(self) -> int:
        ...
    @torsion_b_atom_3_global_index.setter
    def torsion_b_atom_3_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def torsion_b_atom_4_global_index(self) -> int:
        ...
    @torsion_b_atom_4_global_index.setter
    def torsion_b_atom_4_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class ChargedAtomTypeIndex:
    def __init__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
class CompoundAtomIndex:
    def __init__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
class Context:
    def __init__(self, arg0: str, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: RunType, arg4: typing.SupportsInt | typing.SupportsIndex, arg5: typing.SupportsInt | typing.SupportsIndex, arg6: bool) -> None:
        ...
    def addReplica(self) -> None:
        """
        Add an empty replica to the context.
        """
    def addThermodynamicState(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: collections.abc.Sequence[AcceptRejectMode], arg2: VectorInt, arg3: collections.abc.Sequence[str], arg4: VectorInt, arg5: VectorInt, arg6: collections.abc.Sequence[IntegratorType], arg7: VectorInt, arg8: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], arg9: VectorInt) -> None:
        """
        Add an empty themodynamic state to the context.
        """
    def add_world(self, fixman_torque: bool, samples_per_round: typing.SupportsInt | typing.SupportsIndex, roll_flexibilities: collections.abc.Sequence[collections.abc.Sequence[BondFlexibility]], want_spatial_force_history: bool) -> None:
        """
                    Add an empty world.
        
                    Args:
                        roll_flexibilities: A list of lists of BondFlexibility objects, 
                                            e.g., [[rb.BondFlexibility(), ...], [...]]
                        want_spatial_force_history: A boolean indicating whether to track spatial force history.
        """
    def calculate_openmm_energy(self, worldIndex: typing.SupportsInt | typing.SupportsIndex) -> float:
        """
        Calculate the OpenMM energy of the current state for a specific world index.
        """
    def getAtomNameByPrmtopIndex(self, prmtopIndex: typing.SupportsInt | typing.SupportsIndex) -> str:
        """
        Get the unique atom name for a given prmtop index.
        """
    def getWorld(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> World:
        ...
    def getWorlds(self) -> list[World]:
        ...
    def initialize_openmm(self) -> bool:
        """
        Load an OpenMM system from components.
        """
    def loadAmberSystem(self, arg0: SystemTopology, arg1: ForceFieldParams, arg2: SimulationSettings, arg3: collections.abc.Sequence[ZMatrixRow]) -> None:
        """
        Load an AMBER system.
        """
    def run_rex(self, num_equilibration_rounds: typing.SupportsInt | typing.SupportsIndex, num_production_rounds: typing.SupportsInt | typing.SupportsIndex, write_frequency: typing.SupportsInt | typing.SupportsIndex, write_to_stdio: bool) -> None:
        """
        Run replica exchange.
        """
    def setGBSAOptions(self, arg0: bool, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set GBSA-OBC2 options.
        """
    def setNonbonded(self, arg0: NonbondedMethod, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Set nonbonded method and cutoff.
        """
    def setPdbRestartFreq(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Set the PDB restart frequency.
        """
    def setVerbose(self, arg0: bool) -> None:
        """
        Control if you want extraneous output to cout.
        """
    def set_world_temperatures(self, world_temperatures_in_k: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        """
        Set the temperatures of the worlds in the context.
        """
class Exclusion:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: typing.SupportsInt | typing.SupportsIndex, atom_2_global_index: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_1_global_index(self) -> int:
        ...
    @atom_1_global_index.setter
    def atom_1_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_2_global_index(self) -> int:
        ...
    @atom_2_global_index.setter
    def atom_2_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class ForceFieldParams:
    has_nbfix: bool
    nonbonded_method: NonbondedMethod
    use_gbsaobc2: bool
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, has_nbfix: bool = False, num_types: typing.SupportsInt | typing.SupportsIndex = 0, a_coef: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] | None = None, b_coef: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] | None = None, use_gbsaobc2: bool = True, gbsa_solvent_dielectric: typing.SupportsFloat | typing.SupportsIndex = 78.5, gbsa_solute_dielectric: typing.SupportsFloat | typing.SupportsIndex = 1.0, nonbonded_method: robo_bindings.NonbondedMethod | None = None, nonbonded_cutoff_in_nm: typing.SupportsFloat | typing.SupportsIndex = 1.2) -> None:
        ...
    @property
    def a_coef(self) -> list[float]:
        ...
    @a_coef.setter
    def a_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def b_coef(self) -> list[float]:
        ...
    @b_coef.setter
    def b_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def gbsa_solute_dielectric(self) -> float:
        ...
    @gbsa_solute_dielectric.setter
    def gbsa_solute_dielectric(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def gbsa_solvent_dielectric(self) -> float:
        ...
    @gbsa_solvent_dielectric.setter
    def gbsa_solvent_dielectric(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def nonbonded_cutoff_in_nm(self) -> float:
        ...
    @nonbonded_cutoff_in_nm.setter
    def nonbonded_cutoff_in_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def num_types(self) -> int:
        ...
    @num_types.setter
    def num_types(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: IntegratorType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: IntegratorType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: NonbondedMethod) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: NonbondedMethod) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.SupportsInt | typing.SupportsIndex) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __init__(self, zero: typing.SupportsInt | typing.SupportsIndex, one: typing.SupportsInt | typing.SupportsIndex, two: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def one(self) -> int:
        ...
    @one.setter
    def one(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def two(self) -> int:
        ...
    @two.setter
    def two(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def zero(self) -> int:
        ...
    @zero.setter
    def zero(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class RoboAngle:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"], prmtop_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, stiffness_in_kj_per_rad_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_angle_in_deg: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[CompoundAtomIndex], "FixedSize(3)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[CompoundAtomIndex], "FixedSize(3)"]) -> None:
        ...
    @property
    def global_indices(self) -> typing.Annotated[list[int], "FixedSize(3)"]:
        ...
    @global_indices.setter
    def global_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nominal_angle_in_deg(self) -> float:
        ...
    @nominal_angle_in_deg.setter
    def nominal_angle_in_deg(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def prmtop_indices(self) -> typing.Annotated[list[int], "FixedSize(3)"]:
        ...
    @prmtop_indices.setter
    def prmtop_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
    @property
    def stiffness_in_kj_per_rad_sq(self) -> float:
        ...
    @stiffness_in_kj_per_rad_sq.setter
    def stiffness_in_kj_per_rad_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboAtom:
    connectivity: RoboAtomConnectivity
    element_info: RoboAtomElement
    identity: RoboAtomIdentity
    physics: RoboAtomPhysics
    position: Vec3
    def __init__(self, identity: RoboAtomIdentity, element_info: RoboAtomElement, physics: RoboAtomPhysics, connectivity: RoboAtomConnectivity, position: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
class RoboAtomConnectivity:
    neighbors_global_indices: VectorInt
    root: bool
    def __init__(self, neighbors_global_indices: VectorInt, root: bool) -> None:
        ...
class RoboAtomElement:
    elementName: str
    elementSymbol: str
    def __init__(self, element_name: str, element_symbol: str, atomic_number: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def atomicNumber(self) -> int:
        ...
    @atomicNumber.setter
    def atomicNumber(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class RoboAtomIdentity:
    atom_class_index: AtomClassIndex
    atom_class_name: str
    charged_atom_type_index: ChargedAtomTypeIndex
    charged_atom_type_name: str
    compound_atom_index: CompoundAtomIndex
    residue_name: str
    unique_name: str
    def __init__(self, unique_name: str, residue_name: str, atom_class_name: str, charged_atom_type_name: str, global_index: typing.SupportsInt | typing.SupportsIndex, prmtop_index: typing.SupportsInt | typing.SupportsIndex, molecule_index: typing.SupportsInt | typing.SupportsIndex, residue_index: typing.SupportsInt | typing.SupportsIndex, nonbonded_index: typing.SupportsInt | typing.SupportsIndex, compound_atom_index: typing.SupportsInt | typing.SupportsIndex, atom_class_index: typing.SupportsInt | typing.SupportsIndex, charged_atom_type_index: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def global_index(self) -> int:
        ...
    @global_index.setter
    def global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nonbonded_index(self) -> int:
        ...
    @nonbonded_index.setter
    def nonbonded_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def prmtop_index(self) -> int:
        ...
    @prmtop_index.setter
    def prmtop_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def residue_index(self) -> int:
        ...
    @residue_index.setter
    def residue_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class RoboAtomPhysics:
    def __init__(self, charge_e: typing.SupportsFloat | typing.SupportsIndex, mass_daltons: typing.SupportsFloat | typing.SupportsIndex, vdw_radius_nm: typing.SupportsFloat | typing.SupportsIndex, vdw_well_depth_kj: typing.SupportsFloat | typing.SupportsIndex, sigma_nm: typing.SupportsFloat | typing.SupportsIndex, solvent_radius_nm: typing.SupportsFloat | typing.SupportsIndex, screen: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def charge_e(self) -> float:
        ...
    @charge_e.setter
    def charge_e(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def mass_daltons(self) -> float:
        ...
    @mass_daltons.setter
    def mass_daltons(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def screen(self) -> float:
        ...
    @screen.setter
    def screen(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def sigma_nm(self) -> float:
        ...
    @sigma_nm.setter
    def sigma_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def solvent_radius_nm(self) -> float:
        ...
    @solvent_radius_nm.setter
    def solvent_radius_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def vdw_radius_nm(self) -> float:
        ...
    @vdw_radius_nm.setter
    def vdw_radius_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def vdw_well_depth_kj(self) -> float:
        ...
    @vdw_well_depth_kj.setter
    def vdw_well_depth_kj(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboBond:
    dihedral_type: str
    ring_closing: bool
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"], prmtop_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"], stiffness_in_kj_per_nm_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_length_in_nm: typing.SupportsFloat | typing.SupportsIndex, molecule_index: typing.SupportsInt | typing.SupportsIndex, ring_closing: bool, dihedral_type: str) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[CompoundAtomIndex], "FixedSize(2)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[CompoundAtomIndex], "FixedSize(2)"]) -> None:
        ...
    @property
    def global_indices(self) -> typing.Annotated[list[int], "FixedSize(2)"]:
        ...
    @global_indices.setter
    def global_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"]) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nominal_length_in_nm(self) -> float:
        ...
    @nominal_length_in_nm.setter
    def nominal_length_in_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def prmtop_indices(self) -> typing.Annotated[list[int], "FixedSize(2)"]:
        ...
    @prmtop_indices.setter
    def prmtop_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"]) -> None:
        ...
    @property
    def stiffness_in_kj_per_nm_sq(self) -> float:
        ...
    @stiffness_in_kj_per_nm_sq.setter
    def stiffness_in_kj_per_nm_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboHarmonicImproperTorsion:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], prmtop_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, stiffness_in_kj_per_rad_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_angle_in_rad: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[CompoundAtomIndex], "FixedSize(4)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[CompoundAtomIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def global_indices(self) -> typing.Annotated[list[int], "FixedSize(4)"]:
        ...
    @global_indices.setter
    def global_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nominal_angle_in_rad(self) -> float:
        ...
    @nominal_angle_in_rad.setter
    def nominal_angle_in_rad(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def prmtop_indices(self) -> typing.Annotated[list[int], "FixedSize(4)"]:
        ...
    @prmtop_indices.setter
    def prmtop_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def stiffness_in_kj_per_rad_sq(self) -> float:
        ...
    @stiffness_in_kj_per_rad_sq.setter
    def stiffness_in_kj_per_rad_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboPeriodicTorsion:
    improper: bool
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], prmtop_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, improper: bool, terms: collections.abc.Sequence[RoboPeriodicTorsionTerm]) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[CompoundAtomIndex], "FixedSize(4)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[CompoundAtomIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def global_indices(self) -> typing.Annotated[list[int], "FixedSize(4)"]:
        ...
    @global_indices.setter
    def global_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def prmtop_indices(self) -> typing.Annotated[list[int], "FixedSize(4)"]:
        ...
    @prmtop_indices.setter
    def prmtop_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def terms(self) -> typing.Annotated[list[RoboPeriodicTorsionTerm], "FixedSize(5)"]:
        ...
    @terms.setter
    def terms(self, arg0: typing.Annotated[collections.abc.Sequence[RoboPeriodicTorsionTerm], "FixedSize(5)"]) -> None:
        ...
class RoboPeriodicTorsionTerm:
    def __init__(self, amplitude_kj: typing.SupportsFloat | typing.SupportsIndex, phase_deg: typing.SupportsFloat | typing.SupportsIndex, periodicity: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def amplitude_kj(self) -> float:
        ...
    @amplitude_kj.setter
    def amplitude_kj(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def periodicity(self) -> int:
        ...
    @periodicity.setter
    def periodicity(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def phase_deg(self) -> float:
        ...
    @phase_deg.setter
    def phase_deg(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: RootMobility) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: RootMobility) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: RunType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: RunType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: SamplerName) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: SamplerName) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: typing.SupportsInt | typing.SupportsIndex, atom_4_global_index: typing.SupportsInt | typing.SupportsIndex, charge_product: typing.SupportsFloat | typing.SupportsIndex, epsilon: typing.SupportsFloat | typing.SupportsIndex, sigma: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_1_global_index(self) -> int:
        ...
    @atom_1_global_index.setter
    def atom_1_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_4_global_index(self) -> int:
        ...
    @atom_4_global_index.setter
    def atom_4_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def charge_product(self) -> float:
        ...
    @charge_product.setter
    def charge_product(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def epsilon(self) -> float:
        ...
    @epsilon.setter
    def epsilon(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def sigma(self) -> float:
        ...
    @sigma.setter
    def sigma(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class SimulationSettings:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, thermostat_temperature_in_k: typing.SupportsFloat | typing.SupportsIndex = 300.0, collision_frequency: typing.SupportsFloat | typing.SupportsIndex = 1.0, seed: typing.SupportsInt | typing.SupportsIndex = 0) -> None:
        ...
    @property
    def collision_frequency(self) -> float:
        ...
    @collision_frequency.setter
    def collision_frequency(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def seed(self) -> int:
        ...
    @seed.setter
    def seed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def thermostat_temperature_in_k(self) -> float:
        ...
    @thermostat_temperature_in_k.setter
    def thermostat_temperature_in_k(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
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
    def __init__(self, root_atom_global_indices: robo_bindings.VectorInt | None = None, topology_ranges: robo_bindings.VectorTopologyRange | None = None, atoms: robo_bindings.VectorRoboAtom | None = None, bonds: robo_bindings.VectorRoboBond | None = None, angles: robo_bindings.VectorRoboAngle | None = None, periodic_torsions: robo_bindings.VectorRoboPeriodicTorsion | None = None, harmonic_improper_torsions: robo_bindings.VectorRoboHarmonicImproperTorsion | None = None, cmap_grids: robo_bindings.VectorCMAPGrid | None = None, cmap_torsions: robo_bindings.VectorCMAPTorsion | None = None, urey_bradleys: robo_bindings.VectorUreyBradley | None = None, scaling14s: robo_bindings.VectorScaling14 | None = None, exclusions: robo_bindings.VectorExclusion | None = None, root_mobilities: robo_bindings.VectorRootMobility | None = None) -> None:
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
    @typing.overload
    def __eq__(self, other: ThermostatName) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: ThermostatName) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __eq__(self, other: TopologyRangeType) -> bool:
        ...
    @typing.overload
    def __eq__(self, other: typing.Any) -> bool:
        ...
    def __getstate__(self) -> int:
        ...
    def __hash__(self) -> int:
        ...
    def __index__(self) -> int:
        ...
    def __init__(self, value: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def __int__(self) -> int:
        ...
    @typing.overload
    def __ne__(self, other: TopologyRangeType) -> bool:
        ...
    @typing.overload
    def __ne__(self, other: typing.Any) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    def __setstate__(self, state: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, atom_1_global_index: typing.SupportsInt | typing.SupportsIndex, atom_3_global_index: typing.SupportsInt | typing.SupportsIndex, stiffness_in_kj_per_nm_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_length_in_nm: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_1_global_index(self) -> int:
        ...
    @atom_1_global_index.setter
    def atom_1_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def atom_3_global_index(self) -> int:
        ...
    @atom_3_global_index.setter
    def atom_3_global_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def nominal_length_in_nm(self) -> float:
        ...
    @nominal_length_in_nm.setter
    def nominal_length_in_nm(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def stiffness_in_kj_per_nm_sq(self) -> float:
        ...
    @stiffness_in_kj_per_nm_sq.setter
    def stiffness_in_kj_per_nm_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class Vec3:
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> float:
        ...
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __init__(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    def __repr__(self) -> str:
        ...
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class VectorCMAPGrid:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> CMAPGrid:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[CMAPGrid]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: CMAPGrid) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: CMAPGrid) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> CMAPGrid:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> CMAPGrid:
        """
        Remove and return the item at index ``i``
        """
class VectorCMAPTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> CMAPTorsion:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[CMAPTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: CMAPTorsion) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: CMAPTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> CMAPTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> CMAPTorsion:
        """
        Remove and return the item at index ``i``
        """
class VectorExclusion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> Exclusion:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[Exclusion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: Exclusion) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: Exclusion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> Exclusion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> Exclusion:
        """
        Remove and return the item at index ``i``
        """
class VectorInt:
    __hash__: typing.ClassVar[None] = None
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    def __contains__(self, x: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Return true the container contains ``x``
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> int:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[int]:
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
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @typing.overload
    def __setitem__(self, arg0: slice, arg1: VectorInt) -> None:
        """
        Assign list elements using a slice object
        """
    def append(self, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Add an item to the end of the list
        """
    def clear(self) -> None:
        """
        Clear the contents
        """
    def count(self, x: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Return the number of times ``x`` appears in the list
        """
    @typing.overload
    def extend(self, L: VectorInt) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    @typing.overload
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> int:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> int:
        """
        Remove and return the item at index ``i``
        """
    def remove(self, x: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Remove the first item from the list whose value is x. It is an error if there is no such item.
        """
class VectorRoboAngle:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RoboAngle:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RoboAngle]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RoboAngle) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RoboAngle) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboAngle:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RoboAngle:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboAtom:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RoboAtom:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RoboAtom]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RoboAtom) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RoboAtom) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboAtom:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RoboAtom:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboBond:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RoboBond:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RoboBond]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RoboBond) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RoboBond) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboBond:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RoboBond:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboHarmonicImproperTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RoboHarmonicImproperTorsion:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RoboHarmonicImproperTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RoboHarmonicImproperTorsion) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RoboHarmonicImproperTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboHarmonicImproperTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RoboHarmonicImproperTorsion:
        """
        Remove and return the item at index ``i``
        """
class VectorRoboPeriodicTorsion:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RoboPeriodicTorsion:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RoboPeriodicTorsion]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RoboPeriodicTorsion) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RoboPeriodicTorsion) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RoboPeriodicTorsion:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RoboPeriodicTorsion:
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
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> RootMobility:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[RootMobility]:
        ...
    def __len__(self) -> int:
        ...
    def __ne__(self, arg0: VectorRootMobility) -> bool:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: RootMobility) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: RootMobility) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> RootMobility:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> RootMobility:
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
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> Scaling14:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[Scaling14]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: Scaling14) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: Scaling14) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> Scaling14:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> Scaling14:
        """
        Remove and return the item at index ``i``
        """
class VectorTopologyRange:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> TopologyRange:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[TopologyRange]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: TopologyRange) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: TopologyRange) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> TopologyRange:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> TopologyRange:
        """
        Remove and return the item at index ``i``
        """
class VectorUreyBradley:
    def __bool__(self) -> bool:
        """
        Check whether the list is nonempty
        """
    @typing.overload
    def __delitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __getitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> UreyBradley:
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
    def __init__(self, arg0: collections.abc.Iterable) -> None:
        ...
    def __iter__(self) -> collections.abc.Iterator[UreyBradley]:
        ...
    def __len__(self) -> int:
        ...
    def __repr__(self) -> str:
        ...
    @typing.overload
    def __setitem__(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: UreyBradley) -> None:
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
    def extend(self, L: collections.abc.Iterable) -> None:
        """
        Extend the list by appending all the items in the given list
        """
    def insert(self, i: typing.SupportsInt | typing.SupportsIndex, x: UreyBradley) -> None:
        """
        Insert an item at a given position.
        """
    @typing.overload
    def pop(self) -> UreyBradley:
        """
        Remove and return the last item
        """
    @typing.overload
    def pop(self, i: typing.SupportsInt | typing.SupportsIndex) -> UreyBradley:
        """
        Remove and return the item at index ``i``
        """
class World:
    def add_sampler(self, sampler_name: SamplerName, integrator_type: IntegratorType, thermostat_name: ThermostatName, use_fixman_potential: bool, use_nuts: bool) -> bool:
        """
        Add a sampler to the world.
        """
class ZMatrixRow:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], molecule_index: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[CompoundAtomIndex], "FixedSize(4)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[CompoundAtomIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def global_indices(self) -> typing.Annotated[list[int], "FixedSize(4)"]:
        ...
    @global_indices.setter
    def global_indices(self, arg0: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"]) -> None:
        ...
    @property
    def molecule_index(self) -> int:
        ...
    @molecule_index.setter
    def molecule_index(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
def align_flip_and_translate_frame_along_x_axis(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> numpy.typing.NDArray[numpy.float64]:
    ...
def calculate_angle_in_rad(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> float:
    ...
def calculate_dihedral_in_rad(arg0: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg1: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg2: typing.Annotated[numpy.typing.ArrayLike, numpy.float64], arg3: typing.Annotated[numpy.typing.ArrayLike, numpy.float64]) -> float:
    ...
def calculate_log_sum_exp2(arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: typing.SupportsFloat | typing.SupportsIndex) -> float:
    ...
def calculate_mag_sq(arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> float:
    ...
def multiply_by_scalar(arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], arg1: typing.SupportsFloat | typing.SupportsIndex, arg2: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
    ...
def normalize_in_place(arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
    ...
def safe_log_sine_sqr(arg0: typing.SupportsFloat | typing.SupportsIndex) -> float:
    ...
CutoffNonPeriodic: NonbondedMethod  # value = <NonbondedMethod.CutoffNonPeriodic: 1>
NoCutoff: NonbondedMethod  # value = <NonbondedMethod.NoCutoff: 0>
