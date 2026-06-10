"""
Robosample bindings
"""
from __future__ import annotations
import collections.abc
import numpy
import numpy.typing
import typing
__all__: list[str] = ['AcceptRejectMode', 'AtomClassIndex', 'BondCenterChirality', 'BondFlexibility', 'BondMobility', 'CMAPGrid', 'CMAPTorsion', 'ChargedAtomTypeIndex', 'CompoundAtomIndex', 'Context', 'CoordinateTransferError', 'CutoffNonPeriodic', 'Exclusion', 'ForceFieldParams', 'IntegratorType', 'NoCutoff', 'NonbondedMethod', 'ReferenceIndices', 'RoboAngle', 'RoboAtom', 'RoboAtomConnectivity', 'RoboAtomElement', 'RoboAtomIdentity', 'RoboAtomPhysics', 'RoboBond', 'RoboHarmonicImproperTorsion', 'RoboPeriodicTorsion', 'RoboPeriodicTorsionTerm', 'RootMobility', 'RunType', 'SamplerName', 'Scaling14', 'SimulationSettings', 'SystemTopology', 'SystemTopology_NEW', 'ThermostatName', 'TopologyRange', 'TopologyRangeType', 'UreyBradley', 'Vec3', 'VectorCMAPGrid', 'VectorCMAPTorsion', 'VectorExclusion', 'VectorInt', 'VectorRoboAngle', 'VectorRoboAtom', 'VectorRoboBond', 'VectorRoboHarmonicImproperTorsion', 'VectorRoboPeriodicTorsion', 'VectorRootMobility', 'VectorScaling14', 'VectorTopologyRange', 'VectorUreyBradley', 'World', 'ZMatrixRow', 'align_flip_and_translate_frame_along_x_axis', 'calculate_angle_in_rad', 'calculate_dihedral_in_rad', 'calculate_log_sum_exp2', 'calculate_mag_sq', 'multiply_by_scalar', 'normalize_in_place', 'safe_log_sine_sqr']
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
    def __init__(self, base_name: str, seed: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def run_rex(self, run_type: RunType, num_equilibration_rounds: typing.SupportsInt | typing.SupportsIndex, num_production_rounds: typing.SupportsInt | typing.SupportsIndex, write_frequency: typing.SupportsInt | typing.SupportsIndex, write_to_stdio: bool) -> None:
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
    def validate_context(self) -> bool:
        """
        Validates all worlds and replicas in the context.
        """
class CoordinateTransferError:
    @property
    def angles(self) -> float:
        ...
    @angles.setter
    def angles(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def anglesMax(self) -> float:
        ...
    @anglesMax.setter
    def anglesMax(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def bonds(self) -> float:
        ...
    @bonds.setter
    def bonds(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def bondsMax(self) -> float:
        ...
    @bondsMax.setter
    def bondsMax(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def cartesian(self) -> float:
        ...
    @cartesian.setter
    def cartesian(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def cartesianMax(self) -> float:
        ...
    @cartesianMax.setter
    def cartesianMax(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def improperDihedrals(self) -> float:
        ...
    @improperDihedrals.setter
    def improperDihedrals(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def improperDihedralsMax(self) -> float:
        ...
    @improperDihedralsMax.setter
    def improperDihedralsMax(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def matchResiduals(self) -> list[float]:
        ...
    @matchResiduals.setter
    def matchResiduals(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def properDihedrals(self) -> float:
        ...
    @properDihedrals.setter
    def properDihedrals(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def properDihedralsMax(self) -> float:
        ...
    @properDihedralsMax.setter
    def properDihedralsMax(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
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
class SystemTopology_NEW:
    def __init__(self) -> None:
        ...
    @property
    def a_coef(self) -> list[float]:
        """
        Lennard-Jones A coefficients (repulsive term) for each type pair, stored as a flat upper-triangular matrix in kJ/mol*nm^12.
        """
    @a_coef.setter
    def a_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_begin(self) -> VectorInt:
        """
        Begin index into angles arrays for each molecule.
        """
    @angles_begin.setter
    def angles_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def angles_end(self) -> VectorInt:
        """
        End index (exclusive) into angles arrays for each molecule.
        """
    @angles_end.setter
    def angles_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def angles_equilibrium(self) -> list[float]:
        """
        Equilibrium bond angle in radians [rad].
        """
    @angles_equilibrium.setter
    def angles_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_i(self) -> VectorInt:
        """
        Global index of angle atom 1 (outer).
        """
    @angles_i.setter
    def angles_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def angles_j(self) -> VectorInt:
        """
        Global index of angle atom 2 (central).
        """
    @angles_j.setter
    def angles_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def angles_k(self) -> VectorInt:
        """
        Global index of angle atom 3 (outer).
        """
    @angles_k.setter
    def angles_k(self, arg0: VectorInt) -> None:
        ...
    @property
    def angles_stiffness(self) -> list[float]:
        """
        Harmonic angle force constant in kilojoules per mole per radian squared [kJ/mol/rad^2].
        """
    @angles_stiffness.setter
    def angles_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_begin(self) -> VectorInt:
        """
        Begin index into atoms arrays for each molecule. Length equals num_molecules.
        """
    @atoms_begin.setter
    def atoms_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_charge(self) -> list[float]:
        """
        Partial charge in units of the proton charge [e].
        """
    @atoms_charge.setter
    def atoms_charge(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_charged_atom_type_index(self) -> VectorInt:
        """
        Index into the nonbonded parameter table for this atom's charged type.
        """
    @atoms_charged_atom_type_index.setter
    def atoms_charged_atom_type_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_charged_type_names(self) -> list[str]:
        """
        Charged atom type name for each atom.
        """
    @atoms_charged_type_names.setter
    def atoms_charged_type_names(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_class_index(self) -> VectorInt:
        """
        Index into the nonbonded parameter table for this atom's class.
        """
    @atoms_class_index.setter
    def atoms_class_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_class_names(self) -> list[str]:
        """
        Force-field atom class name for each atom.
        """
    @atoms_class_names.setter
    def atoms_class_names(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_compound_atom_index(self) -> VectorInt:
        """
        Index of the atom within its compound (residue) for each atom.
        """
    @atoms_compound_atom_index.setter
    def atoms_compound_atom_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_element_name(self) -> list[str]:
        """
        Full element name string (e.g. "Carbon").
        """
    @atoms_element_name.setter
    def atoms_element_name(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_element_symbol(self) -> list[str]:
        """
        Element symbol string (e.g. "C").
        """
    @atoms_element_symbol.setter
    def atoms_element_symbol(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_end(self) -> VectorInt:
        """
        End index (exclusive) into atoms arrays for each molecule.
        """
    @atoms_end.setter
    def atoms_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_epsilon(self) -> list[float]:
        """
        Lennard-Jones epsilon (well depth) in kilojoules per mole [kJ/mol].
        """
    @atoms_epsilon.setter
    def atoms_epsilon(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_mass(self) -> list[float]:
        """
        Mass in daltons [Da].
        """
    @atoms_mass.setter
    def atoms_mass(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_nonbonded_index(self) -> VectorInt:
        """
        Index into the nonbonded parameter table for this atom.
        """
    @atoms_nonbonded_index.setter
    def atoms_nonbonded_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_radius(self) -> list[float]:
        """
        GBSA implicit-solvent radius in nanometers [nm].
        """
    @atoms_radius.setter
    def atoms_radius(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_root_index(self) -> VectorInt:
        """
        Per-molecule index of the root atom for Z-matrix tree traversal.
        """
    @atoms_root_index.setter
    def atoms_root_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def atoms_screen(self) -> list[float]:
        """
        OBC screening factor (dimensionless).
        """
    @atoms_screen.setter
    def atoms_screen(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_sigma(self) -> list[float]:
        """
        Lennard-Jones sigma (van der Waals radius) in nanometers [nm].
        """
    @atoms_sigma.setter
    def atoms_sigma(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_unique_name(self) -> list[str]:
        """
        Unique atom name string (e.g. "ALA1_CA_3").
        """
    @atoms_unique_name.setter
    def atoms_unique_name(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_x(self) -> list[float]:
        """
        x coordinate of the reference structure in nanometers [nm].
        """
    @atoms_x.setter
    def atoms_x(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_y(self) -> list[float]:
        """
        y coordinate of the reference structure in nanometers [nm].
        """
    @atoms_y.setter
    def atoms_y(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_z(self) -> list[float]:
        """
        z coordinate of the reference structure in nanometers [nm].
        """
    @atoms_z.setter
    def atoms_z(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def b_coef(self) -> list[float]:
        """
        Lennard-Jones B coefficients (attractive term) for each type pair, stored as a flat upper-triangular matrix in kJ/mol*nm^6.
        """
    @b_coef.setter
    def b_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_begin(self) -> VectorInt:
        """
        Begin index into bonds arrays for each molecule.
        """
    @bonds_begin.setter
    def bonds_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def bonds_end(self) -> VectorInt:
        """
        End index (exclusive) into bonds arrays for each molecule.
        """
    @bonds_end.setter
    def bonds_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def bonds_equilibrium(self) -> list[float]:
        """
        Equilibrium bond length in nanometers [nm].
        """
    @bonds_equilibrium.setter
    def bonds_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_i(self) -> VectorInt:
        """
        BFS index of bond endpoint atom 1.
        """
    @bonds_i.setter
    def bonds_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def bonds_j(self) -> VectorInt:
        """
        BFS index of bond endpoint atom 2.
        """
    @bonds_j.setter
    def bonds_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def bonds_molecule_index(self) -> VectorInt:
        """
        Index of the molecule this bond belongs to.
        """
    @bonds_molecule_index.setter
    def bonds_molecule_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def bonds_ring_closing(self) -> list[bool]:
        """
        True if this bond closes a ring (i.e. is not a tree bond).
        """
    @bonds_ring_closing.setter
    def bonds_ring_closing(self, arg0: collections.abc.Sequence[bool]) -> None:
        ...
    @property
    def bonds_stiffness(self) -> list[float]:
        """
        Harmonic bond force constant in kilojoules per mole per nanometer squared [kJ/mol/nm^2].
        """
    @bonds_stiffness.setter
    def bonds_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_grid_energy(self) -> list[float]:
        """
        Flattened CMAP energy grid in row-major order in kilojoules per mole [kJ/mol].
        """
    @cmap_grid_energy.setter
    def cmap_grid_energy(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_grid_size(self) -> int:
        """
        Dimension of each CMAP grid (grid has cmap_grid_size x cmap_grid_size points).
        """
    @cmap_grid_size.setter
    def cmap_grid_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def cmap_torsion_a1(self) -> VectorInt:
        """
        Global index of torsion A atom 1.
        """
    @cmap_torsion_a1.setter
    def cmap_torsion_a1(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_a2(self) -> VectorInt:
        """
        Global index of torsion A atom 2.
        """
    @cmap_torsion_a2.setter
    def cmap_torsion_a2(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_a3(self) -> VectorInt:
        """
        Global index of torsion A atom 3.
        """
    @cmap_torsion_a3.setter
    def cmap_torsion_a3(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_a4(self) -> VectorInt:
        """
        Global index of torsion A atom 4.
        """
    @cmap_torsion_a4.setter
    def cmap_torsion_a4(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_b1(self) -> VectorInt:
        """
        Global index of torsion B atom 1.
        """
    @cmap_torsion_b1.setter
    def cmap_torsion_b1(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_b2(self) -> VectorInt:
        """
        Global index of torsion B atom 2.
        """
    @cmap_torsion_b2.setter
    def cmap_torsion_b2(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_b3(self) -> VectorInt:
        """
        Global index of torsion B atom 3.
        """
    @cmap_torsion_b3.setter
    def cmap_torsion_b3(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_b4(self) -> VectorInt:
        """
        Global index of torsion B atom 4.
        """
    @cmap_torsion_b4.setter
    def cmap_torsion_b4(self, arg0: VectorInt) -> None:
        ...
    @property
    def cmap_torsion_map_index(self) -> VectorInt:
        """
        Index into the CMAP grid table for each torsion pair.
        """
    @cmap_torsion_map_index.setter
    def cmap_torsion_map_index(self, arg0: VectorInt) -> None:
        ...
    @property
    def collision_frequency(self) -> float:
        """
        Langevin collision frequency (friction coefficient) in inverse picoseconds [ps^-1].
        """
    @collision_frequency.setter
    def collision_frequency(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def exclusion_begin(self) -> VectorInt:
        """
        Begin index into exclusions arrays for each molecule.
        """
    @exclusion_begin.setter
    def exclusion_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def exclusion_end(self) -> VectorInt:
        """
        End index (exclusive) into exclusions arrays for each molecule.
        """
    @exclusion_end.setter
    def exclusion_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def exclusion_i(self) -> VectorInt:
        """
        Global index of atom 1 of the excluded pair.
        """
    @exclusion_i.setter
    def exclusion_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def exclusion_j(self) -> VectorInt:
        """
        Global index of atom 2 of the excluded pair.
        """
    @exclusion_j.setter
    def exclusion_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def gbsa_solute_dielectric(self) -> float:
        """
        Relative dielectric constant of the solute interior (dimensionless). Default 1.0.
        """
    @gbsa_solute_dielectric.setter
    def gbsa_solute_dielectric(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def gbsa_solvent_dielectric(self) -> float:
        """
        Relative dielectric constant of the solvent (dimensionless). Default 78.5 (water at 298 K).
        """
    @gbsa_solvent_dielectric.setter
    def gbsa_solvent_dielectric(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def harmonic_torsions_begin(self) -> VectorInt:
        """
        Begin index into harmonic torsions arrays for each molecule.
        """
    @harmonic_torsions_begin.setter
    def harmonic_torsions_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_end(self) -> VectorInt:
        """
        End index (exclusive) into harmonic torsions arrays for each molecule.
        """
    @harmonic_torsions_end.setter
    def harmonic_torsions_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_i(self) -> VectorInt:
        """
        Global index of torsion atom 1.
        """
    @harmonic_torsions_i.setter
    def harmonic_torsions_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_j(self) -> VectorInt:
        """
        Global index of torsion atom 2.
        """
    @harmonic_torsions_j.setter
    def harmonic_torsions_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_k(self) -> VectorInt:
        """
        Global index of torsion atom 3.
        """
    @harmonic_torsions_k.setter
    def harmonic_torsions_k(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_l(self) -> VectorInt:
        """
        Global index of torsion atom 4.
        """
    @harmonic_torsions_l.setter
    def harmonic_torsions_l(self, arg0: VectorInt) -> None:
        ...
    @property
    def harmonic_torsions_phase(self) -> list[float]:
        """
        Equilibrium dihedral angle in radians [rad].
        """
    @harmonic_torsions_phase.setter
    def harmonic_torsions_phase(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_stiffness(self) -> list[float]:
        """
        Force constant in kilojoules per mole [kJ/mol].
        """
    @harmonic_torsions_stiffness.setter
    def harmonic_torsions_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def has_nbfix(self) -> bool:
        """
        True if NBfix pairwise corrections are present.
        """
    @has_nbfix.setter
    def has_nbfix(self, arg0: bool) -> None:
        ...
    @property
    def nonbonded_cutoff(self) -> float:
        """
        Distance cutoff for nonbonded interactions in nanometers [nm].
        """
    @nonbonded_cutoff.setter
    def nonbonded_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def nonbonded_method(self) -> NonbondedMethod:
        """
        Algorithm used to evaluate nonbonded interactions.
        """
    @nonbonded_method.setter
    def nonbonded_method(self, arg0: NonbondedMethod) -> None:
        ...
    @property
    def num_angles(self) -> int:
        """
        Total number of angles.
        """
    @num_angles.setter
    def num_angles(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_atoms(self) -> int:
        """
        Total number of atoms.
        """
    @num_atoms.setter
    def num_atoms(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_bonds(self) -> int:
        """
        Total number of bonds (tree bonds plus ring-closing bonds).
        """
    @num_bonds.setter
    def num_bonds(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_exclusions(self) -> int:
        """
        Total number of non-bonded exclusion pairs.
        """
    @num_exclusions.setter
    def num_exclusions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_harmonic_torsions(self) -> int:
        """
        Total number of harmonic torsion terms.
        """
    @num_harmonic_torsions.setter
    def num_harmonic_torsions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_molecules(self) -> int:
        """
        Total number of molecules.
        """
    @num_molecules.setter
    def num_molecules(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_nb_types(self) -> int:
        """
        Number of distinct nonbonded atom types.
        """
    @num_nb_types.setter
    def num_nb_types(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_periodic_torsions(self) -> int:
        """
        Total number of periodic torsion terms.
        """
    @num_periodic_torsions.setter
    def num_periodic_torsions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_scaling14(self) -> int:
        """
        Total number of 1-4 pair scaling terms.
        """
    @num_scaling14.setter
    def num_scaling14(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_urey_bradley(self) -> int:
        """
        Total number of Urey-Bradley 1-3 terms.
        """
    @num_urey_bradley.setter
    def num_urey_bradley(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_z_matrix_rows(self) -> int:
        """
        Number of Z-matrix rows (equals number of atoms).
        """
    @num_z_matrix_rows.setter
    def num_z_matrix_rows(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def periodic_torsions_begin(self) -> VectorInt:
        """
        Begin index into periodic torsions arrays for each molecule.
        """
    @periodic_torsions_begin.setter
    def periodic_torsions_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_end(self) -> VectorInt:
        """
        End index (exclusive) into periodic torsions arrays for each molecule.
        """
    @periodic_torsions_end.setter
    def periodic_torsions_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_i(self) -> VectorInt:
        """
        Global index of torsion atom 1.
        """
    @periodic_torsions_i.setter
    def periodic_torsions_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_improper(self) -> list[bool]:
        """
        True if the torsion is an improper (out-of-plane) term.
        """
    @periodic_torsions_improper.setter
    def periodic_torsions_improper(self, arg0: collections.abc.Sequence[bool]) -> None:
        ...
    @property
    def periodic_torsions_j(self) -> VectorInt:
        """
        Global index of torsion atom 2.
        """
    @periodic_torsions_j.setter
    def periodic_torsions_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_k(self) -> VectorInt:
        """
        Global index of torsion atom 3.
        """
    @periodic_torsions_k.setter
    def periodic_torsions_k(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_l(self) -> VectorInt:
        """
        Global index of torsion atom 4.
        """
    @periodic_torsions_l.setter
    def periodic_torsions_l(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_n(self) -> VectorInt:
        """
        Torsion periodicity (integer, dimensionless).
        """
    @periodic_torsions_n.setter
    def periodic_torsions_n(self, arg0: VectorInt) -> None:
        ...
    @property
    def periodic_torsions_phase(self) -> list[float]:
        """
        Phase offset in radians [rad].
        """
    @periodic_torsions_phase.setter
    def periodic_torsions_phase(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_stiffness(self) -> list[float]:
        """
        Force constant in kilojoules per mole [kJ/mol].
        """
    @periodic_torsions_stiffness.setter
    def periodic_torsions_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_begin(self) -> VectorInt:
        """
        Begin index into 1-4 scaling arrays for each molecule.
        """
    @scaling14_begin.setter
    def scaling14_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def scaling14_charge_product(self) -> list[float]:
        """
        q1 times q4 pre-scaled by the 1-4 electrostatic factor, in units of proton charge squared [e^2].
        """
    @scaling14_charge_product.setter
    def scaling14_charge_product(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_end(self) -> VectorInt:
        """
        End index (exclusive) into 1-4 scaling arrays for each molecule.
        """
    @scaling14_end.setter
    def scaling14_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def scaling14_epsilon(self) -> list[float]:
        """
        Combined Lennard-Jones well depth in kilojoules per mole [kJ/mol].
        """
    @scaling14_epsilon.setter
    def scaling14_epsilon(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_i(self) -> VectorInt:
        """
        Global index of atom 1 (first atom of dihedral i-j-k-l).
        """
    @scaling14_i.setter
    def scaling14_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def scaling14_l(self) -> VectorInt:
        """
        Global index of atom 4 (last atom of dihedral i-j-k-l).
        """
    @scaling14_l.setter
    def scaling14_l(self, arg0: VectorInt) -> None:
        ...
    @property
    def scaling14_sigma(self) -> list[float]:
        """
        Combined Lennard-Jones radius in nanometers [nm].
        """
    @scaling14_sigma.setter
    def scaling14_sigma(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def seed(self) -> int:
        """
        Random number seed for the integrator. 0 means system-generated.
        """
    @seed.setter
    def seed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def thermostat_temperature(self) -> float:
        """
        Target temperature for the Langevin thermostat in Kelvin [K].
        """
    @thermostat_temperature.setter
    def thermostat_temperature(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def urey_bradley_begin(self) -> VectorInt:
        """
        Begin index into Urey-Bradley arrays for each molecule.
        """
    @urey_bradley_begin.setter
    def urey_bradley_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def urey_bradley_end(self) -> VectorInt:
        """
        End index (exclusive) into Urey-Bradley arrays for each molecule.
        """
    @urey_bradley_end.setter
    def urey_bradley_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def urey_bradley_equilibrium(self) -> list[float]:
        """
        Nominal 1-3 distance in nanometers [nm].
        """
    @urey_bradley_equilibrium.setter
    def urey_bradley_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_i(self) -> VectorInt:
        """
        Global index of atom 1 (outer atom of angle i-j-k).
        """
    @urey_bradley_i.setter
    def urey_bradley_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def urey_bradley_k(self) -> VectorInt:
        """
        Global index of atom 3 (outer atom of angle i-j-k).
        """
    @urey_bradley_k.setter
    def urey_bradley_k(self, arg0: VectorInt) -> None:
        ...
    @property
    def urey_bradley_stiffness(self) -> list[float]:
        """
        Harmonic force constant in kilojoules per mole per nanometer squared [kJ/mol/nm^2].
        """
    @urey_bradley_stiffness.setter
    def urey_bradley_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def use_gbsa_obc2(self) -> bool:
        """
        True if the GBSA OBC2 implicit solvent model is active.
        """
    @use_gbsa_obc2.setter
    def use_gbsa_obc2(self, arg0: bool) -> None:
        ...
    @property
    def z_matrix_begin(self) -> VectorInt:
        """
        Begin index into Z-matrix arrays for each molecule.
        """
    @z_matrix_begin.setter
    def z_matrix_begin(self, arg0: VectorInt) -> None:
        ...
    @property
    def z_matrix_end(self) -> VectorInt:
        """
        End index (exclusive) into Z-matrix arrays for each molecule.
        """
    @z_matrix_end.setter
    def z_matrix_end(self, arg0: VectorInt) -> None:
        ...
    @property
    def z_matrix_i(self) -> VectorInt:
        """
        Global atom index of the atom placed at row r.
        """
    @z_matrix_i.setter
    def z_matrix_i(self, arg0: VectorInt) -> None:
        ...
    @property
    def z_matrix_j(self) -> VectorInt:
        """
        Bond-length reference atom at row r. Row 0 holds sentinel -1 (root has no bond reference).
        """
    @z_matrix_j.setter
    def z_matrix_j(self, arg0: VectorInt) -> None:
        ...
    @property
    def z_matrix_k(self) -> VectorInt:
        """
        Bond-angle reference atom at row r. Rows 0-1 hold sentinel -1.
        """
    @z_matrix_k.setter
    def z_matrix_k(self, arg0: VectorInt) -> None:
        ...
    @property
    def z_matrix_l(self) -> VectorInt:
        """
        Dihedral reference atom at row r. Rows 0-2 hold sentinel -1.
        """
    @z_matrix_l.setter
    def z_matrix_l(self, arg0: VectorInt) -> None:
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
    def get_coordinate_transfer_errors(self) -> list[CoordinateTransferError]:
        """
        Get the coordinate transfer errors for all samplers in the world.
        """
    def has_rigid_body_violations(self, timeStep: typing.SupportsFloat | typing.SupportsIndex, numSteps: typing.SupportsInt | typing.SupportsIndex) -> bool:
        """
        Checks for rigid body violations.
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
