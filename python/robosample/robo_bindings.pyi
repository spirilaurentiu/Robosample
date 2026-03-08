"""
Robosample bindings
"""
from __future__ import annotations
import collections.abc
import typing
__all__: list[str] = ['AcceptRejectMode', 'BondFlexibility', 'BondMobility', 'CMAPGrid', 'CMAPTorsion', 'Context', 'Exclusion', 'IntegratorType', 'NonbondedMethod', 'RoboAngle', 'RoboAtom', 'RoboAtomConnectivity', 'RoboAtomElement', 'RoboAtomIdentity', 'RoboAtomPhysics', 'RoboBond', 'RoboHarmonicImproperTorsion', 'RoboPeriodicTorsion', 'RoboPeriodicTorsionTerm', 'RootMobility', 'RunType', 'SamplerName', 'Scaling14', 'ThermostatName', 'TopologyRange', 'TopologyRangeType', 'UreyBradley', 'World']
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
    def __init__(self) -> None:
        ...
    @property
    def a1(self) -> int:
        ...
    @a1.setter
    def a1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a2(self) -> int:
        ...
    @a2.setter
    def a2(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a3(self) -> int:
        ...
    @a3.setter
    def a3(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a4(self) -> int:
        ...
    @a4.setter
    def a4(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def b1(self) -> int:
        ...
    @b1.setter
    def b1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def b2(self) -> int:
        ...
    @b2.setter
    def b2(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def b3(self) -> int:
        ...
    @b3.setter
    def b3(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def b4(self) -> int:
        ...
    @b4.setter
    def b4(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def mapIndex(self) -> int:
        ...
    @mapIndex.setter
    def mapIndex(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
class Context:
    def Initialize(self) -> None:
        """
        Initializes the context after all worlds and replicas have been set.
        """
    def RunREX(self, arg0: typing.SupportsInt | typing.SupportsIndex, arg1: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Run replica exchange.
        """
    def __init__(self, arg0: str, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: typing.SupportsInt | typing.SupportsIndex, arg3: typing.SupportsInt | typing.SupportsIndex, arg4: RunType, arg5: typing.SupportsInt | typing.SupportsIndex, arg6: typing.SupportsInt | typing.SupportsIndex, arg7: bool) -> None:
        ...
    def addReplica(self) -> None:
        """
        Add an empty replica to the context.
        """
    def addThermodynamicState(self, arg0: typing.SupportsFloat | typing.SupportsIndex, arg1: collections.abc.Sequence[AcceptRejectMode], arg2: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg3: collections.abc.Sequence[str], arg4: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg5: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg6: collections.abc.Sequence[IntegratorType], arg7: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg8: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], arg9: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        """
        Add an empty themodynamic state to the context.
        """
    def addWorld(self, arg0: bool, arg1: typing.SupportsInt | typing.SupportsIndex, arg2: RootMobility, arg3: collections.abc.Sequence[collections.abc.Sequence[BondFlexibility]], arg4: bool, arg5: bool, arg6: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Add an empty world.
        """
    def getWorld(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> World:
        ...
    def getWorlds(self) -> list[World]:
        ...
    def initialize_openmm(self, arg0: collections.abc.Sequence[RoboAtom], arg1: collections.abc.Sequence[RoboBond], arg2: collections.abc.Sequence[RoboAngle], arg3: collections.abc.Sequence[RoboPeriodicTorsion], arg4: collections.abc.Sequence[RoboHarmonicImproperTorsion], arg5: collections.abc.Sequence[CMAPGrid], arg6: collections.abc.Sequence[CMAPTorsion], arg7: collections.abc.Sequence[UreyBradley], arg8: bool, arg9: typing.SupportsInt | typing.SupportsIndex, arg10: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], arg11: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], arg12: collections.abc.Sequence[Exclusion], arg13: collections.abc.Sequence[Scaling14]) -> dict[str, float]:
        """
        Load an OpenMM system from components.
        """
    def loadAmberSystem(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], arg1: collections.abc.Sequence[RoboAtom], arg2: collections.abc.Sequence[RoboBond], arg3: collections.abc.Sequence[RoboAngle], arg4: collections.abc.Sequence[RoboPeriodicTorsion], arg5: collections.abc.Sequence[RoboHarmonicImproperTorsion], arg6: collections.abc.Sequence[TopologyRange]) -> None:
        """
        Load an AMBER system.
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
    def setPrintFreq(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Set the print frequency.
        """
    def setVerbose(self, arg0: bool) -> None:
        """
        Control if you want extraneous output to cout.
        """
class Exclusion:
    @typing.overload
    def __init__(self) -> None:
        ...
    @typing.overload
    def __init__(self, a1: typing.SupportsInt | typing.SupportsIndex, a2: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a1(self) -> int:
        ...
    @a1.setter
    def a1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a2(self) -> int:
        ...
    @a2.setter
    def a2(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
class RoboAngle:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(3)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, stiffness_in_kj_per_rad_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_angle_in_deg: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[...], "FixedSize(3)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[...], "FixedSize(3)"]) -> None:
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
    def __init__(self, identity: RoboAtomIdentity, element_info: RoboAtomElement, physics: RoboAtomPhysics, connectivity: RoboAtomConnectivity, position: typing.Annotated[collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex], "FixedSize(3)"]) -> None:
        ...
    @property
    def position(self) -> ...:
        ...
    @position.setter
    def position(*args, **kwargs):
        """
        (self: robo_bindings.RoboAtom, arg0: SimTK::Vec<3, double, 1>) -> None
        """
class RoboAtomConnectivity:
    root: bool
    def __init__(self, neighbors_global_indices: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], root: bool) -> None:
        ...
    @property
    def neighbors_global_indices(self) -> list[int]:
        ...
    @neighbors_global_indices.setter
    def neighbors_global_indices(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
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
    atom_class_index: ...
    atom_class_name: str
    charged_atom_type_index: ...
    charged_atom_type_name: str
    compound_atom_index: ...
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
    ring_closing: bool
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(2)"], stiffness_in_kj_per_nm_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_length_in_nm: typing.SupportsFloat | typing.SupportsIndex, molecule_index: typing.SupportsInt | typing.SupportsIndex, ring_closing: bool) -> None:
        ...
    @property
    def compound_atom_indices(self) -> typing.Annotated[list[...], "FixedSize(2)"]:
        ...
    @compound_atom_indices.setter
    def compound_atom_indices(self, arg0: typing.Annotated[collections.abc.Sequence[...], "FixedSize(2)"]) -> None:
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
    def stiffness_in_kj_per_nm_sq(self) -> float:
        ...
    @stiffness_in_kj_per_nm_sq.setter
    def stiffness_in_kj_per_nm_sq(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboHarmonicImproperTorsion:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, stiffness_in_kj_per_rad_sq: typing.SupportsFloat | typing.SupportsIndex, nominal_angle_in_rad: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
class RoboPeriodicTorsion:
    def __init__(self, global_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], compound_atom_indices: typing.Annotated[collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex], "FixedSize(4)"], molecule_index: typing.SupportsInt | typing.SupportsIndex, improper: bool, terms: collections.abc.Sequence[RoboPeriodicTorsionTerm]) -> None:
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
    def __init__(self, a1: typing.SupportsInt | typing.SupportsIndex, a4: typing.SupportsInt | typing.SupportsIndex, charge_product: typing.SupportsFloat | typing.SupportsIndex, epsilon: typing.SupportsFloat | typing.SupportsIndex, sigma: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def a1(self) -> int:
        ...
    @a1.setter
    def a1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a4(self) -> int:
        ...
    @a4.setter
    def a4(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
    def __init__(self, startCounts: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    def close(self, endCounts: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
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
    def __init__(self) -> None:
        ...
    @property
    def a1(self) -> int:
        ...
    @a1.setter
    def a1(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def a3(self) -> int:
        ...
    @a3.setter
    def a3(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
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
class World:
    def addSampler(self, arg0: SamplerName, arg1: IntegratorType, arg2: ThermostatName, arg3: bool) -> bool:
        """
        Add a sampler to the world.
        """
    def getAcceptanceRMSD(self) -> list[tuple[bool, float]]:
        """
        Get RMSD values for accepted/rejected moves during testing.
        """
    def getCumulativeAngleDisplacements(self) -> list[float]:
        """
        Get cumulative angle displacements.
        """
    def getCumulativeBondDisplacements(self) -> list[float]:
        """
        Get cumulative bond displacements.
        """
    def getCumulativeCartesianDisplacements(self) -> list[float]:
        """
        Get cumulative Cartesian displacements.
        """
    def getCumulativeTorsionDisplacements(self) -> list[float]:
        """
        Get cumulative torsion displacements.
        """
    def getMatchAtomTargetLocationsResiduals(self) -> list[list[float]]:
        """
        Get residuals from matchAtomTargetLocations during testing.
        """
    def getOpenMMCumulativeAngleDisplacements(self) -> list[float]:
        """
        Get cumulative angle displacements between OpenMM and SimTK during testing.
        """
    def getOpenMMCumulativeBondDisplacements(self) -> list[float]:
        """
        Get cumulative bond displacements between OpenMM and SimTK during testing.
        """
    def getOpenMMCumulativeCartesianDisplacements(self) -> list[float]:
        """
        Get cumulative Cartesian displacements between OpenMM and SimTK during testing.
        """
    def getOpenMMCumulativeTorsionDisplacements(self) -> list[float]:
        """
        Get cumulative torsion displacements between OpenMM and SimTK during testing.
        """
    def getRigidBodyAngleDriftInRad(self) -> list[float]:
        """
        Get rigid body angle drift in rad during testing.
        """
    def getRigidBodyBondRMSDInNm(self) -> list[float]:
        """
        Get rigid body bond RMSD in nm during testing.
        """
    def getRigidBodyImproperTorsionDriftInRad(self) -> list[float]:
        """
        Get rigid body improper torsion drift in rad during testing.
        """
    def getRigidBodyProperTorsionDriftInRad(self) -> list[float]:
        """
        Get rigid body proper torsion drift in rad during testing.
        """
    def getRollFlexibilities(self) -> bool:
        """
        Get rolling sequence.
        """
    def setFlexibilities(self, arg0: collections.abc.Sequence[BondFlexibility]) -> None:
        """
        Set the flexibilities of the bonds.
        """
    def setRollFlexibilities(self, arg0: bool) -> None:
        """
        Set rolling sequence.
        """
