"""
Robosample C++ bindings (robo_bindings)
"""
from __future__ import annotations
import collections.abc
import typing
__all__: list[str] = ['AcceptRejectMode', 'BondMobility', 'Context', 'ForceGroupEnergy', 'JointType', 'MoveType', 'NonbondedMethod', 'RootMobility', 'Selection', 'SystemTopology', 'World']
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
class BondMobility:
    """
    Members:
    
      Rigid
    
      Torsion
    
      Free
    
      Ball
    
      Pin
    
      Slider
    
      Cylinder
    
      BendStretch
    """
    Ball: typing.ClassVar[BondMobility]  # value = <BondMobility.Ball: 3>
    BendStretch: typing.ClassVar[BondMobility]  # value = <BondMobility.BendStretch: 7>
    Cylinder: typing.ClassVar[BondMobility]  # value = <BondMobility.Cylinder: 6>
    Free: typing.ClassVar[BondMobility]  # value = <BondMobility.Free: 2>
    Pin: typing.ClassVar[BondMobility]  # value = <BondMobility.Pin: 4>
    Rigid: typing.ClassVar[BondMobility]  # value = <BondMobility.Rigid: 0>
    Slider: typing.ClassVar[BondMobility]  # value = <BondMobility.Slider: 5>
    Torsion: typing.ClassVar[BondMobility]  # value = <BondMobility.Torsion: 1>
    __members__: typing.ClassVar[dict[str, BondMobility]]  # value = {'Rigid': <BondMobility.Rigid: 0>, 'Torsion': <BondMobility.Torsion: 1>, 'Free': <BondMobility.Free: 2>, 'Ball': <BondMobility.Ball: 3>, 'Pin': <BondMobility.Pin: 4>, 'Slider': <BondMobility.Slider: 5>, 'Cylinder': <BondMobility.Cylinder: 6>, 'BendStretch': <BondMobility.BendStretch: 7>}
    @typing.overload
    def __eq__(self, other: BondMobility) -> bool:
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
class Context:
    system_topology: SystemTopology
    def __init__(self, base_name: str, seed: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    def add_cartesian_world(self) -> World:
        """
        Add a Cartesian (OpenMM-MD) world; returns it for .add_sampler(...).
        """
    def add_docking_world(self, ligand_molecule_indices: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> World:
        """
        Add a rigid-body docking world: the listed ligand molecules get Free roots, every other molecule is welded/rigid. Chain .add_sampler(sphere_radius=...).
        """
    def add_robotic_world(self, selection: Selection) -> World:
        """
        Add an internal-coordinate (torsional) world for the given selection.
        """
    def add_torsional_world(self, selection: Selection) -> World:
        """
        Alias of add_robotic_world.
        """
    def build_flexibilities(self, bonds: collections.abc.Sequence[tuple[typing.SupportsInt | typing.SupportsIndex, typing.SupportsInt | typing.SupportsIndex]] | None, mobility: BondMobility, flag: bool) -> Selection:
        """
        Build a per-bond mobility selection (bonds=None => all eligible).
        """
    def calc_openmm_potential_energy(self) -> float:
        """
        Set the reference coordinates and return the OpenMM potential energy [kJ/mol].
        """
    def calc_openmm_potential_energy_by_group(self) -> tuple[float, list[ForceGroupEnergy]]:
        """
        Compute potential energy by OpenMM force group, returning (group, name, energy) tuples.
        """
    def initialize(self, temperatures: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex] = []) -> None:
        """
        Build OpenMM, set the replica temperature ladder, seed coordinates.
        """
    def initialize_openmm(self) -> bool:
        """
        Build the single OpenMM System/Context from system_topology.
        """
    def run_rex(self, equil_rounds: typing.SupportsInt | typing.SupportsIndex, prod_rounds: typing.SupportsInt | typing.SupportsIndex, write_freq: typing.SupportsInt | typing.SupportsIndex, verbose: bool) -> None:
        """
        Run replica exchange: Gibbs sweep over worlds + adjacent swaps.
        """
    def set_enforce_periodic_box(self, enabled: bool) -> None:
        """
        Whether OpenMM wraps coordinates into the primary box when state is pulled back. MUST stay False (the default) under explicit solvent so the robot engine receives whole molecules; energies/forces are unaffected (minimum image is always applied internally).
        """
    def set_mts(self, enabled: bool, inner_substeps: typing.SupportsInt | typing.SupportsIndex = 4) -> None:
        """
        Enable r-RESPA multiple-timestep OpenMM MD (Cartesian world): slow forces once per outer step, fast bonded forces inner_substeps times. Call before initialize().
        """
    def set_root_mobility(self, molecule_index: typing.SupportsInt | typing.SupportsIndex, mobility: RootMobility) -> None:
        """
        Override a molecule's root attachment to Ground.
        """
    def set_separate_force_groups(self, enabled: bool) -> None:
        """
        Enable/disable separate OpenMM force groups for each Force.
        """
class ForceGroupEnergy:
    @property
    def energy(self) -> float:
        ...
    @property
    def group(self) -> int:
        ...
    @property
    def name(self) -> str:
        ...
class JointType:
    """
    Members:
    
      Weld
    
      Pin
    
      Slider
    
      Cylinder
    
      BendStretch
    
      Translation
    
      Ball
    
      SphericalCoords
    
      FreeLine
    
      Free
    """
    Ball: typing.ClassVar[JointType]  # value = <JointType.Ball: 6>
    BendStretch: typing.ClassVar[JointType]  # value = <JointType.BendStretch: 4>
    Cylinder: typing.ClassVar[JointType]  # value = <JointType.Cylinder: 3>
    Free: typing.ClassVar[JointType]  # value = <JointType.Free: 9>
    FreeLine: typing.ClassVar[JointType]  # value = <JointType.FreeLine: 8>
    Pin: typing.ClassVar[JointType]  # value = <JointType.Pin: 1>
    Slider: typing.ClassVar[JointType]  # value = <JointType.Slider: 2>
    SphericalCoords: typing.ClassVar[JointType]  # value = <JointType.SphericalCoords: 7>
    Translation: typing.ClassVar[JointType]  # value = <JointType.Translation: 5>
    Weld: typing.ClassVar[JointType]  # value = <JointType.Weld: 0>
    __members__: typing.ClassVar[dict[str, JointType]]  # value = {'Weld': <JointType.Weld: 0>, 'Pin': <JointType.Pin: 1>, 'Slider': <JointType.Slider: 2>, 'Cylinder': <JointType.Cylinder: 3>, 'BendStretch': <JointType.BendStretch: 4>, 'Translation': <JointType.Translation: 5>, 'Ball': <JointType.Ball: 6>, 'SphericalCoords': <JointType.SphericalCoords: 7>, 'FreeLine': <JointType.FreeLine: 8>, 'Free': <JointType.Free: 9>}
    @typing.overload
    def __eq__(self, other: JointType) -> bool:
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
    def __ne__(self, other: JointType) -> bool:
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
class MoveType:
    """
    Members:
    
      MdHmc
    
      RigidKick
    
      NcmcSwitch
    """
    MdHmc: typing.ClassVar[MoveType]  # value = <MoveType.MdHmc: 0>
    NcmcSwitch: typing.ClassVar[MoveType]  # value = <MoveType.NcmcSwitch: 2>
    RigidKick: typing.ClassVar[MoveType]  # value = <MoveType.RigidKick: 1>
    __members__: typing.ClassVar[dict[str, MoveType]]  # value = {'MdHmc': <MoveType.MdHmc: 0>, 'RigidKick': <MoveType.RigidKick: 1>, 'NcmcSwitch': <MoveType.NcmcSwitch: 2>}
    @typing.overload
    def __eq__(self, other: MoveType) -> bool:
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
    def __ne__(self, other: MoveType) -> bool:
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
    
      CutoffPeriodic
    
      Ewald
    
      PME
    """
    CutoffNonPeriodic: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.CutoffNonPeriodic: 1>
    CutoffPeriodic: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.CutoffPeriodic: 2>
    Ewald: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.Ewald: 3>
    NoCutoff: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.NoCutoff: 0>
    PME: typing.ClassVar[NonbondedMethod]  # value = <NonbondedMethod.PME: 4>
    __members__: typing.ClassVar[dict[str, NonbondedMethod]]  # value = {'NoCutoff': <NonbondedMethod.NoCutoff: 0>, 'CutoffNonPeriodic': <NonbondedMethod.CutoffNonPeriodic: 1>, 'CutoffPeriodic': <NonbondedMethod.CutoffPeriodic: 2>, 'Ewald': <NonbondedMethod.Ewald: 3>, 'PME': <NonbondedMethod.PME: 4>}
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
class Selection:
    def __init__(self) -> None:
        ...
class SystemTopology:
    has_nb_fix: bool
    nonbonded_method: NonbondedMethod
    use_gbsa_obc2: bool
    def __init__(self) -> None:
        ...
    @property
    def a_coef(self) -> list[float]:
        ...
    @a_coef.setter
    def a_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_begin(self) -> list[int]:
        ...
    @angles_begin.setter
    def angles_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_end(self) -> list[int]:
        ...
    @angles_end.setter
    def angles_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_equilibrium(self) -> list[float]:
        ...
    @angles_equilibrium.setter
    def angles_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_i(self) -> list[int]:
        ...
    @angles_i.setter
    def angles_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_j(self) -> list[int]:
        ...
    @angles_j.setter
    def angles_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_k(self) -> list[int]:
        ...
    @angles_k.setter
    def angles_k(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def angles_stiffness(self) -> list[float]:
        ...
    @angles_stiffness.setter
    def angles_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_atomic_number(self) -> list[int]:
        ...
    @atoms_atomic_number.setter
    def atoms_atomic_number(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_begin(self) -> list[int]:
        ...
    @atoms_begin.setter
    def atoms_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_charge(self) -> list[float]:
        ...
    @atoms_charge.setter
    def atoms_charge(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_element_name(self) -> list[str]:
        ...
    @atoms_element_name.setter
    def atoms_element_name(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_element_symbol(self) -> list[str]:
        ...
    @atoms_element_symbol.setter
    def atoms_element_symbol(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_end(self) -> list[int]:
        ...
    @atoms_end.setter
    def atoms_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_epsilon(self) -> list[float]:
        ...
    @atoms_epsilon.setter
    def atoms_epsilon(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_mass(self) -> list[float]:
        ...
    @atoms_mass.setter
    def atoms_mass(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_nonbonded_index(self) -> list[int]:
        ...
    @atoms_nonbonded_index.setter
    def atoms_nonbonded_index(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_num_bonds_involved(self) -> list[int]:
        ...
    @atoms_num_bonds_involved.setter
    def atoms_num_bonds_involved(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_prmtop_index(self) -> list[int]:
        ...
    @atoms_prmtop_index.setter
    def atoms_prmtop_index(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_radius(self) -> list[float]:
        ...
    @atoms_radius.setter
    def atoms_radius(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_root_index(self) -> list[int]:
        ...
    @atoms_root_index.setter
    def atoms_root_index(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_screen(self) -> list[float]:
        ...
    @atoms_screen.setter
    def atoms_screen(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_sigma(self) -> list[float]:
        ...
    @atoms_sigma.setter
    def atoms_sigma(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_unique_name(self) -> list[str]:
        ...
    @atoms_unique_name.setter
    def atoms_unique_name(self, arg0: collections.abc.Sequence[str]) -> None:
        ...
    @property
    def atoms_x(self) -> list[float]:
        ...
    @atoms_x.setter
    def atoms_x(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_y(self) -> list[float]:
        ...
    @atoms_y.setter
    def atoms_y(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def atoms_z(self) -> list[float]:
        ...
    @atoms_z.setter
    def atoms_z(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def b_coef(self) -> list[float]:
        ...
    @b_coef.setter
    def b_coef(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_begin(self) -> list[int]:
        ...
    @bonds_begin.setter
    def bonds_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_end(self) -> list[int]:
        ...
    @bonds_end.setter
    def bonds_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_equilibrium(self) -> list[float]:
        ...
    @bonds_equilibrium.setter
    def bonds_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_i(self) -> list[int]:
        ...
    @bonds_i.setter
    def bonds_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_j(self) -> list[int]:
        ...
    @bonds_j.setter
    def bonds_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_molecule_index(self) -> list[int]:
        ...
    @bonds_molecule_index.setter
    def bonds_molecule_index(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def bonds_ring_closing(self) -> list[bool]:
        ...
    @bonds_ring_closing.setter
    def bonds_ring_closing(self, arg0: collections.abc.Sequence[bool]) -> None:
        ...
    @property
    def bonds_stiffness(self) -> list[float]:
        ...
    @bonds_stiffness.setter
    def bonds_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def box_vectors(self) -> list[float]:
        ...
    @box_vectors.setter
    def box_vectors(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_grid_energy(self) -> list[float]:
        ...
    @cmap_grid_energy.setter
    def cmap_grid_energy(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_grid_size(self) -> int:
        ...
    @cmap_grid_size.setter
    def cmap_grid_size(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def cmap_torsion_a1(self) -> list[int]:
        ...
    @cmap_torsion_a1.setter
    def cmap_torsion_a1(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_a2(self) -> list[int]:
        ...
    @cmap_torsion_a2.setter
    def cmap_torsion_a2(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_a3(self) -> list[int]:
        ...
    @cmap_torsion_a3.setter
    def cmap_torsion_a3(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_a4(self) -> list[int]:
        ...
    @cmap_torsion_a4.setter
    def cmap_torsion_a4(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_b1(self) -> list[int]:
        ...
    @cmap_torsion_b1.setter
    def cmap_torsion_b1(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_b2(self) -> list[int]:
        ...
    @cmap_torsion_b2.setter
    def cmap_torsion_b2(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_b3(self) -> list[int]:
        ...
    @cmap_torsion_b3.setter
    def cmap_torsion_b3(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_b4(self) -> list[int]:
        ...
    @cmap_torsion_b4.setter
    def cmap_torsion_b4(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def cmap_torsion_map_index(self) -> list[int]:
        ...
    @cmap_torsion_map_index.setter
    def cmap_torsion_map_index(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def collision_frequency(self) -> float:
        ...
    @collision_frequency.setter
    def collision_frequency(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def ewald_error_tolerance(self) -> float:
        ...
    @ewald_error_tolerance.setter
    def ewald_error_tolerance(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def exclusion_begin(self) -> list[int]:
        ...
    @exclusion_begin.setter
    def exclusion_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def exclusion_end(self) -> list[int]:
        ...
    @exclusion_end.setter
    def exclusion_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def exclusion_i(self) -> list[int]:
        ...
    @exclusion_i.setter
    def exclusion_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def exclusion_j(self) -> list[int]:
        ...
    @exclusion_j.setter
    def exclusion_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
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
    def harmonic_torsions_begin(self) -> list[int]:
        ...
    @harmonic_torsions_begin.setter
    def harmonic_torsions_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_end(self) -> list[int]:
        ...
    @harmonic_torsions_end.setter
    def harmonic_torsions_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_i(self) -> list[int]:
        ...
    @harmonic_torsions_i.setter
    def harmonic_torsions_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_j(self) -> list[int]:
        ...
    @harmonic_torsions_j.setter
    def harmonic_torsions_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_k(self) -> list[int]:
        ...
    @harmonic_torsions_k.setter
    def harmonic_torsions_k(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_l(self) -> list[int]:
        ...
    @harmonic_torsions_l.setter
    def harmonic_torsions_l(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_phase(self) -> list[float]:
        ...
    @harmonic_torsions_phase.setter
    def harmonic_torsions_phase(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def harmonic_torsions_stiffness(self) -> list[float]:
        ...
    @harmonic_torsions_stiffness.setter
    def harmonic_torsions_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def nonbonded_cutoff(self) -> float:
        ...
    @nonbonded_cutoff.setter
    def nonbonded_cutoff(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def num_angles(self) -> int:
        ...
    @num_angles.setter
    def num_angles(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_atoms(self) -> int:
        ...
    @num_atoms.setter
    def num_atoms(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_bonds(self) -> int:
        ...
    @num_bonds.setter
    def num_bonds(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_exclusions(self) -> int:
        ...
    @num_exclusions.setter
    def num_exclusions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_harmonic_torsions(self) -> int:
        ...
    @num_harmonic_torsions.setter
    def num_harmonic_torsions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_molecules(self) -> int:
        ...
    @num_molecules.setter
    def num_molecules(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_nb_types(self) -> int:
        ...
    @num_nb_types.setter
    def num_nb_types(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_periodic_torsions(self) -> int:
        ...
    @num_periodic_torsions.setter
    def num_periodic_torsions(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_scaling14(self) -> int:
        ...
    @num_scaling14.setter
    def num_scaling14(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_urey_bradley(self) -> int:
        ...
    @num_urey_bradley.setter
    def num_urey_bradley(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_virtual_sites(self) -> int:
        ...
    @num_virtual_sites.setter
    def num_virtual_sites(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def num_z_matrix_rows(self) -> int:
        ...
    @num_z_matrix_rows.setter
    def num_z_matrix_rows(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def periodic_torsions_begin(self) -> list[int]:
        ...
    @periodic_torsions_begin.setter
    def periodic_torsions_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_end(self) -> list[int]:
        ...
    @periodic_torsions_end.setter
    def periodic_torsions_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_i(self) -> list[int]:
        ...
    @periodic_torsions_i.setter
    def periodic_torsions_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_improper(self) -> list[bool]:
        ...
    @periodic_torsions_improper.setter
    def periodic_torsions_improper(self, arg0: collections.abc.Sequence[bool]) -> None:
        ...
    @property
    def periodic_torsions_j(self) -> list[int]:
        ...
    @periodic_torsions_j.setter
    def periodic_torsions_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_k(self) -> list[int]:
        ...
    @periodic_torsions_k.setter
    def periodic_torsions_k(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_l(self) -> list[int]:
        ...
    @periodic_torsions_l.setter
    def periodic_torsions_l(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_n(self) -> list[int]:
        ...
    @periodic_torsions_n.setter
    def periodic_torsions_n(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_phase(self) -> list[float]:
        ...
    @periodic_torsions_phase.setter
    def periodic_torsions_phase(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def periodic_torsions_stiffness(self) -> list[float]:
        ...
    @periodic_torsions_stiffness.setter
    def periodic_torsions_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def root_mobilities(self) -> list[RootMobility]:
        ...
    @root_mobilities.setter
    def root_mobilities(self, arg0: collections.abc.Sequence[RootMobility]) -> None:
        ...
    @property
    def scaling14_begin(self) -> list[int]:
        ...
    @scaling14_begin.setter
    def scaling14_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_charge_product(self) -> list[float]:
        ...
    @scaling14_charge_product.setter
    def scaling14_charge_product(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_end(self) -> list[int]:
        ...
    @scaling14_end.setter
    def scaling14_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_epsilon(self) -> list[float]:
        ...
    @scaling14_epsilon.setter
    def scaling14_epsilon(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_i(self) -> list[int]:
        ...
    @scaling14_i.setter
    def scaling14_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_l(self) -> list[int]:
        ...
    @scaling14_l.setter
    def scaling14_l(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def scaling14_sigma(self) -> list[float]:
        ...
    @scaling14_sigma.setter
    def scaling14_sigma(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def seed(self) -> int:
        ...
    @seed.setter
    def seed(self, arg0: typing.SupportsInt | typing.SupportsIndex) -> None:
        ...
    @property
    def thermostat_temperature(self) -> float:
        ...
    @thermostat_temperature.setter
    def thermostat_temperature(self, arg0: typing.SupportsFloat | typing.SupportsIndex) -> None:
        ...
    @property
    def urey_bradley_begin(self) -> list[int]:
        ...
    @urey_bradley_begin.setter
    def urey_bradley_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_end(self) -> list[int]:
        ...
    @urey_bradley_end.setter
    def urey_bradley_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_equilibrium(self) -> list[float]:
        ...
    @urey_bradley_equilibrium.setter
    def urey_bradley_equilibrium(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_i(self) -> list[int]:
        ...
    @urey_bradley_i.setter
    def urey_bradley_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_k(self) -> list[int]:
        ...
    @urey_bradley_k.setter
    def urey_bradley_k(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def urey_bradley_stiffness(self) -> list[float]:
        ...
    @urey_bradley_stiffness.setter
    def urey_bradley_stiffness(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_atom1(self) -> list[int]:
        ...
    @vs_atom1.setter
    def vs_atom1(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_atom2(self) -> list[int]:
        ...
    @vs_atom2.setter
    def vs_atom2(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_atom3(self) -> list[int]:
        ...
    @vs_atom3.setter
    def vs_atom3(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_site(self) -> list[int]:
        ...
    @vs_site.setter
    def vs_site(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_weight1(self) -> list[float]:
        ...
    @vs_weight1.setter
    def vs_weight1(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_weight2(self) -> list[float]:
        ...
    @vs_weight2.setter
    def vs_weight2(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def vs_weight3(self) -> list[float]:
        ...
    @vs_weight3.setter
    def vs_weight3(self, arg0: collections.abc.Sequence[typing.SupportsFloat | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_begin(self) -> list[int]:
        ...
    @z_matrix_begin.setter
    def z_matrix_begin(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_end(self) -> list[int]:
        ...
    @z_matrix_end.setter
    def z_matrix_end(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_i(self) -> list[int]:
        ...
    @z_matrix_i.setter
    def z_matrix_i(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_j(self) -> list[int]:
        ...
    @z_matrix_j.setter
    def z_matrix_j(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_k(self) -> list[int]:
        ...
    @z_matrix_k.setter
    def z_matrix_k(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
    @property
    def z_matrix_l(self) -> list[int]:
        ...
    @z_matrix_l.setter
    def z_matrix_l(self, arg0: collections.abc.Sequence[typing.SupportsInt | typing.SupportsIndex]) -> None:
        ...
class World:
    def add_sampler(self, timeStep: typing.SupportsFloat | typing.SupportsIndex, mdSteps: typing.SupportsInt | typing.SupportsIndex, acceptRejectMode: AcceptRejectMode, use_nuts: bool, sphere_factor: typing.SupportsFloat | typing.SupportsIndex = 1.0, use_fixman: bool | None = None, always_kick: bool = False, clash_threshold: typing.SupportsFloat | typing.SupportsIndex = 10.0, max_initial_kick_tries: typing.SupportsInt | typing.SupportsIndex = 0) -> World:
        """
        Configure this world's sampler; returns the world for chaining. sphere_factor scales the auto-sized per-ligand binding sphere (R = R_receptor + sphere_factor*R_ligand). The docking kick relocates a ligand only when its COM leaves the sphere (always_kick=True perturbs every round). A proposed pose is rejected -- in ALL modes, including AlwaysAccept -- if its potential energy is non-finite or |PE| exceeds clash_threshold, so overlapping geometry never passes. use_fixman=None auto-enables Fixman+logSineSqr on non-Cartesian worlds. max_initial_kick_tries>0 enables a pre-round-0 retry loop that keeps drawing random placements until a clash-free starting pose is found (dPE <= maxStartPE), or raises RuntimeError after the budget is exhausted.
        """
    def configure_ncmc(self, atom_begin: typing.SupportsInt | typing.SupportsIndex, atom_end: typing.SupportsInt | typing.SupportsIndex, ncmc_steps: typing.SupportsInt | typing.SupportsIndex, hold_fraction: typing.SupportsFloat | typing.SupportsIndex = 0.0) -> None:
        """
        Make this a per-molecule NCMC world: soften [atom_begin,atom_end) x rest nonbonded during a lambda:1->0->1 switch. Call AFTER add_sampler.
        """
    def set_body_mass_scale(self, body: typing.SupportsInt | typing.SupportsIndex, scale: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Per-body form of set_mass_scale_by_joint (scale=1.0 is physical/off).
        """
    def set_mass_scale_by_joint(self, joint_type: JointType, scale: typing.SupportsFloat | typing.SupportsIndex) -> None:
        """
        Inflate the spatial inertia used ONLY in the proposal (momentum draw, KE, Fixman ln det M) for every body of the given JointType, raising the stable dt ~sqrt(scale) with zero configurational bias. scale=1.0 is physical (off). Typical: set_mass_scale_by_joint(JointType.Free, 16.0) on a solvent world to tame water libration.
        """
    def set_reversibility_check(self, interval: typing.SupportsInt | typing.SupportsIndex) -> None:
        """
        Low-level setter for the periodic reversibility probe. Prefer the reversibility_check_every=N argument to context.add_*_world(), which forwards here. interval<=0 disables it (default); interval>0 runs RobotEngine::checkReversibility every `interval` rounds (starting at round 0, so it doubles as a startup check), integrating mdSteps forward+back at this world's timestep from the freshly seeded state and logging the relative round-trip residual (~1e-12 ideal; a large or non-finite value means the timestep is too large for the current geometry). Non-destructive; it is a smoke test for the current configuration only -- the always-on guard is the per-step corrector throw in the integrator. See THEORY 5.7.
        """
    @property
    def index(self) -> int:
        ...
    @property
    def is_cartesian(self) -> bool:
        ...
    @property
    def is_docking(self) -> bool:
        ...
