"""amber_loader.py

Direct (parmed-free) reader for AMBER formatted coordinate files (``.rst7`` /
``.inpcrd``) plus the molecule-instance/prototype dedup. Part of the
fast-loader rewrite (``docs/specs/fast-amber-loader.md``):

* component (A) -- :func:`read_amber_coordinates`: replaces
  ``parm.atoms[i].xx/.xy/.xz`` with a single vectorized numpy read.
* component (B) -- :func:`partition_molecules`: replaces ``parm.split()``
  with ``scipy.sparse.csgraph.connected_components`` over the bond graph
  plus a conservative per-instance fingerprint.

File format (AMBER restart / inpcrd, ASCII)
--------------------------------------------
Mirrors ``openmm.app.internal.amber_file_parser.AmberAsciiRestart`` (the
reference this module is checked against -- see
``docs/specs/loader-feature-matrix.md``)::

    line 0:  title
    line 1:  NATOM [TIME]           -- whitespace-separated (width varies
                                        across writers; NOT a fixed I5/I6
                                        field -- some restarts pad it, some
                                        don't)
    coords:  ceil(NATOM/2) lines, 2 atoms (6 floats, F12.7) per line, last
             line may hold only 3 floats (odd NATOM)
    [vels]:  same shape as coords, present iff the file was written with
             velocities
    [box]:   one line, 6 floats (F12.7): a, b, c, alpha, beta, gamma

Whether velocities/box are present is inferred from the total line count
(the format carries no explicit flag for this), exactly as OpenMM's
``AmberAsciiRestart._parse`` does, including its documented ``NATOM in (1, 2)``
disambiguation edge case.
"""

from __future__ import annotations

import math
import os
from dataclasses import dataclass

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

from . import prmtop_reader
from .units import ANG_TO_NM, DEG_TO_RAD

__all__ = [
    "read_amber_coordinates",
    "partition_molecules",
    "MoleculePartition",
    "PrototypeTopology",
    "UnsupportedTopologyFeature",
    "box_vectors_from_lengths_angles",
    "VirtualSiteRecord",
    "extract_virtual_sites",
]


class UnsupportedTopologyFeature(Exception):
    """A prmtop feature was parsed but has no consumer yet (fast-loader
    spec, docs/specs/fast-amber-loader.md §4a "parsed-but-raises" policy).

    Raised only when the offending feature is actually present in the input
    (never unconditionally); the message names the feature and, where
    relevant, the OpenMM behavior being deferred.
    """


def _read_fixed_width_floats(line: str, n: int, width: int = 12) -> list[float]:
    return [float(line[i : i + width]) for i in range(0, n * width, width)]


def read_amber_coordinates(
    path: str | os.PathLike[str],
) -> tuple[np.ndarray, np.ndarray | None]:
    """Read atom coordinates (and, if present, the box) from an AMBER
    restart/inpcrd file, without parmed.

    Parameters
    ----------
    path : str or os.PathLike[str]
        Path to the ``.rst7``/``.inpcrd`` file.

    Returns
    -------
    coords_nm : (N, 3) float64 ndarray
        Atom positions in **prmtop atom order**, converted to nm.
    box : (6,) float64 ndarray or None
        ``[a, b, c, alpha, beta, gamma]`` -- box lengths in nm, angles in
        radians -- or ``None`` if the file carries no box line. NOTE: this is
        the raw box geometry, not reduced lattice vectors; ``Context.
        load_amber`` still sources ``system_topology.box_vectors`` from
        ParmEd's already-reduced vectors (see ``docs/specs/
        fast-amber-loader.md`` §3/§6 Step 2 -- box-vector construction is
        deferred to a later step of the rewrite).

    Raises
    ------
    ValueError
        If the file cannot be parsed as a well-formed AMBER restart (bad
        atom count, wrong number of lines, malformed box line).
    """
    with open(path, "r") as fh:
        lines = fh.readlines()
    while lines and not lines[-1].strip():
        lines.pop()

    if len(lines) < 2:
        raise ValueError(f"'{path}': too short to be an AMBER restart/inpcrd file.")

    try:
        natom = int(lines[1].split()[0])
    except (IndexError, ValueError) as exc:
        raise ValueError(f"'{path}': could not parse NATOM from line 2.") from exc

    n_coord_lines = int(math.ceil(natom / 2.0))

    # Disambiguate whether velocities and/or a box line follow the
    # coordinates, from the total line count alone (the format has no
    # explicit flag for this) -- mirrors
    # openmm.app.internal.amber_file_parser.AmberAsciiRestart._parse.
    if len(lines) == n_coord_lines + 2:
        has_box = has_vels = False
    elif natom in (1, 2) and len(lines) == 4:
        line = lines[3]
        if natom == 1:
            tmp = [line[i : i + 12] for i in range(0, 72, 12) if line[i : i + 12].strip()]
            if len(tmp) == 3:
                has_vels, has_box = True, False
            elif len(tmp) == 6:
                has_box, has_vels = True, False
            else:
                raise ValueError(f"'{path}': unrecognized line in restart file.")
        else:
            # Ambiguous case: velocities (scaled by ~20.445) have ~0% chance
            # of exceeding 60.0, so a value >= 60 in the last line means box
            # lengths/angles, not velocities.
            is_box_like = any(
                float(line[i : i + 12]) >= 60.0 for i in range(0, 72, 12)
            )
            has_box, has_vels = (True, False) if is_box_like else (False, True)
    elif len(lines) == n_coord_lines + 3:
        has_box, has_vels = True, False
    elif len(lines) == 2 * n_coord_lines + 2:
        has_box, has_vels = False, True
    elif len(lines) == 2 * n_coord_lines + 3:
        has_box = has_vels = True
    else:
        raise ValueError(
            f"'{path}': badly formatted restart file -- {len(lines)} lines for {natom} atoms."
        )

    coord_lines = lines[2 : 2 + n_coord_lines]
    # The coordinate block is a fixed-width Fortran 6F12.7 record: fields abut
    # with no separator, so a value filling all 12 columns (a negative magnitude
    # >= 100, or a positive >= 1000 -- both normal for large unwrapped systems)
    # has no space before its neighbour. Parse by 12-column slices, NOT
    # str.split() (which fuses such neighbours into one bad token), mirroring
    # openmm.app.internal.amber_file_parser.AmberAsciiRestart and the box line
    # below (which already used _read_fixed_width_floats).
    coords_list: list[float] = []
    for coord_line in coord_lines:
        stripped = coord_line.rstrip("\n")
        for i in range(0, len(stripped), 12):
            field = stripped[i : i + 12]
            if field.strip():
                coords_list.append(float(field))
    coords_flat = np.array(coords_list, dtype=np.float64)
    if coords_flat.size != natom * 3:
        raise ValueError(
            f"'{path}': parsed {coords_flat.size} coordinate values, expected {natom * 3}."
        )
    coords_nm = coords_flat.reshape(natom, 3) * ANG_TO_NM

    box: np.ndarray | None = None
    if has_box:
        box_line_idx = 2 + n_coord_lines + (n_coord_lines if has_vels else 0)
        box_vals = _read_fixed_width_floats(lines[box_line_idx], 6)
        lengths_nm = [v * ANG_TO_NM for v in box_vals[:3]]
        angles_rad = [v * DEG_TO_RAD for v in box_vals[3:]]
        box = np.array(lengths_nm + angles_rad, dtype=np.float64)

    return coords_nm, box


# ====================================================================== #
# Component (B): molecule-instance / prototype dedup -- replaces
# ``parm.split()`` (docs/specs/fast-amber-loader.md §3(B)/§5).
# ====================================================================== #


@dataclass(frozen=True)
class MoleculePartition:
    """Molecule-instance / prototype grouping (component (B) of the
    fast-loader spec).

    Replaces ParmEd's ``Structure.split()``: partitions the whole-system atom
    set into connected-component molecule *instances* (bond-graph
    connectivity only, ring-closing bonds included) and groups instances that
    are provably identical -- same per-atom identity/parameters AND the same
    intra-molecule bond graph -- into *prototypes*.

    Attributes
    ----------
    instance_atoms : list[np.ndarray]
        One entry per molecule instance: ascending 0-based prmtop atom
        indices, list position == ``instance_idx``. Instances are numbered by
        first-appearing atom in prmtop order, matching ParmEd's
        ``tag_molecules``.
    prototype_of_instance : list[int]
        ``prototype_of_instance[instance_idx]`` -> prototype id, parallel to
        ``instance_atoms``.
    prototype_representative_atoms : list[np.ndarray]
        One entry per UNIQUE prototype, ordered by first appearance (the
        ``instance_idx`` of the first instance assigned to it) -- matching
        the order ParmEd's ``split()`` builds its ``structs`` list.  Ascending
        0-based prmtop atom indices of the *representative* (first-occurring)
        instance of that prototype: the atom set to slice the full ``parm``
        structure at (``parm[mask]``) to build that prototype's
        ``pmd.Structure``.
    """

    instance_atoms: list[np.ndarray]
    prototype_of_instance: list[int]
    prototype_representative_atoms: list[np.ndarray]


def _lj_diagonal_per_atom(raw_data: dict) -> tuple[np.ndarray, np.ndarray]:
    """Per-atom ``(Rmin/2, epsilon)`` from the diagonal AMBER LJ self-terms.

    Mirrors ParmEd's ``Atom.rmin`` / ``Atom.epsilon`` -- both are per LJ
    *type*, not per pair: ``Rmin_i = (2*A_ii/B_ii)^(1/6)``,
    ``epsilon_i = B_ii^2 / (4*A_ii)`` from the diagonal of the
    ``NONBONDED_PARM_INDEX`` table (the same source
    ``prmtop_reader.has_nbfix_fast`` uses for its combining-rule check).

    Used only for the dedup fingerprint (:func:`partition_molecules`) -- NOT
    the per-atom sigma/epsilon fed to ``MoleculePrototype`` (those still come
    from ParmEd on the sliced prototype structure, unchanged).
    """
    num_types = int(raw_data["POINTERS"][1])
    atom_type_index = np.asarray(raw_data["ATOM_TYPE_INDEX"], dtype=np.int64) - 1
    nb_index = np.asarray(raw_data["NONBONDED_PARM_INDEX"], dtype=np.int64).reshape(
        num_types, num_types
    )
    acoef = np.asarray(raw_data["LENNARD_JONES_ACOEF"], dtype=np.float64)
    bcoef = np.asarray(raw_data["LENNARD_JONES_BCOEF"], dtype=np.float64)

    diag_idx = nb_index.diagonal() - 1  # 0-based into ACOEF/BCOEF, per LJ type
    a_ii = acoef[diag_idx]
    b_ii = bcoef[diag_idx]
    # Zero out near-zero self-interaction coefficients (dummy/virtual-site
    # types) using the SAME 1e-10 absolute threshold as ParmEd's
    # AmberParm.fill_LJ (parmed/amber/_amberparm.py) -- not just exact-zero /
    # non-finite -- so this reproduces parmed's Atom.rmin / Atom.epsilon
    # (and, transitively, Atom.sigma) exactly for the PrototypeTopology atom
    # builder (amber_loader.PrototypeTopology), not merely "close".
    zero_mask = (a_ii < 1e-10) | (b_ii < 1e-10)
    with np.errstate(divide="ignore", invalid="ignore"):
        rmin_full = (2.0 * a_ii / b_ii) ** (1.0 / 6.0)
        eps = 0.25 * b_ii**2 / a_ii
    rmin_half = np.where(zero_mask, 0.0, rmin_full / 2.0)
    eps = np.where(zero_mask, 0.0, eps)

    return rmin_half[atom_type_index], eps[atom_type_index]


def partition_molecules(raw_data: dict, natom: int) -> MoleculePartition:
    """Partition the whole-system atom set into molecule instances + prototypes.

    Replaces ``parm.split()``. Two steps:

    1. **Instances.** Build the bond graph from ``BONDS_INC_HYDROGEN`` +
       ``BONDS_WITHOUT_HYDROGEN`` (connectivity only; ring-closing bonds
       included) and find its connected components with
       ``scipy.sparse.csgraph.connected_components`` -- one component per
       molecule instance. Instances are renumbered by first-appearing atom in
       prmtop order, matching ParmEd's ``tag_molecules``.
    2. **Prototypes.** Group instances by a conservative fingerprint: the
       per-atom ``(residue_name, atom_name, charge, rmin, epsilon)`` tuple
       (ParmEd's single-residue ``split()`` key) PLUS the intra-instance bond
       set (component-local indices). Identical fingerprint implies
       identical molecule; adding the bond set only ever makes the
       fingerprint *finer* than ParmEd's (safe to over-split -- never merges
       non-identical molecules -- see ``docs/specs/fast-amber-loader.md``
       §5).

    Parameters
    ----------
    raw_data : dict
        The ``raw_data`` dict from
        :func:`robosample.prmtop_reader.parse_prmtop`.
    natom : int
        Total atom count (``len(raw_data["ATOM_NAME"])``).

    Returns
    -------
    MoleculePartition
    """
    # ---- 1. Bond graph -> connected components ----------------------------
    bonds_raw = np.concatenate(
        [
            np.asarray(raw_data.get("BONDS_INC_HYDROGEN", []), dtype=np.int64),
            np.asarray(raw_data.get("BONDS_WITHOUT_HYDROGEN", []), dtype=np.int64),
        ]
    ).reshape(-1, 3)
    # Each triplet is (3*atom_i, 3*atom_j, bond_type_idx); atom indices are
    # 0-based once divided by 3 (same ÷3 decoding as the DIHEDRALS pointers in
    # prmtop_reader.load_nonbonded_exceptions).
    bond_i = bonds_raw[:, 0] // 3
    bond_j = bonds_raw[:, 1] // 3

    if bond_i.size:
        data = np.ones(bond_i.size, dtype=np.int8)
        graph = coo_matrix((data, (bond_i, bond_j)), shape=(natom, natom))
    else:
        graph = coo_matrix((natom, natom), dtype=np.int8)

    n_instances, raw_labels = connected_components(graph, directed=False)

    # Renumber components by first-appearing atom (ascending atom index),
    # matching ParmEd's tag_molecules numbering.
    first_atom_of_label = np.full(n_instances, natom, dtype=np.int64)
    np.minimum.at(first_atom_of_label, raw_labels, np.arange(natom, dtype=np.int64))
    instance_of_atom = np.argsort(np.argsort(first_atom_of_label))[raw_labels]

    # Group atom indices per instance, ascending within each instance (stable
    # sort by instance id preserves ascending original atom order for ties,
    # since the pre-sort order is 0..natom-1).
    sort_idx = np.argsort(instance_of_atom, kind="stable")
    sorted_instance = instance_of_atom[sort_idx]
    boundaries = np.searchsorted(sorted_instance, np.arange(n_instances + 1))
    instance_atoms = [
        sort_idx[boundaries[i] : boundaries[i + 1]] for i in range(n_instances)
    ]

    # ---- 2. Per-atom fingerprint fields ------------------------------------
    residue_starts0 = np.asarray(raw_data["RESIDUE_POINTER"], dtype=np.int64) - 1
    residue_of_atom = (
        np.searchsorted(residue_starts0, np.arange(natom), side="right") - 1
    )
    resname_of_atom = np.asarray(raw_data["RESIDUE_LABEL"])[residue_of_atom]
    atom_name_of_atom = np.asarray(raw_data["ATOM_NAME"])
    charge_of_atom = np.asarray(raw_data["CHARGE"], dtype=np.float64)
    rmin_of_atom, epsilon_of_atom = _lj_diagonal_per_atom(raw_data)

    # Local (within-instance, ascending-atom-order) rank of every atom --
    # used to remap bonds to component-local indices for the fingerprint.
    pos_in_sorted = np.empty(natom, dtype=np.int64)
    pos_in_sorted[sort_idx] = np.arange(natom, dtype=np.int64)
    local_rank_of_atom = pos_in_sorted - boundaries[instance_of_atom]

    # Bonds grouped by instance (both endpoints of any bond edge are
    # guaranteed to share a connected component, hence the same instance id).
    edge_instance = instance_of_atom[bond_i]
    edge_order = np.argsort(edge_instance, kind="stable")
    edge_instance_sorted = edge_instance[edge_order]
    edge_local_i = local_rank_of_atom[bond_i][edge_order]
    edge_local_j = local_rank_of_atom[bond_j][edge_order]
    edge_boundaries = np.searchsorted(edge_instance_sorted, np.arange(n_instances + 1))

    def _fingerprint(instance_idx: int) -> tuple:
        atoms = instance_atoms[instance_idx]
        atom_key = tuple(
            (
                str(resname_of_atom[a]),
                str(atom_name_of_atom[a]),
                round(float(charge_of_atom[a]), 6),
                round(float(rmin_of_atom[a]), 6),
                round(float(epsilon_of_atom[a]), 6),
            )
            for a in atoms
        )
        e0, e1 = edge_boundaries[instance_idx], edge_boundaries[instance_idx + 1]
        bond_key = tuple(
            sorted(
                (int(min(i, j)), int(max(i, j)))
                for i, j in zip(
                    edge_local_i[e0:e1].tolist(), edge_local_j[e0:e1].tolist()
                )
            )
        )
        return (atom_key, bond_key)

    # ---- 3. Group instances into prototypes (first-seen order) ------------
    prototype_of_fingerprint: dict[tuple, int] = {}
    prototype_of_instance: list[int] = [0] * n_instances
    prototype_representative_atoms: list[np.ndarray] = []
    for i in range(n_instances):
        key = _fingerprint(i)
        pid = prototype_of_fingerprint.get(key)
        if pid is None:
            pid = len(prototype_representative_atoms)
            prototype_of_fingerprint[key] = pid
            prototype_representative_atoms.append(instance_atoms[i])
        prototype_of_instance[i] = pid

    return MoleculePartition(
        instance_atoms=instance_atoms,
        prototype_of_instance=prototype_of_instance,
        prototype_representative_atoms=prototype_representative_atoms,
    )


# ====================================================================== #
# Component (C): PrototypeTopology shim -- replaces the parmed ``Structure``
# built per prototype via ``parm[mask]`` (docs/specs/fast-amber-loader.md
# §3(C)). A lightweight, numpy/raw-array-backed object exposing exactly the
# ParmEd attribute surface consumed by ``MoleculePrototype``,
# ``acyclic_graph``, ``z_matrix``, and ``amber_dihedral_classifier``, so those
# modules run unchanged (module docstrings there still say "pmd.Structure" /
# "pmd.Atom" for the type hints -- those are non-binding under
# ``from __future__ import annotations`` and are left as-is; the objects
# handed to them at runtime are the shims defined below, never a real ParmEd
# object).
#
# Units and index-space conventions exactly mirror what a ParmEd-sliced
# ``Structure`` would expose (verified empirically against ParmEd for several
# example systems, incl. CHAMBER -- see docs/specs/fast-amber-loader.md):
#   * atom.sigma / atom.epsilon   : Å / kcal/mol (raw AMBER units; the SAME
#                                    formula ParmEd's AmberParm.fill_LJ uses --
#                                    sigma = (rmin/2)*2^(-1/6)*2, from the
#                                    diagonal LJ A/B self-terms).
#   * atom.solvent_radius/.screen : Å / dimensionless (raw prmtop RADII/SCREEN).
#   * atom.xx/.xy/.xz             : Å (raw prmtop/inpcrd coordinates).
#   * bond.type.k/.req            : kcal/mol/Å^2, Å            (raw prmtop).
#   * angle.type.k/.theteq        : kcal/mol/rad^2, DEGREES    (ParmEd converts
#                                    the raw prmtop radians to degrees when
#                                    building AngleType; MoleculePrototype
#                                    converts back with DEG_TO_RAD).
#   * dihedral.type.phi_k/.per/.phase/.scee/.scnb :
#                                    kcal/mol, unitless, DEGREES, unitless,
#                                    unitless (phase likewise converted
#                                    rad->deg to match ParmEd's DihedralType).
#   * improper.type.psi_k/.psi_eq : kcal/mol/rad^2, DEGREES (CHARMM_IMPROPER_
#                                    PHASE is radians in the prmtop; ParmEd's
#                                    ChamberParm._load_improper_info converts
#                                    to degrees -- replicated here).
#   * urey_bradley.type.k/.req    : kcal/mol/Å^2, Å (raw prmtop, like bonds).
# AMBER/CHAMBER-loaded dihedrals are NEVER grouped into a ParmEd
# ``DihedralTypeList`` (verified: ``AmberParm``/``ChamberParm`` never call
# ``Structure.join_dihedrals()``) -- multi-term torsions appear as several
# separate single-term ``Dihedral`` records sharing the same 4 atoms, exactly
# as the raw ``DIHEDRALS_*_HYDROGEN`` sections store them. ``.type`` is
# therefore always a single value-holder object here, never a list.
# ====================================================================== #

# Periodic-table symbols indexed by atomic number (index 0 == "EP", ParmEd's
# placeholder for extra points / virtual sites). A static copy of
# ``parmed.periodic_table.Element`` kept here so ``atoms_element_name`` /
# ``atoms_element_symbol`` stay byte-identical to the pre-rewrite (ParmEd
# ``Atom.element_name``) values without importing parmed into this
# (deliberately parmed-free) module.
_ELEMENT_SYMBOLS: tuple[str, ...] = (
    "EP", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg",
    "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn",
    "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb",
    "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm",
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta",
    "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At",
    "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og",
)

# sigma = (Rmin/2) * 2^(-1/6) * 2  -- ParmEd's Atom.sigma formula (Atom.rmin
# IS Rmin/2, the per-atom LJ radius already returned by _lj_diagonal_per_atom).
_SIGMA_FULL_SCALE: float = 2.0 ** (-1.0 / 6.0) * 2.0


@dataclass(frozen=True, slots=True)
class _BondTypeShim:
    k: float
    req: float


@dataclass(frozen=True, slots=True)
class _AngleTypeShim:
    k: float
    theteq: float
    """Degrees (matches ParmEd's AngleType.theteq convention)."""


@dataclass(frozen=True, slots=True)
class _DihedralTypeShim:
    phi_k: float
    per: int
    """ParmEd's DihedralType.__init__ does ``self.per = int(per)`` -- the raw
    prmtop DIHEDRAL_PERIODICITY section stores it as a float (%FORMAT 5E16.8)
    even though the value is always integral; replicate the int() cast so
    ``periodic_torsions_n`` stays an int array, matching the pybind11
    SystemTopology field signature (Sequence[SupportsInt])."""
    phase: float
    """Degrees (matches ParmEd's DihedralType.phase convention)."""
    scee: float
    scnb: float


@dataclass(frozen=True, slots=True)
class _ImproperTypeShim:
    psi_k: float
    psi_eq: float
    """Degrees (matches ParmEd's ImproperType.psi_eq convention)."""


class _ResidueShim:
    __slots__ = ("idx", "name")

    def __init__(self, idx: int, name: str) -> None:
        self.idx = idx
        self.name = name


class _AtomShim:
    """Duck-types the ParmEd ``Atom`` attribute surface consumed downstream
    (``MoleculePrototype``, ``acyclic_graph``, ``z_matrix``,
    ``amber_dihedral_classifier``, ``bond_util.is_rigid_bond``)."""

    __slots__ = (
        "idx", "name", "mass", "charge", "sigma", "epsilon", "solvent_radius",
        "screen", "nb_idx", "atomic_number", "element_name", "type",
        "residue", "xx", "xy", "xz", "bond_partners", "bonds", "dihedrals",
    )

    def __init__(
        self,
        idx: int,
        name: str,
        mass: float,
        charge: float,
        sigma: float,
        epsilon: float,
        solvent_radius: float,
        screen: float,
        nb_idx: int,
        atomic_number: int,
        element_name: str,
        type: str,  # noqa: A002 (matches parmed.Atom.type's name)
        residue: _ResidueShim,
        xx: float,
        xy: float,
        xz: float,
    ) -> None:
        self.idx = idx
        self.name = name
        self.mass = mass
        self.charge = charge
        self.sigma = sigma
        self.epsilon = epsilon
        self.solvent_radius = solvent_radius
        self.screen = screen
        self.nb_idx = nb_idx
        self.atomic_number = atomic_number
        self.element_name = element_name
        self.type = type
        self.residue = residue
        self.xx = xx
        self.xy = xy
        self.xz = xz
        # Populated once bonds are built (see _extract_bonds): bond_partners
        # is every OTHER atom directly bonded to this one (both hydrogen and
        # non-hydrogen sections, ring-closing bonds included -- the FULL bond
        # graph degree), ascending by local idx (matches ParmEd's
        # Atom.bond_partners, which sorts by Atom.idx).
        self.bond_partners: list[_AtomShim] = []
        # Bond objects touching this atom (needed by bond_util.is_rigid_bond's
        # `_find_bond` fast path; AMBER-loaded bonds always carry order=1.0,
        # so that fast path never actually triggers, but the attribute must
        # exist and behave like ParmEd's).
        self.bonds: list[_BondShim] = []
        # Never populated (nothing downstream reads atom.dihedrals -- only the
        # Structure-level molecule.dihedrals collection is consumed); kept so
        # any accidental future access fails on an empty list, not AttributeError.
        self.dihedrals: list = []


class _BondShim:
    __slots__ = ("atom1", "atom2", "type", "order")

    def __init__(self, atom1: _AtomShim, atom2: _AtomShim, type: _BondTypeShim) -> None:  # noqa: A002
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type
        self.order = 1.0  # AMBER-loaded bonds never carry a resolved bond order.


class _AngleShim:
    __slots__ = ("atom1", "atom2", "atom3", "type")

    def __init__(self, atom1, atom2, atom3, type: _AngleTypeShim) -> None:  # noqa: A002
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.type = type


class _DihedralShim:
    __slots__ = ("atom1", "atom2", "atom3", "atom4", "type", "improper", "ignore_end")

    def __init__(
        self, atom1, atom2, atom3, atom4, type: _DihedralTypeShim, improper: bool, ignore_end: bool  # noqa: A002
    ) -> None:
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.type = type
        self.improper = improper
        self.ignore_end = ignore_end


class _ImproperShim:
    __slots__ = ("atom1", "atom2", "atom3", "atom4", "type")

    def __init__(self, atom1, atom2, atom3, atom4, type: _ImproperTypeShim) -> None:  # noqa: A002
        self.atom1 = atom1
        self.atom2 = atom2
        self.atom3 = atom3
        self.atom4 = atom4
        self.type = type


class _UreyBradleyShim:
    __slots__ = ("atom1", "atom2", "type")

    def __init__(self, atom1, atom2, type: _BondTypeShim) -> None:  # noqa: A002
        self.atom1 = atom1
        self.atom2 = atom2
        self.type = type


def _extract_bonds(
    raw_data: dict, local_of_global: np.ndarray, atoms: list[_AtomShim]
) -> list[_BondShim]:
    force_k = np.asarray(raw_data["BOND_FORCE_CONSTANT"], dtype=np.float64)
    equil = np.asarray(raw_data["BOND_EQUIL_VALUE"], dtype=np.float64)
    bonds: list[_BondShim] = []
    # WITHOUT_HYDROGEN before INC_HYDROGEN -- matches ParmEd's
    # AmberParm._load_bond_info build order for Structure.bonds.
    for section in ("BONDS_WITHOUT_HYDROGEN", "BONDS_INC_HYDROGEN"):
        raw = np.asarray(raw_data.get(section, []), dtype=np.int64).reshape(-1, 3)
        for gi3, gj3, t in raw.tolist():
            gi, gj = gi3 // 3, gj3 // 3
            li = int(local_of_global[gi])
            if li < 0:
                continue
            lj = int(local_of_global[gj])
            btype = _BondTypeShim(k=float(force_k[t - 1]), req=float(equil[t - 1]))
            bond = _BondShim(atoms[li], atoms[lj], btype)
            atoms[li].bonds.append(bond)
            atoms[lj].bonds.append(bond)
            bonds.append(bond)
    # bond_partners: every atom directly bonded, ascending local idx.
    partner_ids: list[set[int]] = [set() for _ in atoms]
    for b in bonds:
        partner_ids[b.atom1.idx].add(b.atom2.idx)
        partner_ids[b.atom2.idx].add(b.atom1.idx)
    for a in atoms:
        a.bond_partners = [atoms[j] for j in sorted(partner_ids[a.idx])]
    return bonds


def _extract_angles(
    raw_data: dict, local_of_global: np.ndarray, atoms: list[_AtomShim]
) -> list[_AngleShim]:
    force_k = np.asarray(raw_data["ANGLE_FORCE_CONSTANT"], dtype=np.float64)
    equil_rad = np.asarray(raw_data["ANGLE_EQUIL_VALUE"], dtype=np.float64)
    angles: list[_AngleShim] = []
    for section in ("ANGLES_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN"):
        raw = np.asarray(raw_data.get(section, []), dtype=np.int64).reshape(-1, 4)
        for gi3, gj3, gk3, t in raw.tolist():
            gi, gj, gk = gi3 // 3, gj3 // 3, gk3 // 3
            li = int(local_of_global[gi])
            if li < 0:
                continue
            lj, lk = int(local_of_global[gj]), int(local_of_global[gk])
            atype = _AngleTypeShim(
                k=float(force_k[t - 1]), theteq=math.degrees(float(equil_rad[t - 1]))
            )
            angles.append(_AngleShim(atoms[li], atoms[lj], atoms[lk], atype))
    return angles


def _extract_dihedrals(
    raw_data: dict, local_of_global: np.ndarray, atoms: list[_AtomShim]
) -> list[_DihedralShim]:
    force_k = np.asarray(raw_data["DIHEDRAL_FORCE_CONSTANT"], dtype=np.float64)
    per = np.asarray(raw_data["DIHEDRAL_PERIODICITY"], dtype=np.float64)
    phase_rad = np.asarray(raw_data["DIHEDRAL_PHASE"], dtype=np.float64)
    scee = raw_data.get("SCEE_SCALE_FACTOR")
    scee = np.full(force_k.shape, 1.2) if scee is None else np.asarray(scee, dtype=np.float64)
    scnb = raw_data.get("SCNB_SCALE_FACTOR")
    scnb = np.full(force_k.shape, 2.0) if scnb is None else np.asarray(scnb, dtype=np.float64)

    dihedrals: list[_DihedralShim] = []
    # WITHOUT_HYDROGEN before INC_HYDROGEN -- matches ParmEd's
    # AmberParm._load_dihedral_info build order for Structure.dihedrals (the
    # ENERGY-term consumer; see prmtop_reader.load_nonbonded_exceptions for
    # the DIFFERENT INC-before-WITHOUT order used by the 1-4/exclusion table,
    # which mirrors that pre-existing function's own convention instead).
    for section in ("DIHEDRALS_WITHOUT_HYDROGEN", "DIHEDRALS_INC_HYDROGEN"):
        raw = np.asarray(raw_data.get(section, []), dtype=np.int64).reshape(-1, 5)
        for gi3, gj3, k_raw, l_raw, t in raw.tolist():
            gi, gj = gi3 // 3, gj3 // 3
            li = int(local_of_global[gi])
            if li < 0:
                continue
            lj = int(local_of_global[gj])
            ignore_end = k_raw < 0
            improper = l_raw < 0
            gk, gl = abs(k_raw) // 3, abs(l_raw) // 3
            lk, ll = int(local_of_global[gk]), int(local_of_global[gl])
            dtype = _DihedralTypeShim(
                phi_k=float(force_k[t - 1]),
                per=int(per[t - 1]),
                phase=math.degrees(float(phase_rad[t - 1])),
                scee=float(scee[t - 1]),
                scnb=float(scnb[t - 1]),
            )
            dihedrals.append(
                _DihedralShim(
                    atoms[li], atoms[lj], atoms[lk], atoms[ll], dtype, bool(improper), bool(ignore_end)
                )
            )
    return dihedrals


def _extract_impropers(
    raw_data: dict, local_of_global: np.ndarray, atoms: list[_AtomShim]
) -> list[_ImproperShim]:
    """CHAMBER harmonic (CHARMM-style) impropers -- CHARMM_IMPROPERS section.

    Distinct from the periodic (cosine) improper torsions carried in
    DIHEDRALS_*_HYDROGEN with a negative 4th pointer (those become
    ``_DihedralShim(..., improper=True)`` records, handled by
    ``_extract_dihedrals``).
    """
    if "CHARMM_IMPROPERS" not in raw_data:
        return []
    force_k = np.asarray(raw_data["CHARMM_IMPROPER_FORCE_CONSTANT"], dtype=np.float64)
    eq_rad = np.asarray(raw_data["CHARMM_IMPROPER_PHASE"], dtype=np.float64)
    raw = np.asarray(raw_data["CHARMM_IMPROPERS"], dtype=np.int64).reshape(-1, 5)
    two_pi = 2.0 * math.pi
    impropers: list[_ImproperShim] = []
    for i, j, k, l, t in raw.tolist():  # noqa: E741 (matches AMBER field names)
        gi, gj, gk, gl = i - 1, j - 1, k - 1, l - 1  # plain 1-based, NOT x3-encoded
        li = int(local_of_global[gi])
        if li < 0:
            continue
        lj, lk, ll = (
            int(local_of_global[gj]),
            int(local_of_global[gk]),
            int(local_of_global[gl]),
        )
        eq = float(eq_rad[t - 1])
        # ParmEd's ChamberParm._load_improper_info heuristic: prmtop stores
        # CHARMM_IMPROPER_PHASE in radians; convert to degrees (the branch
        # that would keep it as-is only fires for already-degree legacy files
        # with |eq| > 2*pi, which a freshly parsed prmtop never is).
        eq_deg = math.degrees(eq) if abs(eq) <= two_pi else eq
        itype = _ImproperTypeShim(psi_k=float(force_k[t - 1]), psi_eq=eq_deg)
        impropers.append(_ImproperShim(atoms[li], atoms[lj], atoms[lk], atoms[ll], itype))
    return impropers


def _extract_urey_bradleys(
    raw_data: dict, local_of_global: np.ndarray, atoms: list[_AtomShim]
) -> list[_UreyBradleyShim]:
    """CHAMBER Urey-Bradley (1,3-pair) terms -- CHARMM_UREY_BRADLEY section."""
    if "CHARMM_UREY_BRADLEY" not in raw_data:
        return []
    force_k = np.asarray(raw_data["CHARMM_UREY_BRADLEY_FORCE_CONSTANT"], dtype=np.float64)
    equil = np.asarray(raw_data["CHARMM_UREY_BRADLEY_EQUIL_VALUE"], dtype=np.float64)
    raw = np.asarray(raw_data["CHARMM_UREY_BRADLEY"], dtype=np.int64).reshape(-1, 3)
    ubs: list[_UreyBradleyShim] = []
    for i, j, t in raw.tolist():
        gi, gj = i - 1, j - 1  # plain 1-based, NOT x3-encoded
        li = int(local_of_global[gi])
        if li < 0:
            continue
        lj = int(local_of_global[gj])
        ubtype = _BondTypeShim(k=float(force_k[t - 1]), req=float(equil[t - 1]))
        ubs.append(_UreyBradleyShim(atoms[li], atoms[lj], ubtype))
    return ubs


class PrototypeTopology:
    """Numpy/raw-array-backed replacement for the ParmEd ``Structure`` a
    prototype used to be built from via ``parm[mask]``
    (docs/specs/fast-amber-loader.md §3(C)).

    Exposes ``atoms``, ``bonds``, ``angles``, ``dihedrals``, ``impropers``,
    ``urey_bradleys``, ``residues``, and ``nonbonded_tables`` -- exactly the
    surface ``MoleculePrototype``/``acyclic_graph``/``z_matrix``/
    ``amber_dihedral_classifier`` consume (see module-level docstring section
    above for the full unit/index-space contract). Built once per UNIQUE
    prototype from the WHOLE-SYSTEM ``raw_data`` plus one representative
    instance's GLOBAL atom-index set -- no ParmEd object is built or sliced.
    """

    def __init__(
        self,
        raw_data: dict,
        atom_indices: np.ndarray,
        coords_nm: np.ndarray,
        residue_of_atom: np.ndarray,
        resname_by_residue: np.ndarray,
    ) -> None:
        """
        Parameters
        ----------
        raw_data : dict
            Whole-system ``raw_data`` from :func:`prmtop_reader.parse_prmtop`.
        atom_indices : np.ndarray
            Ascending 0-based GLOBAL prmtop atom indices of the representative
            instance this prototype is built from (see
            :class:`MoleculePartition`).
        coords_nm : np.ndarray
            Whole-system ``(N, 3)`` coordinate array, nm (from
            :func:`read_amber_coordinates`); sliced at *atom_indices* and
            converted back to Å here (``atom.xx/.xy/.xz`` -- ParmEd convention).
        residue_of_atom : np.ndarray
            Whole-system, 0-based GLOBAL residue index per atom.
        resname_by_residue : np.ndarray
            Whole-system residue-name array, indexed by GLOBAL residue index.
        """
        natom_total = len(raw_data["ATOM_NAME"])
        local_of_global = np.full(natom_total, -1, dtype=np.int64)
        local_of_global[atom_indices] = np.arange(len(atom_indices), dtype=np.int64)

        charge = np.asarray(raw_data["CHARGE"], dtype=np.float64)
        mass = np.asarray(raw_data["MASS"], dtype=np.float64)
        nb_idx_arr = np.asarray(raw_data["ATOM_TYPE_INDEX"], dtype=np.int64)
        atom_name_arr = np.asarray(raw_data["ATOM_NAME"])
        amber_type_arr = np.asarray(raw_data["AMBER_ATOM_TYPE"])
        if "ATOMIC_NUMBER" not in raw_data:
            raise ValueError(
                "prmtop is missing the ATOMIC_NUMBER section; the fast loader "
                "requires it (legacy mass-based element inference is not "
                "implemented)."
            )
        atomic_number_arr = np.asarray(raw_data["ATOMIC_NUMBER"], dtype=np.int64)
        radii_arr = np.asarray(
            raw_data.get("RADII", np.zeros(natom_total)), dtype=np.float64
        )
        screen_arr = np.asarray(
            raw_data.get("SCREEN", np.zeros(natom_total)), dtype=np.float64
        )
        rmin_half, eps = _lj_diagonal_per_atom(raw_data)  # whole-system, per-atom

        # ---- Residues: local idx = rank of first appearance in ascending
        # atom order (matches ParmEd's per-prototype residue renumbering; safe
        # because AMBER residues are always contiguous atom ranges, so global
        # residue ids along atom_indices are already non-decreasing and
        # np.unique's SORTED order equals first-appearance order). ------------
        global_res = residue_of_atom[atom_indices]
        uniq_res = np.unique(global_res)
        local_res_of_global_res = {int(g): i for i, g in enumerate(uniq_res.tolist())}
        residues = [
            _ResidueShim(idx=i, name=str(resname_by_residue[int(g)]))
            for i, g in enumerate(uniq_res.tolist())
        ]

        atoms: list[_AtomShim] = []
        coords_ang = coords_nm[atom_indices] / ANG_TO_NM
        for p, g in enumerate(atom_indices.tolist()):
            residue = residues[local_res_of_global_res[int(global_res[p])]]
            z = int(atomic_number_arr[g])
            atoms.append(
                _AtomShim(
                    idx=p,
                    name=str(atom_name_arr[g]),
                    mass=float(mass[g]),
                    charge=float(charge[g]),
                    sigma=float(rmin_half[g]) * _SIGMA_FULL_SCALE,
                    epsilon=float(eps[g]),
                    solvent_radius=float(radii_arr[g]),
                    screen=float(screen_arr[g]),
                    nb_idx=int(nb_idx_arr[g]),
                    atomic_number=z,
                    element_name=_ELEMENT_SYMBOLS[z] if 0 <= z < len(_ELEMENT_SYMBOLS) else "EP",
                    type=str(amber_type_arr[g]),
                    residue=residue,
                    xx=float(coords_ang[p, 0]),
                    xy=float(coords_ang[p, 1]),
                    xz=float(coords_ang[p, 2]),
                )
            )

        self.atoms: list[_AtomShim] = atoms
        self.residues: list[_ResidueShim] = residues
        self.bonds: list[_BondShim] = _extract_bonds(raw_data, local_of_global, atoms)
        self.angles: list[_AngleShim] = _extract_angles(raw_data, local_of_global, atoms)
        self.dihedrals: list[_DihedralShim] = _extract_dihedrals(
            raw_data, local_of_global, atoms
        )
        self.impropers: list[_ImproperShim] = _extract_impropers(
            raw_data, local_of_global, atoms
        )
        self.urey_bradleys: list[_UreyBradleyShim] = _extract_urey_bradleys(
            raw_data, local_of_global, atoms
        )
        self.nonbonded_tables = prmtop_reader.load_nonbonded_exceptions(
            raw_data, atom_indices
        )


# ====================================================================== #
# Component (D) additions -- fast-loader Step 4b (docs/specs/
# fast-amber-loader.md): box/periodicity from IFBOX + the rst7 box line
# (replaces ``parm.box_vectors``), and virtual-site (extra point) frame
# extraction from raw bond/angle arrays (replaces ``ExtraPoint``/
# ``ThreeParticleExtraPointFrame`` attribute access on a loaded ``parm``).
# ====================================================================== #

# ParmEd's ``parmed.constants.TINY`` (1e-8 Angstrom) used to snap near-zero
# box-vector components to exactly 0.0, scaled to nm.
_BOX_TINY_NM: float = 1e-9


def box_vectors_from_lengths_angles(
    lengths_nm: np.ndarray, angles_rad: np.ndarray
) -> np.ndarray:
    """Reduce ``(a, b, c, alpha, beta, gamma)`` box geometry to 3 lattice
    vectors, nm.

    Replicates ``parmed.geometry.box_lengths_and_angles_to_vectors`` exactly
    (verified bit-for-bit against it for orthorhombic and synthetic
    truncated-octahedral geometry; see docs/specs/fast-amber-loader.md Step
    4b): the standard lower-triangular reduction (``a`` along x; ``b`` in the
    xy-plane; ``c`` completes the parallelepiped), with near-zero components
    snapped to exactly ``0.0``. The formula depends only on the lengths and
    angles -- not on ``IFBOX`` -- so it handles any valid box, including
    truncated-octahedral (``IFBOX == 2``, ``alpha == beta == gamma ~=
    109.4712 deg``) with no special-casing.

    Parameters
    ----------
    lengths_nm : (3,) array-like
        ``[a, b, c]`` box vector lengths, nm.
    angles_rad : (3,) array-like
        ``[alpha, beta, gamma]``, radians.

    Returns
    -------
    (3, 3) ndarray
        Rows are the 3 reduced lattice vectors ``[a_vec; b_vec; c_vec]``, nm.
    """
    a, b, c = (float(x) for x in lengths_nm)
    alpha, beta, gamma = (float(x) for x in angles_rad)

    bx = b * math.cos(gamma)
    by = b * math.sin(gamma)
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(c * c - cx * cx - cy * cy)

    if abs(bx) < _BOX_TINY_NM:
        bx = 0.0
    if abs(by) < _BOX_TINY_NM:
        by = 0.0
    if abs(cx) < _BOX_TINY_NM:
        cx = 0.0
    if abs(cy) < _BOX_TINY_NM:
        cy = 0.0
    if abs(cz) < _BOX_TINY_NM:
        cz = 0.0

    return np.array([[a, 0.0, 0.0], [bx, by, 0.0], [cx, cy, cz]], dtype=np.float64)


@dataclass(frozen=True, slots=True)
class VirtualSiteRecord:
    """One 3-particle-average virtual site (extra point), 0-based GLOBAL
    prmtop atom indices."""

    site: int
    """The extra point's own prmtop atom index."""
    atom1: int
    """Parent (heavy) atom the EP is bonded to."""
    atom2: int
    atom3: int
    """The other two atoms bonded to the parent (the frame's 2nd/3rd points)."""
    weight1: float
    weight2: float
    weight3: float


def extract_virtual_sites(
    raw_data: dict, coords_nm: np.ndarray
) -> list[VirtualSiteRecord]:
    """Extract 3-particle-average virtual sites (4-point-water extra points:
    OPC/TIP4P/TIP4P-Ew) directly from raw prmtop bond/angle arrays + the
    instance coordinates -- replaces ParmEd's ``ExtraPoint.frame_type`` /
    ``ThreeParticleExtraPointFrame.get_weights()`` (docs/specs/
    fast-amber-loader.md Step 4b).

    An atom is a virtual site iff ``ATOMIC_NUMBER == 0`` (ParmEd's own
    ``AmberParm._load_atoms_and_residues`` criterion --
    ``parmed/amber/_amberparm.py``). AMBER represents a rigid 4-point water
    with 3 "bonds" off the oxygen (O-H1, O-H2, O-EP) and, because the whole
    O/H1/H2 triangle is rigid, usually no explicit H1-O-H2 ``ANGLE`` record
    -- the H1-H2 distance is instead carried as a third "bond". This mirrors
    ParmEd's ``ThreeParticleExtraPointFrame.get_weights()`` exactly:

    * If an ``a1-parent-a2`` angle IS present, the 2-3 distance is derived
      from it via the law of cosines.
    * Otherwise, a direct ``a1-a2`` bond supplies the 2-3 distance (the
      rigid-triangle representation above).
    * ``weight = req(parent-EP) / sqrt(req(parent-a1) * req(parent-a2) -
      0.25 * req23**2)``; "inside" (EP between the frame atoms, e.g. TIP4P)
      vs "outside" is read from the instance geometry (the angle between
      ``a1 - parent`` and ``EP - parent`` at the loaded coordinates), exactly
      as ParmEd's ``ExtraPoint.frame_type`` does when coordinates are
      available (always true here -- ``load_amber`` always has the rst7).

    Any other bonding pattern (the EP bonded to something other than exactly
    one parent, or a parent with a bond count other than 3 -- e.g. a
    2-particle frame or a 5-point out-of-plane frame like TIP5P) raises
    :class:`UnsupportedTopologyFeature` naming the pattern, per the
    fast-loader §4a parsed-but-raises policy -- it is never silently
    approximated.

    Parameters
    ----------
    raw_data : dict
        Whole-system ``raw_data`` from :func:`prmtop_reader.parse_prmtop`.
    coords_nm : np.ndarray
        Whole-system ``(N, 3)`` coordinate array, nm, prmtop atom order (from
        :func:`read_amber_coordinates`).

    Returns
    -------
    list[VirtualSiteRecord]
        One record per virtual site, in ascending prmtop atom-index order (0
        or more; empty for systems with no extra points).
    """
    atomic_number = np.asarray(raw_data.get("ATOMIC_NUMBER", []), dtype=np.int64)
    ep_indices = np.flatnonzero(atomic_number == 0)
    if ep_indices.size == 0:
        return []  # fast path: the overwhelming majority of systems have no EPs

    natom = atomic_number.shape[0]
    equil = np.asarray(raw_data["BOND_EQUIL_VALUE"], dtype=np.float64)  # raw (A)

    # Per-atom bond adjacency, built in ParmEd's own ``Atom.bonds`` insertion
    # order: BONDS_WITHOUT_HYDROGEN rows first, then BONDS_INC_HYDROGEN rows,
    # each section in file (prmtop) row order -- matches ``_extract_bonds``'s
    # convention above and is required to reproduce ParmEd's "first bond not
    # containing the EP" geometry probe below exactly.
    bonds_by_atom: list[list[tuple[int, float]]] = [[] for _ in range(natom)]
    for section in ("BONDS_WITHOUT_HYDROGEN", "BONDS_INC_HYDROGEN"):
        raw = np.asarray(raw_data.get(section, []), dtype=np.int64).reshape(-1, 3)
        for gi3, gj3, t in raw.tolist():
            gi, gj = gi3 // 3, gj3 // 3
            req = float(equil[t - 1])
            bonds_by_atom[gi].append((gj, req))
            bonds_by_atom[gj].append((gi, req))

    # Angle lookup: (vertex_atom, frozenset({outer1, outer2})) -> equilibrium
    # angle, radians (raw prmtop units -- ANGLE_EQUIL_VALUE is not degree-
    # converted by parse_prmtop, unlike CHARMM_IMPROPER_PHASE elsewhere).
    angle_theta: dict[tuple[int, frozenset], float] = {}
    ang_equil = np.asarray(raw_data.get("ANGLE_EQUIL_VALUE", []), dtype=np.float64)
    for section in ("ANGLES_WITHOUT_HYDROGEN", "ANGLES_INC_HYDROGEN"):
        raw = np.asarray(raw_data.get(section, []), dtype=np.int64).reshape(-1, 4)
        for gi3, gj3, gk3, t in raw.tolist():
            gi, gj, gk = gi3 // 3, gj3 // 3, gk3 // 3  # gj: central (vertex) atom
            angle_theta[(gj, frozenset((gi, gk)))] = float(ang_equil[t - 1])

    records: list[VirtualSiteRecord] = []
    for ep in ep_indices.tolist():
        ep_bonds = bonds_by_atom[ep]
        if len(ep_bonds) != 1:
            raise UnsupportedTopologyFeature(
                f"Virtual site at prmtop atom index {ep} (0-based) has "
                f"{len(ep_bonds)} bonds; only a single-bond (one parent) "
                "extra point is supported."
            )
        parent, req_ep = ep_bonds[0]
        parent_bonds = bonds_by_atom[parent]
        if len(parent_bonds) != 3:
            raise UnsupportedTopologyFeature(
                f"Virtual site at prmtop atom index {ep} (0-based): parent "
                f"atom {parent} has {len(parent_bonds)} bonds (frame type "
                f"needs {len(parent_bonds)} particles). Only the 3-bond "
                "in-plane average frame (TIP4P/OPC-style rigid 4-point "
                "water) is supported -- 2-particle and out-of-plane "
                "(e.g. TIP5P, 5 particles) virtual-site frames are not "
                "implemented."
            )

        others = [(a, r) for (a, r) in parent_bonds if a != ep]
        (a1, req_pa1), (a2, req_pa2) = others
        other_atom_for_inside = others[0][0]

        theta = angle_theta.get((parent, frozenset((a1, a2))))
        if theta is not None:
            req23 = math.sqrt(
                req_pa1 * req_pa1 + req_pa2 * req_pa2
                - 2.0 * req_pa1 * req_pa2 * math.cos(theta)
            )
        else:
            req23 = None
            for nbr, r in bonds_by_atom[a1]:
                if nbr == a2:
                    req23 = r
                    break
            if req23 is None:
                raise UnsupportedTopologyFeature(
                    f"Virtual site at prmtop atom index {ep} (0-based): "
                    f"cannot determine frame geometry -- no angle "
                    f"{a1}-{parent}-{a2} and no direct bond {a1}-{a2} in "
                    "the topology."
                )

        weight = req_ep / math.sqrt(req_pa1 * req_pa2 - 0.25 * req23 * req23)

        p = coords_nm[parent]
        o = coords_nm[other_atom_for_inside]
        e = coords_nm[ep]
        v1 = o - p
        v2 = e - p
        cos_ang = float(np.dot(v1, v2)) / (
            float(np.linalg.norm(v1)) * float(np.linalg.norm(v2))
        )
        inside = math.acos(max(-1.0, min(1.0, cos_ang))) < (math.pi / 2.0)

        if inside:
            w_parent, w_other = 1.0 - weight, weight / 2.0
        else:
            w_parent, w_other = 1.0 + weight, -weight / 2.0

        records.append(
            VirtualSiteRecord(
                site=ep,
                atom1=parent,
                atom2=a1,
                atom3=a2,
                weight1=w_parent,
                weight2=w_other,
                weight3=w_other,
            )
        )

    return records
