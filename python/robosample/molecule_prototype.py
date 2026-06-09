"""molecule_prototype.py

Pre-processed, simulation-ready snapshot of a single ParmEd molecule.

Design notes
------------
* All physical quantities are converted to SI-adjacent units on construction
  (lengths: nm, energies: kJ/mol, angles: rad).
* Atom, bond, angle, and torsion arrays follow *ParmEd order* -- no BFS
  reordering.  The Z-matrix is the only quantity stored in spanning-tree
  traversal order.
* Dihedral types are classified once per bond and cached in
  ``_dihedral_type_cache``.  Every subsequent lookup (``bonds_dihedral_type``,
  ``dihedral_bond_records``, ``_tree_neighbors``) reads from this cache.
* The Z-matrix lists are *uniform length n* (number of atoms).  Rows that
  have fewer than four references store the sentinel value ``-1`` in the
  inapplicable positions.
"""

from __future__ import annotations

import logging
from collections import deque
from functools import cached_property
from typing import Any

import mdtraj as md
import networkx as nx
import numpy as np
import parmed as pmd

from .amber_dihedral_classifier import AmberDihedralClassifier
from .amber_dihedral_types import DihedralType
from .bond_util import is_rigid_bond
from .secondary_structure import BondDSSP, DSSPCode

logger: logging.Logger = logging.getLogger(__name__)

# Sentinel stored in z_matrix_j/k/l when a reference atom does not exist.
_ZM_SENTINEL: int = -1


class MoleculePrototype:
    """
    Immutable pre-processed snapshot of a single ParmEd molecule.

    All mutable state is set in ``__init__``.  Expensive derived quantities
    (``dihedral_bond_records``) are computed lazily via ``cached_property``.
    """

    # ------------------------------------------------------------------
    # Unit conversion factors (class-level constants)
    # ------------------------------------------------------------------

    KCAL_TO_KJ: float = pmd.unit.kilocalories_per_mole.conversion_factor_to(
        pmd.unit.kilojoules_per_mole
    )
    ANG_TO_NM: float = pmd.unit.angstrom.conversion_factor_to(pmd.unit.nanometer)
    DEG_TO_RAD: float = pmd.unit.degree.conversion_factor_to(pmd.unit.radian)

    # Converts the AMBER r_min (equilibrium pair distance) to the LJ sigma:
    #   sigma = r_min / 2^(1/6)   =>   SIGMA_SCALE = 2^(-1/6) ~= 0.8909
    SIGMA_SCALE: float = 2.0 ** (-1.0 / 6.0)

    # ------------------------------------------------------------------
    # Attribute declarations (values set in __init__)
    # ------------------------------------------------------------------

    num_residues: int
    """Number of residues in the molecule."""

    num_atoms: int
    """Total number of atoms."""

    atoms_root_index: int
    """
    ParmEd atom index of the Z-matrix root (heaviest terminal atom, or the
    sole atom for n == 1).
    """

    # Atom arrays -- ParmEd order, unit-converted
    atoms_mass: list[float]
    """Mass in daltons [Da]."""

    atoms_charge: list[float]
    """Partial charge in units of the proton charge."""

    atoms_sigma: list[float]
    """Lennard-Jones sigma (van der Waals radius) in nm."""

    atoms_epsilon: list[float]
    """Lennard-Jones epsilon (well depth) in kJ/mol."""

    atoms_radius: list[float]
    """GBSA solvent radius in nm."""

    atoms_screen: list[float]
    """OBC screening factor (dimensionless)."""

    atoms_x: list[float]
    """x coordinate from the reference structure in nm."""

    atoms_y: list[float]
    """y coordinate from the reference structure in nm."""

    atoms_z: list[float]
    """z coordinate from the reference structure in nm."""

    # Bond arrays -- ParmEd order
    num_bonds: int
    """Total number of bonds (tree bonds + ring-closing bonds)."""

    bonds_i: list[int]
    """ParmEd atom index of bond endpoint atom1."""

    bonds_j: list[int]
    """ParmEd atom index of bond endpoint atom2."""

    bonds_stiffness: list[float]
    """Harmonic force constant in kJ/mol/nm^2."""

    bonds_equilibrium: list[float]
    """Equilibrium bond length in nm."""

    bonds_is_ring_closing: list[bool]
    """True when the bond was removed from the spanning tree (ring closure)."""

    bonds_secondary_structure: list[DSSPCode]
    """DSSP secondary structure code of the bond's first endpoint atom."""

    bonds_dihedral_type: list[DihedralType]
    """
    Best DihedralType found for any dihedral whose *middle bond* is this bond.
    ``DihedralType.UNKNOWN`` when no dihedral uses this bond as its middle bond
    or when the classifier returned UNKNOWN for all candidates.
    """

    # Angle arrays -- ParmEd order
    num_angles: int
    angles_i: list[int]
    angles_j: list[int]
    angles_k: list[int]
    angles_equilibrium: list[float]
    """Equilibrium angle in rad."""
    angles_stiffness: list[float]
    """Harmonic force constant in kJ/mol/rad^2."""

    # Periodic torsion arrays -- ParmEd order, expanded per term
    num_periodic_torsions: int
    periodic_torsions_improper: list[bool]
    periodic_torsions_i: list[int]
    periodic_torsions_j: list[int]
    periodic_torsions_k: list[int]
    periodic_torsions_l: list[int]
    periodic_torsions_n: list[int]
    """Periodicity."""
    periodic_torsions_phase: list[float]
    """Phase offset in rad."""
    periodic_torsions_stiffness: list[float]
    """Force constant in kJ/mol."""

    # Harmonic (CHARMM-style improper) torsion arrays -- ParmEd order
    num_harmonic_torsions: int
    harmonic_torsions_i: list[int]
    harmonic_torsions_j: list[int]
    harmonic_torsions_k: list[int]
    harmonic_torsions_l: list[int]
    harmonic_torsions_stiffness: list[float]
    """Force constant in kJ/mol."""
    harmonic_torsions_phase: list[float]
    """Equilibrium angle in rad."""

    # Urey-Bradley terms -- ParmEd order
    num_urey_bradley: int
    """Number of Urey-Bradley (1,3-distance) terms."""

    urey_bradley_i: list[int]
    """Atom index of the first 1,3-end atom of each Urey-Bradley term."""

    urey_bradley_k: list[int]
    """Atom index of the second 1,3-end atom of each Urey-Bradley term."""

    urey_bradley_stiffness: list[float]
    """Harmonic force constant in kJ/mol/nm^2."""

    urey_bradley_equilibrium: list[float]
    """Equilibrium 1,3-distance in nm."""

    # 1-4 non-bonded scaling pairs -- parsed from AMBER prmtop dihedral pointers
    num_scaling14: int
    """Number of unique 1-4 non-bonded pair records."""

    scaling14_i: list[int]
    """Atom index of the first atom (i, position 1 in the originating dihedral)."""

    scaling14_l: list[int]
    """Atom index of the fourth atom (l, position 4 in the originating dihedral)."""

    scaling14_charge_product: list[float]
    """
    Electrostatic charge product q_i * q_l / scee in units of the proton charge
    squared, where scee is the per-dihedral-type SCEE_SCALE_FACTOR from the
    prmtop.  The SCEE factor is already absorbed; no further scaling is needed.
    """

    scaling14_epsilon: list[float]
    """
    Combined LJ well depth for the pair in kJ/mol, derived from the prmtop
    LENNARD_JONES_14_ACOEF / BCOEF tables as Bcoef^2 / (4*Acoef), then divided
    by the per-dihedral-type SCNB_SCALE_FACTOR.  The SCNB factor is already
    absorbed.
    """

    scaling14_sigma: list[float]
    """
    Combined LJ sigma for the pair in nm, derived from the prmtop LJ coefficients
    as (2*Acoef/Bcoef)^(1/6) * SIGMA_SCALE, where SIGMA_SCALE = 2^(-1/6)
    converts the AMBER r_min to the LJ sigma.
    """

    # Exclusions -- parsed from AMBER prmtop EXCLUDED_ATOMS_LIST
    num_exclusions: int
    """Number of explicit non-bonded exclusion pairs (not already in scaling14)."""

    exclusion_i: list[int]
    """Atom index of the first atom of each exclusion pair."""

    exclusion_j: list[int]
    """Atom index of the second atom of each exclusion pair."""

    # Z-matrix arrays -- spanning-tree traversal order, length == num_atoms
    z_matrix_i: list[int]
    """
    Global atom index of the atom placed at row r.
    Length n.  Rows are in tree-traversal order rooted at atoms_root_index.
    """

    z_matrix_j: list[int]
    """
    Bond-length reference atom index at row r.
    Length n.  Row 0 stores the sentinel -1 (root has no bond reference).
    """

    z_matrix_k: list[int]
    """
    Bond-angle reference atom index at row r.
    Length n.  Rows 0-1 store the sentinel -1.
    """

    z_matrix_l: list[int]
    """
    Dihedral reference atom index at row r.
    Length n.  Rows 0-2 store the sentinel -1.
    """

    # ------------------------------------------------------------------
    # Internal attributes
    # ------------------------------------------------------------------

    acyclic_graph: nx.Graph
    """Spanning forest of the molecule (ring-closing bonds removed)."""

    _atom_by_idx: dict[int, pmd.Atom]
    """O(1) atom lookup by ParmEd atom index.  Built once in __init__."""

    _dihedral_type_cache: dict[frozenset[int], DihedralType]
    """
    Dihedral type keyed by the frozenset of the two *middle-bond* atom indices.
    Populated by ``_build_acyclic_graph``; read by ``bonds_dihedral_type``,
    ``dihedral_bond_records``.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(
        self,
        molecule: pmd.Structure,
        dihedral_classifier: AmberDihedralClassifier,
    ) -> None:
        self.molecule = molecule
        self.dihedral_classifier = dihedral_classifier

        self.num_atoms = len(molecule.atoms)
        self.num_residues = max((a.residue.idx + 1 for a in molecule.atoms), default=0)

        # O(n) lookup; built before any method that needs it
        self._atom_by_idx = {a.idx: a for a in molecule.atoms}

        # Dihedral type cache is populated as a side-effect of _build_acyclic_graph
        self._dihedral_type_cache = {}
        self.acyclic_graph = self._build_acyclic_graph()

        # Root: heaviest terminal atom (or the sole atom for n == 1)
        if self.num_atoms == 1:
            root_atom: pmd.Atom = molecule.atoms[0]
        else:
            terminals: list[pmd.Atom] = [
                a for a in molecule.atoms if len(a.bond_partners) == 1
            ]
            root_atom = self._sort_atoms_by_mass(terminals)[0]
        self.atoms_root_index = root_atom.idx

        logger.debug(
            "Root atom: idx=%d  name=%s  mass=%.3f",
            root_atom.idx,
            getattr(root_atom, "name", "?"),
            getattr(root_atom, "mass", 0.0),
        )

        # -- Atoms (ParmEd order, unit-converted) --------------------------
        self.atoms_mass = [a.mass for a in molecule.atoms]
        self.atoms_charge = [a.charge for a in molecule.atoms]
        self.atoms_sigma = [a.sigma * self.ANG_TO_NM for a in molecule.atoms]
        self.atoms_epsilon = [a.epsilon * self.KCAL_TO_KJ for a in molecule.atoms]
        self.atoms_radius = [a.solvent_radius * self.ANG_TO_NM for a in molecule.atoms]
        self.atoms_screen = [a.screen for a in molecule.atoms]
        self.atoms_x = [a.xx * self.ANG_TO_NM for a in molecule.atoms]
        self.atoms_y = [a.xy * self.ANG_TO_NM for a in molecule.atoms]
        self.atoms_z = [a.xz * self.ANG_TO_NM for a in molecule.atoms]

        # -- Bonds (ParmEd order) ------------------------------------------

        self.num_bonds = len(molecule.bonds)
        self.bonds_i = [b.atom1.idx for b in molecule.bonds]
        self.bonds_j = [b.atom2.idx for b in molecule.bonds]
        self.bonds_stiffness = [
            b.type.k * self.KCAL_TO_KJ / self.ANG_TO_NM**2 for b in molecule.bonds
        ]
        self.bonds_equilibrium = [b.type.req * self.ANG_TO_NM for b in molecule.bonds]

        acyclic_edges: set[frozenset[int]] = {
            frozenset(e) for e in self.acyclic_graph.edges()
        }
        self.bonds_is_ring_closing = [
            frozenset((b.atom1.idx, b.atom2.idx)) not in acyclic_edges
            for b in molecule.bonds
        ]
        self.bonds_dihedral_type = [
            self._dihedral_type_cache.get(
                frozenset((b.atom1.idx, b.atom2.idx)), DihedralType.UNKNOWN
            )
            for b in molecule.bonds
        ]

        # -- Bonds DSSP (ParmEd order) -------------------------------------
        traj = md.Trajectory(
            xyz=np.array([self.atoms_x, self.atoms_y, self.atoms_z]).T,
            topology=md.Topology.from_openmm(molecule.topology),
        )

        dssp = md.compute_dssp(traj, simplified=False)
        dssp = dssp[0]
        self.bonds_secondary_structure = [
            BondDSSP(
                atom1_code=DSSPCode.from_mdtraj(dssp[b.atom1.residue.idx]),
                atom2_code=DSSPCode.from_mdtraj(dssp[b.atom2.residue.idx]),
            ).resolve(strategy="priority")
            for b in molecule.bonds
        ]

        # -- Angles (ParmEd order) -----------------------------------------
        self.num_angles = len(molecule.angles)
        self.angles_i = [a.atom1.idx for a in molecule.angles]
        self.angles_j = [a.atom2.idx for a in molecule.angles]
        self.angles_k = [a.atom3.idx for a in molecule.angles]
        self.angles_equilibrium = [
            a.type.theteq * self.DEG_TO_RAD for a in molecule.angles
        ]
        self.angles_stiffness = [a.type.k * self.KCAL_TO_KJ for a in molecule.angles]

        # -- Periodic torsions (ParmEd order, one record per term) ---------
        periodic_records: list[tuple[pmd.Dihedral, pmd.DihedralType]] = []
        for d in molecule.dihedrals:
            if isinstance(d.type, pmd.DihedralTypeList):
                periodic_records.extend((d, dt) for dt in d.type)
            elif isinstance(d.type, pmd.DihedralType):
                periodic_records.append((d, d.type))

        self.num_periodic_torsions = len(periodic_records)
        self.periodic_torsions_improper = [d.improper for d, dt in periodic_records]
        self.periodic_torsions_i = [d.atom1.idx for d, dt in periodic_records]
        self.periodic_torsions_j = [d.atom2.idx for d, dt in periodic_records]
        self.periodic_torsions_k = [d.atom3.idx for d, dt in periodic_records]
        self.periodic_torsions_l = [d.atom4.idx for d, dt in periodic_records]
        self.periodic_torsions_n = [dt.per for d, dt in periodic_records]
        self.periodic_torsions_phase = [
            dt.phase * self.DEG_TO_RAD for d, dt in periodic_records
        ]
        self.periodic_torsions_stiffness = [
            dt.phi_k * self.KCAL_TO_KJ for d, dt in periodic_records
        ]

        # -- Harmonic (CHARMM-style improper) torsions (ParmEd order) ------
        self.num_harmonic_torsions = len(molecule.impropers)
        self.harmonic_torsions_i = [imp.atom1.idx for imp in molecule.impropers]
        self.harmonic_torsions_j = [imp.atom2.idx for imp in molecule.impropers]
        self.harmonic_torsions_k = [imp.atom3.idx for imp in molecule.impropers]
        self.harmonic_torsions_l = [imp.atom4.idx for imp in molecule.impropers]
        self.harmonic_torsions_stiffness = [
            imp.type.psi_k * self.KCAL_TO_KJ for imp in molecule.impropers
        ]
        self.harmonic_torsions_phase = [
            imp.type.psi_eq * self.DEG_TO_RAD for imp in molecule.impropers
        ]

        # -- Urey-Bradley terms (ParmEd order) --------------------------------
        # Each UreyBradley connects the two end-atoms of an angle (1,3 pair).
        # The type carries a BondType-style k and req (equilibrium 1,3 distance).
        self.num_urey_bradley = len(molecule.urey_bradleys)
        self.urey_bradley_i = [ub.atom1.idx for ub in molecule.urey_bradleys]
        self.urey_bradley_k = [ub.atom2.idx for ub in molecule.urey_bradleys]
        self.urey_bradley_stiffness = [
            ub.type.k * self.KCAL_TO_KJ / self.ANG_TO_NM**2
            for ub in molecule.urey_bradleys
        ]
        self.urey_bradley_equilibrium = [
            ub.type.req * self.ANG_TO_NM for ub in molecule.urey_bradleys
        ]

        # -- 1-4 scaling pairs and exclusions ---------------------------------
        # These are not available on the ParmEd Structure API; they must be
        # parsed from the raw AMBER prmtop section data.  The arrays are
        # initialised empty here and populated by _parse_amber_prmtop_nonbonded
        # when the molecule carries parm_data (i.e. was loaded from a .prmtop).
        self.num_scaling14 = 0
        self.scaling14_i = []
        self.scaling14_l = []
        self.scaling14_charge_product = []
        self.scaling14_epsilon = []
        self.scaling14_sigma = []

        self.num_exclusions = 0
        self.exclusion_i = []
        self.exclusion_j = []

        if hasattr(molecule, "parm_data"):
            self._parse_amber_prmtop_nonbonded()
        else:
            logger.warning(
                "molecule has no parm_data attribute; scaling14 and exclusions "
                "will be empty.  Load the molecule from an AMBER .prmtop file "
                "to populate these arrays."
            )

        # -- Z-matrix -------------------------------------------------------
        self._build_z_matrix(root_atom)

    # ------------------------------------------------------------------
    # Public helpers
    # ------------------------------------------------------------------

    @cached_property
    def dihedral_bond_records(self) -> list[dict[str, Any]]:
        """
        One record per unique proper-torsion central bond.

        Ring-closing central bonds and improper torsions are excluded.
        Atom indices are raw ParmEd indices (no BFS reordering).

        Each record contains:
          j             -- atom index of the first  central-bond atom
          k             -- atom index of the second central-bond atom
          dihedral_type -- ``DihedralType`` enum value
          resid         -- residue index of atom j
          resname       -- residue name  of atom j
        """
        ring_closing_bonds: set[frozenset[int]] = {
            frozenset((self.bonds_i[n], self.bonds_j[n]))
            for n in range(self.num_bonds)
            if self.bonds_is_ring_closing[n]
        }

        seen: dict[frozenset[int], dict[str, Any]] = {}
        for d in self.molecule.dihedrals:
            if d.improper:
                continue
            j_idx: int = d.atom2.idx
            k_idx: int = d.atom3.idx
            key: frozenset[int] = frozenset((j_idx, k_idx))
            if key in ring_closing_bonds:
                continue
            dtype: DihedralType = self._dihedral_type_cache.get(
                key, DihedralType.UNKNOWN
            )
            if key not in seen or (
                seen[key]["dihedral_type"] is DihedralType.UNKNOWN
                and dtype is not DihedralType.UNKNOWN
            ):
                seen[key] = {
                    "j": j_idx,
                    "k": k_idx,
                    "dihedral_type": dtype,
                    "resid": d.atom2.residue.idx,
                    "resname": d.atom2.residue.name,
                }
        return list(seen.values())

    @staticmethod
    def _sort_atoms_by_mass(atoms: list[pmd.Atom]) -> list[pmd.Atom]:
        """
        Return a new list sorted by descending mass.  Atom index breaks ties
        (ascending) to guarantee a reproducible ordering.
        """
        return sorted(atoms, key=lambda a: (-a.mass, a.idx))

    def _tree_neighbors(self, atom: pmd.Atom) -> list[pmd.Atom]:
        """
        Return the tree-neighbors of *atom* (atoms bonded via non-ring-closing
        bonds).  Uses ``acyclic_graph`` for O(degree) lookup.
        """
        return [self._atom_by_idx[nb] for nb in self.acyclic_graph.neighbors(atom.idx)]

    # ------------------------------------------------------------------
    # AMBER prmtop non-bonded section parser
    # ------------------------------------------------------------------

    def _parse_amber_prmtop_nonbonded(self) -> None:
        """
        Populate ``scaling14_*`` and ``exclusions_*`` arrays from raw AMBER
        prmtop section data stored in ``self.molecule.parm_data``.

        This method must only be called when ``hasattr(molecule, 'parm_data')``
        is True (i.e. the molecule is a ``pmd.amber.AmberParm``).

        1-4 scaling pairs
        -----------------
        Parsed from ``DIHEDRALS_INC_HYDROGEN`` + ``DIHEDRALS_WITHOUT_HYDROGEN``
        (groups of five integers: i, j, k, l, dihedral_type_index).

        * Entries with ``k < 0`` signal that the 1-4 non-bonded interaction for
          this dihedral is excluded (the pair is already handled elsewhere); they
          are skipped.
        * Entries with ``l < 0`` are improper dihedrals; they are skipped.
        * LJ parameters are taken from ``LENNARD_JONES_14_ACOEF / BCOEF`` when
          those sections exist (AMBER ff14SB and later), falling back to
          ``LENNARD_JONES_ACOEF / BCOEF`` otherwise.
        * ``SCEE_SCALE_FACTOR`` and ``SCNB_SCALE_FACTOR`` (per dihedral type)
          are absorbed into ``scaling14_charge_product`` and ``scaling14_epsilon``
          respectively so that callers can use the values directly.
        * Duplicate pairs (same atom pair reached via different dihedrals) are
          deduplicated; the first occurrence wins.

        Exclusions
        ----------
        Parsed from ``NUMBER_EXCLUDED_ATOMS`` + ``EXCLUDED_ATOMS_LIST``.
        Pairs that are already in the scaling14 set are excluded here so that
        a given atom pair appears in exactly one of the two lists.
        """
        molecule = self.molecule
        parm_data: dict = molecule.parm_data  # type: ignore[attr-defined]

        num_types: int = molecule.ptr("NTYPES")  # type: ignore[attr-defined]

        # Prefer 1-4-specific LJ tables; fall back to standard tables.
        lj14_a: list[float] = parm_data.get(
            "LENNARD_JONES_14_ACOEF", parm_data["LENNARD_JONES_ACOEF"]
        )
        lj14_b: list[float] = parm_data.get(
            "LENNARD_JONES_14_BCOEF", parm_data["LENNARD_JONES_BCOEF"]
        )
        nb_index: list[int] = parm_data["NONBONDED_PARM_INDEX"]
        charges: list[float] = parm_data["CHARGE"]  # AMBER internal units
        scee_factors: list[float] = parm_data["SCEE_SCALE_FACTOR"]
        scnb_factors: list[float] = parm_data["SCNB_SCALE_FACTOR"]

        # AMBER stores charges as q * 18.2223 (sqrt(332.0636) in kcal*Ang/e^2).
        # ParmEd exposes atom.charge already converted to proton charge units,
        # but the parm_data["CHARGE"] array retains the raw AMBER units.
        # We derive charge_product from atom.charge (already in proton charges)
        # so no manual conversion is needed.

        seen_14: set[tuple[int, int]] = set()

        dihedral_ptrs: list[int] = (
            parm_data["DIHEDRALS_INC_HYDROGEN"]
            + parm_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        )

        for ii in range(0, len(dihedral_ptrs), 5):
            i_raw, _j_raw, k_raw, l_raw, dtype_idx = dihedral_ptrs[ii : ii + 5]

            # k < 0: 1-4 interaction suppressed (ring / already excluded)
            if k_raw < 0:
                continue
            # l < 0: improper dihedral
            if l_raw < 0:
                continue

            atom_i_idx: int = i_raw // 3
            atom_l_idx: int = l_raw // 3

            # Safety: atoms must belong to this molecule
            if (
                atom_i_idx not in self._atom_by_idx
                or atom_l_idx not in self._atom_by_idx
            ):
                continue

            atom_i: pmd.Atom = self._atom_by_idx[atom_i_idx]
            atom_l: pmd.Atom = self._atom_by_idx[atom_l_idx]

            # Canonical key: (min, max) so (i,l) and (l,i) are the same pair
            key: tuple[int, int] = (
                min(atom_i_idx, atom_l_idx),
                max(atom_i_idx, atom_l_idx),
            )
            if key in seen_14:
                continue

            # LJ pair index (0-based) from the NONBONDED_PARM_INDEX table
            nb_i: int = atom_i.nb_idx - 1  # ParmEd nb_idx is 1-based
            nb_l: int = atom_l.nb_idx - 1
            pair_idx: int = nb_index[nb_i * num_types + nb_l] - 1  # 0-based

            if pair_idx < 0:
                continue

            acoef: float = lj14_a[pair_idx]
            bcoef: float = lj14_b[pair_idx]

            if acoef != 0.0 and bcoef != 0.0:
                # epsilon = Bcoef^2 / (4 * Acoef)   [kcal/mol]
                epsilon_kcal: float = (bcoef**2) / (4.0 * acoef)
                # r_min = (2*Acoef/Bcoef)^(1/6)     [Angstrom]
                r_min_ang: float = (2.0 * acoef / bcoef) ** (1.0 / 6.0)
                epsilon_kj: float = epsilon_kcal * self.KCAL_TO_KJ
                sigma_nm: float = r_min_ang * self.ANG_TO_NM * self.SIGMA_SCALE
            else:
                # Zero LJ coefficients: no vdW interaction for this pair
                epsilon_kj = 0.0
                sigma_nm = 1.0 * self.ANG_TO_NM  # placeholder; epsilon is 0

            scee: float = scee_factors[dtype_idx - 1]
            scnb: float = scnb_factors[dtype_idx - 1]

            # atom.charge is already in proton-charge units (ParmEd converts it)
            charge_product: float = atom_i.charge * atom_l.charge / scee
            epsilon_scaled: float = epsilon_kj / scnb

            seen_14.add(key)
            self.scaling14_i.append(key[0])
            self.scaling14_l.append(key[1])
            self.scaling14_charge_product.append(charge_product)
            self.scaling14_epsilon.append(epsilon_scaled)
            self.scaling14_sigma.append(sigma_nm)

        self.num_scaling14 = len(self.scaling14_i)
        logger.debug("Parsed %d 1-4 scaling pairs from prmtop.", self.num_scaling14)

        # ---- Exclusions -------------------------------------------------------
        n_excluded_list: list[int] = parm_data["NUMBER_EXCLUDED_ATOMS"]
        excluded_atoms: list[int] = parm_data["EXCLUDED_ATOMS_LIST"]

        # Start with the 1-4 pairs already handled above so the two lists are
        # mutually exclusive.
        seen_excl: set[tuple[int, int]] = set(seen_14)

        offset: int = 0
        for i_atom in range(len(molecule.atoms)):
            n: int = int(n_excluded_list[i_atom])
            for j_atom_1based in excluded_atoms[offset : offset + n]:
                j: int = int(j_atom_1based)
                if j <= 0:
                    # prmtop uses j=0 as a placeholder for atoms with no exclusions
                    continue
                j_idx: int = j - 1  # prmtop is 1-based

                if j_idx not in self._atom_by_idx:
                    continue

                key = (min(i_atom, j_idx), max(i_atom, j_idx))
                if key in seen_excl:
                    continue
                seen_excl.add(key)
                self.exclusion_i.append(key[0])
                self.exclusion_j.append(key[1])

            offset += n

        self.num_exclusions = len(self.exclusion_i)
        logger.debug("Parsed %d exclusions from prmtop.", self.num_exclusions)

    def _build_acyclic_graph(self) -> nx.Graph:
        """
        Build a spanning forest by removing ring-closing (cotree) bonds.

        Pipeline
        --------
        1. Build the full molecular bond graph and classify every bond's
           dihedral type, populating ``_dihedral_type_cache``.  Bonds whose
           best classification is ``PROTEIN_RING_DIHEDRAL`` are recorded as
           *candidate* ring-closing bonds.
        2. First pass -- cut each candidate ring-closing bond, but only where
           doing so is safe:
             * A bond whose removal would disconnect the molecule (a *bridge*,
               e.g. an inter-chain disulfide that is the only covalent link
               between two chains) is kept in the spanning tree and a warning
               is logged.
             * At most one closure is anchored per atom in this pass.  This is
               a *soft* preference: bonds skipped because their endpoint is
               already used are simply deferred to step 3 rather than forced.
        3. Reduce any residual cycles to a spanning forest by choosing the
           ring-closing bonds as the *complement of a maximum-weight spanning
           tree*.  ``keep_weight`` encodes (lexicographically, via separated
           magnitudes) how desirable a bond is to KEEP as a normal tree edge:

               forbidden (SG--SG) > one-closure-per-atom > flexible-over-rigid
               > near-a-forbidden-atom > higher-degree

           A spanning tree always exists for a connected graph, so this never
           dead-ends on heavily fused polycyclic systems (e.g. cucurbiturils)
           the way greedy cycle-basis breaking does, and it cuts exactly the
           minimal number of bonds (the circuit rank).  Edges that are
           topological bridges are guaranteed to remain in the tree, so
           connectivity is preserved automatically.
        4. Validate the final ring-closing bond set via
           ``_validate_ring_closing_bonds``.
        5. Run the connectivity and acyclicity checks as a final sanity gate.

        Difference from the legacy heuristic
        -------------------------------------
        Disulfide (SG--SG) bonds are *no longer* forced to be ring-closing.
        They are treated as "forbidden" bonds that are preferentially KEPT in
        the spanning tree; a disulfide is only cut when it genuinely lies on a
        cycle and cutting it is the least-bad option.  This prevents a bridging
        inter-chain disulfide from being severed (which previously disconnected
        the molecule and raised in ``_check_disconnected_graph``).

        Side effect
        -----------
        Populates ``self._dihedral_type_cache`` unconditionally for every bond.

        Returns
        -------
        nx.Graph
            Acyclic (forest) graph whose nodes are ParmEd atom indices.

        Raises
        ------
        ValueError
            If the result is still cyclic or disconnected after all removal
            steps, or if ring-closing bond validation fails.
        """
        # ----------------------------------------------------------------
        # Step 1: build full molecular bond graph + classify + cache
        # ----------------------------------------------------------------
        full_graph: nx.Graph = nx.Graph()
        full_graph.add_nodes_from(a.idx for a in self.molecule.atoms)

        # Bonds whose dihedral classification marks them as ring dihedrals.
        # These are *candidates* for cutting; the actual cut decision is made
        # in steps 2-3 so that bridges are never severed.
        ring_dihedral_bonds: set[frozenset[int]] = set()

        for bond in self.molecule.bonds:
            # Convention kept from original: atom2 is the "parent" side.
            parent: pmd.Atom = bond.atom2
            child: pmd.Atom = bond.atom1
            key: frozenset[int] = frozenset((parent.idx, child.idx))

            full_graph.add_edge(parent.idx, child.idx)

            best_type: DihedralType = DihedralType.UNKNOWN
            is_ring_dihedral: bool = False

            for grandparent in parent.bond_partners:
                if grandparent is child:
                    continue
                for gchild in child.bond_partners:
                    if gchild is parent:
                        continue
                    candidate = pmd.Dihedral(grandparent, parent, child, gchild)
                    dtype: DihedralType = self.dihedral_classifier.classify(candidate)

                    if dtype == DihedralType.PROTEIN_RING_DIHEDRAL:
                        is_ring_dihedral = True
                        best_type = dtype
                        break
                    if dtype != DihedralType.UNKNOWN:
                        best_type = dtype
                if is_ring_dihedral:
                    break

            self._dihedral_type_cache[key] = best_type

            if is_ring_dihedral:
                ring_dihedral_bonds.add(key)

        # ----------------------------------------------------------------
        # Step 2: cut classified ring dihedrals where it is safe to do so
        # ----------------------------------------------------------------
        g: nx.Graph = full_graph.copy()

        # SG--SG (disulfide) bonds are "forbidden": preferentially KEPT in the
        # spanning tree.  They carry no torsional information worth modelling,
        # and an inter-chain disulfide may be the sole covalent link between two
        # chains, so cutting it would disconnect the molecule.
        forbidden_edges: set[frozenset[int]] = self._find_forbidden_bonds(
            [("SG", "SG")]
        )

        # ring_closing_set: bonds removed from the spanning tree (cotree bonds).
        # Populated by the first pass (step 2) and the spanning step (step 3).
        ring_closing_set: set[frozenset[int]] = set()

        # Soft "one ring-closing bond per atom" bookkeeping, shared by both
        # the first pass and the maximum-spanning-tree weighting below.
        closures_per_atom: dict[int, int] = {}

        def _can_cut(u: int, v: int) -> bool:
            """
            Remove edge (u, v) from *g* iff the graph stays connected.

            Returns True and leaves the edge removed when cutting is safe (the
            edge lies on a cycle).  Returns False and restores the edge when it
            is a bridge.
            """
            g.remove_edge(u, v)
            if nx.has_path(g, u, v):
                return True
            g.add_edge(u, v)  # bridge: put it back
            return False

        # Iterate bonds in ParmEd order (not over the set) for reproducibility.
        for bond in self.molecule.bonds:
            u: int = bond.atom2.idx
            v: int = bond.atom1.idx
            key = frozenset((u, v))
            if key not in ring_dihedral_bonds:
                continue

            # Soft one-closure-per-atom: prefer not to anchor two closures on
            # the same atom.  In fused polycyclic systems this is sometimes
            # unavoidable, so DEFER rather than force -- leave the bond in the
            # graph and let the spanning-tree step decide.
            if closures_per_atom.get(u, 0) >= 1 or closures_per_atom.get(v, 0) >= 1:
                continue

            # Only cut if both endpoints remain connected through the rest of
            # the graph (i.e. the bond lies on a cycle, not a bridge).
            if not _can_cut(u, v):
                logger.warning(
                    "Bond %d -- %d is classified as ring-closing but is a bridge "
                    "(its removal would disconnect the molecule); keeping it in the "
                    "spanning tree to preserve connectivity.",
                    u,
                    v,
                )
                continue

            ring_closing_set.add(key)
            for idx in (u, v):
                closures_per_atom[idx] = closures_per_atom.get(idx, 0) + 1

        # ----------------------------------------------------------------
        # Step 3: reduce residual cycles via a maximum-weight spanning tree
        # ----------------------------------------------------------------
        forbidden_nodes: set[int] = {n for edge in forbidden_edges for n in tuple(edge)}
        # Distance from each node to the nearest forbidden atom (0 when there
        # are no forbidden atoms).  Computed on the post-first-pass graph.
        forbidden_dist: dict[int, int] = (
            nx.multi_source_shortest_path_length(g, forbidden_nodes)
            if forbidden_nodes
            else {}
        )

        def _keep_weight(u: int, v: int) -> float:
            """
            Soft preference for KEEPING bond (u, v) as a normal tree edge
            (higher == keep).  The lowest-weight edges become ring closures.
            Magnitudes are separated so the ordering is effectively
            lexicographic.
            """
            weight: float = 0.0

            # Forbidden bonds (e.g. SG--SG) should stay in the tree if possible.
            if frozenset((u, v)) in forbidden_edges:
                weight += 1e9

            # Soft "one ring-closing bond per atom" preference.
            if closures_per_atom.get(u, 0) or closures_per_atom.get(v, 0):
                weight += 1e6

            # Prefer to cut rigid bonds, i.e. keep flexible ones in the tree.
            if not is_rigid_bond(self._atom_by_idx[u], self._atom_by_idx[v]):
                weight += 1e3

            # Prefer to cut bonds far from forbidden atoms, i.e. keep near ones.
            if forbidden_dist:
                dist = max(forbidden_dist.get(u, 0), forbidden_dist.get(v, 0))
                weight += -0.1 * float(dist)

            # Mild preference to keep higher-degree bonds in the tree.
            weight += 1e-3 * (g.degree[u] + g.degree[v])

            return weight

        for u, v in g.edges:
            g[u][v]["keep_weight"] = _keep_weight(u, v)

        # Maximum spanning forest: kept edges form the tree, the rest close rings.
        spanning: nx.Graph = nx.maximum_spanning_tree(g, weight="keep_weight")
        kept_edges: set[frozenset[int]] = {frozenset(e) for e in spanning.edges}

        for u, v in list(g.edges):
            if frozenset((u, v)) in kept_edges:
                continue
            g.remove_edge(u, v)
            key = frozenset((u, v))
            ring_closing_set.add(key)
            for idx in (u, v):
                closures_per_atom[idx] = closures_per_atom.get(idx, 0) + 1

        # Report (do not fail) when fused rings force an atom into >1 closure.
        overloaded: list[int] = sorted(
            idx for idx, count in closures_per_atom.items() if count > 1
        )
        if overloaded:
            shown = ", ".join(self._atom_label(idx) for idx in overloaded[:10])
            logger.warning(
                "%d atom(s) participate in more than one ring-closing bond. This is "
                "expected for fused polycyclic systems (e.g. cucurbituril) and is "
                "treated as a soft preference rather than an error. Atoms: %s%s.",
                len(overloaded),
                shown,
                ", ..." if len(overloaded) > 10 else "",
            )

        # ----------------------------------------------------------------
        # Step 4: validate ring-closing bond set
        # ----------------------------------------------------------------
        if full_graph.number_of_nodes() != self.num_atoms:
            raise ValueError(
                "Acyclic graph node count %d != atom count %d."
                % (full_graph.number_of_nodes(), self.num_atoms)
            )
        self._validate_ring_closing_bonds(full_graph, g, ring_closing_set)

        # ----------------------------------------------------------------
        # Step 5: final sanity checks
        # ----------------------------------------------------------------
        self._check_residual_cycles(g)
        self._check_disconnected_graph(g)

        return g

    # ------------------------------------------------------------------
    # Ring-closing bond helpers
    # ------------------------------------------------------------------

    def _find_forbidden_bonds(
        self, atom_type_pairs: list[tuple[str, str]]
    ) -> set[frozenset[int]]:
        """
        Identify bonds whose two endpoint atom *names* match any forbidden
        atom-type pair.  Matching is case-insensitive and order-independent.

        Forbidden bonds are ones we prefer to KEEP in the spanning tree rather
        than turn into ring closures -- e.g. disulfides (``("SG", "SG")``),
        which carry no useful torsional degree of freedom and may be the only
        covalent link between two chains.

        Parameters
        ----------
        atom_type_pairs : list[tuple[str, str]]
            Forbidden atom-name pairs, e.g. ``[("SG", "SG")]``.

        Returns
        -------
        set[frozenset[int]]
            Forbidden bonds as frozensets of the two ParmEd atom indices.
        """
        forbidden_pairs: set[frozenset[str]] = {
            frozenset((a.lower(), b.lower())) for a, b in atom_type_pairs
        }

        forbidden_edges: set[frozenset[int]] = set()
        for bond in self.molecule.bonds:
            a1, a2 = bond.atom1, bond.atom2
            if frozenset((a1.name.lower(), a2.name.lower())) in forbidden_pairs:
                forbidden_edges.add(frozenset((a1.idx, a2.idx)))
        return forbidden_edges

    def _validate_ring_closing_bonds(
        self,
        full_graph: nx.Graph,
        spanning_tree: nx.Graph,
        ring_closing_set: set[frozenset[int]],
    ) -> None:
        """
        Validate that the selected ring-closing bonds are consistent.

        Checks
        ------
        1. **Cyclomatic number** -- ``|ring_closing_set|`` must equal
           ``|E_full| - |V| + |connected_components|``.  Too few means some
           cycle has no ring-closing bond; too many means some non-cyclic bond
           was incorrectly removed.
        2. **Each ring-closing bond closes exactly one cycle** -- for each
           ring-closing bond ``(u, v)``, ``u`` and ``v`` must be connected in
           the spanning tree (otherwise the bond was a bridge, not part of any
           cycle, and its removal disconnects the graph).
        3. **No shared ring-closing bonds** -- no ring-closing bond may appear
           as an edge inside another ring-closing bond's fundamental cycle.
           A fundamental cycle is: the ring-closing bond ``(u, v)`` plus the
           unique path from ``u`` to ``v`` through the spanning tree.

           This check is theoretically guaranteed by spanning-tree construction
           (ring-closing bonds are chords; chords are never spanning-tree edges;
           fundamental cycles consist only of spanning-tree edges + the chord).
           The explicit check here serves as a defence against bugs in the
           classifier or bond-type annotation.

        Raises
        ------
        ValueError
            On any of the three failures above, with a diagnostic message.
        """
        # ---- Check 1: cyclomatic number ----
        n_comp: int = nx.number_connected_components(spanning_tree)
        cyclomatic: int = (
            full_graph.number_of_edges() - full_graph.number_of_nodes() + n_comp
        )
        if len(ring_closing_set) != cyclomatic:
            raise ValueError(
                "Ring-closing bond count (%d) != cyclomatic number (%d). "
                "Either a cycle has no ring-closing bond, or a non-cyclic "
                "(bridge) bond was incorrectly removed."
                % (len(ring_closing_set), cyclomatic)
            )

        # ---- Checks 2 & 3: per ring-closing bond ----
        # Compute fundamental cycles: for each chord (u,v), the fundamental
        # cycle is the set of spanning-tree edges on the path u -> v, plus
        # the chord itself.
        fundamental_cycle_edges: dict[frozenset[int], frozenset[frozenset[int]]] = {}

        for rc_bond in ring_closing_set:
            u, v = tuple(rc_bond)

            # Check 2: endpoints must be connected in the spanning tree
            if not nx.has_path(spanning_tree, u, v):
                raise ValueError(
                    "Ring-closing bond %d -- %d has no path in the spanning "
                    "tree.  The bond may be a bridge (not part of any cycle), "
                    "or the spanning tree is disconnected at this node." % (u, v)
                )

            path: list[int] = nx.shortest_path(spanning_tree, u, v)
            tree_edges: frozenset[frozenset[int]] = frozenset(
                frozenset((path[i], path[i + 1])) for i in range(len(path) - 1)
            )
            fundamental_cycle_edges[rc_bond] = tree_edges | {rc_bond}

        # Check 3: no ring-closing bond appears inside another's fundamental cycle
        # (the fundamental cycle uses only spanning-tree edges + the chord itself;
        # ring-closing bonds are never spanning-tree edges, so this should always
        # pass -- any failure signals a bug upstream)
        for rc_bond, cycle_edges in fundamental_cycle_edges.items():
            for other_rc in ring_closing_set:
                if other_rc == rc_bond:
                    continue
                if other_rc in cycle_edges:
                    u1, v1 = tuple(rc_bond)
                    u2, v2 = tuple(other_rc)
                    raise ValueError(
                        "Ring-closing bond %d -- %d appears inside the "
                        "fundamental cycle of ring-closing bond %d -- %d.  "
                        "Ring-closing bonds must be independent (one per "
                        "fundamental cycle, no sharing).  This indicates a "
                        "bug in the classifier or bond-removal logic."
                        % (u2, v2, u1, v1)
                    )

        logger.debug(
            "Ring-closing bond validation passed: %d bond(s), cyclomatic number %d.",
            len(ring_closing_set),
            cyclomatic,
        )

    def _atom_label(self, idx: int) -> str:
        """Return a human-readable label for atom *idx* (used in log messages)."""
        a = self.molecule.atoms[idx]
        return "%s%d_%s_%d" % (a.residue.name, a.residue.idx + 1, a.name, idx + 1)

    def _check_residual_cycles(self, g: nx.Graph) -> None:
        """
        Raise ``ValueError`` if *g* contains cycles, logging each cycle.
        Called as a final sanity gate after all ring-closing removal steps.
        """
        if nx.is_forest(g):
            return

        cycles = nx.cycle_basis(g)
        logger.error("Graph contains %d residual cycle(s):", len(cycles))
        for ci, cycle in enumerate(cycles, start=1):
            logger.error("  Cycle %d:", ci)
            for j in range(len(cycle)):
                n1 = self._atom_label(cycle[j])
                n2 = self._atom_label(cycle[(j + 1) % len(cycle)])
                logger.error("    %s -- %s", n1, n2)
        raise ValueError(
            "Acyclic graph still contains cycles after ring-closure removal."
        )

    def _check_disconnected_graph(self, g: nx.Graph) -> None:
        """
        Raise ``ValueError`` if *g* is disconnected, logging each component.
        """
        if nx.is_connected(g):
            return

        components = list(nx.connected_components(g))
        logger.error(
            "Graph is disconnected (%d components) after ring-closure removal.",
            len(components),
        )
        for ci, comp in enumerate(components, start=1):
            logger.error("  Component %d (%d atoms):", ci, len(comp))
            for idx in comp:
                logger.error("    - %s", self._atom_label(idx))
        raise ValueError("Acyclic graph is disconnected after ring-closure removal.")

    # ------------------------------------------------------------------
    # Z-matrix
    # ------------------------------------------------------------------

    def _build_z_matrix(self, root: pmd.Atom) -> None:
        """
        Build a Z-matrix (BAT-style internal-coordinate tree) rooted at *root*.

        All four output lists have length ``num_atoms``.  Rows that have fewer
        than four defined references store the sentinel ``-1``:

        +---------+---+---------+---------+---------+
        | Row     | i | j       | k       | l       |
        +=========+===+=========+=========+=========+
        | 0 (root)| * | -1      | -1      | -1      |
        +---------+---+---------+---------+---------+
        | 1       | * | valid   | -1      | -1      |
        +---------+---+---------+---------+---------+
        | 2       | * | valid   | valid   | -1      |
        +---------+---+---------+---------+---------+
        | 3+      | * | valid   | valid   | valid   |
        +---------+---+---------+---------+---------+

        Reference-atom constraints (rows 3+):
          j -- selected tree-neighbor that "discovers" atom i
          k -- selected, non-terminal (full-graph degree > 1) tree-neighbor of j
          l -- any selected tree-neighbor of k, excluding j

        Ring-closing bonds are excluded via ``self.acyclic_graph``.

        Root-triplet selection heuristic
        ---------------------------------
        initial  = root (caller-specified)
        second   = heaviest non-terminal tree-neighbor of root; falls back to
                   heaviest if all neighbors are terminal
        third    = heaviest non-terminal tree-neighbor of second (excl. root);
                   falls back to heaviest; for n == 3 the non-terminal constraint
                   is dropped entirely

        Algorithm
        ---------
        A work-queue / deferred-retry scheme replaces the original O(n^2) scan.
        Complexity: O(n * d) where d is the maximum tree degree.

        Corner cases
        ------------
        n == 1   -- all four lists get one entry: [root.idx], [-1], [-1], [-1]
        n == 2   -- each list has two entries; j[0] = k[0] = k[1] = l[0] = l[1] = -1
        n == 3   -- l is all -1; k[0] = k[1] = -1
        Isolated root    -- ValueError
        Unroutable triplet -- ValueError
        Disconnected graph -- ValueError
        Stalled traversal  -- ValueError (degenerate topology or wrong root)

        Parameters
        ----------
        root : pmd.Atom
            Root atom.  Must be a node in ``self.acyclic_graph``.

        Raises
        ------
        ValueError
            See corner cases above.
        """
        n: int = self.num_atoms

        if n == 0:
            raise ValueError("Molecule contains no atoms; cannot build Z-matrix.")
        if root.idx not in self.acyclic_graph:
            raise ValueError(
                "Root atom (idx=%d, name=%r) is not a node in the acyclic graph."
                % (root.idx, getattr(root, "name", "?"))
            )

        self.z_matrix_i = []
        self.z_matrix_j = []
        self.z_matrix_k = []
        self.z_matrix_l = []

        # Pre-compute sorted adjacency lists once; used throughout.
        # Using acyclic_graph ensures ring-closing bonds are excluded.
        tree_adj: dict[int, list[pmd.Atom]] = {
            a.idx: self._sort_atoms_by_mass(
                [self._atom_by_idx[nb] for nb in self.acyclic_graph.neighbors(a.idx)]
            )
            for a in self.molecule.atoms
        }

        # -- n == 1 ----------------------------------------------------------
        if n == 1:
            self.z_matrix_i.append(root.idx)
            self.z_matrix_j.append(_ZM_SENTINEL)
            self.z_matrix_k.append(_ZM_SENTINEL)
            self.z_matrix_l.append(_ZM_SENTINEL)
            return

        root_nbrs: list[pmd.Atom] = tree_adj[root.idx]
        if not root_nbrs:
            raise ValueError(
                "Root atom (idx=%d) has no tree-neighbors in a %d-atom molecule. "
                "The acyclic graph may be disconnected at this node." % (root.idx, n)
            )

        # -- n == 2 ----------------------------------------------------------
        if n == 2:
            second: pmd.Atom = root_nbrs[0]
            self.z_matrix_i.extend([root.idx, second.idx])
            self.z_matrix_j.extend([_ZM_SENTINEL, root.idx])
            self.z_matrix_k.extend([_ZM_SENTINEL, _ZM_SENTINEL])
            self.z_matrix_l.extend([_ZM_SENTINEL, _ZM_SENTINEL])
            return

        # -- Root triplet for n >= 3 -----------------------------------------
        initial_atom: pmd.Atom = root

        non_term_root: list[pmd.Atom] = [
            a for a in root_nbrs if len(a.bond_partners) > 1
        ]
        second_atom: pmd.Atom = (self._sort_atoms_by_mass(non_term_root) or root_nbrs)[
            0
        ]

        second_excl: list[pmd.Atom] = [
            a for a in tree_adj[second_atom.idx] if a.idx != initial_atom.idx
        ]
        if not second_excl:
            raise ValueError(
                "Cannot build Z-matrix root triplet: second atom (idx=%d) has "
                "no tree-neighbors other than root (idx=%d) in a %d-atom molecule. "
                "Choose a different root atom." % (second_atom.idx, initial_atom.idx, n)
            )

        if n > 3:
            non_term_excl: list[pmd.Atom] = [
                a for a in second_excl if len(a.bond_partners) > 1
            ]
            third_atom: pmd.Atom = (
                self._sort_atoms_by_mass(non_term_excl)
                or self._sort_atoms_by_mass(second_excl)
            )[0]
        else:
            # n == 3: drop the non-terminal constraint for the third atom
            third_atom = self._sort_atoms_by_mass(second_excl)[0]

        # Write root triplet rows 0, 1, 2
        self.z_matrix_i.append(initial_atom.idx)
        self.z_matrix_j.append(_ZM_SENTINEL)
        self.z_matrix_k.append(_ZM_SENTINEL)
        self.z_matrix_l.append(_ZM_SENTINEL)

        self.z_matrix_i.append(second_atom.idx)
        self.z_matrix_j.append(initial_atom.idx)
        self.z_matrix_k.append(_ZM_SENTINEL)
        self.z_matrix_l.append(_ZM_SENTINEL)

        self.z_matrix_i.append(third_atom.idx)
        self.z_matrix_j.append(second_atom.idx)
        self.z_matrix_k.append(initial_atom.idx)
        self.z_matrix_l.append(_ZM_SENTINEL)

        if n == 3:
            return

        # -- General case: n >= 4 --------------------------------------------
        #
        # Work-queue:  each entry is (a0, a1) where a0 is an unselected atom
        # discovered via the already-selected atom a1.
        #
        # Deferred:    pairs whose (a2, a3) search failed in this batch; retried
        # after each successful placement batch when new selected atoms may
        # provide the missing a2/a3 references.
        #
        # Stall detection:  if the deferred retry produces zero placements, the
        # construction cannot proceed and a ValueError is raised.

        selected: set[int] = {initial_atom.idx, second_atom.idx, third_atom.idx}

        # Seed in insertion order (initial, second, third) -- same ordering as
        # the original algorithm's `selected_atoms` list for determinism.
        work: deque[tuple[pmd.Atom, pmd.Atom]] = deque()
        for seed in (initial_atom, second_atom, third_atom):
            for nb in tree_adj[seed.idx]:
                if nb.idx not in selected:
                    work.append((nb, seed))

        deferred: list[tuple[pmd.Atom, pmd.Atom]] = []
        # Number of atoms placed the last time we flushed deferred -> work.
        # If this does not grow before the next flush, the traversal is stuck.
        placed_at_last_retry: int = 3

        while len(self.z_matrix_i) < n:
            if not work:
                if not deferred:
                    missing: list[int] = [
                        idx for idx in self.acyclic_graph.nodes if idx not in selected
                    ]
                    raise ValueError(
                        "Z-matrix construction terminated at %d/%d atoms with "
                        "empty queues. Disconnected acyclic graph? "
                        "Unplaced atom indices: %s."
                        % (len(self.z_matrix_i), n, missing)
                    )
                if len(self.z_matrix_i) == placed_at_last_retry:
                    raise ValueError(
                        "Z-matrix stalled at %d/%d atoms. %d deferred pair(s) "
                        "could not be resolved: %s. "
                        "Likely cause: degenerate topology (e.g. all neighbors "
                        "of the central atom are terminal in the full molecular "
                        "graph) or an unsuitable root atom."
                        % (
                            len(self.z_matrix_i),
                            n,
                            len(deferred),
                            [(a0.idx, a1.idx) for a0, a1 in deferred],
                        )
                    )
                placed_at_last_retry = len(self.z_matrix_i)
                work.extend(deferred)
                deferred.clear()

            next_work: deque[tuple[pmd.Atom, pmd.Atom]] = deque()

            while work:
                a0, a1 = work.popleft()

                # Guard: duplicate queue entry (cannot happen in a true tree,
                # but prevents silent double-placement if seeds overlap)
                if a0.idx in selected:
                    continue

                placed: bool = False
                for a2 in tree_adj[a1.idx]:
                    if a2.idx == a0.idx or a2.idx not in selected:
                        continue
                    if len(a2.bond_partners) <= 1:
                        # a2 must be non-terminal in the *full* molecular graph
                        continue
                    for a3 in tree_adj[a2.idx]:
                        if a3.idx == a1.idx or a3.idx not in selected:
                            continue
                        # Valid (a2, a3) pair found -- emit row
                        self.z_matrix_i.append(a0.idx)
                        self.z_matrix_j.append(a1.idx)
                        self.z_matrix_k.append(a2.idx)
                        self.z_matrix_l.append(a3.idx)
                        selected.add(a0.idx)
                        placed = True
                        for nb in tree_adj[a0.idx]:
                            if nb.idx not in selected:
                                next_work.append((nb, a0))
                        break
                    if placed:
                        break

                if not placed:
                    deferred.append((a0, a1))

            work = next_work
            # If work is now empty and deferred is non-empty, the top of the
            # loop will either retry (if progress was made) or raise.
