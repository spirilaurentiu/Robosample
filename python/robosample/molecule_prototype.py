"""molecule_prototype.py

Pre-processed, simulation-ready snapshot of a single ParmEd molecule.

Index spaces
------------
Two atom-index spaces are used throughout (see ``indexing.CompoundIndex``):

* **local**    -- ParmEd / AMBER prmtop atom index.
* **compound** -- BFS index from the spanning-forest root; index 0 is the root,
  and a smaller index means closer to the root.

Every topology array is exposed in **compound** indices (the canonical,
BFS-ordered representation) with a parallel ``*_local`` array carrying the
ParmEd indices for the same elements.  Atom-level arrays are stored in
*compound order* (position == compound index).

Orientation
-----------
Path-like elements (bonds, angles, proper torsions, Urey-Bradley, 1-4 pairs,
exclusions) are oriented so the first atom is closest to the root and the last
is farthest.  The whole tuple is reversed when needed, so the central atom of
an angle stays central and proper-dihedral energies are preserved (a proper
dihedral is invariant under full atom-order reversal).

* **Improper torsions are NOT reoriented** -- reversing a 4-atom improper flips
  the dihedral sign, changing ``k*(phi - phi0)^2`` when ``phi0 != 0``.  They keep
  their ParmEd atom order and are only index-remapped.
* **The Z-matrix is NOT reoriented** -- its rows define each atom relative to
  earlier-placed atoms (back toward the root) by construction.

Unit conventions: lengths nm, energies kJ/mol, angles rad.
"""

from __future__ import annotations

import logging
from functools import cached_property
from typing import Any

import mdtraj as md
import networkx as nx
import numpy as np
import parmed as pmd

from . import acyclic_graph as ag
from . import prmtop_reader
from . import z_matrix as zm
from .amber_dihedral_classifier import AmberDihedralClassifier
from .amber_dihedral_types import DihedralType
from .indexing import CompoundIndex
from .secondary_structure import BondDSSP, DSSPCode
from .units import ANG_TO_NM, DEG_TO_RAD, KCAL_TO_KJ

logger: logging.Logger = logging.getLogger(__name__)

# Sentinel for absent Z-matrix references.
_ZM_SENTINEL: int = zm.ZM_SENTINEL


class MoleculePrototype:
    """Immutable pre-processed snapshot of a single ParmEd molecule."""

    def __init__(
        self,
        molecule: pmd.Structure,
        dihedral_classifier: AmberDihedralClassifier,
    ) -> None:
        self.molecule = molecule
        self.dihedral_classifier = dihedral_classifier

        self.num_atoms: int = len(molecule.atoms)
        self.num_residues: int = max(
            (a.residue.idx + 1 for a in molecule.atoms), default=0
        )

        self._atom_by_idx: dict[int, pmd.Atom] = {a.idx: a for a in molecule.atoms}

        # -- Spanning forest + dihedral-type cache (local indices) ---------
        self.acyclic_graph, self._dihedral_type_cache = ag.build_acyclic_graph(
            molecule, dihedral_classifier
        )

        # -- Root: heaviest terminal atom (or the sole atom for n == 1) ----
        if self.num_atoms == 1:
            root_atom: pmd.Atom = molecule.atoms[0]
        else:
            terminals = [a for a in molecule.atoms if len(a.bond_partners) == 1]
            root_atom = self._sort_atoms_by_mass(terminals)[0]
        self.atoms_root_index: int = root_atom.idx  # local index of root
        self.atoms_root_compound_index: int = 0  # root is always compound 0

        # -- BFS compound indexing -----------------------------------------
        self.compound_index = CompoundIndex(self.acyclic_graph, self.atoms_root_index)
        self.compound_to_local: list[int] = self.compound_index.compound_to_local
        self.local_to_compound: dict[int, int] = self.compound_index.local_to_compound

        logger.debug(
            "Root atom: local=%d compound=0 name=%s mass=%.3f",
            root_atom.idx,
            getattr(root_atom, "name", "?"),
            getattr(root_atom, "mass", 0.0),
        )

        self._build_atom_arrays()
        dssp = self._compute_dssp()  # per-residue codes (order-independent)
        self._build_bond_arrays(dssp)
        self._build_angle_arrays()
        self._build_periodic_torsion_arrays()
        self._build_harmonic_torsion_arrays()
        self._build_urey_bradley_arrays()
        self._build_nonbonded_arrays()
        self._build_z_matrix_arrays()

    # ------------------------------------------------------------------
    # Atoms (compound order)
    # ------------------------------------------------------------------

    def _build_atom_arrays(self) -> None:
        order = self.compound_to_local
        atoms = [self._atom_by_idx[local] for local in order]

        # Local index of each compound atom (compound c -> ParmEd idx).
        self.atoms_local_index: list[int] = list(order)

        # SimTK Compound::AtomIndex: the atom's index within its compound.
        # Atoms are stored in compound order, so this is simply 0..n-1 and is
        # local to the compound (no system atom offset is applied downstream).
        self.atoms_compound_atom_index: list[int] = list(range(self.num_atoms))

        self.atoms_mass = [a.mass for a in atoms]
        self.atoms_charge = [a.charge for a in atoms]
        self.atoms_sigma = [a.sigma * ANG_TO_NM for a in atoms]
        self.atoms_epsilon = [a.epsilon * KCAL_TO_KJ for a in atoms]
        self.atoms_radius = [a.solvent_radius * ANG_TO_NM for a in atoms]
        self.atoms_screen = [a.screen for a in atoms]
        self.atoms_x = [a.xx * ANG_TO_NM for a in atoms]
        self.atoms_y = [a.xy * ANG_TO_NM for a in atoms]
        self.atoms_z = [a.xz * ANG_TO_NM for a in atoms]

        # ---- Identity / type-table fields (compound order) -----------------
        # Stored in compound (BFS) order like every other atom array.  None of
        # these is a *system atom index*, so downstream flattening applies no
        # atom offset to them.
        #
        # Element name (e.g. "C", "N").  element_symbol mirrors element_name:
        # the consumer only needs *a* stable per-atom string here, not a
        # canonical symbol, so the exact value is not significant.
        self.atoms_element_name: list[str] = [a.element_name for a in atoms]
        self.atoms_element_symbol: list[str] = [a.element_name for a in atoms]
        # Lennard-Jones (nonbonded) parameter-table index, 0-based.  ParmEd's
        # nb_idx is 1-based, hence the -1.  This indexes the LJ *type* table, not
        # the atom list, so it is never shifted by the atom offset.
        self.atoms_nonbonded_index: list[int] = [a.nb_idx - 1 for a in atoms]

        # NOTE: atoms_unique_name is deliberately NOT built here.  It embeds
        # GLOBAL, whole-system 1-based residue and prmtop atom numbers
        # (e.g. "ALA1_N_4", where the trailing 4 is the atom's index in the
        # .prmtop file).  A prototype is shared across every instance of its
        # molecule type and has no knowledge of its global position, so a
        # prototype-local index would be wrong.  The name is assembled in
        # Context.load_amber, where the per-instance atom/residue offsets are
        # known.

    # ------------------------------------------------------------------
    # DSSP (per residue; independent of atom-index space)
    # ------------------------------------------------------------------

    def _compute_dssp(self):
        # The trajectory must use *local*-order coordinates so it matches the
        # local-order OpenMM topology.  DSSP is returned per residue, so the
        # result is independent of atom (compound) reordering.
        xyz_local = (
            np.array(
                [
                    [a.xx for a in self.molecule.atoms],
                    [a.xy for a in self.molecule.atoms],
                    [a.xz for a in self.molecule.atoms],
                ]
            ).T
            * ANG_TO_NM
        )
        traj = md.Trajectory(
            xyz=xyz_local,
            topology=md.Topology.from_openmm(self.molecule.topology),
        )
        return md.compute_dssp(traj, simplified=False)[0]

    # ------------------------------------------------------------------
    # Bonds (BFS order; oriented parent -> child, i.e. closest -> farthest)
    # ------------------------------------------------------------------

    def _build_bond_arrays(self, dssp) -> None:
        bond_lookup = {
            frozenset((b.atom1.idx, b.atom2.idx)): b for b in self.molecule.bonds
        }

        # Tree bonds in BFS order: bfs_edges yields (parent, child) with the
        # parent already closer to the root (smaller compound index).
        tree_edges = list(
            nx.bfs_edges(self.acyclic_graph, source=self.atoms_root_index)
        )
        tree_edge_set = {frozenset(e) for e in tree_edges}

        # ordered_local: (i_local, j_local, is_ring_closing), i closest to root.
        ordered_local: list[tuple[int, int, bool]] = []
        for parent, child in tree_edges:
            i_l, j_l = self.compound_index.orient_local((parent, child))
            ordered_local.append((i_l, j_l, False))
        for b in self.molecule.bonds:
            key = frozenset((b.atom1.idx, b.atom2.idx))
            if key in tree_edge_set:
                continue
            i_l, j_l = self.compound_index.orient_local((b.atom1.idx, b.atom2.idx))
            ordered_local.append((i_l, j_l, True))

        if len(ordered_local) != len(self.molecule.bonds):
            raise ValueError(
                "Bond bookkeeping mismatch: ordered %d bonds but molecule has %d."
                % (len(ordered_local), len(self.molecule.bonds))
            )

        self.num_bonds = len(ordered_local)
        self.bonds_i_local = [t[0] for t in ordered_local]
        self.bonds_j_local = [t[1] for t in ordered_local]
        self.bonds_i = [self.compound_index.to_compound(t[0]) for t in ordered_local]
        self.bonds_j = [self.compound_index.to_compound(t[1]) for t in ordered_local]
        self.bonds_is_ring_closing = [t[2] for t in ordered_local]

        self.bonds_stiffness = []
        self.bonds_equilibrium = []
        self.bonds_dihedral_type = []
        self.bonds_secondary_structure = []
        for i_l, j_l, _rc in ordered_local:
            b = bond_lookup[frozenset((i_l, j_l))]
            self.bonds_stiffness.append(b.type.k * KCAL_TO_KJ / ANG_TO_NM**2)
            self.bonds_equilibrium.append(b.type.req * ANG_TO_NM)
            self.bonds_dihedral_type.append(
                self._dihedral_type_cache.get(
                    frozenset((i_l, j_l)), DihedralType.UNKNOWN
                )
            )
            self.bonds_secondary_structure.append(
                BondDSSP(
                    atom1_code=DSSPCode.from_mdtraj(
                        dssp[self._atom_by_idx[i_l].residue.idx]
                    ),
                    atom2_code=DSSPCode.from_mdtraj(
                        dssp[self._atom_by_idx[j_l].residue.idx]
                    ),
                ).resolve(strategy="priority")
            )

    # ------------------------------------------------------------------
    # Angles (oriented i -> k, closest -> farthest; central atom j fixed)
    # ------------------------------------------------------------------

    def _build_angle_arrays(self) -> None:
        self.num_angles = len(self.molecule.angles)
        self.angles_i_local, self.angles_j_local, self.angles_k_local = [], [], []
        self.angles_i, self.angles_j, self.angles_k = [], [], []
        self.angles_equilibrium, self.angles_stiffness = [], []

        for ang in self.molecule.angles:
            i_l, j_l, k_l = self.compound_index.orient_local(
                (ang.atom1.idx, ang.atom2.idx, ang.atom3.idx)
            )
            self.angles_i_local.append(i_l)
            self.angles_j_local.append(j_l)
            self.angles_k_local.append(k_l)
            self.angles_i.append(self.compound_index.to_compound(i_l))
            self.angles_j.append(self.compound_index.to_compound(j_l))
            self.angles_k.append(self.compound_index.to_compound(k_l))
            self.angles_equilibrium.append(ang.type.theteq * DEG_TO_RAD)
            self.angles_stiffness.append(ang.type.k * KCAL_TO_KJ)

    # ------------------------------------------------------------------
    # Periodic torsions (one record per term).
    # Proper torsions oriented i -> l; impropers keep ParmEd order.
    # ------------------------------------------------------------------

    def _build_periodic_torsion_arrays(self) -> None:
        records: list[tuple[pmd.Dihedral, pmd.DihedralType]] = []
        for d in self.molecule.dihedrals:
            if isinstance(d.type, pmd.DihedralTypeList):
                records.extend((d, dt) for dt in d.type)
            elif isinstance(d.type, pmd.DihedralType):
                records.append((d, d.type))

        self.num_periodic_torsions = len(records)
        self.periodic_torsions_improper = []
        self.periodic_torsions_i_local, self.periodic_torsions_j_local = [], []
        self.periodic_torsions_k_local, self.periodic_torsions_l_local = [], []
        self.periodic_torsions_i, self.periodic_torsions_j = [], []
        self.periodic_torsions_k, self.periodic_torsions_l = [], []
        self.periodic_torsions_n = []
        self.periodic_torsions_phase = []
        self.periodic_torsions_stiffness = []

        for d, dt in records:
            quad = (d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx)
            # Reorient only proper torsions; impropers keep ParmEd order.
            oriented = quad if d.improper else self.compound_index.orient_local(quad)
            i_l, j_l, k_l, l_l = oriented

            self.periodic_torsions_improper.append(d.improper)
            self.periodic_torsions_i_local.append(i_l)
            self.periodic_torsions_j_local.append(j_l)
            self.periodic_torsions_k_local.append(k_l)
            self.periodic_torsions_l_local.append(l_l)
            self.periodic_torsions_i.append(self.compound_index.to_compound(i_l))
            self.periodic_torsions_j.append(self.compound_index.to_compound(j_l))
            self.periodic_torsions_k.append(self.compound_index.to_compound(k_l))
            self.periodic_torsions_l.append(self.compound_index.to_compound(l_l))
            self.periodic_torsions_n.append(dt.per)
            self.periodic_torsions_phase.append(dt.phase * DEG_TO_RAD)
            self.periodic_torsions_stiffness.append(dt.phi_k * KCAL_TO_KJ)

    # ------------------------------------------------------------------
    # Harmonic (CHARMM) improper torsions -- impropers, ParmEd order kept.
    # ------------------------------------------------------------------

    def _build_harmonic_torsion_arrays(self) -> None:
        self.num_harmonic_torsions = len(self.molecule.impropers)
        self.harmonic_torsions_i_local, self.harmonic_torsions_j_local = [], []
        self.harmonic_torsions_k_local, self.harmonic_torsions_l_local = [], []
        self.harmonic_torsions_i, self.harmonic_torsions_j = [], []
        self.harmonic_torsions_k, self.harmonic_torsions_l = [], []
        self.harmonic_torsions_stiffness, self.harmonic_torsions_phase = [], []

        for imp in self.molecule.impropers:
            quad = (imp.atom1.idx, imp.atom2.idx, imp.atom3.idx, imp.atom4.idx)
            i_l, j_l, k_l, l_l = quad  # impropers: order is meaningful, keep it
            self.harmonic_torsions_i_local.append(i_l)
            self.harmonic_torsions_j_local.append(j_l)
            self.harmonic_torsions_k_local.append(k_l)
            self.harmonic_torsions_l_local.append(l_l)
            self.harmonic_torsions_i.append(self.compound_index.to_compound(i_l))
            self.harmonic_torsions_j.append(self.compound_index.to_compound(j_l))
            self.harmonic_torsions_k.append(self.compound_index.to_compound(k_l))
            self.harmonic_torsions_l.append(self.compound_index.to_compound(l_l))
            self.harmonic_torsions_stiffness.append(imp.type.psi_k * KCAL_TO_KJ)
            self.harmonic_torsions_phase.append(imp.type.psi_eq * DEG_TO_RAD)

    # ------------------------------------------------------------------
    # Urey-Bradley (1,3 pair; oriented i -> k closest -> farthest)
    # ------------------------------------------------------------------

    def _build_urey_bradley_arrays(self) -> None:
        self.num_urey_bradley = len(self.molecule.urey_bradleys)
        self.urey_bradley_i_local, self.urey_bradley_k_local = [], []
        self.urey_bradley_i, self.urey_bradley_k = [], []
        self.urey_bradley_stiffness, self.urey_bradley_equilibrium = [], []

        for ub in self.molecule.urey_bradleys:
            i_l, k_l = self.compound_index.orient_local((ub.atom1.idx, ub.atom2.idx))
            self.urey_bradley_i_local.append(i_l)
            self.urey_bradley_k_local.append(k_l)
            self.urey_bradley_i.append(self.compound_index.to_compound(i_l))
            self.urey_bradley_k.append(self.compound_index.to_compound(k_l))
            self.urey_bradley_stiffness.append(ub.type.k * KCAL_TO_KJ / ANG_TO_NM**2)
            self.urey_bradley_equilibrium.append(ub.type.req * ANG_TO_NM)

    # ------------------------------------------------------------------
    # 1-4 scaling pairs and exclusions (parsed from raw prmtop, then oriented)
    # ------------------------------------------------------------------

    def _build_nonbonded_arrays(self) -> None:
        self.num_scaling14 = 0
        self.scaling14_i_local, self.scaling14_l_local = [], []
        self.scaling14_i, self.scaling14_l = [], []
        self.scaling14_charge_product, self.scaling14_epsilon = [], []
        self.scaling14_sigma = []

        self.num_exclusions = 0
        self.exclusion_i_local, self.exclusion_j_local = [], []
        self.exclusion_i, self.exclusion_j = [], []

        if not hasattr(self.molecule, "parm_data"):
            logger.warning(
                "molecule has no parm_data attribute; scaling14 and exclusions "
                "will be empty.  Load the molecule from an AMBER .prmtop file "
                "to populate these arrays."
            )
            return

        tables = prmtop_reader.load_nonbonded_exceptions(
            self.molecule.parm_data, self.molecule
        )

        for rec in tables.scaling14:
            # Pair quantities are symmetric; orientation only fixes which is i/l.
            i_l, l_l = self.compound_index.orient_local((rec.i_local, rec.l_local))
            self.scaling14_i_local.append(i_l)
            self.scaling14_l_local.append(l_l)
            self.scaling14_i.append(self.compound_index.to_compound(i_l))
            self.scaling14_l.append(self.compound_index.to_compound(l_l))
            self.scaling14_charge_product.append(rec.charge_product)
            self.scaling14_epsilon.append(rec.epsilon)
            self.scaling14_sigma.append(rec.sigma)
        self.num_scaling14 = len(self.scaling14_i)

        for rec in tables.exclusions:
            i_l, j_l = self.compound_index.orient_local((rec.i_local, rec.j_local))
            self.exclusion_i_local.append(i_l)
            self.exclusion_j_local.append(j_l)
            self.exclusion_i.append(self.compound_index.to_compound(i_l))
            self.exclusion_j.append(self.compound_index.to_compound(j_l))
        self.num_exclusions = len(self.exclusion_i)

    # ------------------------------------------------------------------
    # Z-matrix (own traversal order; references remapped to compound)
    # ------------------------------------------------------------------

    def _build_z_matrix_arrays(self) -> None:
        z_i_l, z_j_l, z_k_l, z_l_l = zm.build_z_matrix(
            self.molecule, self.acyclic_graph, self.atoms_root_index
        )
        self.z_matrix_i_local = z_i_l
        self.z_matrix_j_local = z_j_l
        self.z_matrix_k_local = z_k_l
        self.z_matrix_l_local = z_l_l

        def remap(values: list[int]) -> list[int]:
            return [
                self.compound_index.to_compound(v)
                if v != _ZM_SENTINEL
                else _ZM_SENTINEL
                for v in values
            ]

        self.z_matrix_i = remap(z_i_l)
        self.z_matrix_j = remap(z_j_l)
        self.z_matrix_k = remap(z_k_l)
        self.z_matrix_l = remap(z_l_l)

    # ------------------------------------------------------------------
    # Derived / lazy
    # ------------------------------------------------------------------

    @cached_property
    def dihedral_bond_records(self) -> list[dict[str, Any]]:
        """
        One record per unique proper-torsion central bond, in compound indices.

        Ring-closing central bonds and improper torsions are excluded.  Each
        record carries the central-bond atoms oriented closest -> farthest
        (``j`` closer to root than ``k``), in both compound and local indices.
        """
        ring_closing_bonds: set[frozenset[int]] = {
            frozenset((self.bonds_i_local[n], self.bonds_j_local[n]))
            for n in range(self.num_bonds)
            if self.bonds_is_ring_closing[n]
        }

        seen: dict[frozenset[int], dict[str, Any]] = {}
        for d in self.molecule.dihedrals:
            if d.improper:
                continue
            key = frozenset((d.atom2.idx, d.atom3.idx))
            if key in ring_closing_bonds:
                continue
            dtype = self._dihedral_type_cache.get(key, DihedralType.UNKNOWN)
            if key not in seen or (
                seen[key]["dihedral_type"] is DihedralType.UNKNOWN
                and dtype is not DihedralType.UNKNOWN
            ):
                j_l, k_l = self.compound_index.orient_local((d.atom2.idx, d.atom3.idx))
                seen[key] = {
                    "j": self.compound_index.to_compound(j_l),
                    "k": self.compound_index.to_compound(k_l),
                    "j_local": j_l,
                    "k_local": k_l,
                    "dihedral_type": dtype,
                    "resid": self._atom_by_idx[j_l].residue.idx,
                    "resname": self._atom_by_idx[j_l].residue.name,
                }
        return list(seen.values())

    @staticmethod
    def _sort_atoms_by_mass(atoms: list[pmd.Atom]) -> list[pmd.Atom]:
        """Sort by descending mass; ascending atom index breaks ties."""
        return sorted(atoms, key=lambda a: (-a.mass, a.idx))
