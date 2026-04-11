from dataclasses import dataclass, field
from typing import Dict, FrozenSet, Iterable, List, Optional, Set, Tuple

import networkx as nx
import parmed as pmd

from . import robo_bindings as rb
from .amber_dihedral import DihedralClassifier
from .bond_util import is_rigid_bond
from .graph_utils import GraphTraversalUtils


@dataclass(frozen=True)
class BondProperties:
    parent_atom_prmtop_index: int
    child_atom_prmtop_index: int
    is_rigid: bool
    sg_sg_bond_distance: int | None
    sum_of_degrees: int


@dataclass(slots=True)
class AtomParams:
    local_index: int
    compound_atom_index: int
    element_name: str
    element_symbol: str
    atomic_number: int
    charge_e: float
    mass_daltons: float
    vdw_radius_nm: float
    vdw_well_depth_kj: float
    sigma_nm: float
    solvent_radius_nm: float
    screen: float
    neighbors_local_indices: List[int]
    root: bool


@dataclass(slots=True)
class BondParams:
    parent_local_index: int
    child_local_index: int
    parent_compound_atom_index: int
    child_compound_atom_index: int
    stiffness_in_kj_per_nm_sq: float
    nominal_length_in_nm: float
    is_ring_closing: bool
    dihedral_type: str
    gparent_local_index: Optional[int] = field(default=None)
    nephew_local_index: Optional[int] = field(default=None)


@dataclass(slots=True)
class AngleParams:
    local_indices: Tuple[int, int, int]
    compound_atom_indices: Tuple[int, int, int]
    stiffness_in_kj_per_rad_sq: float
    nominal_angle_in_deg: float


@dataclass(slots=True)
class PeriodicTorsionTerm:
    amplitude_in_kj: float
    phase_in_deg: float
    periodicity: int


@dataclass(slots=True)
class PeriodicTorsionParams:
    local_indices: Tuple[int, int, int, int]
    compound_atom_indices: Tuple[int, int, int, int]
    is_improper: bool
    terms: List[PeriodicTorsionTerm]


@dataclass(slots=True)
class HarmonicImproperTorsionParams:
    local_indices: Tuple[int, int, int, int]
    compound_atom_indices: Tuple[int, int, int, int]
    stiffness_in_kj_per_rad_sq: float
    nominal_angle_in_rad: float


@dataclass(slots=True)
class ZMatrixRowAtomGroups:
    i_global_atom_index: int | None = field(default=None)
    j_global_atom_index: int | None = field(default=None)
    k_global_atom_index: int | None = field(default=None)
    l_global_atom_index: int | None = field(default=None)

    i_compound_atom_index: int | None = field(default=None)
    j_compound_atom_index: int | None = field(default=None)
    k_compound_atom_index: int | None = field(default=None)
    l_compound_atom_index: int | None = field(default=None)


class MoleculePrototype:
    def __init__(self, molecule: pmd.Structure):
        self.molecule = molecule
        self.dihedral_classifier = DihedralClassifier()

        # Find bonds from macrocycles
        # Their order and dihedral labeling is irrelevant
        self.macrocycle_bonds = self._build_bonds_macrocycles()

        # Root is the heaviest terminal atom with the smallest index
        if len(self.molecule.atoms) == 1:
            terminal_atoms = self.molecule.atoms
        else:
            terminal_atoms = [
                a for a in self.molecule.atoms if len(a.bond_partners) == 1
            ]
            terminal_atoms = self._sort_atoms_by_mass(terminal_atoms)
        root = terminal_atoms[0].idx

        # Find (parent, child) bonds in BFS order starting from root
        self.acyclic_graph, self.bonds = self._build_bonds_bfs(root)

        # Create a mapping between the original indices and the BFS-explored indices
        # Again, Molmodel adds atoms via Molmodel via bondAtom(idx1, idx2)
        # Essentially, we create a mapping between the order in which atoms were added to  and their original indices
        self.nodes = [root]
        for parent_local_index, child_local_index, dihedral_type, resid in self.bonds:
            if "ring" not in dihedral_type:
                self.nodes.append(child_local_index)

        self.local_to_compound_atom_index_map = {}
        for compound_atom_index, prmtop_index in enumerate(self.nodes):
            self.local_to_compound_atom_index_map[prmtop_index] = compound_atom_index

        # Generate Z matrix atom indices
        self.z_matrix = self._build_z_matrix()

        # Build atom parameters
        self.atom_params: list[AtomParams] = []
        self.num_residues = 0

        for local_index in self.nodes:
            a = self.molecule.atoms[local_index]
            self.atom_params.append(
                AtomParams(
                    local_index=local_index,
                    compound_atom_index=self._local_to_compound_atom_index(local_index),
                    element_name=a.element_name,
                    element_symbol=a.element_name,  # TODO
                    atomic_number=a.atomic_number,
                    charge_e=a.ucharge.value_in_unit(pmd.unit.elementary_charge),
                    mass_daltons=a.umass.value_in_unit(pmd.unit.dalton),
                    vdw_radius_nm=a.urmin.value_in_unit(pmd.unit.nanometer) * 2,
                    vdw_well_depth_kj=a.uepsilon.value_in_unit(
                        pmd.unit.kilojoule_per_mole
                    ),
                    sigma_nm=a.usigma.value_in_unit(pmd.unit.nanometer),
                    solvent_radius_nm=a.usolvent_radius.value_in_unit(
                        pmd.unit.nanometer
                    ),
                    screen=a.screen,
                    neighbors_local_indices=[node.idx for node in a.bond_partners],
                    root=(local_index == root),
                )
            )
            self.num_residues = max(self.num_residues, a.residue.idx + 1)

        # Now we create bond parameters
        # We must canonicalize the bond indices to make sure that (i, j) and (j, i) map to the same parameter set
        # This is because we reorder bond direction (parent, child) during BFS traversal
        # We begin by creating a lookup table
        bond_lookup = {
            frozenset((b.atom1.idx, b.atom2.idx)): b for b in self.molecule.bonds
        }

        # Build bond parametersBondParams
        self.bond_params: List[BondParams] = []
        for parent_local_index, child_local_index, dihedral_type, resid in self.bonds:
            bond = bond_lookup.get(frozenset((parent_local_index, child_local_index)))

            if bond is None:
                raise RuntimeError(
                    f"Internal Error: Bond between atoms {parent_local_index} and {child_local_index} not found in bond lookup."
                )

            self.bond_params.append(
                BondParams(
                    parent_local_index=parent_local_index,
                    child_local_index=child_local_index,
                    parent_compound_atom_index=self._local_to_compound_atom_index(
                        parent_local_index
                    ),
                    child_compound_atom_index=self._local_to_compound_atom_index(
                        child_local_index
                    ),
                    stiffness_in_kj_per_nm_sq=bond.type.uk.value_in_unit(
                        pmd.unit.kilojoule_per_mole / pmd.unit.nanometer**2
                    ),
                    nominal_length_in_nm=bond.type.ureq.value_in_unit(
                        pmd.unit.nanometer
                    ),
                    is_ring_closing=("ring" in dihedral_type),
                    dihedral_type=dihedral_type,
                )
            )

        # Build angle parameters
        self.angle_params: List[AngleParams] = []
        for angle in self.molecule.angles:
            self.angle_params.append(
                AngleParams(
                    local_indices=(angle.atom1.idx, angle.atom2.idx, angle.atom3.idx),
                    compound_atom_indices=(
                        self._local_to_compound_atom_index(angle.atom1.idx),
                        self._local_to_compound_atom_index(angle.atom2.idx),
                        self._local_to_compound_atom_index(angle.atom3.idx),
                    ),
                    stiffness_in_kj_per_rad_sq=angle.type.uk.value_in_unit(
                        pmd.unit.kilojoule_per_mole / pmd.unit.radian**2
                    ),
                    nominal_angle_in_deg=angle.type.utheteq.value_in_unit(
                        pmd.unit.degree
                    ),
                )
            )

        # Build AMBER proper and improper periodic torsions
        # They have multiple terms, each with their own amplitude, phase and periodicity
        periodic_torsion_terms: dict[
            tuple[int, int, int, int], list[rb.RoboPeriodicTorsionTerm]
        ] = {}
        for d in self.molecule.dihedrals:
            # Phase is sensitive to floating point errors
            EPS = 1e-3
            phase_deg = d.type.uphase.value_in_unit(pmd.unit.degree)
            phase_deg = ((phase_deg + 180) % 360) - 180
            if phase_deg > 180.0 + EPS:
                raise ValueError(
                    f"Torsion phase {phase_deg:.6f}° exceeds 180 beyond tolerance ({EPS}°)"
                )
            else:
                phase_deg = 180.0 if phase_deg > 180.0 else phase_deg

            key = d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx
            term = rb.RoboPeriodicTorsionTerm(
                amplitude_kj=d.type.uphi_k.value_in_unit(pmd.unit.kilojoule_per_mole),
                phase_deg=phase_deg,
                periodicity=d.type.per,
            )

            if key not in periodic_torsion_terms:
                periodic_torsion_terms[key] = {"improper": d.improper, "terms": [term]}
            else:
                # Check that we don't mix proper and improper torsions for the same set of atoms
                if periodic_torsion_terms[key]["improper"] != d.improper:
                    raise ValueError(f"Mixed proper/improper torsion for key {key}")

                # Check that all periodicities are unique for this torsion
                existing_periodicities = {
                    t.periodicity for t in periodic_torsion_terms[key]["terms"]
                }
                if term.periodicity in existing_periodicities:
                    raise ValueError(
                        f"Duplicate periodicity {term.periodicity} found for torsion between atoms {key[0]}, {key[1]}, {key[2]}, {key[3]}"
                    )

                periodic_torsion_terms[key]["terms"].append(term)

            # TODO cannot be 0 periodicity?
            # TODO we also have amplitude 0 terms that we should ignore perhaps?
            # TODO from DuMM defineBondTorsion: Pay particular attention to the amplitude -- in our convention it is really a half-amplitude since the sinusoids range from -1 to 1 making the full energy and torque "excursion" twice the amplitude.
            # TODO check the other functions as well, they have some cautions about amplitude signs etc
            if d.type.per == 0:
                raise ValueError(
                    f"Torsion term with 0 periodicity found for torsion between atoms {key[0]}, {key[1]}, {key[2]}, {key[3]}"
                )

        # Now that we have collected all periodic torsion terms, build the parameter objects
        self.periodic_torsion_params: List[PeriodicTorsionParams] = []
        for key, value in periodic_torsion_terms.items():
            (
                atom1_local_index,
                atom2_local_index,
                atom3_local_index,
                atom4_local_index,
            ) = key

            self.periodic_torsion_params.append(
                PeriodicTorsionParams(
                    local_indices=(
                        atom1_local_index,
                        atom2_local_index,
                        atom3_local_index,
                        atom4_local_index,
                    ),
                    compound_atom_indices=(
                        self._local_to_compound_atom_index(atom1_local_index),
                        self._local_to_compound_atom_index(atom2_local_index),
                        self._local_to_compound_atom_index(atom3_local_index),
                        self._local_to_compound_atom_index(atom4_local_index),
                    ),
                    is_improper=value["improper"],
                    terms=value["terms"],
                )
            )

        # Build CHARMM improper harmonic torsions
        # They are absent in AMBER prmtop files
        self.improper_harmonic_torsion_terms: List[HarmonicImproperTorsionParams] = []
        for imp in self.molecule.impropers:
            self.improper_harmonic_torsion_terms.append(
                HarmonicImproperTorsionParams(
                    local_indices=(
                        imp.atom1.idx,
                        imp.atom2.idx,
                        imp.atom3.idx,
                        imp.atom4.idx,
                    ),
                    compound_atom_indices=(
                        self._local_to_compound_atom_index(imp.atom1.idx),
                        self._local_to_compound_atom_index(imp.atom2.idx),
                        self._local_to_compound_atom_index(imp.atom3.idx),
                        self._local_to_compound_atom_index(imp.atom4.idx),
                    ),
                    stiffness_in_kj_per_rad_sq=imp.type.upsi_k.value_in_unit(
                        pmd.unit.kilojoule_per_mole / pmd.unit.radian**2
                    ),
                    nominal_angle_in_rad=imp.type.upsi_eq.value_in_unit(
                        pmd.unit.radian
                    ),
                )
            )

    def _local_to_compound_atom_index(self, local_index: int) -> int:
        """
        Map an Amber prmtop atom index [a, b] to the internal compound representation atom index [0, n-1] where n is the number of atoms in the compound.

        Parameters
        ----------
        prmtop_index : int
            The prmtop 0-based atom index to map.

        Returns
        -------
        int
            The corresponding compound atom index in range [0, n-1].
        """
        try:
            return self.local_to_compound_atom_index_map[local_index]

        # At the time of writing this function, the map is a dictionary, but we catch other types for future-proofing
        except (KeyError, IndexError, TypeError) as e:
            # Determine the physical bounds of your map for the error report
            lower_bound = min(self.local_to_compound_atom_index_map.keys())
            upper_bound = max(self.local_to_compound_atom_index_map.keys())
            total_keys = len(self.local_to_compound_atom_index_map)

            raise IndexError(
                f"\n--- TOPOLOGY MAPPING FATAL ERROR ---\n"
                f"Input Index:      {local_index} ({type(local_index).__name__})\n"
                f"Expected Domain:  [{lower_bound}, {upper_bound}], total keys: {total_keys}\n"
                f"Internal Error:   {type(e).__name__}: {e}\n"
                f"Resolution:       Verify that your prmtop file matches the current \n"
                f"                  system state and that you aren't using 1-based \n"
                f"                  indexing for a 0-based lookup.\n"
                f"------------------------------------"
            ) from e

    def _build_bonds_macrocycles(self) -> List[Tuple[int, int]]:
        # Create a full, undirected graph of the molecule
        full_graph = nx.Graph()
        for bond in self.molecule.bonds:
            full_graph.add_edge(bond.atom1.idx, bond.atom2.idx)

        bonds = []
        for cycle in nx.cycle_basis(full_graph):
            # Macrocycles are defined as cycles of length >= 12
            if len(cycle) < 12:
                continue

            # Assign a standardized dihedral type to each bond in the cycle based on its neighboring atoms
            cycle_bonds: List[Tuple[int, int, str]] = []
            for i in range(len(cycle)):
                # Get the current atom and the next one (using modulo for the wrap-around)
                atom_a = cycle[i]
                atom_b = cycle[(i + 1) % len(cycle)]
                dihedral_type = self._get_standardized_dihedral_type(
                    self.molecule[atom_a], self.molecule[atom_b]
                )

                # Append the pair to your bonds list
                cycle_bonds.append([atom_a, atom_b, dihedral_type])
            bonds.append(cycle_bonds)

        return bonds

    def _build_bonds_bfs(
        self,
        root: pmd.Atom,
        forbidden_atom_type_pairs: Iterable[Tuple[str, str]] = [("SG", "SG")],
    ) -> List[Tuple[int, int, bool]]:
        """
        Build an acyclic graph of the molecule and order bonds for Z-matrix generation.

        Non-ring-closing bonds appear first in BFS traversal order. Ring-closing
        bonds are added afterward, respecting forbidden atom-type constraints.

        Parameters
        ----------
        forbidden_atom_type_pairs : iterable of (str, str)
            Atom types between which ring-closing bonds must be avoided.

        Returns
        -------
        List of tuples:
            (parent_local_index, child_local_index, is_ring_closing)
        """
        # Acyclic graph
        acyclic_graph = nx.Graph()
        acyclic_graph.add_nodes_from(atom.idx for atom in self.molecule.atoms)

        # Count how many ring closing bonds each atom is involved in to enforce the constraint that no atom can be involved in more than one ring closing bond
        ring_closing_bonds_involved: Dict[int, int] = {}

        # Set of ring-closing bonds for quick lookup
        ring_closing_bonds: Set[Tuple[int, int]] = set()

        # First pass: identify standardized ring-closing dihedrals
        for bond in self.molecule.bonds:
            parent, child = bond.atom2, bond.atom1
            parent_idx, child_idx = parent.idx, child.idx

            # Check if this is the middle bond of a standardized dihedral
            # Non-standardized dihedral are None and we treat as non-ring-closing for now
            dihedral_type = self._get_standardized_dihedral_type(parent, child)

            if "ring" in dihedral_type:
                # Check rigid bond constraints
                for idx in (parent_idx, child_idx):
                    if ring_closing_bonds_involved.get(idx, 0) >= 1:
                        raise ValueError(
                            f"Atom {idx} involved in multiple rigid bonds."
                        )
                    ring_closing_bonds_involved[idx] = (
                        ring_closing_bonds_involved.get(idx, 0) + 1
                    )
                ring_closing_bonds.add((parent_idx, child_idx))
            else:
                acyclic_graph.add_edge(parent_idx, child_idx)

        # Compute forbidden bonds and distances
        forbidden_edges = self._find_forbidden_bonds(forbidden_atom_type_pairs)
        forbidden_nodes = GraphTraversalUtils.bond_edges_to_nodes(forbidden_edges)
        forbidden_atom_dist = (
            GraphTraversalUtils.nodes_to_distances(acyclic_graph, forbidden_nodes)
            if forbidden_nodes
            else {}
        )

        # Handle residual cycles
        for cycle in sorted(nx.cycle_basis(acyclic_graph), key=len):
            candidate_bonds: List[BondProperties] = []
            n = len(cycle)
            for i in range(n):
                u, v = cycle[i], cycle[(i + 1) % n]
                if frozenset((u, v)) in forbidden_edges:
                    continue

                atom_u, atom_v = self.molecule[u], self.molecule[v]
                degree_sum = acyclic_graph.degree[u] + acyclic_graph.degree[v]
                is_rigid = is_rigid_bond(atom_u, atom_v)

                forbidden_distance = (
                    max(forbidden_atom_dist.get(u, 0), forbidden_atom_dist.get(v, 0))
                    if forbidden_atom_dist
                    else None
                )

                candidate_bonds.append(
                    BondProperties(u, v, is_rigid, forbidden_distance, degree_sum)
                )

            # Sort bonds: rigid first, farthest from forbidden atoms, lowest degree sum
            candidate_bonds.sort(
                key=lambda b: (
                    not b.is_rigid,
                    -(b.sg_sg_bond_distance or -1),
                    b.sum_of_degrees,
                )
            )

            # Pick first eligible bond
            for bond in candidate_bonds:
                u, v = bond.parent_atom_prmtop_index, bond.child_atom_prmtop_index
                if (u, v) not in ring_closing_bonds and (
                    v,
                    u,
                ) not in ring_closing_bonds:
                    break
            else:
                continue  # no bond eligible

            # Check rigid bond constraints
            for idx in (u, v):
                if ring_closing_bonds_involved.get(idx, 0) >= 1:
                    raise ValueError(f"Atom {idx} involved in multiple rigid bonds.")
                ring_closing_bonds_involved[idx] = (
                    ring_closing_bonds_involved.get(idx, 0) + 1
                )

            acyclic_graph.remove_edge(u, v)
            ring_closing_bonds.add((u, v))

        # Validation
        if acyclic_graph.number_of_nodes() != len(self.molecule.atoms):
            raise ValueError("Graph contains unknown nodes.")
        if not nx.is_forest(acyclic_graph):
            raise ValueError("Graph contains residual cycles.")
        if not nx.is_connected(acyclic_graph):
            raise ValueError("Graph disconnected after ring closure removal.")

        # Generate BFS order for non-ring-closing bonds
        bond_to_resid = {
            frozenset((b.atom1.idx, b.atom2.idx)): b.atom1.residue.idx
            for b in self.molecule.bonds
        }
        edges_with_flags = [
            (
                u,
                v,
                self._get_standardized_dihedral_type(
                    self.molecule[u], self.molecule[v]
                ),
                bond_to_resid.get(frozenset((u, v)), -1),
            )
            for u, v in nx.bfs_edges(acyclic_graph, source=root)
        ]

        # Append ring-closing bonds (order irrelevant)
        edges_with_flags += [
            (u, v, "ring-closing", bond_to_resid.get(frozenset((u, v)), -1))
            for u, v in ring_closing_bonds
        ]

        # Make sure we have all bonds accounted for
        if len(edges_with_flags) != len(self.molecule.bonds):
            raise ValueError("Mismatch in number of bonds after processing.")

        # We are implicitly constructing a spanning tree used to define internal coordinates
        # We need to validate parents have already been visited and children have not been visited yet
        GraphTraversalUtils.validate_bfs_parent_child_edges(
            acyclic_graph, root, edges_with_flags
        )

        return acyclic_graph, edges_with_flags

    @staticmethod
    def _sort_atoms_by_mass(atoms: list[pmd.Atom]) -> list[pmd.Atom]:
        r"""Sorts a list of atoms in descending order by mass.
        The atom index is used as a tiebreaker so that the ordering is reproducible.
        If two atoms have the same mass, the one with the lower index comes first.

        Parameters
        ----------
        atoms : list of pmd.Atom
            List of atoms to sort.

        Returns
        -------
        sorted_atoms : list of pmd.Atom
            Sorted list.
        """

        # Negate mass to sort it descending, keep idx as is for ascending
        return sorted(atoms, key=lambda a: (-a.mass, a.idx))

    def _get_standardized_dihedral_type(
        self, parent_atom: pmd.Atom, child_atom: pmd.Atom
    ) -> str:
        """
        Determine the standardized dihedral type for a given parent-child atom pair.

        This function generates candidate dihedrals by combining each bond partner
        of the parent with each bond partner of the child, excluding trivial overlaps.
        Each candidate is a `parmed.topologyobjects.Dihedral` and is classified via
        `self.dihedral_classifier`. The first match found is returned.

        Parameters
        ----------
        parent_atom : pmd.Atom
            Atom serving as the "parent" in the dihedral.
        child_atom : pmd.Atom
            Atom serving as the "child" in the dihedral.

        Returns
        -------
        str
            The standardized dihedral type.
        """
        candidates = [  # noqa
            pmd.Dihedral(grandparent, parent_atom, child_atom, gchild)
            for grandparent in parent_atom.bond_partners
            if grandparent != child_atom
            for gchild in child_atom.bond_partners
            if gchild != parent_atom
        ]

        for gparent in parent_atom.bond_partners:
            for gchild in child_atom.bond_partners:
                if gparent == child_atom or gchild == parent_atom:
                    continue
                dihedral_type = self.dihedral_classifier.classify(  # noqa
                    gparent, parent_atom, child_atom, gchild
                )

        return "non-standard"

    def _find_forbidden_bonds(
        self, atom_type_pairs: Iterable[Tuple[str, str]]
    ) -> Set[FrozenSet[int]]:
        """
        Identify all bonds in a molecule whose endpoint atom names match any of
        the specified forbidden atom-type pairs.

        Parameters
        ----------
        atom_type_pairs : iterable of (str, str)
            Forbidden atom-type pairs (e.g. [('SG', 'SG'), ('ZN', 'SG')]).
            Matching is case-insensitive and order-independent.

        Returns
        -------
        forbidden_edges : set of frozenset(int, int)
            Set of forbidden bonds represented as frozensets of atom indices.
        """
        # Normalize pairs for case-insensitive, order-independent matching
        forbidden_pairs = {
            frozenset((a.lower(), b.lower())) for a, b in atom_type_pairs
        }

        forbidden_edges: Set[FrozenSet[int]] = set()

        for bond in self.molecule.bonds:
            a1, a2 = bond.atom1, bond.atom2
            atom_pair = frozenset((a1.name.lower(), a2.name.lower()))

            if atom_pair in forbidden_pairs:
                forbidden_edges.add(frozenset((a1.idx, a2.idx)))

        return forbidden_edges

    def _build_z_matrix(self) -> list[ZMatrixRowAtomGroups]:
        """
        Construct a Z-matrix (BAT-style internal coordinate tree) from an MDAnalysis Universe molecule.
        Cycles are handled by explicitly excluding ring-closing bonds, turning the molecular graph into a spanning tree.
        This function expects that the molecular graph is acyclic (a forest) after removing the pre-defined ring-closing bonds.
        Here, we try to select root atoms that are heavy and non-terminal to improve numerical stability:
        1. The first atom is chosen as the heaviest terminal atom.
        2. The second atom is the its bonded neighbor.
        3. The third atom is the heaviest non-terminal bonded neighbor of the second atom that is not part of a ring-closing bond and is not collinear with the first two atoms.

        The result is torsional spanning tree while explicitly removing redundant ring closing bonds.

        Parameters
        ----------
        ag_o : list of MDAnalysis Atoms
            List to sort
        reverse : bool
            Atoms will be in descending order

        Returns
        -------
        ag_n : list of Atoms
            an ordered, loop-free internal coordinate definition
        """

        # The molecular graph must be acyclic for the Z-matrix to be built correctly
        # This is ensured by removing ring-closing bonds beforehand (see above)
        # Basically, we check here that our ring closing bond detection worked correctly and we did not miss any cycles
        # A forest is a graph with no undirected cycles
        if not nx.is_forest(self.acyclic_graph):
            raise ValueError("Molecule contains cycles; cannot build Z-matrix.")

        # Build the Z matrix
        z_matrix: list[ZMatrixRowAtomGroups] = []

        # We begin by excluding degenerate cases with 1, 2, or 3 atoms
        if len(self.molecule.atoms) == 0:
            raise ValueError("Molecule contains no atoms; cannot build Z-matrix.")

        if len(self.molecule.atoms) == 1:
            initial_atom = self.molecule.atoms[0]
            z_matrix.append(
                ZMatrixRowAtomGroups(
                    i_global_atom_index=initial_atom.idx,
                    i_compound_atom_index=self._local_to_compound_atom_index(
                        initial_atom.idx
                    ),
                )
            )
            return z_matrix

        if len(self.molecule.atoms) == 2:
            initial_atom = self.molecule.atoms[0]
            z_matrix.append(
                ZMatrixRowAtomGroups(
                    i_global_atom_index=initial_atom.idx,
                    i_compound_atom_index=self._local_to_compound_atom_index(
                        initial_atom.idx
                    ),
                )
            )

            second_atom = self.molecule.atoms[1]
            z_matrix.append(
                ZMatrixRowAtomGroups(
                    i_global_atom_index=second_atom.idx,
                    i_compound_atom_index=self._local_to_compound_atom_index(
                        second_atom.idx
                    ),
                    j_global_atom_index=initial_atom.idx,
                    j_compound_atom_index=self._local_to_compound_atom_index(
                        initial_atom.idx
                    ),
                )
            )

            return z_matrix

        # We select the initial atom from terminal atoms (atoms with only one bond)
        # TODO Macrocycles with substituents may give terminals that lie on side chains, producing pathological Z-matrices.
        # TODO Root selection must be graph-theoretic, not terminal-based. Use: degree-1 atoms if they exist, otherwise fall back to highest-degree or highest-mass atom not in a ring closure set.
        terminal_atoms = [a for a in self.molecule.atoms if len(a.bonds) == 1]
        terminal_atoms = self._sort_atoms_by_mass(terminal_atoms)

        # Select the heaviest root atom from the heaviest terminal atoms
        initial_atom = terminal_atoms[0]
        z_matrix.append(
            ZMatrixRowAtomGroups(
                i_global_atom_index=initial_atom.idx,
                i_compound_atom_index=self._local_to_compound_atom_index(
                    initial_atom.idx
                ),
            )
        )

        # The next atom in the root is bonded to the initial atom
        # Since the initial atom is a terminal atom, there is only one bonded atom
        second_atom = initial_atom.bond_partners[0]
        z_matrix.append(
            ZMatrixRowAtomGroups(
                i_global_atom_index=second_atom.idx,
                i_compound_atom_index=self._local_to_compound_atom_index(
                    second_atom.idx
                ),
                j_global_atom_index=initial_atom.idx,
                j_compound_atom_index=self._local_to_compound_atom_index(
                    initial_atom.idx
                ),
            )
        )

        # The last atom in the root is the heaviest atom bonded to the second atom
        # If there are more than three atoms, then the last atom cannot be a terminal atom
        # Moreover, the three atoms must not be collinear
        if len(self.molecule.atoms) != 3:
            third_atom_candidates = []
            for a in self._tree_neighbors(second_atom):
                if (
                    (a != initial_atom) and (a not in terminal_atoms)
                ):  # and not self._are_collinear(initial_atom.position, second_atom.position, a.position)
                    third_atom_candidates.append(a)
            third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]
        else:
            third_atom_candidates = []
            for a in self._tree_neighbors(second_atom):
                if a != initial_atom:
                    third_atom_candidates.append(a)
            third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]

        third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]
        z_matrix.append(
            ZMatrixRowAtomGroups(
                i_global_atom_index=third_atom.idx,
                i_compound_atom_index=self._local_to_compound_atom_index(
                    third_atom.idx
                ),
                j_global_atom_index=second_atom.idx,
                j_compound_atom_index=self._local_to_compound_atom_index(
                    second_atom.idx
                ),
                k_global_atom_index=initial_atom.idx,
                k_compound_atom_index=self._local_to_compound_atom_index(
                    initial_atom.idx
                ),
            )
        )

        # Root triplet
        root = [initial_atom, second_atom, third_atom]

        for atom_i, atom_j, atom_k, atom_l in self._find_torsions(root):
            z_matrix.append(
                ZMatrixRowAtomGroups(
                    i_global_atom_index=atom_i.idx,
                    i_compound_atom_index=self._local_to_compound_atom_index(
                        atom_i.idx
                    ),
                    j_global_atom_index=atom_j.idx,
                    j_compound_atom_index=self._local_to_compound_atom_index(
                        atom_j.idx
                    ),
                    k_global_atom_index=atom_k.idx,
                    k_compound_atom_index=self._local_to_compound_atom_index(
                        atom_k.idx
                    ),
                    l_global_atom_index=atom_l.idx,
                    l_compound_atom_index=self._local_to_compound_atom_index(
                        atom_l.idx
                    ),
                )
            )

        return z_matrix

    def _tree_neighbors(self, atom: pmd.Atom) -> list[pmd.Atom]:
        """
        Builds the list atoms bonded to the given atom, excluding those involved in ring-closing bonds.

        Parameters
        ----------
        atom : MDAnalysis Atom
            The atom for which to find tree neighbors.

        Returns
        -------
        neighbors : list of MDAnalysis Atom
            List of bonded atoms not involved in ring-closing bonds.
        """
        return [
            a
            for a in atom.bond_partners
            if "ring" not in self._get_standardized_dihedral_type(atom, a)
        ]

    def _find_torsions(
        self, root: list[pmd.Atom]
    ) -> list[tuple[pmd.Atom, pmd.Atom, pmd.Atom, pmd.Atom]]:
        """
        Constructs a list of torsion angles.

        Returns
        -------
        torsions : list of AtomGroup
            list of AtomGroup objects that define torsion angles
        """
        torsions: list[tuple[pmd.Atom, pmd.Atom, pmd.Atom, pmd.Atom]] = []
        selected_atoms = list(root)
        selected_atoms_set = set([root[0].idx, root[1].idx, root[2].idx])

        # Build tree neighbors
        tree_adj: dict[int, list[pmd.Atom]] = {
            # atom.idx: self._sort_atoms_by_mass([a for a in atom.bond_partners if not self.is_ring_closing(atom, a)]) for atom in self.molecule.atoms
            atom.idx: self._sort_atoms_by_mass([a for a in atom.bond_partners])
            for atom in self.molecule.atoms
        }

        while len(selected_atoms) < len(self.molecule.atoms):
            torsionAdded = False
            for a1 in selected_atoms:
                # Find a0, which is a new atom connected to the selected atom
                a0_list = [
                    a for a in tree_adj[a1.idx] if a.idx not in selected_atoms_set
                ]
                for a0 in a0_list:
                    # Find a2, which is connected to a1, is not a terminal atom and has been selected
                    a2_list = [
                        a
                        for a in tree_adj[a1.idx]
                        if (a != a0)
                        and len(a.bond_partners) > 1
                        and (a.idx in selected_atoms_set)
                    ]
                    for a2 in a2_list:
                        # Find a3, which is connected to a2, has been selected, and is not a1
                        a3_list = [
                            a
                            for a in tree_adj[a2.idx]
                            if (a != a1) and (a.idx in selected_atoms_set)
                        ]
                        for a3 in a3_list:
                            # Add the torsion to the list of torsions
                            torsions.append((a0, a1, a2, a3))

                            # Add the new atom to selected_atoms which extends the loop
                            selected_atoms.append(a0)
                            selected_atoms_set.add(a0.idx)
                            torsionAdded = True
                            break
                        break

            if torsionAdded is False:
                print("Selected atoms:")
                print([a.idx + 1 for a in selected_atoms])
                print("Torsions found:")
                print([list(t.indices + 1) for t in torsions])
                raise ValueError("Additional torsions not found.")

        return torsions
