import parmed as pmd
from collections import deque
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import List, Iterable, Tuple, Set, FrozenSet, Optional, Dict, Hashable
import robo_bindings as rb
import numpy as np
import networkx as nx

class GraphTraversalUtils:
    @staticmethod
    def validate_bfs_parent_child_edges(graph: nx.Graph, root: Hashable, bfs_edges: Iterable[Tuple[Hashable, Hashable]]) -> None:
        """
        Validate that a sequence of edges corresponds to a valid BFS tree
        traversal starting from a given root.

        This function enforces the invariant that for every edge (u, v):
        - u has already been visited when the edge is produced
        - v has not been visited before (i.e., v is discovered via u)
        - (u, v) is an actual edge in the graph

        Parameters
        ----------
        graph : networkx.Graph
            The graph on which BFS is assumed to have been performed.
        root : hashable
            The BFS root node.
        bfs_edges : iterable of (hashable, hashable)
            Edges produced by a BFS traversal, typically from
            ``networkx.bfs_edges``.

        Raises
        ------
        ValueError
            If any edge violates BFS parent-child semantics.
        """
        visited = {root}

        for step, (u, v) in enumerate(bfs_edges):
            if u not in visited:
                raise ValueError(f"BFS invariant violated at step {step}: parent node {u} has not been visited yet.")
            if v in visited:
                raise ValueError(f"BFS invariant violated at step {step}: child node {v} was already visited.")
            if not graph.has_edge(u, v):
                raise ValueError(f"BFS invariant violated at step {step}: edge ({u}, {v}) does not exist in graph.")
            visited.add(v)


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
    local_indices: Tuple[int, int]
    compound_atom_indices: Tuple[int, int]
    stiffness_in_kj_per_nm_sq: float
    nominal_length_in_nm: float
    is_ring_closing: bool

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

# Hold compound atom indices, not prmtop indices
@dataclass
class ZMatrixRowAtomGroups:
    i: pmd.Atom
    j: pmd.Atom | None
    k: pmd.Atom | None
    l: pmd.Atom | None

# def bond_ag(row: ZMatrixRowAtomGroups) -> mda.AtomGroup | None:
#     return None if row.j is None else mda.AtomGroup([row.i, row.j])

# def angle_ag(row: ZMatrixRowAtomGroups, u: mda.Universe) -> mda.AtomGroup | None:
#     return None if row.k is None else u.atoms[[row.i, row.j, row.k]]

# def dihedral_ag(row: ZMatrixRowAtomGroups, u: mda.Universe) -> mda.AtomGroup | None:
#     return None if row.l is None else u.atoms[[row.i, row.j, row.k, row.l]]




class Molecule(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def prmtop_to_compound_atom_index(self, prmtop_index: int) -> int:
        pass

    @abstractmethod
    def get_nodes_as_prmtop_indices(self) -> List[int]:
        pass

    @abstractmethod
    def get_bonds(self) -> List[Tuple[int, int, bool]]:
        pass







# For pretty images: https://emleddin.github.io/comp-chem-website/AMBERguide-AAs-DNA-RNA.html
PROTEIN_DIHEDRAL_SELECTION = {
    # N-terminus ACE
    'ACE': {},

    # C-terminus NME
    'NME': {},

    # Alanine
    'ALA': {
        'chi1': ('N', 'CA', 'CB', 'HB1'),
    },

    # Arginine as a terminal large, resonance-stabilized, planar structure with a diffuse positive charge guanidinium group
    # I.e CZ connected to to two nitrogen atoms (NH1, NH2) which are connected to two hydrogens each (HH11, HH12 and HH21, HH22)
    'ARG': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'NE'),
        'chi4': ('CG', 'CD', 'NE', 'CZ'),
        # 'chi5': ('CD', 'NE', 'CZ', 'NH1'),
    },

    # Asparagine
    'ASN': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND2'),
        'chi3': ('CB', 'CG', 'ND2', 'HD21'),
    },

    # Aspartate has a terminal rigid planar carboxylate group: one carbon (CG) connected to two oxygens (OD1, OD2)
    # It can be unprotonated (ASP, charge -1) or protonated (ASH, charge 0) in which there is one hydrogen atom connected to one of the oxygens
    # We treat the terminus as rigid, similar to ARG
    'ASP': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'OD1'),
    },
    'ASH': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'OD1'),
    },

    # Cysteine (CYS) has a terminal thiol group (HG connected to SG)
    # If HG is lost, it can become CYX and form disulfide bonds or CYM and coordinate metals
    'CYS': {
        'chi1': ('N', 'CA', 'CB', 'SG'),
        'chi2': ('CA', 'CB', 'SG', 'HG'),
    },
    'CYX': {
        'chi1': ('N', 'CA', 'CB', 'SG'),
    },
    'CYM': {
        # TODO check if metal coordination induces rigid dihedrals
        'chi1': ('N', 'CA', 'CB', 'SG'),
    },

    # Glutamine
    'GLN': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'NE2'),
        'chi4': ('CG', 'CD', 'NE2', 'HE21'),
    },

    # Glutamate has a terminal rigid planar carboxylate group: one carbon (CD) connected to two oxygens (OE1, OE2)
    # It can be unprotonated (GLU, charge -1) or protonated (GLH, charge 0) in which there is one hydrogen atom connected to one of the oxygens
    # We treat the terminus as rigid, similar to ARG
    'GLU': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'OE1'),
    },
    'GLH': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'OE1'),
    },

    # Glycine has no side chain
    'GLY': {},

    # Imidazole (HID/HIE/HIP) is aromatic and essentially planar
    # TODO does it matter from a kinematic point of view where we place the ring closing bond?
    'HIP': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND1'),
        'ring_closing_1': ('CG', 'ND1', 'CE1', 'NE2'),
    },
    'HIE': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND1'),
        'ring_closing_1': ('CG', 'ND1', 'CE1', 'NE2'),
    },
    'HID': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND1'),
        'ring_closing_1': ('CG', 'ND1', 'CE1', 'NE2'),
    },

    # Isoleucine
    'ILE': {
        'chi1': ('N', 'CA', 'CB', 'CG1'),
        'chi2.1': ('CA', 'CB', 'CG1', 'CD1'),
        'chi2.2': ('CA', 'CB', 'CG2', 'HG21'),
        'chi3': ('CB', 'CG1', 'CD1', 'HD11'),
    },

    # Leucine
    'LEU': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD1'),
        'chi3.1': ('CB', 'CG', 'CD1', 'HD11'),
        'chi3.2': ('CB', 'CG', 'CD2', 'HD21'),
    },

    # Lysine has a long flexible side chain ending in a positively charged amino group (NZ connected to HZ1, HZ2, HZ3)
    # If one of the hydrogens is lost, it can become LYN (neutral)
    # This does not affect the chi angles
    'LYS': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'CE'),
        'chi4': ('CG', 'CD', 'CE', 'NZ'),
        'chi5': ('CD', 'CE', 'NZ', 'HZ1'),
    },
    'LYN': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'CE'),
        'chi4': ('CG', 'CD', 'CE', 'NZ'),
        'chi5': ('CD', 'CE', 'NZ', 'HZ1'),
    },

    # Methionine
    'MET': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'SD'),
        'chi3': ('CB', 'CG', 'SD', 'CE'),
        'chi4': ('CG', 'SD', 'CE', 'HE1'),
    },

    # Phenylalanine
    'PHE': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD1'),
        'ring_closing_1': ('CG', 'CD1', 'CE1', 'CZ'), # Same as TYR
    },

    # In prooline, there are 5 chi angles
    # Chi1 and Chi2 are the most important for puckering
    # Chi5 is equivalent to phi
    # Chi4 is defined as ring closing
    'PRO': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD'),
        'chi3': ('CB', 'CG', 'CD', 'N'),
        'ring_closing_1': ('CG', 'CD', 'N', 'CA'),
    },

    # Serine
    'SER': {
        'chi1': ('N', 'CA', 'CB', 'OG'),
        'chi2': ('CA', 'CB', 'OG', 'HG'),
    },

    # Threonine
    'THR': {
        'chi1': ('N', 'CA', 'CB', 'OG1'),
        'chi2.1': ('CA', 'CB', 'OG1', 'HG1'),
        'chi2.2': ('CA', 'CB', 'CG2', 'HG21'),
    },

    # Tryptophan's side chain is a bicyclic structure called the indole group, which is composed of two fused rings: benzene and pyrrole
    'TRP': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD1'),
        'ring_closing_1': ('CG', 'CD1', 'NE1', 'CE2'), # Pyrrole ring
        'ring_closing_2': ('CE2', 'CZ2', 'CH2', 'CZ3'), # Benzene ring
    },

    # Tyrosine
    'TYR': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD1'),
        'chi3': ('CE1', 'CZ', 'OH', 'HH'),
        'ring_closing_1': ('CG', 'CD1', 'CE1', 'CZ'), # Same as PHE
    },

    # Valine
    'VAL': {
        'chi1': ('N', 'CA', 'CB', 'CG1'),
        'chi2.1': ('CA', 'CB', 'CG1', 'HG11'),
        'chi2.2': ('CA', 'CB', 'CG2', 'HG21'),
    }
}

AROMATIC_TYPES = {
    # GAFF aromatic carbons / heteroatoms
    "ca", "cp", "cq", "cc", "cd", "ce", "cf",
    "na", "nb", "nc", "nd", "ne", "nf",

    # ff19SB aromatic
    "c", "ca", "cb", "cc", "cd", "ce", "cf", "cw", "cr"
}

# sp / sp2 carbons and heteroatoms (very conservative)
MULTIPLE_BOND_TYPES = {
    # carbonyls, sp2/sp
    "c", "c1", "c2", "ce", "cf", "cg",
    "o", "o2", "os",
    "n", "n2", "n1",
}

# canonical AMBER amide pattern
AMIDE_C = {"c"}      # carbonyl carbon
AMIDE_N = {"n", "nh"}  # amide nitrogens

def is_aromatic(a1 : str, a2: str) -> bool:
    return (a1 in AROMATIC_TYPES) and (a2 in AROMATIC_TYPES)

def is_amide(a1: str, a2: str) -> bool:
    return ((a1 in AMIDE_C and a2 in AMIDE_N) or
            (a2 in AMIDE_C and a1 in AMIDE_N))

def is_multiple_like(a1: str, a2: str) -> bool:
    return (a1 in MULTIPLE_BOND_TYPES and
            a2 in MULTIPLE_BOND_TYPES)

def is_rigid_bond(parent_atom: pmd.Atom, child_atom: pmd.Atom) -> bool:
	a1_type = parent_atom.type.lower()
	a2_type = child_atom.type.lower()

	if is_aromatic(a1_type, a2_type):
		return True
	if is_amide(a1_type, a2_type):
		return True
	if is_multiple_like(a1_type, a2_type):
		return True
	if parent_atom.name == 'SG' and child_atom.name == 'SG':
		return True

	return False







class MoleculeTest:
    _COLLINEARITY_TOL = 1e-12
    _EPS = 1e-15

    def __init__(self,
                 molecule: pmd.Structure,
                 include_omega: bool = False,
                 include_chi1: bool = True,
                 include_chi2: bool = True,
                 include_chi3: bool = True,
                 include_chi4: bool = True,
                 include_chi5: bool = True,
                 want_n_terminus_phi_rigid: bool = True,
                 want_c_terminus_psi_rigid: bool = True):
        
        self.molecule = molecule
        
        # Build a list of standardized ring-closing dihedrals to exclude
        standardized_ring_closing_dihedral_types = set()
        for residue, dihedral_list in PROTEIN_DIHEDRAL_SELECTION.items():
            for dihedral_name, atom_types in dihedral_list.items():
                if 'ring_closing' in dihedral_name:
                    standardized_ring_closing_dihedral_types.add(atom_types)
        standardized_ring_closing_dihedral_types = list(standardized_ring_closing_dihedral_types)
        self.standardized_ring_closing_dihedral_ref_set = {min(t, t[::-1]) for t in standardized_ring_closing_dihedral_types}

        # Build a list of standardized dihedrals to include
        self.standardized_dihedrals = {
            # Protein backbone
            ("C", "N", "CA", "C") : 'omega',
            ("N", "CA", "C", "N") : 'phi',
            ("CA", "C", "N", "CA") : 'psi',
        }
        for residue, dihedral_list in PROTEIN_DIHEDRAL_SELECTION.items():
            for dihedral_name, atom_types in dihedral_list.items():
                if 'ring_closing' in dihedral_name:
                    continue
                self.standardized_dihedrals[atom_types] = dihedral_name
        self.standardized_dihedrals_ref_set = {min(t, t[::-1]) for t in self.standardized_dihedrals.keys()}

        # Save ring closing bonds using local prmtop indices
        # We are guaranteed to have no duplicates
        # Set is useful because we need to check if any bond is ring closing multiple times and this is O(1)
        self.ring_closing_bonds, acyclic_graph = self._find_ring_closing_bonds()

        # Create the full graph of the molecule (including rings)
        full_graph = nx.Graph()
        full_graph.add_nodes_from(range(len(self.molecule.atoms)))
        for bond in self.molecule.bonds:
            full_graph.add_edge(bond.atom2.idx, bond.atom1.idx)

        # Create a list of bonds which are in a rigid cycle
        rigid_cycles_bonds = set()

        for cycle in nx.cycle_basis(full_graph):
            edges = list(zip(cycle, cycle[1:] + cycle[:1]))
            if all(is_rigid_bond(self.molecule[a], self.molecule[b]) for a, b in edges):
                rigid_cycles_bonds.update(edges)

        # # Build Z-matrix
        # self.z_matrix = self.build_z_matrix()

        # # # Get root
        # # self.root_prmtop_index = self.z_matrix[0].i.index

        # # Get atoms as visited in the Z-matrix and map their prmtop indices to compound atom indices
        # self.nodes: list[int] = []
        # self.prmtop_to_compound_atom_index_map: dict[int, int] = {}
        
        # for i, row in enumerate(self.z_matrix):
        #     self.nodes.append(row.i.idx)
        #     self.prmtop_to_compound_atom_index_map[row.i.idx] = i

        # # Get bonds as (parent_prmtop_index, child_prmtop_index, is_ring_closing)
        # self.all_bonds: list[Tuple[int, int, bool]] = []

        # for row in self.z_matrix:
        #     if not row.j:
        #         continue

        #     parent_atom = row.j
        #     child_atom = row.i

        #     parent_atom_prmtop_index = parent_atom.idx
        #     child_atom_prmtop_index = child_atom.idx

        #     parent_atom_unique_name = parent_atom.residue.name + str(parent_atom.residue.idx+1) + '_' + parent_atom.name + '_' + str(parent_atom.idx+1)
        #     child_atom_unique_name = child_atom.residue.name + str(child_atom.residue.idx+1) + '_' + child_atom.name + '_' + str(child_atom.idx+1)

        #     # Make sure we did not add ring closing bonds in the Z-matrix traversal
        #     if (parent_atom_prmtop_index, child_atom_prmtop_index) in self.ring_closing_bonds:
        #         raise RuntimeError(f"Internal Error: Ring closing bond ({parent_atom_unique_name}, {child_atom_unique_name}) found in Z-matrix traversal.")
        #     if (child_atom_prmtop_index, parent_atom_prmtop_index) in self.ring_closing_bonds:
        #         raise RuntimeError(f"Internal Error: Ring closing bond ({child_atom_unique_name}, {parent_atom_unique_name}) found in Z-matrix traversal.")

        #     self.all_bonds.append((parent_atom_prmtop_index, child_atom_prmtop_index, False))

        # # Add ring closing bonds now
        # for (rcb_parent_prmtop_index, rcb_child_prmtop_index) in self.ring_closing_bonds:
        #     self.all_bonds.append((rcb_parent_prmtop_index, rcb_child_prmtop_index, True))

        # # Find non-redundant bonds (central bonds in a well defined Z-matrix dihedral)
        # self.non_redundant_bonds : set[tuple[int, int]] = set()
        # for row in self.z_matrix:
        #     if not row.j or not row.k or not row.l:
        #         continue
        #     if is_rigid_bond(row.j, row.k):
        #         continue

        #     child_atom_prmtop_index = row.j.idx
        #     parent_atom_prmtop_index = row.k.idx
        #     self.non_redundant_bonds.add((parent_atom_prmtop_index, child_atom_prmtop_index))

        # Root is the heaviest terminal atom with the smallest prmtop index
        terminal_atoms = [a for a in self.molecule.atoms if len(a.bond_partners) == 1]
        terminal_atoms = self._sort_atoms_by_mass(terminal_atoms)
        root = terminal_atoms[0].idx

        # nx.bfs_edges yields only edges that lead to previously unvisited nodes, in the order they are discovered during BFS
        # Back-edges or cross-edges are not yielded (these are ring closing bonds that we removed earlier)
        # Nodes in the edges are expressed as local prmtop indices (same as in the original prmtop file, but starting from 0 for each molecule)
        acyclic_graph_edges = nx.bfs_edges(acyclic_graph, source=root)
        acyclic_graph_edges = list(acyclic_graph_edges)

        # We are implicitly constructing a spanning tree used to define internal coordinates
        # Although networkx.bfs_edges guarantees BFS traversal, we still need to validate the parent-child relationships
        GraphTraversalUtils.validate_bfs_parent_child_edges(acyclic_graph, root, acyclic_graph_edges)

        # Create a mapping between the original indices and the BFS-explored indices
        # Again, Molmodel adds atoms via Molmodel via bondAtom(idx1, idx2)
        # Essentially, we create a mapping between the order in which atoms were added to  and their original indices
        acyclic_graph_nodes = [root] + [child_local_index for parent_local_index, child_local_index in acyclic_graph_edges]
        self.local_to_compound_atom_index_map = {}
        for compound_atom_index, prmtop_index in enumerate(acyclic_graph_nodes):
            self.local_to_compound_atom_index_map[prmtop_index] = compound_atom_index

        # Make sure all nodes have been visited via BFS
        if len(acyclic_graph_nodes) != len(self.molecule.atoms):
            missing_nodes = set(range(len(self.molecule.atoms))) - set(acyclic_graph_nodes)
            missing_atom_names = [self.molecule[idx].name for idx in missing_nodes]
            raise RuntimeError(f"Internal Error: Not all atoms were visited during BFS traversal of the molecular graph. Missing atoms: {missing_atom_names}")
        
        # Get bonds as explored by BFS starting from the root atom
        self.all_bonds = [(parent, child, False) for parent, child in acyclic_graph_edges]

        # BFS doesn't yield edges if they extend back to an already explored node - a ring closing bond
        # At this point we don't want to add ring closing bonds, so we need to filter for them
        # Molmodel adds atoms by calling bondAtom(idx1,idx2) which expects idx1 to be already bonded via bondAtom(idx0,idx1)
        # Ring closing bonds add bonds to non-bonded atoms - bondAtom(idx2,idx3) - which yields errors
        self.all_bonds += [(bond[0], bond[1], True) for bond in self.ring_closing_bonds]

        # Make sure all bonds are present
        if len(self.all_bonds) != len(self.molecule.bonds):
            raise RuntimeError("Internal Error: Number of bonds in the molecule does not match the number of bonds found during BFS traversal and ring-closing bond addition.")

        # for parent_prmtop_index, child_prmtop_index, ring_closing in self.all_bonds:
        #     parent_atom = self.molecule[parent_prmtop_index]
        #     child_atom = self.molecule[child_prmtop_index]
        #     parent_atom_unique_name = parent_atom.residue.name + str(parent_atom.residue.idx+1) + '_' + parent_atom.name + '_' + str(parent_atom.idx+1)
        #     child_atom_unique_name = child_atom.residue.name + str(child_atom.residue.idx+1) + '_' + child_atom.name + '_' + str(child_atom.idx+1)
        #     print(f"Bond between {parent_atom_unique_name} (type: {parent_atom.type}, name: {parent_atom.name}) and {child_atom_unique_name} (type: {child_atom.type}, name: {child_atom.name}). Is ring closing: {ring_closing}")

        # Find non-redundant bonds (central bonds in a well defined Z-matrix dihedral)
        # A non-redundant bond is a bond that can and should be simulated
        self.non_redundant_bonds : set[tuple[int, int]] = set()
        for child_atom_prmtop_index, parent_atom_prmtop_index in acyclic_graph_edges:
            child_atom = self.molecule[child_atom_prmtop_index]
            parent_atom = self.molecule[parent_atom_prmtop_index]

            # Terminal atoms are fixed
            if len(child_atom.bond_partners) == 1 or len(parent_atom.bond_partners) == 1:
                continue

            # We also exclude bonds involved in rigid cycles
            if (parent_atom_prmtop_index, child_atom_prmtop_index) in rigid_cycles_bonds or (child_atom_prmtop_index, parent_atom_prmtop_index) in rigid_cycles_bonds:
                continue

            # parent_atom_unique_name = parent_atom.residue.name + str(parent_atom.residue.idx+1) + '_' + parent_atom.name + '_' + str(parent_atom.idx+1)
            # child_atom_unique_name = child_atom.residue.name + str(child_atom.residue.idx+1) + '_' + child_atom.name + '_' + str(child_atom.idx+1)
            # print(f"NRB between {parent_atom_unique_name} (type: {parent_atom.type}, name: {parent_atom.name}) and {child_atom_unique_name} (type: {child_atom.type}, name: {child_atom.name}). All rigid: {all_rigid}, Overconstrained: {overconstrained}.")
            # print(f"\tRigid bonds involved for parent atom: {rigid_bonds_involved.get(parent_atom_prmtop_index, 0)}, child atom: {rigid_bonds_involved.get(child_atom_prmtop_index, 0)}")

            # Check if this bond is connected to a ring closing bond
            overconstrained = any(parent_atom_prmtop_index in t for t in self.ring_closing_bonds)
            if not overconstrained:
                overconstrained = any(child_atom_prmtop_index in t for t in self.ring_closing_bonds)

            # We can safely ignore rigid bonds that are not connected to any ring closing bond
            # Otherwise, for an atom, we risk putting two rigid bonds
            # First is from the ring closing bond (Molmodel sets it to rigid)
            # Second is from the actual bond rigidity
            if is_rigid_bond(parent_atom, child_atom) and not overconstrained:
                continue

            self.non_redundant_bonds.add((parent_atom_prmtop_index, child_atom_prmtop_index))
        # print(f"Total non-redundant bonds found: {len(self.non_redundant_bonds)}")

        # Find backbone dihedrals from the standardized list
        # Residue index -> {"phi": (a2.idx, a3.idx), "psi": (a2.idx, a3.idx)} - local prmtop indices of the central bond defining the dihedral
        from collections import defaultdict
        backbone_by_residue = defaultdict(dict)

        for dih in self.molecule.dihedrals:
            atom_types = (dih.atom1.name, dih.atom2.name, dih.atom3.name, dih.atom4.name)
            dihedral_type = (
                self.standardized_dihedrals.get(atom_types)
                or self.standardized_dihedrals.get(atom_types[::-1])
            )

            if dihedral_type not in {"phi", "psi"}:
                continue

            # Central bond defines residue ownership
            residue_idx = min(dih.atom2.residue.idx, dih.atom3.residue.idx)
            backbone_by_residue[residue_idx][dihedral_type] = (dih.atom2.idx, dih.atom3.idx)

        # Final structure: residue_idx -> ((phi2, phi3), (psi2, psi3))
        self.backbone_dihedral_bonds = []
        for res_idx, d in backbone_by_residue.items():
            if "phi" in d and "psi" in d:
                self.backbone_dihedral_bonds.append([d["phi"], d["psi"]])

        # Build atom parameters
        self.atom_params: list[AtomParams] = []
        self.num_residues = 0

        for local_index in acyclic_graph_nodes:
            a = self.molecule.atoms[local_index]
            self.atom_params.append(AtomParams(
                local_index=local_index,
                compound_atom_index=self.local_to_compound_atom_index_map[local_index],
                element_name=a.element_name,
                element_symbol=a.element_name, # TODO
                atomic_number=a.atomic_number,
                charge_e=a.ucharge.value_in_unit(pmd.unit.elementary_charge),
                mass_daltons=a.umass.value_in_unit(pmd.unit.dalton),
                vdw_radius_nm=a.urmin.value_in_unit(pmd.unit.nanometer) * 2,
                vdw_well_depth_kj=a.uepsilon.value_in_unit(pmd.unit.kilojoule_per_mole),
                sigma_nm=a.usigma.value_in_unit(pmd.unit.nanometer),
                solvent_radius_nm=a.usolvent_radius.value_in_unit(pmd.unit.nanometer),
                screen=a.screen,
                neighbors_local_indices=[node.idx for node in a.bond_partners],
                root=(local_index == root),
            ))
            self.num_residues = max(self.num_residues, a.residue.idx + 1)

        # Now we create bond parameters
        # We must canonicalize the bond indices to make sure that (i, j) and (j, i) map to the same parameter set
        # This is because we reorder bond direction (parent, child) during BFS traversal
        # We begin by creating a lookup table
        bond_lookup = {frozenset((b.atom1.idx, b.atom2.idx)): b for b in self.molecule.bonds}

        # Build bond parameters
        self.bond_params: List[BondParams] = []
        for parent_local_index, child_local_index, is_ring_closing in self.all_bonds:
            bond = bond_lookup.get(frozenset((parent_local_index, child_local_index)))
            
            if bond is None:
                raise RuntimeError(f"Internal Error: Bond between atoms {parent_local_index} and {child_local_index} not found in bond lookup.")
            
            self.bond_params.append(BondParams(
                local_indices=(
                    parent_local_index,
                    child_local_index
                ),
                compound_atom_indices=(
                    self.local_to_compound_atom_index_map[parent_local_index],
                    self.local_to_compound_atom_index_map[child_local_index]
                ),
                stiffness_in_kj_per_nm_sq=bond.type.uk.value_in_unit(pmd.unit.kilojoule_per_mole / pmd.unit.nanometer**2),
                nominal_length_in_nm=bond.type.ureq.value_in_unit(pmd.unit.nanometer),
                is_ring_closing=is_ring_closing
            ))

        # Build angle parameters
        self.angle_params: List[AngleParams] = []
        for angle in self.molecule.angles:
            self.angle_params.append(AngleParams(
                local_indices=(
                    angle.atom1.idx,
                    angle.atom2.idx,
                    angle.atom3.idx
                ),
                compound_atom_indices=(
                    self.local_to_compound_atom_index_map[angle.atom1.idx],
                    self.local_to_compound_atom_index_map[angle.atom2.idx],
                    self.local_to_compound_atom_index_map[angle.atom3.idx]
                ),
                stiffness_in_kj_per_rad_sq=angle.type.uk.value_in_unit(pmd.unit.kilojoule_per_mole / pmd.unit.radian**2),
                nominal_angle_in_deg=angle.type.utheteq.value_in_unit(pmd.unit.degree)
            ))


        # Build AMBER proper and improper periodic torsions
        # They have multiple terms, each with their own amplitude, phase and periodicity
        periodic_torsion_terms: dict[tuple[int, int, int, int], list[rb.RoboPeriodicTorsionTerm]] = {}
        for d in self.molecule.dihedrals:

            # TODO cannot be 0 periodicity?
            # TODO we also have amplitude 0 terms that we should ignore perhaps?
            # TODO from DuMM defineBondTorsion: Pay particular attention to the amplitude -- in our convention it is really a half-amplitude since the sinusoids range from -1 to 1 making the full energy and torque "excursion" twice the amplitude.
            # TODO check the other functions as well, they have some cautions about amplitude signs etc
            if d.type.per == 0:
                raise ValueError(f"Torsion term with 0 periodicity found for torsion between atoms {key[0]}, {key[1]}, {key[2]}, {key[3]}")
            
            # Phase is sensitive to floating point errors
            EPS = 1e-3
            phase_deg=d.type.uphase.value_in_unit(pmd.unit.degree)
            phase_deg = ((phase_deg + 180) % 360) - 180
            if phase_deg > 180.0 + EPS:
                raise ValueError(f"Torsion phase {phase_deg:.6f}° exceeds 180 beyond tolerance ({EPS}°)")
            else:
                phase_deg = 180.0 if phase_deg > 180.0 else phase_deg

            key = d.atom1.idx, d.atom2.idx, d.atom3.idx, d.atom4.idx
            term = rb.RoboPeriodicTorsionTerm(
                amplitude_kj=d.type.uphi_k.value_in_unit(pmd.unit.kilojoule_per_mole),
                phase_deg=phase_deg,
                periodicity=d.type.per
            )

            if key not in periodic_torsion_terms:
                periodic_torsion_terms[key] = { "improper": d.improper, "terms": [term] }
            else:
                # Check that we don't mix proper and improper torsions for the same set of atoms
                if periodic_torsion_terms[key]["improper"] != d.improper:
                    raise ValueError(f"Mixed proper/improper torsion for key {key}")
                
                # Check that all periodicities are unique for this torsion
                existing_periodicities = {t.periodicity for t in periodic_torsion_terms[key]["terms"]}
                if term.periodicity in existing_periodicities:
                    raise ValueError(f"Duplicate periodicity {term.periodicity} found for torsion between atoms {key[0]}, {key[1]}, {key[2]}, {key[3]}")

                periodic_torsion_terms[key]["terms"].append(term)

        # Now that we have collected all periodic torsion terms, build the parameter objects
        self.periodic_torsion_params: List[PeriodicTorsionParams] = []
        for key, value in periodic_torsion_terms.items():
            atom1_local_index, atom2_local_index, atom3_local_index, atom4_local_index = key

            self.periodic_torsion_params.append(PeriodicTorsionParams(
                local_indices=(
                    atom1_local_index,
                    atom2_local_index,
                    atom3_local_index,
                    atom4_local_index
                ),
                compound_atom_indices=(
                    self.local_to_compound_atom_index_map[atom1_local_index],
                    self.local_to_compound_atom_index_map[atom2_local_index],
                    self.local_to_compound_atom_index_map[atom3_local_index],
                    self.local_to_compound_atom_index_map[atom4_local_index]
                ),
                is_improper=value["improper"],
                terms=value["terms"]
            ))

        # Build CHARMM improper harmonic torsions
        # They are absent in AMBER prmtop files
        self.improper_harmonic_torsion_terms: List[HarmonicImproperTorsionParams] = []
        for imp in self.molecule.impropers:
            self.improper_harmonic_torsion_terms.append(HarmonicImproperTorsionParams(
                local_indices=(
                    imp.atom1.idx,
                    imp.atom2.idx,
                    imp.atom3.idx,
                    imp.atom4.idx
                ),
                compound_atom_indices=(
                    self.local_to_compound_atom_index_map[imp.atom1.idx],
                    self.local_to_compound_atom_index_map[imp.atom2.idx],
                    self.local_to_compound_atom_index_map[imp.atom3.idx],
                    self.local_to_compound_atom_index_map[imp.atom4.idx]
                ),
                stiffness_in_kj_per_rad_sq=imp.type.upsi_k.value_in_unit(pmd.unit.kilojoule_per_mole / pmd.unit.radian**2),
                nominal_angle_in_rad=imp.type.upsi_eq.value_in_unit(pmd.unit.radian)
            ))
        
    def _find_forbidden_bonds(self, atom_type_pairs: Iterable[Tuple[str, str]]) -> Set[FrozenSet[int]]:
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
            frozenset((a.lower(), b.lower()))
            for a, b in atom_type_pairs
        }

        forbidden_edges: Set[FrozenSet[int]] = set()

        for bond in self.molecule.bonds:
            a1, a2 = bond.atom1, bond.atom2
            atom_pair = frozenset((a1.name.lower(), a2.name.lower()))

            if atom_pair in forbidden_pairs:
                forbidden_edges.add(frozenset((a1.idx, a2.idx)))

        return forbidden_edges

    @staticmethod
    def _bond_edges_to_source_nodes(edges: Iterable[FrozenSet[int]]) -> Set[int]:
        """
        Convert a collection of bond edges into a set of source nodes.

        Parameters
        ----------
        edges : iterable of frozenset(int, int)
            Bond edges.

        Returns
        -------
        nodes : set of int
            Union of all atom indices participating in the bonds.
        """
        nodes: Set[int] = set()
        for edge in edges:
            nodes.update(edge)
        return nodes
    
    @staticmethod
    def _distances_to_nodes(graph: nx.Graph, source_nodes: Iterable[int]) -> Dict[int, Optional[int]]:
        """
        Compute the shortest path distance (in number of bonds) from every node
        in a graph to the nearest source node.

        Parameters
        ----------
        source_nodes : iterable of int
            Nodes that act as distance-0 sources.

        Returns
        -------
        distances : dict[int, int or None]
            Mapping from node index to shortest distance.
            Nodes unreachable from any source have value None.
        """
        dist = {node: None for node in graph.nodes}
        q = deque()

        for node in source_nodes:
            if node in graph:
                dist[node] = 0
                q.append(node)

        while q:
            u = q.popleft()
            for v in graph.neighbors(u):
                if dist[v] is None:
                    dist[v] = dist[u] + 1
                    q.append(v)

        return dist

    def _find_ring_closing_bonds(self) -> tuple[set[tuple[int, int]], nx.Graph]:
        # Build the molecular graph
        # It may contain cycles (proline, aromatic rings, macrocycles between disuflide bonds, etc)
		# The algorithm that computes the Z-matrix requires a loop-free structure
		# Therefore, we start by removing standardized ring-closing dihedrals from the graph (proline, tryptophan, etc)
        acyclic_graph = nx.Graph()

        # Atoms in non-redundant bonds must be involved in at most 1 rigid bond
        # Otherwise, the system becomes overconstrained
        rigid_bonds_involved = {}
        ring_closing_bonds = set()

        # Iterate over all bonds in the molecule
        for i, bond in enumerate(self.molecule.bonds):
            # Get local prmtop indices of parent and child atoms
            # A local prmtop index is the index of the atom in the molecule structure
            # For the 1st molecule, it is equal to the prmtop index in the original file
            # For other molecules, it is offset by the number of atoms in previous molecules
            parent_atom = bond.atom2
            parent_atom_prmtop_index = parent_atom.idx

            child_atom = bond.atom1
            child_atom_prmtop_index = child_atom.idx

            # Build the possible dihedrals for this bond
            # This means getting all parents of parent atom and all children of child atom
            # We are only interested in the atom types
            candidate_dihedral_atom_types = self.get_possible_dihedrals_atom_types(parent_atom, child_atom)

            # If any of the candidate dihedrals is a standardized ring-closing dihedral, skip adding this bond to the graph
            if self.is_standardized_ring_closing_dihedral(candidate_dihedral_atom_types):
                if parent_atom_prmtop_index in rigid_bonds_involved and rigid_bonds_involved[parent_atom_prmtop_index] >= 1:
                    raise ValueError(f"Atom with prmtop index {parent_atom_prmtop_index} is involved in multiple rigid bonds.")
                if child_atom_prmtop_index in rigid_bonds_involved and rigid_bonds_involved[child_atom_prmtop_index] >= 1:
                    raise ValueError(f"Atom with prmtop index {child_atom_prmtop_index} is involved in multiple rigid bonds.")
                
                # Store the ring closing bond and remove it from the graph
                rigid_bonds_involved[parent_atom_prmtop_index] = rigid_bonds_involved.get(parent_atom_prmtop_index, 0) + 1
                rigid_bonds_involved[child_atom_prmtop_index] = rigid_bonds_involved.get(child_atom_prmtop_index, 0) + 1

                ring_closing_bonds.add((parent_atom_prmtop_index, child_atom_prmtop_index))
            else:
                acyclic_graph.add_edge(parent_atom_prmtop_index, child_atom_prmtop_index)

        # Find forbidden bonds
        forbidden_edges = self._find_forbidden_bonds(atom_type_pairs=[("SG", "SG")])
        forbidden_nodes = self._bond_edges_to_source_nodes(forbidden_edges)

        if forbidden_nodes:
            forbidden_atom_dist = self._distances_to_nodes(acyclic_graph, forbidden_nodes)
        else:
            forbidden_atom_dist = {}

        # Check if we have any residual ring closing bonds
        # At this point, they can be disulfide bonds or rings not present in our database
        # A basis for cycles of a network is a minimal collection of cycles such that any cycle in the network can be written as a sum of cycles in the basis
        # This means that a bond may be part of multiple cycles

        # Travese the shortest cycles first
        candidate_cycles = nx.cycle_basis(acyclic_graph)
        for candidate_cycle in sorted(candidate_cycles, key=len):
            n = len(candidate_cycle)
            candidate_bonds: list[BondProperties] = []

            for i in range(n):
                u = candidate_cycle[i]
                v = candidate_cycle[(i + 1) % n]

                atom_u = self.molecule[u]
                atom_v = self.molecule[v]

                # Skip forbidden bonds entirely
                if frozenset((u, v)) in forbidden_edges:
                    continue

                degree_sum = acyclic_graph.degree[u] + acyclic_graph.degree[v]
                is_rigid = is_rigid_bond(atom_u, atom_v)

                if forbidden_atom_dist:
                    # Distance to a bond = max distance of its two endpoint atoms
                    forbidden_bond_distance = max(
                        forbidden_atom_dist[u],
                        forbidden_atom_dist[v]
                    )
                else:
                    forbidden_bond_distance = None

                candidate_bonds.append(BondProperties(
                    parent_atom_prmtop_index=u,
                    child_atom_prmtop_index=v,
                    is_rigid=is_rigid,
                    sg_sg_bond_distance=forbidden_bond_distance,
                    sum_of_degrees=degree_sum
                ))

            # Rigid first, then bonds farthest from forbidden bonds, then lowest degree first
            candidate_bonds.sort(key=lambda bond: (
                not bond.is_rigid,
                -bond.sg_sg_bond_distance if bond.sg_sg_bond_distance is not None else -1,
                bond.sum_of_degrees
            ))

            for bond in candidate_bonds:
                u, v = bond.parent_atom_prmtop_index, bond.child_atom_prmtop_index
                if (u, v) in ring_closing_bonds or (v, u) in ring_closing_bonds:
                    continue
                selected_bond = bond
                break

            # One atom cannot be involved in multiple rigid bonds because that would create an overconstrained system
            u, v = selected_bond.parent_atom_prmtop_index, selected_bond.child_atom_prmtop_index
            if u in rigid_bonds_involved and rigid_bonds_involved[u] >= 1:
                raise ValueError(f"Atom with prmtop index {u} is involved in multiple rigid bonds.")
            if v in rigid_bonds_involved and rigid_bonds_involved[v] >= 1:
                raise ValueError(f"Atom with prmtop index {v} is involved in multiple rigid bonds.")
            
            # Store the ring closing bond and remove it from the graph
            rigid_bonds_involved[u] = rigid_bonds_involved.get(u, 0) + 1
            rigid_bonds_involved[v] = rigid_bonds_involved.get(v, 0) + 1
            
            acyclic_graph.remove_edge(u, v)
            ring_closing_bonds.add((u, v))

        if acyclic_graph.number_of_nodes() != len(self.molecule.atoms):
            raise ValueError("Molecular graph contains nodes that are not in the molecule's atom list.")
        if not nx.is_forest(acyclic_graph):
            raise ValueError("Molecular graph contains cycles: not enough ring-closing bonds were removed.")
        if not nx.is_connected(acyclic_graph):
            raise ValueError("Molecular graph is disconnected: too many ring-closing bonds were removed.")
        
        return ring_closing_bonds, acyclic_graph
    
    def get_ring_closing_bonds(self) -> List[Tuple[int, int]]:
        return list(self.ring_closing_bonds)

    def get_possible_dihedrals_atom_types(self, parent_atom: pmd.Atom, child_atom: pmd.Atom) -> List[Tuple[str, str, str, str]]:
        dihedral_candidates = []
        for grandparent_atom in parent_atom.bond_partners:
            if grandparent_atom == child_atom:
                continue
            for nephew_atom in child_atom.bond_partners:
                if nephew_atom == parent_atom:
                    continue
                dihedral_candidates.append((grandparent_atom.name, parent_atom.name, child_atom.name, nephew_atom.name))
        return dihedral_candidates

    def is_standardized_ring_closing_dihedral(self, candidate_dihedral_atom_types: List[Tuple[str, str, str, str]]) -> bool:
        return any(min(t, t[::-1]) in self.standardized_ring_closing_dihedral_ref_set for t in candidate_dihedral_atom_types)
    
    def is_ring_closing(self, a: pmd.Atom, b: pmd.Atom) -> bool:
        i, j = a.idx, b.idx
        return (i, j) in self.ring_closing_bonds or (j, i) in self.ring_closing_bonds
    
    def local_to_compound_atom_index(self, local_index: int) -> int:
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
    
    def get_nodes_as_prmtop_indices(self) -> list[int]:
        return list(self.nodes)
    
    def get_nonredundant_bonds(self) -> list[tuple[int, int]]:
        return list(self.non_redundant_bonds)

    def get_bonds_as_local_indices(self) -> list[tuple[int, int, bool]]:
        """
        Get local indices of bonds in the orderd they were explored.
        Note that local indices are different from compound atom indices.

        Returns
        -------
        list[tuple[int, int, bool]]
            A list of bonds as tuples of (parent_local_index, child_local_index, is_ring_closing)
        """
        return self.all_bonds
    
    def get_angles_as_local_indices(self) -> list[tuple[int, int, int]]:
        """
        Get local indices of angles in the molecule.
        Indices are not canonicalized ie there is no guarantee that local_index_1 < local_index_3.
        Note that local indices are different from compound atom indices.

        Returns
        -------
        list[tuple[int, int, int]]
            A list of angles as tuples of (local_index_1, local_index_2, local_index_3)
        """
        return self.angles

    def get_periodic_dihedrals_as_local_indices(self) -> list[tuple[int, int, int, int]]:
        """
        Get local indices of periodic dihedrals in the molecule.
        Note that local indices are different from compound atom indices.

        Returns
        -------
        list[tuple[int, int, int, int]]
            A list of periodic dihedrals as tuples of (local_index_1, local_index_2, local_index_3, local_index_4)
        """
        return list(self.periodic_dihedrals)
    
    def get_harmonic_impropers_as_local_indices(self) -> list[tuple[int, int, int, int]]:
        """
        Get local indices of harmonic impropers in the molecule.
        Note that local indices are different from compound atom indices.

        Returns
        -------
        list[tuple[int, int, int, int]]
            A list of harmonic impropers as tuples of (local_index_1, local_index_2, local_index_3, local_index_4)
        """
        return self.harmonic_impropers

    
    def get_standardized_dihedral_type(self, parent_atom: pmd.Atom, child_atom: pmd.Atom) -> str | None:
        # Build a list of dihedral candidates by looking at bonded atoms to parent and child
        dihedral_type = None
        for grandparent_atom in parent_atom.bond_partners:
            for nephew_atom in child_atom.bond_partners:
                if grandparent_atom == child_atom:
                    continue
                if nephew_atom == parent_atom:
                    continue

                # Try to match into our database of standardized dihedrals
                dihedral_atom_types = (grandparent_atom.name, parent_atom.name, child_atom.name, nephew_atom.name)
                dihedral_type = self.standardized_dihedrals.get(dihedral_atom_types, None)
                if dihedral_type is None:
                    dihedral_type = self.standardized_dihedrals.get(tuple(reversed(dihedral_atom_types)), None)
                if dihedral_type is not None:
                    break

            if dihedral_type is not None:
                break

        return dihedral_type
    
    # def get_ring_closing_bonds(self, graph: nx.Graph) -> List[Tuple[int, int]]:
    #     # Find all cycles in the graph
    #     cycles = nx.cycle_basis(graph)

    #     ring_closing_bonds = []
    #     for cycle in cycles:
    #         # For each cycle, get the edges that form the cycle
    #         cycle_edges = []
    #         for i in range(len(cycle)):
    #             j = (i + 1) % len(cycle)
    #             edge = (min(cycle[i], cycle[j]), max(cycle[i], cycle[j]))
    #             cycle_edges.append(edge)

    #         # Add the edges to the ring closing bonds list
    #         ring_closing_bonds.extend(cycle_edges)

    #     return ring_closing_bonds
    
    # def build_z_matrix(self) -> list[ZMatrixRowAtomGroups]:
    #     """
    #     Construct a Z-matrix (BAT-style internal coordinate tree) from an MDAnalysis Universe molecule.
    #     Cycles are handled by explicitly excluding ring-closing bonds, turning the molecular graph into a spanning tree.
    #     This function expects that the molecular graph is acyclic (a forest) after removing the pre-defined ring-closing bonds.
    #     Here, we try to select root atoms that are heavy and non-terminal to improve numerical stability:
    #     1. The first atom is chosen as the heaviest terminal atom.
    #     2. The second atom is the its bonded neighbor.
    #     3. The third atom is the heaviest non-terminal bonded neighbor of the second atom that is not part of a ring-closing bond and is not collinear with the first two atoms.

    #     The result is torsional spanning tree while explicitly removing redundant ring closing bonds.

    #     Parameters
    #     ----------
    #     ag_o : list of MDAnalysis Atoms
    #         List to sort
    #     reverse : bool
    #         Atoms will be in descending order

    #     Returns
    #     -------
    #     ag_n : list of Atoms
    #         an ordered, loop-free internal coordinate definition
    #     """

    #     # The molecular graph must be acyclic for the Z-matrix to be built correctly
    #     # This is ensured by removing ring-closing bonds beforehand (see above)
    #     # Basically, we check here that our ring closing bond detection worked correctly and we did not miss any cycles
    #     # A forest is a graph with no undirected cycles
    #     if not nx.is_forest(self.graph):
    #         raise ValueError("Molecule contains cycles; cannot build Z-matrix.")
        
    #     # Build the Z matrix
    #     z_matrix: list[ZMatrixRowAtomGroups] = []
        
    #     # We begin by excluding degenerate cases with 1, 2, or 3 atoms
    #     if len(self.molecule.atoms) == 0:
    #         raise ValueError("Molecule contains no atoms; cannot build Z-matrix.")
    #     if len(self.molecule.atoms) == 1:
    #         initial_atom = self.molecule.atoms[0]
    #         z_matrix.append(ZMatrixRowAtomGroups(i=initial_atom, j=None, k=None, l=None))
    #         return z_matrix
    #     if len(self.molecule.atoms) == 2:
    #         initial_atom = self.molecule.atoms[0]
    #         second_atom = self.molecule.atoms[1]
    #         z_matrix.append(ZMatrixRowAtomGroups(i=initial_atom, j=None, k=None, l=None))
    #         z_matrix.append(ZMatrixRowAtomGroups(i=second_atom, j=initial_atom, k=None, l=None))
    #         return z_matrix
        
    #     # We select the initial atom from terminal atoms (atoms with only one bond)
    #     # TODO Macrocycles with substituents may give terminals that lie on side chains, producing pathological Z-matrices.
    #     # TODO Root selection must be graph-theoretic, not terminal-based. Use: degree-1 atoms if they exist, otherwise fall back to highest-degree or highest-mass atom not in a ring closure set.
    #     terminal_atoms = [a for a in self.molecule.atoms if len(a.bonds) == 1]
    #     terminal_atoms = self._sort_atoms_by_mass(terminal_atoms)

    #     # Select the heaviest root atom from the heaviest terminal atoms
    #     initial_atom = terminal_atoms[0]
    #     z_matrix.append(ZMatrixRowAtomGroups(i=initial_atom, j=None, k=None, l=None))

    #     # The next atom in the root is bonded to the initial atom
    #     # Since the initial atom is a terminal atom, there is only one bonded atom
    #     second_atom = initial_atom.bond_partners[0]
    #     z_matrix.append(ZMatrixRowAtomGroups(i=second_atom, j=initial_atom, k=None, l=None))

    #     # The last atom in the root is the heaviest atom bonded to the second atom
    #     # If there are more than three atoms, then the last atom cannot be a terminal atom
    #     # Moreover, the three atoms must not be collinear
    #     if len(self.molecule.atoms) != 3:
    #         third_atom_candidates = []
    #         for a in self._tree_neighbors(second_atom):
    #             if (a != initial_atom) and (a not in terminal_atoms): # and not self._are_collinear(initial_atom.position, second_atom.position, a.position)
    #                 third_atom_candidates.append(a)
    #         third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]
    #     else:
    #         third_atom_candidates = []
    #         for a in self._tree_neighbors(second_atom):
    #             if a != initial_atom:
    #                 third_atom_candidates.append(a)
    #         third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]

    #     third_atom = self._sort_atoms_by_mass(third_atom_candidates)[0]
    #     z_matrix.append(ZMatrixRowAtomGroups(i=third_atom, j=second_atom, k=initial_atom, l=None))

    #     # Root triplet
    #     root = [initial_atom, second_atom, third_atom]

    #     for i, j, k, l in self._find_torsions(root):
    #         z_matrix.append(ZMatrixRowAtomGroups(i=i, j=j, k=k, l=l))

    #     # # Print unique atom names for each row in the Z-matrix
    #     # for row in z_matrix:
    #     #     atom_name_i = row.i.residue.name + row.i.resid.__str__() + '_' + row.i.name + '_' + row.i.index.__str__()
    #     #     atom_name_j = None if row.j is None else row.j.resname + row.j.resid.__str__() + '_' + row.j.name + '_' + row.j.index.__str__()
    #     #     atom_name_k = None if row.k is None else row.k.resname + row.k.resid.__str__() + '_' + row.k.name + '_' + row.k.index.__str__()
    #     #     atom_name_l = None if row.l is None else row.l.resname + row.l.resid.__str__() + '_' + row.l.name + '_' + row.l.index.__str__()
    #     #     sele = 'sele id ' + row.i.index.__str__() + '+' + (row.j.index.__str__() if row.j is not None else '') + '+' + (row.k.index.__str__() if row.k is not None else '') + '+' + (row.l.index.__str__() if row.l is not None else '')
    #     #     print(f"Z-matrix row: i={atom_name_i}, j={atom_name_j}, k={atom_name_k}, l={atom_name_l}", sele)

    #     return z_matrix
    
    # def _tree_neighbors(self, atom: pmd.Atom) -> list[pmd.Atom]:
    #     """
    #     Builds the list atoms bonded to the given atom, excluding those involved in ring-closing bonds.

    #     Parameters
    #     ----------
    #     atom : MDAnalysis Atom
    #         The atom for which to find tree neighbors.

    #     Returns
    #     -------
    #     neighbors : list of MDAnalysis Atom
    #         List of bonded atoms not involved in ring-closing bonds.        
    #     """
    #     return [a for a in atom.bond_partners if not self.is_ring_closing(atom, a)]
    
    # def _find_torsions(self, root: list[pmd.Atom]) -> list[tuple[pmd.Atom, pmd.Atom, pmd.Atom, pmd.Atom]]:
    #     """
    #     Constructs a list of torsion angles.

    #     Returns
    #     -------
    #     torsions : list of AtomGroup
    #         list of AtomGroup objects that define torsion angles
    #     """
    #     torsions: list[tuple[pmd.Atom, pmd.Atom, pmd.Atom, pmd.Atom]] = []
    #     selected_atoms = list(root)
    #     selected_atoms_set = set([root[0].idx, root[1].idx, root[2].idx])

    #     # Build tree neighbors
    #     tree_adj: dict[int, list[pmd.Atom]] = {
    #         # atom.idx: self._sort_atoms_by_mass([a for a in atom.bond_partners if not self.is_ring_closing(atom, a)]) for atom in self.molecule.atoms
    #         atom.idx: self._sort_atoms_by_mass([a for a in atom.bond_partners]) for atom in self.molecule.atoms
    #     }

    #     while len(selected_atoms) < len(self.molecule.atoms):
    #         torsionAdded = False
    #         for a1 in selected_atoms:
    #             # Find a0, which is a new atom connected to the selected atom
    #             a0_list = [a for a in tree_adj[a1.idx] if a.idx not in selected_atoms_set]
    #             for a0 in a0_list:
    #                 # Find a2, which is connected to a1, is not a terminal atom and has been selected
    #                 a2_list = [a for a in tree_adj[a1.idx] if (a != a0) and len(a.bond_partners) > 1 and (a.idx in selected_atoms_set)]
    #                 for a2 in a2_list:
    #                     # Find a3, which is connected to a2, has been selected, and is not a1
    #                     a3_list = [a for a in tree_adj[a2.idx] if (a != a1) and (a.idx in selected_atoms_set)]
    #                     for a3 in a3_list:
    #                         # Add the torsion to the list of torsions
    #                         torsions.append((a0, a1, a2, a3))

    #                         # Add the new atom to selected_atoms which extends the loop
    #                         selected_atoms.append(a0)
    #                         selected_atoms_set.add(a0.idx)
    #                         torsionAdded = True
    #                         break
    #                     break

    #         if torsionAdded is False:
    #             print("Selected atoms:")
    #             print([a.idx + 1 for a in selected_atoms])
    #             print("Torsions found:")
    #             print([list(t.indices + 1) for t in torsions])
    #             raise ValueError("Additional torsions not found.")

    #     return torsions
    
    @staticmethod
    def _sort_atoms_by_mass(atoms: list[pmd.Atom]) -> list[pmd.Atom]:
        r"""Sorts a list of atoms in descending order by mass.
        The atom index is used as a tiebreaker so that the ordering is reproducible.
        If two atoms have the same mass, the one with the lower index comes first.

        Parameters
        ----------
        ag_o : list of MDAnalysis Atoms
            List of atoms to sort

        Returns
        -------
        ag_n : list of Atoms
            Sorted list
        """

        # Negate mass to sort it descending, keep idx as is for ascending
        return sorted(atoms, key=lambda a: (-a.mass, a.idx))
    
    # @staticmethod
    # def _are_collinear(p1: np.ndarray, p2: np.ndarray, p3: np.ndarray, eps=_EPS, tol=_COLLINEARITY_TOL) -> bool:
    #     """Determines if three points are collinear using a scale-normalized cross-product area check.

    #     Parameters
    #     ----------
    #     p1 : np.ndarray
    #         First point.
    #     p2 : np.ndarray
    #         Second point.
    #     p3 : np.ndarray
    #         Third point.
    #     eps : float
    #         Tolerance for considering points coincident.
    #     tol : float
    #         Tolerance for considering points collinear.

    #     Returns
    #     -------
    #     bool
    #         True if the points are collinear, False otherwise.
    #     """
    #     # Create vectors relative to p1
    #     v1 = p2 - p1
    #     v2 = p3 - p1

    #     # Compute the magnitude of the cross product (Area of parallelogram)
    #     # For N-dimensions, we use the fact that ||a x b|| = sqrt(||a||^2 ||b||^2 - (a.b)^2)
    #     # This is more robust than a standard 3D cross product call
    #     norm_v1_sq = np.dot(v1, v1)
    #     norm_v2_sq = np.dot(v2, v2)
    #     dot_v1_v2 = np.dot(v1, v2)
        
    #     # Area squared (Lagrange's Identity)
    #     area_sq = norm_v1_sq * norm_v2_sq - dot_v1_v2**2
        
    #     # Handle coincident points: if any two points are the same, they are collinear
    #     if np.isclose(norm_v1_sq, 0, atol=eps) or np.isclose(norm_v2_sq, 0, atol=eps):
    #         return True

    #     # Scale-Invariant Normalization
    #     # We divide the squared area by the product of the squared lengths
    #     # This is equivalent to sin^2(theta)
    #     # Unlike cos(theta), sin(theta) is highly sensitive near 0
    #     relative_error_sq = area_sq / (norm_v1_sq * norm_v2_sq)

    #     # We use a very small tolerance here
    #     # Since we are dealing with squared values, 1e-12 corresponds to  a sine of 1e-6
    #     return abs(relative_error_sq) < tol