import MDAnalysis as mda
import MDAnalysis.analysis.bat as bat

from enum import IntEnum
from abc import ABC, abstractmethod
from typing import List, Tuple
from dataclasses import dataclass

# Hold compound atom indices, not prmtop indices
@dataclass
class ZMatrixRowAtomGroups:
    i: mda.core.groups.Atom
    j: mda.core.groups.Atom | None
    k: mda.core.groups.Atom | None
    l: mda.core.groups.Atom | None

def bond_ag(row: ZMatrixRowAtomGroups) -> mda.AtomGroup | None:
    return None if row.j is None else mda.AtomGroup([row.i, row.j])

# def angle_ag(row: ZMatrixRowAtomGroups, u: mda.Universe) -> mda.AtomGroup | None:
#     return None if row.k is None else u.atoms[[row.i, row.j, row.k]]

# def dihedral_ag(row: ZMatrixRowAtomGroups, u: mda.Universe) -> mda.AtomGroup | None:
#     return None if row.l is None else u.atoms[[row.i, row.j, row.k, row.l]]

def build_z_matrix(molecule: mda.Universe) -> list[ZMatrixRowAtomGroups]:
    # Builds the torsional spanning tree ie the subgraph connecting all vertices with the fewest possible edges, creating a loop-free structure

    # Find terminal atoms (atoms with only one bond)
    terminal_atoms = bat._sort_atoms_by_mass(
        [a for a in molecule.atoms if len(a.bonds) == 1], reverse=True
    )

    # Select the heaviest root atom from the heaviest terminal atoms
    initial_atom = terminal_atoms[0]

    # The next atom in the root is bonded to the initial atom
    # Since the initial atom is a terminal atom, there is only one bonded atom
    second_atom = initial_atom.bonded_atoms[0]

    # The last atom in the root is the heaviest atom bonded to the second atom
    # If there are more than three atoms, then the last atom cannot be a terminal atom
    if molecule.n_atoms != 3:
        third_atom = bat._sort_atoms_by_mass(
            [
                a
                for a in second_atom.bonded_atoms
                if (a in molecule)
                and (a != initial_atom)
                and (a not in terminal_atoms)
            ],
            reverse=True,
        )[0]
    else:
        third_atom = bat._sort_atoms_by_mass(
            [
                a
                for a in second_atom.bonded_atoms
                if (a in molecule) and (a != initial_atom)
            ],
            reverse=True,
        )[0]

        # # Assert that the three atoms are not collinear
        # vec1 = second_atom.position - initial_atom.position
        # vec2 = third_atom.position - second_atom.position
        # cos_angle = bat._cosine_of_angle_between_vectors(vec1, vec2)
        # assert abs(cos_angle) < 0.999, "The three root atoms are collinear."

    # Root triplet
    root = mda.AtomGroup([initial_atom, second_atom, third_atom])

    # Build the Z matrix
    z_matrix: list[ZMatrixRowAtomGroups] = []

    z_matrix.append(ZMatrixRowAtomGroups(i=initial_atom, j=None, k=None, l=None))
    z_matrix.append(ZMatrixRowAtomGroups(i=second_atom, j=initial_atom, k=None, l=None))
    z_matrix.append(ZMatrixRowAtomGroups(i=third_atom, j=second_atom, k=initial_atom, l=None))

    for i, j, k, l in bat._find_torsions(root, molecule):
        z_matrix.append(ZMatrixRowAtomGroups(i=i, j=j, k=k, l=l))

    return z_matrix

class MoleculeType(IntEnum):
    PROTEIN = 0
    NUCLEIC_ACID = 1
    LIGAND = 2
    WATER = 3
    ION = 4

def get_molecule_type(molecule: mda.Universe) -> MoleculeType:
    return MoleculeType.PROTEIN

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
