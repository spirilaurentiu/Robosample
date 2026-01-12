import MDAnalysis as mda
import networkx as nx
from typing import Tuple, List, Set

import molecule

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
        'chi5': ('CD', 'NE', 'CZ', 'NH1'),
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
        'ring_closing_1': ('CB', 'CG', 'ND1', 'CE1'),
    },
    'HIE': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND1'),
        'ring_closing_1': ('CB', 'CG', 'ND1', 'CE1'),
    },
    'HID': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'ND1'),
        'ring_closing_1': ('CB', 'CG', 'ND1', 'CE1'),
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
        'ring_closing_1': ('CB', 'CG', 'CD1', 'CE1'), # Same as TYR
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
        'ring_closing_1': ('CB', 'CG', 'CD1', 'NE1'), # Pyrrole ring
        'ring_closing_2': ('CG', 'CD2', 'CE3', 'CZ3'), # Benzene ring
    },

    # Tyrosine
    'TYR': {
        'chi1': ('N', 'CA', 'CB', 'CG'),
        'chi2': ('CA', 'CB', 'CG', 'CD1'),
        'chi3': ('CE1', 'CZ', 'OH', 'HH'),
        'ring_closing_1': ('CB', 'CG', 'CD1', 'CE1'), # Same as PHE
    },

    # Valine
    'VAL': {
        'chi1': ('N', 'CA', 'CB', 'CG1'),
        'chi2.1': ('CA', 'CB', 'CG1', 'HG11'),
        'chi2.2': ('CA', 'CB', 'CG2', 'HG21'),
    }
}

class Protein(molecule.Molecule):
    def __init__(self,
                 universe: mda.Universe,
                 parm,
                 include_omega: bool = False,
                 include_chi1: bool = True,
                 include_chi2: bool = True,
                 include_chi3: bool = True,
                 include_chi4: bool = True,
                 include_chi5: bool = True,
                 want_n_terminus_phi_rigid: bool = True,
                 want_c_terminus_psi_rigid: bool = True):
        
        super().__init__()
        self.universe = universe
        
        # Create the dihedral selections based on input flags
        self.dihedral_sele = PROTEIN_DIHEDRAL_SELECTION
        for res in self.dihedral_sele:
            if 'chi1' in self.dihedral_sele[res] and not include_chi1:
                del self.dihedral_sele[res]['chi1']

            if 'chi2' in self.dihedral_sele[res] and not include_chi2:
                del self.dihedral_sele[res]['chi2']
            if 'chi2.1' in self.dihedral_sele[res] and not include_chi2:
                del self.dihedral_sele[res]['chi2.1']
            if 'chi2.2' in self.dihedral_sele[res] and not include_chi2:
                del self.dihedral_sele[res]['chi2.2']

            if 'chi3' in self.dihedral_sele[res] and not include_chi3:
                del self.dihedral_sele[res]['chi3']
            if 'chi3.1' in self.dihedral_sele[res] and not include_chi3:
                del self.dihedral_sele[res]['chi3.1']
            if 'chi3.2' in self.dihedral_sele[res] and not include_chi3:
                del self.dihedral_sele[res]['chi3.2']
                
            if 'chi4' in self.dihedral_sele[res] and not include_chi4:
                del self.dihedral_sele[res]['chi4']

            if 'chi5' in self.dihedral_sele[res] and not include_chi5:
                del self.dihedral_sele[res]['chi5']

        # Phi in undefined for the first residue
        self.phi_dihedrals = [phi for res in self.universe.residues if (phi := res.phi_selection())]

        # Psi is undefined for the last residue
        self.psi_dihedrals = [psi for res in self.universe.residues if (psi := res.psi_selection())]

        # Omega dihedrals are defined for all residues
        self.omega_dihedrals = [omega for res in self.universe.residues if (omega := res.omega_selection())]

        # Build chi dihedrals from user-defined selections
        self.chi_dihedrals = []
        for res in self.universe.residues:
            for chi_name, chi_atoms in self.dihedral_sele[res.resname].items():
                chi_atom_ix = [self.universe.select_atoms(f"resid {res.resid} and name {atom_name}").ix[0] for atom_name in chi_atoms]
                chi_dihedral = mda.AtomGroup(chi_atom_ix, self.universe)

                assert len(chi_dihedral) == 4, f"Dihedral {chi_name} in residue {res.resname} {res.resid} does not have 4 atoms: {chi_atoms}, found {chi_dihedral}"
                self.chi_dihedrals.append((chi_dihedral[0].index, chi_dihedral[1].index, chi_dihedral[2].index, chi_dihedral[3].index))

        # Find prmtop index of root atom
        root_residue = self.universe.residues[0]
        # root_index_prmtop = self.universe.select_atoms(f"resid {root_residue.resid} and name N")[0].index
        root_index_prmtop = self.universe.atoms[0].index

        dihh = molecule.find_dihedrals(self.universe)
        for d in dihh:
            atom_0 = d[0].index
            atom_1 = d[1].index
            atom_2 = d[2].index
            atom_3 = d[3].index

            a0 = parm.atoms[atom_0]
            a1 = parm.atoms[atom_1]
            a2 = parm.atoms[atom_2]
            a3 = parm.atoms[atom_3]

            atom_0_unique_name = self.universe.atoms[atom_0].resname + str(self.universe.atoms[atom_0].resid) + '_' + self.universe.atoms[atom_0].name + '_' + str(self.universe.atoms[atom_0].id)
            atom_1_unique_name = self.universe.atoms[atom_1].resname + str(self.universe.atoms[atom_1].resid) + '_' + self.universe.atoms[atom_1].name + '_' + str(self.universe.atoms[atom_1].id)
            atom_2_unique_name = self.universe.atoms[atom_2].resname + str(self.universe.atoms[atom_2].resid) + '_' + self.universe.atoms[atom_2].name + '_' + str(self.universe.atoms[atom_2].id)
            atom_3_unique_name = self.universe.atoms[atom_3].resname + str(self.universe.atoms[atom_3].resid) + '_' + self.universe.atoms[atom_3].name + '_' + str(self.universe.atoms[atom_3].id)

            atom_names = (self.universe.atoms[atom_0].name, self.universe.atoms[atom_1].name, self.universe.atoms[atom_2].name, self.universe.atoms[atom_3].name)

            # ix = None
            # for i, d in enumerate(self.universe.dihedrals):
            #     dihh_indices = (d.atoms[0].index, d.atoms[1].index, d.atoms[2].index, d.atoms[3].index)
            #     univ_indices = (atom_0, atom_1, atom_2, atom_3)
            #     if dihh_indices == univ_indices or dihh_indices == univ_indices[::-1]:
            #         ix = i
            #         break

            print(atom_3_unique_name, atom_2_unique_name, atom_1_unique_name, atom_0_unique_name, "aromaticity:", a0.aromatic, a1.aromatic, a2.aromatic, a3.aromatic)

            # # find the type
            # for res in PROTEIN_DIHEDRAL_SELECTION:
            #      for dihedral_type, dihedral_atoms in PROTEIN_DIHEDRAL_SELECTION[res].items():
            #         if list(atom_names) == dihedral_atoms or list(reversed(atom_names)) == dihedral_atoms:
            #             print("Found dihedral type:", dihedral_type, "between atoms:", atom_0_unique_name, atom_1_unique_name, atom_2_unique_name, atom_3_unique_name)

        exit()

        # Initialize the undirected molecular graph
        self.graph = nx.Graph()
        for bond in self.universe.bonds:
            i = bond.atoms[0].index
            j = bond.atoms[1].index
            self.graph.add_edge(i, j)

        # Build ring-closing bonds and remove them from the graph
        self.ring_closing_bonds = self.build_kinematic_ring_closing_bonds()
        for rcb in self.ring_closing_bonds:
            self.graph.remove_edge(rcb[0], rcb[1])

        # Perform BFS on the new graph
        # Nodes in the edges are expressed as prmtop indices
        self.edges = nx.bfs_edges(self.graph, source=root_index_prmtop)
        self.edges = list(self.edges)

        # Create a mapping between the original indices and the BFS-explored indices
        # Again, Molmodel adds atoms via Molmodel via bondAtom(idx1, idx2)
        # Essentially, we create a mapping between the order in which atoms were added to  and their original indices
        self.nodes = [root_index_prmtop] + [v for u, v in self.edges]
        self.prmtop_to_compound_atom_index_map = {}
        for compound_atom_index, prmtop_index in enumerate(self.nodes):
            self.prmtop_to_compound_atom_index_map[prmtop_index] = compound_atom_index

        # Get bonds as explored by BFS starting from the root atom
        # Bonds are expressed with prmtop indices
        self.all_bonds = [(parent, child, False) for parent, child in self.edges]

        # Make sure all nodes have been visited via BFS
        # Nodes are expressed as prmtop indices
        visited_nodes: Set[int] = set()
        for parent, child in self.edges:
            visited_nodes.add(parent)
            visited_nodes.add(child)
        assert len(visited_nodes) == len(self.universe.atoms), "Not all atoms were visited during BFS traversal of the molecular graph"

        # BFS doesn't yield edges if they extend back to an already explored node - a ring closing bond
        # At this point we don't want to add ring closing bonds, so we need to filter for them
        # Molmodel adds atoms by calling bondAtom(idx1,idx2) which expects idx1 to be already bonded via bondAtom(idx0,idx1)
        # Ring closing bonds add bonds to non-bonded atoms - bondAtom(idx2,idx3) - which yields errors
        self.all_bonds += [(bond[0], bond[1], True) for bond in self.ring_closing_bonds]

        assert len(self.all_bonds) == len(self.universe.bonds), "Number of bonds in the molecule does not match the number of bonds in the graph"

        # Build non-redundant bonds: phi, psi and user-specified chi (excluding omega and ring-closing)
        self.non_redundant_bonds = []
        for phi in self.phi_dihedrals:
            self.non_redundant_bonds.append((phi[1].index, phi[2].index))
        for psi in self.psi_dihedrals:
            self.non_redundant_bonds.append((psi[1].index, psi[2].index))
        for chi in self.chi_dihedrals:
            self.non_redundant_bonds.append((chi[1], chi[2]))

        # First residue index (1-based), starts at whatever residue index there is in the prmtop file
        first_resid_prmtop = self.universe.residues[0].resid
        last_resid_prmtop = self.universe.residues[-1].resid

        # MDAnalysis uses 0-based resid indexing internally and starts from 0 for each molecule
        first_resid_universe = 0
        last_resid_universe = len(self.universe.residues) - 1

        if not want_n_terminus_phi_rigid:
            if self.universe.residues[first_resid_universe].resname == 'ACE':
                h1 = self.universe.select_atoms(f"resid {first_resid_prmtop} and name H1")
                ch3 = self.universe.select_atoms(f"resid {first_resid_prmtop} and name CH3")
                c = self.universe.select_atoms(f"resid {first_resid_prmtop} and name C")
                n = self.universe.select_atoms(f"resid {first_resid_prmtop+1} and name N")
                phi = mda.AtomGroup([h1.ix[0], ch3.ix[0], c.ix[0], n.ix[0]], self.universe)
                assert phi is not None, "Could not build N-terminus phi dihedral for ACE capping group"
            else:
                h1 = self.universe.select_atoms(f"resid {first_resid_prmtop} and name H1")
                n = self.universe.select_atoms(f"resid {first_resid_prmtop} and name N")
                ca = self.universe.select_atoms(f"resid {first_resid_prmtop} and name CA")
                c = self.universe.select_atoms(f"resid {first_resid_prmtop} and name C")
                phi = mda.AtomGroup([h1.ix[0], n.ix[0], ca.ix[0], c.ix[0]], self.universe)
                assert phi is not None, "Could not build N-terminus phi dihedral for residue " + self.universe.residues[first_resid_universe].resname + str(first_resid_universe)
            self.non_redundant_bonds.append((phi.ix[1], phi.ix[2]))

        if not want_c_terminus_psi_rigid:
            if self.universe.residues[last_resid_universe].resname == 'NME':
                c_prev = self.universe.select_atoms(f"resid {last_resid_prmtop-1} and name C")
                n = self.universe.select_atoms(f"resid {last_resid_prmtop} and name N")
                c = self.universe.select_atoms(f"resid {last_resid_prmtop} and name C")
                h1 = self.universe.select_atoms(f"resid {last_resid_prmtop} and name H1")
                psi = mda.AtomGroup([c_prev.ix[0], n.ix[0], c.ix[0], h1.ix[0]], self.universe)
            else:
                n = self.universe.select_atoms(f"resid {last_resid_prmtop} and name N")
                ca = self.universe.select_atoms(f"resid {last_resid_prmtop} and name CA")
                c = self.universe.select_atoms(f"resid {last_resid_prmtop} and name C")
                oxt = self.universe.select_atoms(f"resid {last_resid_prmtop} and name OXT")
                psi = mda.AtomGroup([n.ix[0], ca.ix[0], c.ix[0], oxt.ix[0]], self.universe)
                assert psi is not None, "Could not build C-terminus psi dihedral for residue " + self.universe.residues[last_resid_universe].resname + str(last_resid_universe)
            self.non_redundant_bonds.append((psi.ix[1], psi.ix[2]))

    def prmtop_to_compound_atom_index(self, prmtop_index: int) -> int:
        return self.prmtop_to_compound_atom_index_map[prmtop_index]
    
    def get_nodes_as_prmtop_indices(self) -> List[int]:
        return self.nodes

    def get_bonds(self) -> List[Tuple[int, int, bool]]:
        return self.all_bonds

    def build_kinematic_ring_closing_bonds(self):
        # Ring-closing bonds should be chosen only from backbone mobilizers (phi, psi, omega) that lie on the kinematic loop
        # Never from side-chain or intra-residue dihedrals, including those in cyclic residues (Pro, Phe, Tyr, Trp, His)
        kinematic_cycles = []

        # Check if the node is part of multiple cycles
        # One cycle needs exactly one unique ring closing bond (can't have the same bond closing multiple cycles)
        is_node_available = {}

        for candidate_cycle in nx.cycle_basis(self.graph):
            residues_in_cycle = set()
            for atom_index in candidate_cycle:
                atom = self.universe.atoms[atom_index]
                is_node_available[atom_index] = True
                residues_in_cycle.add(atom.residue.resid)

            cycle_ids = []
            for atom_index in candidate_cycle:
                a = self.universe.atoms[atom_index]
                cycle_ids.append(a.index+1)
            
            if len(residues_in_cycle) > 1:
                kinematic_cycles.append(candidate_cycle)
        
        # filter nodes that first appear in the current kinematic cycle
        filtered_kinematic_cycles = []
        for cycle in kinematic_cycles:
            filtered_cycle = []
            for atom_index in cycle:
                if is_node_available[atom_index]:
                    filtered_cycle.append(atom_index)
                    is_node_available[atom_index] = False
            filtered_kinematic_cycles.append(filtered_cycle)

        # Build the list of ring-closing bonds
        ring_closing_bonds: List[Tuple[int, int]] = []

        # Start by looking at backbone omega dihedrals
        for cycle in filtered_kinematic_cycles:
            # Find all bonds in the cycle
            cycle_bonds = []
            for i in range(len(cycle)):
                atom1_index = cycle[i]
                atom2_index = cycle[(i + 1) % len(cycle)]

                # Cycles ends may not be directly bonded
                if self.graph.has_edge(atom1_index, atom2_index):
                    cycle_bonds.append((atom1_index, atom2_index))
            
            # select the first omega bond found in the cycle
            selected_bond = None
            for bond in cycle_bonds:
                for omega in self.omega_dihedrals:
                    bond_indices = tuple(sorted([bond[0], bond[1]]))
                    omega_indices = tuple(sorted([omega[1].index, omega[2].index]))

                    if bond_indices == omega_indices:
                        selected_bond = bond
                        break
                if selected_bond is not None:
                    break

            assert selected_bond is not None, "No suitable ring-closing bond found in the kinematic cycle"
            # if selected_bond is None:
            #     selected_bond = cycle_bonds[0]
            
            ring_closing_bonds.append(selected_bond)

        # Look for user-specified ring-closing dihedrals
        for res in self.universe.residues:
            for chi_name, chi_atoms in self.dihedral_sele[res.resname].items():
                chi_atom_ix = [self.universe.select_atoms(f"resid {res.resid} and name {atom_name}").ix[0] for atom_name in chi_atoms]
                chi_dihedral = mda.AtomGroup(chi_atom_ix, self.universe)

                assert len(chi_dihedral) == 4, f"Dihedral {chi_name} in residue {res.resname} {res.resid} does not have 4 atoms: {chi_atoms}, found {chi_dihedral}"

                if 'closing' in chi_name:
                    ring_closing_bonds.append((chi_dihedral[1].index, chi_dihedral[2].index))
                    continue

        return ring_closing_bonds