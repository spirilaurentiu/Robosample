import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals


class BATCorrelations:
    def __init__(self, dcd_files, prmtop_file):

        dihedral_sele = {
            'ALA': {
                'chi1': ['N', 'CA', 'CB', 'HB1'],
            },
            'ARG': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD'],
                'chi3': ['CB', 'CG', 'CD', 'NE'],
                'chi4': ['CG', 'CD', 'NE', 'CZ'],
                'chi5': ['CD', 'NE', 'CZ', 'NH1'],
            },
            'ASN': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'ND2'],
                'chi3': ['CB', 'CG', 'ND2', 'HD21'],
            },
            'ASP': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'OD1'],
            },
            'CYS': {
                'chi1': ['N', 'CA', 'CB', 'SG'],
                'chi2': ['CA', 'CB', 'SG', 'HG'],
            },
            'CYX': {
                'chi1': ['N', 'CA', 'CB', 'SG'],
            },
            'GLN': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD'],
                'chi3': ['CB', 'CG', 'CD', 'NE2'],
                'chi4': ['CG', 'CD', 'NE2', 'HE21'],
            },
            'GLU': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD'],
                'chi3': ['CB', 'CG', 'CD', 'OE1'],
            },
            'GLY': {},
            'HIP': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'ND1'],
            },
            'HIE': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'ND1'],
            },
            'HID': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'ND1'],
            },
            'ILE': {
                'chi1': ['N', 'CA', 'CB', 'CG1'],
                'chi2.1': ['CA', 'CB', 'CG1', 'CD1'],
                'chi2.2': ['CA', 'CB', 'CG2', 'HG21'],
                'chi3': ['CB', 'CG1', 'CD1', 'HD11'],
            },
            'LEU': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD1'],
                'chi3.1': ['CB', 'CG', 'CD1', 'HD11'],
                'chi3.2': ['CB', 'CG', 'CD2', 'HD21'],
            },
            'LYS': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD'],
                'chi3': ['CB', 'CG', 'CD', 'CE'],
                'chi4': ['CG', 'CD', 'CE', 'NZ'],
                'chi5': ['CD', 'CE', 'NZ', 'HZ1']
            },
            'MET': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'SD'],
                'chi3': ['CB', 'CG', 'SD', 'CE'],
                'chi4': ['CG', 'SD', 'CE', 'HE1'],
            },
            'PHE': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD1'],
            },
            'PRO': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD'],
                'chi3': ['CB', 'CG', 'CD', 'N'],
            },
            'SER': {
                'chi1': ['N', 'CA', 'CB', 'OG'],
                'chi2': ['CA', 'CB', 'OG', 'HG'],
            },
            'THR': {
                'chi1': ['N', 'CA', 'CB', 'OG1'],
                'chi2.1': ['CA', 'CB', 'OG1', 'HG1'],
                'chi2.2': ['CA', 'CB', 'CG2', 'HG21'],
            },
            'TRP': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD1'],
            },
            'TYR': {
                'chi1': ['N', 'CA', 'CB', 'CG'],
                'chi2': ['CA', 'CB', 'CG', 'CD1'],
                'chi3': ['CE1', 'CZ', 'OH', 'HH'],
            },
            'VAL': {
                'chi1': ['N', 'CA', 'CB', 'CG1'],
                'chi2.1': ['CA', 'CB', 'CG1', 'HG11'],
                'chi2.2': ['CA', 'CB', 'CG2', 'HG21'],
            }
        }

        # @TODO asp glu arg lys - protonated

        # will hold (dihedral_type, atom_group, atom_indices, residue_name)
        self.dihedral_types = []
        self.atom_indices = []
        self.residue_names = []
        self.residue_ids = []
        self.dihedral_values = []
        populated = False

        for dcd in dcd_files:

            universe = mda.Universe(prmtop_file, dcd)
            if not populated:
                # Build and store the full molecular graph
                self.universe = universe
                self.base_graph = self._build_graph()
                self.atom_masses = np.array([atom.mass for atom in self.universe.atoms])
            
            # Get the dihedrals
            atom_groups = []
            for res in universe.residues:

                # The phi angle of the first residue is not defined
                phi = res.phi_selection()
                if not phi:
                    h1 = universe.select_atoms(f"resid {res.resid} and name H1")
                    n = universe.select_atoms(f"resid {res.resid} and name N")
                    ca = universe.select_atoms(f"resid {res.resid} and name CA")
                    c = universe.select_atoms(f"resid {res.resid} and name C")
                    phi = mda.AtomGroup([h1.ix[0], n.ix[0], ca.ix[0], c.ix[0]], universe)

                atom_groups.append(phi)
                if not populated:
                    self.dihedral_types.append('phi')
                    self.atom_indices.append(phi.indices)
                    self.residue_names.append(res.resname)
                    self.residue_ids.append(res.resid)

                # The psi angle of the last residue is not defined
                psi = res.psi_selection()
                if not psi:
                    n = universe.select_atoms(f"resid {res.resid} and name N")
                    ca = universe.select_atoms(f"resid {res.resid} and name CA")
                    c = universe.select_atoms(f"resid {res.resid} and name C")
                    oxt = universe.select_atoms(f"resid {res.resid} and name OXT")
                    psi = mda.AtomGroup([n.ix[0], ca.ix[0], c.ix[0], oxt.ix[0]], universe)

                atom_groups.append(psi)
                if not populated:
                    self.dihedral_types.append('psi')
                    self.atom_indices.append(psi.indices)
                    self.residue_names.append(res.resname)
                    self.residue_ids.append(res.resid)

                # Present in all residues
                omega = res.omega_selection()
                if omega:
                    atom_groups.append(omega)
                    if not populated:
                        self.dihedral_types.append('omega')
                        self.atom_indices.append(omega.indices)
                        self.residue_names.append(res.resname)
                        self.residue_ids.append(res.resid)

                # Chi angles
                for chi_name, chi_atoms in dihedral_sele[res.resname].items():
                    chi_atom_ix = [universe.select_atoms(f"resid {res.resid} and name {atom_name}").ix[0] for atom_name in chi_atoms]
                    chi_dihedral = mda.AtomGroup(chi_atom_ix, universe)
                    atom_groups.append(chi_dihedral)
                    if not populated:
                        self.dihedral_types.append(chi_name)
                        self.atom_indices.append(chi_dihedral.indices)
                        self.residue_names.append(res.resname)
                        self.residue_ids.append(res.resid)

            # Find the disulfide bonds
            disulfide_bonds = set()
            for sg in self.universe.select_atoms("resname CYX and name SG"):
                for b in sg.bonds:
                    if b.atoms[0].name == 'SG' and b.atoms[1].name == 'SG':
                        disulfide_bonds.add(tuple(sorted(b.indices)))
            disulfide_bonds = list(disulfide_bonds)

            # print the atoms in the disulfide bonds
            for bond in disulfide_bonds:
                sg0_atom = self.universe.atoms[bond[0]]
                sg1_atom = self.universe.atoms[bond[1]]

                sg_atom_pairs = [(sg0_atom, sg1_atom), (sg1_atom, sg0_atom)]
                for sg_atom_0, sg_atom_1 in sg_atom_pairs:
                    sg0_selection = self.universe.select_atoms(f"resid {sg_atom_0.resid} and name SG")
                    sg1_selection = self.universe.select_atoms(f"resid {sg_atom_1.resid} and name SG")
                    cb_selection = self.universe.select_atoms(f"resid {sg_atom_1.resid} and name CB")
                    ca_selection = self.universe.select_atoms(f"resid {sg_atom_1.resid} and name CA")

                    atom_group = mda.AtomGroup([sg0_selection.atoms[0], sg1_selection.atoms[0], cb_selection.atoms[0], ca_selection.atoms[0]])
                    atom_groups.append(atom_group)

                    if not populated:
                        self.dihedral_types.append('chi2')
                        self.atom_indices.append(atom_group.indices)
                        self.residue_names.append('CYX')
                        self.residue_ids.append(sg_atom_1.resid)
            populated = True

            # Compute dihedrals
            values = dihedrals.Dihedral(atom_groups).run().angles + 180
            values = np.deg2rad(values)
            values = values.astype(np.float32)

            # [num_dcds, num_dihedrals, num_frames]
            self.dihedral_values.append(values)

        self.numDihedrals = len(self.atom_indices)

    def get_dihedral_atom_indices(self):
        indices = []
        for aix0, aix1, aix2, aix3 in self.atom_indices:
            indices.append((aix1, aix2))
        return indices

    def calculate_correlation(self, dihs, pair):
        """Calculates circular correlation coefficient for a single pair of dihedrals."""
        ix, jx = pair

        corr = self.dihedralsCorrelation_claude(dihs[:, ix], dihs[:, jx])
        return (ix, jx, corr)


    def compute_correlations(self):
                
        correlation_matrices = []

        # Compute the correlations for each DCD
        for dihedrals in self.dihedral_values:

            # Create argument list (pairs of dihedral indices)
            arg_list = []
            for ix in range(self.numDihedrals):
                for jx in range(ix + 1, self.numDihedrals):
                    arg_list.append((ix, jx))

            correlation_matrix = np.ones((self.numDihedrals, self.numDihedrals))
            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.calculate_correlation, 
                dihedrals, pair) for pair in arg_list]
                for future in as_completed(futures):
                    ix, jx, corr = future.result()
                    correlation_matrix[ix][jx] = corr
                    correlation_matrix[jx][ix] = corr

            correlation_matrices.append(correlation_matrix)

        correlation_matrices = np.array(correlation_matrices)
        return correlation_matrices

    def cluster(self, correlation_matrix, n_clusters=5):

        # Perform hierarchical clustering
        linkage = sch.linkage(correlation_matrix, method='ward', optimal_ordering=False)
        
        # Perform clustering
        cluster_labels = sch.fcluster(linkage, t=n_clusters, criterion='maxclust')

        filtered_clusters = [[] for _ in range(n_clusters)]
        for i, c in enumerate(cluster_labels):
            if self.dihedral_types[i] == 'omega':
                continue

            aix1 = self.atom_indices[i][1]
            aix2 = self.atom_indices[i][2]
            dihedral_type = self.dihedral_types[i]

            filtered_clusters[c-1].append((aix1, aix2, dihedral_type))

        return filtered_clusters

    def _build_graph(self):
        """Precomputes the molecular graph without any cuts."""
        G = nx.Graph()
        for atom in self.universe.atoms:
            G.add_node(atom.index)
        for bond in self.universe.bonds:
            idx1, idx2 = bond[0].index, bond[1].index
            G.add_edge(idx1, idx2)
        return G

    def compute_protein_cut_masses(self, cut_bonds, num_components=None):
        # Create a copy of the base graph (shallow copy, preserving structure)
        G = self.base_graph = self._build_graph()

        # Remove cut bonds
        if cut_bonds:
            G.remove_edges_from(cut_bonds)

        # Find connected components
        parts = list(nx.connected_components(G))
        if num_components and len(parts) != num_components:
            # print(f"Warning: found {len(parts)} components, expected {num_components}")
            return None
        else:
            # Compute masses using NumPy
            # print(f"Found {len(parts)} components")
            return np.array([self.atom_masses[list(part)].sum() for part in parts])

    def chose_random_bonds(self, n=10):
        random_dihedrals = np.random.choice(len(self.atom_indices), n, replace=False)
        random_bonds = []
        for r in random_dihedrals:
            aix1 = self.atom_indices[r][1]
            aix2 = self.atom_indices[r][2]
            random_bonds.append((min(aix1, aix2), max(aix1, aix2)))

        # sort in ascending order by first index
        random_bonds = sorted(random_bonds, key=lambda x: x[0])
        return random_bonds
    
    def dihedralsCorrelation_claude(self, x, y, variance_threshold=1e-5):
        """Circular correlation between two dihedral series"""

        # Unwrap angles to handle discontinuities
        x = np.unwrap(x)
        y = np.unwrap(y)
        
        # Calculate circular variance for each series
        x_var = circstats.circvar(x)
        y_var = circstats.circvar(y)
        
        # Filter out series with low variance (high noise)
        if x_var < variance_threshold or y_var < variance_threshold:
            return 0  # Not enough meaningful variation
        
        # Compute circular correlation
        x_centered = x - np.mean(x)
        y_centered = y - np.mean(y)
        
        numerator = np.sum(np.sin(x_centered) * np.sin(y_centered))
        denominator = np.sqrt(np.sum(np.sin(x_centered)**2) * np.sum(np.sin(y_centered)**2))
        
        return numerator / denominator

    def chose_decoy_bonds(self, ref_bonds, steps=3000):
        ref_masses = self.compute_protein_cut_masses(ref_bonds)
        np.ndarray.sort(ref_masses)

        best_bonds = None
        best_mass_diff = 1e10
        best_index_diff = 0

        for i in range(steps):

            # Check that we have the same number of masses
            random_bonds = self.chose_random_bonds(len(ref_bonds))
            random_masses = self.compute_protein_cut_masses(random_bonds, len(ref_masses))
            while random_masses is None:
                random_bonds = self.chose_random_bonds(len(ref_bonds))
                random_masses = self.compute_protein_cut_masses(random_bonds, len(ref_masses))
            np.ndarray.sort(random_masses)

            # Check that we have the same minimum mass
            if np.min(ref_masses) != np.min(random_masses):
                continue

            # Calculate total mass difference
            total_mass_diff = np.sum(np.abs(ref_masses - random_masses))
            if total_mass_diff >= best_mass_diff:
                continue

            # Calculate index difference
            index_diff = 0
            for r, b in zip(ref_bonds, random_bonds):
                index_diff += abs(r[0] - b[0])
            if index_diff < best_index_diff:
                continue

            best_mass_diff = total_mass_diff
            best_bonds = random_bonds
            best_index_diff = index_diff

        return best_bonds

# if __name__ == '__main__':
#     # for pdbid in ['1A5E', '1APQ', '1BLA', '2JR6', '2KB1', '2KCU', '2L29', '2L3B', '2LGJ']:
#     for pdbid in ['1A5E']:
#         test_prmtop = f"data-raw/{pdbid}.prmtop"
#         test_dcds = [f"md/{pdbid}_0.dcd", f"md/{pdbid}_1.dcd", f"md/{pdbid}_2.dcd", f"md/{pdbid}_3.dcd", f"md/{pdbid}_4.dcd"]
#         flexor = BATCorrelations(test_dcds, test_prmtop)

#         # corr = flexor.compute_correlations()
#         corr = np.load(f"results/{pdbid}_correlations.npy")
#         abs_corr = np.abs(corr)
#         mean_corr = np.mean(abs_corr, axis=0)
#         max_corr = np.max(abs_corr, axis=0)

#         clusters = flexor.cluster(max_corr)
#         for c in clusters:
#             bond_list = []
#             for aix1, aix2, dihedral_type in c:
#                 bond_list.append((aix1, aix2))