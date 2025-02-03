import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals


class BATCorrelations:
    def __init__(self, dcd_files, prmtop_file, inpcrd_file):

        dihedral_sele = {
            'ALA': {
                'chi1': '',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'ARG': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD',
                'chi3': 'name CB or name CG or name CD or name NE',
                'chi4': 'name CG or name CD or name NE or name CZ',
                'chi5': 'name CD or name NE or name CZ or name NH1'
            },
            'ASN': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name OD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'ASP': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name OD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'CYS': {
                'chi1': 'name N or name CA or name CB or name SG',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'CYX': {
                'chi1': 'name N or name CA or name CB or name SG',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'GLN': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD',
                'chi3': 'name CB or name CG or name CD or name OE1',
                'chi4': '',
                'chi5': ''
            },
            'GLU': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD',
                'chi3': 'name CB or name CG or name CD or name OE1',
                'chi4': '',
                'chi5': ''
            },
            'GLY': {
                'chi1': '',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'HIP': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name ND1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'HIE': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name ND1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'HID': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name ND1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'ILE': {
                'chi1': 'name N or name CA or name CB or name CG1',
                'chi2': 'name CA or name CB or name CG1 or name CD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'LEU': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'LYS': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD',
                'chi3': 'name CB or name CG or name CD or name CE',
                'chi4': 'name CG or name CD or name CE or name NZ',
                'chi5': ''
            },
            'MET': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name SD',
                'chi3': 'name CB or name CG or name SD or name CE',
                'chi4': '',
                'chi5': ''
            },
            'PHE': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'PRO': {
                'chi1': '',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'SER': {
                'chi1': 'name N or name CA or name CB or name OG',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'THR': {
                'chi1': 'name N or name CA or name CB or name OG1',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'TRP': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'TYR': {
                'chi1': 'name N or name CA or name CB or name CG',
                'chi2': 'name CA or name CB or name CG or name CD1',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            },
            'VAL': {
                'chi1': 'name N or name CA or name CB or name CG1',
                'chi2': '',
                'chi3': '',
                'chi4': '',
                'chi5': ''
            }
        }

        # Load the PDB file once
        self.universe = mda.Universe(prmtop_file, inpcrd_file)

        # Build and store the full molecular graph
        self.base_graph = self._build_graph()
        self.atom_masses = np.array([atom.mass for atom in self.universe.atoms])

            


        # will hold (dihedral_type, atom_group, atom_indices, residue_name)
        self.dihedral_types = []
        self.atom_indices = []
        self.residue_names = []
        self.residue_ids = []
        self.dihedral_values = []
        populated = False

        for dcd in dcd_files:
            universe = mda.Universe(prmtop_file, dcd)
            atom_groups = []
            
            # Get the dihedrals
            for res in universe.residues:
                # Not present on N-terminus (first residue)
                phi = res.phi_selection()
                if phi:
                    atom_groups.append(phi)
                    if not populated:
                        self.dihedral_types.append('phi')
                        self.atom_indices.append(phi.indices)
                        self.residue_names.append(res.resname)
                        self.residue_ids.append(res.resid)

                # Not present on C-terminus (last residue)
                psi = res.psi_selection()
                if psi:
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
                for chi in ['chi1', 'chi2', 'chi3', 'chi4', 'chi5']:
                    chi_names = dihedral_sele[res.resname][chi]
                    if not chi_names:
                        continue

                    sele = f"resid {res.resid} and ({chi_names})"
                    chi_dihedral = universe.select_atoms(sele)
                    atom_groups.append(chi_dihedral)
                    if not populated:
                        self.dihedral_types.append(chi)
                        self.atom_indices.append(chi_dihedral.indices)
                        self.residue_names.append(res.resname)
                        self.residue_ids.append(res.resid)
            
            populated = True

            # Compute dihedrals
            values = dihedrals.Dihedral(atom_groups).run().angles + 180
            values = np.deg2rad(values)
            values = values.astype(np.float32)

            # [num_dcds, num_dihedrals, num_frames]
            self.dihedral_values.append(values)

        self.numDihedrals = len(self.atom_indices)

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