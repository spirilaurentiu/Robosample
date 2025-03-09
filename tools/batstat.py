import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import scipy.stats as stats
from scipy import linalg

class BATCorrelations:
    def __init__(self, dcd_files, prmtop_file, inpcrd_file):

        self.tol = 1e-6

        dihedral_sele = {
            'ALA': {
                'chi1': ['N', 'CA', 'CB', 'HB1'],
            },
            'ARG': { # Has a terminal guanidinium group that will not move
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

        # Load the PDB file once
        self.universe = mda.Universe(prmtop_file) #, inpcrd_file

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
            # universe = mda.Universe(prmtop_file, inpcrd_file)
            atom_groups = []
            
            # Get the dihedrals
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

            if not populated:
                self.numDihedrals = len(self.atom_indices)
            populated = True

            # Compute dihedrals
            values = dihedrals.Dihedral(atom_groups).run().angles # + 180
            values = np.deg2rad(values)
            values = values.astype(np.float32)

            # [num_dcds, num_dihedrals, num_frames]
            self.dihedral_values.append(values)

            # break

            # # print atom indices and dihedral values
            # for i in range(len(self.atom_indices)):
            #     aix0 = self.atom_indices[i][0]
            #     aix1 = self.atom_indices[i][1]
            #     aix2 = self.atom_indices[i][2]
            #     aix3 = self.atom_indices[i][3]
            #     val = np.rad2deg(values[0][i]) - 180
            #     t = self.dihedral_types[i]
            #     print(f"measure dihed {{{aix0} {aix1} {aix2} {aix3}}} \t {val} \t {t}")

        # self.dihedrals_to_pdb()
        # exit()

    def get_dihedral_atom_indices(self):
        indices = []
        for aix0, aix1, aix2, aix3 in self.atom_indices:
            indices.append((aix1, aix2))
        return indices

    def dihedrals_to_pdb(self):
        with open('dihedrals.pdb', 'w') as f:

            hetatm = []
            conect = []
            
            for aix, dihedral_type in zip(self.atom_indices, self.dihedral_types):
                aix1 = aix[1]
                aix2 = aix[2]

                atom1 = self.universe.atoms[aix1]
                atom2 = self.universe.atoms[aix2]

                hetatm.append(f"HETATM{aix1+1:5d} {dihedral_type[:4]:4s} {atom1.resname:3s} {atom1.resid:4d}    {atom1.position[0]:8.3f}{atom1.position[1]:8.3f}{atom1.position[2]:8.3f}  1.00  0.00          {atom1.name:2s}")
                hetatm.append(f"HETATM{aix2+2:5d} {dihedral_type[:4]:4s} {atom2.resname:3s} {atom2.resid:4d}    {atom2.position[0]:8.3f}{atom2.position[1]:8.3f}{atom2.position[2]:8.3f}  1.00  0.00          {atom2.name:2s}")
                conect.append(f"CONECT{aix1+1:5d}{aix2+2:5d}")
            
            f.write("MODEL\n")
            for line in hetatm:
                f.write(line + '\n')
            for line in conect:
                f.write(line + '\n')
            f.write("ENDMDL\n")

    def circcorr(self, dihs, pair):
        """Calculates circular correlation coefficient for a single pair of dihedrals."""
        ix, jx = pair

        x = dihs[:, ix]
        y = dihs[:, jx]

        # Unwrap angles to handle discontinuities
        x = np.unwrap(dihs[:, ix])
        y = np.unwrap(dihs[:, jx])
        
        # Calculate circular variance for each series
        x_var = circstats.circvar(x)
        y_var = circstats.circvar(y)
        
        # Filter out series with low variance (high noise)
        if x_var < self.tol or y_var < self.tol:
            return ix, jx, 0
        
        # Compute circular correlation
        x_centered = x - circstats.circmean(x)
        y_centered = y - circstats.circmean(y)
        
        numerator = np.sum(np.sin(x_centered) * np.sin(y_centered))
        denominator = np.sqrt(np.sum(np.sin(x_centered)**2) * np.sum(np.sin(y_centered)**2))
        correlation = numerator / denominator
        
        return ix, jx, correlation
    
    def validate_covariance_matrix(self, covariance_matrix):
        if covariance_matrix.shape[0] != covariance_matrix.shape[1]:
            raise ValueError("Covariance matrix is not square.")
        
        if not np.allclose(covariance_matrix, covariance_matrix.T, atol=self.tol):
            raise ValueError("Covariance matrix is not symmetric.")
        
        if not np.all(np.diag(covariance_matrix) >= self.tol):
            raise ValueError("Covariance matrix has non-positive diagonal elements.")
        
        try:
            eigenvalues = np.linalg.eigvalsh(covariance_matrix)
            linalg.cholesky(covariance_matrix)
        except linalg.LinAlgError:
            raise ValueError("Cholesky decomposition failed: covariance matrix is not positive definite.")
        

    def compute_correlations(self):

        correlation_matrices = []

        # [num_dcds, num_frames, num_dihedrals]
        # Compute the correlations for each DCD
        for dihedrals in self.dihedral_values:

            # Create argument list (pairs of dihedral indices)
            arg_list = []
            for ix in range(self.numDihedrals):
                for jx in range(ix + 1, self.numDihedrals):
                    arg_list.append((ix, jx))
                    # print(f"Computing correlation for dihedrals {ix} and {jx}: {self.circcorr(dihedrals, (ix, jx))}")

            correlation_matrix = np.ones((self.numDihedrals, self.numDihedrals))
            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.circcorr, dihedrals, pair) for pair in arg_list]
                for future in as_completed(futures):
                    ix, jx, correlation = future.result()
                    correlation_matrix[ix][jx] = correlation
                    correlation_matrix[jx][ix] = correlation

            self.validate_covariance_matrix(correlation_matrix)
            correlation_matrices.append(correlation_matrix)

        # [num_dcds, num_dihedrals, num_dihedrals]
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
    
    def consensus_clustering(self, heatmaps, num_clusters):
        """
        Perform consensus clustering from multiple correlation heatmaps.
        
        Parameters:
            heatmaps (list of np.array): List of NxN correlation matrices.
            num_clusters (int): Number of clusters to form.
        
        Returns:
            final_labels (np.array): Cluster labels for each datapoint.
        """
        n = heatmaps[0].shape[0]  # Number of datapoints
        num_matrices = len(heatmaps)
        co_occurrence = np.zeros((n, n))
        
        for heatmap in heatmaps:
            # Convert similarity (correlation) into a distance metric
            distance_matrix = 1 - heatmap
            
            # Perform hierarchical clustering
            linkage_matrix = sch.linkage(distance_matrix, method='average')
            labels = sch.fcluster(linkage_matrix, num_clusters, criterion='maxclust')
            
            # Construct co-occurrence matrix
            for i in range(n):
                for j in range(n):
                    if labels[i] == labels[j]:
                        co_occurrence[i, j] += 1
        
        # Normalize co-occurrence matrix
        co_occurrence /= num_matrices
        
        # Perform final clustering on co-occurrence matrix
        final_distance_matrix = 1 - co_occurrence
        final_linkage = sch.linkage(final_distance_matrix, method='average')
        final_labels = sch.fcluster(final_linkage, num_clusters, criterion='maxclust')

        filtered_clusters = [[] for _ in range(num_clusters)]
        for i, c in enumerate(final_labels):
            if self.dihedral_types[i] == 'omega':
                continue

            aix1 = self.atom_indices[i][1]
            aix2 = self.atom_indices[i][2]
            dihedral_type = self.dihedral_types[i]

            filtered_clusters[c-1].append((aix1, aix2, dihedral_type))

        return filtered_clusters, co_occurrence
    
    def cluster_dihedrals(self, correlation_matrix, alpha=0.5, n_clusters=5):
        from sklearn.cluster import SpectralClustering

        """
        Partitions the dihedrals into blocks based on both marginal and partial correlations.
        
        Parameters:
        - correlation_matrix: NxN numpy array (assumed symmetric and invertible)
        - threshold: distance threshold for forming clusters (lower => more clusters)
        - alpha: weight for marginal correlations (0 <= alpha <= 1). (1-alpha) is the weight for partial correlations.
        - n_clusters: number of clusters to form
        
        Returns:
        - clusters: list of clusters, each cluster is a list of indices
        - distance_matrix: the distance matrix used for clustering
        - partial_corr: computed partial correlation matrix
        - combined_corr: the weighted combination of marginal and partial correlations
        """

        # Compute precision matrix (inverse of correlation matrix)
        try:
            # Try regularized inverse first
            from sklearn.covariance import GraphicalLassoCV
            model = GraphicalLassoCV(alphas=[alpha])
            model.fit(correlation_matrix)
            precision = model.precision_
            print(f"Graphical lasso succeeded with alpha={alpha}")
        except:
            # Fallback to pseudo-inverse if graphical lasso fails
            precision = np.linalg.inv(correlation_matrix)
            print(f"Graphical lasso failed, using pseudo-inverse")

        # Compute partial correlations: 
        diag = np.diag(precision)
        partial_corr = -precision / np.sqrt(np.outer(diag, diag))
        np.fill_diagonal(partial_corr, 1.0)  # Set diagonal to 1

        # Combine marginal and partial correlations using geometric mean
        # This linear combination lets you balance overall dependency and direct coupling.
        combined_corr = np.sign(correlation_matrix) * np.sqrt(np.abs(correlation_matrix * partial_corr))
        
        # Define a distance: higher combined correlation -> shorter distance.
        distance_matrix = 1 - combined_corr

        def optimal_cluster_count(distance_matrix, max_clusters=100):
            from sklearn.metrics import silhouette_score
            scores = []
            for k in range(2, max_clusters + 1):
                clustering = SpectralClustering(n_clusters=k, affinity='precomputed', random_state=42)
                labels = clustering.fit_predict(distance_matrix)
                score = silhouette_score(distance_matrix, labels, metric='precomputed')
                scores.append(score)

            # smooth out scores
            kernel_size = int(max_clusters / 5)
            scores = np.array(scores)
            scores = np.convolve(scores, np.ones(kernel_size) / kernel_size, mode='valid')

            from kneed import KneeLocator
            kn = KneeLocator(range(2, scores.shape[0] + 2), scores, curve='convex', direction='decreasing')
            best_k = kn.knee

            # import matplotlib.pyplot as plt
            # plt.plot(range(2, scores.shape[0] + 2), scores)
            # plt.plot(best_k, scores[best_k-2], 'ro')
            # plt.show()

            # best_k = np.argmax(scores) + 2  # because range starts from 2
            return best_k
        
        best_k = optimal_cluster_count(distance_matrix)

        # Perform spectral clustering
        n_clusters = best_k
        spectral = SpectralClustering(n_clusters=n_clusters, affinity='precomputed', random_state=42)
        cluster_labels = spectral.fit_predict(distance_matrix)

        # # Hierarchical clustering using Ward's method
        # from sklearn.cluster import AgglomerativeClustering
        # clustering = AgglomerativeClustering(
        #     n_clusters=n_clusters,
        #     # metric='precomputed',
        #     linkage='ward',
        #     # distance_threshold=0.3
        # )
        # cluster_labels = clustering.fit_predict(distance_matrix)


        # Create a list of clusters, each cluster is a list of indices
        clusters = [[] for _ in range(cluster_labels.max() + 1)]
        for idx, label in enumerate(cluster_labels):
            aix1 = self.atom_indices[idx][1]
            aix2 = self.atom_indices[idx][2]
            dihedral_type = self.dihedral_types[idx]
            clusters[label].append((aix1, aix2, dihedral_type))

        return clusters
    
    def greedy_collapse(self, corr_matrix, degrees, num_vars, alpha, gamma, threshold):
        """
        Heuristic greedy collapsing based on Venugopal & Gogate's paper.
        Select variables for collapsing based on correlation and edge cost.
        """
        collapsed = set()
        edges_added = 0

        # Score each variable's "collapsibility"
        scores = np.zeros(num_vars)
        for i in range(num_vars):
            scores[i] = np.sum(corr_matrix[i]) / num_vars  # average correlation with all others

        while True:
            # Select the most collapsible variable that doesn't violate alpha-degree constraint
            candidates = [i for i in range(num_vars) if i not in collapsed and degrees[i] <= alpha]
            if not candidates:
                break

            best_var = max(candidates, key=lambda x: scores[x])

            # Estimate edges added if collapsing this variable
            neighbors = np.where(corr_matrix[best_var] > threshold)[0]
            new_edges = len(neighbors) * (len(neighbors) - 1) // 2  # clique formation

            if edges_added + new_edges > gamma:
                break

            # Collapse this variable
            collapsed.add(best_var)
            edges_added += new_edges

            # Update degrees (removing best_var reduces degrees of its neighbors)
            for neighbor in neighbors:
                degrees[neighbor] -= 1

        return collapsed

    def greedy_block(self, corr_matrix, collapsed, beta, threshold):
        """
        Greedy heuristic to form blocks of variables to be sampled jointly.
        """
        num_vars = corr_matrix.shape[0]
        blocks = [{i} for i in range(num_vars) if i not in collapsed]

        def block_score(block):
            indices = list(block)
            if len(indices) < 2:
                return 0
            pairwise_corrs = [corr_matrix[i, j] for i in block for j in block if i < j]
            mean_corr = np.mean(pairwise_corrs) if pairwise_corrs else 0
            if mean_corr < threshold:
                return -np.inf  # Reject low-cohesion blocks
            return np.sum(pairwise_corrs)  # Original sum score


        merged = True
        while merged:
            merged = False
            best_merge = None
            best_score = -np.inf

            for i, block1 in enumerate(blocks):
                for j, block2 in enumerate(blocks):
                    if i >= j:
                        continue

                    # Check if merge respects treewidth limit (here simplified as block size)
                    if len(block1 | block2) > beta:
                        continue

                    merged_score = block_score(block1 | block2)
                    if merged_score > best_score:
                        best_score = merged_score
                        best_merge = (i, j)

            if best_merge is None:
                break

            # Merge blocks
            i, j = best_merge
            blocks[i] = blocks[i] | blocks[j]
            del blocks[j]
            merged = True

        return blocks

    def refine_blocks(self, blocks, collapsed, corr_matrix, threshold):
        """
        After greedy block merging, move 'dead' blocks (low internal correlation)
        into the collapsed set.
        """
        refined_blocks = []
        for block in blocks:
            if len(block) < 5:
                # Single dihedral block - candidate for collapsing.
                collapsed.update(block)
                continue

            # Check average intra-block correlation
            pairwise_corrs = [corr_matrix[i, j] for i in block for j in block if i < j]
            mean_corr = np.mean(pairwise_corrs) if pairwise_corrs else 0

            if mean_corr < threshold:
                # Weak block - better to collapse.
                collapsed.update(block)
            else:
                refined_blocks.append(block)

        return refined_blocks, collapsed

    def print_block_correlations(self, blocks, collapsed, corr_matrix):
        def mean_corr(block):
            indices = list(block)
            if len(indices) < 2:
                return 0
            pairwise_corrs = [corr_matrix[i, j] for i in block for j in block if i < j]
            return np.mean(pairwise_corrs) if pairwise_corrs else 0

        print("\nBlock intra-correlations:")
        for idx, block in enumerate(blocks):
            print(f"  Block {idx}: mean intra-correlation = {mean_corr(block):.4f} with {len(block)} variables")

        collapsed_list = list(collapsed)
        if len(collapsed_list) > 1:
            pairwise_corrs = [corr_matrix[i, j] for i in collapsed for j in collapsed if i < j]
            collapsed_intra_corr = np.mean(pairwise_corrs) if pairwise_corrs else 0
        else:
            collapsed_intra_corr = 0

        print(f"\nCollapsed intra-correlation: {collapsed_intra_corr:.4f} with {len(collapsed)} variables")
    
    def dynamic_partitioning(self, corr_matrix):
        """
        Full pipeline to partition dihedrals into blocks and collapsed variables.
        """

        THRESHOLD = 0.1

        num_vars = corr_matrix.shape[0]
        degrees = np.zeros(num_vars, dtype=int)
        for i in range(num_vars):
            degrees[i] = np.count_nonzero(corr_matrix[i] > THRESHOLD)
        average_degree = np.mean(degrees)

        ALPHA = int(0.5 * average_degree)  # Allow collapsing of hubs that aren't "too large"
        BETA = int(0.5 * average_degree)   # Blocks can tolerate a bit more density
        GAMMA = 50 * ALPHA                 # Penalty grows with allowed collapse degree

        collapsed = self.greedy_collapse(corr_matrix, degrees, num_vars, ALPHA, GAMMA, THRESHOLD)
        blocks = self.greedy_block(corr_matrix, collapsed, BETA, THRESHOLD)
        blocks, collapsed = self.refine_blocks(blocks, collapsed, corr_matrix, THRESHOLD)

        # check if all dihedrals are in blocks
        all_vars = set(range(num_vars))
        for block in blocks:
            all_vars -= block
        all_vars -= collapsed
        assert not all_vars, f"Error: {len(all_vars)} variables not included in blocks or collapsed set"

        self.print_block_correlations(blocks, collapsed, corr_matrix)

        blocks_atom_list = []
        for block in blocks:
            atom_list = []
            for dihedral_index in list(block):
                aix1 = self.atom_indices[dihedral_index][1]
                aix2 = self.atom_indices[dihedral_index][2]
                atom_list.append((aix1, aix2))
            blocks_atom_list.append(atom_list)

        collapsed_atom_list = []
        for dihedral_index in list(collapsed):
            aix1 = self.atom_indices[dihedral_index][1]
            aix2 = self.atom_indices[dihedral_index][2]
            collapsed_atom_list.append((aix1, aix2))

        return blocks_atom_list, collapsed_atom_list

    def chose_decoy_bonds(self, ref_bonds, steps=3000):
        ref_masses = self.compute_protein_cut_masses(ref_bonds)
        np.ndarray.sort(ref_masses)

        best_bonds = None
        best_mass_diff = 1e10
        best_index_diff = 0

        for i in range(steps):

            # Check that we have the same number of masses
            random_bonds = stats.chose_random_bonds(len(ref_bonds))
            random_masses = stats.compute_protein_cut_masses(random_bonds, len(ref_masses))
            while random_masses is None:
                random_bonds = stats.chose_random_bonds(len(ref_bonds))
                random_masses = stats.compute_protein_cut_masses(random_bonds, len(ref_masses))
            np.ndarray.sort(random_masses)

            # Check that we have the same minimum mass
            if np.min(ref_masses) != np.min(random_masses):
                continue

            # Calculate total mass difference
            total_mass_diff = np.sum(np.abs(ref_masses - random_masses))
            if total_mass_diff >= best_mass_diff:
                continue

            # Calculate index difference
            index_diff = np.sum(np.abs(ref_bonds[:, 0] - random_bonds[:, 0]))
            if index_diff < best_index_diff:
                continue

            best_mass_diff = total_mass_diff
            best_bonds = random_bonds
            best_index_diff = index_diff

        return best_bonds