import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import scipy.stats as stats
from scipy import linalg
import community as community_louvain

class BATCorrelations:
    def __init__(self, prmtop_file, inpcrd_file, seed):
        self.tol = 1e-6
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        # Load the PDB file once
        self.prmtop_file = prmtop_file
        self.inpcrd_file = inpcrd_file

        self.universe = mda.Universe(prmtop_file, inpcrd_file)

        # Build and store the full molecular graph
        self.base_graph = self._build_graph()
        self.atom_masses = np.array([atom.mass for atom in self.universe.atoms])

        # will hold (dihedral_type, atom_group, atom_indices, residue_name)
        atom_groups, dihedral_types, atom_indices, residue_names, residue_ids, num_dihedrals = self.build_dihedral_atom_groups(self.universe)
        self.dihedral_types = dihedral_types
        self.atom_indices = atom_indices
        self.residue_names = residue_names
        self.residue_ids = residue_ids
        self.num_dihedrals = num_dihedrals
        self.dihedral_values = []

    def compute_dihedrals_from_dcd(self, dcd_files, start, stop=None, step=1):
        for dcd in dcd_files:
            universe = mda.Universe(self.prmtop_file, dcd)
            atom_groups, dihedral_types, atom_indices, residue_names, residue_ids, num_dihedrals = self.build_dihedral_atom_groups(universe)

            # Compute dihedrals
            stop = stop if stop else universe.trajectory.n_frames
            values = dihedrals.Dihedral(atom_groups).run(start, stop, step).angles
            values = np.deg2rad(values)
            values = values.astype(np.float32)

            # [num_dcds, num_dihedrals, num_frames]
            self.dihedral_values.append(values)

        # self.dihedrals_to_pdb()

    def build_dihedral_atom_groups(self, universe):
        atom_groups = []
        dihedral_types = []
        atom_indices = []
        residue_names = []
        residue_ids = []

        for res in universe.residues:

            # The phi angle of the first residue is not defined
            phi = res.phi_selection()
            if not phi:
                h2 = universe.select_atoms(f"resid {res.resid} and name H2")
                n = universe.select_atoms(f"resid {res.resid} and name N")
                ca = universe.select_atoms(f"resid {res.resid} and name CA")
                c = universe.select_atoms(f"resid {res.resid} and name C")
                phi = mda.AtomGroup([h2.ix[0], n.ix[0], ca.ix[0], c.ix[0]], universe)

            atom_groups.append(phi)
            dihedral_types.append('phi')
            atom_indices.append(phi.indices)
            residue_names.append(res.resname)
            residue_ids.append(res.resid)

            # The psi angle of the last residue is not defined
            psi = res.psi_selection()
            if not psi:
                n = universe.select_atoms(f"resid {res.resid} and name N")
                ca = universe.select_atoms(f"resid {res.resid} and name CA")
                c = universe.select_atoms(f"resid {res.resid} and name C")
                oxt = universe.select_atoms(f"resid {res.resid} and name OXT")
                psi = mda.AtomGroup([n.ix[0], ca.ix[0], c.ix[0], oxt.ix[0]], universe)

            atom_groups.append(psi)
            dihedral_types.append('psi')
            atom_indices.append(psi.indices)
            residue_names.append(res.resname)
            residue_ids.append(res.resid)

            # Present in all residues
            if self.include_omega:
                omega = res.omega_selection()
                if omega:
                    atom_groups.append(omega)
                    dihedral_types.append('omega')
                    atom_indices.append(omega.indices)
                    residue_names.append(res.resname)
                    residue_ids.append(res.resid)

            # Chi angles
            for chi_name, chi_atoms in self.dihedral_sele[res.resname].items():
                chi_atom_ix = [universe.select_atoms(f"resid {res.resid} and name {atom_name}").ix[0] for atom_name in chi_atoms]
                chi_dihedral = mda.AtomGroup(chi_atom_ix, universe)

                atom_groups.append(chi_dihedral)
                dihedral_types.append(chi_name)
                atom_indices.append(chi_dihedral.indices)
                residue_names.append(res.resname)
                residue_ids.append(res.resid)

                print(f"Found dihedral {chi_name} in residue {res.resname} {res.resid} with atoms {chi_atoms}")

        # Find the disulfide bonds
        disulfide_bonds = set()
        for sg in universe.select_atoms("resname CYX and name SG"):
            for b in sg.bonds:
                if b.atoms[0].name == 'SG' and b.atoms[1].name == 'SG':
                    disulfide_bonds.add(tuple(sorted(b.indices)))
        disulfide_bonds = list(disulfide_bonds)

        # print the atoms in the disulfide bonds
        for bond in disulfide_bonds:
            sg0_atom = universe.atoms[bond[0]]
            sg1_atom = universe.atoms[bond[1]]

            sg_atom_pairs = [(sg0_atom, sg1_atom), (sg1_atom, sg0_atom)]
            for sg_atom_0, sg_atom_1 in sg_atom_pairs:
                sg0_selection = universe.select_atoms(f"resid {sg_atom_0.resid} and name SG")
                sg1_selection = universe.select_atoms(f"resid {sg_atom_1.resid} and name SG")
                cb_selection = universe.select_atoms(f"resid {sg_atom_1.resid} and name CB")
                ca_selection = universe.select_atoms(f"resid {sg_atom_1.resid} and name CA")

                atom_group = mda.AtomGroup([sg0_selection.atoms[0], sg1_selection.atoms[0], cb_selection.atoms[0], ca_selection.atoms[0]])
                atom_groups.append(atom_group)
                dihedral_types.append('chi2')
                atom_indices.append(atom_group.indices)
                residue_names.append('CYX')
                residue_ids.append(sg_atom_1.resid)

        num_dihedrals = len(atom_indices)

        return atom_groups, dihedral_types, atom_indices, residue_names, residue_ids, num_dihedrals

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

    def circcorr(self, dihs, pair, shuffle=False):
        """Calculates circular correlation coefficient for a single pair of dihedrals.
        
        If `shuffle=True`, independently shuffles each dihedral to break temporal correlation (null model).
        """
        ix, jx = pair

        x = dihs[:, ix]
        y = dihs[:, jx]

        # Optional shuffling for null model
        if shuffle:
            x = x.copy()
            y = y.copy()

            x = self.rng.permutation(x)
            y = self.rng.permutation(y)

        # Unwrap angles to handle discontinuities
        x = np.unwrap(x)
        y = np.unwrap(y)
        
        # Calculate circular variance for each series
        x_var = circstats.circvar(x)
        y_var = circstats.circvar(y)
        
        # Filter out series with low variance (high noise)
        if x_var < self.tol or y_var < self.tol:
            return ix, jx, 0.0
        
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
        

    def compute_correlations(self, null=False):

        correlation_matrices = []

        # [num_dcds, num_frames, num_dihedrals]
        # Compute the correlations for each DCD
        for dihedrals in self.dihedral_values:

            # Create argument list (pairs of dihedral indices)
            arg_list = []
            for ix in range(self.num_dihedrals):
                for jx in range(ix + 1, self.num_dihedrals):
                    arg_list.append((ix, jx))
                    # print(f"Computing correlation for dihedrals {ix} and {jx}: {self.circcorr(dihedrals, (ix, jx))}")

            correlation_matrix = np.ones((self.num_dihedrals, self.num_dihedrals))

            with ThreadPoolExecutor() as executor:
                futures = [executor.submit(self.circcorr, dihedrals, pair, null) for pair in arg_list]
                for future in as_completed(futures):
                    ix, jx, correlation = future.result()
                    correlation_matrix[ix][jx] = correlation
                    correlation_matrix[jx][ix] = correlation

            # self.validate_covariance_matrix(correlation_matrix)
            correlation_matrices.append(correlation_matrix)

        # [num_dcds, num_dihedrals, num_dihedrals]
        correlation_matrices = np.array(correlation_matrices)
        return correlation_matrices

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

    def chose_random_bonds(self, n):
        random_dihedrals = self.rng.choice(len(self.atom_indices), n, replace=False)
        random_bonds = []
        for r in random_dihedrals:
            aix1 = self.atom_indices[r][1]
            aix2 = self.atom_indices[r][2]
            random_bonds.append((min(aix1, aix2), max(aix1, aix2)))

        # sort in ascending order by first index
        random_bonds = sorted(random_bonds, key=lambda x: x[0])
        return random_bonds

    def chose_decoy_bonds(self, ref_bonds, steps, seed):
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
            index_diff = np.sum(np.abs(ref_bonds[:, 0] - random_bonds[:, 0]))
            if index_diff < best_index_diff:
                continue

            best_mass_diff = total_mass_diff
            best_bonds = random_bonds
            best_index_diff = index_diff

        return best_bonds
    
    def correlation_to_graph(self, corr):
        G = nx.Graph()
        n = corr.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                weight = abs(corr[i, j])
                G.add_edge(i, j, weight=weight)
        return G

    def chose_correlated_bonds(self, corr):
        G = self.correlation_to_graph(corr)
        partition = community_louvain.best_partition(G, weight='weight', resolution=1.0, random_state=self.seed)
        labels = np.array([partition[i] for i in range(len(partition))])
        modularity = community_louvain.modularity(partition, G)

        blocks = []
        for i in range(len(np.unique(labels))):
            block = np.where(labels == i)[0]

            block_corr = np.mean(corr[block][:, block], axis=0)
            percentile = np.percentile(block_corr, 50)
            high_group = block[np.where(block_corr >= percentile)[0]]
            low_group = block[np.where(block_corr < percentile)[0]]

            blocks.append(high_group)

        return blocks, modularity
