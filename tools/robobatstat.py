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
    def __init__(self, prmtop_file, inpcrd_file, seed, include_omega=False, include_chi1=False, include_chi2=False, include_chi3=False, include_chi4=False, include_chi5=False):
        self.dihedral_sele = {
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

        # remove unwanted dihedrals
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
                
        # @TODO asp glu arg lys - protonated

        self.include_omega = include_omega
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
    
if __name__ == '__main__':

    # test()
    # exit()

    # prmtop_file = 'data-raw/1A5E.prmtop'
    # dcd_files = ['md/1A5E_0.dcd']

    # stats = BATCorrelations(dcd_files, prmtop_file)

    # exit()

    # ref_bonds = stats.chose_random_bonds(n=100)
    # ref_masses = stats.compute_protein_cut_masses(ref_bonds)
    # np.ndarray.sort(ref_masses)

    # print(f"Reference masses: {ref_masses}")

    # best_bonds = None
    # best_masses = None
    # best_mass_diff = 1e10
    # best_index_diff = 0

    # import tqdm
    # for i in tqdm.tqdm(range(10000)):

    #     while True:
    #         # Check that we have the same number of masses
    #         random_bonds = stats.chose_random_bonds(len(ref_bonds))
    #         random_masses = stats.compute_protein_cut_masses(random_bonds, len(ref_masses))
    #         while random_masses is None:
    #             random_bonds = stats.chose_random_bonds(len(ref_bonds))
    #             random_masses = stats.compute_protein_cut_masses(random_bonds, len(ref_masses))
    #         np.ndarray.sort(random_masses)

    #         # # Check that we have different bonds
    #         # if np.any(np.isin(ref_bonds, random_bonds)):
    #         #     continue

    #         break

    #     # Check that we have the same minimum mass
    #     if np.min(ref_masses) != np.min(random_masses):
    #         continue

    #     # Calculate total mass difference
    #     total_mass_diff = np.sum(np.abs(ref_masses - random_masses))
    #     if total_mass_diff >= best_mass_diff:
    #         continue

    #     # Calculate index difference
    #     index_diff = np.sum(np.abs(ref_bonds[:, 0] - random_bonds[:, 0]))
    #     if index_diff < best_index_diff:
    #         continue

    #     best_mass_diff = total_mass_diff
    #     best_bonds = random_bonds
    #     best_masses = random_masses
    #     best_index_diff = index_diff
    #     print(f"New best mass difference: {best_mass_diff}")

    # print(f"Best mass difference: {best_mass_diff}")
    # print(f"Reference masses: {ref_masses}")
    # print(f"Best masses: {best_masses}")
    # print("Reference indices:", ref_bonds)
    # print("Best indices:", best_bonds)

    # exit()

    # for pdbid in ['1A5E', '1APQ', '1BLA', '2JR6', '2KB1', '2KCU', '2L29', '2L3B', '2LGJ']:
    for pdbid in ['1APQ']:
    # for pdbid in ['2KTA', '2KQ8', '2HHI', '2KNQ', '2EZA', '2LGJ', '6V88', '2KCU', '2LT2', '2KPU', '1L1I', '1BLA', '2L29', '2L3B', '2JOZ', '2LTD', '2JR6', '1APQ', '1A5E', '2N9B', '2KB1']:
        test_prmtop = f"data-raw/{pdbid}.prmtop"
        test_inpcrd = f"data-raw/{pdbid}_min.inpcrd"
        # test_dcds = [f"1APQ.gbsa.all/1APQ_tdnr_6002.repl0.dcd"]
        test_dcds = [
            # '1APQ.test.0/1APQ_tdnr_6000.repl0.dcd',
            # '1APQ.test.0/1APQ_tdnr_6001.repl0.dcd',
            # '1APQ.test.0/1APQ_tdnr_6004.repl0.dcd',
            # '1APQ.test.0/1APQ_tdnr_6005.repl0.dcd',
            # '1APQ.test.0/1APQ_tdnr_6007.repl0.dcd',
            '1APQ.test.0/1APQ_tdnr_6008.repl0.dcd',
        ]

        # flexor = BATCorrelations(test_prmtop, test_inpcrd)
        # flexor.compute_dihedrals_from_dcd(test_dcds, 750, 2650)
        # # corr_begin = flexor.compute_correlations()
        # # np.save(f"results/{pdbid}_corr_begin.npy", corr_begin)
        # corr_begin = np.load(f"results/{pdbid}_corr_begin.npy", allow_pickle=True)[0]

        # flexor = BATCorrelations(test_prmtop, test_inpcrd)
        # flexor.compute_dihedrals_from_dcd(test_dcds, 2680)
        # # corr_end = flexor.compute_correlations()
        # # np.save(f"results/{pdbid}_corr_end.npy", corr_end)
        # corr_end = np.load(f"results/{pdbid}_corr_end.npy", allow_pickle=True)[0]

        # # max_begin = np.max(corr_begin)
        # # max_end = np.max(corr_end)
        # # threshold = np.max([max_begin, max_end]) * 0.333

        # # corr_begin = np.abs(corr_begin)
        # # corr_begin[corr_begin < threshold] = 0
        # # corr_begin[corr_begin >= threshold] = 1

        # # corr_end = np.abs(corr_end)
        # # corr_end[corr_end < threshold] = 0
        # # corr_end[corr_end >= threshold] = 1

        # # Assuming `flexor.residue_ids` is your array of values
        # residue_ids = flexor.residue_ids
        # first_index_8 = 27
        # last_index_15 = 55

        # # # difference between the two correlation matrices
        # # union = np.logical_xor(corr_begin, corr_end)
        # # difference = np.logical_and(corr_begin, np.logical_not(corr_end))

        # # Plot both correlation matrices at the same time
        # import matplotlib.pyplot as plt
        # fig, axes = plt.subplots(1, 2, figsize=(20, 5))
        # axes[0].imshow(corr_begin, cmap='coolwarm', vmin=-1, vmax=1)
        # axes[0].set_title(f"Frames 750 to 2650")
        # axes[0].axvline(x=first_index_8, color='blue', linestyle='--', label='First index of 8')
        # axes[0].axhline(y=first_index_8, color='blue', linestyle='--')
        # axes[0].axvline(x=last_index_15, color='blue', linestyle='--', label='First index of 8')
        # axes[0].axhline(y=last_index_15, color='blue', linestyle='--')

        # axes[1].imshow(corr_end, cmap='coolwarm', vmin=-1, vmax=1)
        # axes[1].set_title(f"Frames 2680 to end")
        # axes[1].axvline(x=first_index_8, color='red', linestyle='--', label='Last index of 15')
        # axes[1].axhline(y=first_index_8, color='red', linestyle='--')
        # axes[1].axvline(x=last_index_15, color='blue', linestyle='--', label='First index of 8')
        # axes[1].axhline(y=last_index_15, color='blue', linestyle='--')

        # plt.tight_layout()
        # plt.show()


        # exit()

        # import matplotlib.pyplot as plt
        # # Assuming `flexor.dihedral_types` is a list of dihedral types
        # dihedral_types = flexor.dihedral_types

        # # plot all matrices
        # fig, axes = plt.subplots(1, len(corr), figsize=(20, 5))
        # for i, c in enumerate(corr):
        #     ax = axes[i]
        #     ax.imshow(c, cmap='coolwarm', vmin=-1, vmax=1)
        #     ax.set_title(f"Matrix {i+1}")
        #     ax.set_xticks(range(len(dihedral_types)))
        #     ax.set_yticks(range(len(dihedral_types)))
        #     ax.set_xticklabels(dihedral_types, rotation=90)
        #     ax.set_yticklabels(dihedral_types)

        # # # # Plot the heatmap
        # # # plt.figure(figsize=(10, 8))
        # # # plt.imshow(np.abs(corr[1]), cmap='coolwarm', vmin=0, vmax=1)
        # # # plt.colorbar()

        # # # # Set the dihedral types as labels on the x and y axes
        # # # plt.xticks(ticks=np.arange(len(dihedral_types)), labels=dihedral_types, rotation=90)
        # # # plt.yticks(ticks=np.arange(len(dihedral_types)), labels=dihedral_types)

        # # # # Add labels and title
        # # # plt.xlabel("Dihedral Types")
        # # # plt.ylabel("Dihedral Types")
        # # # plt.title("Heatmap of Corrcorrelations")

        # # plt.tight_layout()
        # # plt.show()

        # # # plot the stdev for the mean of the correlation matrices
        # # stdev = np.max(np.abs(corr), axis=0)
        # # plt.figure(figsize=(10, 8))
        # # plt.imshow(stdev, cmap='coolwarm')
        # # plt.colorbar()
        # # plt.xticks(ticks=np.arange(len(dihedral_types)), labels=dihedral_types, rotation=90)
        # # plt.yticks(ticks=np.arange(len(dihedral_types)), labels=dihedral_types)
        # # plt.xlabel("Dihedral Types")
        # # plt.ylabel("Dihedral Types")
        # # plt.title("Heatmap of Correlation Standard Deviations")
        # plt.tight_layout()
        # plt.show()

        # # abs_corr = np.abs(corr)
        # # mean_corr = np.mean(abs_corr, axis=0)
        # # max_corr = np.max(abs_corr, axis=0)

        # # mean = np.mean(np.abs(corr), axis=0)
        # # flexor.validate_covariance_matrix(mean)
        # # clusters = flexor.cluster_dihedrals(mean)

        # # alpha=50
        # # beta=40
        # # gamma=600
        # # correlation_matrix = np.mean(np.abs(corr), axis=0)

        # # collapsed_set = greedy_collapse(correlation_matrix, alpha, gamma)
        # # blocks = greedy_block(correlation_matrix, collapsed_set, beta)

        flexor = BATCorrelations(test_prmtop, test_inpcrd)
        flexor.compute_dihedrals_from_dcd(test_dcds, 750, 2650)

        # corr = flexor.compute_correlations()
        # np.save(f"results/{pdbid}_corr.npy", corr)
        corr = np.load(f"results/{pdbid}_corr.npy", allow_pickle=True)[0]
        corr = np.abs(corr)

        # # corr_null = []
        # # for _ in range(10):
        # #     corr_null.append(flexor.compute_correlations(True)) # [10, 1, 112, 112]
        # # np.save(f"results/{pdbid}_corr_null.npy", corr_null)
        # corr_null = np.load(f"results/{pdbid}_corr_null.npy", allow_pickle=True)
        # corr_null = corr_null.reshape(10, 112, 112)
        # corr_null = np.abs(corr_null)
        # # corr_null = corr_null[0]
        # corr_null = np.mean(corr_null, axis=0)

        # # fill the diagnol with 0s
        # np.fill_diagonal(corr_null, 0)

        # threshold = flexor.extract_threshold(corr_null)
        threshold = 0.09
        print(f"Threshold: {threshold}")


        # # plot corr and corr_null on the same figure
        # import matplotlib.pyplot as plt
        # fig, axes = plt.subplots(1, 3, figsize=(20, 5))
        # axes[0].imshow(corr, cmap='Reds', vmin=0, vmax=1)
        # axes[0].set_title(f"Corr")
        # axes[1].imshow(corr_null, cmap='Reds', vmin=0, vmax=1)
        # axes[1].set_title(f"Corr Null")
        # axes[2].imshow(corr>threshold, cmap='Reds', vmin=0, vmax=1)
        # axes[2].set_title(f"Corr Null")
        # plt.tight_layout()
        # plt.show()

        ####################################################################
        import numpy as np
        import matplotlib.pyplot as plt
        import networkx as nx
        import community as community_louvain  # pip install python-louvain

        # Create a NetworkX graph from the adjacency matrix
        adj_matrix = np.abs(corr) > threshold  # Threshold to create edges
        np.fill_diagonal(adj_matrix, 0) # Remove self-loops for community detection

        G = nx.from_numpy_array(adj_matrix)
        partition = community_louvain.best_partition(G, random_state=42) # TODO add seed from simulation
        labels = np.array([partition[i] for i in range(len(partition))])
        num_nodes = len(labels)
        modularity = community_louvain.modularity(partition, G)

        blocks = []
        for i in range(len(np.unique(labels))):
            block = np.where(labels == i)[0]
            blocks.append(block)

        # make sure all blocks were selected
        selected = []
        for b in blocks:
            for d in b:
                assert d not in selected
                selected.append(d)
        assert len(selected) == num_nodes

        print(f"Modularity: {modularity:.4f}")
        # exit()

        # Plot original adjacency matrix
        plt.figure(figsize=(8, 8))
        plt.imshow(adj_matrix, cmap='viridis')
        plt.title("Adjacency Matrix with Louvain Communities")
        
        # Determine community boundaries in original order
        prev_label = labels[0]
        for i, label in enumerate(labels[1:], start=1):
            if label != prev_label:
                prev_label = label
                pos = i - 0.5
                plt.axvline(pos)
                plt.axhline(pos)

        # Colored community labels as bars on top and left
        unique_labels = np.unique(labels)
        color_map = plt.get_cmap('tab10', len(unique_labels))
        colors = [color_map(l) for l in labels]
        for i, color in enumerate(colors):
            plt.gca().add_patch(plt.Rectangle((-0.5, i - 0.5), -0.5, 1, color=color))  # left
            plt.gca().add_patch(plt.Rectangle((i - 0.5, -0.5), 1, -0.5, color=color))  # top

        plt.xlim([-1, num_nodes])
        plt.ylim([num_nodes, -1])
        plt.show()



        intra_corrs = []
        for block in blocks:
            submatrix = corr[np.ix_(block, block)]
            # Take upper triangle without diagonal
            triu_vals = submatrix[np.triu_indices(len(block), k=1)]
            intra_corrs.append(triu_vals.mean() if len(triu_vals) else 0)
        print(f"Intra-block correlations: {intra_corrs}")

        

        num_blocks = len(blocks)
        inter_corr_matrix = np.zeros((num_blocks, num_blocks))

        for i in range(num_blocks):
            for j in range(num_blocks):
                if i == j:
                    submatrix = corr[np.ix_(blocks[i], blocks[i])]
                    triu_vals = submatrix[np.triu_indices(len(blocks[i]), k=1)]
                    inter_corr_matrix[i, i] = triu_vals.mean() if len(triu_vals) else 0
                elif i < j:
                    submatrix = corr[np.ix_(blocks[i], blocks[j])]
                    mean_val = submatrix.mean()
                    inter_corr_matrix[i, j] = mean_val
                    inter_corr_matrix[j, i] = mean_val  # symmetry


        # Plot the inter-block correlation matrix
        plt.figure(figsize=(8, 6))
        plt.imshow(inter_corr_matrix, cmap='Reds', vmin=0, vmax=1)
        plt.colorbar(label='Mean Correlation')
        plt.xticks(ticks=np.arange(num_blocks), labels=[f"Block {i+1}" for i in range(num_blocks)], rotation=90)
        plt.yticks(ticks=np.arange(num_blocks), labels=[f"Block {i+1}" for i in range(num_blocks)])
        plt.title("Inter-Block Correlation Matrix")
        plt.xlabel("Blocks")
        plt.ylabel("Blocks")
        plt.tight_layout()
        plt.show()


        oversample_blocks = []
        roll_block = []
        samples_per_round = []

        # For each block, split the dihedrals into two groups: high and low average correlation (with respect to the others in the same group)
        for block in blocks:
            block_corr = np.mean(corr[block][:, block], axis=0)
            # if len(block) < 2:
            #     continue

            # # histogram of block correlation
            # plt.hist(block_corr)
            # plt.title(f"Block {block}")
            # plt.xlabel("Correlation")
            # plt.ylabel("Frequency")
            # plt.axvline(x=np.mean(block_corr), color='r', linestyle='--', label='Mean')
            # plt.axvline(x=np.percentile(block_corr, 25), color='g', linestyle='--', label='25th Percentile')
            # plt.axvline(x=np.percentile(block_corr, 50), color='g', linestyle='--', label='50th Percentile')
            # plt.axvline(x=np.percentile(block_corr, 75), color='g', linestyle='--', label='75th Percentile')
            # plt.axvline(x=np.percentile(block_corr, 95), color='b', linestyle='--', label='95th Percentile')
            # plt.axvline(x=np.percentile(block_corr, 99), color='y', linestyle='--', label='99th Percentile')
            # plt.legend()
            # plt.show()

            # # K-means clustering for two groups
            # from sklearn.cluster import KMeans
            # kmeans = KMeans(n_clusters=2, random_state=0).fit(block_corr.reshape(-1, 1))
            # labels = kmeans.labels_
            # high_group = block[np.where(labels == 0)[0]]
            # low_group = block[np.where(labels == 1)[0]]

            # 95% percentile is the threshold for high correlation
            threshold = np.percentile(block_corr, 50)
            high_group = block[np.where(block_corr >= threshold)[0]]
            low_group = block[np.where(block_corr < threshold)[0]]

            # high_group = block

            # Print the groups and their average correlations
            print(f"Block {block}:")
            print(f"  High group: {high_group}")
            print(f"  Low group: {low_group}")

            oversample_blocks.append(high_group)
            roll_block += list(low_group)

            # print mean block correlation
            block_corr = np.mean(corr[high_group][:, high_group])
            samples_per_round.append(block_corr * len(high_group))
            print(f"  Block correlation: {block_corr:.4f}")


            # # Print each dihedral in the group, its correlation with the others and its assigned group (high or low correlation)
            # for dihedral in high_group:
            #     dihedral_corr = np.mean(corr[dihedral][:, block])
            #     print(f"    Dihedral {dihedral}: {dihedral_corr:.4f} (High)")


        # Compute mean correlation for uncorrelated dihedrals block
        submatrix = corr[np.ix_(roll_block, roll_block)]
        triu_vals = submatrix[np.triu_indices(len(roll_block), k=1)]
        intra_corrs.append(triu_vals.mean() if len(triu_vals) else 0)
        print(f"Intra-correlation for uncorrelated dihedrals block: {intra_corrs[-1]:.4f}")



        # compute samples per round
        def scale_list(data, lower, upper):
            min_val = min(data)
            max_val = max(data)
            if max_val == min_val:
                return [lower] * len(data)  # Avoid division by zero
            return [lower + (x - min_val) * (upper - lower) / (max_val - min_val) for x in data]

        samples_per_round = scale_list(samples_per_round, 1, 10)
        samples_per_round = [int(round(x)) for x in samples_per_round]
        samples_per_round = samples_per_round + [1]  # Add at the end for roll
        print(f"Samples per round: {samples_per_round}")
        # exit()

        

        # print blocks as a python list, comma separated and all
        print("Oversample blocks:")
        for block in oversample_blocks:
            print('[', end='')
            print(', '.join([str(i) for i in block]), end='')
            print('],')
        print("Roll block:")
        print('[', end='')
        print(', '.join([str(i) for i in roll_block]), end='')
        print('],')

        # 19 20 25 28 33 38 40
        # 10 11 13 15 17 20 21
        # sele resi 10+11+13+15+17+20+21
        # sele resi 33+34+36+38+40+43+44
        ####################################################################

        plt.imshow(np.abs(corr), cmap='Reds', vmin=0, vmax=1)
        plt.show()

        # from matplotlib import pyplot as plt
        # plt.imshow(np.abs(corr), cmap='Reds', vmin=0, vmax=1)
        # plt.show()

        # blocks, collapsed = flexor.dynamic_partitioning(np.abs(corr))
        # blocks.append(list(collapsed))

        # # make sure all dihedrals were selected
        # selected = []
        # for b in blocks:
        #     for d in b:
        #         assert d not in selected
        #         selected.append(d)
        # assert len(selected) == flexor.num_dihedrals

        # write out the blocks
        final_blocks = oversample_blocks + [roll_block]

        for i in range(5):
            print("blocks", flexor.chose_decoy_blocks(final_blocks))

        with open(f"results/{pdbid}_cluster.pdb", 'w') as f:
            for b in final_blocks:
                hetatm = []
                conect = []

                for c in b:
                    aix1 = flexor.atom_indices[c][1]
                    aix2 = flexor.atom_indices[c][2]
                    atom1 = flexor.universe.atoms[aix1]
                    atom2 = flexor.universe.atoms[aix2]
                    dihedral_type = 'A' # flexor.dihedral_types[c]

                    hetatm.append(f"HETATM{aix1+1:5d} {dihedral_type[:4]:4s} {atom1.resname:3s} {atom1.resid:4d}    {atom1.position[0]:8.3f}{atom1.position[1]:8.3f}{atom1.position[2]:8.3f}  1.00  0.00          {atom1.name:2s}")
                    hetatm.append(f"HETATM{aix2+2:5d} {dihedral_type[:4]:4s} {atom2.resname:3s} {atom2.resid:4d}    {atom2.position[0]:8.3f}{atom2.position[1]:8.3f}{atom2.position[2]:8.3f}  1.00  0.00          {atom2.name:2s}")
                    conect.append(f"CONECT{aix1+1:5d}{aix2+2:5d}")
                
                f.write("MODEL\n")
                for line in hetatm:
                    f.write(line + '\n')
                for line in conect:
                    f.write(line + '\n')
                f.write("ENDMDL\n")

        exit()

        # import matplotlib.pyplot as plt
        # import seaborn as sns

        # clusters, co_occurrence_matrix = flexor.consensus_clustering(abs_corr, num_clusters=5)
    
        # # Plot co-occurrence matrix
        # plt.figure(figsize=(8, 6))
        # plt.imshow(co_occurrence_matrix, cmap='coolwarm')
        # plt.show()

        # plt.imshow(mean_corr, cmap='coolwarm', vmin=0, vmax=1)
        # plt.colorbar()
        # plt.show()

        # plt.imshow(max_corr, cmap='coolwarm', vmin=0, vmax=1)
        # plt.colorbar()
        # plt.show()


        # # clusters = flexor.cluster(mean_corr, n_clusters=5)
        # with open(f"results/{pdbid}_cluster.pdb", 'w') as f:
        #     for c in clusters:
        #         hetatm = []
        #         conect = []

        #         for aix1, aix2, dihedral_type in c:
        #             atom1 = flexor.universe.atoms[aix1]
        #             atom2 = flexor.universe.atoms[aix2]

        #             hetatm.append(f"HETATM{aix1+1:5d} {dihedral_type[:4]:4s} {atom1.resname:3s} {atom1.resid:4d}    {atom1.position[0]:8.3f}{atom1.position[1]:8.3f}{atom1.position[2]:8.3f}  1.00  0.00          {atom1.name:2s}")
        #             hetatm.append(f"HETATM{aix2+2:5d} {dihedral_type[:4]:4s} {atom2.resname:3s} {atom2.resid:4d}    {atom2.position[0]:8.3f}{atom2.position[1]:8.3f}{atom2.position[2]:8.3f}  1.00  0.00          {atom2.name:2s}")
        #             conect.append(f"CONECT{aix1+1:5d}{aix2+2:5d}")
            
        #         f.write("MODEL\n")
        #         for line in hetatm:
        #             f.write(line + '\n')
        #         for line in conect:
        #             f.write(line + '\n')
        #         f.write("ENDMDL\n")






# Block intra-correlations:
# Block 0: mean intra-correlation = 0.1220 with 19 variables # BEFORE MAX
# Block 1: mean intra-correlation = 0.1332 with 19 variables
# Block 2: mean intra-correlation = 0.1998 with 19 variables
# Block 3: mean intra-correlation = 0.2332 with 19 variables # MAX
# Block 4: mean intra-correlation = 0.1939 with 19 variables
# Block 5: mean intra-correlation = 0.1073 with 19 variables
# Block 6: mean intra-correlation = 0.1169 with 19 variables # AFTER MAX
# Block 7: mean intra-correlation = 0.1791 with 19 variables
# Block 8: mean intra-correlation = 0.1039 with 5 variables
# Block 9: mean intra-correlation = 0.1026 with 12 variables
# Block 10: mean intra-correlation = 0.1007 with 19 variables
# Block 11: mean intra-correlation = 0.1661 with 19 variables


# redo simulations without omega
# plot dihedral stdev
# why are distant chi angles selected?
# all vs all correlation matrices - islam says that we should oversample blocks that are highly correlated with their neighbours
# does the block order matter?
