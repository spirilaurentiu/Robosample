from dataclasses import dataclass
from typing import Self, Iterable, Tuple, List, Set
from collections import deque

import parmed as pmd
import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
from MDAnalysis.core.universe import Merge
import scipy.stats as stats
from scipy import linalg


import robo_bindings as rb
import protein

# @dataclass
# class AtomDefinition(rb.RoboAtomDefinition):
#     globalIndex : int
#     parentAtomGlobalIndex : int
#     molecule_index : int
#     residueIndex : int
#     atomClassIndex : int
#     chargedAtomTypeIndex : int
#     atomName : str
#     residueName : str
#     atomClassName : str
#     chargedAtomName : str
#     neighborsGlobalIndices : list[int]
#     bondsInvolvedGlobalIndex : list[int]
#     availableBonds : int
#     root : bool
#     atomicNumber : int
#     charge : float
#     mass : float
#     vdwRadiusInNm : float
#     vdwWellDepthInKJ : float
#     x : float
#     y : float
#     z : float

#     def __init__(self, globalIndex : int, molIx : int, atomicNumber : int, charge : float, mass : float, vdw : float, lj : float, resName : str, resIx : int, x : float, y : float, z : float, name : str, root : bool):
#         self.globalIndex = globalIndex
#         moleculeIx = molIx
#         self.atomicNumber = atomicNumber
#         self.charge = charge
#         self.mass = mass
#         self.vdw = vdw
#         self.lj = lj
#         self.resName = resName
#         self.resIx = resIx
#         self.x = x
#         self.y = y
#         self.z = z
#         self.name = name
#         self.root = root
#         super().__init__(globalIndex, molIx, atomicNumber, charge, mass, vdw, lj, resName, resIx, x, y, z, name, root)

# @dataclass
# class Bond(rb.BondLink):
#     parentAtomGlobalIndex : int
#     childAtomGlobalIndex : int
#     bondGlobalIndex : int
#     molecule_index : int
#     ringClosing : bool
#     forceK : float
#     forceEquil : float

#     def __init__(self, parentAtomGlobalIndex : int, childAtomGlobalIndex : int, bondGlobalIndex : int, molecule_index : int, ringClosing : bool, forceK : float, forceEquil : float):
#         self.parentAtomGlobalIndex = parentAtomGlobalIndex
#         self.childAtomGlobalIndex = childAtomGlobalIndex
#         self.bondGlobalIndex = bondGlobalIndex
#         self.molecule_index = molecule_index
#         self.ringClosing = ringClosing
#         self.forceK = forceK
#         self.forceEquil = forceEquil
#         super().__init__(parentAtomGlobalIndex, childAtomGlobalIndex, bondGlobalIndex, molecule_index, ringClosing, forceK, forceEquil)

@dataclass
class AtomClassIndexPair:
    universeIndex: int
    parmIndices: List[int]

@dataclass
class Sampler:
    samplerName: rb.SamplerName
    integratorType: rb.IntegratorType
    thermostatName: rb.ThermostatName
    useFixmanPotential: bool
    timeStep: float
    mdSteps: int
    boostMDSteps: int
    acceptRejectMode: rb.AcceptRejectMode
    distortOption: int
    distortArgs: str
    flow: int

@dataclass
class World:
    fixmanTorque: bool
    samplesPerRound: int
    rootMobility: rb.RootMobility
    flexibilities: list[rb.BondFlexibility]
    useOpenMM: bool
    visual: bool
    visualizerFrequency: float
    isCartesian: bool
    samplers: list[Sampler]

    def addSampler(self,
                   timeStep: float,
                   mdSteps: int,
                   boostMDSteps: int,
                   samplerName: rb.SamplerName = rb.SamplerName.HMC,
                   integratorType: rb.IntegratorType = rb.IntegratorType.VERLET,
                   thermostatName: rb.ThermostatName = rb.ThermostatName.ANDERSEN,
                   acceptRejectMode: rb.AcceptRejectMode = rb.AcceptRejectMode.MetropolisHastings,
                   useFixmanPotential: bool = True,
                   distortOption: int = 0,
                   distortArgs: str = "0",
                   flow: int = 0) -> Self:
        assert len(self.samplers) == 0, "Can only add samplers before initializing the context"

        if self.isCartesian:
            integratorType = rb.IntegratorType.OMMVV
            useFixmanPotential = False

        s = Sampler(samplerName=samplerName,
                    integratorType=integratorType,
                    thermostatName=thermostatName,
                    useFixmanPotential=useFixmanPotential,
                    timeStep=timeStep,
                    mdSteps=mdSteps,
                    boostMDSteps=boostMDSteps,
                    acceptRejectMode=acceptRejectMode,
                    distortOption=distortOption,
                    distortArgs=distortArgs,
                    flow=flow)
        self.samplers.append(s)
        return self
    
class Context(rb.Context):
    def __init__(self,
                 name: str,
                 seed: int,
                 prmtop: str,
                 inpcrd: str,
                 write_freq: int,
                 runType: rb.RunType = rb.RunType.REMC,
                 replicaSwapFreq: int = 1,
                 fixmanSwapFreq: int = 0,
                 pdb_restart_freq: int = 0,
                 threads: int = 0,
                 nofRoundsTillReblock: int = 1,
                 nonbonded_method: int = 0,
                 nonbonded_cutoff: float = 1,
                 gbsa: bool = True,
                 verbose: bool = False,
                 testing: bool = False,
                 include_omega: bool = False,
                 include_chi1: bool = True,
                 include_chi2: bool = True,
                 include_chi3: bool = True,
                 include_chi4: bool = True,
                 include_chi5: bool = True,
                 want_n_terminus_phi_rigid: bool = False,
                 want_c_terminus_psi_rigid: bool = False):
        
        super().__init__(name, seed, threads, nofRoundsTillReblock, runType, replicaSwapFreq, fixmanSwapFreq, testing)
        self.setPdbRestartFreq(pdb_restart_freq)
        self.setPrintFreq(write_freq)
        self.setNonbonded(nonbonded_method, nonbonded_cutoff)
        self.setVerbose(verbose)

        # TODO i think it doesn't work
        if gbsa:
            self.setGBSA(1)
        else:
            self.setGBSA(0)

        self.tol = 1e-6
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        # We split world and sampler creation into two parts
        # 1: Add info about them to the context
        # 2: Actually add them to the robo_bindings context in initialize()
        # This is because inside the cpp context, we store worlds in a vector and adding multiple elements will invalidate the samplers (also a vector)
        self.worlds = list[World]()

        # Read the raw data files
        self.parm = pmd.load_file(prmtop, xyz=inpcrd)
        self.universe = mda.Universe(prmtop, inpcrd)
        self.prmtop = prmtop
        self.inpcrd = inpcrd

        # The unit for charge in the prmtop file is e * 18.2223
        # This conversion is automatically handled by ParmEd when loading the file
        # The resulting charge is in units of elementary charge (e) as requested by DuMM in Robosample
        # We still need to correct floaing point errors and normalize to 4 digits
        for a in self.parm.atoms:
            a.charge = round(a.charge, 4)

        # DuMM atom classes are defined by their atom type (XC, C8, N3 etc), not atom name (N, CA, C, O etc)
        atom_classes = set([a.type for a in self.parm.atoms])
        atom_classes = sorted(atom_classes) # Sort to ensure consistent ordering, set() does not guarantee order
        self.atom_class_indices = {atom_type: i for i, atom_type in enumerate(atom_classes)}

        # DuMM charged atom types are AMBER atom types plus their partial charge
        # DuMMForceFieldSubsystemRep::setBiotypeChargedAtomType - there is 1:1 correspondence between biotype and charged atom type
        charged_atom_types = set([self.get_charged_atom_type_name(a) for a in self.parm.atoms])
        charged_atom_types = sorted(charged_atom_types) # Sort to ensure consistent ordering, set() does not guarantee order
        self.charged_atom_type_indices = {atom_type: i for i, atom_type in enumerate(charged_atom_types)}

        self.root_indices = list[int]()
        self.atoms = list[rb.RoboAtom]()
        self.bond_stretches = list[rb.RoboBondStretch]()
        self.bond_bends = list[rb.RoboBondBend]()
        self.bond_torsions = list[rb.RoboBondTorsion]()

        self.atom_ranges: list[tuple[int, int]] = []
        self.bond_stretch_ranges: list[tuple[int, int]] = []
        self.bond_bend_ranges: list[tuple[int, int]] = []
        self.bond_torsion_ranges: list[tuple[int, int]] = []
        
        # Build global index look-up tables
        self.bond_indices: dict[tuple[int, int], AtomClassIndexPair] = {}
        self.angle_indices: dict[tuple[int, int, int], AtomClassIndexPair] = {}
        self.dihedral_indices: dict[tuple[int, int, int, int], AtomClassIndexPair] = {}
        self.improper_indices: dict[tuple[int, int, int, int], AtomClassIndexPair] = {}
        self.create_index_mappings()
        
        # Map original prmtop atom indices to the order in which atoms are added to the Robosample context
        # This is the order in which atoms are explored via BFS starting from the root atom of each molecule offset by the number of atoms in previous molecules
        self.prmtop_to_global_index = {}

        # Store non-redundant bonds
        self.non_redundant_bonds = list[tuple[int, int]]()

        # Parse all the molecules
        for molecule_index, mda_molecule in enumerate(self.universe.atoms.fragments):

            molecule = protein.Protein(mda_molecule,
                                       self.parm,
                                       include_omega=include_omega,
                                       include_chi1=include_chi1,
                                       include_chi2=include_chi2,
                                       include_chi3=include_chi3,
                                       include_chi4=include_chi4,
                                       include_chi5=include_chi5,
                                       want_n_terminus_phi_rigid=want_n_terminus_phi_rigid,
                                       want_c_terminus_psi_rigid=want_c_terminus_psi_rigid)
            
            # Create new atom, bond, angle and dihedral ranges
            atom_range = rb.IteratorPair()
            atom_range.begin = len(self.atoms)
            atom_range.last = len(self.atoms) + len(mda_molecule.atoms)
            self.atom_ranges.append(atom_range)

            bond_range = rb.IteratorPair()
            bond_range.begin = len(self.bond_stretches)
            bond_range.last = len(self.bond_stretches) + len(mda_molecule.bonds)
            self.bond_stretch_ranges.append(bond_range)

            angle_range = rb.IteratorPair()
            angle_range.begin = len(self.bond_bends)
            angle_range.last = len(self.bond_bends) + len(mda_molecule.angles)
            self.bond_bend_ranges.append(angle_range)

            torsion_range = rb.IteratorPair()
            torsion_range.begin = len(self.bond_torsions)
            torsion_range.last = len(self.bond_torsions) + len(mda_molecule.dihedrals) + len(mda_molecule.impropers)
            self.bond_torsion_ranges.append(torsion_range)

            print('Atom range for molecule ', molecule_index, ':(', self.atom_ranges[-1].begin, ',', self.atom_ranges[-1].last, ']')
            print('Bond stretch range for molecule ', molecule_index, ':(', self.bond_stretch_ranges[-1].begin, ',', self.bond_stretch_ranges[-1].last, ']')
            print('Bond bend range for molecule ', molecule_index, ':(', self.bond_bend_ranges[-1].begin, ',', self.bond_bend_ranges[-1].last, ']')
            print('Bond torsion range for molecule ', molecule_index, ':(', self.bond_torsion_ranges[-1].begin, ',', self.bond_torsion_ranges[-1].last, ']')
            
            # Collect prmtop atom indices of non-redundant bonds
            self.non_redundant_bonds.extend(molecule.non_redundant_bonds)
        
            # Get the atoms' prmtop indices in the order they are explored via BFS
            for i, node_prmtop_index in enumerate(molecule.get_nodes_as_prmtop_indices()):

                # Map prmtop index to global index (BFS order)
                self.prmtop_to_global_index[node_prmtop_index] = len(self.prmtop_to_global_index)

                # First atom is always the root
                is_root = False
                if i == 0:
                    is_root = True
                    self.root_indices.append(self.prmtop_to_global_index[node_prmtop_index])

                # Actual atom definition
                a = self.parm.atoms[node_prmtop_index]
                spec = rb.RoboAtomDefinition()

                spec.global_index = self.prmtop_to_global_index[node_prmtop_index]
                spec.prmtop_index = a.idx
                spec.compound_atom_index = molecule.prmtop_to_compound_atom_index(node_prmtop_index)
                spec.molecule_index = molecule_index
                spec.residue_index = a.residue.idx
                spec.atom_class_name = a.type
                spec.atom_class_index = self.atom_class_indices[a.type]
                spec.charged_atom_type_name = self.get_charged_atom_type_name(a)
                spec.charged_atom_type_index = self.charged_atom_type_indices[self.get_charged_atom_type_name(a)]
                spec.residue_name = a.residue.name # ALA
                spec.neighbors_global_indices = [n.idx for n in a.bond_partners] # TODO make it global indices
                spec.atomic_number = a.atomic_number
                spec.charge_in_e = a.charge # The partial atomic charge of this atom in fractions of an electron
                spec.mass_in_daltons = a.mass # The atomic mass of this atom in daltons
                spec.vdw_radius_nm = a.rmin * 2 / 10 # rmin is actually radius/2. We then convert from Angstroms to nm
                spec.sigma_nm = (a.rmin * 2 / 10) / (2 ** (1/6))  # Convert rmin/2 to sigma in nm
                spec.vdw_well_depth_kj = a.epsilon * 4.184 # Convert from kcal/mol to kJ/mol
                spec.solventRadiusInNm = a.solvent_radius / 10 # Convert from A to nm
                spec.screen = a.screen
                spec.x_nm = a.xx / 10 # Convert from A to nm
                spec.y_nm = a.xy / 10 # Convert from A to nm
                spec.z_nm = a.xz / 10 # Convert from A to nm

                # Notice that residue and atom indices in unique atom name are 1-based, not 0-based
                spec.unique_atom_name = a.residue.name + str(a.residue.idx+1) + '_' + a.name + '_' + str(a.idx+1) # e.g. ALA1_N_4
                if is_root:
                    spec.unique_atom_name += '_ROOT' # e.g. ALA1_N_4_ROOT
                spec.root = is_root

                # print('Adding atom', spec.unique_atom_name, 'with global index', spec.global_index, 'and prmtop index', spec.prmtop_index, 'in molecule', molecule_index)

                self.atoms.append(rb.RoboAtom(spec))

            # # Iterate non-redundant bonds
            # for bond in molecule.non_redundant_bonds:
            #     parent_atom = self.prmtop_to_global_index[bond[0]]
            #     child_atom = self.prmtop_to_global_index[bond[1]]
            #     print('Nnon redundant bond between', parent_atom, 'and', child_atom, 'in molecule', molecule_index)
                    
            # Add bonds
            for parent_prmtop_index, child_prmtop_index, is_ring in molecule.get_bonds():
                pair = self.bond_indices[(parent_prmtop_index, child_prmtop_index)]
                assert(len(pair.parmIndices) == 1), f"Multiple bond parameters found for bond between atoms {parent_prmtop_index} and {child_prmtop_index}"
                bond = self.parm.bonds[pair.parmIndices[0]]

                spec = rb.RoboBondStretchDefinition()
                spec.parentAtomGlobalIndex = self.prmtop_to_global_index[parent_prmtop_index]
                spec.childAtomGlobalIndex = self.prmtop_to_global_index[child_prmtop_index]
                spec.parentCompoundAtomIndex = molecule.prmtop_to_compound_atom_index(parent_prmtop_index)
                spec.childCompoundAtomIndex = molecule.prmtop_to_compound_atom_index(child_prmtop_index)
                spec.bondGlobalIndex = pair.universeIndex
                spec.moleculeIndex = molecule_index
                spec.ringClosing = is_ring
                spec.stiffnessInKJPerNmSq = bond.type.k * 4.184 * (10 ** 2) # Convert from kcal/(mol*A^2) to kJ/(mol*nm^2)
                spec.nominalLengthInNm = bond.type.req / 10.0 # Convert from Angstroms to nm

                self.bond_stretches.append(rb.RoboBondStretch(spec))

            # Find angles in the molecule
            for mol_angle in mda_molecule.angles:
                pair = self.angle_indices[(mol_angle.atoms[0].index, mol_angle.atoms[1].index, mol_angle.atoms[2].index)]
                assert(len(pair.parmIndices) == 1), f"Multiple angle parameters found for angle between atoms {mol_angle.atoms[0].index}, {mol_angle.atoms[1].index}, {mol_angle.atoms[2].index}"
                angle = self.parm.angles[pair.parmIndices[0]]
                
                spec = rb.RoboBondBendDefinition()
                spec.globalIndex1 = self.prmtop_to_global_index[mol_angle.atoms[0].index]
                spec.globalIndex2 = self.prmtop_to_global_index[mol_angle.atoms[1].index]
                spec.globalIndex3 = self.prmtop_to_global_index[mol_angle.atoms[2].index]
                spec.compoundAtomIndex1 = molecule.prmtop_to_compound_atom_index(mol_angle.atoms[0].index)
                spec.compoundAtomIndex2 = molecule.prmtop_to_compound_atom_index(mol_angle.atoms[1].index)
                spec.compoundAtomIndex3 = molecule.prmtop_to_compound_atom_index(mol_angle.atoms[2].index)
                spec.moleculeIndex = molecule_index
                spec.stiffnessInKJPerRadSq = angle.type.k * 4.184
                spec.nominalAngleInDeg = angle.type.theteq
                self.bond_bends.append(rb.RoboBondBend(spec))

            # Find torsions in the molecule
            propers_and_impropers = [(dihedral, True) for dihedral in mda_molecule.impropers] + [(dihedral, False) for dihedral in mda_molecule.dihedrals]
            for mol_torsion, is_improper in propers_and_impropers:
                if is_improper:
                    pair = self.improper_indices[(mol_torsion.atoms[0].index, mol_torsion.atoms[1].index, mol_torsion.atoms[2].index, mol_torsion.atoms[3].index)]
                    assert (len(pair.parmIndices) > 0), f"No improper dihedral parameters found for improper between atoms {mol_torsion.atoms[0].index}, {mol_torsion.atoms[1].index}, {mol_torsion.atoms[2].index}, {mol_torsion.atoms[3].index}"
                else:
                    pair = self.dihedral_indices[(mol_torsion.atoms[0].index, mol_torsion.atoms[1].index, mol_torsion.atoms[2].index, mol_torsion.atoms[3].index)]
                    assert(len(pair.parmIndices) > 0), f"No dihedral parameters found for dihedral between atoms {mol_torsion.atoms[0].index}, {mol_torsion.atoms[1].index}, {mol_torsion.atoms[2].index}, {mol_torsion.atoms[3].index}"

                spec = rb.RoboBondTorsionDefinition()
                spec.globalIndex1 = self.prmtop_to_global_index[mol_torsion.atoms[0].index]
                spec.globalIndex2 = self.prmtop_to_global_index[mol_torsion.atoms[1].index]
                spec.globalIndex3 = self.prmtop_to_global_index[mol_torsion.atoms[2].index]
                spec.globalIndex4 = self.prmtop_to_global_index[mol_torsion.atoms[3].index]
                spec.compoundAtomIndex1 = molecule.prmtop_to_compound_atom_index(mol_torsion.atoms[0].index)
                spec.compoundAtomIndex2 = molecule.prmtop_to_compound_atom_index(mol_torsion.atoms[1].index)
                spec.compoundAtomIndex3 = molecule.prmtop_to_compound_atom_index(mol_torsion.atoms[2].index)
                spec.compoundAtomIndex4 = molecule.prmtop_to_compound_atom_index(mol_torsion.atoms[3].index)
                spec.moleculeIndex = molecule_index
                spec.improper = is_improper

                num_terms = len(pair.parmIndices)
                if (num_terms > 0):
                    torsion = self.parm.dihedrals[pair.parmIndices[0]]
                    spec.ampInKJ_1 = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                    spec.phaseInDegrees_1 = torsion.type.phase
                    spec.periodicity_1 = torsion.type.per
                if (num_terms > 1):
                    torsion = self.parm.dihedrals[pair.parmIndices[1]]
                    spec.ampInKJ_2 = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                    spec.phaseInDegrees_2 = torsion.type.phase
                    spec.periodicity_2 = torsion.type.per
                if (num_terms > 2):
                    torsion = self.parm.dihedrals[pair.parmIndices[2]]
                    spec.ampInKJ_3 = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                    spec.phaseInDegrees_3 = torsion.type.phase
                    spec.periodicity_3 = torsion.type.per
                if (num_terms > 3):
                    torsion = self.parm.dihedrals[pair.parmIndices[3]]
                    spec.ampInKJ_4 = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                    spec.phaseInDegrees_4 = torsion.type.phase
                    spec.periodicity_4 = torsion.type.per
                if (num_terms > 4):
                    torsion = self.parm.dihedrals[pair.parmIndices[4]]
                    spec.ampInKJ_5 = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                    spec.phaseInDegrees_5 = torsion.type.phase
                    spec.periodicity_5 = torsion.type.per

                self.bond_torsions.append(rb.RoboBondTorsion(spec))

    def get_cyclomatic_number(self, G: nx.Graph) -> int:
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        num_components = nx.number_connected_components(G)
        cyclomatic_number = num_edges - num_nodes + num_components
        return cyclomatic_number
    
    def getNonRundantBonds(self) -> list[rb.BondFlexibility]:
        flexiblities = []
        for bond in self.non_redundant_bonds:
            flexiblities.append(rb.BondFlexibility(bond[0], bond[1], rb.BondMobility.Torsion))
        return flexiblities

    def addCartesianWorld(self, samplesPerRound: int = 1) -> World:
        # Create a list of all bonds in the system with full flexibility
        flex = []
        for b in self.universe.bonds:
            bat_parent = self.prmtop_to_global_index[b.atoms[0].index]
            bat_child = self.prmtop_to_global_index[b.atoms[1].index]
            flex.append(rb.BondFlexibility(bat_parent, bat_child, rb.BondMobility.Translation))

        w = World(fixmanTorque=False,
                  samplesPerRound=samplesPerRound,
                  rootMobility=rb.RootMobility.WELD,
                  flexibilities=flex,
                  useOpenMM=True,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=True,
                  samplers=list[Sampler]())
        self.worlds.append(w)
        return self.worlds[-1]
    
    def addTorsionalWorld(self, torsional_bonds: list[rb.BondFlexibility], samplesperRound: int = 1):
        """ # !!!!!!!!!!!!!!!!!!!!! torsional bonds are in prmtop order, we reorder them in bat coordinates inside this function !!!!!!!!!!!!!!!!!!!!! """

        # Reorder torsional bonds to match the global indices used in Robosample
        torsional_bonds_reordered = []

        for b in torsional_bonds:
            assert b.mobility == rb.BondMobility.Torsion, "Only torsional flexibilities are supported in torsional world. Rigid bonds are implicitly defined by the absence of a flexibility."

            # This map holds (a,b) and (b,a)
            id = self.bond_indices.get((b.i, b.j))
            assert id is not None, f"Bond between atoms {b.i} and {b.j} not found in bond indices"
            torsional_bonds_reordered.append(rb.BondFlexibility(self.prmtop_to_global_index[b.i], self.prmtop_to_global_index[b.j], b.mobility))

        w = World(fixmanTorque=True,
                  samplesPerRound=samplesperRound,
                  rootMobility=rb.RootMobility.WELD,
                  flexibilities=torsional_bonds_reordered,
                  useOpenMM=True,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=False,
                  samplers=list[Sampler]())
        self.worlds.append(w)

        return self.worlds[-1]
    
    def initialize(self, replicaTemperatures: list[float]):

        # Load the system into Robosample
        super().loadAmberSystem(self.root_indices,
                                self.atoms,
                                self.bond_stretches,
                                self.bond_bends,
                                self.bond_torsions,
                                self.atom_ranges,
                                self.bond_stretch_ranges,
                                self.bond_bend_ranges,
                                self.bond_torsion_ranges)
        
        # Extract thermodynamics state info from all samplers
        accept_reject_modes = list[rb.AcceptRejectMode]()
        distort_options = list[int]()
        distort_args = list[str]()
        flow = list[int]()
        work = list[int]()
        integrators = []
        timesteps = list[float]()
        worldIndexes = list[int]()
        mdsteps = list[int]()
        boost_md_steps = list[int]()
        
        for w in self.worlds:
            assert len(w.samplers) == 1, "Currently only one sampler per world is supported"
            s = w.samplers[0]

            accept_reject_modes.append(s.acceptRejectMode)
            distort_options.append(s.distortOption)
            distort_args.append(s.distortArgs)
            flow.append(s.flow)
            work.append(0)
            integrators.append(s.integratorType)
            timesteps.append(s.timeStep)
            worldIndexes.append(len(worldIndexes))
            mdsteps.append(s.mdSteps)
            boost_md_steps.append(s.boostMDSteps)

        # Add worlds
        # PyBind11 does not allow this binding to use keyword arguments (i.e flexiblities=flex), so we have to pass all arguments in order
        # This applies to all calls to robo_bindings functions
        for world in self.worlds:
            print("Adding world with the following parameters: ",
                  f"fixmanTorque={world.fixmanTorque}, ",
                  f"samplesPerRound={world.samplesPerRound}, ",
                  f"rootMobility={world.rootMobility}, ",
                  f"flexibilities length={len(world.flexibilities)}, ",
                  f"useOpenMM={world.useOpenMM}, ",
                  f"visual={world.visual}, ",
                  f"visualizerFrequency={world.visualizerFrequency}")
            super().addWorld(world.fixmanTorque, world.samplesPerRound, world.rootMobility, world.flexibilities, world.useOpenMM, world.visual, world.visualizerFrequency)

        print("Added worlds to context.")

        # Add samplers
        for i, w in enumerate(self.worlds):
            s = w.samplers[0]
            super().getWorld(i).addSampler(s.samplerName, s.integratorType, s.thermostatName, s.useFixmanPotential)

        print("Initializing context with the following parameters:")

        # Add replicas and thermodynamic states
        for i, T in enumerate(replicaTemperatures):
            super().addReplica(i)
            super().addThermodynamicState(i, T, accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

        # Initialize the context
        super().Initialize()
        
    def create_index_mappings(self):
        # Fill MD Analysis Universe bond indices
        for i, bond in enumerate(self.universe.bonds):
            a1, a2 = bond.atoms[0].index, bond.atoms[1].index
            for key in ((a1, a2), (a2, a1)):
                self.bond_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndices=[])

        for i, bond in enumerate(self.parm.bonds):
            a1, a2 = bond.atom1.idx, bond.atom2.idx
            for key in ((a1, a2), (a2, a1)):
                assert key in self.bond_indices, f"Bond between atoms {a1} and {a2} not found in universe bonds"
                self.bond_indices[key].parmIndices.append(i)

        # Fill MD Analysis Universe angle indices
        for i, angle in enumerate(self.universe.angles):
            a1, a2, a3 = angle.atoms[0].index, angle.atoms[1].index, angle.atoms[2].index
            for key in ((a1, a2, a3), (a3, a2, a1)):
                self.angle_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndices=[])
                
        for i, angle in enumerate(self.parm.angles):
            a1, a2, a3 = angle.atom1.idx, angle.atom2.idx, angle.atom3.idx
            for key in ((a1, a2, a3), (a3, a2, a1)):
                assert key in self.angle_indices, f"Angle between atoms {a1}, {a2}, {a3} not found in universe angles"
                self.angle_indices[key].parmIndices.append(i)
                
        # Fill MD Analysis Universe dihedral indices
        for i, dihedral in enumerate(self.universe.dihedrals):
            a1, a2, a3, a4 = dihedral.atoms[0].index, dihedral.atoms[1].index, dihedral.atoms[2].index, dihedral.atoms[3].index
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                self.dihedral_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndices=[])
                
        # Fill MD Analysis Universe improper dihedral indices
        for i, dihedral in enumerate(self.universe.impropers):
            a1, a2, a3, a4 = dihedral.atoms[0].index, dihedral.atoms[1].index, dihedral.atoms[2].index, dihedral.atoms[3].index
            # print('Improper in universe:', a1, a2, a3, a4)
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                self.improper_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndices=[])

        for i, dihedral in enumerate(self.parm.dihedrals):
            a1, a2, a3, a4 = dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                # skip if this dihedral is actually an improper
                if key in self.dihedral_indices:
                    self.dihedral_indices[key].parmIndices.append(i)
                elif key in self.improper_indices:
                    self.improper_indices[key].parmIndices.append(i)
                else:
                    assert False, f"Dihedral between atoms {a1}, {a2}, {a3}, {a4} not found in universe dihedrals or impropers"

    def get_charged_atom_type_name(self, a: pmd.Atom) -> str:
        return a.type + ':' + str(len(a.bond_partners)) + ':' + str(a.charge)
    
        # # This is how Molmodel does it in TinkerAmber99.cpp
        # return a.residue.name + '_' + a.name

        # # This will always work
        # return a.name + ':' + str(len(a.bond_partners)) + ':' + str(a.charge) + ':' + str(a.idx)
    