from dataclasses import dataclass
from typing import Self

import numpy as np
import networkx as nx
import astropy.stats.circstats as circstats
import scipy.cluster.hierarchy as sch
from concurrent.futures import ThreadPoolExecutor, as_completed
import MDAnalysis as mda
from MDAnalysis.analysis import dihedrals
import scipy.stats as stats
from scipy import linalg

import robo_bindings as rb



DIHEDRAL_SELECTIONS = {
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
} # @TODO asp glu arg lys - protonated

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
                 nonbonded_cutoff: float = 1.2,
                 gbsa: float = 1,
                 verbose: bool = False,
                 include_omega: bool = False,
                 include_chi1: bool = True,
                 include_chi2: bool = True,
                 include_chi3: bool = True,
                 include_chi4: bool = True,
                 include_chi5: bool = True):
        
        super().__init__(name, seed, threads, nofRoundsTillReblock, runType, replicaSwapFreq, fixmanSwapFreq)
        self.setPdbRestartFreq(pdb_restart_freq)
        self.setPrintFreq(write_freq)
        self.setNonbonded(nonbonded_method, nonbonded_cutoff) # ce pula mea e metoda asta
        self.setGBSA(gbsa) # true or false in pula mea
        self.setVerbose(verbose)

        self.tol = 1e-6
        self.seed = seed

        # Create the dihedral selections based on input flags
        self.dihedral_sele = DIHEDRAL_SELECTIONS
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
                
        self.include_omega = include_omega
        self.tol = 1e-6
        self.seed = seed
        self.rng = np.random.default_rng(seed)
        self.universe = mda.Universe(prmtop, inpcrd)

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

        # We split world and sampler creation into two parts
        # 1: Add info about them to the context
        # 2: Actually add them to the robo_bindings context in initialize()
        # This is because inside the cpp context, we store worlds in a vector and adding multiple elements will invalidate the samplers (also a vector)
        self.worlds = list[World]()

        # Load the system into Robosample
        super().loadAmberSystem(prmtop, inpcrd)

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
    
    def _build_graph(self):
        """Precomputes the molecular graph without any cuts."""
        G = nx.Graph()
        for atom in self.universe.atoms:
            G.add_node(atom.index)
        for bond in self.universe.bonds:
            idx1, idx2 = bond[0].index, bond[1].index
            G.add_edge(idx1, idx2)
        return G

    def addCartesianWorld(self, samplesPerRound: int = 1) -> World:

        # Create a list of all bonds in the system with full flexibility
        flex = []
        for b in self.universe.bonds:
            flex.append(rb.BondFlexibility(b.atoms[0].index, b.atoms[1].index, rb.BondMobility.Torsion))

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
    
    def addTorsionalWorld(self, flex, samplesperRound: int = 1):
        w = World(fixmanTorque=True,
                  samplesPerRound=samplesperRound,
                  rootMobility=rb.RootMobility.WELD,
                  flexibilities=flex,
                  useOpenMM=False,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=False,
                  samplers=list[Sampler]())
        self.worlds.append(w)
        return self.worlds[-1]

    # bonds: full (tbd), non-redundant (get_dihedral_atom_indices)

    # def getWorld(self, index: int) -> rb.World:
    #     return super().getWorld(index)
    
    def initialize(self, replicaTemperatures: list[float]):

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

            accept_reject_modes.append(w.samplers[0].acceptRejectMode)
            distort_options.append(w.samplers[0].distortOption)
            distort_args.append(w.samplers[0].distortArgs)
            flow.append(w.samplers[0].flow)
            work.append(0)
            integrators.append(w.samplers[0].integratorType)
            timesteps.append(w.samplers[0].timeStep)
            worldIndexes.append(len(worldIndexes))
            mdsteps.append(w.samplers[0].mdSteps)
            boost_md_steps.append(w.samplers[0].boostMDSteps)

        # Add worlds
        # PyBind11 does not allow this binding to use keyword arguments (i.e flexiblities=flex), so we have to pass all arguments in order
        # This applies to all calls to robo_bindings functions
        for world in self.worlds:
            super().addWorld(world.fixmanTorque, world.samplesPerRound, world.rootMobility, world.flexibilities, world.useOpenMM, world.visual, world.visualizerFrequency)

        # Add samplers
        for i, w in enumerate(self.worlds):
            super().getWorld(i).addSampler(w.samplers[0].samplerName, w.samplers[0].integratorType, w.samplers[0].thermostatName, w.samplers[0].useFixmanPotential)

        # Add replicas and thermodynamic states
        for i, T in enumerate(replicaTemperatures):
            super().addReplica(i)
            super().addThermodynamicState(i, T, accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

        # Initialize the context
        super().Initialize()