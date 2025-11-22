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

DIHEDRAL_SELECTIONS = {
    'ALA': {
        'chi1': ['N', 'CA', 'CB', 'HB1'],
    },

    # Has a terminal large, resonance-stabilized, planar structure with a diffuse positive charge guanidinium group
    # Excluded dihedrals: ['N', 'CZ', 'NH1', 'HH11'] and ['N', 'CZ', 'NH2', 'HH21']
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

    # Imidazole (HID/HIE/HIP) is aromatic and essentially planar — there’s little biologically relevant puckering (only small out-of-plane distortions that correlate with protonation/H-bonding)
    'HIP': {
        'chi1': ['N', 'CA', 'CB', 'CG'],
        'chi2': ['CA', 'CB', 'CG', 'ND1'],
        'ring_closing_1': ['CB', 'CG', 'ND1', 'CE1'],
    },
    'HIE': {
        'chi1': ['N', 'CA', 'CB', 'CG'],
        'chi2': ['CA', 'CB', 'CG', 'ND1'],
        'ring_closing_1': ['CB', 'CG', 'ND1', 'CE1'],
    },
    # 'HID': {
    #     'chi1': ['N', 'CA', 'CB', 'CG'],
    #     'chi2': ['CA', 'CB', 'CG', 'ND1'],
    # },

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
        'ring_closing_1': ['CB', 'CG', 'CD1', 'CE1'], # Same as TYR
    },
    # There are 5 chi angles. Chi1 and Chi2 are the most important for puckering. Chi5 is equivalent to phi. Chi4 is defined as ring closing
    'PRO': {
        'chi1': ['N', 'CA', 'CB', 'CG'],
        'chi2': ['CA', 'CB', 'CG', 'CD'],
        'chi3': ['CB', 'CG', 'CD', 'N'],
        'ring_closing_1': ['CG', 'CD', 'N', 'CA'],
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
    # Tryptophan's side chain is a bicyclic structure called the indole group, which is composed of two fused rings: benzene and pyrrole
    'TRP': {
        'chi1': ['N', 'CA', 'CB', 'CG'],
        'chi2': ['CA', 'CB', 'CG', 'CD1'],
        'ring_closing_1': ['CB', 'CG', 'CD1', 'NE1'], # Pyrrole ring
        'ring_closing_2': ['CG', 'CD2', 'CE3', 'CZ3'], # Benzene ring
    },
    'TYR': {
        'chi1': ['N', 'CA', 'CB', 'CG'],
        'chi2': ['CA', 'CB', 'CG', 'CD1'],
        'chi3': ['CE1', 'CZ', 'OH', 'HH'],
        'ring_closing_1': ['CB', 'CG', 'CD1', 'CE1'], # Same as PHE
    },
    'VAL': {
        'chi1': ['N', 'CA', 'CB', 'CG1'],
        'chi2.1': ['CA', 'CB', 'CG1', 'HG11'],
        'chi2.2': ['CA', 'CB', 'CG2', 'HG21'],
    }
} # @TODO asp glu arg lys - protonated

# @dataclass
# class AtomDefinition(rb.AtomDefinition):
#     globalIndex : int
#     parentAtomGlobalIndex : int
#     moleculeIndex : int
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
#         self.molIx = molIx
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
#     moleculeIndex : int
#     ringClosing : bool
#     forceK : float
#     forceEquil : float

#     def __init__(self, parentAtomGlobalIndex : int, childAtomGlobalIndex : int, bondGlobalIndex : int, moleculeIndex : int, ringClosing : bool, forceK : float, forceEquil : float):
#         self.parentAtomGlobalIndex = parentAtomGlobalIndex
#         self.childAtomGlobalIndex = childAtomGlobalIndex
#         self.bondGlobalIndex = bondGlobalIndex
#         self.moleculeIndex = moleculeIndex
#         self.ringClosing = ringClosing
#         self.forceK = forceK
#         self.forceEquil = forceEquil
#         super().__init__(parentAtomGlobalIndex, childAtomGlobalIndex, bondGlobalIndex, moleculeIndex, ringClosing, forceK, forceEquil)

@dataclass
class AtomClassIndexPair:
    universeIndex: int
    parmIndex: int

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
                 gbsa: bool = True,
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
        self.setNonbonded(nonbonded_method, nonbonded_cutoff)
        self.setVerbose(verbose)

        if gbsa:
            self.setGBSA(1)
        else:
            self.setGBSA(0)

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

        # Read the raw data files
        self.parm = pmd.load_file(prmtop, inpcrd)

        # The unit for charge in the prmtop file is e * 18.2223
        # This conversion is automatically handled by ParmEd when loading the file
        # The resulting charge is in units of elementary charge (e) as requested by DuMM in Robosample
        # We still need to correct floaing point errors and normalize to 4 digits
        for a in self.parm.atoms:
            a.charge = round(a.charge, 4)

        # DuMM atom classes are defined by their atom type (XC, C8, N3 etc), not atom name (N, CA, C, O etc)
        atom_types = set([a.type for a in self.parm.atoms])
        atom_types = sorted(atom_types) # Sort to ensure consistent ordering, set() does not guarantee order
        atom_types_indices = {atom_type: i for i, atom_type in enumerate(atom_types)}

        # DuMM charged atom types are AMBER atom types plus their partial charge
        # DuMMForceFieldSubsystemRep::setBiotypeChargedAtomType - there is 1:1 correspondence between biotype and charged atom type
        charged_atom_types = set([a.name + ':' + str(len(a.bond_partners)) + ':' + str(a.charge) for a in self.parm.atoms])
        charged_atom_types = sorted(charged_atom_types) # Sort to ensure consistent ordering, set() does not guarantee order
        charged_atom_type_indices = {atom_type: i for i, atom_type in enumerate(charged_atom_types)}

        # Iterate all molecules
        self.dihedral_types = []
        self.atom_indices = []
        self.residue_names = []
        self.residue_ids = []
        self.num_dihedrals = 0
        self.dihedral_values = []

        self.root_indices = list[int]()
        self.atoms = list[rb.Atom]()
        self.bond_stretches = list[rb.BondStretch]()
        self.bond_bends = list[rb.BondBend]()
        self.bond_torsions = list[rb.BondTorsion]()
        
        # Build global index look-up tables
        self.bond_indices: dict[tuple[int, int], AtomClassIndexPair] = {}
        self.angle_indices: dict[tuple[int, int, int], AtomClassIndexPair] = {}
        self.dihedral_indices: dict[tuple[int, int, int, int], AtomClassIndexPair] = {}
        self.improper_indices: dict[tuple[int, int, int, int], AtomClassIndexPair] = {}
        self.create_index_mappings()
            
        # Parse all the molecules
        for moleculeIndex, molecule in enumerate(self.universe.atoms.fragments):
            parent = Merge(molecule)
            G = self.build_molecular_graph(parent)

            # The root is the first N atom of the protein
            root_index = molecule.select_atoms("resid 1 and name N")[0].index
            self.root_indices.append(root_index)

            # Add atoms of this molecule
            for i, atom in enumerate(molecule.atoms):
                a = self.parm.atoms[atom.index]
                spec = rb.AtomDefinition()
                spec.globalIndex = a.idx
                spec.moleculeIndex = moleculeIndex
                spec.residueIndex = a.residue.idx
                spec.atomClassIndex = atom_types_indices[a.type]
                spec.chargedAtomTypeIndex = charged_atom_type_indices[a.name + ':' + str(len(a.bond_partners)) + ':' + str(a.charge)]
                spec.biotypeAtomName = a.name + ':' + str(len(a.bond_partners)) + ':' + str(a.charge)  # AMBER atom name + ':' + valence (number of actual bonds, not typical valence)
                spec.atomClassName = a.type
                spec.chargedAtomName = a.name + ':' + str(len(a.bond_partners)) + ':' + str(a.charge)
                spec.residueName = a.residue.name
                spec.uniqueAtomName = a.residue.name + str(a.residue.idx) + '_' + a.name + '_' + str(a.idx)
                spec.neighborsGlobalIndices = [n.idx for n in a.bond_partners]
                spec.root = root_index == a.idx
                spec.atomicNumber = a.atomic_number
                spec.chargeInE = a.charge # The partial atomic charge of this atom in fractions of an electron
                spec.massInDaltons = a.mass # The atomic mass of this atom in daltons
                spec.vdwRadiusInNm = a.rmin * 2 / 10 # rmin is actually radius/2. We then convert from Angstroms to nm
                spec.sigmaInNm = (a.rmin * 2 / 10) / (2 ** (1/6))  # Convert rmin/2 to sigma in nm
                spec.vdwWellDepthInKJ = a.epsilon * 4.184 # Convert from kcal/mol to kJ/mol
                spec.x = a.xx
                spec.y = a.xy
                spec.z = a.xz
                self.atoms.append(rb.Atom(spec))
                
                print(f"Assigned to atom name {a.name} (global index {a.idx}) atom class index {spec.atomClassIndex} ({spec.atomClassName}) and charged atom class index {spec.chargedAtomTypeIndex} ({spec.chargedAtomName})", flush=True)

            # Find redundant (ring closing) dihedrals and build non-redundant set
            atom_groups, dihedral_types, atom_indices, residue_names, residue_ids, num_dihedrals, ring_closing_dihedrals = self.buildNonRedundantTorsions(self.universe)
            pre_ring_closing_bonds = [(rc.indices[1], rc.indices[2]) for rc in ring_closing_dihedrals]
            self.dihedral_types.extend(dihedral_types)
            self.atom_indices.extend(atom_indices)
            self.residue_names.extend(residue_names)
            self.residue_ids.extend(residue_ids)
            self.num_dihedrals += num_dihedrals

            self.non_redundant_bonds = [(ag.indices[1], ag.indices[2]) for ag in atom_groups]

            cyclomatic_number = self.get_cyclomatic_number(G)
            assert cyclomatic_number == len(pre_ring_closing_bonds), f"Cyclomatic number {cyclomatic_number} doesn't match number of ring closing dihedrals {len(pre_ring_closing_bonds)}"

            # Build bonds for the molecule using BFS to ensure a spanning tree
            # This is a custom function that will return all bonds, including ring closing ones
            # There is logic further down to mark ring closing bonds
            tree_edges, ring_only = self.bfs_and_ring_edges(G, root=0, ring_edges=pre_ring_closing_bonds)
            all_edges = [(parent, child, False) for parent, child in tree_edges] + [(parent, child, True)  for parent, child in ring_only]

            for parent, child, is_ring in all_edges:
                pair = self.bond_indices[(parent, child)]
                bond = self.parm.bonds[pair.parmIndex]
                
                spec = rb.BondStretchDefinition()
                spec.parentAtomGlobalIndex = parent
                spec.childAtomGlobalIndex = child
                spec.bondGlobalIndex = pair.universeIndex
                spec.moleculeIndex = moleculeIndex
                spec.ringClosing = is_ring
                spec.stiffnessInKJPerNmSq = bond.type.k * 4.184 * (10 ** 2) # Convert from kcal/(mol*A^2) to kJ/(mol*nm^2)
                spec.nominalLengthInNm = bond.type.req / 10.0 # Convert from Angstroms to nm
                self.bond_stretches.append(rb.BondStretch(spec))

            # Find angles in the molecule
            for mol_angle in molecule.angles:
                pair = self.angle_indices[(mol_angle.atoms[0].index, mol_angle.atoms[1].index, mol_angle.atoms[2].index)]
                angle = self.parm.angles[pair.parmIndex]
                
                spec = rb.BondBendDefinition()
                spec.globalIndex1 = angle.atom1.idx
                spec.globalIndex2 = angle.atom2.idx
                spec.globalIndex3 = angle.atom3.idx
                spec.stiffnessInKJPerRadSq = angle.type.k * 4.184
                spec.nominalAngleInDeg = angle.type.theteq
                self.bond_bends.append(rb.BondBend(spec))

            # Find torsions in the molecule
            for mol_torsion in molecule.dihedrals:
                pair = self.dihedral_indices[(mol_torsion.atoms[0].index, mol_torsion.atoms[1].index, mol_torsion.atoms[2].index, mol_torsion.atoms[3].index)]
                torsion = self.parm.dihedrals[pair.parmIndex]

                spec = rb.BondTorsionDefinition()
                spec.globalIndex1 = torsion.atom1.idx
                spec.globalIndex2 = torsion.atom2.idx
                spec.globalIndex3 = torsion.atom3.idx
                spec.globalIndex4 = torsion.atom4.idx
                spec.ampInKJ = torsion.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                spec.phaseInDegrees = torsion.type.phase
                spec.periodicity = torsion.type.per
                spec.improper = False
                self.bond_torsions.append(rb.BondTorsion(spec))
                
            for mol_improper in molecule.impropers:
                pair = self.improper_indices[(mol_improper.atoms[0].index, mol_improper.atoms[1].index, mol_improper.atoms[2].index, mol_improper.atoms[3].index)]
                improper = self.parm.dihedrals[pair.parmIndex]

                spec = rb.BondTorsionDefinition()
                spec.globalIndex1 = improper.atom1.idx
                spec.globalIndex2 = improper.atom2.idx
                spec.globalIndex3 = improper.atom3.idx
                spec.globalIndex4 = improper.atom4.idx
                spec.ampInKJ = improper.type.phi_k * 4.184 # Convert from kcal/mol to kJ/mol
                spec.phaseInDegrees = improper.type.phase
                spec.periodicity = improper.type.per
                spec.improper = True
                self.bond_torsions.append(rb.BondTorsion(spec))

            # assert len(self.ordered_bonds) == G.number_of_edges(), f"Number of BFS bonds {len(self.ordered_bonds)} doesn't match number of edges {G.number_of_edges()}"
            # self.ordered_bonds.append(self.ordered_bonds)

        # # Print roots and bonds
        # for r in self.root_indices:
        #     print(f"Root atom index: {r} ({self.universe.atoms[r].resname} {self.universe.atoms[r].resid} {self.universe.atoms[r].name})")
        # for i, bonds in enumerate(self.ordered_bonds):
        #     print(f"Molecule {i} bonds:")
        #     for j, b in enumerate(bonds):
        #         print(f"  Bond {j} between atoms {b.parent} and {b.child}, index {b.index}, ring closing: {b.ring_closing}")

        # exit()


        # We split world and sampler creation into two parts
        # 1: Add info about them to the context
        # 2: Actually add them to the robo_bindings context in initialize()
        # This is because inside the cpp context, we store worlds in a vector and adding multiple elements will invalidate the samplers (also a vector)
        self.worlds = list[World]()

        self.prmtop = prmtop
        self.inpcrd = inpcrd

    def buildNonRedundantTorsions(self, universe):
        """
        Build a list of all non-redundant torsions in the system.
        This includes phi, psi and chi angles.
        It removes one torsion per independent cycle (the cyclomatic number) and any other chemically fixed dihedrals (e.g. omega).
        The specific chemical nature of the torsion you remove doesn't matter for topological redundancy
        if the cyclomatic number is greater than the number of removed omega torsions, then additional torsions need to be removed from the remaining independent cycles that weren't "satisfied" by the omega removals.

        DOESN'T CYX REMOVAL GUARANTEE THAT ALL DISULFIDE CYCLES ARE REMOVED?
        """
        atom_groups = []
        dihedral_types = []
        atom_indices = []
        residue_names = []
        residue_ids = []

        ring_closing_dihedrals = []

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

                assert len(chi_dihedral) == 4, f"Dihedral {chi_name} in residue {res.resname} {res.resid} does not have 4 atoms: {chi_atoms}, found {chi_dihedral}"

                if 'closing' in chi_name:
                    ring_closing_dihedrals.append(chi_dihedral)
                    continue

                atom_groups.append(chi_dihedral)
                dihedral_types.append(chi_name)
                atom_indices.append(chi_dihedral.indices)
                residue_names.append(res.resname)
                residue_ids.append(res.resid)

                # print(f"Found dihedral {chi_name} in residue {res.resname} {res.resid} with atoms {chi_atoms}")

        # Find the disulfide bonds
        disulfide_bonds = set()
        for sg in universe.select_atoms("resname CYX and name SG"):
            for b in sg.bonds:
                if b.atoms[0].name == 'SG' and b.atoms[1].name == 'SG':
                    disulfide_bonds.add(tuple(sorted(b.indices)))
        disulfide_bonds = list(disulfide_bonds)

        # print the atoms in the disulfide bonds
        for bond in disulfide_bonds:
            sg1_atom = universe.atoms[bond[0]]
            sg2_atom = universe.atoms[bond[1]]

            sg_atom_pairs = [(sg1_atom, sg2_atom), (sg2_atom, sg1_atom)]
            for sg_atom_1, sg_atom_2 in sg_atom_pairs:
                sg1_selection = universe.select_atoms(f"resid {sg_atom_1.resid} and name SG")
                sg2_selection = universe.select_atoms(f"resid {sg_atom_2.resid} and name SG")
                cb_selection = universe.select_atoms(f"resid {sg_atom_2.resid} and name CB")
                ca_selection = universe.select_atoms(f"resid {sg_atom_2.resid} and name CA")

                atom_group = mda.AtomGroup([sg1_selection.atoms[0], sg2_selection.atoms[0], cb_selection.atoms[0], ca_selection.atoms[0]])
                atom_groups.append(atom_group)
                dihedral_types.append('chi2')
                atom_indices.append(atom_group.indices)
                residue_names.append('CYX')
                residue_ids.append(sg_atom_2.resid)

            # Add ring closing on the actual disulfide bond
            cb1_selection = universe.select_atoms(f"resid {sg1_atom.resid} and name CB")
            sg1_selection = universe.select_atoms(f"resid {sg1_atom.resid} and name SG")
            sg2_selection = universe.select_atoms(f"resid {sg2_atom.resid} and name SG")
            cb2_selection = universe.select_atoms(f"resid {sg2_atom.resid} and name CB")

            atom_group = mda.AtomGroup([cb1_selection.atoms[0], sg1_selection.atoms[0], sg2_selection.atoms[0], cb2_selection.atoms[0]])
            ring_closing_dihedrals.append(atom_group)

        num_dihedrals = len(atom_indices)

        return atom_groups, dihedral_types, atom_indices, residue_names, residue_ids, num_dihedrals, ring_closing_dihedrals
    
    def build_molecular_graph(self, parent: mda.Universe) -> nx.Graph:
        """Precomputes the molecular graph without any cuts."""
        G = nx.Graph()
        for atom in parent.atoms:
            G.add_node(atom.index)
        for bond in parent.bonds:
            idx1, idx2 = bond[0].index, bond[1].index
            G.add_edge(idx1, idx2)
        return G
    
    def get_cyclomatic_number(self, G: nx.Graph) -> int:
        num_nodes = G.number_of_nodes()
        num_edges = G.number_of_edges()
        num_components = nx.number_connected_components(G)
        cyclomatic_number = num_edges - num_nodes + num_components
        return cyclomatic_number
    
    def bfs_and_ring_edges(self, G: nx.Graph,root: int,ring_edges: Iterable[Tuple[int, int]]) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
        """
        Return BFS-tree edges and ring-closing edges separately.

        Parameters
        ----------
        G : nx.Graph
            Molecular graph (undirected).
        root : int
            Starting heavy atom.
        ring_edges : iterable of (int,int)
            Undirected ring-closing bonds.

        Returns
        -------
        bfs_edges : list of (int,int)
            Edges discovered by BFS, excluding ring bonds.
            For each (u,v), u is guaranteed to have been visited first.
        ring_only : list of (int,int)
            Unique ring-closing edges (undirected, sorted tuple).
        """
        ring_set: Set[Tuple[int, int]] = {tuple(sorted(e)) for e in ring_edges}
        bfs_edges: List[Tuple[int, int]] = []
        seen_nodes = {root}
        q = deque([root])

        while q:
            src = q.popleft()
            for nbr in G.neighbors(src):
                edge_key = tuple(sorted((src, nbr)))
                if edge_key in ring_set:
                    # ring bonds are collected separately
                    continue
                if nbr not in seen_nodes:
                    seen_nodes.add(nbr)
                    q.append(nbr)
                    bfs_edges.append((src, nbr))

        # return the ring bonds as a unique, order-independent list
        ring_only = list(ring_set)
        return bfs_edges, ring_only
    
    def get_bond(self, molecule, parent: int, child: int):
        """Return (index, bond) for the bond connecting the two atoms, or (-1, None) if not found."""
        for i, b in enumerate(self.molecule.bonds):
            a0, a1 = b.atoms
            if (a0.index == parent and a1.index == child) or (a0.index == child and a1.index == parent):
                return i, b
        return -1, None

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
                  useOpenMM=True,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=False,
                  samplers=list[Sampler]())
        self.worlds.append(w)
        return self.worlds[-1]
    
    def getNonRundantBonds(self) -> list[Tuple[int, int, rb.BondMobility]]:
        flexibilities = []
        for bond in self.non_redundant_bonds:
            flexibilities.append(rb.BondFlexibility(bond[0], bond[1], rb.BondMobility.Torsion))
        return flexibilities

    def initialize(self, replicaTemperatures: list[float]):

        # Load the system into Robosample
        super().loadAmberSystem(self.root_indices, self.atoms, self.bond_stretches, self.bond_bends, self.bond_torsions)

        print('Loaded system into context.')

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
        # assert len(self.universe.bonds) == len(self.parm.bonds), f"Number of bonds in universe ({len(self.universe.bonds)}) does not match number of bonds in parm ({len(self.parm.bonds)})"
        # assert len(self.universe.angles) == len(self.parm.angles), f"Number of angles in universe ({len(self.universe.angles)}) does not match number of angles in parm ({len(self.parm.angles)})"
        
        
        # print(len({tuple(sorted((dih.atom1.idx, dih.atom2.idx, dih.atom3.idx, dih.atom4.idx))) for dih in self.parm.dihedrals}))
        # print(len(self.parm.dihedrals), len(self.parm.impropers))
        # print(len(self.universe.dihedrals), len(self.universe.impropers))
        
        # exit()
        
        # assert len(self.universe.dihedrals) + len(self.universe.impropers) == len(self.parm.dihedrals), f"Number of dihedrals in universe ({len(self.universe.dihedrals) + len(self.universe.impropers)}) does not match number of dihedrals in parm ({len(self.parm.dihedrals)})"
        
        # Fill MD Analysis Universe bond indices
        for i, bond in enumerate(self.universe.bonds):
            a1, a2 = bond.atoms[0].index, bond.atoms[1].index
            for key in ((a1, a2), (a2, a1)):
                self.bond_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndex=-1)

        for i, bond in enumerate(self.parm.bonds):
            a1, a2 = bond.atom1.idx, bond.atom2.idx
            for key in ((a1, a2), (a2, a1)):
                assert key in self.bond_indices, f"Bond between atoms {a1} and {a2} not found in universe bonds"
                self.bond_indices[key].parmIndex = i
                    
        # Fill MD Analysis Universe angle indices
        for i, angle in enumerate(self.universe.angles):
            a1, a2, a3 = angle.atoms[0].index, angle.atoms[1].index, angle.atoms[2].index
            for key in ((a1, a2, a3), (a3, a2, a1)):
                self.angle_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndex=-1)
                
        for i, angle in enumerate(self.parm.angles):
            a1, a2, a3 = angle.atom1.idx, angle.atom2.idx, angle.atom3.idx
            for key in ((a1, a2, a3), (a3, a2, a1)):
                assert key in self.angle_indices, f"Angle between atoms {a1}, {a2}, {a3} not found in universe angles"
                self.angle_indices[key].parmIndex = i
                    
        # Fill MD Analysis Universe dihedral indices
        for i, dihedral in enumerate(self.universe.dihedrals):
            a1, a2, a3, a4 = dihedral.atoms[0].index, dihedral.atoms[1].index, dihedral.atoms[2].index, dihedral.atoms[3].index
            print(f"MDanalysis defined dihedral ({a1}, {a2}, {a3}, {a4}) or ({a4}, {a3}, {a2}, {a1}) at index {i}", flush=True)
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                self.dihedral_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndex=-1)
                
        for i, dihedral in enumerate(self.parm.dihedrals):
            a1, a2, a3, a4 = dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx
            print(f"Parmed defined dihedral ({a1}, {a2}, {a3}, {a4}) or ({a4}, {a3}, {a2}, {a1}) at index {i}", flush=True)
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                # skip if this dihedral is actually an improper
                if key in self.dihedral_indices:
                    self.dihedral_indices[key].parmIndex = i
                
        # Fill MD Analysis Universe improper dihedral indices
        for i, dihedral in enumerate(self.universe.impropers):
            a1, a2, a3, a4 = dihedral.atoms[0].index, dihedral.atoms[1].index, dihedral.atoms[2].index, dihedral.atoms[3].index
            print(f"MDanalysis defined improper dihedral ({a1}, {a2}, {a3}, {a4}) or ({a4}, {a3}, {a2}, {a1}) at index {i}", flush=True)
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                self.improper_indices[key] = AtomClassIndexPair(universeIndex=i, parmIndex=-1)
                
        for i, dihedral in enumerate(self.parm.impropers):
            a1, a2, a3, a4 = dihedral.atom1.idx, dihedral.atom2.idx, dihedral.atom3.idx, dihedral.atom4.idx
            print(f"Parmed defined improper dihedral ({a1}, {a2}, {a3}, {a4}) or ({a4}, {a3}, {a2}, {a1}) at index {i}", flush=True)
            for key in ((a1, a2, a3, a4), (a4, a3, a2, a1)):
                assert key in self.improper_indices, f"Improper dihedral between atoms {a1}, {a2}, {a3}, {a4} not found in universe dihedrals"
                self.improper_indices[key].parmIndex = i
                    
