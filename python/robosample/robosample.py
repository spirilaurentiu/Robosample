from dataclasses import dataclass
from typing import Self, Iterable, Tuple, List, Set, Dict
from collections import deque, defaultdict
from enum import IntEnum, unique
import mdtraj as md
from contextlib import contextmanager
from matplotlib import units
import parmed as pmd
from parmed import unit as u
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
from molecule_prototype import MoleculePrototype
import prmtop_reader

@unique
class NonbondedMethod(IntEnum):
    """
    Nonbonded interaction treatment for OpenMM NonbondedForce.

    This enum is a high-level, typed, IDE-friendly facade over
    ``rb.NonbondedMethod`` (pybind11 binding of
    ``OpenMM::NonbondedForce::NonbondedMethod``).

    Members map 1:1 to the underlying OpenMM values and can be
    passed transparently to C++ bindings.
    """

    NoCutoff = rb.NonbondedMethod.NoCutoff
    """
    No cutoff is applied to nonbonded interactions.
    The full set of N^2 interactions is computed exactly.
    This necessarily means that periodic boundary conditions cannot be used.
    """

    CutoffNonPeriodic = rb.NonbondedMethod.CutoffNonPeriodic
    """
    Interactions beyond the cutoff distance are ignored.
    Coulomb interactions closer than the cutoff distance are modified using the reaction field method.
    """

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
    SIGMA_SCALE = 2**(-1./6.)

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
                 use_gbsa_obc2: bool = True,
                 gbsa_solvent_dielectric: float = 78.5,
                 gbsa_solute_dielectric: float = 1.0,
                 nonbonded_method: NonbondedMethod = NonbondedMethod.CutoffNonPeriodic,
                 nonbonded_cutoff_in_nm: float = 1.2,
                 verbose: bool = False,
                 testing: bool = False,
                 rigid_protein_phi: bool = False,
                 rigid_protein_psi: bool = False,
                 rigid_protein_omega: bool = True,
                 rigid_protein_chi1: bool = False,
                 rigid_protein_chi2: bool = False,
                 rigid_protein_chi3: bool = False,
                 rigid_protein_chi4: bool = False,
                 rigid_protein_chi5: bool = False):
        
        super().__init__(name, seed, threads, nofRoundsTillReblock, runType, replicaSwapFreq, fixmanSwapFreq, testing)
        self.setPdbRestartFreq(pdb_restart_freq)
        self.setPrintFreq(write_freq)
        self.setNonbonded(rb.NonbondedMethod.CutoffNonPeriodic, nonbonded_cutoff_in_nm)
        self.setVerbose(verbose)
        self.setGBSAOptions(use_gbsa_obc2, gbsa_solvent_dielectric, gbsa_solute_dielectric)

        self.tol = 1e-3
        self.seed = seed
        self.rng = np.random.default_rng(seed)

        # We split world and sampler creation into two parts
        # 1: Add info about them to the context
        # 2: Actually add them to the robo_bindings context in initialize()
        # This is because inside the cpp context, we store worlds in a vector and adding multiple elements will invalidate the samplers (also a vector)
        self.worlds = list[World]()

        # Read the raw data files
        self.parm = pmd.load_file(prmtop, xyz=inpcrd)
        self.prmtop = prmtop
        self.inpcrd = inpcrd
        self.num_types = self.parm.pointers['NTYPES']

        # parmed does nasty rounding when loading and loses some precision that adds up to a few kj
        # prmtop files hold more decimal places than can be stored via Python float64 (IEEE 754 double) has ~16 decimal digits of precision
        # prmtop holds more that 16, so this function will lose a few digits (fewer than parmed)
        parm_file = prmtop_reader.parse_prmtop_numpy(prmtop)
        parm_data = parm_file['raw_data']

        # Nonbonded fix (NBFIX) is a technique that replaces standard Lennard-Jones (LJ) interaction parameters (epsilon and sigma) between specific atom pairs
        # This overrides default combination rules to fix overbinding artifacts, particularly between cations/anions and protein/lipid functional groups
        # It is commonly used in CHARMM force fields to improve hydration and binding accuracy
        self.has_nbfix = prmtop_reader.has_nbfix_fast(parm_data['NONBONDED_PARM_INDEX'], self.num_types, parm_data['LENNARD_JONES_ACOEF'], parm_data['LENNARD_JONES_BCOEF'])

        ene_conv = pmd.unit.kilocalories_per_mole.conversion_factor_to(pmd.unit.kilojoules_per_mole)
        length_conv = pmd.unit.angstroms.conversion_factor_to(pmd.unit.nanometers)
        afac = np.sqrt(ene_conv) * length_conv**6
        bfac = ene_conv * length_conv**6

        self.acoef = [0 for _ in range(self.num_types * self.num_types)]
        self.bcoef = [0 for _ in range(self.num_types * self.num_types)]

        for i in range(self.num_types):
            for j in range(self.num_types):
                idx = parm_data['NONBONDED_PARM_INDEX'][i * self.num_types + j] - 1
                if idx < 0:
                    raise ValueError(f"Invalid nonbonded index for atom types {i} and {j}")
                self.acoef[i * self.num_types + j] = np.sqrt(parm_data['LENNARD_JONES_ACOEF'][idx]) * afac
                self.bcoef[i * self.num_types + j] = parm_data['LENNARD_JONES_BCOEF'][idx] * bfac

        # DuMM atom classes are defined by their atom type (XC, C8, N3 etc), not atom name (N, CA, C, O etc)
        atom_classes = self.generate_synthetic_atom_classes()
        unique_atom_classes = set(atom_classes.values())
        unique_atom_classes = sorted(unique_atom_classes)
        self.atom_class_indices = {key: i+1 for i, key in enumerate(unique_atom_classes)}

        # DuMM charged atom types are AMBER atom types plus their partial charge
        # DuMMForceFieldSubsystemRep::setBiotypeChargedAtomType - there is 1:1 correspondence between biotype and charged atom type
        charged_atom_types = set([atom_classes[a] + ':' + str(a.charge) for a in self.parm.atoms])
        charged_atom_types = sorted(charged_atom_types) # Sort to ensure consistent ordering, set() does not guarantee order
        self.charged_atom_type_indices = {atom_type: i for i, atom_type in enumerate(charged_atom_types)}

        self.root_indices = list[int]()
        self.topology_ranges: list[rb.TopologyRange] = []
        self.atoms: list[rb.RoboAtom] = []
        self.bond_stretches: list[rb.RoboBond] = []
        self.bond_bends: list[rb.RoboAngle] = []
        self.proper_periodic_torsions: list[rb.RoboPeriodicTorsion] = []
        self.improper_harmonic_torsions: list[rb.RoboHarmonicImproperTorsion] = []
        self.z_matrix = list[tuple[int, int, int, int]]()
        self.num_atom_offset = 0
        self.num_residues_offset = 0

        # Map original prmtop atom indices to the order in which atoms are added to the Robosample context
        # This is the order in which atoms are explored via BFS starting from the root atom of each molecule offset by the number of atoms in previous molecules
        self.prmtop_to_global_index = {}

        # Store non-redundant bonds
        self.standard_dihedral_bonds = list[list[tuple[int, int]]]()
        self.macrocycle_dihedral_bonds = list[list[list[tuple[int, int]]]]()

        # Split system into unique molecule prototypes and their occurrences
        # Each entry: (prototype_structure, instance_indices)
        # Example: parm_prototypes[0] holds the first molecule type found in the system
        # parm_prototypes[0][0] is the parmed Structure of that molecule type
        # parm_prototypes[0][1] is a list of all instance indices of that molecule type in the full system eg [0, 3, 5] if that molecule occurs at those indices
        parm_prototypes = self.parm.split()

        # Parsing atom, bond, angle and torsion definitions is expensive, so we cache for each prototype type
        molecule_prototypes = [MoleculePrototype(mol_struct) for mol_struct, _ in parm_prototypes]

        # [(instance_index, prototype_index)]
        molecules: list[tuple[int, int]] = []
        for protoptype_index in range(len(parm_prototypes)):
            for instance_index in parm_prototypes[protoptype_index][1]:
                molecules.append((instance_index, protoptype_index))
        molecules.sort(key=lambda x: x[0])

        # Map each molecule instance index to its prototype index
        # instance_to_prototype[i] = j -> molecule instance i in the full system is of prototype parm_prototypes[j] and its structure is parm_prototypes[j][0]
        # Molecule instanc_index is of type prototype_index
        # instance_to_prototype: dict[int, int] = {}
        for instance_index, prototype_index in molecules:
            # Record begin and end indices for this molecule's topology definitions
            with self.record_topology(molecule_prototypes[prototype_index]) as topo_range:
                # Add atoms
                for atom in molecule_prototypes[prototype_index].atom_params:
                    prmtop_index = atom.local_index + self.num_atom_offset
                    a = self.parm.atoms[prmtop_index]
                    if a.idx != prmtop_index:
                        raise ValueError(f"Atom index mismatch: expected {prmtop_index}, got {a.idx}")
                    if prmtop_index in self.prmtop_to_global_index:
                        raise ValueError(f"Duplicate prmtop index found in mapping: {prmtop_index}")
                        
                    # Map prmtop index to global index (BFS order)
                    global_index = len(self.prmtop_to_global_index)
                    self.prmtop_to_global_index[prmtop_index] = global_index
                        
                    # Get atom class and charged atom type
                    atom_class_name = atom_classes[a]
                    charged_atom_type_name = atom_classes[a] + ':' + str(a.charge)

                    # Residue and atom indices in unique atom name are 1-based, not 0-based
                    unique_atom_name = a.residue.name + str(a.residue.idx+1) + '_' + a.name + '_' + str(prmtop_index+1) # e.g. ALA1_N_4

                    # First atom is always the root and holds magic properties
                    if atom.root:
                        unique_atom_name += '_ROOT' # e.g. ALA1_N_4_ROOT
                        self.root_indices.append(global_index)

                    # Get Lennard-Jones parameters for this atom
                    # Nonbonded indices are 1-based
                    idx = (a.nb_idx - 1) * self.num_types + (a.nb_idx - 1)
                    if idx < 0:
                        raise ValueError(f"Invalid nonbonded index for atom {a.idx}")
                    nb_parm_idx = parm_data['NONBONDED_PARM_INDEX'][idx] - 1
                    acoef = parm_data['LENNARD_JONES_ACOEF'][nb_parm_idx]
                    bcoef = parm_data['LENNARD_JONES_BCOEF'][nb_parm_idx]

                    # Parameters may be undefined for some atoms, typicall hydrogen atoms
                    if acoef != 0.0 and bcoef != 0.0:
                        r_min = (2*acoef/bcoef)**(1/6.0)
                        epsilon = 0.25*bcoef*bcoef/acoef
                    else:
                        r_min = 1.0
                        epsilon = 0.0

                    # Convert to sigma and epsilon
                    lengthConversionFactor = pmd.unit.angstrom.conversion_factor_to(pmd.unit.nanometer)
                    energyConversionFactor = pmd.unit.kilocalorie_per_mole.conversion_factor_to(pmd.unit.kilojoule_per_mole)

                    rVdw = r_min / 2.0 * lengthConversionFactor
                    sigma = rVdw * 2.0 * self.SIGMA_SCALE
                    epsilon = epsilon * energyConversionFactor

                    # Define the atom
                    self.atoms.append(
                        rb.RoboAtom(
                            identity=rb.RoboAtomIdentity(
                                unique_name=unique_atom_name,
                                residue_name=a.residue.name,
                                atom_class_name=atom_class_name,
                                charged_atom_type_name=charged_atom_type_name,
                                global_index=global_index,
                                prmtop_index=prmtop_index,
                                molecule_index=instance_index,
                                residue_index=a.residue.idx,
                                nonbonded_index=a.nb_idx-1,
                                compound_atom_index=atom.compound_atom_index,
                                atom_class_index=self.atom_class_indices[atom_class_name],
                                charged_atom_type_index=self.charged_atom_type_indices[charged_atom_type_name]
                            ),
                            element_info=rb.RoboAtomElement(
                                element_name=atom.element_name,
                                element_symbol=atom.element_symbol,
                                atomic_number=atom.atomic_number
                            ),
                            physics=rb.RoboAtomPhysics(
                                charge_e=parm_data['CHARGE'][prmtop_index],
                                mass_daltons=parm_data['MASS'][prmtop_index],
                                vdw_radius_nm=atom.vdw_radius_nm,
                                vdw_well_depth_kj=epsilon,
                                sigma_nm=sigma,
                                solvent_radius_nm=parm_data['RADII'][prmtop_index]/10,
                                screen=parm_data['SCREEN'][prmtop_index],
                            ),
                            connectivity=rb.RoboAtomConnectivity(
                                neighbors_global_indices=[n.idx for n in a.bond_partners],
                                root=atom.root
                            ),
                            position=[a.xx / 10, a.xy / 10, a.xz / 10]
                        )
                    )

                # Add bonds
                for bond in molecule_prototypes[prototype_index].bond_params:
                    prmtop_indices = [l+self.num_atom_offset for l in bond.local_indices]
                    for prmtop_index in prmtop_indices:
                        atom = self.parm.atoms[prmtop_index]
                        if atom.idx != prmtop_index:
                            raise ValueError(f"Bond atom index mismatch: expected {prmtop_index}, got {atom.idx}")
                        
                    self.bond_stretches.append(
                        rb.RoboBond(
                            global_indices=tuple(self.prmtop_to_global_index[p] for p in prmtop_indices),
                            compound_atom_indices=bond.compound_atom_indices,
                            stiffness_in_kj_per_nm_sq=bond.stiffness_in_kj_per_nm_sq,
                            nominal_length_in_nm=bond.nominal_length_in_nm,
                            molecule_index=instance_index,
                            ring_closing=bond.is_ring_closing,
                            dihedral_type=bond.dihedral_type
                        )
                    )

                # import mdtraj as md
                # traj = md.load(self.inpcrd, top=self.prmtop)
                # ss = md.compute_dssp(traj)
                


                # # # # # Convert parmed to mdtraj for DSSP secondary structure assignment to determine which dihedral bonds are rotatable
                # # # # top = molecule_prototypes[prototype_index].molecule.topology
                # # # # xyz = molecule_prototypes[prototype_index].molecule.coordinates / 10.0

                # # # # traj = md.Trajectory(xyz=[xyz], topology=top)
                # # # # ss = md.compute_dssp(traj)

                # Add bonds which are the middle bond of a standard dihedral
                # We don't convert to global indices and must keep original prmtop ones
                for (p1, p2, dihedral_type, resid) in molecule_prototypes[prototype_index].bonds:
                    atom1_prmtop = p1+self.num_atom_offset
                    if self.parm.atoms[atom1_prmtop].idx != atom1_prmtop:
                        raise ValueError(f"Standard dihedral bond atom index mismatch: expected {atom1_prmtop}, got {self.parm.atoms[atom1_prmtop].idx}")
                    atom2_prmtop = p2+self.num_atom_offset
                    if self.parm.atoms[atom2_prmtop].idx != atom2_prmtop:
                        raise ValueError(f"Standard dihedral bond atom index mismatch: expected {atom2_prmtop}, got {self.parm.atoms[atom2_prmtop].idx}")
                    
                    # Mandatory
                    if 'ring' in dihedral_type:
                        continue

                    # if ss[0][resid] != 'C':
                    #     continue

                    if 'pro-' in dihedral_type:
                        continue

                    # Optional
                    if dihedral_type == 'protein-phi' and not rigid_protein_phi:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-psi' and not rigid_protein_psi:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-omega' and not rigid_protein_omega:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-chi1' and not rigid_protein_chi1:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-chi2' and not rigid_protein_chi2:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-chi3' and not rigid_protein_chi3:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-chi4' and not rigid_protein_chi4:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if dihedral_type == 'protein-chi5' and not rigid_protein_chi5:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                    if 'mandatory' in dihedral_type:
                        self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))

                # # Add macrocycle dihedral bonds
                # for macrocycle in molecule_prototypes[prototype_index].macrocycle_bonds:
                #     macrocycle_bonds = list[list[tuple[int, int]]]()

                #     for (p1, p2, dihedral_type, resid) in macrocycle:
                #         atom1_prmtop = p1+self.num_atom_offset
                #         if self.parm.atoms[atom1_prmtop].idx != atom1_prmtop:
                #             raise ValueError(f"Macrocycle dihedral bond atom index mismatch: expected {atom1_prmtop}, got {self.parm.atoms[atom1_prmtop].idx}")
                #         atom2_prmtop = p2+self.num_atom_offset
                #         if self.parm.atoms[atom2_prmtop].idx != atom2_prmtop:
                #             raise ValueError(f"Macrocycle dihedral bond atom index mismatch: expected {atom2_prmtop}, got {self.parm.atoms[atom2_prmtop].idx}")
                        
                #         # Mandatory
                #         if 'ring' in dihedral_type:
                #             continue

                #         # Optional
                #         if dihedral_type == 'protein-phi' and not rigid_protein_phi:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-psi' and not rigid_protein_psi:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-omega' and not rigid_protein_omega:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-chi1' and not rigid_protein_chi1:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-chi2' and not rigid_protein_chi2:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-chi3' and not rigid_protein_chi3:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-chi4' and not rigid_protein_chi4:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))
                #         if dihedral_type == 'protein-chi5' and not rigid_protein_chi5:
                #             self.standard_dihedral_bonds.append((atom1_prmtop, atom2_prmtop))

                #     self.macrocycle_dihedral_bonds.append(macrocycle_bonds)

                # Add angles
                for angle in molecule_prototypes[prototype_index].angle_params:
                    self.bond_bends.append(
                        rb.RoboAngle(
                            global_indices=tuple(self.prmtop_to_global_index[l+self.num_atom_offset] for l in angle.local_indices),
                            compound_atom_indices=angle.compound_atom_indices,
                            stiffness_in_kj_per_rad_sq=angle.stiffness_in_kj_per_rad_sq,
                            nominal_angle_in_deg=angle.nominal_angle_in_deg,
                            molecule_index=instance_index,
                        )
                    )

                # Add periodic torsions
                for periodic_torsion in molecule_prototypes[prototype_index].periodic_torsion_params:
                    self.proper_periodic_torsions.append(
                        rb.RoboPeriodicTorsion(
                            global_indices=tuple(self.prmtop_to_global_index[l+self.num_atom_offset] for l in periodic_torsion.local_indices),
                            compound_atom_indices=periodic_torsion.compound_atom_indices,
                            molecule_index=instance_index,
                            improper=periodic_torsion.is_improper,
                            terms=periodic_torsion.terms
                        )
                    )

                # Add harmonic improper torsions
                for harmonic_improper_torsion in molecule_prototypes[prototype_index].improper_harmonic_torsion_terms:
                    self.improper_harmonic_torsions.append(
                        rb.RoboHarmonicImproperTorsion(
                            global_indices=tuple(self.prmtop_to_global_index[l+self.num_atom_offset] for l in harmonic_improper_torsion.local_indices),
                            compound_atom_indices=harmonic_improper_torsion.compound_atom_indices,
                            stiffness_in_kj_per_rad_sq=harmonic_improper_torsion.stiffness_in_kj_per_rad_sq,
                            nominal_angle_in_rad=harmonic_improper_torsion.nominal_angle_in_rad,
                            molecule_index=instance_index,
                        )
                    )

                # Add Z matrix atom indices for this molecule
                for row in molecule_prototypes[prototype_index].z_matrix:
                    global_indices = (
                        self.prmtop_to_global_index[row.i_global_atom_index+self.num_atom_offset],
                        self.prmtop_to_global_index[row.j_global_atom_index+self.num_atom_offset] if row.j_global_atom_index is not None else -1,
                        self.prmtop_to_global_index[row.k_global_atom_index+self.num_atom_offset] if row.k_global_atom_index is not None else -1,
                        self.prmtop_to_global_index[row.l_global_atom_index+self.num_atom_offset] if row.l_global_atom_index is not None else -1,
                    )
                    compound_atom_indices = (
                        row.i_compound_atom_index,
                        row.j_compound_atom_index if row.j_compound_atom_index is not None else -1,
                        row.k_compound_atom_index if row.k_compound_atom_index is not None else -1,
                        row.l_compound_atom_index if row.l_compound_atom_index is not None else -1,
                    )
                    row = rb.ZMatrixRow(global_indices=global_indices, compound_atom_indices=compound_atom_indices, molecule_index=instance_index)
                    self.z_matrix.append(row)

        self.scaling14s = list[rb.Scaling14]()
        excluded = set()
        
        length_conv = pmd.unit.angstrom.conversion_factor_to(pmd.unit.nanometers)
        ene_conv = pmd.unit.kilocalories_per_mole.conversion_factor_to(pmd.unit.kilojoules_per_mole)

        dihedral_pointers = self.parm.parm_data["DIHEDRALS_INC_HYDROGEN"] + self.parm.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        for ii in range(0, len(dihedral_pointers), 5):
            i, j, k, l, dihedral_type_index = dihedral_pointers[ii:ii+5]

            # If negative, the 1-4 non-bonded interactions for this specific dihedral are not calculated
            # This prevents double-counting if the atoms are already part of a ring or another excluded group
            if k < 0:
                continue

            # If negative, this identifies the dihedral as an improper dihedral
            if l < 0:
                continue

            atom1_global_index = self.prmtop_to_global_index[i//3]
            atom1 = self.atoms[atom1_global_index]
            atom1_charge = parm_data['CHARGE'][i//3]

            atom4_global_index = self.prmtop_to_global_index[l//3]
            atom4 = self.atoms[atom4_global_index]
            atom4_charge = parm_data['CHARGE'][l//3]

            idx = parm_data['NONBONDED_PARM_INDEX'][atom1.identity.nonbonded_index * self.num_types + atom4.identity.nonbonded_index] - 1
            if idx < 0:
                continue

            if len(parm_data.get('LENNARD_JONES_14_ACOEF', [])) > 0:
                acoef = parm_data['LENNARD_JONES_14_ACOEF'][idx]
            else:
                acoef = parm_data['LENNARD_JONES_ACOEF'][idx]

            if len(parm_data.get('LENNARD_JONES_14_BCOEF', [])) > 0:
                bcoef = parm_data['LENNARD_JONES_14_BCOEF'][idx]
            else:
                bcoef = parm_data['LENNARD_JONES_BCOEF'][idx]

            if acoef != 0.0 and bcoef != 0.0:
                epsilon = (bcoef * bcoef) / (4 * acoef) * ene_conv
                r_min = (2*acoef/bcoef)**(1/6.0) * length_conv
            else:
                epsilon = 0.0
                r_min = 1.0

            charge_product = (atom1_charge * atom4_charge) / parm_data['SCEE_SCALE_FACTOR'][dihedral_type_index-1]
            epsilon /= parm_data['SCNB_SCALE_FACTOR'][dihedral_type_index-1]
            sigma = r_min * self.SIGMA_SCALE

            key = (min(atom1_global_index, atom4_global_index), max(atom1_global_index, atom4_global_index))
            if key in excluded:
                continue

            excluded.add(key)
            self.scaling14s.append(rb.Scaling14(a1=key[0], a4=key[1], charge_product=charge_product, epsilon=epsilon, sigma=sigma))

        numExcludedAtomsList=self.parm.parm_data["NUMBER_EXCLUDED_ATOMS"]
        excludedAtomsList=self.parm.parm_data["EXCLUDED_ATOMS_LIST"]
        self.exclusions = list[rb.Exclusion]()
        total=0
        for iAtom in range(self.parm.ptr('NATOM')):
            index0=total
            n=int(numExcludedAtomsList[iAtom])
            total+=n
            index1=total
            for jAtom in excludedAtomsList[index0:index1]:
                j=int(jAtom)
                if j>0:
                    global_i = self.prmtop_to_global_index[iAtom]
                    global_j = self.prmtop_to_global_index[j-1]
                    key = (min(global_i, global_j), max(global_i, global_j))
                    if key not in excluded:
                        excluded.add(key)
                        self.exclusions.append(rb.Exclusion(a1=key[0], a2=key[1]))

        # Now that we have mapped all atom indices, we need to update the neighbor indices to this mapping
        mapping = self.prmtop_to_global_index.get # Use .get for safety, or just the dict
        for atom in self.atoms:
            conn = atom.connectivity
            conn.neighbors_global_indices = list(map(mapping, conn.neighbors_global_indices))

        # Now parse CMAPs
        self.cmap_torsions: list[rb.CMAPTorsion] = []
        self.cmap_grids: list[rb.CMAPGrid] = []
        
        cmap_resolution = self.parm.parm_data.get('CMAP_RESOLUTION', [])
        num_cmap_grids = len(cmap_resolution)
        for i in range(num_cmap_grids):
            res = cmap_resolution[i]

            key = f'CMAP_PARAMETER_{i+1:02d}'
            cmap = self.parm.parm_data[key]

            # A CMAP defines a 2D correction surface E(phi, psi) for one specific pair of backbone dihedrals (typically phi and psi)
            # Typically, size = 24, so phi and psi are each discretized into 24 bins over [-pi, pi), [-pi, pi), giving 24*24=576 grid points
            # Those 576 values define one reusable energy surface, not per-atom data
            # For each residue assigned a CMAP, OpenMM computes the current phi and psi angles from four atoms per dihedral (8 atoms total, with overlap)
            # The angles are mapped onto the 2D grid.
            # Energy and forces are obtained by periodic bicubic interpolation of that same 576-point map.
            # The same map is reused for every residue that references it
            grid = rb.CMAPGrid()
            grid.size = res
            
            # Create a 1D list for OpenMM (phi varies fastest)
            # OpenMM_Index = phi_idx + res * psi_idx
            new_energy = [0.0] * (res * res)
            
            # We need to shift by 180 degrees (half the resolution)
            half = res // 2

            for phi_omm in range(res):
                for psi_omm in range(res):
                    # Map OpenMM index (starts at 0) back to Amber index (starts at -180)
                    # OpenMM 0 deg is Amber Index 12 (if res=24)
                    phi_amber = (phi_omm + half) % res
                    psi_amber = (psi_omm + half) % res
                    
                    # Amber/ParmEd index: psi changes fastest
                    old_index = phi_amber * res + psi_amber
                    
                    # OpenMM index: phi changes fastest
                    new_index = phi_omm + res * psi_omm
                    
                    # Convert from kcal/mol to kJ/mol
                    new_energy[new_index] = cmap[old_index] * 4.184

            grid.energy = new_energy
            self.cmap_grids.append(grid)

        # Add torsions that need correction from CMAPs
        cmap_index = self.parm.parm_data.get('CMAP_INDEX', [])
        for i in range(0, len(cmap_index), 6):
            map_index = cmap_index[i+5]
            if map_index < 1 or map_index > num_cmap_grids:
                raise ValueError(f"CMAP index out of range: {map_index}, number of CMAP grids: {num_cmap_grids}")
            
            torsion = rb.CMAPTorsion()
            torsion.mapIndex = map_index - 1

            torsion.a1 = self.prmtop_to_global_index[cmap_index[i+0]-1]
            torsion.a2 = self.prmtop_to_global_index[cmap_index[i+1]-1]
            torsion.a3 = self.prmtop_to_global_index[cmap_index[i+2]-1]
            torsion.a4 = self.prmtop_to_global_index[cmap_index[i+3]-1]

            torsion.b1 = self.prmtop_to_global_index[cmap_index[i+1]-1]
            torsion.b2 = self.prmtop_to_global_index[cmap_index[i+2]-1]
            torsion.b3 = self.prmtop_to_global_index[cmap_index[i+3]-1]
            torsion.b4 = self.prmtop_to_global_index[cmap_index[i+4]-1]

            self.cmap_torsions.append(torsion)

        # Add Urey-Bradley terms
        ub_terms = self.parm.parm_data.get('CHARMM_UREY_BRADLEY', [])
        ub_k = self.parm.parm_data.get('CHARMM_UREY_BRADLEY_FORCE_CONSTANT', [])
        ub_eq = self.parm.parm_data.get('CHARMM_UREY_BRADLEY_EQUIL_VALUE', [])

        self.urey_bradleys: list[rb.UreyBradley] = []
        for i in range(0, len(ub_terms), 3):
            ub = rb.UreyBradley()
            ub.a1 = self.prmtop_to_global_index[ub_terms[i]-1]
            ub.a3 = self.prmtop_to_global_index[ub_terms[i+1]-1]
            ub.stiffness_in_kj_per_nm_sq = ub_k[ub_terms[i+2]-1] * 4.184 * 100
            ub.nominal_length_in_nm = ub_eq[ub_terms[i+2]-1] / 10
            self.urey_bradleys.append(ub)

        print(f"Loaded system with {len(self.atoms)} atoms, {len(self.bond_stretches)} bonds, {len(self.bond_bends)} angles, {len(self.proper_periodic_torsions)} proper torsions, {len(self.improper_harmonic_torsions)} improper torsions, {len(self.cmap_torsions)} CMAP torsions, and {len(self.urey_bradleys)} Urey-Bradley terms.")

    @contextmanager
    def record_topology(self, molecule: MoleculePrototype) -> Iterable[rb.TopologyRange]:
        start_counts = self._get_current_counts()
        t_range = rb.TopologyRange(start_counts)

        yield t_range
    
        t_range.close(self._get_current_counts())
        self.topology_ranges.append(t_range)
        self.num_atom_offset += len(molecule.atom_params)
        self.num_residues_offset += molecule.num_residues
    
    def _get_current_counts(self) -> Tuple[int, int, int, int, int]:
        return (len(self.atoms),
                len(self.bond_stretches),
                len(self.bond_bends),
                len(self.proper_periodic_torsions),
                len(self.improper_harmonic_torsions))

    def getDefaultBonds(self, bonds_type) -> list[rb.BondFlexibility]:
        flexibilities = []

        if bonds_type == 'rama':
            raise ValueError("not supported lol")
            source = self.backbone_dihedral_bonds
        elif bonds_type == 'standard':
            source = [self.standard_dihedral_bonds]
        elif bonds_type == 'macrocycle':
            source = self.macrocycle_dihedral_bonds
        else:
            raise ValueError(f"Unsupported bonds_type: {bonds_type}")

        # original prmtop indices
        for sublist in source:
            flexibilities.append([])
            for bond in sublist:
                flex = rb.BondFlexibility()

                atom1 = self.parm.atoms[bond[0]]
                atom2 = self.parm.atoms[bond[1]]

                flex.globalIndex1 = self.prmtop_to_global_index[atom1.idx]
                flex.globalIndex2 = self.prmtop_to_global_index[atom2.idx]

                flex.uniqueAtomName1 = atom1.residue.name + str(atom1.residue.idx+1) + '_' + atom1.name + '_' + str(atom1.idx+1)
                flex.uniqueAtomName2 = atom2.residue.name + str(atom2.residue.idx+1) + '_' + atom2.name + '_' + str(atom2.idx+1)

                flex.mobility = rb.BondMobility.Torsion

                flexibilities[-1].append(flex)

        return flexibilities
    
    def selectBonds(self, query: str, excludeTerminal: bool = True) -> list[rb.BondFlexibility]:

        # Load files and run query
        traj = md.load(self.inpcrd, top=self.prmtop)
        sel = traj.topology.select(query)
        sub_traj = traj.atom_slice(sel)

        # Parse all selected bonds
        flexibilities = []
        for bond in sub_traj.topology.bonds:
            if excludeTerminal:
                if bond.atom1.n_bonds == 1 or bond.atom2.n_bonds == 1:
                    continue
            
            # Convert from prmtop indices to global indices
            global_index_1 = self.prmtop_to_global_index[bond.atom1.index]
            global_index_2 = self.prmtop_to_global_index[bond.atom2.index]

            # Create flexibility
            flex = rb.BondFlexibility()
            flex.globalIndex1 = global_index_1
            flex.globalIndex2 = global_index_2
            flex.uniqueAtomName1 = self.atoms[global_index_1].identity.unique_name
            flex.uniqueAtomName2 = self.atoms[global_index_2].identity.unique_name
            flex.mobility = rb.BondMobility.Torsion

            flexibilities.append(flex)

        return flexibilities
    
    def addCartesianWorld(self, samplesPerRound: int = 1) -> World:
        flexibilities = []
        for bond in self.bond_stretches:
            if bond.ring_closing:
                continue
            flex = rb.BondFlexibility()
            flex.globalIndex1 = bond.global_indices[0]
            flex.globalIndex2 = bond.global_indices[1]
            flex.uniqueAtomName1 = self.atoms[bond.global_indices[0]].identity.unique_name
            flex.uniqueAtomName2 = self.atoms[bond.global_indices[1]].identity.unique_name
            flex.mobility = rb.BondMobility.Translation
            flexibilities.append(flex)

        w = World(fixmanTorque=False,
                  samplesPerRound=samplesPerRound,
                  rootMobility=rb.RootMobility.WELD,
                  flexibilities=[flexibilities],
                  useOpenMM=True,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=True,
                  samplers=list[Sampler]())
        self.worlds.append(w)
        return self.worlds[-1]
    
    def addTorsionalWorld(self, torsional_bonds: list[rb.BondFlexibility], samplesperRound: int = 1):
        """ # !!!!!!!!!!!!!!!!!!!!! torsional bonds are in prmtop order, we reorder them in bat coordinates inside this function !!!!!!!!!!!!!!!!!!!!! """

        # # Reorder torsional bonds to match the global indices used in Robosample
        # torsional_bonds_reordered = []

        # for b in torsional_bonds:
        #     assert b.mobility == rb.BondMobility.Torsion, "Only torsional flexibilities are supported in torsional world. Rigid bonds are implicitly defined by the absence of a flexibility."

        #     # This map holds (a,b) and (b,a)
        #     id = bond_params.get((b.i, b.j))
        #     assert id is not None, f"Bond between atoms {b.i} and {b.j} not found in bond indices"
        #     torsional_bonds_reordered.append(rb.BondFlexibility(self.prmtop_to_global_index[b.i], self.prmtop_to_global_index[b.j], b.mobility))

        w = World(fixmanTorque=True,
                  samplesPerRound=samplesperRound,
                  rootMobility=rb.RootMobility.WELD,
                  flexibilities=torsional_bonds,
                  useOpenMM=True,
                  visual=False,
                  visualizerFrequency=0,
                  isCartesian=False,
                  samplers=list[Sampler]())
        self.worlds.append(w)

        return self.worlds[-1]
    
    def initialize(self, replicaTemperatures: list[float]):

        # Load the system into Robosample
        super().loadAmberSystem(
            self.root_indices,
            self.atoms,
            self.bond_stretches,
            self.bond_bends,
            self.proper_periodic_torsions,
            self.improper_harmonic_torsions,
            self.topology_ranges,
            self.z_matrix
        )
        
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

        # OpenMM must be initialized before adding samplers since they want to calculate energies when initializing
        super().initialize_openmm(
            self.atoms,
            self.bond_stretches,
            self.bond_bends,
            self.proper_periodic_torsions,
            self.improper_harmonic_torsions,
            self.cmap_grids,
            self.cmap_torsions,
            self.urey_bradleys,
            self.has_nbfix,
            self.num_types,
            self.acoef,
            self.bcoef,
            self.exclusions,
            self.scaling14s
        )

        # Add samplers
        for i, w in enumerate(self.worlds):
            s = w.samplers[0]
            super().getWorld(i).addSampler(s.samplerName, s.integratorType, s.thermostatName, s.useFixmanPotential)

        # Add replicas and thermodynamic states
        for temp in replicaTemperatures:
            super().addReplica()
            super().addThermodynamicState(temp, accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

        if not super().validate_context():
            raise ValueError("Invalid context.")

    def generate_synthetic_atom_classes(self):
        # Signatures store the "parameter environment" of each atom
        signatures = defaultdict(lambda: {"type": None, "mass": None, "bonds": [], "angles": [], "dihedrals": []})

        # Capture base types
        for atom in self.parm.atoms:
            signatures[atom]["type"] = atom.type
            signatures[atom]["mass"] = atom.mass

        # Capture Bond Environments
        for bond in self.parm.bonds:
            # Use id(bond.type) to ensure we distinguish between two 
            # BondType objects that might have different parameters
            bt_id = id(bond.type)
            signatures[bond.atom1]["bonds"].append(bt_id)
            signatures[bond.atom2]["bonds"].append(bt_id)

        # Capture Angle Environments
        for angle in self.parm.angles:
            at_id = id(angle.type)
            signatures[angle.atom1]["angles"].append(at_id)
            signatures[angle.atom2]["angles"].append(at_id)
            signatures[angle.atom3]["angles"].append(at_id)

        # Capture Dihedral Environments
        for dihed in self.parm.dihedrals:
            dt_id = id(dihed.type)
            signatures[dihed.atom1]["dihedrals"].append(dt_id)
            signatures[dihed.atom2]["dihedrals"].append(dt_id)
            signatures[dihed.atom3]["dihedrals"].append(dt_id)
            signatures[dihed.atom4]["dihedrals"].append(dt_id)

        # Collapse signatures into unique class names
        unique_sig_to_name = {}
        atom_to_class_name = {}
        type_counters = defaultdict(int)

        for atom in self.parm.atoms:
            # Create a hashable tuple representing the unique environment
            # Sorting is critical so that atom order in a bond doesn't change the signature
            sig = (
                signatures[atom]["type"],
                signatures[atom]["mass"],
                tuple(sorted(signatures[atom]["bonds"])),
                tuple(sorted(signatures[atom]["angles"])),
                tuple(sorted(signatures[atom]["dihedrals"]))
            )

            if sig not in unique_sig_to_name:
                base_type = signatures[atom]["type"]
                type_counters[base_type] += 1
                unique_sig_to_name[sig] = f"{base_type}_{type_counters[base_type]}"
            
            atom_to_class_name[atom] = unique_sig_to_name[sig]

        return atom_to_class_name
    