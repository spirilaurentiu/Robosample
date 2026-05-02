import warnings
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Iterable, Self, Tuple

import community as community_louvain
import MDAnalysis as mda
import mdtraj as md
import networkx as nx
import numpy as np
import pandas as pd
import parmed as pmd
import scipy.linalg as linalg
from MDAnalysis.analysis import dihedrals

from . import robo_bindings as rb
from .molecule_prototype import MoleculePrototype
from .prmtop_reader import has_nbfix_fast, parse_prmtop_numpy


@dataclass
class Sampler:
    samplerName: rb.SamplerName
    integratorType: rb.IntegratorType
    thermostatName: rb.ThermostatName
    useFixmanPotential: bool
    use_nuts: bool
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
    isCartesian: bool
    samplers: list[Sampler]

    def add_sampler(
        self,
        timeStep: float,
        mdSteps: int,
        boostMDSteps: int,
        samplerName: rb.SamplerName = rb.SamplerName.HMC,
        integratorType: rb.IntegratorType = rb.IntegratorType.VERLET,
        thermostatName: rb.ThermostatName = rb.ThermostatName.ANDERSEN,
        acceptRejectMode: rb.AcceptRejectMode = rb.AcceptRejectMode.MetropolisHastings,
        useFixmanPotential: bool = True,
        use_nuts: bool = False,
        distortOption: int = 0,
        distortArgs: str = "0",
        flow: int = 0,
    ) -> Self:
        assert len(self.samplers) == 0, (
            "Can only add samplers before initializing the context"
        )

        if self.isCartesian:
            integratorType = rb.IntegratorType.OMMVV
            useFixmanPotential = False

        s = Sampler(
            samplerName=samplerName,
            integratorType=integratorType,
            thermostatName=thermostatName,
            useFixmanPotential=useFixmanPotential,
            use_nuts=use_nuts,
            timeStep=timeStep,
            mdSteps=mdSteps,
            boostMDSteps=boostMDSteps,
            acceptRejectMode=acceptRejectMode,
            distortOption=distortOption,
            distortArgs=distortArgs,
            flow=flow,
        )
        self.samplers.append(s)
        return self


class Context(rb.Context):
    SIGMA_SCALE = 2 ** (-1.0 / 6.0)

    def __init__(
        self,
        name: str,
        seed: int,
        prmtop: str,
        inpcrd: str,
        write_freq: int,
        runType: rb.RunType = rb.RunType.REMC,
        replicaSwapFreq: int = 1,
        fixmanSwapFreq: int = 0,
        pdb_restart_freq: int = 0,
        nofRoundsTillReblock: int = 1,
        use_gbsa_obc2: bool = True,
        gbsa_solvent_dielectric: float = 78.5,
        gbsa_solute_dielectric: float = 1.0,
        nonbonded_method: rb.NonbondedMethod = rb.NonbondedMethod.CutoffNonPeriodic,
        nonbonded_cutoff_in_nm: float = 1.2,
        verbose: bool = False,
        testing: bool = False,
    ):

        super().__init__(
            name,
            seed,
            nofRoundsTillReblock,
            runType,
            replicaSwapFreq,
            fixmanSwapFreq,
            testing,
        )
        self.setPdbRestartFreq(pdb_restart_freq)
        self.setNonbonded(rb.NonbondedMethod.CutoffNonPeriodic, nonbonded_cutoff_in_nm)
        self.setVerbose(verbose)
        self.setGBSAOptions(
            use_gbsa_obc2, gbsa_solvent_dielectric, gbsa_solute_dielectric
        )

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

        self.system_topology = rb.SystemTopology()

        self.ff_params = rb.ForceFieldParams()
        self.ff_params.use_gbsaobc2 = use_gbsa_obc2
        self.ff_params.gbsa_solvent_dielectric = gbsa_solvent_dielectric
        self.ff_params.gbsa_solute_dielectric = gbsa_solute_dielectric
        self.ff_params.nonbonded_method = nonbonded_method
        self.ff_params.nonbonded_cutoff_in_nm = nonbonded_cutoff_in_nm

        self.sim_settings = rb.SimulationSettings(
            seed=seed, thermostat_temperature_in_k=300, collision_frequency=1.0
        )

        # parmed does nasty rounding when loading and loses some precision that adds up to a few kj
        # prmtop files hold more decimal places than can be stored via Python float64 (IEEE 754 double) has ~16 decimal digits of precision
        # prmtop holds more that 16, so this function will lose a few digits (fewer than parmed)
        parm_file = parse_prmtop_numpy(prmtop)
        parm_data = parm_file["raw_data"]

        # Nonbonded fix (NBFIX) is a technique that replaces standard Lennard-Jones (LJ) interaction parameters (epsilon and sigma) between specific atom pairs
        # This overrides default combination rules to fix overbinding artifacts, particularly between cations/anions and protein/lipid functional groups
        # It is commonly used in CHARMM force fields to improve hydration and binding accuracy
        self.ff_params.num_types = self.parm.pointers["NTYPES"]
        self.ff_params.has_nbfix = has_nbfix_fast(
            parm_data["NONBONDED_PARM_INDEX"],
            self.ff_params.num_types,
            parm_data["LENNARD_JONES_ACOEF"],
            parm_data["LENNARD_JONES_BCOEF"],
        )

        ene_conv = pmd.unit.kilocalories_per_mole.conversion_factor_to(
            pmd.unit.kilojoules_per_mole
        )
        length_conv = pmd.unit.angstroms.conversion_factor_to(pmd.unit.nanometers)
        afac = np.sqrt(ene_conv) * length_conv**6
        bfac = ene_conv * length_conv**6

        self.ff_params.a_coef = [
            0 for _ in range(self.ff_params.num_types * self.ff_params.num_types)
        ]
        self.ff_params.b_coef = [
            0 for _ in range(self.ff_params.num_types * self.ff_params.num_types)
        ]

        for i in range(self.ff_params.num_types):
            for j in range(self.ff_params.num_types):
                idx = (
                    parm_data["NONBONDED_PARM_INDEX"][i * self.ff_params.num_types + j]
                    - 1
                )
                if idx < 0:
                    raise ValueError(
                        f"Invalid nonbonded index for atom types {i} and {j}"
                    )
                self.ff_params.a_coef[i * self.ff_params.num_types + j] = (
                    np.sqrt(parm_data["LENNARD_JONES_ACOEF"][idx]) * afac
                )
                self.ff_params.b_coef[i * self.ff_params.num_types + j] = (
                    parm_data["LENNARD_JONES_BCOEF"][idx] * bfac
                )

        # DuMM atom classes are defined by their atom type (XC, C8, N3 etc), not atom name (N, CA, C, O etc)
        atom_classes = self.generate_synthetic_atom_classes()
        unique_atom_classes = set(atom_classes.values())
        unique_atom_classes = sorted(unique_atom_classes)
        self.atom_class_indices = {
            key: i + 1 for i, key in enumerate(unique_atom_classes)
        }

        # DuMM charged atom types are AMBER atom types plus their partial charge
        # DuMMForceFieldSubsystemRep::setBiotypeChargedAtomType - there is 1:1 correspondence between biotype and charged atom type
        charged_atom_types = set(
            [atom_classes[a] + ":" + str(a.charge) for a in self.parm.atoms]
        )

        # Sort to ensure consistent ordering, set() does not guarantee order
        charged_atom_types = sorted(charged_atom_types)

        self.charged_atom_type_indices = {
            atom_type: i for i, atom_type in enumerate(charged_atom_types)
        }

        self.z_matrix = list[tuple[int, int, int, int]]()
        self.num_atom_offset = 0
        self.num_residues_offset = 0

        # Map original prmtop atom indices to the order in which atoms are added to the Robosample context
        # This is the order in which atoms are explored via BFS starting from the root atom of each molecule offset by the number of atoms in previous molecules
        self.prmtop_to_global_index = {}

        # Compute secondary structure for each residue
        traj = md.load(self.inpcrd, top=self.prmtop)
        secondary_structure = md.compute_dssp(traj)

        # Store bonds that correspond to standardized dihedrals (protein backbone phi, psi, sidechain chi, etc)
        standard_dihedral_bond_columns = [
            "atom1_prmtop_index",
            "atom2_prmtop_index",
            "dihedral_type",
            "resname",
            "resid",
            "dss",
            "molecule_index",
        ]
        self.standard_dihedral_bonds = pd.DataFrame(
            columns=standard_dihedral_bond_columns
        )
        self.standard_dihedral_atom_groups = []

        # Split system into unique molecule prototypes and their occurrences
        # Each entry: (prototype_structure, instance_indices)
        # Example: parm_prototypes[0] holds the first molecule type found in the system
        # parm_prototypes[0][0] is the parmed Structure of that molecule type
        # parm_prototypes[0][1] is a list of all instance indices of that molecule type in the full system eg [0, 3, 5] if that molecule occurs at those indices
        parm_prototypes = self.parm.split()

        # Parsing atom, bond, angle and torsion definitions is expensive, so we cache for each prototype type
        molecule_prototypes = [
            MoleculePrototype(mol_struct) for mol_struct, _ in parm_prototypes
        ]

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
            with self.record_topology(molecule_prototypes[prototype_index]) as _:
                # Add atoms
                for atom in molecule_prototypes[prototype_index].atom_params:
                    prmtop_index = atom.local_index + self.num_atom_offset
                    a = self.parm.atoms[prmtop_index]
                    if a.idx != prmtop_index:
                        raise ValueError(
                            f"Atom index mismatch: expected {prmtop_index}, got {a.idx}"
                        )
                    if prmtop_index in self.prmtop_to_global_index:
                        raise ValueError(
                            f"Duplicate prmtop index found in mapping: {prmtop_index}"
                        )

                    # Map prmtop index to global index (BFS order)
                    global_index = len(self.prmtop_to_global_index)
                    self.prmtop_to_global_index[prmtop_index] = global_index

                    # Get atom class and charged atom type
                    atom_class_name = atom_classes[a]
                    charged_atom_type_name = atom_classes[a] + ":" + str(a.charge)

                    # Residue and atom indices in unique atom name are 1-based, not 0-based
                    unique_atom_name = (
                        a.residue.name
                        + str(a.residue.idx + 1)
                        + "_"
                        + a.name
                        + "_"
                        + str(prmtop_index + 1)
                    )  # e.g. ALA1_N_4

                    # First atom is always the root and holds magic properties
                    if atom.root:
                        unique_atom_name += "_ROOT"  # e.g. ALA1_N_4_ROOT
                        self.system_topology.root_atom_global_indices.append(
                            global_index
                        )

                    # Get Lennard-Jones parameters for this atom
                    # Nonbonded indices are 1-based
                    idx = (a.nb_idx - 1) * self.ff_params.num_types + (a.nb_idx - 1)
                    if idx < 0:
                        raise ValueError(f"Invalid nonbonded index for atom {a.idx}")
                    nb_parm_idx = parm_data["NONBONDED_PARM_INDEX"][idx] - 1
                    acoef = parm_data["LENNARD_JONES_ACOEF"][nb_parm_idx]
                    bcoef = parm_data["LENNARD_JONES_BCOEF"][nb_parm_idx]

                    # Parameters may be undefined for some atoms, typical hydrogen atoms
                    if acoef != 0.0 and bcoef != 0.0:
                        r_min = (2 * acoef / bcoef) ** (1 / 6.0)
                        epsilon = 0.25 * bcoef * bcoef / acoef
                    else:
                        r_min = 1.0
                        epsilon = 0.0

                    # Convert to sigma and epsilon
                    lengthConversionFactor = pmd.unit.angstrom.conversion_factor_to(
                        pmd.unit.nanometer
                    )
                    energyConversionFactor = (
                        pmd.unit.kilocalorie_per_mole.conversion_factor_to(
                            pmd.unit.kilojoule_per_mole
                        )
                    )

                    rVdw = r_min / 2.0 * lengthConversionFactor
                    sigma = rVdw * 2.0 * self.SIGMA_SCALE
                    epsilon = epsilon * energyConversionFactor

                    # Define the atom
                    self.system_topology.atoms.append(
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
                                nonbonded_index=a.nb_idx - 1,
                                compound_atom_index=atom.compound_atom_index,
                                atom_class_index=self.atom_class_indices[
                                    atom_class_name
                                ],
                                charged_atom_type_index=self.charged_atom_type_indices[
                                    charged_atom_type_name
                                ],
                            ),
                            element_info=rb.RoboAtomElement(
                                element_name=atom.element_name,
                                element_symbol=atom.element_symbol,
                                atomic_number=atom.atomic_number,
                            ),
                            physics=rb.RoboAtomPhysics(
                                charge_e=parm_data["CHARGE"][prmtop_index],
                                mass_daltons=parm_data["MASS"][prmtop_index],
                                vdw_radius_nm=atom.vdw_radius_nm,
                                vdw_well_depth_kj=epsilon,
                                sigma_nm=sigma,
                                solvent_radius_nm=parm_data["RADII"][prmtop_index] / 10,
                                screen=parm_data["SCREEN"][prmtop_index],
                            ),
                            connectivity=rb.RoboAtomConnectivity(
                                neighbors_global_indices=[
                                    n.idx for n in a.bond_partners
                                ],
                                root=atom.root,
                            ),
                            position=[a.xx / 10, a.xy / 10, a.xz / 10],
                        )
                    )

                # Add bonds
                for bond in molecule_prototypes[prototype_index].bond_params:
                    prmtop_indices = [
                        bond.parent_local_index + self.num_atom_offset,
                        bond.child_local_index + self.num_atom_offset,
                    ]
                    for prmtop_index in prmtop_indices:
                        atom = self.parm.atoms[prmtop_index]
                        if atom.idx != prmtop_index:
                            raise ValueError(
                                f"Bond atom index mismatch: expected {prmtop_index}, got {atom.idx}"
                            )

                    self.system_topology.bonds.append(
                        rb.RoboBond(
                            global_indices=tuple(
                                self.prmtop_to_global_index[p] for p in prmtop_indices
                            ),
                            prmtop_indices=tuple(prmtop_indices),
                            compound_atom_indices=(
                                bond.parent_compound_atom_index,
                                bond.child_compound_atom_index,
                            ),
                            stiffness_in_kj_per_nm_sq=bond.stiffness_in_kj_per_nm_sq,
                            nominal_length_in_nm=bond.nominal_length_in_nm,
                            molecule_index=instance_index,
                            ring_closing=bond.is_ring_closing,
                            dihedral_type=bond.dihedral_type,
                        )
                    )

                # Add bonds which are the middle bond of a standard dihedral
                # We don't convert to global indices and must keep original prmtop ones
                for p1, p2, dihedral_type, resid in molecule_prototypes[
                    prototype_index
                ].bonds:
                    atom1_prmtop = p1 + self.num_atom_offset
                    if self.parm.atoms[atom1_prmtop].idx != atom1_prmtop:
                        raise ValueError(
                            f"Standard dihedral bond atom index mismatch: expected {atom1_prmtop}, got {self.parm.atoms[atom1_prmtop].idx}"
                        )
                    atom2_prmtop = p2 + self.num_atom_offset
                    if self.parm.atoms[atom2_prmtop].idx != atom2_prmtop:
                        raise ValueError(
                            f"Standard dihedral bond atom index mismatch: expected {atom2_prmtop}, got {self.parm.atoms[atom2_prmtop].idx}"
                        )

                    # Always skip ring closing bonds
                    if "ring" in dihedral_type:
                        continue

                    # Define standardized dihedral names
                    if dihedral_type != "non-standard":
                        new_dihedral_bond = {
                            "atom1_prmtop_index": atom1_prmtop,
                            "atom2_prmtop_index": atom2_prmtop,
                            "dihedral_type": dihedral_type,
                            "resname": self.parm.residues[resid].name,
                            "resid": resid,
                            "dss": secondary_structure[0][resid],
                            "molecule_index": instance_index,
                        }
                        self.standard_dihedral_bonds = pd.concat(
                            [
                                self.standard_dihedral_bonds,
                                pd.DataFrame([new_dihedral_bond]),
                            ],
                            ignore_index=True,
                        )

                        gparent, parent, child, gchild, _ = molecule_prototypes[
                            prototype_index
                        ].get_standardized_dihedral_type_2(
                            self.parm.atoms[atom1_prmtop], self.parm.atoms[atom2_prmtop]
                        )
                        self.standard_dihedral_atom_groups.append(
                            (gparent.idx, parent.idx, child.idx, gchild.idx),
                        )

                # Add angles
                for angle in molecule_prototypes[prototype_index].angle_params:
                    prmtop_indices = [
                        local_index + self.num_atom_offset
                        for local_index in angle.local_indices
                    ]
                    for prmtop_index in prmtop_indices:
                        atom = self.parm.atoms[prmtop_index]
                        if atom.idx != prmtop_index:
                            raise ValueError(
                                f"Bond atom index mismatch: expected {prmtop_index}, got {atom.idx}"
                            )

                    self.system_topology.angles.append(
                        rb.RoboAngle(
                            global_indices=tuple(
                                self.prmtop_to_global_index[
                                    local_index + self.num_atom_offset
                                ]
                                for local_index in angle.local_indices
                            ),
                            prmtop_indices=tuple(prmtop_indices),
                            compound_atom_indices=angle.compound_atom_indices,
                            stiffness_in_kj_per_rad_sq=angle.stiffness_in_kj_per_rad_sq,
                            nominal_angle_in_deg=angle.nominal_angle_in_deg,
                            molecule_index=instance_index,
                        )
                    )

                # Add periodic torsions
                for periodic_torsion in molecule_prototypes[
                    prototype_index
                ].periodic_torsion_params:
                    prmtop_indices = [
                        local_index + self.num_atom_offset
                        for local_index in periodic_torsion.local_indices
                    ]
                    for prmtop_index in prmtop_indices:
                        atom = self.parm.atoms[prmtop_index]
                        if atom.idx != prmtop_index:
                            raise ValueError(
                                f"Bond atom index mismatch: expected {prmtop_index}, got {atom.idx}"
                            )

                    self.system_topology.periodic_torsions.append(
                        rb.RoboPeriodicTorsion(
                            global_indices=tuple(
                                self.prmtop_to_global_index[
                                    local_index + self.num_atom_offset
                                ]
                                for local_index in periodic_torsion.local_indices
                            ),
                            prmtop_indices=tuple(prmtop_indices),
                            compound_atom_indices=periodic_torsion.compound_atom_indices,
                            molecule_index=instance_index,
                            improper=periodic_torsion.is_improper,
                            terms=periodic_torsion.terms,
                        )
                    )

                # Add harmonic improper torsions
                for harmonic_improper_torsion in molecule_prototypes[
                    prototype_index
                ].improper_harmonic_torsion_terms:
                    prmtop_indices = [
                        local_index + self.num_atom_offset
                        for local_index in harmonic_improper_torsion.local_indices
                    ]
                    for prmtop_index in prmtop_indices:
                        atom = self.parm.atoms[prmtop_index]
                        if atom.idx != prmtop_index:
                            raise ValueError(
                                f"Bond atom index mismatch: expected {prmtop_index}, got {atom.idx}"
                            )

                    self.system_topology.harmonic_improper_torsions.append(
                        rb.RoboHarmonicImproperTorsion(
                            global_indices=tuple(
                                self.prmtop_to_global_index[
                                    local_index + self.num_atom_offset
                                ]
                                for local_index in harmonic_improper_torsion.local_indices
                            ),
                            prmtop_indices=tuple(prmtop_indices),
                            compound_atom_indices=harmonic_improper_torsion.compound_atom_indices,
                            stiffness_in_kj_per_rad_sq=harmonic_improper_torsion.stiffness_in_kj_per_rad_sq,
                            nominal_angle_in_rad=harmonic_improper_torsion.nominal_angle_in_rad,
                            molecule_index=instance_index,
                        )
                    )

                # Add Z matrix atom indices for this molecule
                for row in molecule_prototypes[prototype_index].z_matrix:
                    global_indices = (
                        self.prmtop_to_global_index[
                            row.i_global_atom_index + self.num_atom_offset
                        ],
                        self.prmtop_to_global_index[
                            row.j_global_atom_index + self.num_atom_offset
                        ]
                        if row.j_global_atom_index is not None
                        else -1,
                        self.prmtop_to_global_index[
                            row.k_global_atom_index + self.num_atom_offset
                        ]
                        if row.k_global_atom_index is not None
                        else -1,
                        self.prmtop_to_global_index[
                            row.l_global_atom_index + self.num_atom_offset
                        ]
                        if row.l_global_atom_index is not None
                        else -1,
                    )
                    compound_atom_indices = (
                        row.i_compound_atom_index,
                        row.j_compound_atom_index
                        if row.j_compound_atom_index is not None
                        else -1,
                        row.k_compound_atom_index
                        if row.k_compound_atom_index is not None
                        else -1,
                        row.l_compound_atom_index
                        if row.l_compound_atom_index is not None
                        else -1,
                    )
                    row = rb.ZMatrixRow(
                        global_indices=global_indices,
                        compound_atom_indices=compound_atom_indices,
                        molecule_index=instance_index,
                    )
                    self.z_matrix.append(row)

        excluded = set()
        length_conv = pmd.unit.angstrom.conversion_factor_to(pmd.unit.nanometers)
        ene_conv = pmd.unit.kilocalories_per_mole.conversion_factor_to(
            pmd.unit.kilojoules_per_mole
        )

        dihedral_pointers = (
            self.parm.parm_data["DIHEDRALS_INC_HYDROGEN"]
            + self.parm.parm_data["DIHEDRALS_WITHOUT_HYDROGEN"]
        )
        for ii in range(0, len(dihedral_pointers), 5):
            index_i, index_j, index_k, index_l, dihedral_type_index = dihedral_pointers[
                ii : ii + 5
            ]

            # If negative, the 1-4 non-bonded interactions for this specific dihedral are not calculated
            # This prevents double-counting if the atoms are already part of a ring or another excluded group
            if index_k < 0:
                continue

            # If negative, this identifies the dihedral as an improper dihedral
            if index_l < 0:
                continue

            atom1_global_index = self.prmtop_to_global_index[index_i // 3]
            atom1 = self.system_topology.atoms[atom1_global_index]
            atom1_charge = parm_data["CHARGE"][index_i // 3]

            atom4_global_index = self.prmtop_to_global_index[index_l // 3]
            atom4 = self.system_topology.atoms[atom4_global_index]
            atom4_charge = parm_data["CHARGE"][index_l // 3]

            idx = (
                parm_data["NONBONDED_PARM_INDEX"][
                    atom1.identity.nonbonded_index * self.ff_params.num_types
                    + atom4.identity.nonbonded_index
                ]
                - 1
            )
            if idx < 0:
                continue

            if len(parm_data.get("LENNARD_JONES_14_ACOEF", [])) > 0:
                acoef = parm_data["LENNARD_JONES_14_ACOEF"][idx]
            else:
                acoef = parm_data["LENNARD_JONES_ACOEF"][idx]

            if len(parm_data.get("LENNARD_JONES_14_BCOEF", [])) > 0:
                bcoef = parm_data["LENNARD_JONES_14_BCOEF"][idx]
            else:
                bcoef = parm_data["LENNARD_JONES_BCOEF"][idx]

            if acoef != 0.0 and bcoef != 0.0:
                epsilon = (bcoef * bcoef) / (4 * acoef) * ene_conv
                r_min = (2 * acoef / bcoef) ** (1 / 6.0) * length_conv
            else:
                epsilon = 0.0
                r_min = 1.0

            charge_product = (atom1_charge * atom4_charge) / parm_data[
                "SCEE_SCALE_FACTOR"
            ][dihedral_type_index - 1]
            epsilon /= parm_data["SCNB_SCALE_FACTOR"][dihedral_type_index - 1]
            sigma = r_min * self.SIGMA_SCALE

            key = (
                min(atom1_global_index, atom4_global_index),
                max(atom1_global_index, atom4_global_index),
            )
            if key in excluded:
                continue

            excluded.add(key)
            self.system_topology.scaling14s.append(
                rb.Scaling14(
                    atom_1_global_index=key[0],
                    atom_4_global_index=key[1],
                    charge_product=charge_product,
                    epsilon=epsilon,
                    sigma=sigma,
                )
            )

        numExcludedAtomsList = self.parm.parm_data["NUMBER_EXCLUDED_ATOMS"]
        excludedAtomsList = self.parm.parm_data["EXCLUDED_ATOMS_LIST"]
        total = 0
        for iAtom in range(self.parm.ptr("NATOM")):
            index0 = total
            n = int(numExcludedAtomsList[iAtom])
            total += n
            index1 = total
            for jAtom in excludedAtomsList[index0:index1]:
                j = int(jAtom)
                if j > 0:
                    global_i = self.prmtop_to_global_index[iAtom]
                    global_j = self.prmtop_to_global_index[j - 1]
                    key = (min(global_i, global_j), max(global_i, global_j))
                    if key not in excluded:
                        excluded.add(key)
                        self.system_topology.exclusions.append(
                            rb.Exclusion(
                                atom_1_global_index=key[0], atom_2_global_index=key[1]
                            )
                        )

        # Now that we have mapped all atom indices, we need to update the neighbor indices to this mapping
        mapping = (
            self.prmtop_to_global_index.get
        )  # Use .get for safety, or just the dict
        for atom in self.system_topology.atoms:
            conn = atom.connectivity
            conn.neighbors_global_indices = list(
                map(mapping, conn.neighbors_global_indices)
            )

        # Now parse CMAPs
        cmap_resolution = self.parm.parm_data.get("CMAP_RESOLUTION", [])
        num_cmap_grids = len(cmap_resolution)
        for i in range(num_cmap_grids):
            res = cmap_resolution[i]

            key = f"CMAP_PARAMETER_{i + 1:02d}"
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
            self.system_topology.cmap_grids.append(grid)

        # Add torsions that need correction from CMAPs
        cmap_index = self.parm.parm_data.get("CMAP_INDEX", [])
        for i in range(0, len(cmap_index), 6):
            map_index = cmap_index[i + 5]
            if map_index < 1 or map_index > num_cmap_grids:
                raise ValueError(
                    f"CMAP index out of range: {map_index}, number of CMAP grids: {num_cmap_grids}"
                )

            torsion = rb.CMAPTorsion()
            torsion.mapIndex = map_index - 1

            torsion.torsion_a_atom_1_global_index = self.prmtop_to_global_index[
                cmap_index[i + 0] - 1
            ]
            torsion.torsion_a_atom_2_global_index = self.prmtop_to_global_index[
                cmap_index[i + 1] - 1
            ]
            torsion.torsion_a_atom_3_global_index = self.prmtop_to_global_index[
                cmap_index[i + 2] - 1
            ]
            torsion.torsion_a_atom_4_global_index = self.prmtop_to_global_index[
                cmap_index[i + 3] - 1
            ]

            torsion.torsion_b_atom_1_global_index = self.prmtop_to_global_index[
                cmap_index[i + 1] - 1
            ]
            torsion.torsion_b_atom_2_global_index = self.prmtop_to_global_index[
                cmap_index[i + 2] - 1
            ]
            torsion.torsion_b_atom_3_global_index = self.prmtop_to_global_index[
                cmap_index[i + 3] - 1
            ]
            torsion.torsion_b_atom_4_global_index = self.prmtop_to_global_index[
                cmap_index[i + 4] - 1
            ]

            self.system_topology.cmap_torsions.append(torsion)

        # Add Urey-Bradley terms
        ub_terms = self.parm.parm_data.get("CHARMM_UREY_BRADLEY", [])
        ub_k = self.parm.parm_data.get("CHARMM_UREY_BRADLEY_FORCE_CONSTANT", [])
        ub_eq = self.parm.parm_data.get("CHARMM_UREY_BRADLEY_EQUIL_VALUE", [])

        for i in range(0, len(ub_terms), 3):
            ub = rb.UreyBradley()
            ub.atom_1_global_index = self.prmtop_to_global_index[ub_terms[i] - 1]
            ub.atom_3_global_index = self.prmtop_to_global_index[ub_terms[i + 1] - 1]
            ub.stiffness_in_kj_per_nm_sq = ub_k[ub_terms[i + 2] - 1] * 4.184 * 100
            ub.nominal_length_in_nm = ub_eq[ub_terms[i + 2] - 1] / 10
            self.system_topology.urey_bradleys.append(ub)

        print(
            f"Loaded system with {len(self.system_topology.atoms)} atoms, {len(self.system_topology.bonds)} bonds, {len(self.system_topology.angles)} angles, {len(self.system_topology.periodic_torsions)} proper torsions, {len(self.system_topology.harmonic_improper_torsions)} improper torsions, {len(self.system_topology.cmap_torsions)} CMAP torsions, and {len(self.system_topology.urey_bradleys)} Urey-Bradley terms."
        )

    @contextmanager
    def record_topology(
        self, molecule: MoleculePrototype
    ) -> Iterable[rb.TopologyRange]:
        start_counts = self._get_current_counts()
        t_range = rb.TopologyRange(start_counts)

        yield t_range

        t_range.close(self._get_current_counts())
        self.system_topology.topology_ranges.append(t_range)
        self.num_atom_offset += len(molecule.atom_params)
        self.num_residues_offset += molecule.num_residues

    def _get_current_counts(self) -> Tuple[int, int, int, int, int]:
        return (
            len(self.system_topology.atoms),
            len(self.system_topology.bonds),
            len(self.system_topology.angles),
            len(self.system_topology.periodic_torsions),
            len(self.system_topology.harmonic_improper_torsions),
        )

    def build_flexibilities(self, source: pd.DataFrame) -> list[rb.BondFlexibility]:
        """
        [[rb.BondFlexibility(), rb.BondFlexibility(), ...], [...], ...]
        """
        flexibilities = []
        for _, bond in source.iterrows():
            flex = rb.BondFlexibility()

            atom1 = self.parm.atoms[bond["atom1_prmtop_index"]]
            atom2 = self.parm.atoms[bond["atom2_prmtop_index"]]

            flex.globalIndex1 = self.prmtop_to_global_index[atom1.idx]
            flex.globalIndex2 = self.prmtop_to_global_index[atom2.idx]

            flex.uniqueAtomName1 = self.system_topology.atoms[
                flex.globalIndex1
            ].identity.unique_name
            flex.uniqueAtomName2 = self.system_topology.atoms[
                flex.globalIndex2
            ].identity.unique_name

            flex.mobility = rb.BondMobility.Torsion

            flexibilities.append(flex)

        return [flexibilities]

    def create_torsional_bonds(
        self, bond_indices: list[tuple[int, int]]
    ) -> list[rb.BondFlexibility]:
        flexibilities = []
        for bond in bond_indices:
            flex = rb.BondFlexibility()

            atom1 = self.parm.atoms[bond[0]]
            atom2 = self.parm.atoms[bond[1]]

            flex.globalIndex1 = self.prmtop_to_global_index[atom1.idx]
            flex.globalIndex2 = self.prmtop_to_global_index[atom2.idx]

            flex.uniqueAtomName1 = self.system_topology.atoms[
                flex.globalIndex1
            ].identity.unique_name
            flex.uniqueAtomName2 = self.system_topology.atoms[
                flex.globalIndex2
            ].identity.unique_name

            flex.mobility = rb.BondMobility.Torsion

            flexibilities.append(flex)

        return [flexibilities]

    def selectBonds(
        self, query: str, excludeTerminal: bool = True
    ) -> list[rb.BondFlexibility]:

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
            flex.uniqueAtomName1 = self.system_topology.atoms[
                global_index_1
            ].identity.unique_name
            flex.uniqueAtomName2 = self.system_topology.atoms[
                global_index_2
            ].identity.unique_name
            flex.mobility = rb.BondMobility.Torsion

            flexibilities.append(flex)

        return flexibilities

    def addCartesianWorld(self, samplesPerRound: int = 1) -> World:
        flexibilities = []
        for bond in self.system_topology.bonds:
            if bond.ring_closing:
                continue
            flex = rb.BondFlexibility()
            flex.globalIndex1 = bond.global_indices[0]
            flex.globalIndex2 = bond.global_indices[1]
            flex.uniqueAtomName1 = self.system_topology.atoms[
                bond.global_indices[0]
            ].identity.unique_name
            flex.uniqueAtomName2 = self.system_topology.atoms[
                bond.global_indices[1]
            ].identity.unique_name
            flex.mobility = rb.BondMobility.Translation
            flexibilities.append(flex)

        w = World(
            fixmanTorque=False,
            samplesPerRound=samplesPerRound,
            rootMobility=rb.RootMobility.WELD,
            flexibilities=[flexibilities],
            isCartesian=True,
            samplers=list[Sampler](),
        )
        self.worlds.append(w)
        return self.worlds[-1]

    def addTorsionalWorld(
        self, torsional_bonds: list[rb.BondFlexibility], samplesperRound: int = 1
    ):
        """# !!!!!!!!!!!!!!!!!!!!! torsional bonds are in prmtop order, we reorder them in bat coordinates inside this function !!!!!!!!!!!!!!!!!!!!!"""

        # # Reorder torsional bonds to match the global indices used in Robosample
        # torsional_bonds_reordered = []

        # for b in torsional_bonds:
        #     assert b.mobility == rb.BondMobility.Torsion, "Only torsional flexibilities are supported in torsional world. Rigid bonds are implicitly defined by the absence of a flexibility."

        #     # This map holds (a,b) and (b,a)
        #     id = bond_params.get((b.i, b.j))
        #     assert id is not None, f"Bond between atoms {b.i} and {b.j} not found in bond indices"
        #     torsional_bonds_reordered.append(rb.BondFlexibility(self.prmtop_to_global_index[b.i], self.prmtop_to_global_index[b.j], b.mobility))

        w = World(
            fixmanTorque=True,
            samplesPerRound=samplesperRound,
            rootMobility=rb.RootMobility.WELD,
            flexibilities=torsional_bonds,
            isCartesian=False,
            samplers=list[Sampler](),
        )
        self.worlds.append(w)

        return self.worlds[-1]

    def initialize(self, replicaTemperatures: list[float]):

        # Load the system into Robosample
        super().loadAmberSystem(
            self.system_topology,
            self.ff_params,
            self.sim_settings,
            self.z_matrix,
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
            assert len(w.samplers) == 1, (
                "Currently only one sampler per world is supported"
            )
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
            print(
                "Adding world with the following parameters: ",
                f"fixmanTorque={world.fixmanTorque}, ",
                f"samplesPerRound={world.samplesPerRound}, ",
                f"rootMobility={world.rootMobility}, ",
                f"flexibilities length={len(world.flexibilities)}, ",
            )
            super().add_world(
                world.fixmanTorque,
                world.samplesPerRound,
                world.rootMobility,
                world.flexibilities,
            )

        # OpenMM must be initialized before adding samplers since they want to calculate energies when initializing
        ok = super().initialize_openmm()
        if not ok:
            raise ValueError("Failed to initialize OpenMM system in Robosample.")

        # Add samplers
        for i, w in enumerate(self.worlds):
            s = w.samplers[0]
            super().getWorld(i).add_sampler(
                s.samplerName,
                s.integratorType,
                s.thermostatName,
                s.useFixmanPotential,
                s.use_nuts,
            )

        # Add replicas and thermodynamic states
        for temp in replicaTemperatures:
            super().addReplica()
            super().addThermodynamicState(
                temp,
                accept_reject_modes,
                distort_options,
                distort_args,
                flow,
                work,
                integrators,
                worldIndexes,
                timesteps,
                mdsteps,
            )

        # if not super().validate_context():
        #     raise ValueError("Invalid context.")

        print(
            "Context initialized successfully with the following parameters: ",
            f"\n\tNumber of atoms: {len(self.system_topology.atoms)}, ",
            f"\n\tNumber of bond stretches: {len(self.system_topology.bonds)}, ",
            f"\n\tNumber of bond bends: {len(self.system_topology.angles)}, ",
            f"\n\tNumber of periodic torsions: {len(self.system_topology.periodic_torsions)}, ",
            f"\n\tNumber of improper harmonic torsions: {len(self.system_topology.harmonic_improper_torsions)}, ",
            f"\n\tNumber of CMAP torsions: {len(self.system_topology.cmap_torsions)}, ",
            f"\n\tNumber of Urey-Bradley terms: {len(self.system_topology.urey_bradleys)}, ",
            f"\n\tNumber of exclusions: {len(self.system_topology.exclusions)}, ",
            f"\n\tNumber of 1-4 scalings: {len(self.system_topology.scaling14s)}, ",
            f"\n\tNumber of worlds: {len(self.worlds)}, ",
            f"\n\tNumber of replicas: {len(replicaTemperatures)}, ",
            f"\n\tNumber of thermodynamic states: {len(replicaTemperatures)}",
        )

    def generate_synthetic_atom_classes(self):
        # Signatures store the "parameter environment" of each atom
        signatures = defaultdict(
            lambda: {
                "type": None,
                "mass": None,
                "bonds": [],
                "angles": [],
                "dihedrals": [],
            }
        )

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
                tuple(sorted(signatures[atom]["dihedrals"])),
            )

            if sig not in unique_sig_to_name:
                base_type = signatures[atom]["type"]
                type_counters[base_type] += 1
                unique_sig_to_name[sig] = f"{base_type}_{type_counters[base_type]}"

            atom_to_class_name[atom] = unique_sig_to_name[sig]

        return atom_to_class_name

    def _circmean_vectorised(self, angles: np.ndarray, axis: int = 0) -> np.ndarray:
        return np.arctan2(
            np.mean(np.sin(angles), axis=axis),
            np.mean(np.cos(angles), axis=axis),
        )

    def _circvar_vectorised(self, angles: np.ndarray, axis: int = 0) -> np.ndarray:
        R_bar = np.sqrt(
            np.mean(np.cos(angles), axis=axis) ** 2
            + np.mean(np.sin(angles), axis=axis) ** 2
        )
        return 1.0 - R_bar

    def circcorr(
        self,
        dihs: np.ndarray,
        pair: tuple,
        shuffle: bool = False,
    ) -> tuple:
        ix, jx = pair

        # Extract the two time series.
        # .copy() is mandatory before in-place shuffling to avoid mutating the
        # shared `dihs` array, which would corrupt concurrent calculations.
        x = dihs[:, ix].copy() if shuffle else dihs[:, ix]
        y = dihs[:, jx].copy() if shuffle else dihs[:, jx]

        if shuffle:
            # Independent permutations preserve marginals but destroy all
            # temporal structure and cross-series coupling simultaneously.
            x = self.rng.permutation(x)
            y = self.rng.permutation(y)

        x_var = float(self._circvar_vectorised(x))
        y_var = float(self._circvar_vectorised(y))
        if x_var < self.tol or y_var < self.tol:
            return ix, jx, 0.0

        mu_x = float(self._circmean_vectorised(x))
        mu_y = float(self._circmean_vectorised(y))

        sin_x = np.sin(x - mu_x)
        sin_y = np.sin(y - mu_y)

        numerator = np.sum(sin_x * sin_y)
        denom_sq = np.sum(sin_x**2) * np.sum(sin_y**2)

        if denom_sq < np.finfo(np.float64).eps:
            return ix, jx, 0.0

        correlation = numerator / np.sqrt(denom_sq)
        correlation = float(np.clip(correlation, -1.0, 1.0))

        return ix, jx, correlation

    def compute_dihedral_correlation_matrix(self, dcd_file: str) -> np.ndarray:
        universe = mda.Universe(self.prmtop, dcd_file)
        num_dihedrals = len(self.standard_dihedral_atom_groups)

        atom_groups = []
        for gparent, parent, child, gchild in self.standard_dihedral_atom_groups:
            atom_groups.append(
                mda.AtomGroup(universe.atoms[[gparent, parent, child, gchild]])
            )

        values = dihedrals.Dihedral(atom_groups).run().angles
        values = np.deg2rad(values).astype(np.float64)  # Shape: (T, N)
        print(type(values), values.shape)

        circ_means = self._circmean_vectorised(values, axis=0)
        circ_vars = self._circvar_vectorised(values, axis=0)
        low_var_mask = circ_vars < self.tol  # (N,) boolean

        if low_var_mask.sum() > 0:
            warnings.warn(
                f"{low_var_mask.sum()} dihedral(s) have circular variance < {self.tol} "
                f"and will be treated as uncorrelated (indices: "
                f"{np.where(low_var_mask)[0].tolist()}). "
                "These correspond to near-rigid bonds where correlation is undefined.",
                UserWarning,
                stacklevel=2,
            )

        # Centred sine matrix: S[t, j] = sin(θ_j(t) − μ_j)
        # Shape: (T, N)
        S = np.sin(values - circ_means[np.newaxis, :])

        # Numerator matrix via single BLAS dgemm: C = Sᵀ S
        # C[i, j] = Σₜ sin(θᵢ(t)−μᵢ)·sin(θⱼ(t)−μⱼ) — the JS numerator for pair (i,j)
        # Shape: (N, N)
        C = S.T @ S

        # Column 2-norms: d_j = √(C[j,j]) = √(Σₜ S[t,j]²)
        # Shape: (N,)
        col_norms = np.sqrt(np.diag(C))

        # Denominator matrix: D[i,j] = d_i · d_j (outer product)
        # Shape: (N, N)
        D = np.outer(col_norms, col_norms)

        # Elementwise division. D[i,j] = 0 only if dihedral i or j has a
        # constant sin-centred projection, caught by the variance filter below.
        # We suppress the divide-by-zero warning here; those entries are zeroed out.
        with np.errstate(invalid="ignore", divide="ignore"):
            correlation_matrix = np.where(D > 0.0, C / D, 0.0)

        # --- Apply low-variance mask -------------------------------------------
        # Zero all correlations involving near-rigid dihedrals (both row and
        # column). These entries would otherwise be 0/0 = NaN or numerically
        # large from floating-point noise.
        correlation_matrix[low_var_mask, :] = 0.0
        correlation_matrix[:, low_var_mask] = 0.0

        # Force exact unit diagonal: ρ(X,X) = 1 analytically; float64 arithmetic
        # may yield 0.999...9 or 1.000...1 due to rounding in the C/D division.
        np.fill_diagonal(correlation_matrix, 1.0)

        # Clamp all entries to [−1, 1] to absorb floating-point noise.
        np.clip(correlation_matrix, -1.0, 1.0, out=correlation_matrix)

        # --- Validity checks ---------------------------------------------------

        # Check (1): Square
        if correlation_matrix.shape[0] != correlation_matrix.shape[1]:
            raise ValueError(
                f"Correlation matrix is not square: shape is {correlation_matrix.shape}. "
                "This should be impossible with the vectorised implementation and "
                "indicates a bug in trajectory loading or atom group construction."
            )

        # Check (2): Symmetry
        # The vectorised formula produces a symmetric matrix by construction
        # (Sᵀ S is symmetric). Asymmetry here implies NaN propagation upstream.
        max_asymmetry = np.max(np.abs(correlation_matrix - correlation_matrix.T))
        if not np.allclose(correlation_matrix, correlation_matrix.T, atol=self.tol):
            raise ValueError(
                f"Correlation matrix is not symmetric (max |M − Mᵀ| = {max_asymmetry:.3e}). "
                "Check for NaN or Inf values in the input trajectory angles, which can "
                "arise from clashing atoms or corrupted DCD frames."
            )

        # Check (3): Unit diagonal
        min_diag = np.min(np.diag(correlation_matrix))
        if min_diag < 1.0 - self.tol:
            bad_indices = np.where(np.diag(correlation_matrix) < 1.0 - self.tol)[0]
            raise ValueError(
                f"Diagonal entries below 1 at dihedral indices: {bad_indices.tolist()}. "
                f"Minimum diagonal value: {min_diag:.6f}. "
                "These dihedrals have near-zero sin-centred projections that bypassed "
                "the circular variance filter. Consider increasing self.tol."
            )

        # Check (4): Positive definiteness via Cholesky.
        # A valid circular correlation matrix is PSD. Strict PD is required for
        # the downstream Gibbs sampler to define a proper conditional distribution.
        # Near-PSD matrices (all eigenvalues ≥ 0 but some < machine epsilon) arise
        # from near-collinear dihedral blocks and are handled by Tikhonov
        # regularisation: M ← M + εI.
        # Tikhonov regularisation is equivalent to assuming each dihedral carries
        # an independent isotropic noise floor of amplitude √ε — statistically
        # defensible and standard practice in covariance estimation.
        try:
            linalg.cholesky(correlation_matrix, lower=True)
        except linalg.LinAlgError:
            eps = 1e-6
            regularised = correlation_matrix + eps * np.eye(num_dihedrals)
            try:
                linalg.cholesky(regularised, lower=True)
                warnings.warn(
                    f"Correlation matrix is not strictly positive definite. "
                    f"Applied Tikhonov regularisation ε = {eps} (M ← M + εI). "
                    "This is equivalent to assuming an isotropic noise floor of "
                    f"√{eps:.0e} ≈ {np.sqrt(eps):.4f} rad on each dihedral. "
                    "Likely cause: near-collinear dihedral blocks or too few "
                    "independent frames relative to the number of dihedrals.",
                    UserWarning,
                    stacklevel=2,
                )
                correlation_matrix = regularised
            except linalg.LinAlgError:
                eigvals = np.linalg.eigvalsh(correlation_matrix)
                raise ValueError(
                    "Cholesky decomposition failed even after Tikhonov regularisation. "
                    f"Minimum eigenvalue: {eigvals.min():.4e}. "
                    "Possible causes: "
                    "(1) Duplicate dihedral atom groups producing rank-deficient rows/columns. "
                    "(2) Excessive low-variance zeroing creating a rank-deficient submatrix. "
                    "(3) Severely corrupted trajectory frames with NaN angles. "
                    "Inspect the matrix eigenspectrum and dihedral variance array."
                )

        return correlation_matrix

    def correlation_to_graph(self, correlation_matrix):
        G = nx.Graph()
        n = correlation_matrix.shape[0]
        for i in range(n):
            for j in range(i + 1, n):
                weight = abs(correlation_matrix[i, j])
                G.add_edge(i, j, weight=weight)
        return G

    def block_modularity_contributions(self, G, partition):
        # Computes the standard Newman–Girvan modularity contribution for community
        m2 = sum(w for _, _, w in G.edges(data="weight", default=1)) * 2
        communities = {}
        for node, comm in partition.items():
            communities.setdefault(comm, []).append(node)

        contributions = {}
        for comm_id, nodes in communities.items():
            node_set = set(nodes)
            # internal edges (e_c)
            e_c = sum(
                d.get("weight", 1)
                for u, v, d in G.edges(data=True)
                if u in node_set and v in node_set
            )
            # total degree of community nodes (a_c)
            a_c = sum(
                sum(d.get("weight", 1) for _, _, d in G.edges(u, data=True))
                for u in nodes
            )
            contributions[comm_id] = (e_c / m2) * 2 - (a_c / m2) ** 2

        return contributions

    def chose_correlated_bonds(self, correlation_matrix, rogue_corr_threshold=0.1):
        G = self.correlation_to_graph(correlation_matrix)
        partition = community_louvain.best_partition(
            G, weight="weight", resolution=1.0, random_state=self.seed
        )
        labels = np.array([partition[i] for i in range(len(partition))])
        modularity = self.block_modularity_contributions(G, partition)

        strong_blocks = []
        weak_blocks = []
        rogue_blocks = []

        for i in np.unique(labels):
            block = np.where(labels == i)[0]

            if len(block) == 1:
                rogue_blocks.extend(block.tolist())
                continue

            # Mean correlation of each dihedral to others in its block
            block_corr = np.mean(correlation_matrix[np.ix_(block, block)], axis=0)
            percentile = np.percentile(block_corr, 50)

            high_mask = block_corr >= percentile
            low_mask = ~high_mask

            high_group = block[high_mask]
            low_group = block[low_mask]

            # High group: only promote to Gibbs if genuinely correlated
            if len(high_group) > 1:
                mean_intra = np.mean(correlation_matrix[np.ix_(high_group, high_group)])
                if mean_intra >= rogue_corr_threshold:
                    strong_blocks.append(high_group)
                else:
                    # Correlation too weak even in the high half
                    rogue_blocks.extend(high_group.tolist())
            else:
                rogue_blocks.extend(high_group.tolist())

            # Low group: keep as weak block or demote to torsional
            if len(low_group) > 1:
                weak_blocks.append(low_group)
            else:
                rogue_blocks.extend(low_group.tolist())

        # Sanity check - all dihedrals must be accounted for
        all_assigned = (
            [d for b in strong_blocks for d in b]
            + [d for b in weak_blocks for d in b]
            + rogue_blocks
        )
        assert len(all_assigned) == correlation_matrix.shape[0], (
            f"Dihedral count mismatch: {len(all_assigned)} vs {correlation_matrix.shape[0]}"
        )
        assert len(set(all_assigned)) == len(all_assigned), (
            "Duplicate dihedral assignment detected"
        )

        return strong_blocks, weak_blocks, rogue_blocks, modularity

    def build_gibbs_blocks_from_trajectory(self, dcd_file):
        correlation_matrix = self.compute_dihedral_correlation_matrix(dcd_file)
        strong_blocks, weak_blocks, rogue_blocks, modularity = (
            self.chose_correlated_bonds(correlation_matrix)
        )
        strong_blocks_correlation = [
            np.mean(correlation_matrix[np.ix_(block, block)]) for block in strong_blocks
        ]
        weak_blocks_correlation = [
            np.mean(correlation_matrix[np.ix_(block, block)]) for block in weak_blocks
        ]

        print("=== Gibbs Blocks (high intra-correlation) ===")
        for i, block in enumerate(strong_blocks):
            mean_corr = strong_blocks_correlation[i]
            print(
                f"  Block {i}: {len(block)} dihedrals | "
                f"mean intra-corr: {mean_corr:.3f} | "
                f"modularity contrib: {modularity.get(i, float('nan')):.4f}"
            )

        print("\n=== Weak Blocks (low intra-correlation, still grouped) ===")
        for i, block in enumerate(weak_blocks):
            mean_corr = weak_blocks_correlation[i]
            print(
                f"  Weak block {i}: {len(block)} dihedrals | "
                f"mean intra-corr: {mean_corr:.3f}"
            )

        print("\n=== Rogue Blocks (low correlation, singletons) ===")
        print(f"  {len(rogue_blocks)} dihedrals")
        print(
            f"\nTotal: {sum(len(b) for b in strong_blocks)} Gibbs + "
            f"{sum(len(b) for b in weak_blocks)} weak + "
            f"{len(rogue_blocks)} rogue = "
            f"{sum(len(b) for b in strong_blocks) + sum(len(b) for b in weak_blocks) + len(rogue_blocks)}"
        )

        return (
            (strong_blocks, strong_blocks_correlation),
            (weak_blocks, weak_blocks_correlation),
            rogue_blocks,
        )
