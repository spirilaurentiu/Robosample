import argparse
import time
from sys import stdout

import mdtraj as md
import parmed as pmd

import openmm as mm
import robosample
from openmm import app, unit

# python3 python/robosample/run_ffar1.py ffar1 examples/ffar1.prmtop examples/ffar1.rst7 6000 0 20 1

# Create the parser
parser = argparse.ArgumentParser(description="Process PDB code and seed.")

# Add the arguments
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="The number of equilibration steps.")
parser.add_argument("prod_steps", type=int, help="The number of production steps.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")

# Parse the arguments
args = parser.parse_args()


def run_openmm_equilibration():
    prmtop = app.AmberPrmtopFile(args.prmtop)
    inpcrd = app.AmberInpcrdFile(args.inpcrd)
    system = prmtop.createSystem(
        nonbondedMethod=app.CutoffNonPeriodic,
        nonbondedCutoff=1.2 * unit.nanometer,
        implicitSolvent=app.OBC2,
        constraints=None,
        removeCMMotion=False,
    )

    for force in system.getForces():
        if isinstance(force, mm.HarmonicBondForce):
            force.setForceGroup(0)  # fast — evaluated every inner step
        else:
            force.setForceGroup(1)  # slow — evaluated every outer step

    # 4 inner bond steps per 1 outer step → effective bond dt = 0.25 fs
    timestep = 1 * unit.femtoseconds
    integrator = mm.MTSIntegrator(
        timestep,  # outer timestep = 1 fs
        [(1, 1), (0, 4)],  # group 1 once, group 0 four times
    )

    platform = mm.Platform.getPlatformByName("CUDA")
    simulation = app.Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)

    simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

    simulation.minimizeEnergy()

    # ---- PRINT ENERGIES BEFORE STEP ----
    state = simulation.context.getState(getEnergy=True)
    print(f"Kinetic Energy: {state.getKineticEnergy()}")
    print(f"Potential Energy: {state.getPotentialEnergy()}")
    print(f"Total Energy: {state.getKineticEnergy() + state.getPotentialEnergy()}")

    simulation.step(1)

    # ---- PRINT ENERGIES BEFORE STEP ----
    state = simulation.context.getState(getEnergy=True)
    print(f"Kinetic Energy: {state.getKineticEnergy()}")
    print(f"Potential Energy: {state.getPotentialEnergy()}")
    print(f"Total Energy: {state.getKineticEnergy() + state.getPotentialEnergy()}")

    # Reporters
    steps_per_frame = 1000
    output_dcd = f"{args.name}_equil_openmm.dcd"
    simulation.reporters.append(app.DCDReporter(output_dcd, steps_per_frame))
    simulation.reporters.append(
        app.StateDataReporter(
            stdout,
            steps_per_frame,
            step=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
        )
    )

    simulation.step(10000)

    traj = md.load(output_dcd, top=args.prmtop)
    last = traj[-1]

    parm = pmd.load_file(args.prmtop)
    parm.coordinates = last.xyz[0] * 10.0
    if last.unitcell_lengths is not None:
        parm.box = list(last.unitcell_lengths[0] * 10.0) + list(last.unitcell_angles[0])
    parm.save("last_frame.rst7", format="rst7", overwrite=True)


# run_openmm_equilibration()

T0 = 300.0
T_MAX = 1000.0
NOF_REPLICAS = 1
R = 1 if NOF_REPLICAS == 1 else (T_MAX / T0) ** (1.0 / (NOF_REPLICAS - 1))

TIMESTEP_CARTESIAN = 0.001
MDSTEPS_CARTESIAN = 500  # 500 fs

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# [[rb.BondFlexibility(), rb.BondFlexibility(), ...], [...], ...]

# context.addCartesianWorld().add_sampler(
#     timeStep=TIMESTEP_CARTESIAN,
#     mdSteps=MDSTEPS_CARTESIAN,
#     boostMDSteps=MDSTEPS_CARTESIAN,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=False,
# )

bonds = context.find_bonds(
    [
        (3168, 3171),
        (3780, 3770),
        (2566, 2586),
        (1140, 1138),
        (1668, 1653),
        (533, 548),
    ]
)
sele = context.build_flexibilities(bonds)
context.addTorsionalWorld(sele).add_sampler(
    timeStep=0.1,
    mdSteps=100,  # ignored if using NUTS
    boostMDSteps=100,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    use_nuts=False,
)

# mask_phi = context.standard_dihedral_bonds["dihedral_type"] == "phi"
# mask_psi = context.standard_dihedral_bonds["dihedral_type"] == "psi"
# mask_target_mol_ix = context.standard_dihedral_bonds["molecule_index"] == 0
# mask_base = mask_target_mol_ix

# # N terminus
# mask_n_terminus = context.standard_dihedral_bonds["resid"].between(290, 299)
# bonds = context.standard_dihedral_bonds[mask_n_terminus & mask_target_mol_ix]
# sele = context.build_flexibilities(bonds)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=0.01,
#     mdSteps=10,  # ignored if using NUTS
#     boostMDSteps=10,  # ignored if using NUTS
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=True,
# )

# # C terminus
# mask_c_terminus = context.standard_dihedral_bonds["resid"].between(290, 299)
# bonds = context.standard_dihedral_bonds[mask_c_terminus]
# sele = context.build_flexibilities(bonds)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=0.01,
#     mdSteps=2000,  # ignored if using NUTS
#     boostMDSteps=2000,  # ignored if using NUTS
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=True,
# )

##################### SIDECHAIN ##############################
# chi_types = ["chi1", "chi2", "chi3", "chi4", "chi5"]
# mask_target_chis = context.standard_dihedral_bonds["dihedral_type"].isin(chi_types)
# mask_not_cyx = context.standard_dihedral_bonds["resname"] != "CYX"
# mask_sidechain = mask_target_chis & mask_not_cyx

# bonds = context.standard_dihedral_bonds[mask_sidechain]
# sele = context.build_flexibilities(bonds)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=0.01,
#     mdSteps=1000,  # ignored if using NUTS
#     boostMDSteps=1000,  # ignored if using NUTS
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=True,
# )


# chi_types = ["chi1", "chi2", "chi3", "chi4", "chi5"]
# df = context.standard_dihedral_bonds

# mask_chi = df["dihedral_type"].isin(chi_types)
# mask_not_cyx = df["resname"] != "CYX"
# mask_base = mask_chi & mask_not_cyx

# # Light terminals: small rotating group, low barrier, fast oscillation
# WORLD1_RESIDUES = {"SER", "THR", "CYS", "CYM"}
# WORLD1_STEP = 0.005  # Tight wells on OG/OG1/SG need small steps

# # Flexible chains: long polar/charged side chains, medium barriers
# WORLD2_RESIDUES = {
#     "ASN",
#     "ASP",
#     "ASH",
#     "GLN",
#     "GLU",
#     "GLH",
#     "LYS",
#     "LYN",
#     "ARG",
#     "MET",
#     "PRO",
# }
# WORLD2_STEP = 0.01

# # Branched/aromatic: heavy rotating groups, high barriers (~3-5 kcal/mol)
# WORLD3_RESIDUES = {"VAL", "LEU", "ILE", "PHE", "TYR", "TRP", "HID", "HIE", "HIP"}
# WORLD3_STEP = 0.03

# # Side chains must appear in one side chain world
# mask_w1 = mask_base & df["resname"].isin(WORLD1_RESIDUES)
# mask_w2 = mask_base & df["resname"].isin(WORLD2_RESIDUES)
# mask_w3 = mask_base & df["resname"].isin(WORLD3_RESIDUES)

# uncovered = mask_base & ~(mask_w1 | mask_w2 | mask_w3)
# if uncovered.any():
#     missing = df[uncovered]["resname"].unique().tolist()
#     raise ValueError(f"Residues not assigned to any torsional world: {missing}")

# overlap = (mask_w1 & mask_w2) | (mask_w1 & mask_w3) | (mask_w2 & mask_w3)
# if overlap.any():
#     raise ValueError("Residue appears in more than one torsional world")

# # Add side chain worlds
# for mask, step, label in [
#     (mask_w1, WORLD1_STEP, "light-terminals"),
#     (mask_w2, WORLD2_STEP, "flexible-chains"),
#     (mask_w3, WORLD3_STEP, "branched-aromatic"),
# ]:
#     bonds = df[mask]
#     sele = context.build_flexibilities(bonds)
#     context.addTorsionalWorld(sele).add_sampler(
#         timeStep=step,
#         mdSteps=1000,  # Ignored if using NUTS
#         boostMDSteps=1000,  # Ignored if using NUTS
#         acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#         use_nuts=True,
#     )
##################### SIDECHAIN ##############################

# # Add torsional world with standardized dihedrals
# sele = context.build_flexibilities(context.standard_dihedral_bonds)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=TIMESTEP_TD,
#     mdSteps=MDSTEPS_TD,
#     boostMDSteps=MDSTEPS_TD,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
# )

# mask_phi = context.standard_dihedral_bonds["dihedral_type"] == "phi"
# mask_psi = context.standard_dihedral_bonds["dihedral_type"] == "psi"
# mask_disorder = context.standard_dihedral_bonds["dss"] == "C"
# mask_target_mol_ix = context.standard_dihedral_bonds["molecule_index"] == 0
# coil_bonds = context.standard_dihedral_bonds[
#     (mask_phi | mask_psi) & mask_disorder & mask_target_mol_ix
# ]

# mask_proline = context.standard_dihedral_bonds["resname"] == "PRO"
# proline_bonds = context.standard_dihedral_bonds[
#     (mask_phi | mask_psi) & mask_proline & mask_target_mol_ix
# ]

# # mask_glycine = context.standard_dihedral_bonds["resname"] == "GLY"
# # glycine_bonds = context.standard_dihedral_bonds[
# #     (mask_phi | mask_psi) & mask_glycine & mask_target_mol_ix
# # ]

# mask_alpha = context.standard_dihedral_bonds["dihedral_type"] == "alpha"
# headgroup_bonds = context.standard_dihedral_bonds[mask_alpha]

# mask_sn1 = context.standard_dihedral_bonds["dihedral_type"] == "sn1"
# mask_sn2 = context.standard_dihedral_bonds["dihedral_type"] == "sn2"
# ester_bonds = context.standard_dihedral_bonds[mask_sn1 | mask_sn2]

# mask_lipid_chi2 = context.standard_dihedral_bonds["dihedral_type"] == "lipid-chi2"
# mask_lipid_chi7 = context.standard_dihedral_bonds["dihedral_type"] == "lipid-chi7"
# mask_lipid_chi12 = context.standard_dihedral_bonds["dihedral_type"] == "lipid-chi12"
# mask_lipid_chi17 = context.standard_dihedral_bonds["dihedral_type"] == "lipid-chi17"
# acyl_chain_bonds = context.standard_dihedral_bonds[
#     mask_lipid_chi2 | mask_lipid_chi7 | mask_lipid_chi12 | mask_lipid_chi17
# ]

# combined_selection = pd.concat(
#     [coil_bonds, proline_bonds, headgroup_bonds, ester_bonds, acyl_chain_bonds], axis=0
# ).drop_duplicates()

# sele = context.build_flexibilities(combined_selection)
# context.addTorsionalWorld(sele).add_sampler(
#     timeStep=TIMESTEP_TD,
#     mdSteps=MDSTEPS_TD,
#     boostMDSteps=MDSTEPS_TD,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=False,
# )

# Add replicas (geometric temperature ladder)
temperatures = []
for i in range(NOF_REPLICAS):
    temperatures.append(T0 * (R**i))

context.initialize(temperatures)

# Run the simulation
start_time = time.perf_counter()

context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)

end_time = time.perf_counter()
duration = end_time - start_time

print(f"run_rex() took {duration:.4f} seconds")

# # Calculate RMSD
# traj = md.load(f"{args.name}_{args.seed}.repl0.dcd", top=args.prmtop)
# ref = md.load("last_frame.rst7", top=args.prmtop)
# rmsd = md.rmsd(traj, ref)
# print(f"RMSD to reference: {np.max(rmsd)}")

"""
source leaprc.protein.ff19SB
mol1 = sequence { ACE ALA ALA NME }
mol2 = sequence { ACE ALA ALA NME }
translate mol2 { 15.0 0.0 0.0 }
system = combine { mol1 mol2 }
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.inpcrd
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.rst7
savepdb system 2ala.2ala.pdb
quit
"""

"""
source leaprc.protein.ff19SB
system = sequence { ACE GLY ALA VAL LEU ILE SER THR CYS MET PHE TYR TRP PRO ASP GLU ASN GLN HIS LYS ARG NME }
saveamberparm system all.prmtop all.inpcrd
saveamberparm system all.prmtop all.rst7
savepdb system all.pdb
quit
"""
