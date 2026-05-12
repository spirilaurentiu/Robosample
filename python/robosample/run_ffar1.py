import argparse
import time

import robosample

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

# def run_openmm_equilibration():
#     prmtop = app.AmberPrmtopFile(args.prmtop)
#     inpcrd = app.AmberInpcrdFile(args.inpcrd)
#     system = prmtop.createSystem(
#         nonbondedMethod=app.CutoffNonPeriodic,
#         nonbondedCutoff=1.2 * unit.nanometer,
#         implicitSolvent=app.OBC2,
#         constraints=None,
#         removeCMMotion=False,
#     )

#     for force in system.getForces():
#         if isinstance(force, mm.HarmonicBondForce):
#             force.setForceGroup(0)  # fast — evaluated every inner step
#         else:
#             force.setForceGroup(1)  # slow — evaluated every outer step

#     # 4 inner bond steps per 1 outer step → effective bond dt = 0.25 fs
#     timestep = 1 * unit.femtoseconds
#     integrator = mm.MTSIntegrator(
#         timestep,  # outer timestep = 1 fs
#         [(1, 1), (0, 4)],  # group 1 once, group 0 four times
#     )

#     platform = mm.Platform.getPlatformByName("CUDA")
#     simulation = app.Simulation(prmtop.topology, system, integrator, platform)
#     simulation.context.setPositions(inpcrd.positions)

#     simulation.context.setVelocitiesToTemperature(300 * unit.kelvin)

#     simulation.minimizeEnergy()

#     # ---- PRINT ENERGIES BEFORE STEP ----
#     state = simulation.context.getState(getEnergy=True)
#     print(f"Kinetic Energy: {state.getKineticEnergy()}")
#     print(f"Potential Energy: {state.getPotentialEnergy()}")
#     print(f"Total Energy: {state.getKineticEnergy() + state.getPotentialEnergy()}")

#     simulation.step(1)

#     # ---- PRINT ENERGIES BEFORE STEP ----
#     state = simulation.context.getState(getEnergy=True)
#     print(f"Kinetic Energy: {state.getKineticEnergy()}")
#     print(f"Potential Energy: {state.getPotentialEnergy()}")
#     print(f"Total Energy: {state.getKineticEnergy() + state.getPotentialEnergy()}")

#     # Reporters
#     steps_per_frame = 1000
#     output_dcd = f"{args.name}_equil_openmm.dcd"
#     simulation.reporters.append(app.DCDReporter(output_dcd, steps_per_frame))
#     simulation.reporters.append(
#         app.StateDataReporter(
#             stdout,
#             steps_per_frame,
#             step=True,
#             potentialEnergy=True,
#             kineticEnergy=True,
#             totalEnergy=True,
#         )
#     )

#     simulation.step(10000)

#     traj = md.load(output_dcd, top=args.prmtop)
#     last = traj[-1]

#     parm = pmd.load_file(args.prmtop)
#     parm.coordinates = last.xyz[0] * 10.0
#     if last.unitcell_lengths is not None:
#         parm.box = list(last.unitcell_lengths[0] * 10.0) + list(last.unitcell_angles[0])
#     parm.save("last_frame.rst7", format="rst7", overwrite=True)


# run_openmm_equilibration()

T0 = 300.0
T_MAX = 1000.0
NOF_REPLICAS = 1
R = 1 if NOF_REPLICAS == 1 else (T_MAX / T0) ** (1.0 / (NOF_REPLICAS - 1))

import mdtraj as md
import parmed as pmd

traj = md.load("vmd/ffar1_6000.repl0.dcd", top=args.prmtop)
last = traj[-1]
parm = pmd.load_file(args.prmtop)
parm.coordinates = last.xyz[0] * 10.0
if last.unitcell_lengths is not None:
    parm.box = list(last.unitcell_lengths[0] * 10.0) + list(last.unitcell_angles[0])
parm.save("last_frame.rst7", format="rst7", overwrite=True)

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd="last_frame.rst7",
    write_freq=args.write_freq,
    testing=False,
)

# First molecule is the target
context.set_root_mobility(0, robosample.rb.RootMobility.Free)

# Next two molecules are the outer proteins
context.set_root_mobility(1, robosample.rb.RootMobility.Weld)
context.set_root_mobility(2, robosample.rb.RootMobility.Weld)

# All other molecules are lipids and solvent
for mol_ix in range(3, context.get_num_molecules()):
    context.set_root_mobility(mol_ix, robosample.rb.RootMobility.Free)

# ############# CARTESIAN WORLD #############
# context.add_cartesian_world().add_sampler(
#     timeStep=0.001,
#     mdSteps=50_000,
#     boostMDSteps=50_000,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
#     use_nuts=False,
# )
# ############### CARTESIAN WORLD #############

############### World 2a ###############
resid = [103, 136, 173, 182, 190, 201, 232, 236, 243, 257, 275]
dihs = ["chi1", "chi2", "chi3", "chi4", "chi5"]

bonds = context.standard_dihedral_bonds.loc[
    (context.standard_dihedral_bonds["resid"].isin(resid))
    & (context.standard_dihedral_bonds["dihedral_type"].isin(dihs))
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
]

sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.005,  # 5 fs
    mdSteps=500,
    boostMDSteps=500,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
############### World 2a ###############

############### World 2b ###############
side_dihedral_types = [
    "chi1",
    "chi2",
    "chi3",
    "chi4",
    "chi5",
    "alpha",
    "beta",
    "gamma",
    "delta",
    "epsilon",
    "zeta",
    "sn1",
    "sn2",
    "lipid-chi1",
    "lipid-chi1_prime",
    "lipid-chi2",
    "lipid-chi3",
    "lipid-chi4",
    "lipid-chi5",
    "lipid-chi6",
    "lipid-chi7",
    "lipid-chi8",
    "lipid-chi9",
    "lipid-chi10",
    "lipid-chi11",
    "lipid-chi12",
    "lipid-chi13",
    "lipid-chi14",
    "lipid-chi15",
    "lipid-chi16",
    "lipid-chi17",
    "lipid-chi18",
]

bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(side_dihedral_types)
]
sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.01,  # 10 fs
    mdSteps=500,  # ignored if using NUTS
    boostMDSteps=500,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
############### World 2b ###############

################## World 3a ###############
bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"])
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
    & (context.standard_dihedral_bonds["resid"].isin(range(146, 182)))
]
sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.02,  # 20 fs
    mdSteps=300,  # ignored if using NUTS
    boostMDSteps=300,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
################## World 3a ###############

################### World 3b ###############
bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"])
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
    & (context.standard_dihedral_bonds["resid"].isin(range(209, 223)))
]
sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.02,  # 20 fs
    mdSteps=300,  # ignored if using NUTS
    boostMDSteps=300,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
################### World 3b ###############

################### World 3c ###############
bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"])
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
    & (context.standard_dihedral_bonds["resid"].isin(range(235, 242)))
]
sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.02,  # 20 fs
    mdSteps=300,  # ignored if using NUTS
    boostMDSteps=300,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
################### World 3c ###############

################### World 4 ###############
helix_residues = set(range(14, 37))  # TM1
helix_residues |= set(range(41, 65))  # TM2
helix_residues |= set(range(72, 110))  # TM3
helix_residues |= set(range(130, 147))  # TM4
helix_residues |= set(range(182, 210))  # TM5
helix_residues |= set(range(223, 248))  # TM6
helix_residues |= set(range(257, 290))  # TM7

bonds = context.standard_dihedral_bonds.loc[
    (context.standard_dihedral_bonds["dihedral_type"].isin(["phi", "psi"]))
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
    & (~context.standard_dihedral_bonds["resid"].isin(helix_residues))
]
sele = context.build_flexibilities(bonds)
context.add_torsional_world(sele).add_sampler(
    timeStep=0.01,  # 1 fs
    mdSteps=1000,  # ignored if using NUTS
    boostMDSteps=1000,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
################### World 4 ###############

############### World D - Diagnostic / Structural Reporters ###############
diag_resid = [14, 37, 41, 62, 86, 110, 130, 146, 182, 209, 223, 238, 247, 257, 290]

bonds = context.standard_dihedral_bonds.loc[
    (context.standard_dihedral_bonds["resid"].isin(diag_resid))
    & (context.standard_dihedral_bonds["dihedral_type"] == "phi")
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
]

sele = context.build_flexibilities(bonds)
world = context.add_torsional_world(sele, want_spatial_force_history=True).add_sampler(
    timeStep=0.02,  # 20 fs
    mdSteps=250,
    boostMDSteps=250,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)
############### World D - Diagnostic / Structural Reporters ###############

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
