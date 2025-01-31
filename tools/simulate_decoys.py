import flexor
import mdtraj as md
import argparse
import robosample as robosample
import batstat
import numpy as np
import os

# Create the parser
parser = argparse.ArgumentParser(description='Process PDB code and seed.')

# Add the arguments
parser.add_argument('name', type=str, help='Name of the simulation.')
parser.add_argument('prmtop', type=str, help='Relative path to the .prmtop file.')
parser.add_argument('rst7', type=str, help='Relative path to the .rst7 file.')
parser.add_argument('seed', type=int, help='The seed.')
parser.add_argument('equil_steps', type=int, help='The number of equilibration steps.')
parser.add_argument('prod_steps', type=int, help='The number of production steps.')
parser.add_argument('write_freq', type=int, help='CSV and DCD write frequency.')
parser.add_argument('temperature_init', type=int, help='Temperature of the first replica.')
parser.add_argument('pdbid', type=str, help='pdbid')

# Parse the arguments
args = parser.parse_args()

# prepare flexor generator
mdtrajObj = md.load(args.rst7, top=args.prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
# PDBID: 1APQ -> 1APQ_6000_0.dcd, 1APQ_6000_1.dcd, 1APQ_6000_2.dcd
context = robosample.Context(args.name, args.seed, 0, 1, robosample.RunType.REMC, 1, 0)
context.setPdbRestartFreq(0) # WRITE_PDBS
context.setPrintFreq(args.write_freq) # PRINT_FREQ
context.setNonbonded(0, 1.2)
context.setGBSA(1)
context.setVerbose(False)

# load system
context.loadAmberSystem(args.prmtop, args.rst7)

# openmm cartesian
flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
context.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# Cluster the previous simulations
dcd_files = [f"{args.pdbid}_{i}.dcd" for i in range(5)]
stats = batstat.BATCorrelations(dcd_files, args.prmtop)

corr_path = f"{args.pdbid}_correlations.npy"
if os.path.exists(corr_path):
	corr = np.load(corr_path)
else:
	corr = stats.compute_correlations()
	np.save(corr_path, corr)

abs_corr = np.abs(corr)
max_corr = np.max(abs_corr, axis=0)
clusters = stats.cluster(max_corr)
clusters = [clusters[0]] # only one cluster
for cluster in clusters:
	bond_list = []
	for aix1, aix2, dihedral_type in cluster:
		bond_list.append((aix1, aix2))
        
	decoy_bonds = stats.chose_decoy_bonds(bond_list, steps=3000)
	flex = flexorObj.create_from_list(decoy_bonds, robosample.BondMobility.Torsion)
	context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)
	
# samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
for i in range(len(clusters)):
	context.getWorld(i + 1).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

nof_replicas = 1
temperature = args.temperature_init
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] + [robosample.AcceptRejectMode.MetropolisHastings] * len(clusters)
timesteps = [0.0007] + [0.005] * len(clusters)
worldIndexes = [0] + [i+1 for i in range(len(clusters))]
world_indexes = [0] + [i+1 for i in range(len(clusters))]
mdsteps = [1429] + [10] * len(clusters) # 14286 - 1 ps instead of 10 ps
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV] + [robosample.IntegratorType.VERLET] * len(clusters)

distort_options = [0] + [0] * len(clusters)
distort_args = ["0"] + ["0"] * len(clusters)
flow = [0] + [0] * len(clusters)
work = [0] + [0] * len(clusters)

for i in range(nof_replicas):
    context.addReplica(i)
    context.addThermodynamicState(i,
		temperatures[i],
		accept_reject_modes,
		distort_options,
		distort_args,
		flow,
		work,
        integrators,
        worldIndexes,
		timesteps,
		mdsteps)

# initialize the simulation
context.Initialize()

# start the simulation
context.RunREX(args.equil_steps, args.prod_steps)
