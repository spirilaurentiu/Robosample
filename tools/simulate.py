import flexor
import mdtraj as md
import argparse
import robosample
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
stats = batstat.BATCorrelations(dcd_files, args.prmtop, args.rst7)
correlation_file = f"{args.pdbid}_6000_correlation.npy"

print(correlation_file)

if os.path.exists(correlation_file):
	corr = np.load(correlation_file)
else:
	corr = stats.compute_correlations()
	np.save(correlation_file, corr)

blocks, collapsed = stats.dynamic_partitioning(np.mean(np.abs(corr), axis=0))
blocks.append(collapsed)

samples_per_round = [10] * len(blocks)
samples_per_round.append(1)

print(len(blocks))

for block, num_samples in zip(blocks, samples_per_round):
	bond_list = []
	print(block)
	for (aix1, aix2) in block:
		bond_list.append((aix1, aix2))
        
	flex = flexorObj.create_from_list(bond_list, robosample.BondMobility.Torsion)
	context.addWorld(True, num_samples, robosample.RootMobility.WELD, flex, True, False, 0)

# # sidechains pins
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
# c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# # ramachandran pins
# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
# c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)

# c.getWorld(1).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
# c.getWorld(2).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

for i in range(len(blocks) + 1):
	context.getWorld(i).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

nof_replicas = 1
temperature = args.temperature_init
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.0007] + [0.0075] * len(blocks)
worldIndexes = range(len(blocks) + 1)
world_indexes = range(len(blocks) + 1)
mdsteps = [1429] + [10] * len(blocks) # 14286 - 1 ps instead of 10 ps
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV, robosample.IntegratorType.VERLET, robosample.IntegratorType.VERLET]

distort_options = [0] * (len(blocks) + 1)
distort_args = ["0"] * (len(blocks) + 1)
flow = [0] * (len(blocks) + 1)
work = [0] * (len(blocks) + 1)

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


# python3 simulate.py 1APQ_clusters_0 ../../robocath/data-raw/1APQ.prmtop ../../robocath/data-raw/1APQ_min.rst7 6000 10 100 10 300 1APQ