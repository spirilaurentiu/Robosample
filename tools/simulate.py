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
parser.add_argument('inpcrd', type=str, help='Relative path to the .inpcrd file.')
parser.add_argument('seed', type=int, help='The seed.')
parser.add_argument('equil_steps', type=int, help='The number of equilibration steps.')
parser.add_argument('prod_steps', type=int, help='The number of production steps.')
parser.add_argument('write_freq', type=int, help='CSV and DCD write frequency.')
parser.add_argument('temperature_init', type=int, help='Temperature of the first replica.')
parser.add_argument('pdbid', type=str, help='pdbid')
parser.add_argument('type', type=str, help='type of simulation')

# Parse the arguments
args = parser.parse_args()

# prepare flexor generator
mdtrajObj = md.load(args.inpcrd, top=args.prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
context = robosample.Context(args.name, args.seed, 0, 1, robosample.RunType.REMC, 1, 0)
context.setPdbRestartFreq(0) # WRITE_PDBS
context.setPrintFreq(args.write_freq) # PRINT_FREQ
context.setNonbonded(0, 1.2)
context.setGBSA(1)
context.setVerbose(False)

# load system
context.loadAmberSystem(args.prmtop, args.inpcrd)

if args.type == 'tdnr':
	# Do torsional dynamics - non-redundant dihedrals
	stats = batstat.BATCorrelations(args.prmtop, args.inpcrd)
	atom_indices = stats.get_dihedral_atom_indices()

	flex = flexorObj.create_from_list(atom_indices, robosample.BondMobility.Torsion)
	context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)
	context.getWorld(0).addSampler(robosample.SamplerName.HMC, robosample.IntegratorType.OMMVV, robosample.ThermostatName.ANDERSEN, False)

	nof_replicas = 1
	temperature = args.temperature_init
	temperatures = []
	boost_temperatures = []
	for i in range(nof_replicas):
		temperatures.append(temperature + (i * 10))
		boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

	accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings]
	timesteps = [0.0075]
	worldIndexes = [0]
	world_indexes = [0]
	mdsteps = [10]
	boost_md_steps = mdsteps
	integrators = [robosample.IntegratorType.VERLET]

	distort_options = [0]
	distort_args = ["0"]
	flow = [0]
	work = [0]

	for i in range(nof_replicas):
		context.addReplica(i)
		context.addThermodynamicState(i, temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

	context.Initialize()
	context.RunREX(args.equil_steps, args.prod_steps)

elif args.type == 'tdc':
	# Do torsional dynamics - correlated dihedrals

	# Cluster the previous simulations
	stats = batstat.BATCorrelations(args.prmtop, args.inpcrd)
	correlation_file = f"{args.pdbid}_6000_correlation.npy"

	if os.path.exists(correlation_file):
		corr = np.load(correlation_file)
	else:
		dcd_files = [f"{args.pdbid}_{i}.dcd" for i in range(5)]
		stats.compute_dihedrals_from_dcd(dcd_files)
		corr = stats.compute_correlations()
		np.save(correlation_file, corr)

	# Partition the dihedrals into blocks
	blocks, samples_per_round = stats.dynamic_partitioning(np.mean(np.abs(corr), axis=0))
	print("samples_per_round", samples_per_round)

	# Create the flexors from the blocks
	for block, num_samples in zip(blocks, samples_per_round):
		flex = flexorObj.create_from_list(block, robosample.BondMobility.Torsion)
		context.addWorld(True, num_samples, robosample.RootMobility.WELD, flex, True, False, 0)

	for i in range(len(blocks)):
		context.getWorld(i).addSampler(robosample.SamplerName.HMC, robosample.IntegratorType.VERLET, robosample.ThermostatName.ANDERSEN, True)

	nof_replicas = 1
	temperature = args.temperature_init
	temperatures = []
	boost_temperatures = []
	for i in range(nof_replicas):
		temperatures.append(temperature + (i * 10))
		boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

	accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] * (len(blocks))
	timesteps = [0.0075] * len(blocks)
	worldIndexes = range(len(blocks))
	world_indexes = range(len(blocks))
	mdsteps = [10] * len(blocks)
	boost_md_steps = mdsteps
	integrators = [robosample.IntegratorType.VERLET] * len(blocks)

	distort_options = [0] * (len(blocks))
	distort_args = ["0"] * (len(blocks))
	flow = [0] * (len(blocks))
	work = [0] * (len(blocks))

	for i in range(nof_replicas):
		context.addReplica(i)
		context.addThermodynamicState(i, temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

	context.Initialize()
	context.RunREX(args.equil_steps, args.prod_steps)

else:
	# # sidechains pins
	# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
	# c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

	# # ramachandran pins
	# flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
	# c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)



	# c.getWorld(1).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
	# c.getWorld(2).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

	pass


# python3 simulate.py 1APQ_clusters_0 ../../robocath/data-raw/1APQ.prmtop ../../robocath/data-raw/1APQ_min.inpcrd 6000 10 100 10 300 1APQ tdnr