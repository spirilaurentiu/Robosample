import flexor
import mdtraj as md
import argparse
import robosample
import batstat
import numpy as np

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
c = robosample.Context(args.name, args.seed, 0, 1, robosample.RunType.REMC, 1, 0)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(args.write_freq) # PRINT_FREQ
c.setNonbonded(0, 1.2)
c.setGBSA(1)
c.setVerbose(False)

# load system
c.loadAmberSystem(args.prmtop, args.rst7)

# openmm cartesian
flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
c.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# Cluster the previous simulations
dcd_files = [f"{args.pdbid}_6000_{i}.dcd" for i in range(5)]
stats = batstat.BATCorrelations(args.prmtop, dcd_files)
corr = stats.compute_correlations()

abs_corr = np.abs(corr)
max_corr = np.max(abs_corr, axis=0)
clusters = flexor.cluster(max_corr)
for c in clusters:
	bond_list = []
	for aix1, aix2, dihedral_type in c:
		bond_list.append((aix1, aix2))
        
	flex = flexorObj.create_from_list(bond_list, robosample.BondMobility.Torsion)
	c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

c.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
for i in range(len(clusters)):
	c.getWorld(i).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

nof_replicas = 1
temperature = args.temperature_init
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] + [robosample.AcceptRejectMode.MetropolisHastings] * len(clusters)
timesteps = [0.0007] + [0.02] * len(clusters)
worldIndexes = [i for i in range(len(clusters))]
world_indexes = [i for i in range(len(clusters))]
mdsteps = [1429] + [10] * len(clusters) # 14286 - 1 ps instead of 10 ps
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV] + [robosample.IntegratorType.VERLET] * len(clusters)

distort_options = [0] * len(clusters)
distort_args = ["0"] * len(clusters)
flow = [0] * len(clusters)
work = [0] * len(clusters)

for i in range(nof_replicas):
    c.addReplica(i)
    c.addThermodynamicState(i,
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
c.Initialize()

# start the simulation
c.RunREX(args.equil_steps, args.prod_steps)