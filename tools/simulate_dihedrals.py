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

# Add the dihedral flexor
dcd_files = [f"{args.pdbid}_{i}.dcd" for i in range(5)]
stats = batstat.BATCorrelations(dcd_files, args.prmtop)
flex = flexorObj.create_from_list(stats.get_dihedral_atom_indices(), robosample.BondMobility.Torsion)
context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)
context.getWorld(0).addSampler(robosample.SamplerName.HMC, robosample.IntegratorType.VERLET, robosample.ThermostatName.ANDERSEN, True)

# Add replicas
nof_replicas = 1
temperature = args.temperature_init
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.007]
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

# initialize the simulation
context.Initialize()

# start the simulation
context.RunREX(args.equil_steps, args.prod_steps)
