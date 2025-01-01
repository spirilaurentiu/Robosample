import flexor
import mdtraj as md
import argparse
import robosample

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

# Parse the arguments
args = parser.parse_args()

# prepare flexor generator
mdtrajObj = md.load(args.rst7, top=args.prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
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

# sidechains pins
flex = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# ramachandran pins
flex = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
c.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

# samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

c.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
c.getWorld(1).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
c.getWorld(2).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

nof_replicas = 1
temperature = args.temperature_init
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.0007, 0.035, 0.075]
worldIndexes = [0, 1, 2]
world_indexes = [0, 1, 2]
mdsteps = [1429, 10, 10] # 14286 - 1 ps instead of 10 ps
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV, robosample.IntegratorType.VERLET, robosample.IntegratorType.VERLET]

distort_options = [0, 0, 0]
distort_args = ["0", "0" , "0"]
flow = [0, 0, 0]
work = [0, 0, 0]

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