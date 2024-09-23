import flexor
import mdtraj as md
import argparse
import robosample

# Create the parser
parser = argparse.ArgumentParser(description='')
parser.add_argument('name', type=str, help='The name of the simulation')
parser.add_argument('prmtop', type=str, help='AMBER topology file')
parser.add_argument('rst7', type=str, help='AMBER coordinate file')
parser.add_argument('seed', type=int, help='The seed')
parser.add_argument('equil_steps', type=int, help='The number of MD equilibration steps')
parser.add_argument('prod_steps', type=int, help='The number of MD production steps')
parser.add_argument('write_freq', type=int, help='CSV and DCD write frequency')
parser.add_argument('num_replicas', type=int, help='Total number of replicas')
parser.add_argument('t_min', type=int, help='Initial temperature')
parser.add_argument('t_max', type=int, help='Final temperature')
# parser.add_argument('rigid_mdsteps', type=int, help='Rigid world MD steps')
args = parser.parse_args()

# create robosample context
c = robosample.Context(args.name, args.seed, 0, 1, robosample.RunType.REMC, 1, 0)
c.setPdbRestartFreq(0) # do not write PDBs
c.setPrintFreq(args.write_freq) # how often to write CSV and DCD files
c.setGBSA(1) # water - implicit solvent

# load system
c.loadAmberSystem(args.prmtop, args.rst7)

# prepare flexor generator
mdtrajObj = md.load(args.rst7, top=args.prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

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

c.getWorld(0).addSampler(sampler, robosample.IntegratorName.OMMVV, thermostat, False)
c.getWorld(1).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)
c.getWorld(2).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)

# Sampler data
timesteps = [0.001, 0.005, 0.01]
world_indexes = [0, 1, 2]
mdsteps = boost_md_steps = [25, 10, 10]
distort_options = [0, 0, 0]
distort_args = ["0", "0", "0"]
flow = [0, 0, 0]
work = [0, 0, 0]
accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] * len(world_indexes)

# Calculate the common ratio for geometric progression
if args.num_replicas > 1:
    r = (args.t_max / args.t_min) ** (1 / (args.num_replicas - 1))
else:
    r = 1

# Generate the temperatures
temperatures = [args.t_min * r ** i for i in range(args.num_replicas)]
boost_temperatures = [args.t_min * r ** i for i in range(args.num_replicas)]  # used for openmm velocities

for i in range(args.num_replicas):
    c.addReplica(i)
    c.addThermodynamicState(i, temperatures[i], boost_temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, world_indexes, timesteps, mdsteps, boost_md_steps)

# initialize the simulation
c.initialize()

# start the simulation
c.RunREXNew(args.equil_steps, args.prod_steps)

# cd build/ && python3 simulate.py rage.fps.0 rage.fps.0.prmtop rage.fps.0.min.rst7 42 1 999 1 1 300 600 && python3 simulate.py rage.fps.4 rage.fps.4.prmtop rage.fps.4.min.rst7 42 1 499 1 1 300 600
# cd build/ && python3 simulate.py rage.fps.1 rage.fps.1.prmtop rage.fps.1.min.rst7 42 1 999 1 1 300 600 && python3 simulate.py rage.fps.5 rage.fps.5.prmtop rage.fps.5.min.rst7 42 1 499 1 1 300 600
# cd build/ && python3 simulate.py rage.fps.2 rage.fps.2.prmtop rage.fps.2.min.rst7 42 1 999 1 1 300 600 && python3 simulate.py rage.fps.6 rage.fps.6.prmtop rage.fps.6.min.rst7 42 1 499 1 1 300 600
# cd build/ && python3 simulate.py rage.fps.3 rage.fps.3.prmtop rage.fps.3.min.rst7 42 1 999 1 1 300 600 && python3 simulate.py rage.fps.7 rage.fps.7.prmtop rage.fps.7.min.rst7 42 1 499 1 1 300 600