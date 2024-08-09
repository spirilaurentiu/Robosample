import flexor
import mdtraj as md
import argparse
import robosample

# Create the parser
parser = argparse.ArgumentParser(description='Process PDB code and seed.')

# Add the arguments
# 1i42 1w4k 2btg 2eej 2kes 2kyy 2l39 2oa4 2obu 5xf0
parser.add_argument('prmtop', type=str, help='Prmtop file')
parser.add_argument('rst7', type=str, help='Rst7 file')
parser.add_argument('dcd', type=str, help='Output DCD path')
parser.add_argument('Seed', type=int, help='The seed')
parser.add_argument('EquilSteps', type=int, help='The number of MD equilibration steps')
parser.add_argument('ProdSteps', type=int, help='The number of MD production steps')
parser.add_argument('WriteFreq', type=int, help='CSV and DCD write frequency')

# Parse the arguments
args = parser.parse_args()

# Set the variables
seed = args.Seed

prmtop = args.prmtop
inpcrd = args.rst7
dcd = args.dcd

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
c = robosample.Context("placeholder_name", seed, 0, 1, robosample.RunType.REMC, 1, 0)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(args.WriteFreq) # PRINT_FREQ

# load system
c.loadAmberSystem(prmtop, inpcrd)

# # openmm cartesian
# flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
# c.addWorld(False, 1, robosample.RootMobility.CARTESIAN, flex, True, False, 0)

# rigid all
flex = []
c.addWorld(True, 1, robosample.RootMobility.CARTESIAN, flex, True, False, 0)

# samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

c.getWorld(0).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, False)
# c.getWorld(1).addSampler(sampler, robosample.IntegratorName.Verlet, thermostat, True)

nof_replicas = 1
temperature = 300
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

scale = 1
accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.001, 0.003]
world_indexes = [0, 1]
mdsteps = [25 * scale, 50 * scale]
boost_md_steps = [25 * scale, 50 * scale]
samples_per_round = [1, 1]

distort_options = [0, 0]
distort_args = ["0", "0"]
flow = [0, 0]
work = [0, 0]

# remove last element from all lists above
accept_reject_modes.pop()
timesteps.pop()
world_indexes.pop()
mdsteps.pop()
boost_md_steps.pop()
distort_options.pop()
distort_args.pop()
flow.pop()
work.pop()

for i in range(nof_replicas):
    c.addReplica(i)
    c.addThermodynamicState(i, temperatures[i], boost_temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, world_indexes, timesteps, mdsteps, boost_md_steps)

# initialize the simulation
c.initialize()

# start the simulation
c.RunREXNew(args.EquilSteps, args.ProdSteps)
