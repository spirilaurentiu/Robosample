import flexor
import mdtraj as md
import argparse
import robosample

# Create the parser
parser = argparse.ArgumentParser(description='Process PDB code and seed.')

# Add the arguments
# 1i42 1w4k 2btg 2eej 2kes 2kyy 2l39 2oa4 2obu 5xf0
parser.add_argument('PDBCode', type=str, help='The PDB code')
parser.add_argument('Seed', type=int, help='The seed')
parser.add_argument('EquilSteps', type=int, help='The number of MD equilibration steps')
parser.add_argument('ProdSteps', type=int, help='The number of MD production steps')
parser.add_argument('WriteFreq', type=int, help='CSV and DCD write frequency')

# Parse the arguments
args = parser.parse_args()

# Set the variables
pdb_code = args.PDBCode
seed = args.Seed

prmtop = '../datasets/diverse-cath-nmr/' + pdb_code + ".prmtop"
inpcrd = '../datasets/diverse-cath-nmr/' + pdb_code + "_min.rst7"

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
c = robosample.Context(pdb_code, seed, 0, 1, robosample.RunType.REMC, 1, 0)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(args.WriteFreq) # PRINT_FREQ
c.setNonbonded(0, 1.2)
c.setGBSA(1)

# load system
c.loadAmberSystem(prmtop, inpcrd)

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
temperature = 300
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.0007, 0.035, 0.075]
worldIndexes = [0, 1, 2]
world_indexes = [0, 1, 2]
mdsteps = [14286, 20, 40]
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
c.RunREX(args.EquilSteps, args.ProdSteps)