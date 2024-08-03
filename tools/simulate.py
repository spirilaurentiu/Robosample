import flexor
import mdtraj as md
import argparse
import robosample as robosample

# Create the parser
parser = argparse.ArgumentParser(description='Process PDB code and seed.')

# Add the arguments
# 1i42 1w4k 2btg 2eej 2kes 2kyy 2l39 2oa4 2obu 5xf0
parser.add_argument('PDBCode', type=str, help='The PDB code')
parser.add_argument('Seed', type=int, help='The seed')

# Parse the arguments
args = parser.parse_args()

# Set the variables
pdb_code = args.PDBCode
seed = args.Seed
prmtop = pdb_code + "/" + pdb_code + ".H.capped.prmtop"
inpcrd = pdb_code + "/" + pdb_code + ".H.capped.rst7"
dcd = pdb_code + "_" + str(seed) + ".dcd"

# prepare flexor generator
mdtrajObj = md.load(inpcrd, top=prmtop)
flexorObj = flexor.Flexor(mdtrajObj)

# create robosample context
c = robosample.Context(300, 300, 42, 0, 1, robosample.RunType.REMC, 1, 0)
c.setRequiredNofRounds(2)
c.setPdbRestartFreq(0) # WRITE_PDBS
c.setPrintFreq(1) # PRINT_FREQ

# load system
c.loadAmberSystem(prmtop, inpcrd)
c.appendDCDReporter(dcd)

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

nof_replicas = 1
temperature = 300
temperatures = []
boost_temperatures = []
for i in range(nof_replicas):
    temperatures.append(temperature + (i * 10))
    boost_temperatures.append(temperature + (i * 10))  # used for openmm velocities

scale = 1
accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings, robosample.AcceptRejectMode.MetropolisHastings]
timesteps = [0.001, 0.003, 0.004]
world_indexes = [0, 1, 2]
mdsteps = [25 * scale, 50 * scale, 50 * scale]
boost_md_steps = [25 * scale, 50 * scale, 50 * scale]
samples_per_round = [1, 1, 1]

distort_options = [0, 0, 0]
distort_args = ["0", "0" , "0"]
flow = [0, 0, 0]
work = [0, 0, 0]

for i in range(nof_replicas):
    c.addReplica(i)
    c.addThermodynamicState(i, temperatures[i], boost_temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, world_indexes, timesteps, mdsteps, boost_md_steps)

# initialize the simulation
c.initialize()

# start the simulation
c.Run()
