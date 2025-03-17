import sys
sys.path.append("/home/laurentiu/git6/Robosample/bin/")
#sys.path.append("/home/pcuser/git6/Robosample/bin/")
import flexor
import mdtraj as md
import argparse
import robosample
#import batstat
import numpy as np

# Map the mobility types
mobilityMap = {
	"Cartesian": robosample.BondMobility.Translation,
	"Pin": robosample.BondMobility.Torsion,
	"Torsion": robosample.BondMobility.Torsion,
	"Slider": robosample.BondMobility.Slider,
}

# Read the flexibilities from a file
def getFlexibilitiesFromFile(flexFile):
	"""
	Backward compatibility: Reads the flexibilities from a file.
	:param flexFile: The file containing the flexibilities.
	:return: The list of flexibilities.
	"""
	flexibilities = []
	with open(flexFile, 'r') as f:
		for line in f:
			if line[0] == '#':
				continue
			tokens = line.split()
			if(len(tokens) >= 3):
				if(tokens[2] != "Weld"):
					aIx_1 = int(tokens[0])
					aIx_2 = int(tokens[1])
					mobility = mobilityMap[tokens[2]]
					flexibilities.append( robosample.BondFlexibility(aIx_1, aIx_2, mobility) )
	return flexibilities
#

# Print the flexibilities
def printFlexibilities(flexibilities):
	"""
	Prints the flexibilities.
	:param flexibilities: The list of flexibilities.
	"""
	for flexIx, flex in enumerate(flexibilities):
		print(flexibilities[flexIx].i, flexibilities[flexIx].j, flexibilities[flexIx].mobility)
#

#region: Parse the arguments
# python simulate.py baseName prmtop rst7 equil_steps prod_steps write_freq temperature_init seed
# python simulate.py 1a1p ../data/1a1p/1a1p.prmtop ../data/1a1p/1a1p.rst7 1000 10000 300.00 666
parser = argparse.ArgumentParser(description='Process PDB code and seed.')
parser.add_argument('--name', type=str, help='Name of the simulation.')
parser.add_argument('--top', type=str, help='Relative path to the .prmtop file.')
parser.add_argument('--rst7', type=str, help='Relative path to the .rst7 file.')
parser.add_argument('--equilSteps', type=int, help='The number of equilibration steps.')
parser.add_argument('--prodSteps', type=int, help='The number of production steps.')
parser.add_argument('--writeFreq', type=int, help='CSV and DCD write frequency.')
parser.add_argument('--baseTemperature', type=float, help='Temperature of the first replica.')
parser.add_argument('--runType', type=str, help='Run type: DEFAULT, REMC, RENEMC, RENE.')
parser.add_argument('--seed', type=int, help='The seed.')
parser.add_argument('--flexFNs', type=str, nargs='+', default=[], help='The flexFNs.')
args = parser.parse_args()
#endregion

# Create robosample context
run_type = getattr(robosample.RunType, args.runType)
context = robosample.Context(args.name, args.seed, 0, 1, run_type, 1, 0)

# Set parameters
context.setPdbRestartFreq(0) # WRITE_PDBS
context.setPrintFreq(args.writeFreq) # PRINT_FREQ
context.setNonbonded(0, 1.2)
context.setGBSA(0)
context.setVerbose(True)

# Load system
context.loadAmberSystem(args.top, args.rst7)

# Prepare flexor generator
mdtrajObj = md.load(args.rst7, top=args.top)
flexorObj = flexor.Flexor(mdtrajObj)

# Worlds
# region cpp
# void addWorld(
# 	bool fixmanTorque,
# 	int samplesPerRound,
# 	ROOT_MOBILITY rootMobility,
# 	const std::vector<BOND_FLEXIBILITY>& flexibilities,
# 	bool useOpenMM = true,
# 	bool visual = false,
# 	SimTK::Real visualizerFrequency = 0);
# endregion cpp
nofWorlds = 0

flexes_Cart = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
context.addWorld(False, 1, robosample.RootMobility.CARTESIAN, flexes_Cart, True, False, 0)
nofWorlds += 1

flexes_Rama = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
context.addWorld(False, 1, robosample.RootMobility.WELD, flexes_Rama, True, False, 0)
nofWorlds += 1

context.getWorld(1).setRollFlexibilities(True)

flexes_Rama = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
context.addWorld(False, 1, robosample.RootMobility.WELD, flexes_Rama, True, False, 0)
nofWorlds += 1

# Samplers
# region cpp
# bool addSampler(SamplerName samplerName,
# 	IntegratorType integratorType,
# 	ThermostatName thermostatName,
# 	bool useFixmanPotential);
# endregion cpp
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN
context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)

for wIx in range(1, nofWorlds):
	context.getWorld(wIx).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

# Replica exchange
nof_replicas = 2
temperatures = np.zeros(nof_replicas, dtype=np.float64)
boost_temperatures = np.zeros(nof_replicas, dtype=np.float64)
for replIx in range(nof_replicas):
    temperatures[replIx] = args.baseTemperature + (replIx * 50)
    boost_temperatures[replIx] = args.baseTemperature + (replIx * 50)  # used for openmm velocities

accept_reject_modes = nofWorlds * [robosample.AcceptRejectMode.MetropolisHastings]

timesteps = nofWorlds * [0.0007]
timesteps[1] = 0.1

worldIndexes = range(nofWorlds)
mdsteps = nofWorlds * [10]
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV] + ((nofWorlds - 1) * [robosample.IntegratorType.VERLET])

distort_options = nofWorlds * [0]
distort_args = nofWorlds * ["0"]
flow = nofWorlds * [0]
work = nofWorlds * [0]

for replIx in range(nof_replicas):
    context.addReplica(replIx)
    context.addThermodynamicState(replIx,
		temperatures[replIx],
		accept_reject_modes,
		distort_options,
		distort_args,
		flow,
		work,
		integrators,
		worldIndexes,
		timesteps,
		mdsteps)

# Initialize
context.Initialize()

# Run
context.RunREX(args.equilSteps, args.prodSteps)
