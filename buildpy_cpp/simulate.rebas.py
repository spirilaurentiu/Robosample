import sys
sys.path.append("/home/laurentiu/git6/Robosample/bin")
import flexor
import mdtraj as md
import argparse
import robosample
#import batstat
import numpy as np


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

# Openmm Cartesian
nofWorlds = 0

print("ADD WORLD Adding Cartesian world...", end = ' ', flush=True)
flexes_Cart = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
context.addWorld(False, 1, robosample.RootMobility.CARTESIAN, flexes_Cart, True, False, 0)
nofWorlds += 1
print("    done.", flush=True)


# Add worlds
for flexFNIx, flexFN in enumerate(args.flexFNs):
	print("ADD WORLD Adding world from file: ", flexFN, flush=True)	
	flexibilities = getFlexibilitiesFromFile(flexFN)
	printFlexibilities(flexibilities)
	context.addWorld(False, 1, robosample.RootMobility.WELD, flexibilities, True, False, 0)
	nofWorlds += 1
	print("    done.", flush=True)

# Samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN
print("ADD SAMPLER Adding OMMVV sampler...", end = ' ', flush=True)
context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
print("    done.", flush=True)

for worldIx in range(1, nofWorlds):
	print("ADD SAMPLER Adding Verlet sampler...", end = ' ', flush=True)
	context.getWorld(worldIx).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
	print("    done.", flush=True)

# Replica exchange
nof_replicas = 2
temperatures = np.zeros(nof_replicas, dtype=np.float64)
boost_temperatures = np.zeros(nof_replicas, dtype=np.float64)
for replIx in range(nof_replicas):
    temperatures[replIx] = args.baseTemperature + (replIx * 50)
    boost_temperatures[replIx] = args.baseTemperature + (replIx * 50)  # used for openmm velocities

accept_reject_modes = nofWorlds * [robosample.AcceptRejectMode.MetropolisHastings]

worldIndexes = range(nofWorlds)
timesteps = nofWorlds * [0.0007]
mdsteps = nofWorlds * [10]
boost_md_steps = mdsteps

integrators = [robosample.IntegratorType.OMMVV] + ((nofWorlds - 1) * [robosample.IntegratorType.VERLET])

distort_options = nofWorlds * [0]
distort_args = nofWorlds * ["0"]
flow = nofWorlds * [0]
work = nofWorlds * [0]

distort_options[-1] = -3
accept_reject_modes[-1] = robosample.AcceptRejectMode.AlwaysAccept
mdsteps[-1] = 0


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
