#region: Imports
import sys
sys.path.append("/home/laurentiu/git6/Robosample/bin")
import flexor
import mdtraj as md
import argparse
import robosample
#import batstat
import numpy as np
#endregion

mobilityMap = {
	"Cartesian": robosample.BondMobility.Translation,
	"Pin": robosample.BondMobility.Torsion,
	"Torsion": robosample.BondMobility.Torsion,
	"Slider": robosample.BondMobility.Slider,
	"AnglePin": robosample.BondMobility.AnglePin,
	"BendStretch": robosample.BondMobility.BendStretch,
	"Spherical": robosample.BondMobility.Spherical,
	"Cylinder": robosample.BondMobility.Cylinder,
	"OrthoSpherical": robosample.BondMobility.OrthoSpherical
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
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')
	
parser = argparse.ArgumentParser(description='Process PDB code and seed.')
parser.add_argument('--name', type=str, help='Name of the simulation.')
parser.add_argument('--top', type=str, help='Relative path to the .prmtop file.')
parser.add_argument('--rst7', type=str, help='Relative path to the .rst7 file.')
parser.add_argument('--rstDir', type=str, help='Restart directory.')
parser.add_argument('--equilSteps', type=int, help='The number of equilibration steps.')
parser.add_argument('--prodSteps', type=int, help='The number of production steps.')
parser.add_argument('--writeFreq', type=int, help='CSV and DCD write frequency.')
parser.add_argument('--baseTemperature', default=300.00, type=float, help='Temperature of the first replica.')
parser.add_argument('--baseTdiff', type=float, default=10.0, help='Temperature difference between the first two replicas.')
parser.add_argument('--nofReplicas', type=int, default=1, help='Number of replicas.')
parser.add_argument('--runType', type=str, help='Run type: DEFAULT, REMC, RENEMC, RENE, REBAS.')
parser.add_argument('--seed', type=int, help='The seed.')
parser.add_argument('--flexFNs', type=str, nargs='+', default=[], help='The flexFNs.')
parser.add_argument('--FixmanTorque', type=str2bool, nargs='+', default=[], 
                    help='Enable Fixman Torque per world')
parser.add_argument('--roll', type=str2bool, nargs='+', default=[],
					help='Roll over.')
args = parser.parse_args()
#endregion

# Create robosample context
run_type = getattr(robosample.RunType, args.runType)
context = robosample.Context(args.name, args.seed, 0, 1, run_type, 1, 0)

# Set parameters
context.setPdbRestartFreq(0) # WRITE_PDBS
context.setPrintFreq(args.writeFreq) # PRINT_FREQ
context.setNonbonded(0, 1.2)
context.setGBSA(1)
context.setVerbose(True)

# Load system
context.loadAmberSystem(args.top, args.rst7)

# Prepare flexor generator
mdtrajObj = md.load(args.rst7, top=args.top)
flexorObj = flexor.Flexor(mdtrajObj)

nofWorlds = 1 + len(args.flexFNs)
worldIndexes = range(nofWorlds)

# -----------------------------------------------------------
# ------------------ TIMESTEPS and MDSTEPS ------------------
# -----------------------------------------------------------
# region sensitive parameters

nof_replicas = args.nofReplicas

timesteps = [[0.0007 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
mdsteps = [[10 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
samples_per_round = [[1 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
distort_options = [[0 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
distort_args = [["0" for _ in range(nofWorlds)] for _ in range(nof_replicas)]
flow = [[0 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
work = [[0 for _ in range(nofWorlds)] for _ in range(nof_replicas)]
accept_reject_modes = [[robosample.AcceptRejectMode.MetropolisHastings for _ in range(nofWorlds)] for _ in range(nof_replicas)]

# World 0 # Cartesian
if nofWorlds > 0:
	for replIx in range(nof_replicas):
		timesteps[replIx][0] = 0.0007
		mdsteps[replIx][0] = 1000
		samples_per_round[replIx][0] = 1

# World 1 # Roll
if nofWorlds > 1:
	for replIx in range(nof_replicas):
		timesteps[replIx][1] = 0.05
		mdsteps[replIx][1] = 5
		samples_per_round[replIx][1] = 1

# World 2 # BAT stats
if nofWorlds > 2:
	for replIx in range(nof_replicas):
		timesteps[replIx][2] = 0.0007
		mdsteps[replIx][2] = 0
		samples_per_round[replIx][2] = 1

# World 3 # BAT REBAS
if nofWorlds > 3:
	for replIx in range(nof_replicas):
		timesteps[replIx][3] = 0.0007
		mdsteps[replIx][3] = 0
		samples_per_round[replIx][3] = 1
		accept_reject_modes[replIx][3] = robosample.AcceptRejectMode.AlwaysAccept
		distort_options[replIx][3] = -6
		distort_args[replIx][3] = "deterministic"

#endregion
# -----------------------------------------------------------

#region EXPERIMENT SETUP
print("=== EXPERIMENT SETUP ===")
print("Number of worlds:", nofWorlds)
print("Number of replicas:", nof_replicas)
print("Timesteps:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", timesteps[replIx])
print("MD Steps:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", mdsteps[replIx])
print("Samples per round:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", samples_per_round[replIx])
print("Distort options:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", distort_options[replIx])
print("Distort args:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", distort_args[replIx])
print("Flow:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", flow[replIx])
print("Work:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", work[replIx])
print("Accept reject modes:")
for replIx in range(nof_replicas):
	print(" Replica", replIx, ":", accept_reject_modes[replIx])
print("=========================", flush=True)
#endregion EXPERIMENT SETUP

FIXMAN_TORQ = args.FixmanTorque
ROLL = args.roll
addedSoFar = 0

# Add default Openmm Cartesian world
flexes_Cart = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
context.addWorld(FIXMAN_TORQ[addedSoFar], samples_per_round[0][addedSoFar], robosample.RootMobility.CARTESIAN, flexes_Cart, True, False, 0)
addedSoFar += 1

# Add worlds
for flexFNIx, flexFN in enumerate(args.flexFNs):
	flexibilities = getFlexibilitiesFromFile(flexFN)
	#printFlexibilities(flexibilities)

	if FIXMAN_TORQ[addedSoFar] == False:
		print(" === NO FIXMAN TOQRUE FOR WORLD", addedSoFar," ===", file=sys.stderr)

	context.addWorld(FIXMAN_TORQ[addedSoFar], samples_per_round[0][addedSoFar], robosample.RootMobility.WELD, flexibilities, True, False, 0)
	addedSoFar += 1

	if ROLL[addedSoFar - 1] == True:
		print("Setting Roll for world", addedSoFar - 1, "to", ROLL[addedSoFar - 1], file=sys.stderr)
		context.getWorld(addedSoFar - 1).setRollFlexibilities(flexibilities)
	else:
		print(" === NO ROLL FOR WORLD", addedSoFar - 1," ===", file=sys.stderr)

# Add samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN
context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)

for worldIx in range(1, nofWorlds):
	context.getWorld(worldIx).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

# Replica exchange
temperatures = np.zeros(nof_replicas, dtype=np.float64)
Tratio = (args.baseTemperature + args.baseTdiff) / args.baseTemperature
for replIx in range(nof_replicas):
    temperatures[replIx] = args.baseTemperature * (Tratio**replIx)

# Add replicas
integrators = [robosample.IntegratorType.OMMVV] + ((nofWorlds - 1) * [robosample.IntegratorType.VERLET])

# Add replicas with coordinates
context.addReplicasAndLoadCoordinates(args.top, args.rstDir, nof_replicas)

# Add thermodynamic states
for replIx in range(nof_replicas):
    context.addThermodynamicState(replIx,
		temperatures[replIx],
		accept_reject_modes[replIx],
		distort_options[replIx],
		distort_args[replIx],
		flow[replIx],
		work[replIx],
		integrators,
		worldIndexes,
		timesteps[replIx],
		mdsteps[replIx])

# Initialize
context.Initialize()

# Run
context.RunREX(args.equilSteps, args.prodSteps)
