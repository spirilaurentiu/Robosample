import sys
sys.path.append("/home/laurentiu/git6/Robosample/bin")
import flexor
import mdtraj as md
import argparse
import robosample
#import batstat
import numpy as np

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

# Openmm Cartesian
flexibilities = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
context.addWorld(False, 1, robosample.RootMobility.WELD, flexibilities, True, False, 0)

# Ramachandran pins
flexibilities = flexorObj.create(range="all", distanceCutoff=0, subset=["rama"], jointType="Pin", sasa_value=-1.0)
context.addWorld(False, 1, robosample.RootMobility.WELD, flexibilities, True, False, 0)

# Sidechains pins
flexibilities = flexorObj.create(range="all", distanceCutoff=0, subset=["all"], jointType="Pin", sasa_value=-1.0)
context.addWorld(False, 1, robosample.RootMobility.WELD, flexibilities, True, False, 0)

# Samplers
sampler = robosample.SamplerName.HMC # rename to type
thermostat = robosample.ThermostatName.ANDERSEN

context.getWorld(0).addSampler(sampler, robosample.IntegratorType.OMMVV, thermostat, False)
context.getWorld(1).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)
context.getWorld(2).addSampler(sampler, robosample.IntegratorType.VERLET, thermostat, True)

# Replica exchange
nof_replicas = 2
temperatures = np.zeros(nof_replicas, dtype=np.float64)
boost_temperatures = np.zeros(nof_replicas, dtype=np.float64)
for replIx in range(nof_replicas):
    temperatures[replIx] = args.baseTemperature + (replIx * 10)
    boost_temperatures[replIx] = args.baseTemperature + (replIx * 10)  # used for openmm velocities

accept_reject_modes = [ robosample.AcceptRejectMode.MetropolisHastings,
						robosample.AcceptRejectMode.MetropolisHastings,
						robosample.AcceptRejectMode.MetropolisHastings]

timesteps = [0.0007, 0.002, 0.0075]
worldIndexes = [0, 1, 2]
mdsteps = [10, 10, 10]
boost_md_steps = mdsteps
integrators = [robosample.IntegratorType.OMMVV, robosample.IntegratorType.VERLET, robosample.IntegratorType.VERLET]

distort_options = [0, 0, 0]
distort_args = ["0", "0" , "0"]
flow = [0, 0, 0]
work = [0, 0, 0]

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
