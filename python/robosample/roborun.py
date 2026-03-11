import argparse

# import rdkit
import robosample
import numpy as np
from openmm import app
import openmm as mm
from openmm import unit

# python3 roborun.py 2but ../examples/2but.prmtop ../examples/2but.rst7 6000 0 10 1

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 1 1 1
# python3 roborun.py 2ala.2ala.test ./data-raw/2ala.2ala.prmtop ./data-raw/2ala.2ala.inpcrd 6000 0 10 1

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 0 1000 1
# vmd -f data-raw/2ala.prmtop 2ala_test_6000.repl0.dcd

# python3 roborun.py gfcd ./data-raw/GfcD121A-1min.prmtop ./data-raw/GfcD121A-1min.inpcrd 6000 0 10 1
# python3 roborun.py gfcd ./data-raw/FFAR1_nanodisc.AMBER.prmtop ./data-raw/FFAR1_nanodisc.AMBER.min.rst7 6000 0 10 1

# python3 roborun.py 2KTA_TEST ./data-raw/2KTA.prmtop ./data-raw/2KTA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KQ8_TEST ./data-raw/2KQ8.prmtop ./data-raw/2KQ8_min.inpcrd 6000 0 10 1
# python3 roborun.py 2HHI_TEST ./data-raw/2HHI.prmtop ./data-raw/2HHI_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KNQ_TEST ./data-raw/2KNQ.prmtop ./data-raw/2KNQ_min.inpcrd 6000 0 10 1
# python3 roborun.py 2EZA_TEST ./data-raw/2EZA.prmtop ./data-raw/2EZA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2LGJ_TEST ./data-raw/2LGJ.prmtop ./data-raw/2LGJ_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KCU_TEST ./data-raw/2KCU.prmtop ./data-raw/2KCU_min.inpcrd 6000 0 10 1
# python3 roborun.py 2LT2_TEST ./data-raw/2LT2.prmtop ./data-raw/2LT2_min.inpcrd 6000 0 10 1
# python3 roborun.py 2KPU_TEST ./data-raw/2KPU.prmtop ./data-raw/2KPU_min.inpcrd 6000 0 10 1
# python3 roborun.py 2JOZ_TEST ./data-raw/2JOZ.prmtop ./data-raw/2JOZ_min.inpcrd 6000 0 10 1
# python3 roborun.py 1L1I_TEST ./data-raw/1L1I.prmtop ./data-raw/1L1I_min.inpcrd 6000 0 10 1
# python3 roborun.py 1BLA_TEST ./data-raw/1BLA.prmtop ./data-raw/1BLA_min.inpcrd 6000 0 10 1
# python3 roborun.py 2L3B_TEST ./data-raw/2L3B.prmtop ./data-raw/2L3B_min.inpcrd 6000 0 10 1
# python3 roborun.py 2JR6_TEST ./data-raw/2JR6.prmtop ./data-raw/2JR6_min.inpcrd 6000 0 10 1
# python3 roborun.py 1APQ_TEST ./data-raw/1APQ.prmtop ./data-raw/1APQ_min.inpcrd 6000 0 10 1
# python3 roborun.py 1A5E_TEST ./data-raw/1A5E.prmtop ./data-raw/1A5E_min.inpcrd 6000 0 10 1

# Create the parser
parser = argparse.ArgumentParser(description='Process PDB code and seed.')

# Add the arguments
parser.add_argument('name', type=str, help='Name of the simulation.')
parser.add_argument('prmtop', type=str, help='Relative path to the .prmtop file.')
parser.add_argument('inpcrd', type=str, help='Relative path to the .inpcrd file.')
parser.add_argument('seed', type=int, help='The seed.')
parser.add_argument('equil_steps', type=int, help='The number of equilibration steps.')
parser.add_argument('prod_steps', type=int, help='The number of production steps.')
parser.add_argument('write_freq', type=int, help='CSV and DCD write frequency.')

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
from MDAnalysis.analysis.bat import BAT
import matplotlib.pyplot as plt

# def calculate_l(dcd):
# 	u = mda.Universe(args.prmtop, dcd)
# 	selection = u.select_atoms('protein')

# 	# Build BAT representation
# 	bat = BAT(selection)
# 	bat.run()

# 	coords = bat.results.bat  # shape: (n_frames, n_internal_coords)

# 	# reference frame
# 	ref = coords[0]

# 	# RMSD in internal coordinate space
# 	diff = coords - ref
# 	rmsd = np.sqrt(np.mean(diff**2, axis=1))

# 	# autocorrelation function
# 	rmsd_centered = rmsd - np.mean(rmsd)
# 	acf = np.correlate(rmsd_centered, rmsd_centered, mode='full')
# 	acf = acf[acf.size // 2:]
# 	acf /= acf[0]

# 	# Integrated autocorrelation time
# 	# Stop at first zero crossing to avoid noise in the tail
# 	cutoff = np.argmax(acf < 0) or len(acf)
# 	tau = 1 + 2 * np.sum(acf[1:cutoff])
# 	L = tau / 2

# 	return L

# def get_native_openmm_energy():
# 	prmtop = app.AmberPrmtopFile(args.prmtop)
# 	inpcrd = app.AmberInpcrdFile(args.inpcrd)

# 	system = prmtop.createSystem(
# 		nonbondedMethod=app.CutoffNonPeriodic,
# 		nonbondedCutoff=1.2,
# 		constraints=None,
# 		implicitSolvent=app.OBC2,
# 		removeCMMotion=False
# 	)

# 	thermostat = mm.AndersenThermostat(300 * unit.kelvin, 1.0 / unit.picosecond)
# 	system.addForce(thermostat)

# 	integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
# 	platform = mm.Platform.getPlatformByName('CUDA')
# 	simulation = app.Simulation(prmtop.topology, system, integrator, platform)

# 	simulation.context.setPositions(inpcrd.positions)
# 	simulation.step(1000)
# 	equil_positions = simulation.context.getState(getPositions=True).getPositions()

# 	for L in [5, 10, 20, 50, 100, 200]:
# 		variances = []
# 		for trial in range(10):  # 50 independent short runs
# 			system = prmtop.createSystem(
# 				nonbondedMethod=app.CutoffNonPeriodic,
# 				nonbondedCutoff=1.2,
# 				constraints=None,
# 				implicitSolvent=app.OBC2,
# 				removeCMMotion=False
# 			)

# 			thermostat = mm.AndersenThermostat(300 * unit.kelvin, 1.0 / unit.picosecond)
# 			system.addForce(thermostat)

# 			integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
# 			platform = mm.Platform.getPlatformByName('CUDA')
# 			simulation = app.Simulation(prmtop.topology, system, integrator, platform)
# 			simulation.context.setPositions(equil_positions)

# 			simulation.reporters.append(app.DCDReporter(f'{args.name}_openmm.dcd', 1))
# 			simulation.step(L)

# 			variances.append(calculate_l(f'{args.name}_openmm.dcd'))
# 		print(f"L={L}: mean_var={np.mean(variances):.4f}")

# 	exit()

# 	# Get total energy first
# 	state = simulation.context.getState(getEnergy=True)
# 	return state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

# Parse the arguments
args = parser.parse_args()

# get_native_openmm_energy()

# Temperature replica exchange parameters
T0 = 300.0
T_MAX = 1000.0
NOF_REPLICAS = 1
R = 1 if NOF_REPLICAS == 1 else (T_MAX / T0) ** (1.0 / (NOF_REPLICAS - 1))

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier , ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_TD = 0.015 # Torsional dymaics time step is 10 fs
MDSTEPS_TD = 100 # Torsional dynamics block trajectory length 1 ps

TIMESTEP_CARTESIAN = 0.001
MDSTEPS_CARTESIAN = 15

# create robosample context
context = robosample.Context(name=args.name, seed=args.seed, prmtop=args.prmtop, inpcrd=args.inpcrd, write_freq=args.write_freq, testing=True)

# Add cartesian world (will integrate with OpenMM)
context.addCartesianWorld().addSampler(timeStep=TIMESTEP_CARTESIAN, mdSteps=MDSTEPS_CARTESIAN, boostMDSteps=MDSTEPS_CARTESIAN, acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept)

# Add torsional world with non-redundant dihedrals
# sele = [
# 	context.selectBonds('resid 0'),
# 	# context.selectBonds('resid 1')
# ]
# sele = context.getDefaultBonds('standard')
# context.addTorsionalWorld(sele).addSampler(timeStep=TIMESTEP_TD, mdSteps=MDSTEPS_TD, boostMDSteps=MDSTEPS_TD)

# for flex in context.getDefaultBonds('macrocycle'):
# 	context.addTorsionalWorld([flex]).addSampler(timeStep=TIMESTEP_TD, mdSteps=MDSTEPS_TD, boostMDSteps=MDSTEPS_TD, acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept)

# Add replicas (geometric temperature ladder)
temperatures = []
for i in range(NOF_REPLICAS):
    temperatures.append(T0 * (R ** i))

context.initialize(temperatures)

# print(context.calculate_openmm_energy(0))
# print(context.calculate_openmm_energy(1))

# # Test OpenMM energies
# robosample_openmm = context.calculate_openmm_energy(1)
# native_openmm = get_native_openmm_energy()
# if robosample_openmm != native_openmm:
# 	message = f"OpenMM energy calculated by Robosample ({robosample_openmm}) does not match native OpenMM energy ({native_openmm})."
# 	raise RuntimeError(message)

# Run the simulation
context.RunREX(args.equil_steps, args.prod_steps)

# Test the simulation
TOL = 1e-4
for i, world in enumerate(context.getWorlds()):
	if i == 0:
		openmmCumulativeCartesianDisplacements = world.getOpenMMCumulativeCartesianDisplacements()
		np.testing.assert_allclose(openmmCumulativeCartesianDisplacements, 0, atol=TOL)

		openmmCumulativeBondDisplacements = world.getOpenMMCumulativeBondDisplacements()
		np.testing.assert_allclose(openmmCumulativeBondDisplacements, 0, atol=TOL)

		openmmCumulativeAngleDisplacements = world.getOpenMMCumulativeAngleDisplacements()
		np.testing.assert_allclose(openmmCumulativeAngleDisplacements, 0, atol=TOL)

		openmmCumulativeTorsionDisplacements = world.getOpenMMCumulativeTorsionDisplacements()
		np.testing.assert_allclose(openmmCumulativeTorsionDisplacements, 0, atol=TOL)
	else:
		matchAtomTargetLocationsResiduals = world.getMatchAtomTargetLocationsResiduals()
		np.testing.assert_allclose(matchAtomTargetLocationsResiduals, 0, atol=TOL)

		cumulativeCartesianDisplacements = world.getCumulativeCartesianDisplacements()
		np.testing.assert_allclose(cumulativeCartesianDisplacements, 0, atol=TOL)

		cumulativeBondDisplacements = world.getCumulativeBondDisplacements()
		np.testing.assert_allclose(cumulativeBondDisplacements, 0, atol=TOL)

		cumulativeAngleDisplacements = world.getCumulativeAngleDisplacements()
		np.testing.assert_allclose(cumulativeAngleDisplacements, 0, atol=TOL)

		cumulativeTorsionDisplacements = world.getCumulativeTorsionDisplacements()
		np.testing.assert_allclose(cumulativeTorsionDisplacements, 0, atol=TOL)

	acceptanceRMSD = world.getAcceptanceRMSD()
	for accepted, rmsd in acceptanceRMSD:
		if accepted:
			if rmsd <= TOL:
				raise RuntimeError(f"RMSD for accepted move is too small ({rmsd}): the molecule did not move at all. A possible cause is that equil_steps={args.equil_steps}")
		else:
			if rmsd > TOL:
				raise RuntimeError(f"RMSD for rejected move is too large ({rmsd}): the molecule moved too much.")

	# Bond lengths inside rigid bodies should not change values
	# The world is not perfect however and we need to allow a high margin (up to 1 nm RMSD for the entire molecule)
	rigidBodyBondRMSDInNm = world.getRigidBodyBondRMSDInNm()
	np.testing.assert_allclose(rigidBodyBondRMSDInNm, 0, atol=TOL)

	rigidBodyAngleDriftInRad = world.getRigidBodyAngleDriftInRad()
	np.testing.assert_allclose(rigidBodyAngleDriftInRad, 0, atol=1)

	rigidBodyProperTorsionDriftInRad = world.getRigidBodyProperTorsionDriftInRad()
	np.testing.assert_allclose(rigidBodyProperTorsionDriftInRad, 0, atol=TOL)

	rigidBodyImproperTorsionDriftInRad = world.getRigidBodyImproperTorsionDriftInRad()
	np.testing.assert_allclose(rigidBodyImproperTorsionDriftInRad, 0, atol=TOL)




"""
source leaprc.protein.ff19SB
mol1 = sequence { ACE ALA ALA NME }
mol2 = sequence { ACE ALA ALA NME }
translate mol2 { 15.0 0.0 0.0 }
system = combine { mol1 mol2 }
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.inpcrd
saveamberparm system 2ala.2ala.prmtop 2ala.2ala.rst7
savepdb system 2ala.2ala.pdb
quit
"""

"""
source leaprc.protein.ff19SB
system = sequence { ACE GLY ALA VAL LEU ILE SER THR CYS MET PHE TYR TRP PRO ASP GLU ASN GLN HIS LYS ARG NME }
saveamberparm system all.prmtop all.inpcrd
saveamberparm system all.prmtop all.rst7
savepdb system all.pdb
quit
"""