import argparse

import rdkit
import robosample
import numpy as np
from openmm import app
import openmm as mm
from openmm import unit

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 10 10 1
# python3 roborun.py 2ala.2ala.test ./data-raw/2ala.2ala.prmtop ./data-raw/2ala.2ala.inpcrd 6000 10 10 1

# python3 roborun.py 2ala_test ./data-raw/2ala.prmtop ./data-raw/2ala.inpcrd 6000 0 1000 1
# vmd -f data-raw/2ala.prmtop 2ala_test_6000.repl0.dcd


# python3 roborun.py 2KTA_TEST ./data-raw/2KTA.prmtop ./data-raw/2KTA_min.inpcrd 6000 10 10 1
# python3 roborun.py 2KQ8_TEST ./data-raw/2KQ8.prmtop ./data-raw/2KQ8_min.inpcrd 6000 10 10 1
# python3 roborun.py 2HHI_TEST ./data-raw/2HHI.prmtop ./data-raw/2HHI_min.inpcrd 6000 10 10 1
# python3 roborun.py 2KNQ_TEST ./data-raw/2KNQ.prmtop ./data-raw/2KNQ_min.inpcrd 6000 10 10 1
# python3 roborun.py 2EZA_TEST ./data-raw/2EZA.prmtop ./data-raw/2EZA_min.inpcrd 6000 10 10 1
# python3 roborun.py 2LGJ_TEST ./data-raw/2LGJ.prmtop ./data-raw/2LGJ_min.inpcrd 6000 10 10 1
# python3 roborun.py 2KCU_TEST ./data-raw/2KCU.prmtop ./data-raw/2KCU_min.inpcrd 6000 10 10 1
# python3 roborun.py 2LT2_TEST ./data-raw/2LT2.prmtop ./data-raw/2LT2_min.inpcrd 6000 10 10 1
# python3 roborun.py 2KPU_TEST ./data-raw/2KPU.prmtop ./data-raw/2KPU_min.inpcrd 6000 10 10 1
# python3 roborun.py 2JOZ_TEST ./data-raw/2JOZ.prmtop ./data-raw/2JOZ_min.inpcrd 6000 10 10 1
# python3 roborun.py 1L1I_TEST ./data-raw/1L1I.prmtop ./data-raw/1L1I_min.inpcrd 6000 10 10 1
# python3 roborun.py 1BLA_TEST ./data-raw/1BLA.prmtop ./data-raw/1BLA_min.inpcrd 6000 10 10 1
# python3 roborun.py 2L3B_TEST ./data-raw/2L3B.prmtop ./data-raw/2L3B_min.inpcrd 6000 10 10 1
# python3 roborun.py 2JR6_TEST ./data-raw/2JR6.prmtop ./data-raw/2JR6_min.inpcrd 6000 10 10 1
# python3 roborun.py 1APQ_TEST ./data-raw/1APQ.prmtop ./data-raw/1APQ_min.inpcrd 6000 10 10 1
# python3 roborun.py 1A5E_TEST ./data-raw/1A5E.prmtop ./data-raw/1A5E_min.inpcrd 6000 10 10 1

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

# Parse the arguments
args = parser.parse_args()




import MDAnalysis as mda
import protein


def build_standardized_dihedral_dict():
	dihedrals = {
		# Protein backbone
		("C", "N", "CA", "C") : 'omega',
		("N", "CA", "C", "N") : 'phi',
		("CA", "C", "N", "CA") : 'psi',
	}

	for residue, dihedral_list in protein.PROTEIN_DIHEDRAL_SELECTION.items():
		for dihedral_name, atom_types in dihedral_list.items():
			if 'ring_closing' in dihedral_name:
				continue
			dihedrals[atom_types] = dihedral_name

	return dihedrals

standardized_dihedrals = build_standardized_dihedral_dict()

def get_standardized_dihedral_type(parent_atom: mda.core.groups.Atom, child_atom: mda.core.groups.Atom) -> str | None:
	# Build a list of dihedral candidates by looking at bonded atoms to parent and child
	dihedral_type = None
	for grandparent_atom in parent_atom.bonded_atoms:
		for nephew_atom in child_atom.bonded_atoms:
			if grandparent_atom == child_atom:
				continue
			if nephew_atom == parent_atom:
				continue

			# Try to match into our database of standardized dihedrals
			dihedral_atom_types = (grandparent_atom.name, parent_atom.name, child_atom.name, nephew_atom.name)
			dihedral_type = standardized_dihedrals.get(dihedral_atom_types, None)
			if dihedral_type is None:
				dihedral_type = standardized_dihedrals.get(tuple(reversed(dihedral_atom_types)), None)
			if dihedral_type is not None:
				break

		if dihedral_type is not None:
			break

	return dihedral_type

AROMATIC_TYPES = {
    # GAFF aromatic carbons / heteroatoms
    "ca", "cp", "cq", "cc", "cd", "ce", "cf",
    "na", "nb", "nc", "nd", "ne", "nf",
    # ff19SB aromatic
    "ca", "cb", "cc", "cd", "ce", "cf",
}

# sp / sp2 carbons and heteroatoms (very conservative)
MULTIPLE_BOND_TYPES = {
    # carbonyls, sp2/sp
    "c", "c1", "c2", "ce", "cf", "cg",
    "o", "o2", "os",
    "n", "n2", "n1",
}

# canonical AMBER amide pattern
AMIDE_C = {"c"}      # carbonyl carbon
AMIDE_N = {"n", "nh"}  # amide nitrogens

def is_aromatic(a1 : str, a2: str) -> bool:
    return (a1 in AROMATIC_TYPES) and (a2 in AROMATIC_TYPES)

def is_amide(a1: str, a2: str) -> bool:
    return ((a1 in AMIDE_C and a2 in AMIDE_N) or
            (a2 in AMIDE_C and a1 in AMIDE_N))

def is_multiple_like(a1: str, a2: str) -> bool:
    return (a1 in MULTIPLE_BOND_TYPES and
            a2 in MULTIPLE_BOND_TYPES)

def is_nonrotatable_bond(parent_atom: mda.core.groups.Atom, child_atom: mda.core.groups.Atom) -> bool:
	a1_type = parent_atom.type.lower()
	a2_type = child_atom.type.lower()

	if is_aromatic(a1_type, a2_type):
		return True
	if is_amide(a1_type, a2_type):
		return True
	if is_multiple_like(a1_type, a2_type):
		return True
	return False

import molecule
universe = mda.Universe(args.prmtop, args.inpcrd)

for mol in universe.atoms.fragments:
	# z_matrix = molecule.build_z_matrix(mol)
	# for row in z_matrix:
	# 	# For row (i, j, k, l), parent atom is j and child atom is i
	# 	parent_atom = row.j
	# 	child_atom = row.i

	for bond in mol.bonds:
		parent_atom = bond.atoms[0]
		child_atom = bond.atoms[1]

		# First row does not contain a bond
		if parent_atom is None:
			continue

		# Bonds in which one atom is terminal are non-rotatable
		if len(parent_atom.bonds) == 1 or len(child_atom.bonds) == 1:
			flexible = False
			dihedral_type = "terminal_bond"
			continue
		else:
			flexible = is_nonrotatable_bond(parent_atom, child_atom)
			dihedral_type = get_standardized_dihedral_type(parent_atom, child_atom)
			if dihedral_type is None:
				dihedral_type = "non_standard"

		parent_atom_name = parent_atom.resname + str(parent_atom.resid) + '_' + parent_atom.name + '_' + str(parent_atom.id)
		child_atom_name = child_atom.resname + str(child_atom.resid) + '_' + child_atom.name + '_' + str(child_atom.id)
		
		if not flexible:
			print(f"Non-rotatable bond between parent {parent_atom_name} and child {child_atom_name} of type {dihedral_type}")
		else:
			print(f"Rotatable bond between parent {parent_atom_name} and child {child_atom_name} of type {dihedral_type}")

	break
	print('---')

		

exit()







def getomm_native():
	# Create the system
	prmtop = app.AmberPrmtopFile(args.prmtop)
	inpcrd = app.AmberInpcrdFile(args.inpcrd)

	system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.OBC2)

	# Add Andersen Thermostat
	thermostat = mm.AndersenThermostat(300 * unit.kelvin, 1.0 / unit.picosecond)
	system.addForce(thermostat)

	# Assign each force to a separate group (0, 1, 2, etc.)
	# This must be done before creating the Simulation
	dict_forces = {}
	for i, force in enumerate(system.getForces()):
		force.setForceGroup(i)
		dict_forces[force.getName()] = i

	integrator = mm.VerletIntegrator(0.001 * unit.picoseconds)
	platform = mm.Platform.getPlatformByName('CUDA')
	simulation = app.Simulation(prmtop.topology, system, integrator, platform)
	simulation.context.setPositions(inpcrd.positions)

	# Get total energy first
	state = simulation.context.getState(getEnergy=True)
	total_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)

	# Assign values to each force group
	force_groups = {
		'totalEnergy': total_energy
	}

	# Get individual components by force group
	for name, group_id in dict_forces.items():
		state = simulation.context.getState(getEnergy=True, groups={group_id})
		component_energy = state.getPotentialEnergy().value_in_unit(unit.kilojoules_per_mole)
		force_groups[name] = component_energy

	return force_groups


TOL = 1e-6

# Temperature replica exchange parameters
T0 = 300.0
T_MAX = 600.0
NOF_REPLICAS = 2
R = (T_MAX / T0) ** (1.0 / (NOF_REPLICAS - 1))

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier , ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_TD = 0.005 # Torsional dymaics time step is 5 fs
MDSTEPS_TD = 200 # Torsional dynamics block trajectory length 1 ps

TIMESTEP_CARTESIAN = 0.0007 # Cartesian time step is 0.7 fs since we don't use contraints (e.g. SHAKE)
MDSTEPS_CARTESIAN = 357 # Cartesian block trajectory length 250 fs

# create robosample context
context = robosample.Context(name=args.name, seed=args.seed, prmtop=args.prmtop, inpcrd=args.inpcrd, write_freq=args.write_freq, testing=True)

# Add cartesian world (will integrate with OpenMM)
context.addCartesianWorld().addSampler(timeStep=TIMESTEP_CARTESIAN, mdSteps=MDSTEPS_CARTESIAN, boostMDSteps=MDSTEPS_CARTESIAN, acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept)

# Add torsional world with non-redundant dihedrals
bonds = context.getNonRundantBonds() # [[atom1, atom2, joint_type], [atom3, atom4, joint_type], ...]
context.addTorsionalWorld(bonds).addSampler(timeStep=TIMESTEP_TD, mdSteps=MDSTEPS_TD, boostMDSteps=MDSTEPS_TD)

# Add replicas (geometric temperature ladder)
temperatures = []
for i in range(NOF_REPLICAS):
    temperatures.append(T0 * (R ** i))

context.initialize(temperatures)
context.RunREX(args.equil_steps, args.prod_steps)

# openmm vs dumm

# Test initial OpenMM energy components
omm_native = getomm_native()
omm_robo = context.getInitialOpenMMEnergyComponents()

OPENMM_TOL = 10
assert abs(omm_native['totalEnergy'] - omm_robo.totalEnergy) < OPENMM_TOL,\
	f"Total energy mismatch: Native OpenMM = {omm_native['totalEnergy']}, RoboSample OpenMM = {omm_robo.totalEnergy}"
assert abs(omm_native['HarmonicBondForce'] - omm_robo.harmonicBondForce) < OPENMM_TOL,\
	f"HarmonicBondForce mismatch: Native OpenMM = {omm_native['HarmonicBondForce']}, RoboSample OpenMM = {omm_robo.harmonicBondForce}"
assert abs(omm_native['HarmonicAngleForce'] - omm_robo.harmonicAngleForce) < OPENMM_TOL,\
	f"HarmonicAngleForce mismatch: Native OpenMM = {omm_native['HarmonicAngleForce']}, RoboSample OpenMM = {omm_robo.harmonicAngleForce}"
assert abs(omm_native['PeriodicTorsionForce'] - omm_robo.periodicTorsionForce) < OPENMM_TOL,\
	f"PeriodicTorsionForce mismatch: Native OpenMM = {omm_native['PeriodicTorsionForce']}, RoboSample OpenMM = {omm_robo.periodicTorsionForce}"
assert abs(omm_native['NonbondedForce'] - omm_robo.nonbondedForce) < OPENMM_TOL,\
	f"NonbondedForce mismatch: Native OpenMM = {omm_native['NonbondedForce']}, RoboSample OpenMM = {omm_robo.nonbondedForce}"
assert abs(omm_native['AndersenThermostat'] - omm_robo.andersenThermostat) < OPENMM_TOL,\
	f"AndersenThermostat mismatch: Native OpenMM = {omm_native['AndersenThermostat']}, RoboSample OpenMM = {omm_robo.andersenThermostat}"
assert abs(omm_native['GBSAOBCForce'] - omm_robo.gbsaObcForce) < OPENMM_TOL,\
	f"GBSAOBCForce mismatch: Native OpenMM = {omm_native['GBSAOBCForce']}, RoboSample OpenMM = {omm_robo.gbsaObcForce}"

# Test the simulation
for i, world in enumerate(context.getWorlds()):
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
			assert rmsd > TOL, "RMSD for accepted move is too small: {}: the molecule did not move at all.".format(rmsd)
		else:
			assert rmsd <= TOL, "RMSD for rejected move is too large: {}: the molecule moved too much.".format(rmsd)

	# Bond lengths inside rigid bodies should not change values
	# The world is not perfect however and we need to allow a high margin (up to 1 nm RMSD for the entire molecule)
	rigidBodyBondRMSDInNm = world.getRigidBodyBondRMSDInNm()
	np.testing.assert_allclose(rigidBodyBondRMSDInNm, 0, atol=TOL)

	rigidBodyAngleDriftInRad = world.getRigidBodyAngleDriftInRad()
	np.testing.assert_allclose(rigidBodyAngleDriftInRad, 0, atol=TOL)

	rigidBodyProperTorsionDriftInRad = world.getRigidBodyProperTorsionDriftInRad()
	np.testing.assert_allclose(rigidBodyProperTorsionDriftInRad, 0, atol=TOL)

	rigidBodyImproperTorsionDriftInRad = world.getRigidBodyImproperTorsionDriftInRad()
	np.testing.assert_allclose(rigidBodyImproperTorsionDriftInRad, 0, atol=TOL)

# ALA2_OXT_23 ALA2_C_21 ALA2_CA_15 ALA2_HA_16  : Dihedral index in universe: 43
# ALA2_OXT_23 ALA2_C_21 ALA2_CA_15 ALA2_CB_17  : Dihedral index in universe: 45
# ALA2_HA_16 ALA2_CA_15 ALA2_CB_17 ALA2_HB1_18  : Dihedral index in universe: 39
# ALA2_HA_16 ALA2_CA_15 ALA2_CB_17 ALA2_HB2_19  : Dihedral index in universe: 40
# ALA2_HA_16 ALA2_CA_15 ALA2_CB_17 ALA2_HB3_20  : Dihedral index in universe: 41
# ALA2_HA_16 ALA2_CA_15 ALA2_C_21 ALA2_O_22  : Dihedral index in universe: 42
# ALA2_HB1_18 ALA2_CB_17 ALA2_CA_15 ALA2_N_13  : Dihedral index in universe: 31
# ALA2_HA_16 ALA2_CA_15 ALA2_N_13 ALA2_H_14  : Dihedral index in universe: 36
# ALA2_HA_16 ALA2_CA_15 ALA2_N_13 ALA1_C_11  : Dihedral index in universe: 26
# ALA2_H_14 ALA2_N_13 ALA1_C_11 ALA1_CA_5  : Dihedral index in universe: 14
# ALA2_N_13 ALA1_C_11 ALA1_CA_5 ALA1_HA_6  : Dihedral index in universe: 20
# ALA2_N_13 ALA1_C_11 ALA1_CA_5 ALA1_CB_7  : Dihedral index in universe: 22
# ALA1_HA_6 ALA1_CA_5 ALA1_CB_7 ALA1_HB1_8  : Dihedral index in universe: 16
# ALA1_HA_6 ALA1_CA_5 ALA1_CB_7 ALA1_HB2_9  : Dihedral index in universe: 17
# ALA1_HA_6 ALA1_CA_5 ALA1_CB_7 ALA1_HB3_10  : Dihedral index in universe: 18
# ALA1_HA_6 ALA1_CA_5 ALA1_C_11 ALA1_O_12  : Dihedral index in universe: 19
# ALA1_HB1_8 ALA1_CB_7 ALA1_CA_5 ALA1_N_1  : Dihedral index in universe: 0
# ALA1_HA_6 ALA1_CA_5 ALA1_N_1 ALA1_H1_2  : Dihedral index in universe: 5
# ALA1_HA_6 ALA1_CA_5 ALA1_N_1 ALA1_H2_3  : Dihedral index in universe: 8
# ALA1_HA_6 ALA1_CA_5 ALA1_N_1 ALA1_H3_4  : Dihedral index in universe: 11


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