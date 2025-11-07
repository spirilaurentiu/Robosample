# import flexor
import mdtraj as md
import argparse
import robosample
# import robobatstat
import numpy as np
# import os
import argparse
import json

def parse_nested_list(arg):
    """
    Custom type function to parse a JSON string into a Python nested list.
    """
    try:
        data = json.loads(arg)
        if isinstance(data, list) and all(isinstance(sublist, list) for sublist in data):
            # Optional: Add more specific validation if needed
            # e.g., check if all elements in sublists are integers
            for sublist in data:
                for item in sublist:
                    if not isinstance(item, (int, float)): # Or just int if strictly integers
                        raise argparse.ArgumentTypeError(f"List elements must be numbers: {item}")
            return data
        else:
            raise argparse.ArgumentTypeError("Argument must be a nested list (e.g., [[1, 2], [3, 4]])")
    except json.JSONDecodeError:
        raise argparse.ArgumentTypeError("Argument is not a valid JSON string.")
    except Exception as e:
        raise argparse.ArgumentTypeError(f"Error parsing list: {e}")


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
parser.add_argument('type', type=str, help='type of simulation')
parser.add_argument('blocks', type=parse_nested_list, help='A nested list in JSON format (e.g., "[[1, 2, 3], [4, 5, 6]]")')

# Parse the arguments
args = parser.parse_args()

BURN_IN_PERCENTAGE = 0.1
INITIAL_TEMPERATURE = 300.0
NOF_REPLICAS = 1

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier, ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_TD = 0.005 # Torsional dymaics time step is 5 fs
MDSTEPS_TD = 200 # Torsional dynamics block trajectory length 1 ps

TIMESTEP_CARTESIAN = 0.0007 # Cartesian time step is 0.7 fs since we don't use contraints (e.g. SHAKE)
MDSTEPS_CARTESIAN = 357 # Cartesian block trajectory length 250 fs

# create robosample context
context = robosample.Context(name=args.name, seed=args.seed, prmtop=args.prmtop, inpcrd=args.inpcrd, write_freq=args.write_freq)

if args.type == 'tdnr':
    # # 2: Torsional dynamics for phi and psi dihedrals
	# stats = robobatstat.BATCorrelations(args.prmtop, args.inpcrd, args.seed, include_omega=False, include_chi1=True, include_chi2=True, include_chi3=True, include_chi4=True, include_chi5=True)
	# atom_indices = stats.get_dihedral_atom_indices()
	
    # # Add Cartesian flexors (OpenMM)
	# flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
	# context.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

    # # Add torsional flexors for phi and psi dihedrals
	# flex = flexorObj.create_from_list(atom_indices, robosample.BondMobility.Torsion)
	# context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)
     
	context.addCartesianWorld().addSampler(timeStep=TIMESTEP_CARTESIAN, mdSteps=MDSTEPS_CARTESIAN, boostMDSteps=MDSTEPS_CARTESIAN)
    
	bonds = context.getNonRundantBonds() # [[atom1, atom2, joint_type], [atom3, atom4, joint_type], ...]
	context.addTorsionalWorld(bonds).addSampler(timeStep=TIMESTEP_TD, mdSteps=MDSTEPS_TD, boostMDSteps=MDSTEPS_TD)
     
	temperatures = []
	for i in range(NOF_REPLICAS):
		temperatures.append(INITIAL_TEMPERATURE + (i * 10))

	context.initialize(temperatures)
     
	exit()

	context.RunREX(args.equil_steps, args.prod_steps)
	
elif args.type == 'tdc':
    pass

	# # Do torsional dynamics - correlated dihedrals

	# # # Cluster the previous simulations
	# stats = batstat.BATCorrelations(args.prmtop, args.inpcrd)
	# stats = robobatstat.BATCorrelations(args.prmtop, args.inpcrd, args.seed, include_omega=False, include_chi1=True, include_chi2=True, include_chi3=True, include_chi4=True, include_chi5=True)
	# blocks_as_bond_list = []
	# for block in args.blocks:
	# 	l = []
	# 	for b in block:
	# 		aix1 = stats.atom_indices[b][1]
	# 		aix2 = stats.atom_indices[b][2]
	# 		l.append([aix1, aix2])
	# 	blocks_as_bond_list.append(l)
	# blocks = blocks_as_bond_list

	# # Add Cartesian flexors (OpenMM)
	# flex = flexorObj.create(range="all", subset=["all"], jointType="Cartesian")
	# context.addWorld(False, 1, robosample.RootMobility.WELD, flex, True, False, 0)

	# # Create the flexors from the blocks
	# for block in blocks:
	# 	flex = flexorObj.create_from_list(block, robosample.BondMobility.Torsion)
	# 	context.addWorld(True, 1, robosample.RootMobility.WELD, flex, True, False, 0)

	# # Add samplers
	# context.getWorld(0).addSampler(robosample.SamplerName.HMC, robosample.IntegratorType.OMMVV, robosample.ThermostatName.ANDERSEN, False)
	# for i in range(1, len(blocks) + 1):
	# 	context.getWorld(i).addSampler(robosample.SamplerName.HMC, robosample.IntegratorType.VERLET, robosample.ThermostatName.ANDERSEN, True)

	# temperatures = []
	# boost_temperatures = []
	# for i in range(NOF_REPLICAS):
	# 	temperatures.append(INITIAL_TEMPERATURE + (i * 10))
	# 	boost_temperatures.append(INITIAL_TEMPERATURE + (i * 10))  # used for openmm velocities

	# accept_reject_modes = [robosample.AcceptRejectMode.MetropolisHastings] * (len(blocks) + 1)
	# timesteps = [TIMESTEP_CARTESIAN] + [TIMESTEP_TD] * len(blocks)
	# worldIndexes = range(len(blocks) + 1)
	# mdsteps = [MDSTEPS_CARTESIAN] + [MDSTEPS_TD] * len(blocks)
	# boost_md_steps = mdsteps
	# integrators = [robosample.IntegratorType.OMMVV] + [robosample.IntegratorType.VERLET] * len(blocks)

	# distort_options = [0] * (len(blocks) + 1)
	# distort_args = ["0"] * (len(blocks) + 1)
	# flow = [0] * (len(blocks) + 1)
	# work = [0] * (len(blocks) + 1)

	# for i in range(NOF_REPLICAS):
	# 	context.addReplica(i)
	# 	context.addThermodynamicState(i, temperatures[i], accept_reject_modes, distort_options, distort_args, flow, work, integrators, worldIndexes, timesteps, mdsteps)

	# context.Initialize()
	# context.RunREX(args.equil_steps, args.prod_steps)


