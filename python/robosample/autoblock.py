import argparse

import numpy as np

import robosample

# python3 python/robosample/autoblock.py 1apq examples/1APQ.prmtop examples/1APQ.rst7 6000 0 20000 1

# Create the parser
parser = argparse.ArgumentParser(description="Process PDB code and seed.")

# Add the arguments
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Relative path to the .prmtop file.")
parser.add_argument("inpcrd", type=str, help="Relative path to the .inpcrd file.")
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="The number of equilibration steps.")
parser.add_argument("prod_steps", type=int, help="The number of production steps.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")

# Parse the arguments
args = parser.parse_args()

# Mean first passage time to cross an energy barrier
# 3 kcal/mol - sub-picosecond to picosecond transitions (modest barrier, ~5KbT)
# 6 kcal/mol - tens to hundreds of picoseconds (moderate barrier, ~10KbT)
# 10 kcal/mol - nanoseconds or longer (high barrier , ~16KbT)
# 2 ps of MD is enough to explore shallow wells, but not to cross deep barriers without enhanced sampling (e.g., HMC, replica exchange)
TIMESTEP_CARTESIAN = 0.001
MDSTEPS_CARTESIAN = 2000

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# # Add cartesian world (will integrate with OpenMM)
# context.addCartesianWorld().add_sampler(
#     timeStep=TIMESTEP_CARTESIAN,
#     mdSteps=MDSTEPS_CARTESIAN,
#     boostMDSteps=MDSTEPS_CARTESIAN,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
# )

# # Add replicas (geometric temperature ladder)
# context.initialize([300])
# context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True, False)

TIMESTEP_TD = 0.001
MDSTEPS_TD = 64
dcd_file = "1apq_6000.repl0.dcd"

for cycle in range(10):
    strong_blocks, weak_blocks, rogue_blocks = (
        context.build_gibbs_blocks_from_trajectory(dcd_file)
    )
    all_blocks = strong_blocks + weak_blocks
    all_blocks.append(np.array(rogue_blocks))

    context = robosample.Context(
        name=args.name + f"_cycle{cycle}",
        seed=args.seed,
        prmtop=args.prmtop,
        inpcrd=args.inpcrd,
        write_freq=args.write_freq,
        testing=False,
    )

    # Add cartesian world (will integrate with OpenMM)
    context.addCartesianWorld().add_sampler(
        timeStep=TIMESTEP_CARTESIAN,
        mdSteps=MDSTEPS_CARTESIAN,
        boostMDSteps=MDSTEPS_CARTESIAN,
        acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    )

    for block in all_blocks:
        bonds = context.standard_dihedral_bonds.loc[block]
        sele = context.build_flexibilities(bonds)
        context.addTorsionalWorld(sele).add_sampler(
            timeStep=TIMESTEP_TD,
            mdSteps=MDSTEPS_TD,
            boostMDSteps=MDSTEPS_TD,
            acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
        )

    context.initialize([300])
    context.run_rex(0, 10, args.write_freq, True, True)
    # context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True, False)

    break
