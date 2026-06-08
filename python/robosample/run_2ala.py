import argparse

import robosample

# python3 python/robosample/roborun.py ala-dipeptide examples/ala-dipeptide.prmtop examples/ala-dipeptide.rst7 6000 0 10 1

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

# create robosample context
context = robosample.Context(
    name=args.name,
    seed=args.seed,
    prmtop=args.prmtop,
    inpcrd=args.inpcrd,
    write_freq=args.write_freq,
    testing=False,
)

# context.add_cartesian_world(temperature=300).add_sampler(
#     timeStep=0.001,
#     mdSteps=500,
#     boostMDSteps=500,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=False,
# )

mask_phi = context.standard_dihedral_bonds["dihedral_type"] == "phi"
mask_psi = context.standard_dihedral_bonds["dihedral_type"] == "psi"
bonds = context.standard_dihedral_bonds[mask_phi | mask_psi]
# bonds = bonds.iloc[[0]]
sele = context.build_flexibilities(bonds, robosample.rb.BondMobility.Torsion, False)
context.add_robotic_world(sele, temperature=3000).add_sampler(
    timeStep=0.025,
    mdSteps=20,  # ignored if using NUTS
    boostMDSteps=20,  # ignored if using NUTS
    acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    use_nuts=False,
)

# for row in context.z_matrix:
#     print(row.global_indices)

context.initialize([300])

# # Run the simulation
# start_time = time.perf_counter()

# context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)

# end_time = time.perf_counter()
# duration = end_time - start_time

# print(f"run_rex() took {duration:.4f} seconds")
