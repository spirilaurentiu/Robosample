import argparse

import robosample

# python3 python/robosample/run.py example examples/example.prmtop examples/example.rst7 6000 0 1 1

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

# All molecule roots are `robosample.rb.RootMobility.Weld`
# context.set_root_mobility(0, robosample.rb.RootMobility.Free)

# All other molecules are lipids and solvent
for mol_ix in range(3, context.get_num_molecules()):
    context.set_root_mobility(mol_ix, robosample.rb.RootMobility.Free)

# context.add_cartesian_world().add_sampler(
#     timeStep=0.001,
#     mdSteps=50_000,
#     boostMDSteps=50_000,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
#     use_nuts=False,
# )

############### World 2a ###############
resid = [103, 136, 173, 182, 190, 201, 232, 236, 243, 257, 275]
dihs = ["chi1", "chi2", "chi3", "chi4", "chi5"]

bonds = context.standard_dihedral_bonds.loc[
    (context.standard_dihedral_bonds["resid"].isin(resid))
    & (context.standard_dihedral_bonds["dihedral_type"].isin(dihs))
    & (context.standard_dihedral_bonds["molecule_index"] == 0)
]

sele = context.build_flexibilities(bonds, robosample.rb.BondMobility.Torsion)
context.add_robotic_world(sele).add_sampler(
    timeStep=0.005,  # 5 fs
    mdSteps=500,
    boostMDSteps=500,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)

# Add replicas (geometric temperature ladder)
context.initialize([300])

context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, False)
