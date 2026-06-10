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

# Create robosample context
dih_classifier = robosample.AmberDihedralClassifier()
context = robosample.Context(args.name, args.seed, dih_classifier)
context.load_amber(args.prmtop, args.inpcrd)

# All molecule roots are `robosample.rb.RootMobility.Weld`
# context.set_root_mobility(0, robosample.rb.RootMobility.Free)

# Host is welded to ground
context.set_root_mobility(0, robosample.rb.RootMobility.Weld)

# Guest is free to roam
context.set_root_mobility(1, robosample.rb.RootMobility.Free)

context.add_cartesian_world().add_sampler(
    timeStep=0.001,
    mdSteps=50_000,
    boostMDSteps=50_000,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)

sele = context.build_flexibilities(None, robosample.rb.BondMobility.Torsion, False)
context.add_robotic_world(sele).add_sampler(
    timeStep=0.25,
    mdSteps=10,
    boostMDSteps=10,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
)

# Add replicas (geometric temperature ladder)
context.initialize([300])

context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)
