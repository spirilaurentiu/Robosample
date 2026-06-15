import argparse

# import openmm_validation
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

# platform = "CUDA"
# ok, cpp_pe, ref_pe, per_class = openmm_validation.compare_by_force_group(
#     context, args.prmtop, args.inpcrd, platform_name=platform
# )

# # All molecule roots are `robosample.rb.RootMobility.Weld`
# context.set_root_mobility(0, robosample.rb.RootMobility.WELD)

# ---- Worlds (the Gibbs sweep order = the move schedule) ----------------------

# (1) Docking world: molecule 0 is the ligand (Free root); everything else is
#     welded/rigid. The binding sphere is sized AUTOMATICALLY, per ligand, as
#     R_receptor + sphere_factor * R_ligand. The ligand is repositioned (uniform
#     position + reorientation) only when its COM leaves the sphere; a proposal
#     whose energy is non-finite or |PE| > clash_threshold is rejected (in every
#     mode), so an overlap with the receptor never passes.
context.add_docking_world(ligand_molecule_indices=[0]).add_sampler(
    # 0.1 ps was ~100x too large for all-atom MD: from any mildly strained pose the
    # Verlet step diverges to non-finite within mdSteps, so the move is rejected and
    # the pose is restored unchanged -- a frozen-PE absorbing state. 0.002 ps is in
    # the stable range the engine recommends. mdSteps raised 10 -> 50 to preserve the
    # ~0.1 ps trajectory length (0.002 * 50 = 0.1 ps) at the smaller step.
    timeStep=0,
    mdSteps=40,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
    sphere_factor=0.5,
    clash_threshold=10.0,  # reject if |peNew| > 10*|pePre| (relative; 1 order of magnitude)
    always_kick=False,
    max_initial_kick_tries=1000,
)


# # (2) Torsional relaxation world (phi/psi flexible) WITH the Fixman correction in
# #     the acceptance Hamiltonian -- required for rigorous Boltzmann sampling of a
# #     constrained (internal-coordinate) world.
# dihedrals = [
#     robosample.DihedralType.PROTEIN_PHI.value,
#     robosample.DihedralType.PROTEIN_PSI.value,
# ]
# bonds = context.standard_dihedral_bonds.loc[
#     context.standard_dihedral_bonds["dihedral_type"].isin(dihedrals)
# ]
# sele = context.build_flexibilities(bonds, robosample.rb.BondMobility.Torsion, False)
# context.add_robotic_world(sele).add_sampler(
#     timeStep=0.005,
#     mdSteps=20,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
#     use_nuts=False,
#     use_fixman=True,
# )

# (3) Cartesian all-atom MD world for mixing (uses the MTS integrator).
context.add_cartesian_world().add_sampler(
    timeStep=0.001,
    mdSteps=1000,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,  # Andersen-thermostatted MD
    use_nuts=False,
)

# [rex] round=12175 replica=0 T=300.0 world=0(docking) kick=N acc=rej PE=189.8145 KE=0.0000 Fix=-902.1686 H=-677.8935
# From -2400 -> wrong conf, but cannot detect it's wrong

context.initialize([300])
context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)
