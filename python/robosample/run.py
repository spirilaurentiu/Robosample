import argparse

import openmm_validation

import robosample

# Explicit-solvent (PME / PBC) smoke test.
#
#   python3 python/robosample/run.py myrun system.prmtop system.rst7 1 5 50 1
#
# The prmtop/rst7 MUST carry a periodic box (a solvated system from tleap:
# solvateBox / solvateOct). The rst7 is read for both coordinates AND the box.

parser = argparse.ArgumentParser(description="Explicit-solvent (PME) Robosample run.")
parser.add_argument("name", type=str, help="Name of the simulation.")
parser.add_argument("prmtop", type=str, help="Path to the .prmtop file.")
parser.add_argument(
    "inpcrd", type=str, help="Path to the .inpcrd/.rst7 file (with box)."
)
parser.add_argument("seed", type=int, help="The seed.")
parser.add_argument("equil_steps", type=int, help="Number of equilibration rounds.")
parser.add_argument("prod_steps", type=int, help="Number of production rounds.")
parser.add_argument("write_freq", type=int, help="CSV and DCD write frequency.")
parser.add_argument(
    "validate", type=bool, help="Whether to run the OpenMM validation (0 or 1)."
)
args = parser.parse_args()

# ---- Build the context in EXPLICIT-SOLVENT mode -----------------------------
# explicit_solvent=True:
#   * PME electrostatics, GBSA off,
#   * box read from the rst7 and stored as reduced lattice vectors,
#   * 1.0 nm cutoff (override with nonbonded_cutoff=...),
#   * every molecule (solute + each water) gets a FREE 6-DOF root.
dih_classifier = robosample.AmberDihedralClassifier()
context = robosample.Context(args.name, args.seed, dih_classifier)
context.load_amber(args.prmtop, args.inpcrd, explicit_solvent=True)

# Keep whole molecules across worlds: never let OpenMM wrap the coordinates that
# flow into the robot engine. This is the default; set here explicitly so the
# intent is on the page. (Energies/forces are unaffected -- OpenMM always applies
# the minimum image internally; this only governs returned positions.)
context.set_enforce_periodic_box(False)

# ---- Optional: verify the PME energy matches an OpenMM reference -------------
if args.validate:
    ok, cpp_pe, ref_pe, per_class = openmm_validation.compare_by_force_group(
        context, args.prmtop, args.inpcrd, platform_name="Reference"
    )
    print(f"[validate] C++ PME PE = {cpp_pe:.6f} kJ/mol")
    print(f"[validate] ref PME PE = {ref_pe:.6f} kJ/mol")
    print(f"[validate] per-force-group match: {'OK' if ok else 'MISMATCH'}")
    for name, (c, r) in per_class.items():
        print(f"[validate]   {name:24s} cpp={c:14.4f}  ref={r:14.4f}  d={c - r:+.4e}")

# ---- Worlds (the Gibbs sweep order = the move schedule) ---------------------

# # (1) Cartesian all-atom MD world. This runs on the device through OpenMM with
# #     PME, so PBC is handled entirely by OpenMM. Velocities are reseeded to the
# #     replica temperature each round (Andersen-style), then `mdSteps` of Verlet.
# context.add_cartesian_world().add_sampler(
#     timeStep=0.001,
#     mdSteps=100,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
#     use_nuts=False,
# )

# (2) Optional internal-coordinate world on the solute (molecule 0) phi/psi.
#     The torsional Verlet needs no box -- it moves only internal DOF and the
#     forces arrive minimum-imaged from OpenMM. Every water is a Free-root rigid
#     body here, so it gets a rigid-body HMC move as well.
dihedrals = [
    robosample.DihedralType.PROTEIN_PHI.value,
    robosample.DihedralType.PROTEIN_PSI.value,
]
bonds = context.standard_dihedral_bonds.loc[
    context.standard_dihedral_bonds["dihedral_type"].isin(dihedrals)
]
sele = context.build_flexibilities(bonds, robosample.rb.BondMobility.Torsion, False)
# context.add_robotic_world(sele).add_sampler(
#     timeStep=0.04,
#     mdSteps=25,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
#     use_nuts=False,
#     use_fixman=True,  # required for rigorous Boltzmann sampling of a constrained world
# )

context.add_ncmc_world(
    sele,
    molecule_index=0,
    timestep=0.002,
    ncmc_steps=1000,
    hold_fraction=0.2,
    use_fixman=True,  # rigorous Boltzmann sampling of the constrained world
)


# ---- Run --------------------------------------------------------------------
# initialize() runs an O(N^2) startup clash scan (now minimum-image aware). On a
# large solvent box this takes a moment; if your box is not yet minimized and the
# scan flags clashes, set ROBO_ALLOW_BAD_START=1 in the environment to downgrade
# the hard error to a warning, or minimize first.
context.initialize([300])
context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)
