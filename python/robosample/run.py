import argparse

import openmm_validation

import robosample

# Explicit-solvent (PME / PBC) smoke test.
#
#   python3 python/robosample/run.py 2ala tip3p/2ala.prmtop tip3p/2ala.rst7 6000 0 100 1 true
#
#   python3 python/robosample/run.py 2ala.implicit examples/2ala/2ala.implicit.prmtop examples/2ala/2ala.implicit.rst7 6000 0 0 1 true
#   python3 python/robosample/run.py 2ala.tip3p examples/2ala/2ala.tip3p.prmtop examples/2ala/2ala.tip3p.rst7 6000 0  0 1 true
#
# The prmtop/rst7 MUST carry a periodic box (a solvated system from tleap:
# solvateBox / solvateOct). The rst7 is read for both coordinates AND the box.

parser = argparse.ArgumentParser(description="Explicit-solvent (PME) Robosample run.")
parser.add_argument("--name", type=str, help="Name of the simulation.")
parser.add_argument("--prmtop", type=str, help="Path to the .prmtop file.")
parser.add_argument(
    "inpcrd", type=str, help="Path to the .inpcrd/.rst7 file (with box)."
)
parser.add_argument("--seed", type=int, help="The seed.")
parser.add_argument("--equil_steps", type=int, help="Number of equilibration rounds.")
parser.add_argument("--prod_steps", type=int, help="Number of production rounds.")
parser.add_argument("--write_freq", type=int, help="CSV and DCD write frequency.")
parser.add_argument(
    "--validate", type=bool, help="Whether to run the OpenMM validation (0 or 1)."
)
args = parser.parse_args()

# ---- Build the context ------------------------------------------------------
# The solvent model is auto-detected from the box (OpenMM-style): this rst7
# carries a periodic box, so load_amber selects EXPLICIT solvent automatically:
#   * PME electrostatics, GBSA off,
#   * box read from the rst7 and stored as reduced lattice vectors,
#   * 1.0 nm cutoff (override with nonbonded_cutoff=...),
#   * every molecule (solute + each water) gets a FREE 6-DOF root.
dih_classifier = robosample.AmberDihedralClassifier()
context = robosample.Context(args.name, args.seed, dih_classifier)
context.load_amber(args.prmtop, args.inpcrd)

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
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
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
sele = context.build_flexibilities(bonds, robosample.rb.JointType.Torsion, False)

# w = context.add_robotic_world(sele)
# w.add_sampler(
#     timeStep=0.002,
#     mdSteps=50,
#     acceptRejectMode=robosample.rb.AcceptRejectMode.MetropolisHastings,
#     use_nuts=False,
#     use_fixman=True,  # required for rigorous Boltzmann sampling of a constrained world
#     distort_option=None,  # robosample.rb.DistortOption.NMA
#     nma_bias_scale=0,  # 1
# )

# # minimized Ground-frame coords (nm), same order add_sampler/OpenMM use:
# st = context.system_topology
# pos_flat = []
# for a in range(st.num_atoms):
#     pos_flat += [st.atoms_x[a], st.atoms_y[a], st.atoms_z[a]]

# omega2 = w.set_nma_soft_mode_from_hessian(pos_flat, h=1e-5)
# print("softest internal mode omega^2 =", omega2)

# Root mobility is a per-world property handled inside add_ncmc_world: it detects
# water and WELDS all solvent AND all solute roots to Ground (only the solute's
# internal torsions move, plus the alchemical stride-through). Solvent relaxation
# and overall ergodicity are owned by the separate full-atom MD Gibbs block. No
# manual per-molecule weld loop is needed (Context.set_root_mobility no longer
# exists).
context.add_ncmc_world(
    sele,
    # 2 fs. The internal-coordinate velocity corrector must CONVERGE for the step
    # to sample: at too-large a dt it cannot find its fixed point, the step is
    # taken but pumps energy, and Metropolis then rejects it (coordinates stay
    # put). If you see "[verlet] ... corrector did not converge" warnings with
    # everything rejected, reduce this further; raise it only while moves keep
    # being accepted.
    # 1 fs. Historically (Construction I, endpoint-DeltaH) per-move acceptance was
    # set almost entirely by the integrator drift, which scales ~dt^2 over the
    # freed-water DOF: dH(5fs)=+70, dH(2fs)=+15, dH(1fs)=+2.5 kJ/mol -> acceptance
    # 0%, 0%, ~55%, and ncmc_steps 20->40 @ 1fs drove acceptance back to 0 (bath
    # shadow work is extensive in the propagated-DOF count, docs/specs/
    # ncmc-explicit-solvent/00-diagnosis-and-scaling.md). use_metropolized_inner
    # below (Construction II) removes that bath term from the OUTER acceptance
    # entirely, so this dt is no longer acceptance-critical the same way; kept at
    # 1fs as a conservative default for the inner GHMC's own per-substep accept.
    timestep=0.005,
    ncmc_steps=100,  # protocol length; under Construction II more steps mainly
    # buys smoother reorganization work, not more shadow-work exposure.
    hold_fraction=0.1,
    use_fixman=True,  # rigorous Boltzmann sampling of the constrained world; ALSO
    # required correctness-wise for Construction II's inner GHMC accept (F1/F2:
    # the inner H_lambda must include U_F, docs/specs/ncmc-explicit-solvent/
    # 10-acceptance-construction.md Sec.3 NOTE F2).
    # All solvent is welded here (no Free-rooted shell waters), so there is no water
    # libration to mass-scale away; physical masses (None) are the default. A
    # kinetic-metric mass_scale on the moving (Torsion) joints can still raise the
    # stable dt ~sqrt(scale) with no configurational bias if needed.
    mass_scale=None,
    accept_reject_mode=robosample.rb.AcceptRejectMode.MetropolisHastings,
    # Solvent-relaxing NCMC: the welded waters' atoms are advanced in flat
    # Cartesian space by OpenMM-force velocity-Verlet INSIDE the proposal, so the
    # cage relaxes during the lambda stride instead of being a rigid wall.
    relax_solvent=True,
    # Construction II (Metropolized-dynamics NCMC, docs/specs/ncmc-explicit-solvent/
    # 10-acceptance-construction.md): each fixed-lambda propagate substep
    # (including the relax_solvent Cartesian-Verlet steps) is Metropolized
    # against the full H_lambda, so the bath's shadow work is absorbed into inner
    # rejections instead of crushing the outer acceptance. This is the fix for
    # the near-zero explicit-solvent NCMC acceptance Construction I hit above --
    # relax_solvent is now pure benefit (lowers reorganization work) instead of
    # self-defeating (every relaxation step charged as shadow work).
    use_metropolized_inner=True,
)


# ---- Run --------------------------------------------------------------------
# initialize() runs an O(N^2) startup clash scan (now minimum-image aware). On a
# large solvent box this takes a moment; if your box is not yet minimized and the
# scan flags clashes, set ROBO_ALLOW_BAD_START=1 in the environment to downgrade
# the hard error to a warning, or minimize first.
context.initialize([300])
context.run_rex(args.equil_steps, args.prod_steps, args.write_freq, True)
