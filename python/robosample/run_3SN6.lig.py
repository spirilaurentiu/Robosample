"""Per-body force monitoring for 3SN6 (human beta2-AR) WITH agonist, in a nanodisc.

One rigid body per TM helix (TM1..TM7): weld the whole receptor and free a SINGLE
phi joint in each of the six inter-helix loops -> 7 rigid TM bodies. The reporter
(enable_reaction_reporter) self-derives those bodies and prints one net applied
spatial force (bodyForceG: force+torque about the body origin, in Ground) per body,
per DCD frame, to `<name>.<replica>.reactions.csv`. Identify each body from the
CSV's atom_idx (-> residue -> which TM). See
docs/specs/reaction-force-monitoring.md and
docs/specs/gpcr-world-design/20-monitoring-selections-beta-adrenergic.md (5A).

Receptor-only nanodisc (molecule 0 = receptor; molecules 1,2 = MSP belt; then
lipids; then agonist). Implicit solvent (no box) -> load_amber selects GBSA-OBC2.

Run from the repo root:
    python3 python/robosample/run_3SN6.lig.py
"""

import argparse

import parmed as pmd

import robosample

# ---- Per-model configuration -----------------------------------------------
NAME = "3SN6.lig"
PRMTOP = "examples/febs/3SN6.lig.nanodisc.prmtop"
INPCRD = "examples/febs/3SN6.lig.nanodisc.min.rst7"

# Ordered, contiguous LABELED rigid bodies (model/prmtop numbering; beta2AR / 3SN6).
# ONE dihedral is flexed at the START of each segment after the first (phi; psi at
# prolines, whose phi is ring-locked), so each segment becomes one rigid body.
# TM cores + microswitch anchors are from docs/specs/gpcr-world-design/
# 20-monitoring-selections-beta-adrenergic.md (verified); loop/H8 boundaries are
# approximate -- edit freely. Insert a (label, lo, hi) row to add a sub-helix body.
# Labels print in the run log; tag CSV rows with examples/febs/label_reactions.py.
SEGMENTS = [
    ("TM1", 1, 41),
    ("TM2_EC", 42, 49),
    ("Na_D2.50", 50, 51),  # Na+ pocket microswitch (D2.50)
    ("TM2_IC", 52, 74),
    ("TM3", 75, 100),
    ("DRY_ionic", 101, 120),  # DRY / ionic-lock microswitch (R3.50)
    ("TM4", 121, 168),
    ("TM5", 169, 237),
    ("TM6_EC", 238, 255),  # TM6 extracellular half
    ("CWxP_toggle", 256, 259),  # CWxP toggle + P6.50 kink (W6.48)
    ("TM6_cyto", 260, 275),  # TM6 cytoplasmic half (G-protein swing)
    ("TM7", 276, 292),
    ("NPxxY", 293, 297),  # NPxxY microswitch (Y7.53)
    ("H8", 298, 312),  # amphipathic helix 8
]

# ---- Run parameters (override on the CLI) ----------------------------------
parser = argparse.ArgumentParser(description=f"Per-body force monitoring: {NAME}")
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("--equil", type=int, default=100, help="equilibration rounds")
parser.add_argument(
    "--prod", type=int, default=1000, help="production rounds (CSV only in production)"
)
parser.add_argument(
    "--write_freq", type=int, default=10, help="DCD + force-CSV write frequency"
)
args = parser.parse_args()

# ---- Build the context ------------------------------------------------------
# No periodic box -> load_amber selects IMPLICIT solvent (GBSA-OBC2).
dih_classifier = robosample.AmberDihedralClassifier()
context = robosample.Context(NAME, args.seed, dih_classifier)
context.load_amber(
    PRMTOP,
    INPCRD,
    use_gbsa_obc2=True,
    nonbonded_method=robosample.rb.NonbondedMethod.CutoffNonPeriodic,
    nonbonded_cutoff=1.0,
)
context.set_enforce_periodic_box(False)

# ---- Receptor residues (molecule 0), for junction picking + atom mapping ----
# standard_dihedral_bonds has no residue column (residue_idx is a -1 placeholder),
# so we map residues -> atoms with ParmEd. The receptor is molecule 0, whose
# molecule-local atom index == global prmtop atom index (== ParmEd a.idx).
parm = pmd.load_file(PRMTOP)
receptor = []
_started = False
for _r in parm.residues:
    _is_aa = len(_r.name.strip()) == 3 and _r.name.strip() not in ("ACE", "NME", "NHE")
    if _is_aa:
        _started = True
        receptor.append(_r)
    elif _started:
        break  # receptor chain (molecule 0) ended
resname = {r.number: r.name.strip() for r in receptor}
atoms_by_res = {r.number: [a.idx for a in r.atoms] for r in receptor}

# ---- One dihedral per segment boundary -> len(SEGMENTS) labeled rigid bodies ----
# Flex the phi of each segment's first residue (psi if that residue is proline,
# whose phi is ring-locked). build_flexibilities welds everything else, so each
# segment becomes one rigid body monitored by the reporter.
boundary_resids = [seg[1] for seg in SEGMENTS[1:]]
phi_atoms = set()
psi_atoms = set()
for _b in boundary_resids:
    (psi_atoms if resname.get(_b) == "PRO" else phi_atoms).update(
        atoms_by_res.get(_b, [])
    )

df = context.standard_dihedral_bonds
_PHI = robosample.DihedralType.PROTEIN_PHI.value
_PSI = robosample.DihedralType.PROTEIN_PSI.value
mask_phi = (
    (df["dihedral_type"] == _PHI)
    & df["atom1_idx"].isin(phi_atoms)
    & df["atom2_idx"].isin(phi_atoms)
)
mask_psi = (
    (df["dihedral_type"] == _PSI)
    & df["atom1_idx"].isin(psi_atoms)
    & df["atom2_idx"].isin(psi_atoms)
)
bonds = df.loc[(df["molecule_idx"] == 0) & (mask_phi | mask_psi)]
sele = context.build_flexibilities(bonds, robosample.rb.JointType.Torsion, False)

print(f"[{NAME}] {len(bonds)} boundary joints -> {len(SEGMENTS)} labeled rigid bodies:")
for _lab, _lo, _hi in SEGMENTS:
    print(f"    {_lab:14s} residues {_lo}-{_hi}")

context.add_cartesian_world().add_sampler(
    timeStep=0.001,
    mdSteps=1000,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
    use_fixman=False,
)


# ---- Reporter world ---------------------------------------------------------
# Weld the nanodisc scaffold (MSP = molecules 1,2) to Ground as the fixed frame;
# receptor (0) and lipids (3+) keep Free roots (lipids must stay mobile).
# ORDER MATTERS: set_root_mobilities rebuilds and RENUMBERS the bodies, so it must
# precede add_sampler, and enable_reaction_reporter must run AFTER the final
# rebuild (else the reporter's body indices are stale -- World.hpp:427-431).
mob = [robosample.rb.JointType.Free] * context.system_topology.num_molecules
mob[1] = robosample.rb.JointType.Rigid
mob[2] = robosample.rb.JointType.Rigid

world = context.add_robotic_world(sele)
world.set_root_mobilities(mob)
world.add_sampler(
    timeStep=0.025,
    mdSteps=10,
    acceptRejectMode=robosample.rb.AcceptRejectMode.AlwaysAccept,
    use_nuts=False,
    use_fixman=True,
)
world.enable_reaction_reporter(
    report_free_bodies=False, include_openmm=True, include_reaction=True
)

# ---- Run --------------------------------------------------------------------
context.initialize([300])
context.run_rex(args.equil, args.prod, args.write_freq, True)
