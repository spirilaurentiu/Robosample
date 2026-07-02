# Tier 2 — torsion conformational analysis (pytest)

Scope: does the Fixman-corrected internal (GCHMC) world reproduce the correct
configurational torsion marginal. Molecules: ethane, butane, 2-butanol — fixtures
now present: `examples/{ethane,butane,2butanol}.{prmtop,rst7}`. External reference:
native OpenMM MD (matched FF/T/nonbonded).

## Claim (scoped by DOF-matching, foundations §4)
Fixman-ON GCHMC samples ρ(φ) ∝ exp(−βU(φ; q₀)) at the frozen geometry q₀
(§3, Kandel eq:13). Against fully-flexible Cartesian MD this matches the
torsion-shaping physics (dihedral term + 1-4/nonbonded) up to the residual
rigid-vs-flexible / extrinsic distortion (Kandel eq:18–19, König 2014). So the
comparison is QUALITATIVE (topology, symmetry, populations within a band) unless
DOF-matched (M-full or M-rigid, §4).

### T2.0 — self-consistency (exact, DOF-free) FIRST  [new file]
Before any external comparison, assert the internal sampler matches its OWN target:
- LEMMA: Fixman-ON GCHMC torsion histogram matches exp(−βU_full(φ; q₀)) — a scan
  of the FULL OpenMM potential along φ at frozen bond/angle geometry — at χ²
  alpha=1e-4. This is TestFixmanBoltzmann.cpp generalized to U≠0 with the OpenMM
  potential as the reference density. It uses NO external simulator and is the
  clean regression gate.
- CONTROL: the same histogram REJECTS the bare dihedral-term-only scan
  exp(−βU_dihedral(φ)) whenever 1-4/nonbonded contributes (ethane H···H, butane)
  — demonstrating why the "bare-scan" reference is invalid (framing #4).

### T2.1 — ethane (symmetry; near-exact)  [new file]
- Three staggered minima at φ ≈ ±60°,180°, equipopulated by C3 symmetry.
- INVARIANT: the three basin populations are equal within counting error
  (multinomial 4·stderr on N_eff). Symmetry is FF-independent → a hard oracle,
  robust to the rigid-vs-flexible residual (which is itself 3-fold symmetric).
- SHOULD: GCHMC basin populations match native OpenMM MD within the same band.

### T2.2 — butane (anti/gauche; qualitative)  [new file]
- Anti (180°) global min, two gauche (±60°) minima.
- INVARIANT: GCHMC recovers 3 minima at the right locations and anti > gauche.
- LEMMA: P(anti)/P(gauche) from GCHMC equals the native-OpenMM-MD ratio within a
  stated band (default 15% relative on N_eff), NOT an absolute literature number
  (populations are FF-dependent; the matched native run is the ground truth).

### T2.3 — 2-butanol (asymmetric + coupled O–H)  [new file]
- C–C rotamers inequivalent (asymmetric marginal); O–H rotamer distribution
  depends on the C–C rotamer (coupling).
- INVARIANT: the JOINT P(φ_CC, φ_OH) from GCHMC is asymmetric and shows C–C/O–H
  coupling (mutual dependence detectable vs the product of marginals).
- LEMMA: joint histogram matches native OpenMM MD within band; the coupling
  (difference between joint and product-of-marginals) has the same sign/scale.
  This is the test that a naive per-torsion sampler cannot fake.

## Torsion extraction
z_matrix.py + standard_dihedral_bonds to identify φ_CC, φ_OH; compute dihedral
values from {replica}.dcd (prmtop atom order) via mdtraj. Same extraction on the
native OpenMM trajectory.

## Correctness conditions
- PRECONDITION: GCHMC and native runs share FF/T/nonbonded/constraints (§6).
- PRECONDITION: equality-band tolerances are set from the DOF-mismatch magnitude,
  not from round numbers — calibrate the band on ethane (where the exact answer is
  known by symmetry) and reuse it; document at the call site.
- INVARIANT (T2.0): internal self-consistency vs exp(−βU_full(φ;q₀)); rejects the
  dihedral-only scan.
- INVARIANT (T2.1): ethane 3-fold equipopulation within counting error.
- LEMMA (T2.2/T2.3): population ratios / joint match native MD within band on N_eff.

## Touch list
- Fixtures examples/{ethane,butane,2butanol}.{prmtop,rst7} — ALREADY PRESENT (no
  generation step needed). SHOULD single-point-PE-validate them by registering in
  the CASES table of tests/test_openmm_potential_energy.py and
  tests/loader_differential_lib.py before use.
- New tests/test_torsion_conformational.py.
- Reuse tests/ensemble_stats.py (Tier 1) for N_eff. Conventions at risk: DOF match
  (§4), Fixman ON, φ wrapping, which frozen geometry q₀ defines the self-consistency
  reference, Euler/orientational Jacobian OFF (World.cpp:254 useOrientationJacobian).

## Verification plan
The discriminating oracle is T2.0 + the T2.1 symmetry invariant: T2.0 catches any
Fixman/metric bug against an internal reference (no external simulator noise), and
ethane's exact 3-fold equipopulation catches a broken sampler even when
FF-dependent bands would mask it. T2.3's coupling term is the check a
plausible-but-biased independent-torsion sampler fails. The alanine cis↔trans
motivation (angle+torsion gated by 1-4 clash) is deferred to a follow-up: it
requires the composite OpenMM+GCHMC sampler (Tier 1 rung C), not a pure-torsion
world, because the 1-4 clash relief needs the flexible bond/angle DOF.
