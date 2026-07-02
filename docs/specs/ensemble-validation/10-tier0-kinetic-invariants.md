# Tier 0 — engine-level kinetic invariants (C++ GoogleTest)

Scope: momentum-space only (foundations §5.2–5.3). No trajectory, no accept/reject,
no potential. Fast; always-on gate. Extends TestEquipartition.cpp and
TestEnsembleValidation.cpp.

## Claim
The internal articulated-body engine draws u = √(RT) M(φ)^{-1/2} g so that
2KE = uᵀM(φ)u = RT·gᵀg, giving ⟨2KE⟩ = RT·n_dof and KE ~ Gamma(n_dof/2, kT),
independent of φ and of the joint topology, with n_dof counted per foundations §2.

## What changes vs today
Existing coverage is 6–8 DOF toys (FreeBody, Free+2Torsion, ConstrainedRing).
The gap: no *realistic multi-body internal-coordinate* fixture where M(φ) is dense
and strongly configuration-dependent. Add one.

### T0.1 — extend to a multi-body internal chain  [extends TestEquipartition.cpp]
- SHALL add a fixture of ≥ 8 Torsion joints on a bent chain (reuse the 90°-bend
  coupling idea from TestFixmanBoltzmann.cpp::twoTorsionChain so M is dense and
  off-diagonal), free-6DOF root → n_dof = 6 + n_torsions.
- SHALL assert ⟨2KE⟩/n_dof = kT within 4·stderr (checkEquipartition pattern).
- SHALL assert the T-ratio invariant: run at T1, T2, assert
  ⟨KE(T1)⟩/⟨KE(T2)⟩ = T1/T2 within propagated 4·stderr (DOF-free; foundations §5.3).

### T0.2 — KE-Gamma on the same multi-body fixture  [extends TestEnsembleValidation.cpp]
- SHALL histogram seedMomenta→realizeVelocity→calcKineticEnergy KE and χ²-test
  against Gamma(n_dof/2, kT) at alpha=1e-4, sub-sampling per-bin weights (never
  bin-centre) exactly as KEisMaxwellBoltzmann does.
- SHALL, as a discriminating control, assert the SAME histogram REJECTS
  Gamma(n_dof±1 /2, kT) (a wrong DOF count) at alpha=1e-4 — otherwise the test
  cannot distinguish a correct M^{-1/2} draw from a mis-counted one.

### T0.3 — constraint DOF bookkeeping  [already in ConstrainedRing; keep + assert coupling]
- SHALL keep the n_dof = nu − n_C ring test AND assert the same n_C is what
  calcConstraintLogDet removes, i.e. equipartition uses nu−n_C and Fixman's
  lnDetZ (World.cpp:860) is nonzero for the same ring. (Cross-checks that the two
  DOF accountings agree — foundations §2.)

### T0.4 (optional) — Python-level equipartition smoke  [new pytest; gated on Open Q1]
With `World.n_dof` now exposed to Python (foundations §2, Open Q2 resolved), an
end-to-end smoke check becomes possible:
- SHOULD run a short torsional-world sim on ala-dipeptide, read per-round KE from
  {base}.moves.csv, and assert ⟨2KE⟩/world.n_dof = kT within 4·stderr(N_eff).
- NOTE: subordinate to the authoritative C++ T0.1/T0.2 gates (it cannot test the KE
  *shape* without the C++ seedMomenta hook). BLOCKED until Open Q1 confirms the
  moves.csv KE column is the freshly-resampled KE, not end-of-trajectory — otherwise
  the equipartition target does not apply to the logged value.

## Correctness conditions
- PRECONDITION: KE samples and n_dof come from the SAME engine (§2).
- INVARIANT (T0.1): ⟨2KE⟩/n_dof = kT within 4·stderr, at ≥2 temperatures with the
  ratio equal to T1/T2. Fails if M^{-1/2} draw, mass metric, or DOF count is wrong.
- LEMMA (T0.2): χ²_stat < χ²_crit(nbins−1, 1e-4) against Gamma(n_dof/2, kT);
  and χ²_stat > χ²_crit against Gamma((n_dof−1)/2, kT). Expected χ²≈nbins−1 under H0.
- INVARIANT (T0.3): equipartition denominator and Fixman n_C are identical.

## Touch list
- tests/TestEquipartition.cpp (new multi-body fixture + T-ratio).
- tests/TestEnsembleValidation.cpp (new KE-Gamma fixture + wrong-DOF control).
- No production code. Conventions at risk: n_dof counting (§2), Free vs welded
  root mobility, M(φ) vs M^{-1} in the draw.

## Verification plan
Analytic oracle throughout (Gamma/χ² are closed-form). The discriminating check is
T0.2's wrong-DOF rejection: a plausible-but-biased engine that draws g~N(0,I) but
forgets the M^{-1/2} factor still yields ⟨2KE⟩≈RT·n_dof-ish on some fixtures but
fails the Gamma SHAPE — so the shape test, not the mean, is the real gate.
