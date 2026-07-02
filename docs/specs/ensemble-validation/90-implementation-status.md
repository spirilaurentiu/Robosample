# Ensemble-validation suite — implementation status & handoff

Progress/handoff record (deliberately separate from the tier specs, which state
*what is correct*, not *what is coded yet* — CLAUDE.md spec style). Last updated
2026-07-02. Read this first if you are resuming the suite.

**One-line state:** every test in all three tiers is WRITTEN. Fast/always-on gates
pass. The **slow tier is not green**: Tier-2 conformational tests fail on untuned
test calibration (see §4). The user has DEFERRED the slow re-run; do not launch
hours-long sampling runs without checking in.

---

## 1. Status by tier

| Tier | File(s) | Written | Verified |
|------|---------|---------|----------|
| Tier 0 (C++ kinetic) | `tests/TestEquipartition.cpp`, `tests/TestEnsembleValidation.cpp`, `tests/RobotBuilders.hpp`, `tests/StatTest.hpp` | ✅ | ✅ **PASS** (incl. slow T0.2) |
| Tier 1 (PE ladder) | `tests/test_ensemble_pe_ladder.py`, `tests/ensemble_stats.py` | ✅ | ✅ **PASS** 7/7 (~28 min) |
| Tier 2 (torsion) | `tests/test_torsion_conformational.py` | ✅ | ⚠️ **12 pass / 4 fail** (see §4) |
| Stats self-checks | `tests/ensemble_stats.py` | ✅ | ✅ **PASS** 8/8 |
| Fixture PE parity | `tests/test_openmm_potential_energy.py` CASES (`ethane`/`butane`/`2butanol`) | ✅ | ✅ **PASS** 3/3 |

Not implemented, on purpose: **T0.4** (Python KE smoke) — blocked on Open Q1
(`00-foundations.md` §8, moves.csv KE-column semantics); leave unbuilt until Q1 is
answered.

---

## 2. What was implemented (inventory)

- **T0.1** `Equipartition.MultiBodyBentChain` + `...MultiBodyBentChain_TRatio`
  (`TestEquipartition.cpp`). Dense-M fixture `buildBentTorsionChain` (Free root + 8
  alternating-bend Torsions, n_dof=14, seed `0x70`) added to `RobotBuilders.hpp`,
  shared with T0.2.
- **T0.2** `EnsembleValidation.KEGammaMultiBodyBentChainWithWrongDofControl`
  (`TestEnsembleValidation.cpp`). KE~Gamma(n_dof/2,kT) + the discriminating ±1-DOF
  *rejection* control. Sub-sampled per-bin Gamma weights (`gammaSubsampledWeights`).
- **T0.3** `Equipartition.ConstrainedRingFixmanCoupling` (`TestEquipartition.cpp`):
  ties the equipartition denominator (nu − n_C) to the same n_C the Fixman loop
  term (`calcConstraintLogDet`) uses.
- **Tier 1** `test_ensemble_pe_ladder.py`: rungs A (native OpenMM Langevin) / B
  (Robosample cartesian) / C (OpenMM+GCHMC composite); histogram equality B≡A, C≡B;
  moments; per-rung Shirts slope.
- **Tier 2** `test_torsion_conformational.py`: T2.0 self-consistency + dihedral-only
  control; T2.1 ethane 3-fold equipopulation; T2.2 butane anti/gauche; T2.3 2-butanol
  CC/OH coupling; plus native-MD comparisons and a `_rotate_to_dihedral` vs mdtraj
  geometry self-check.
- **`ensemble_stats.py`**: Shirts slope estimators A (weighted-linear) & B
  (logistic-MLE), N_eff via autocorrelation, moves.csv parsing, + 8 self-check tests.
- **Fixture registration**: `ethane`/`butane`/`2butanol` added to the CASES table of
  `test_openmm_potential_energy.py` for single-point PE parity against native OpenMM
  (the meaningful, external-oracle half of the Tier-2 touch-list SHOULD). They were
  deliberately NOT added to `tests/loader_differential_lib.py`: that gate pins
  *pre-rewrite* loader goldens and `_generate.py` rewrites ALL goldens at once, so a
  new-molecule golden there is either unsafe or self-referential.

---

## 3. How to run

Env: a `robo_cuda*` conda env must be active (`CONDA_PREFIX` set). Build the C++ test
tier once: `cmake --build --preset cuda-tests`.

- **Tier 0 fast (always-on):**
  `ctest --test-dir build/cuda-tests -R 'Equipartition|EnsembleValidation' --output-on-failure`
- **Tier 0 slow (T0.2 control, PE-slope):** prefix `env ROBOSAMPLE_SLOW_TESTS=1`.
  (`EnsembleValidation.PEobeysEnsembleSlope` alone is ~13 min.)
- **Python tiers** (need `ROBOSAMPLE_REQUIRE_OPENMM=1` to hard-import, `ROBOSAMPLE_SLOW_TESTS=1`
  to un-skip sampling): e.g.
  `env ROBOSAMPLE_SLOW_TESTS=1 ROBOSAMPLE_REQUIRE_OPENMM=1 python3 -m pytest tests/test_torsion_conformational.py -q`
  Tier 1 ≈ 28 min; Tier 2 ≈ 2h46m. The stats self-checks
  (`python3 -m pytest tests/ensemble_stats.py`) need neither env var.
- **Authoritative gate:** `nox -s tests` (sets both env vars, builds cuda-tests).

**Environment gotchas (not failures):**
- OpenMM CUDA prints `Error deleting array … CUDA error` at process teardown — harmless.
- pytest `-q` summary is on stdout; the teardown noise is on stderr. Judge by the
  `N passed / N failed` line, not stderr.
- CUDA force reductions are **not** bitwise-reproducible run-to-run, so GCHMC
  trajectories diverge across runs even with a fixed seed (the World RNG itself IS
  deterministic, `World.cpp:205`). Low-N_eff statistical tests are therefore flaky.

---

## 4. Open work — the 4 Tier-2 failures (RESUME HERE)

First full slow run: **12 passed, 4 failed**. All four are conformational-population /
self-consistency checks the (crashed) prior session wrote but never ran, so their
sampling budgets and bands were never tuned against a real trajectory. Diagnoses below
are **UNCONFIRMED** — they need the deferred re-run to settle. Do NOT "loosen bands
until green": that would hollow out the validation. The discriminating experiment is
whether **ethane's exact 3-fold equipopulation** (a symmetry oracle that MUST hold for
a correct sampler) converges to 1/3 with more sampling — if yes, all four are test-side;
if no, there is a real sampler bias worth escalating.

| Test | Observed | dof / N_eff | Leading diagnosis |
|------|----------|-------------|-------------------|
| `test_t2_0_2butanol_self_consistency` | χ²=**9055** vs crit 160 | 99 / 7500 | Coarse-bin reference bias (see below) — the strongest lead; likely a TEST bug, not the sampler |
| `test_t2_0_butane_self_consistency` | χ²=101.8 vs crit 74.9 | 35 / 1154 | Undersampling (slow anti↔gauche crossing) |
| `test_t2_1_ethane_three_basins_equipopulated` | props [0.42, 0.32, 0.27] | — / 750 | Undersampling + CUDA-nondeterminism flakiness |
| `test_t2_2_butane_anti_gauche_ratio_matches_native_md` | GCHMC [0.17,0.05,0.78] vs native [0.33,0.33,0.34] | — | Stuck/undersampled (both sides) |

**Ruled out:** axis-transpose for the 2D 2-butanol case — `mol.torsions` order
(`_molecule_spec`, dict insertion order phi_cc→phi_oh) matches the `["phi_cc","phi_oh"]`
order passed to `_chi2_gof_nd`, so the histogram and the reference grid axes align.

**Proposed fixes when resumed:**
1. **2-butanol χ²=9055 (do this first — cheap, no re-sample of the sampler needed to
   test the hypothesis):** `_chi2_gof_nd` weights the reference by `exp(-β·U)` at the
   **bin CENTRE** (`test_torsion_conformational.py:306`). With coarse 36° 2D bins
   (`_N_BINS_2D=10`) near a steep torsion wall this is badly biased — the Boltzmann
   weight is dominated by the bin's low-energy edge, not its centre. This is exactly
   the bin-centre trap the KE-Gamma test avoids by sub-sampling
   (`00-foundations.md` §5.2). Fix: evaluate `_reference_pe_grid` on a finer sub-grid
   and aggregate the mean Boltzmann weight per coarse bin (mirror
   `TestEnsembleValidation.cpp::gammaSubsampledWeights`). Fine-binned 1D ethane passing
   while coarse-binned 2D 2-butanol fails by 56× is consistent with this being the cause.
2. **Population tests (ethane/butane):** raise N_eff — more `prod_rounds` and/or larger
   `mdSteps` per round for better barrier crossing (`_T2_0_PROD_ROUNDS`,
   `_T2_0_MDSTEPS`, `_NATIVE_BASIN_N_SAMPLES`). Then calibrate the equality bands from
   the symmetry-exact ethane case (spec `30-tier2` PRECONDITION) and reuse.
3. Only after 1–2 converge and ethane equipopulates: if a residual bias persists, it is
   a real finding — escalate, do not paper over.

An ethane 60k-round convergence experiment was started to settle this but cancelled
before finishing (it was over-sized; ~20–25k rounds is enough and finishes in ~15–20
min). See `[[ensemble-validation-suite]]` memory for the same summary.
