# NCMC explicit-solvent acceptance — Reproducer and oracles

## 1. Scope
A minimal reproducer wired to the existing `dH/work/gap` logging that **attributes** a 0-acceptance run to shadow-work vs reorganization vs corrector-nonconvergence, and the checkable oracles. Terms from `00-diagnosis-and-scaling.md` §2.

## 2. Diagnostic decomposition (mission Q4)
Per move, from `ncmcMove`'s log plus one added field, compute:
- `W` = `work` (protocol/reorganization work).
- `ΔU_F = U_F(q_end) − U_F(q_start)` = `state_.energy.fixman` at end minus at start (available: set by `reinitialize` and `currentTotalEnergy`; **the reproducer requires logging `fixman_start`/`fixman_end` and `logSineSqr_start`/`_end` — a diagnostic logging addition, not engine physics**).
- `ΔJ = −½RT·[ln sin²γ₂(q_end) − ln sin²γ₂(q_start)]`.
- **Shadow work** `Q = ΔH − W − ΔU_F − ΔJ` (equivalently `−gap − ΔU_F − ΔJ`). Uses `ΔH` and `work` already logged.
- **Corrector health**: count of `[verlet] ... corrector did not converge` lines during the move, and `checkReversibility(model,state,bridge,constraints, ncmc_steps, dt)` residual at the move's start geometry (non-destructive, `RobotIntegrator.hpp:471`).

**Attribution rule:**
- `Q ≫ W`, reversibility residual small, no warnings → fundamental symplectic shadow (S-a). Fix: Construction II or smaller `dt`/better integrator.
- `Q ≫ W`, warnings fire / residual large → corrector non-convergence (S-b). Fix: corrector (`20-inner-integrator.md`) or lower `dt`.
- `⟨ΔH⟩` grows ~linearly with `ncmc_steps` at fixed `dt` → drift (S-a/S-b/S-c), fixable/Metropolizable; ~flat → bounded symplectic floor.
- `W ≫ Q` → reorganization-dominated. Fix: more stages / more inner solvent relaxation / smaller Region A / minimum-dissipation protocol.

Predicted attribution for the failing config (relax_solvent, ~thousands of water DOF, dt=1fs, n=20): shadow-dominated and drift-positive (`Q ≫ W`, `⟨ΔH⟩` rising with `n`), with an S-b component from the solute corrector.

## 3. Minimal reproducer (system + sweeps)
System: tip3p/2ala (`tip3p/2ala.prmtop`, `.rst7`), the exact `run.py` config. Harness (Python driver + C++ hooks) SHALL run these sweeps, each ≥ a few hundred moves for a stable mean:
- **Sweep A (bath-scaling, the decisive test):** hold Region A and dt fixed; grow the water box (pad shells beyond A's solvation shell). Record `⟨ΔH⟩`, `⟨Q⟩`, `⟨W⟩`, acceptance under Construction I and Construction II.
- **Sweep B (drift signature):** fix dt; vary `ncmc_steps ∈ {10,20,40}`; record `⟨ΔH⟩(n)`.
- **Sweep C (corrector isolation):** toggle `relax_solvent`; and vary corrector tol; record warning counts, reversibility residual, `Q`.
- **Sweep D (Jarzynski/matched-pair):** Construction II with a genuine `π_λ`-invariant `Φ` vs a deliberately non-invariant `Φ` (plain Verlet + accept-on-`W`); record running `⟨e^{−βW}⟩` and a torsional free-energy observable.
- **Sweep E (Construction-I-vs-II marginal equality — the reference-free correctness gate, F1):** run Construction I (endpoint-ΔH, small dt so acceptance is nonzero) and Construction II on the SAME small system to collect the sampled torsional marginal (φ/ψ histogram) of the moving region. Both target `exp(−βH_{λ=1})`, so their marginals SHALL agree within statistics. Run on the OpenMM-free analytic bridge (à la `TestNCMCWork`) so it needs NO external ΔF reference. Include an arm where the Construction-II inner accept deliberately OMITS the `U_F`/pitch terms — that arm SHALL disagree.

## 4. Oracles

**PRECONDITION (runtime guard, not a test).**
- The solute atom block is contiguous (`context.py` already checks). Region A ⊂ that block.
- `λ=1` at both protocol endpoints (`NcmcProtocol.hpp`). Runtime guard.

**INVARIANT INV0 (reference-free correctness discriminator — PRIMARY, F1).** The Construction-I-vs-II marginal-equality test (Sweep E). Both constructions target `exp(−βH_{λ=1})`, so the moving region's sampled torsional marginal SHALL match between I and II within statistics; the `U_F`/pitch-omitting arm SHALL measurably diverge. This is the ONLY reference-free oracle that catches the canonical Fixman/pitch-in-inner-accept omission (L1, L2, INV2 are all blind to it — see the caveats below), and it requires no external ΔF reference. NORMATIVE deliverable alongside the fix.

**LEMMA L1 (flat-λ gate).** With `λ≡1` throughout, `W=0` bitwise and Construction II acceptance reduces to plain torsional-HMC `ΔH`. Expected: `W==0.0` exactly; end-state bit-identical to alchemy-free HMC. (Extends `TestNCMCWork.FlatLambdaOneGivesZeroWorkAndPlainHmcDeltaH`.) CAVEAT (F5): the reference HMC for the Construction-II arm SHALL include the SAME Fixman/pitch terms as `currentTotalEnergy` (`World.cpp:1788`). The existing `TestNCMCWork` reference compares `V+K` only; inheriting that Fixman-less reference would let a `U_F`-omitting inner pass L1 trivially. Specify a Fixman-complete reference.

**LEMMA L2 (Jarzynski cyclic — WEAK self-consistency gate only, F4).** For the cyclic palindromic `λ:1→0→1` protocol, `⟨e^{−βW}⟩ = e^{−βΔF} = 1`. NOTE this does NOT discriminate correct from biased: because `U_F`/pitch are λ-independent, `W` is identical whether or not the inner kernel includes them, so `⟨e^{−βW}⟩=1` holds relative to WHATEVER stationary density the inner preserves — including the wrong Fixman-less one (F1). Additionally the exponential-average estimator is dominated by rare large-negative-`W` tails, so `Var(e^{−βW})` can be effectively undefined at achievable `n_moves` and L2 can pass a biased sampler for lack of statistical power. Treat L2 as a weak consistency check with an honest caveat; put correctness weight on INV0 (reference-free) and INV1 (reference-based).

**INVARIANT INV1 (matched-pair, must-check — the correct-vs-plausible discriminator).** A known-ΔF observable (2ala φ/ψ torsional free-energy difference between two basins, or a decouple-only ΔF against MBAR reference) SHALL be recovered by {Construction II + `π_λ`-invariant `Φ`}. The *same* code with {accept-on-`W` + plain Verlet `Φ`} (invariance deliberately broken) SHALL measurably miss it. A construction that passes both is not actually Metropolizing the inner step. Must exercise angle+torsion coupling together (see INV3).

**INVARIANT INV2 (bath-scaling).** In Sweep A: Construction II `⟨W⟩`, `⟨Q_outer=0⟩`, and acceptance SHALL be ~invariant to bath padding beyond the solvation shell (within statistics); Construction I `⟨ΔH⟩` SHALL grow with `N_water`. This is the property that certifies FFAR1-scale viability; a Construction-II implementation whose acceptance still falls with box padding has a residual bath term (e.g. inner `H_λ` missing a per-water term) and fails.

**INVARIANT INV3 (physics gate — the move must actually cross a real barrier).** The reproducer SHALL include an alanine-dipeptide C7eq↔C7ax (or cis↔trans) transition that is gated by a 1–4 clash requiring **angle + torsion together**: a move that only relaxes torsion without resolving the clash must not spuriously accept. This guards against a construction that "accepts" by sampling a distorted (electrostatics-annihilated, `30-region-and-protocol-policy.md`) internal landscape. Pass: transition rate exceeds plain-Cartesian MD at matched wall-clock while the recovered basin free-energy difference matches the MBAR reference within statistics.

**LEMMA L3 (shadow isolation).** In Sweep C, turning off `relax_solvent` (removing `3·N_water` propagated DOF) SHALL drop `⟨Q⟩` by an amount consistent with the removed-DOF count (`Q` roughly proportional to propagated-DOF count in the drift regime). Confirms `Q` is bath-DOF-driven, not reorganization.

## 5. Pass condition
The reproducer passes when it (a) *attributes* the reported 0-acceptance run via §2 (expected: shadow-dominated, drift-positive), (b) satisfies **INV0** (reference-free), L1, INV2, INV3, L3 for the Construction-II implementation, plus INV1 once the deferred ΔF reference exists, and (c) demonstrates INV2's separation between Construction I (acceptance falls with box padding) and Construction II (acceptance flat).

NOTE (F1, ordering under the deferred ΔF reference): with INV1/INV3 positive arms deferred, **INV0 is the gating correctness oracle** — it is the only reference-free check that fails on the Fixman/pitch-in-inner-accept omission, which L1/L2/INV2 all pass. Do NOT declare the fix correct on L1+L2+INV2 alone; a `U_F`-omitting inner ships green under those three. A Construction-II build is accepted as correct when it passes INV0 (marginal matches Construction I; the `U_F`-omitting arm diverges) and is box-padding-invariant (INV2); INV1/INV3 upgrade the evidence from reference-free-consistency to absolute-ΔF once the reference is generated.

## 6. Touch list (reproducer only)
- Add `fixman_start/end`, `logSineSqr_start/end` to the `ncmcMove` diagnostic log (diagnostic, not physics).
- New test/driver files (Python harness + a C++ `TestNcmcExplicitSolvent` exercising Sweeps A–D against the OpenMM-backed `World`, plus L2/INV1 on the OpenMM-free analytic bridge à la `TestNCMCWork` where possible).
