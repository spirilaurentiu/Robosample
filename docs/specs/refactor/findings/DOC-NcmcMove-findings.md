# DOC-NcmcMove findings

Style: house `//` contract blocks (see DOC-World-findings.md). No code changed.

## State
world/sampler/NcmcMove.cpp (ncmcMove, protocolLambda, ncmcTeleportRoot,
ncmcApplyTroughTeleport, ncmcInnerGhmcStep, configureNcmc) and the World.hpp
declarations (266-276, 604-622) carried full contract documentation from
SPLIT-W4. No residual edit required.

## Verified hypotheses (work accounting -- Critical if wrong)
- Two acceptance constructions, documented as mutually exclusive (CLAIM C1):
  * Construction I (default, useMetropolizedInner=false): accept on the endpoint
    Hend - Hstart (plain Verlet dH); protocol work is DIAGNOSTIC ONLY, not used
    in acceptance (NcmcMove.cpp:336-348, 448-455).
  * Construction II: each fixed-lambda substep is Metropolised against the full
    H_lambda by ncmcInnerGhmcStep, absorbing shadow work internally; the OUTER
    move accepts on the protocol work ALONE, a = min(1, exp(-beta*W)), computed
    by metropolis(0.0, work) (:401-461). Verified: work accumulates only at the
    PERTURB step (lambda change at fixed q, work += V(lam_new) - V(lam_old),
    :380-390); fixed-lambda V drift is heat, not work (:416).
- Inner GHMC reversibility: one fixed-lambda Verlet substep, Metropolis
  accept/reject with momentum flip on reject; non-convergence is an automatic
  reject (F3). Preserves pi_lambda at each lambda (:153-290). Documented.
- Trough teleport: discrete lambda=0 rigid reposition with KE-preserving velocity
  co-rotation; proposal and its correction documented separately (World.hpp
  :606-607, NcmcMove ncmcApplyTroughTeleport).

## Shadow-work sensitivity (context, not a doc claim)
- MEMORY (NCMC explicit-solvent campaign) traced near-zero acceptance to
  integrator shadow work ~ dt^2 * nDOF, not reorganization. Construction II is the
  design response (absorb shadow work in inner rejections). Recorded for the
  reviewer; the doc does not assert this resolves the acceptance problem.

## Coverage gap (ticket-flagged)
- TestNCMCWork (SLOW, tests 9-16) is the work-accounting oracle. TestNcmcTeleport
  (SLOW) covers the teleport. TestNcmcExplicitSolvent is largely stubbed
  (TESTS.md D-T2) and is NOT evidence. The explicit-solvent path is a recorded
  coverage gap.

## Assumed notes
None.
