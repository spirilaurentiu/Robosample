# Refactor Plan - Consistency Review

Status: final, 2026-07-12. Scope: the recovery docs (`../../architecture/`) and
the SPLIT / DOC / TEST / VERIFY ticket set in this directory. This review checks
internal consistency, blatant errors, and parallelizability. It does not re-derive
the physics or re-audit the source; it takes the ticket set as the object under
review.

## Verdict

The plan is sound and unusually rigorous. The core mechanisms - the
two-phase freeze, the comment-stripped-diff oracle, the `nm` public-symbol delta,
the per-ticket restatement of invariants, and the slow-test mutation-power gate -
are correct and mutually reinforcing. Each ticket is self-sufficient with
machine-checkable exit criteria, which is what makes the one-agent-per-ticket
execution model viable.

Three defect classes were found: two same-file line-range collisions that would
corrupt a pure-motion ticket, a set of dangling references, and a sequencing
document that undersold its own parallelism. The collisions and the DAG are fixed
in place (see Disposition); the remaining items are human decisions.

## Findings

Severity per `CLAUDE.md`: Critical = wrong distribution / broken build; High =
regression / API break; Medium = maintainability; Low = style.

| # | Finding | Severity | Status |
|---|---|---|---|
| F1 | `C2` and `C4` both claim `Context.cpp:21` (`kBoltzmann_kJ`). `C2` swept the whole anonymous namespace into `OutputWriter.cpp`; the constant is REX-only and would be unreachable by the REX TUs -> build break. | Critical | Fixed |
| F2 | `I1` and `R4` both extract `RobotIntegrator.hpp:241-258` (quaternion drift) to different targets, and the plan stated the order both ways. Under pure motion this double-moves the same lines. | High | Fixed |
| F3 | Line ranges are Stage-0-absolute and already drift (source is 1455/1334; tickets say 1456/1335). Within a lane each executed ticket invalidates downstream line numbers. | High | Mitigated (DAG mandates symbol re-anchoring) |
| F4 | `INV-10` referenced via `checkInv7AndInv10Guards` but undefined in `ARCHITECTURE section 5`. | Medium | Fixed (defined from the runtime guard) |
| F5 | `SPLIT-DEDUP-CTXFORWARDERS` unsequenced in README/MODULES; a coder running top-to-bottom never executes it. | Medium | Fixed (sequenced as a barrier) |
| F6 | Sequencing presented as a single serial chain; the real structure is four concurrent lanes + leaves + two barriers (~4x parallelism unused). | Medium | Fixed (DAG in README + MODULES section 4) |
| F7 | `SPLIT-W4` (NcmcMove) and `SPLIT-W9` (HmcMove) referenced as dependencies by W2/W5 but not yet written. | Medium | Fixed (both authored; ranges verified disjoint vs siblings; hostile review done - W9 clean, W4's one Blocking finding (missing `RobotIntegrator.hpp` for `stepTo<Bridge>` instantiation) + 2 nits applied) |
| F8 | `DOC-LegacyCoordSwapRex` targets `workflow/rex/LegacyCoordSwapRex.*`, which no SPLIT ticket creates. | Medium | Fixed (`SPLIT-C6` authored as the human-gated producer; DOC ticket + CALIBRATION marked conditional; C4/README/MODULES sequence it) |
| F9 | Dangling external-spec refs: `singular-dof-fixman.md` (R1), `reaction-force-monitoring.md` (W3), `gpu-cartesian-kinematics` (O6) - none exist in the repo. | Medium | Fixed (softened to "campaign, spec not yet in repo; behavior recovered from source"; the tickets restate the invariants in full regardless) |
| F10 | Run-type names diverge: `ARCHITECTURE section 3` and `DOC-SwapAcceptance` said "REMC/RENE/RJMC/REBAS"; the source enum `RUN_TYPE` (`include/ReplicaExchange.hpp:27`) is `REMC/RENEMC/RENE/REBASONTOP`. `RJMC`/`REBAS` do not exist. | Low | Fixed (aligned to the enum) |
| F11 | `C3` edits `runDrivenRound` (driven-REX), which OQ-5 says is never compiled/run and TESTS section 7 says throws - its bit-identity claim escapes the B1/B3 oracle. | Low | Fixed (oracle-coverage note added to `SPLIT-C3`: that site is verified by inspection, not by a passing test) |
| F12 | VERIFY section 5 states the contract-comment gate as universal, but TEST-001/004/005/006 defer it to TEST-007; B5 coverage is SHOULD in section 1 but a hard gate in the TEST tickets. | Low | Fixed (VERIFY section 5 scopes the contract-comment gate to TEST-007 and states B5 is a hard exit for the test refactor) |
| F13 | Inconsistent path/dash notation (`RobotEngine.cpp` vs `src/...`; en/em dashes) defeats the mechanical exit-criteria checks. | Low | Fixed for dashes (docs normalized to ASCII `-`); bare-filename vs `src/` path style left as-is (cosmetic, no machine check reads it) |
| F14 | README ticket inventory lagged the directory during live generation. | Low | Fixed (Phase B `TEST-002..008` enumerated; all SPLIT/DOC tickets confirmed referenced) |

DOC coverage, TEST<->VERIFY section 5 alignment, TEST<->TESTS seam mapping, INV/OQ
resolution (other than F4), and section-shape compliance were all verified clean.

## Parallelizability

The one-agent-per-ticket model is the intended design and works: each ticket
restates its invariants, names exact symbols, and ends in machine-checkable exit
criteria, so a fresh limited-context agent needs only the ticket plus the VERIFY
baseline artifacts - machine state, not conversational context.

Lane-level concurrency is available and now documented as a DAG (README section "Lane
structure", `MODULES.md section 4`): four serial lanes (World / Context / Bridge-OpenMM /
RobotEngine), two independent leaves, and two barriers (`O2`,
`DEDUP-CTXFORWARDERS`) that straddle Lane B and Lane C. Intra-lane parallelism is
NOT available - F3 means same-file tickets SHALL serialize and re-anchor by symbol.

The coordination state a limited-window agent cannot hold, and which therefore
belongs in a persistent orchestrator: the versioned baseline (`DEDUP-FORCEREDUCER`
adds a test to B1 mid-stream), the cross-lane test-file include-edit dedup, and
the per-lane worktree/merge discipline.

## Disposition (fixes applied this pass)

- `ARCHITECTURE.md section 5`: defined `INV-10` (drive/run-type pairing) from the
  `checkInv7AndInv10Guards` source. (F4)
- `SPLIT-C2.md`: narrowed the move to `boxFromReducedVectors` only; `kBoltzmann_kJ`
  and the namespace stay for C4. `SPLIT-C4.md`: reciprocal ownership note. (F1)
- `SPLIT-I1.md` / `SPLIT-R4.md`: decided as I1 -> R4; R4 re-anchors to
  `driftPositions`, no line overlap. (F2)
- `README.md`: added CX-9 (concrete I1/R4 order), CX-10 (C2/C4), the Lane
  structure, and the two barriers. `MODULES.md section 4`: rewrote sequencing as the
  same DAG. (F2, F5, F6)
- Authored `SPLIT-W4` (NcmcMove) and `SPLIT-W9` (HmcMove) - the last two World
  splits, ranges verified disjoint from all siblings. (F7)
- Authored `SPLIT-C6` (LegacyCoordSwapRex, human-gated) as the `runREX`-isolation
  producer; pointed `DOC-LegacyCoordSwapRex` at it, marked it conditional in
  `DOC-CALIBRATION`, sequenced it in Lane C (README/MODULES). (F8)

## Related deliverable - OpenMM test suite

`TEST-008-openmm-suite.md` (infra, outside the freeze): the vendored OpenMM test
suite is now gated behind `BUILD_OPENMM_TESTS` (OFF by default), runs all three
precisions (`single`/`mixed`/`double`) as parallel ctest cases with no GPU mutex,
registers the previously-built-but-unrun non-GPU tests, and is reachable via
`nox -s openmm_tests` (opt-in) while `nox -s tests` excludes them (`-LE openmm`).
`CMakeLists.txt` + `noxfile.py` edited; the `cuda-tests` configure with the flag
on validated clean.

## Human decisions - ACCEPTED (2026-07-13)

Both remaining decisions are accepted; no open human gate blocks Phase A.

1. The three campaign specs (singular-DOF Fixman, reaction-force-monitoring, GPU
   Cartesian-kinematics) are NOT authored as a prerequisite. The tickets proceed on
   their in-line invariant restatements, which are self-sufficient by charter; the
   campaign specs remain optional, deferrable documentation. (F9)
2. The CX-1...CX-8 decisions in `README.md` are accepted as recorded there:
   TU-split-first, `GLOB_RECURSE` + subtree layout, the O2 deletion and CTXFORWARDERS
   additions, the stretch-only BAT dedup, and the `engine_helpers.hpp` home with its
   one predicted Phase-A test edit.

3. OQ-2 (2026-07-13): investigating how to test `jointHDot_FM` found it already
   golden-tested by three independent, passing oracles (Simbody differential
   `TestRoboticsOracle` at 1e-8, the `JointJacobianDot` FD consistency check, and
   the `Integrator` energy-conservation tests). OQ-2 is CLOSED with no residual
   validation debt - the "gap" was a stale source caveat, now corrected. `SPLIT-R4`
   is mergeable. Recorded in `ARCHITECTURE.md` OQ-2.

Build status (2026-07-13): the split tree compiles clean under `cuda-release`
(`robo_bindings.so` links and installs). One cross-lane defect was fixed - the
I1-extracted free function `driftPositions` called `normalizeQuaternions`
unqualified; it now calls `RobotEngine::normalizeQuaternions`.

## Second review pass (2026-07-13)

A follow-up review checked the ticket set against the current source and normalized
language. Applied this pass:

- Corrected the run-type names to the source enum `RUN_TYPE`
  (`REMC/RENEMC/RENE/REBASONTOP`) in `ARCHITECTURE.md`, `DOC-SwapAcceptance`,
  `DOC-DrivenRexDriver`. (F10)
- Corrected the `TestNcmcExplicitSolvent` skip count (3 permanent stubs, 6/8 skip at
  default tier - not 7/8) in `TESTS.md`, `TEST-004`; fixed `runNcmcLoop` site count
  (seven, not five) in `TEST-004`.
- Removed the three non-existent external-spec citations (F9) and the W9/W4
  `RobotIntegrator.hpp` contradiction; added the F11 oracle-coverage note to
  `SPLIT-C3`, the BAT `>1e-9` guard note, the C5 banner-attribution note, and
  scoped the VERIFY contract-comment gate to `TEST-007` (F12).
- Normalized all docs to ASCII, replaced lowercase "must" with SHALL, and removed
  banned jargon (load-bearing, sidecar, fan out) and the "resolved" emphatic marker.
- The source audit found no ticket that introduces new functionality or a hidden
  behavioral change; every non-pure-motion ticket (O2 deletion, DEDUP-BATBENDSTRETCH,
  DEDUP-CTXFORWARDERS, C3, W5/W8 helper hoists, I1/R4) self-declares and is
  human-gated.
