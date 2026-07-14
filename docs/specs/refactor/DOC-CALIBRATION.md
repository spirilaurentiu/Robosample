# DOC-CALIBRATION: style exemplars and the calibration-first rule

Status: draft, 2026-07-12. This ticket runs **first** in the doc phase, before any
other `DOC-###`. Executor: the `documenter` agent (`.claude/agents/documenter.md`
section 13); the documenter loads `styles/reference.md` as it writes. It is a
blocking prerequisite for the rest of the doc phase.

The doc phase converts the existing dense header prose into structured Doxygen and
separates contract from derivation/history (ARCHITECTURE section 8). Scope, per the user's
recorded decision, is **contract-separation on SPLIT-touched files only** - the
modules created or touched by the `SPLIT-###` tickets (MODULES.md section 1: `world/`,
`dynamics/`, `bridge/`, `workflow/`, `model/`, `math/`, `io/`, `util/`), not a
full-repo conversion.

## Why calibrate first

The 40-plus module tickets share one house style, but that style is not yet agreed.
Documenting three representative files first, stopping for human review, and then
attaching the approved versions to every later ticket keeps the whole set coherent
and prevents 40 tickets of drift. This is human gate #3 in the Architect pipeline
(`documenter.md` section 0, section 13). No downstream `DOC-###` ticket starts until the trio
below is approved.

**Rule.** After the trio is approved, its files stay in the documenter's context as
the style exemplars. Each later ticket carries the line "attach the approved
calibration exemplars"; the documenter matches their form. When an exemplar and
`documenter.md` conflict, the spec wins and the conflict is a finding
(`documenter.md` section 8).

## The calibration trio

One representative file per domain (`documenter.md` section 13: one `.cpp`, one CUDA unit,
one test file). The three are chosen to cohere around one invariant (INV-1,
force->wrench reduction) so the exemplars demonstrate contract-vs-mechanism on the
same physics from three angles.

1. **Representative `.cpp` - `bridge/ForceReducer.cpp`** (extracted by
   `SPLIT-DEDUP-FORCEREDUCER`). The single host force->wrench reduction. Small,
   pure, one clear contract (INV-1), real ownership boundary (borrows OpenMM
   per-atom forces, writes `RobotState::bodyForceG`), and a virtual-site skip
   (INV-2) that is pure contract. It exercises: ownership verbs, effect statements,
   an invariant stated as a `@post`, and the "skip massless particles" rule that
   SHALL read as behavior, not as a loop.

2. **CUDA unit - `bridge/GpuKinematics.cpp`** (extracted by `SPLIT-O6`). The
   `gGpuKin` fused nvrtc pipeline: two JIT kernels that push body transforms into
   OpenMM `posq` and reduce forces on-device (`documenter.md` section 5). The engine ships
   no `.cu` file - the CUDA path is nvrtc source strings compiled at runtime, so
   the launch contract (decomposition, launch config, shared memory, stream,
   completion, INV-1 device-side reduction that SHALL match the host reducer) is
   recovered from the host launch sites and the kernel index math together. This is
   the exemplar for the CUDA rules in `documenter.md` section 5; it also anchors the
   host<->device INV-1 parity claim that the ForceReducer exemplar states on the host
   side.

3. **Test file - `tests/TestReactionForces.cpp`** (FAST, contract; TESTS.md section 3).
   A contract test for INV-1 driven through the `AnalyticForceBridge`, so it shows
   the standard gtest form: `@file` block naming module and invariant, fixture-free
   build-in-body prologue (TESTS.md section 4), and one Given/When/Then contract line per
   `TEST` referencing the INV it defends (`documenter.md` section 7). It is the test-side
   view of the same INV-1 the two source exemplars document.

## Deliverables

- Fully documented `bridge/ForceReducer.cpp` (+ its header declarations).
- Fully documented `bridge/GpuKinematics.cpp` (+ header; nvrtc kernel launch
  contracts per `documenter.md` section 5).
- Fully documented `tests/TestReactionForces.cpp` (`documenter.md` section 7 form).
- `findings/DOC-CALIBRATION-findings.md` (present even if empty).
- Any Doxyfile CUDA additions required to build the `GpuKinematics` docs
  warning-free (`documenter.md` section 9) - the one permitted non-comment edit, reported
  in findings.

## Exit criteria

- The three files build **warning-free** under Doxygen (warnings-as-errors).
- Comment-stripped before/after diff of each of the three files is **empty**
  (proves zero code change; VERIFY section 4).
- `findings/DOC-CALIBRATION-findings.md` present; every `@note Assumed:` in the
  diff has a matching findings entry.
- **Human review recorded.** The approved three files become the attached exemplars
  for every subsequent `DOC-###`. Until this approval is recorded, the doc phase is
  blocked.

## Execution order for the module tickets (after approval)

Leaf helpers before public API (`documenter.md` section 0), lower layers before higher, so
a helper's recovered contract is in context when its callers are documented:

1. **Infrastructure / math (leaves):** DOC-RobotMath, DOC-hinge_linalg,
   DOC-robo_debug, DOC-MemoryArena, DOC-Units, DOC-DCDWriter.
2. **Domain data (`model/`):** DOC-PeriodicBox, DOC-TopologyElements,
   DOC-RobotModel, DOC-RobotState.
3. **Bridge:** DOC-ForceReducer\*, DOC-ForceBridge, DOC-GpuKinematics\*,
   DOC-ForceFactory, DOC-AlchemyForceFactory, DOC-MTSIntegrator,
   DOC-OpenMMSystemBuilder, DOC-OpenMMContext.
4. **Algorithms (`dynamics/`):** DOC-JointKernels, DOC-Constraints, DOC-BatScaling,
   DOC-NMA, DOC-NCMCProtocol, DOC-RobotEngine, DOC-RobotIntegrator.
5. **Domain (`world/`):** DOC-ModelBuilder, DOC-GeometryFitter,
   DOC-FixmanCorrection, DOC-ReactionReporter, DOC-CartesianSolvent,
   DOC-VelocityDistortion, DOC-DockingMove, DOC-NcmcMove, DOC-HmcMove, DOC-World.
6. **Workflow:** DOC-StartupValidator, DOC-OutputWriter, DOC-GibbsSweep,
   DOC-Replica, DOC-SwapAcceptance, DOC-ReplicaExchangeDriver, DOC-DrivenRexDriver,
   DOC-LegacyCoordSwapRex (conditional - runs only if SPLIT-C6 ran; see the ticket),
   DOC-Context.

\* Already documented as calibration exemplars; their tickets record that and cover
only residual symbols not touched during calibration.
