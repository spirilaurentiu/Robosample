# Robosample Engine - Target Module Map

Status: recovery draft, 2026-07-12. Companion to [`ARCHITECTURE.md`](ARCHITECTURE.md).
This defines the target directory/file layout, the split proposals, and their
sequencing. It proposes no behavioral change: every split is pure code motion
(no logic edits, no renames beyond file placement). Full `SPLIT-###` /
`DOC-###` tickets are written only after human approval of this map.

Guiding rule (from `.claude/agents/architect.md`): one engineering concept per
header/source pair; split by responsibility count, never by line count.
Header target 100-300 LOC, `.cpp` 200-600 good / 700-1200 investigate / >1200
very likely multiple responsibilities.

---

## 1. Target directory layout

Illustrative, not final; names follow the recovered concepts. Layer order top
(highest) to bottom (lowest).

```
src/, include/robo/
  bindings/
    PyBind11.cpp                     # thin forwards only; prose moves to headers

  workflow/
    Context.{hpp,cpp}                # thin orchestrator: owns Worlds, delegates
    GibbsSweep.{hpp,cpp}             # the one "run a sweep over the schedule" primitive
    StartupValidator.{hpp,cpp}       # checkStartupGeometry
    OutputWriter.{hpp,cpp}           # moves/energy/DCD/reaction CSV + periodic imaging
    rex/
      ReplicaExchangeDriver.{hpp,cpp}   # label-swap RunREX + swap matrices + maps
      SwapAcceptance.{hpp,cpp}          # attemptREXSwap acceptance algebra (4 run-types)
      DrivenRexDriver.{hpp,cpp}         # RENE/REBASONTOP (behind a flag until tested; OQ-5)
      LegacyCoordSwapRex.{hpp,cpp}      # runREX oracle (retire-or-keep; OQ-3)
      Replica.hpp                       # Replica + ThermodynamicState value types

  world/
    World.{hpp,cpp}                  # residual: identity, thermo, schedule setters, composition
    ModelBuilder.{hpp,cpp}          # buildModel + root mobility + DSU + frame graph
    GeometryFitter.{hpp,cpp}        # setAtomsLocationsInGround / recomputeGeometry
    sampler/
      HmcMove.{hpp,cpp}             # reinitialize + metropolis + torsional/Cartesian core
      DockingMove.{hpp,cpp}         # repositionLigands + sphere sampling + findGoodStartingPose
      NcmcMove.{hpp,cpp}            # ncmcMove + inner GHMC + trough teleport
      VelocityDistortion.{hpp,cpp}  # NMA soft-mode + BAT drive momentum coupling
      CartesianSolvent.{hpp,cpp}    # solvent velocity draw / KE / save-restore
    FixmanCorrection.{hpp,cpp}      # calcFixman + logSineSqr + constraint log-det glue
    ReactionReporter.{hpp,cpp}      # per-body spatial-force snapshot + CSV

  dynamics/
    RobotEngine.{hpp,cpp}           # ABA recursions (kinematics/inertia/forward-dynamics)
    RobotEngine_massops.cpp         # mass-matrix operators (same class, cohesive TU split)
    RobotEngine_reaction.cpp        # reaction forces (same class, cohesive TU split)
    RobotIntegrator.{hpp,cpp}       # verletStep/stepTo/checkReversibility (templated on Bridge)
    JointKernels.{hpp,cpp}          # ALL per-joint switches: X_FM, H_FM, HDot_FM, QDot, QDotDot, drift
    Constraints.{hpp,cpp}           # SHAKE/RATTLE/loop-Fixman with shared G-assembly
    BatScaling.{hpp,cpp}            # applyBend + applyStretch + Jacobian
    NMA.hpp                         # Route-B normal-mode analysis
    NCMCProtocol.hpp                # protocolLambda

  bridge/
    ForceBridge.{hpp,cpp}           # robo:: <-> OpenMM adapter; owns nothing OpenMM-typed in .hpp
    ForceReducer.{hpp,cpp}          # THE single force->wrench reduction (host); CUDA kernel validated against it
    OpenMMContext.{hpp,cpp}         # residual: holds System/Context/Integrator, forwards
    OpenMMSystemBuilder.{hpp,cpp}   # initialize() force/particle/box/integrator construction
    ForceFactory.{hpp,cpp}          # create*Force builders + addStandardExclusions
    AlchemyForceFactory.{hpp,cpp}   # alchemy/decoupling forces + enable/setLambda
    GpuKinematics.{hpp,cpp}         # gGpuKin + fused nvrtc pipeline (CUDA-only)
    MTSIntegrator.{hpp,cpp}         # r-RESPA substep integrator

  model/
    RobotModel.hpp                  # immutable topology + joint-fact predicates
    RobotState.hpp                  # per-step mutable cache (stage-validity contract; OQ-1)
    TopologyElements.hpp            # SystemTopology (later: split into sub-structs; OQ deferred)
    PeriodicBox.hpp                 # box value type + reducedBoxVectors (single source)

  math/
    robot_vecmat.hpp                # Vec/Mat/Rotation/Transform/Quat
    robot_spatial.hpp               # SpatialVec/SpatialInertia/ArticulatedInertia/PhiMatrix
    hinge_linalg.{hpp,cpp}          # dense n<=6 solver: invertDense/jacobiSymEig/symSqrt/pseudoLogDet
    robo_debug.hpp                  # ROBO_DEBUG NaN scanner (guarded)

  io/
    DCDWriter.{hpp,cpp}             # already clean; move as-is
  util/
    MemoryArena.hpp  Units.hpp
```

`robot_math.hpp` splitting (into `robot_vecmat` + `robot_spatial`) is optional
and compile-cost-driven; sequence it last or defer. `SystemTopology` splitting
is a larger change gated behind a decision record (deferred).

---

## 2. Split proposals, by source file

Priority reflects debt x independence (self-contained, already-gated concerns
are safest to move first). Sizes are rough targets.

### `World.cpp` (2619) - highest debt, 9 extractions

| # | Extract | From (approx) | Independence | Target |
|---|---|---|---|---|
| W1 | `ModelBuilder` | `buildModel` 601-968, root mobility, DSU, frame graph | high (build-time only) | `world/ModelBuilder` |
| W2 | `DockingMove` | `repositionLigands`, sphere sampling, `findGoodStartingPose` | high (gated by `docking_`) | `world/sampler/DockingMove` |
| W3 | `ReactionReporter` | `captureReactionSnapshot` + 5 `reaction*_` members | high (opt-in, orthogonal) | `world/ReactionReporter` |
| W4 | `NcmcMove` | `ncmcMove`, `ncmcInnerGhmcStep`, trough teleport | medium | `world/sampler/NcmcMove` |
| W5 | `FixmanCorrection` | `calcFixman`, `calcLogSineSqrGamma2`, `lnDetMCartesian_` | high (pure functions) | `world/FixmanCorrection` |
| W6 | `GeometryFitter` | `setAtomsLocationsInGround`, `recomputeGeometry` | medium | `world/GeometryFitter` |
| W7 | `CartesianSolvent` | `setCartesianSolvent`, `drawSolventVelocities`, `calcSolventKE` | medium | `world/sampler/CartesianSolvent` |
| W8 | `VelocityDistortion` | NMA soft-mode + BAT drive momentum coupling | medium | `world/sampler/VelocityDistortion` |
| W9 | `HmcMove` core | `reinitialize`, `metropolis`, `generateSample` branches | low (touches all above) | `world/sampler/HmcMove` |

`generateSample` (272 lines) collapses into a thin dispatch once W2/W4/W9 land;
the three move regimes become a `Move` strategy. W9 is sequenced **last** among
World splits because it depends on the extracted services.

### `Context.cpp` (1334) - REX entanglement, 5 extractions

| # | Extract | From (approx) | Independence | Target |
|---|---|---|---|---|
| C1 | `StartupValidator` | `checkStartupGeometry` 237-395 | high (pure) | `workflow/StartupValidator` |
| C2 | `OutputWriter` | `writeOutputs*`, `writeReactionRows`, DCD scratch, imaging | high | `workflow/OutputWriter` |
| C3 | `GibbsSweep` | the sweep loop duplicated at 457-490 / 977-1021 / 1112-1143 | high (dedup) | `workflow/GibbsSweep` |
| C4 | `ReplicaExchangeDriver` + `SwapAcceptance` | `RunREX`, `attemptREXSwap`, maps, matrices | medium | `workflow/rex/*` |
| C5 | `DrivenRexDriver` | `runDrivenRound`, `driveReplica`, interleave | medium (untested; OQ-5) | `workflow/rex/DrivenRexDriver` |

`runREX` legacy path (C6) is retire-or-keep per OQ-3. OQ-3 is DECIDED "keep as
oracle", so C6 (`SPLIT-C6.md`) is written: it isolates `runREX` into
`workflow/rex/LegacyCoordSwapRex.{hpp,cpp}` after C4, human-gated. Retirement is
the alternative that drops C6.

### `OpenMMContext.cpp` (1165) - 6 extractions

| # | Extract | Independence | Target |
|---|---|---|---|
| O1 | `MTSIntegrator` | high | `bridge/MTSIntegrator` |
| O2 | `PeriodicBox` box math (fold duplicate `computePeriodicBoxVectors_Context`) | high | `model/PeriodicBox` |
| O3 | `ForceFactory` (create*Force + exclusions) | medium | `bridge/ForceFactory` |
| O4 | `AlchemyForceFactory` | medium (NCMC-only) | `bridge/AlchemyForceFactory` |
| O5 | `OpenMMSystemBuilder` (initialize construction body) | medium | `bridge/OpenMMSystemBuilder` |
| O6 | `GpuKinematics` (gGpuKin + fused pipeline) | high (CUDA-only unit) | `bridge/GpuKinematics` |

### `RobotEngine.cpp` (1455) - extract non-ABA concerns

| # | Extract | Independence | Target |
|---|---|---|---|
| R1 | `hinge_linalg` (invertDense/jacobiSymEig/symSqrt/pseudoLogDet, ~450 L) | high (pure numeric) | `math/hinge_linalg` |
| R2 | `robo_debug` NaN scanner (~150 L, `ROBO_DEBUG`) | high | `math/robo_debug` |
| R3 | cohesive TU split of the class (kinematics / massops / reaction) | mechanical | `dynamics/RobotEngine_*.cpp` |
| R4 | consolidate joint switches (`calcQDot`/`calcQDotDot`/drift) into `JointKernels` | medium (touches taxonomy; OQ-2) | `dynamics/JointKernels` |

### `RobotIntegrator.hpp` (638) - split the god-template

| # | Extract | Target |
|---|---|---|
| I1 | `driftPositions`, `velocityCorrector`, `cartesianSolventVerlet` out of `verletStep` (350 L) | privates in `RobotIntegrator` |

### Smaller, self-contained

| Extract | From | Target |
|---|---|---|
| `ForceReducer` (dedup host vs CUDA reduction; INV-1) | `ForceBridge.hpp` + `OpenMMContext.cpp` kernel | `bridge/ForceReducer` |
| `assembleConstraintRow` (dedup 3 G-assembly loops) | `Constraints.cpp` | private in `Constraints` |
| `applyBend`/`applyStretch` + route through `readBatCoord` | `BatScaling.cpp` | privates in `BatScaling` |
| `flatToVec3`/`vec3ToFlat` marshalling | `PyBind11.cpp` + `World` | `util` or `World` overload |
| `Context::setSeparateForceGroups`/`setEnforcePeriodicBox` forwarders | `PyBind11.cpp` singleton reach-through | `Context` |

---

## 3. Dependency-defect resolutions

Each maps to a defect in `ARCHITECTURE.md section 6`. These change include structure,
not behavior; where an API narrows they are flagged for separate human approval,
never bundled into a code-motion split.

1. **Solver->OpenMM inversion (section 6.1):** introduce an abstract force-bridge interface
   (or forward declaration) so `RobotEngine.hpp` no longer transitively includes
   `OpenMMContext.hpp`. `ForceReducer` becomes the one place the reduction lives.
2. **RobotEngine<->Constraints cycle (section 6.2):** keep the forward-declaration
   discipline; codify it with an include-cycle CI check so a future `.hpp`-level
   use fails loudly.
3. **Bindings->singleton (section 6.3) and World hidden global (section 6.4):** add thin
   forwarders on `Context`/`World`; bindings call the object, not `get()`.
4. **Duplicated reduction (section 6.5):** addressed by `ForceReducer`; the CUDA kernel is
   validated against it in a test.
5. **Scattered joint taxonomy (section 6.6):** addressed by R4 (all joint switches in
   `JointKernels`).

---

## 4. Sequencing

The merge hazard is same-file editing, not layer order. Tickets on the same
god-file serialize; tickets on disjoint files run concurrently. The sequencing is
therefore a **DAG of four serial lanes plus leaves and two barriers**, not a
single chain. A flat total order (leaf-first) remains valid - it is any
topological sort of this DAG - but it forgoes the ~4x parallelism the lanes give.

**Lanes (serial within, concurrent across).** Within a lane each ticket shifts
its file's line numbers, so downstream tickets **re-anchor by symbol**; the
line ranges in every ticket are Stage-0 advisory, the symbol name authoritative.

- **Lane E - RobotEngine** (`RobotEngine.cpp` / `.hpp`, `RobotIntegrator.hpp`):
  R1 -> R2 -> R3 -> I1 -> R4. R4 depends on I1 (both edit the quaternion-drift block;
  I1 extracts it into `driftPositions`, R4 lifts it into `JointKernels` - see
  README CX-9); R4 also gated on OQ-2.
- **Lane B - Bridge/OpenMM** (`OpenMMContext.cpp` / `.hpp`, `ForceBridge.hpp`):
  DEDUP-FORCEREDUCER -> O1 -> O6 -> O3 -> O4 -> O5 (SystemBuilder last).
- **Lane W - World** (`World.cpp` / `.hpp`): W1 -> W3 -> W5 -> W2 -> W6 -> W7 -> W8 ->
  W4 -> W9 (HmcMove last). Longest lane, the critical path. W5 precedes W2/W4
  (shared `engine_helpers.hpp`, CX-6).
- **Lane C - Context** (`Context.cpp` / `.hpp`): C1 -> C2 -> C3 -> C4 -> C5. C2 and
  C4 both touch the top-of-file anonymous namespace: C2 owns `boxFromReducedVectors`,
  C4 owns `kBoltzmann_kJ` (README CX-10).
- **Leaves (independent, any time):** DEDUP-CONSTRAINTROW (`Constraints`),
  DEDUP-BATBENDSTRETCH (`BatScaling`), optional `robot_math` split.

**Barriers - edit two god-files at once; block the lanes they straddle.** Run
before both lanes begin, or after both finish, never concurrently:

- O2 (`PeriodicBox` fold): `OpenMMContext.cpp` + `Context.cpp` - Lane B and Lane C.
- DEDUP-CTXFORWARDERS (section 6.3/section 6.4): `Context.hpp` + `OpenMMContext.hpp` +
  `Context.cpp` - Lane B and Lane C; adds two public symbols, so human-approved
  (VERIFY section 3), not pure motion. Recommended slot: after Lane C completes.

Documentation (`DOC-###`) executes **after** the splits so Doxygen lands on final
files. Test motion is a separate phase after all source splits are verified (the
two-phase freeze rule); see [`TESTS.md`](TESTS.md).

---

## 5. What does NOT move

- `DCDWriter`, `Units`, `MemoryArena`, `NCMCProtocol`, `NMA` - already single-concept.
  They get doc tickets and a directory home, no split.
- `RobotModel` / `RobotState` - stay one concept each; `RobotState` gets the
  stage-validity documentation (OQ-1), not a split.
- `SystemTopology` sub-struct split - deferred behind a decision record; too broad
  to bundle now.
- Any public API change - flagged separately for human approval, never inside a
  code-motion ticket. Target: zero exported-symbol delta on public symbols.
