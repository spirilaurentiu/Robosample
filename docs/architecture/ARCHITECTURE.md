# Robosample Engine - Recovered Architecture

Status: recovery draft, 2026-07-12. Scope: the C++/CUDA engine under `src/` and
`include/` only. The nested `Robosample/` directory (the original SimTK/Simbody
sources) and the vendored `openmm/`, `pybind11/`, `pcg-cpp/`, `units/`
submodules are out of scope.

This document recovers the design that already exists in the code. It assumes
the implementation is correct and describes engineering structure only -
ownership, lifetimes, dependencies, execution order, invariants. It proposes no
behavioral change. Module boundaries and split proposals live in
[`MODULES.md`](MODULES.md); the test-suite triage lives in
[`TESTS.md`](TESTS.md).

---

## 1. What the engine is

Robosample samples molecular conformations with blocked Gibbs sampling coupled
to Hamiltonian Monte Carlo. Each molecule is an articulated robot - a kinematic
tree of rigid bodies joined by bond/angle/torsion joints. A *World* exposes one
robot factorization (one choice of active generalized coordinates); a *Context*
owns several Worlds and sweeps them as Gibbs blocks, optionally under replica
exchange. OpenMM supplies per-atom forces and energies; a custom SimTK-free
articulated-body solver turns those into joint-space dynamics for the HMC
integrator.

The engine is ~8.2k lines of `.cpp` and ~5.6k lines of headers across 20
translation units. Four files hold most of the mass and most of the debt:
`World.cpp` (2619), `RobotEngine.cpp` (1455), `Context.cpp` (1334),
`OpenMMContext.cpp` (1165).

**Replica exchange is the intended center (design directive, 2026-07-12).**
Robosample SHALL be organized around replica exchange: a run is `R` replicas on
a temperature ladder, and a single-replica run is the degenerate `R = 1` case,
not the norm. The current *runnable* production path does not reflect this - every
driver script (`run.py`, `roborun.py`, all `run_*.py`) calls the legacy
coordinate-swap `run_rex` at `NOF_REPLICAS = 1`, so what actually ships is
degenerate single-replica sampling. The label-swap driver `run_rex_label_swap`
(`Context::RunREX`) is the real multi-replica REX and is validated: a 4-replica
REMC ladder runs end to end on CUDA (RTX 3090, CUDA 13.0) and produces the correct
alternating-parity neighbour swap matrix.

Scope split. Making `RunREX` REMC the default driver and modelling `R = 1` as its
special case is a behavioral change; it does NOT ride inside the code-motion
refactor. It is tracked separately in `docs/specs/rex-default-driver.md` and runs on
its own (README "Behavioral specs"). What the behavior-preserving split takes from
this directive is only a layout consequence: Context's REX code is treated as core,
not an optional branch, which reshapes where the `Context` split places it
(`MODULES.md section 2`). The split itself moves REX code verbatim and changes no
driver default.

---

## 2. Layer map

Dependencies point downward only. Two upward/lateral edges violate this and are
recorded as defects in section 6.

```
Application    PyBind11.cpp ................ Python module `robo_bindings`
                     |
Workflow       Context / ReplicaExchange ... owns N Worlds, Gibbs sweep, REX, I/O
                     |
Domain         World ...................... one robot factorization + its sampler
                     |
Algorithms     RobotEngine (ABA) * RobotIntegrator (HMC) * Constraints
               * BatScaling * NMA * JointKernels * ForceBridge * NCMCProtocol
                     |
Bridge         OpenMMContext .............. the one OpenMM System/Context/Integrator
                     |
Domain data    RobotModel * RobotState * TopologyElements * PeriodicBox
                     |
Infrastructure robot_math * MemoryArena * Units * DCDWriter
                     |
Platform       OpenMM * CUDA/nvrtc * LAPACK/BLAS * nholthaus units
```

The cleanest boundary in the system is the force bridge: `RobotEngine`'s
integrator is templated on a `Bridge` type (`RobotIntegrator.hpp`), so OpenMM
enters only the `World` translation unit and the tests can substitute an
analytic force bridge. The strongest data-structure property is the
model/state/algorithm split in the solver (section 4).

---

## 3. How a run actually executes

**Construction.** Python (`context.py`) parses AMBER inputs, fills the
`SystemTopology` SoA buffer, and calls `Context(base_name, seed)`. `Context`
adds Worlds (`addCartesianWorld` / `addRoboticWorld` / `addDockingWorld`), each
of which calls `World::buildModel` to turn bonds + per-joint mobilities into a
`RobotModel` (body forest via union-find + BFS, q/u index tables, static joint
frames, mass properties, loop-closure constraints). `Context::initialize`
brings up the single `OpenMMContext` (particles, forces, integrator, platform),
runs `checkStartupGeometry` (an O(N^2) clash/PE sanity scan), and truncates the
output files.

**Inter-world currency.** The only state passed between Worlds is per-atom
Ground-frame coordinates in nm. Worlds hold no persistent per-replica state
(INV-3). At each sweep position `Context` pushes coordinates into the next World
(`setAtomsLocationsInGround`, which re-fits every rigid-body frame and the
`X_PF`/`X_BM` transforms), runs one sampling move, and reads coordinates back
(`getAtomsLocationsInGround`).

**One HMC move** (`World::generateSample` -> `reinitialize` -> the integrator):
1. `reinitialize` seeds momenta from `sqrt(M)` (via `RobotEngine::multiplyBySqrtM`),
   optionally distorts them (NMA soft-mode / BAT drive), draws Cartesian-solvent
   velocities, enforces velocity constraints (RATTLE), evaluates forces, and
   assembles PE/KE/Fixman/Jacobian into the starting Hamiltonian.
2. The integrator (`RobotEngine::stepTo` -> `verletStep`) runs leapfrog in
   generalized coordinates: `realizePosition` -> `realizeVelocity` ->
   `realizeArticulatedBodyInertias` (the ABA inertia recursion) -> `calcUDot`
   (forward dynamics), position drift with quaternion exp-map, SHAKE position
   projection, implicit-trapezoid velocity correction.
3. Forces come from the `ForceBridge`: either the host path (OpenMM returns
   per-atom forces, ForceBridge reduces them to per-body spatial wrenches) or
   the fused CUDA path (two nvrtc kernels push body transforms into OpenMM's
   `posq` and reduce forces on-device).
4. `metropolis` accepts or rejects on the total-energy change; on reject the
   saved `q`/positions are restored.

**A Gibbs sweep** is the loop `for each World in schedule: push coords ->
generateSample -> pull coords`. This exact loop is written three times in
`Context.cpp` (plain REX, driven RENE round, REMC subround).

**Replica exchange** wraps sweeps. Two drivers coexist: `runREX` (legacy
coordinate-swap, retained only as a correctness oracle) and `RunREX`
(label-swap with `Replica`/`ThermodynamicState` and two inverse index maps).
Swap acceptance covers four run-types (REMC / RENEMC / RENE / REBASONTOP) in
`attemptREXSwap`. Output (moves CSV, per-replica energy CSV, per-replica DCD
with whole-molecule periodic imaging, reaction-force CSV) is emitted by
`writeOutputs`.

---

## 4. Ownership and lifetimes

| Object | Owner | Lifetime | Notes |
|---|---|---|---|
| `SystemTopology` | `Context` (public member) | whole run | filled in-place from Python; the Python<->engine payload |
| `World` | `Context` (`vector<unique_ptr<World>>`) | whole run | heap so `World&` stays valid across sweeps; one per factorization |
| `RobotModel` | `World` (by value) | rebuilt on root-mobility change | immutable topology after `buildModel` |
| `RobotState` | `World` (by value) | per-World | one `MemoryArena` slab; hands out raw cache pointers |
| `ForceBridge` | `World` (by value) | per-World | observer of `RobotModel`; adapter to the OpenMM singleton |
| `OpenMMContext` | process singleton (`get()`) | whole run | the one `System`/`Context`/`Integrator`; `gGpuKin` file-static holds CUDA state |
| `ConstraintSet` | `World` (by value) | per-World | loop-closure distance constraints |
| RNG | `World` (`rng_`) and `Context` (`rexRng_`) | respective owner | seeded from the run seed |

`RobotEngine`, `JointKernels`, `NMA::RouteBNMA`, and `ncmc::protocolLambda` own
no state: they are static/free functions over `(const RobotModel&,
RobotState&, ...)`. The engine "class" is a namespace of ~20 static methods, a
deliberate structure-of-arrays driver style rather than object orientation.

---

## 5. Invariants (the "SHALL NOT change" list)

These are recovered from source and corroborated by tests. They feed every
split/doc ticket's "SHALL remain true" section. Items marked *(source-only)*
lack a corroborating test and are hypotheses until one is written.

- **INV-1 Force->wrench convention.** Per-body spatial force in `RobotState::bodyForceG`
  is `(torque about the body origin, net force)` in the Ground frame. Both the
  host reduction (`ForceBridge.hpp`) and the CUDA `reduceForces` kernel SHALL
  produce this identically.
- **INV-2 Virtual sites.** Virtual-site forces are already projected onto parent
  atoms by OpenMM; the reduction skips massless particles to avoid double
  counting.
- **INV-3 Stateless Worlds.** A World carries no persistent per-replica state
  between sweeps; per-atom Ground coordinates (nm) are the sole inter-world
  currency.
- **INV-4 Realization order.** `RobotState` caches are valid only in stage order
  position -> velocity -> articulated-body inertias -> udot. The contract is
  enforced by convention and doc comments, not by types. *(source-only; highest-risk API)*
- **INV-5 Mass metric.** HMC momentum seeding and kinetic energy use the same
  mass matrix operators (`multiplyBySqrtM` / `multiplyByMInv` / `calcKineticEnergy`);
  `calcLogDetM` is the Fixman kinetic term.
- **INV-6 Constraint consistency.** The `G` (constraint Jacobian) assembly is
  identical across SHAKE, RATTLE, and the loop-closure Fixman log-det, so the
  correction matches the projection. `calcConstraintLogDet` returns 0 for
  acyclic molecules.
- **INV-7 BAT map/Jacobian agreement.** `applyBatScaling` and its Cartesian
  log-Jacobian read identical `(r,theta)` geometry; one anchor snapshot per round is
  shared by both swap partners so the paired map is an exact involution (INV-9
  in the BAT spec).
- **INV-8 REX detailed balance.** `runREX` and `RunREX` produce equivalent
  sampling; `runREX` is the differential oracle for `RunREX`. Swap acceptance
  guards reject on NaN/inf.
- **INV-9 Quaternion double cover.** Free/quaternion joints normalize
  quaternions each step; reversibility checks account for the double cover.
- **INV-10 Drive/run-type pairing.** A driven run type pairs with exactly one
  distortion regime: `RENE`/`REBASONTOP` drive a volume-changing BAT-scaling
  world (`distortOption == ScaleBendStretch`, carrying the `lnJac` term);
  `RENEMC` drives a volume-preserving velocity/NMA world (`distortOption == NMA`,
  omitting `lnJac`). At least one driven world SHALL be configured, or the drive
  is silently inert. A mismatched drive biases swap acceptance (INV-1/B9). This
  invariant is orthogonal to the INV-7 Fixman-in-sampler biconditional (the
  source labels that check `INV-7/V9`). Enforced at runtime by
  `Context::checkInv7AndInv10Guards`, which throws before the first round.
  *(source-only; runtime-guarded, no dedicated regression test)*

---

## 6. Dependency defects (recorded, not preserved)

1. **Solver -> OpenMM inversion.** `ForceBridge.hpp` includes `OpenMMContext.hpp`,
   and `RobotEngine.hpp` includes `ForceBridge.hpp`. So the articulated-body
   solver header transitively depends on the OpenMM bridge - an upward edge from
   Algorithms into the Bridge layer. The integrator sidesteps this by templating
   on `Bridge` in a separate header (`RobotIntegrator.hpp`), but the base engine
   header still carries the include. Resolution: forward-declare the bridge
   interface; keep `ForceBridge` as the only place `robo::` meets OpenMM.
2. **Latent RobotEngine <-> Constraints cycle.** `RobotEngine.hpp` forward-declares
   `robo::ConstraintSet`; `RobotEngine.cpp` and `RobotIntegrator.hpp` include
   `Constraints.hpp`. The cycle is managed by discipline (leaf-header rule) and
   re-triggers if any `.hpp`-level use of `ConstraintSet` is added to
   `RobotEngine.hpp`.
3. **Bindings reach the OpenMM singleton.** Two `Context` bindings in
   `PyBind11.cpp` (`set_separate_force_groups`, `set_enforce_periodic_box`) call
   `OpenMMContext::get()` directly and discard their `Context&` argument.
4. **Hidden global in `World`.** `World::nDof()` reaches
   `OpenMMContext::get().getNumDegreesOfFreedom()` - a global coupling not
   visible in the signature.
5. **Force reduction duplicated.** The wrench reduction exists twice (host loop in
   `ForceBridge.hpp`, CUDA kernel in `OpenMMContext.cpp`) with no shared
   function - a divergence risk against INV-1 and not unit-testable in isolation.
6. **Joint taxonomy scattered.** Per-joint-type `switch` statements appear at four
   sites (`JointKernels`, `calcQDot`, `calcQDotDot`, the quaternion-drift branch
   in `verletStep`) instead of being centralized in `JointKernels`. Adding a
   joint touches four files.

---

## 7. Per-module responsibility (recovered)

- **`World`** - *god-object, ~12 responsibilities.* Kinematic-tree construction,
  geometry re-fitting, root-mobility mutation, OpenMM force bridging, torsional
  GC-HMC, Cartesian MD, Fixman/Jacobian corrections, docking moves, NCMC
  alchemical switching, Cartesian-solvent sub-integration, NMA/BAT velocity
  distortion, and reaction-force reporting - plus RNG, thermodynamics, mass
  preconditioning, and telemetry. Extraction seams in [`MODULES.md`](MODULES.md).
- **`Context`** - *god-object.* World lifecycle, OpenMM bring-up, startup
  validation, Gibbs scheduling (inlined three times), two REX drivers, output
  I/O, BAT anchor accumulation. NCMC and NMA are *not* in `Context` - they are
  already modular. REX is ~700 of the 1334 lines and is entangled with itself.
- **`RobotEngine` + `RobotIntegrator`** - one logical class split across a `.cpp`
  and a header (the header carries the two ~350/133-line integrator templates).
  Clean data/algorithm split, but 450 lines of dense `n<=6` linear algebra and
  150 lines of a `ROBO_DEBUG` NaN scanner are buried in the dynamics `.cpp`.
- **`RobotModel` / `RobotState`** - immutable topology and per-step mutable cache.
  `RobotState` exposes ~35 raw cache pointers with only convention-enforced
  stage validity (INV-4) and carries NCMC-solvent fields (concern bleed).
- **`robot_math`** - 867-line header of spatial-algebra value PODs (`Vec3`,
  `Rotation`, `Transform`, `SpatialVec`, `SpatialInertia`, `ArticulatedInertia`,
  `PhiMatrix`, `Quat`). Breadth, not depth; internal consolidation is good.
- **`OpenMMContext`** - the OpenMM adapter singleton. `initialize` (~235 lines)
  mixes particle/box/force/integrator/platform construction. Force builders,
  fused-CUDA kinematics pipeline (`gGpuKin`), and the alchemy force factory are
  separable.
- **`ForceBridge`** - the sole adapter between `robo::` types and OpenMM; holds
  the host force->wrench reduction.
- **`Constraints`** - loop-closure SHAKE/RATTLE + loop Fixman log-det; three
  copies of the `G`-assembly loop.
- **`BatScaling`** - deterministic BAT position-scaling drive for RENE, with a
  112-line `applyBatScaling` god-function (bend + stretch + Jacobian).
- **`NMA`** - header-only Route-B normal-mode analysis over `RobotModel`/`RobotState`.
- **`NCMCProtocol`** - a single `protocolLambda` free function.
- **`TopologyElements`** - the `SystemTopology` SoA god-struct (~110 fields).
- **`DCDWriter`** - self-contained CHARMM DCD writer, zero engine dependency; the
  cleanest module in the set.
- **`PeriodicBox` / `MemoryArena` / `Units`** - small infrastructure value types.
- **`PyBind11`** - the Python binding surface; mostly thin, with reference prose
  mislocated into docstrings and two singleton reach-throughs.

---

## 8. Documentation state

Coverage is *abundant but non-uniform and non-Doxygen*. Headers carry dense
prose contracts (often 10-40 lines, citing spec IDs and bug history); `.cpp`
files carry zero Doxygen blocks. The refactoring risk is not missing docs but
that contract, derivation, and changelog are interleaved in the same comments,
and that transient status ("NOT compiled or run", Stage-2b/2c TODOs) and
spec-ID cross-references will rot. Documentation work is therefore *conversion
and separation* (contract -> structured Doxygen; derivation/history -> module
notes or specs), not authorship from scratch. The genuine gaps:

- `RobotState`'s single-letter cache accessors (`P()`, `Z()`, `G()`, ...) - no
  stated stage-validity contract (INV-4).
- The energy/force evaluation entry points in `OpenMMContext` (the actual
  runtime surface) are the least documented; config toggles are the most.
- `jointHDot_FM` for BendStretch/SphericalCoords/FreeLine is documented as
  "faithful on paper but NOT golden-tested" - an unvalidated guarantee, not a
  doc gap, and a first-class OPEN-QUESTION.

---

## 9. OPEN-QUESTIONS

Recorded for a human, not decided by fiat.

- **OQ-1** INV-4 (realization stage order) is enforced only by convention. Should
  the refactor introduce a type-level stage guard, or is that a behavioral
  change out of scope? Impacts the `RobotState` accessor doc tickets.
- **OQ-2 - CLOSED (2026-07-13): the premise was stale; `jointHDot_FM` IS golden-tested.**
  Investigating how to test `jointHDot_FM` for BendStretch/SphericalCoords/FreeLine
  found it already validated by three independent, passing oracles, so there is no
  gap to carry forward and no new test to write (a fourth would be test-bloat):
  (1) the external Simbody differential `TestRoboticsOracle.{BendStretch,
  SphericalCoords,FreeLine}` (+ FuzzStates) compares the HDot-dependent `udot`/`A_GB`
  against Simbody fixtures at 1e-8 - this is exactly the "Simbody single-body
  reference" the old source caveat demanded; (2) `JointJacobianDot.HDotMatchesDerivativeOfHFM`
  checks HDot_FM == d/dt jointH_FM by central difference (all three types, 150 trials,
  1e-4); (3) `Integrator.UnderValidatedJointsConserveEnergy` +
  `BendStretchIsolatedConservesEnergy` are the energy-conservation dynamics oracle.
  The stale "NOT golden-tested" caveat in `src/JointKernels.cpp` was corrected to cite
  these three. `SPLIT-R4` is mergeable with no residual validation debt.
- **OQ-3 - DECIDED: keep as oracle; add a behavioral oracle.** `runREX` is not
  deprecated junk - it is the coordinate-swap baseline deliberately frozen as the
  INVARIANT-EQUIV differential oracle for the label-swap `RunREX` ("do not
  extend", `Context.hpp:84-90`). It stays. Separately, the REX suite has only
  *algebraic* oracles today; a theory-derived **stationary-distribution /
  detailed-balance** oracle for the REMC/Default chain will be added as its own
  behavioral track before any Context-REX code motion (see `TESTS.md section 7`).
- **OQ-4** `Inertia` vs `SpatialInertia` vs `MassProperties` overlap heavily.
  Consolidate, or is each still used at a call site outside the solver? Needs a
  whole-repo usage scan before any merge.
- **OQ-5** Driven-REX code (RENE/REBASONTOP: `runDrivenRound`, `driveReplica`,
  `runInterleavedRemcSubround`) is documented as reviewed-on-paper but never
  compiled/run. Does it stay in the split scope, or move behind a feature flag
  until it has a test?
- **OQ-6** `RobotState` carrying NCMC-solvent fields: relocate to a separate
  separate auxiliary state, or accept the concern bleed to preserve the single-arena
  allocation?
