# Spec: Robotics oracle — live-Simbody differential for the ABA port

Status: ready for review
Owner: coder (after review)
Gate: every §7 quantity matches the clone-generated reference within §6 tolerances, on the §8 battery — synthetic per-joint + depth/conditioning stress (Scope A) and all three molecular topology classes rigid/regular/cyclic (Scope B, §4.4) — wired into `nox -s tests`

## 1. Problem

The internal-coordinate dynamics engine (`disasm`, in `./src` + `./include`) is a
Simbody-free reimplementation of Featherstone's Articulated Body Algorithm (ABA),
transcribed from the still-vendored SimTK `Simbody01/` tree that now lives only in the
`refactor` clone at `./Robosample/`.

When Simbody left the port build, the historical *"diff every operator against the
live Simbody build"* gate was deleted and **never replaced for robotics**
(`src/RobotEngine.cpp:22-37`). The current robotics oracle-of-record is analytic
forms, finite-difference cross-checks, operator round-trips, statistical ensembles,
**hand-transcribed frozen Simbody constants**, and a live OpenMM diff for potential
energy only.

**The structural hole (confirmed by auditing every robotics test — §3):** the existing
suite is almost entirely *self-referential*. Every deep ABA intermediate
(`P, PPlus, DI, G, Z, zPlus, eps, A_GB, udot`, reaction forces) is validated only
against the port's **own** code — a finite difference of the port's kinematics, an
analytic identity, or a round-trip between two port operators. Self-consistency cannot
catch an internally-coherent transcription error; the port's own docs name one that did
exactly that (the "C-4" quaternion-handedness bug that "hid under self-consistency",
`include/robot_math.hpp:260-269`). No current test pins the *absolute value* of a deep
ABA intermediate against an external ground truth.

Goal: a **live numerical oracle** sourced from the clone's real Simbody that pins those
absolute values, wired into the standing gate — closing the one hole the existing suite
structurally cannot cover, **without duplicating** what analytic/FD checks already cover.

## 2. Why the port is oracle-shaped (background)

Two design choices make this cheap:

1. `RobotState` (`include/RobotState.hpp`) is a SoA mirror of Simbody's
   `SBTree{Position,Velocity,Acceleration}Cache` + `SBArticulatedBodyInertiaCache`.
   Every ABA intermediate is a named array: `q,u,qdot,udot,qdotdot`;
   `X_GB,X_FM,X_PB,Phi,Mk_G,H_FM,H,V_FM,V_PB_G,V_GB,A_GB,gyro,coriolisA,mobCoriolisA`;
   `P,PPlus,G,DI,abCentrifugal,Z,zPlus,eps`; `bodyForceG,mobilityForce`.
2. `RobotEngine` (`include/RobotEngine.hpp`) is static methods over
   `(const RobotModel&, RobotState&)`, each a 1:1 transcription of a named Simbody
   routine (source map at `src/RobotEngine.cpp:12-21`).

Units are identical on both sides — the consistent MD system (nm, dalton, ps, kJ/mol),
`include/Units.hpp:8-14`. **No unit conversion in the oracle.**

The ABA is *exact* (O(n)): for the same model and inputs the two engines must agree to
near machine precision. Two divergences must be understood precisely:

(a) **Floating-point summation order** — bounded, `~n·eps`; a tolerance concern only, never
a blow-up.

(b) **The hinge-inertia inverse method** — Simbody factorizes `D`; the port uses a
cyclic-Jacobi eigensolver with a **hard `1e-12` null-space lock**
(`src/RobotEngine.cpp:318,338,343-352`): an eigendirection contributes either `1/λ` or
**exactly 0**. This is **bounded only AWAY from the lock**. As an eigenvalue `λ` of `D`
approaches the `1e-12` threshold, the two engines diverge *unboundedly*: the port zeroes
the direction (`0`) while Simbody keeps `1/λ` (`~1e12`), or vice-versa across the
threshold — a **step discontinuity**, not a smooth `cond·eps` error. Simbody may also
reject a genuinely singular `D` at realize (it assumes SPD hinge inertia). **Consequence:**
element-wise differential testing of `DI/PPlus/G` is valid **only** when `min-eig(D)` sits
well above the lock on *both* sides; near/at the lock the divergence is by-design (§2
declares the lock legitimate) and must be tested port-only, never diffed. This governs the
conditioning battery (§8.1) and the stage-3 gate (§6).

## 3. Relationship to existing tests (non-redundancy — binding)

A full audit of the robotics tests (`tests/Test{SpatialAlgebra,Inertia,MassMatrix,
Mobilizer,MobilizerKinematics,ReverseMobilizer,JointKernels,ReactionForces,
KineticEnergy,Integrator,FixmanBoltzmann,MassScaleInvariance,Quaternion,Transform,
RotationConstruction,LinearAlgebraOracle,NMALinearAlgebra,AtomTransfer,Transfer,
Geometry}.cpp`) yields the coverage map below. The oracle's scope is defined **by
subtraction**: it only asserts what is currently checked *solely* by self-consistency.

| ABA quantity | Existing check kind | Existing test (file) | Oracle action |
|---|---|---|---|
| `X_GB`, `V_GB` | FD vs port kinematics | MobilizerKinematics, ReverseMobilizer (1e-10) | **ADD** live diff (cheap, anchors stage 1/2; FD is blind to a coherent kinematics bug) |
| `A_GB` | FD (`d/dt V_GB`), 1e-3 | Mobilizer | **ADD** live diff (FD tol is loose at 1e-3; oracle tightens) |
| `X_FM`, `H_FM`, `HDot_FM` per joint | analytic spot-check + FD, 1e-4 | JointKernels | **ADD only for BendStretch/SphericalCoords/FreeLine** (q-dependent H, highest risk); constant-H joints already well-pinned — do **not** re-diff |
| `V_FM` | independent re-derivation, 1e-10 | ReverseMobilizer | **SKIP** (already externally-derived, not self-referential) |
| `P` | **none direct** — only via `M`/`M⁻¹` round-trips | MassMatrix | **ADD** absolute pin via public `getArticulatedBodyInertia(state,mbx)` (singular; the only public ABI accessor) |
| `PPlus, DI, G` | **none direct** | MassMatrix | **ADD** absolute pin — **requires the clone-getter patch (§4.3, now binding)**; there is no public accessor |
| `Z, zPlus, eps` | **none direct** | — | **ADD** absolute pin — **requires the clone-getter/capture patch (§4.3)**; `getArticulatedBodyForces` does **not** exist and Simbody's centrifugal-force accessors are a *different* quantity (velocity-only, no `bodyForceG`/inward accumulation) |
| `udot` | `= M⁻¹τ` at `u=0` (round-trip); determinism | MassMatrix, Mobilizer | **ADD** live diff at **nonzero velocity** (the round-trip only covers `u=0`; the velocity-coupled path is untested absolutely) |
| `M`, `M⁻¹`, `√M⁻¹` | forward-Jacobian dense oracle + round-trips | MassMatrix | **ADD** narrow diff of dense `M` + eig(`M`) (see correctness note below); skip the round-trips |
| `logDetM` | dense-`M` + Fixman-grid, 1e-7 | MassMatrix, FixmanBoltzmann | **ADD** live scalar diff (basis-independent, cheap cross-check) |
| `KE` | analytic rigid-body oracle | KineticEnergy | **SKIP** (already externally-anchored) |
| reaction forces | independent Newton-Euler **re-derived in the same test** | ReactionForces | **ADD** diff vs Simbody `findMobilizerReactionOnBody…` (the re-derivation is still the port's own physics; Simbody is the external check) |
| `Transform`, `Rotation`, quaternion `N`-map, dihedral, pitch, cross product | analytic + frozen constants | Transform, Quaternion, RotationConstruction, Geometry, Inertia | **SKIP** — pure value algebra, already analytic; a live diff adds only maintenance cost |

**Correctness note the oracle specifically closes.** `TestMassMatrix`'s "independent"
forward-Jacobian dense `M` oracle and the ABI `M⁻¹` path both consume `Mk_G` and `V_GB`
from `realizePosition`/`realizeVelocity`. A bug *there* corrupts both and the round-trip
still passes. Simbody computes `X_GB`/`V_GB`/`M` from its own independent code, so the
live diff is the check that is not blind to a shared-input kinematics error.

**Scaffolding reuse (binding — do not create parallel infrastructure):**
- `tests/RobotBuilders.hpp` — `buildForest(specs)`, `buildSingle`, `buildChain`,
  `attachAtoms`, `randomizeState`, `BodySpec`. The oracle builds robots from these.
- `tests/TestHelpers.hpp` — the `Rng` and the tolerance tiers (`kTight/kAlg/kFD/kLoose`)
  and comparators (`NearVec3`, `NearMat33`); §6 tolerances reuse these tiers.
- `tests/RobotLinearAlgebra.hpp` — `jacobiSymEig`, `logDetSymPD` for the dense-`M`/eig
  comparison.
- `tests/AnalyticForceBridge.hpp` — conservative harmonic bridge (no OpenMM) for states
  that need applied forces.
- `tests/HmcDriver.hpp`, `tests/StatTest.hpp` — not needed by the oracle (it is a
  single-state deterministic diff, not a sampler); listed so they are not re-invented.

**Frozen-golden constants become derived artifacts.** The hand-transcribed Simbody
numbers (`TestSpatialAlgebra.cpp:37-59`; `TestInertia.cpp:128-144,203-212`;
`TestJointKernels.cpp:161-209`; `TestKineticEnergy.cpp:141-174`;
`TestQuaternion.cpp:81-94`; `TestRotationConstruction.cpp:44-56`;
`TestGeometry.cpp:68-79`) are kept as fast, dependency-free regression guards, but the
generator (§4) becomes their **authoritative source**: the same driver that produces the
fixtures also emits these constants, so a transcription typo is caught at generation
time instead of living forever as a possibly-wrong literal. This does not change the
existing tests; it changes where their numbers come from.

## 4. Architecture (decided: fixture generated from the live clone)

Neither Python surface exposes the recursion internals (port pybind: topology + energy;
clone pybind: `getAdvancedQs/Us` = q,u only). The deep quantities live in C++ only —
clone via `World::getMatterSubsystem()` (`Robosample/include/World.hpp:514`); port via
gtest access to `RobotState`. So the oracle is C++-driven on both sides:

- **Generator** (`Robosample/tools/gen_robotics_oracle.cpp`), built by the clone
  toolchain (Simbody is already a first-class dep there). Reads the **serialized model
  spec** (§4.1), builds the corresponding Simbody `MultibodySystem`, sets the scripted
  `(q, u, applied-force)` state, realizes through Acceleration, and dumps every §7
  quantity to a versioned fixture. Driver template: `Robosample/Simbody01/Simbody/tests/
  TestMassMatrix.cpp`.
- **Port test** (`tests/TestRoboticsOracle.cpp`), built in the existing `cuda-tests`
  config with **zero new dependencies** (a fixture is data). Reads the same model spec,
  builds a `RobotModel` via `RobotBuilders::buildForest`, replays the same inputs through
  `RobotEngine`, compares array-by-array (§6).

### 4.1 The model spec is serialized, not re-generated per side

The two sides are separate binaries; they **cannot** each RNG their own frames/masses
(`RobotBuilders::buildSingle/buildChain` take an `Rng`). The model must be defined once
and shared. The serialized model spec (source-embedded, §9) therefore stores, per body:
`{parent, jointType, X_PF, X_BM, mass, com_B, inertia_B}`, plus the input state
`{q, u, bodyForceG, mobilityForce}` and a `schema_version`. Both sides build from these
exact values. This removes any "did the two RNGs agree" ambiguity and makes the model the
ground-truth contract. (For Scope B the spec is *derived* from the clone's molecule build,
§4.2, not hand-authored.)

### 4.2 Two scopes (both in-gate), staged in build order

- **Scope A — synthetic robots (build first).** Model specs are hand-authored `BodySpec`
  lists (one `buildSingle` per `JointType`, a mixed chain, a branched forest, plus the
  §8 depth/conditioning stress cases), authored deterministically once and serialized.
  Isolates the ABA math from the molecule-construction layer; this is where the
  least-covered, highest-risk quantities live (§3) and it is the cheapest to build.

  **JointType -> Simbody mobilizer map (verified against the native nodes the port
  transcribes, `src/JointKernels.cpp:6-61`; NOT Custom mobilizers):**
  `Torsion->Pin` (about Z), `Slider->Slider` (on X), `Cylinder->Cylinder`, `Ball->Ball`
  (quaternion), `Cartesian->Translation`, `Free->Free`, `Rigid->Weld`, and the three
  high-risk joints to their **native** counterparts of the same name:
  `BendStretch->MobilizedBody::BendStretch` (`MobilizedBody_BendStretch.h`: rotate about
  Fz==Mz, slide along rotated Mx; matches `JointKernels.cpp:26-31`),
  `SphericalCoords->MobilizedBody::SphericalCoords` **default ctor**
  (`MobilizedBody_SphericalCoords.h:88-92`: azimuth about Fz==Mz, zenith about My, radius
  along Mz, radialAxis=Z, no axis negation, zero offsets -- state these ctor args
  explicitly in the generator; matches `JointKernels.cpp:43-56`),
  `FreeLine->MobilizedBody::FreeLine` (quaternion default, nq=7/nu=5; **not**
  `LineOrientation`, which is 2-DOF orientation-only). Do **not** use
  `Custom`/`FunctionBased`/`LineOrientation` -- their q/u parametrization differs from the
  port's transcription and would false-FAIL stage 3/4 on exactly the priority joints.
  Frame-equality gate (per joint): before any frame-dependent comparison, assert port
  `X_FM` == Simbody `getMobilizerTransform(state)` at `q=0` and one random `q`, to
  machine precision -- this discriminates a native-vs-Custom parametrization mismatch
  immediately.

  **Caveat (a green Scope A is not a green joint-frame build):** the clone normally
  constructs mobilizers through the Molmodel `BondMobility` layer with frame composition
  (`X_to_Z`, spherical `X_to_Y*Y_to_Z`, `Robosample/src/World.cpp:306-371`). Scope A's
  raw `MobilizedBody` construction bypasses that path, so it validates the ABA math
  cleanly but does **not** validate that the port's *molecule-build* frames match
  Molmodel -- that check only arrives in Scope B.
- **Scope B — molecule / pipeline model (in-gate).** Scope B SHALL cover **all three
  molecular topology classes** (§4.4). Concrete build path (resolved):

  - *Clone side (generator):* build each molecule via the clone's own
    `World`+`Compound`+`DuMM`+Molmodel pipeline from the AMBER prmtop → Simbody
    `MultibodySystem`; dump the §7 recursion quantities AND the decomposition (per-body
    atom-index set, joint type, joint frames, q/u layout). NOTE: the heavy clone build
    (Molmodel + AMBER loader + OpenMM) is confirmed working (`cuda-release`, 609 targets,
    Stage 0). NOTE (resolved in Stage 0): the clone has **no C++-only AMBER-prmtop entry
    point** — prmtop parsing is Python-only (ParmEd/mdtraj, `Robosample/python/robosample/
    context.py`), so the generator is NOT a standalone C++ binary like Scope A's. Invert the
    driver: a Python script orchestrates the clone's molecule build via its existing pybind
    `World`/`Context` API, and a small **additive** clone-side pybind C++ function walks
    `getMatterSubsystem()` bodies and dumps the §7 quantities via the §4.3 Oracle getters to
    npz+manifest. Requires an additive `getWorldState()`/`updWorldState()` on
    `Robosample/include/World.hpp` (free rein). Cross-engine atom identity is `prmtopIndex`
    (`SystemTopology.atomsPrmtopIndex` port ↔ `RoboAtomIdentity.prmtopIndex` clone).
  - *Port side (test):* `SystemTopology` (`include/TopologyElements.hpp:32`) is a plain POD
    of vectors and `World::buildModel(SystemTopology, Selection, rootMobilities)` is pure
    C++. So the port test SHALL obtain its `RobotModel` by loading a `SystemTopology` dumped
    by the port's OWN Python AMBER loader (`context.load_amber` → `context.system_topology`,
    the same path `loader_differential` uses; 10ala/1APQ pkls already exist under
    `tests/fixtures/loader_golden/`), then calling `World::buildModel` and running
    `RobotEngine`. A small Python fixture-gen step (analogue of `loader_golden/_generate.py`)
    SHALL serialize the `SystemTopology` into a C++-loadable form (npz/json, reusing the §9
    machinery), per selection: all-Rigid for the rigid class, default-flexible for regular,
    default for cyclic.
  - *Correspondence:* match port bodies ↔ Simbody bodies by identical atom-index set (§5),
    then run the §5.1 topological preflight (T1/T2 on the body graph; T3 port-only for the
    cyclic constraint count) BEFORE the numeric comparison. The atom-set map is what lets
    the two independently-built models be compared; a body-set mismatch is itself a finding.

  This is where the model-construction path — not the recursion — differs sharply between
  the three classes, and the cyclic class also exposes the §4.5 constraint divergence. NOTE:
  verifying the port's `World::buildModel` reproduces the clone's *decomposition* (beyond
  the CGK structural checks) is bounded by whether the two Python loaders agree; the numeric
  recursion differential + the §5.1 invariants are the primary Scope-B gate, and a
  Python-level decomposition differential MAY be a separate follow-up.

### 4.4 The three molecular topology classes (Scope B coverage)

The `World::buildModel` / Molmodel decomposition path behaves differently per topology;
each must appear in the battery (concrete molecules fixed in §8):

1. **Fully rigid (no internal DOF)** — `examples/10ala.*`, all bonds Rigid. The whole
   112-atom peptide is one rigid body. Exercises the degenerate single-body tree and the
   two legal root mobilities separately: a **Weld/Rigid** root (0 total DOF — frozen;
   `X_GB` constant, `udot=0`) and a **Free** root (6-DOF whole-body motion — the classic
   free-rigid-body case, also analytically checkable). Confirms both builds agree on mass
   properties and the root frame.
2. **Regular (linear / branched tree, no ring closures)** — `examples/10ala.*`, default
   flexible mobility (same molecule as case 1 → controls for everything but mobility). A
   spanning tree of flexible mobilizers; exercises the full recursion and the branched
   inward inertia accumulation. Primary Scope-B correctness case.
3. **Cyclic (contains ring-closing bonds)** — `examples/1APQ.*` (prolines + disulfide
   `CYX` cross-links). The ring-closing bond is **never a tree edge** in either engine
   (`RobotModel.hpp:46-49`, `World.cpp:377,410,608`), so the *spanning tree and its
   recursion are identical to the acyclic case* — but the two engines treat the closure
   itself differently (§4.5). This class validates that ring-closure **detection/skip
   agrees** (both build the same tree despite the extra bonds) and that the port's
   constraint layer does **not** contaminate the unconstrained recursion outputs.

### 4.5 Cyclic molecules — the constraint carve-out (correctness-critical)

The port and clone run **different dynamics** for a cyclic molecule, so the differential
must be scoped to the quantities that are still identical:

- **Port (`disasm`):** each ring-closing bond is a constant-distance (Rod) **holonomic
  constraint** enforced by SHAKE (position) + RATTLE (velocity), removing one flexible
  DOF per ring, with its own constraint Jacobian `G`; the port even folds the
  loop-closure term into its Fixman log-det (`include/Constraints.hpp:33,80-83,104`;
  `src/World.cpp:851-852`).
- **Clone (`refactor`):** does **not** enforce the closure as a constraint — the
  ring-closing bond is only a harmonic force term, integrated with a small timestep (it
  *has* constraint machinery but does not use it for this).

Consequence for the oracle:
- **Comparable (differential ON):** the *unconstrained* spanning-tree recursion given the
  same inputs — `X_GB, V_GB` (position/velocity), `P/PPlus/DI/G, Z, A_GB, udot` from the
  raw `RobotEngine::calcUDot` output **before** any SHAKE/RATTLE projection, and the
  *unconstrained* `calcLogDetM` (tree `ln|M_φ|`, no `G` term). These are identical in both
  engines because the closure is not a tree edge and the projection is a separate
  post-step. The applied ring-bond harmonic force, if included, must be identical in
  `bodyForceG` on both sides.
- **NOT comparable (differential OFF — no clone counterpart):** the SHAKE/RATTLE
  projection, the constraint Jacobian `G`, the constraint-corrected (Lagrange-multiplier)
  `udot`, and the Fixman log-det *with* the loop-closure term. These are validated by the
  port's own constraint tests (`TestConstraints.cpp`, `TestConstraintSolver.cpp`) and
  invariants (`C(q)=0`, `G u=0`), never against the clone. The port test must compare the
  **pre-projection** `calcUDot` output for the cyclic case, or it will false-FAIL
  constrained-vs-unconstrained.

Comparability **requires both engines to break the same ring bond** (identical spanning
tree). §4.4-item-3 makes tree-agreement a *validated invariant*: a tree-structure mismatch
(e.g. fused rings where Molmodel breaks a different bond than the port, §8.2 #8) is a
**real build bug to surface**, not a difference to scope out. Verified sound: the clone's
only `SimTK::Constraint::Rod` path (`Robosample/src/World.cpp:1746-1815`) and its caller
(`Robosample/src/Context.cpp:682`) are **dead code** (commented out, docking-specific), so
the clone genuinely leaves ring closures as DuMM harmonic bond-stretch only; and the port's
`calcUDot` `s.G()` is the *articulated* kernel `G=P·H·DI`, not the constraint Jacobian, with
RATTLE running strictly after (`RobotIntegrator.hpp:409`). Use `calcAcceleration-
IgnoringConstraints` on the clone side (§7.1) so the path stays unconstrained even if a
future clone patch re-enables the dead Rod code.

### 4.3 Clone-getter patch (BINDING for the flagship intermediates)

Simbody's public API exposes **only** `getArticulatedBodyInertia(state, mbx)` (singular,
`P`; `SimbodyMatterSubsystem.h:2829`). There is **no** public accessor for `PPlus`,
`DI`, `G`, or for `Z`/`zPlus`/`eps` (`getArticulatedBodyForces` does not exist; the
centrifugal-force accessors are a different, velocity-only quantity). Absolute-pinning
those — the stated primary value of this oracle — therefore **requires a small local
patch to the clone**, which is legitimate because the clone is a test fixture and is not
shipped:

- `PPlus, DI, G`: persisted in `SBArticulatedBodyInertiaCache`
  (`RigidBodyNodeSpec.cpp:274` fills them). Add const forwarders on
  `SimbodyMatterSubsystem`/`Rep` returning them per `MobilizedBodyIndex`.
- `Z, zPlus, eps`: these are **buffers of the acceleration pass**
  (`calcUDotPass1Inward`, `RigidBodyNodeSpec.cpp:382`), not necessarily persisted in a
  cache. The generator captures them by calling the pass with its own output buffers
  (the pass already takes `SpatialVec* allZ, SpatialVec* allZPlus, Real* allEpsilon`),
  rather than adding a getter. The coder must confirm at implementation time whether the
  realized `State` retains them; if not, use the capture path.

Fallback if a given term proves infeasible to expose: validate it **transitively** (via
reactions + `udot` + `M`) and mark it in the port test's field list as *transitively
pinned, not absolute*, so the residual self-consistency is documented, not silently
reintroduced. `P`, `M`, `logDetM`, `udot`, reactions, and all kinematics are absolutely
pinnable **without** any patch.

Out of scope (opt-in follow-up): Architecture A, an in-process dual-link binary behind a
`robotics-oracle` CMake preset for step-by-step bisection when the fixture gate flags a
divergence. Deferred because it re-introduces the Simbody build dependency the port
removed.

## 5. Correspondence — pin it, do not assume it

The port claims index-for-index equivalence with Simbody; that claim is *under test*, so
correspondence is by invariant and any mismatch is itself a finding.

- **Scope A**: the serialized spec IS the correspondence — same body order, same
  `jointType`, so `q`/`u` blocks align by construction; still assert per body that port
  `bodyNQ/bodyNU` (`RobotModel.hpp:74-83`) equal Simbody `getNumQ/getNumU`, and that the
  quaternion block is the first 4 q for Ball/Free/FreeLine (`RobotModel.hpp:567-575`).
- **Scope B**: match port-body ↔ ref-body by **identical atom-index set** (atom index is
  the shared currency: `RobotModel.hpp:134-141`; clone `cAIx`). The generator emits each
  reference body's atom-set + `qIndex/uIndex/nq/nu/jointType`; the port test builds the
  atom-set map and verifies DOF layout before any numeric comparison.

### 5.1 Topological preflight — Chebychev–Grübler–Kutzbach / Euler loop-count (decomposition correctness, engine-independent)

The tree-agreement check (§4.4/§4.5) proves the two engines used the *same* decomposition,
not that it is *correct*. The **atom-level** ring-count is already guarded on load
(`python/robosample/acyclic_graph.py:285-293` fail-louds if `len(ring_closing_set) ≠` the
atom-graph cyclomatic number, and `check_residual_cycles` guards acyclicity), so a pure
Python miscount does not slip through. The **un-guarded** gap is the **C++ atom→body
reduction** (`src/World.cpp`: the DSU rigid-merge 578-585, the BFS tree build 463-490, the
`atomBody` map, and the inter-body constraint filter 673-681) — a bug there is seen by
neither the Python validator nor the numeric port-vs-clone diff. T1/T2 are an
**engine-independent, geometry-independent** structural preflight targeting exactly that
reduction, run before any numeric comparison, on the **body-adjacency graph** (nodes =
rigid bodies, edges = inter-body bonds):

- **Invariant T1 (body-tree acyclicity).** Rebuild the body graph **independently from the
  bond list + `atomBody`** — every bond whose endpoints lie on *distinct* bodies
  (`atomBody[i] ≠ atomBody[j]`) and that is **not** `bondsRingClosing`. Its
  mobilizer-edge subgraph must be a **spanning forest**: cyclomatic `0`, edge count
  `= N_bodies − N_components`. If a loop-closing bond was mis-reduced into a tree edge, this
  subgraph contains a cycle → **fail loud**. Critically, build this from the bond list,
  **NOT from `RobotModel.bodyParent`** — the port's BFS silently drops back-edges
  (`World.cpp:471`), so `bodyParent` is a forest *by construction* and checking it can
  never fail (a Rule-8 vacuous test). The discriminating regression test injects a real
  ring bond with `bondsRingClosing=false` and asserts T1 fails.
- **Invariant T2 (loop count).** The number of **active** ring-closure constraints —
  `ConstraintSet.numConstraints()`, i.e. the **inter-body** closures Rod-constrained at
  `World.cpp:678` — must equal the body-graph cyclomatic number
  `L = E_body − N_bodies + N_components`. **Do not** assert this against the raw
  `bondsRingClosing` count: a ring collapsed entirely inside one rigid body contributes to
  the atom-graph cyclomatic number but creates **no** body-graph edge and **no** constraint,
  so `numConstraints() ≤ #bondsRingClosing` in general (equality only when no ring is
  intra-body-rigid). §8 must report **both** counts and assert only `numConstraints() ==
  bodyGraphL` and `#bondsRingClosing ≥ numConstraints()`.
- **Invariant T3 (CGK mobility bookkeeping — PORT-ONLY, a bound, not an unconditional
  equality).** Tree DOF `nu = Σ fᵢ` over mobilizers (Chebychev–Grübler–Kutzbach
  `M = 6(N−1−J)+Σfᵢ` with `J = N−1` for a spanning tree, so the `6(·)` term vanishes).
  Constrained DOF `= nu − rank(G)`, `G` = the port's loop-closure constraint Jacobian. `G`
  has exactly `numConstraints()` rows, so `rank(G) ≤ L` always; generically each Rod is
  independent so `rank(G) = L` and `M = nu − L` — true for the biomolecular battery (1APQ:
  independent disulfide/proline distance constraints in generic geometry). A strict
  `rank(G) < L` means **redundant/degenerate loop constraints** (physical overconstraint, à
  la Bennett/Bricard) — flag for investigation, **not** an automatic bug. **`rank(G)` is not
  exposed**: assemble `G` (`numConstraints × nu`) row-by-row via
  `Constraints::mapAtomForcesToGeneralizedForces` (`include/Constraints.hpp`) and compute
  rank explicitly (SVD, or `jacobiSymEig` on `GGᵀ`). This validates DOF *bookkeeping*, not
  the recursion, and lives with the port-only constraint tests, not the differential.

Cheap, exact, and complementary to the numeric diff: T1/T2 audit the C++-side *structure*
the numerical oracle then assumes. The §8.3 random-topology fuzz generator reuses the same
CGK identity to validate its generated models (predict `M` from the graph, assert the built
`RobotModel.nu` matches) before baking them.

## 6. Comparison staging and tolerances

Stage the rollout; each layer green before the next, or a position bug cascades and
buries its root cause. All comparisons are tolerance-based (never `==` on floats) and
**fail loud, naming the first offending (body, quantity, index)** (Rule 11). Reuse the
`TestHelpers` tolerance tiers.

Frame discipline (critical): compare **Ground-expressed / frame-invariant** quantities
first. If the port chose an equivalent-but-different internal F/M frame, the *physical*
Ground quantities still match while `X_FM/H_FM/Phi/P/Z` legitimately differ. In Scope A
(serialized shared model spec) the two engines share the body-frame convention, so a
frame-equality assertion between the serialized `X_PF/X_BM` and the generator's mobilizer
frames gates the frame-dependent comparisons.

NOTE (Scope B — resolved by review): across two INDEPENDENTLY-built molecule models the
body-frame *orientation* convention differs — the port's `recomputeGeometry` frame-graph
orients world-root and childless/terminal bodies under-determinedly (roots → identity
`Rotation`; no-first-child → one axis), while Molmodel derives every frame via BondCenter
(`src/World.cpp:867-885` documents that `X_GB.R()` carries only within-block deviation and
must not be read as an absolute orientation). Therefore **`X_GB.R()` is NOT an admissible
cross-engine anchor for Scope B** — the body-frame orientation and `atomStation_B` are built
from the same transform and cancel in `posG = X_GB.p + X_GB.R·station_B`, so a
different-but-self-consistent frame yields identical physics. The Scope-B invariant anchors
SHALL be: **per-atom Ground positions**, `X_GB.p()` (body origin), Ground-expressed `V_GB`/
`A_GB`, `udot`, and the scalars (`KE`, `logDetM`, `eig(M)`). `X_GB.R()` remains a valid
anchor only in Scope A, where the frames are shared by construction.

| Stage | Quantities | Frame class | rtol | Note |
|---|---|---|---|---|
| 1 | `X_GB`, per-atom positions | Ground, invariant | 1e-10 | pure kinematics |
| 2 | `V_GB`, then `qdot` | Ground, invariant | 1e-9 | quaternion: input q stored pre-normalized (§7.1); compare orientation via `R_FM`, not the raw 4-vector — double-cover sign AND magnitude (Simbody re-normalizes in its `State`; the port stores raw q and does not re-orthonormalize `Rotation(Mat33)`, `TestQuaternion.cpp:372-385`) |
| 3 | `P`, `PPlus`, `DI`, `G` | frame-dependent | 1e-8 | **gated on the §6.1 lock check**; element-wise only when `min-eig(D)` is well above the lock on both sides |
| 4 | `Z`, `zPlus`, `eps`, then `udot`, `A_GB` | mixed; `udot`/`A_GB` invariant | 1e-8 | key physical end products |
| 5 | dense `M`, eig(`M`), `logDetM`, reactions at `Bo`/`Mo` in Ground | invariant | 1e-8 (1e-6 eig) | basis-independent cross-check on stage 3 |

### 6.1 Stage-3 lock gate (J1-A — mandatory before any `DI/PPlus/G` element diff)

`DI/PPlus/G` are trustworthy only in the smooth regime of the null-space lock (§2b).
Before diffing any of their elements, per body:
1. assert `min-eig(D)_port ≈ min-eig(D)_Simbody` (both compute it), **and**
2. assert `min-eig(D) > lockTol` (the `1e-12` threshold, scaled).

If (2) fails, **skip the element-wise `DI/PPlus/G` diff for that body** and compare only
the invariant end-products (`udot, A_GB, M, logDetM`), which stay well-conditioned because
near-null directions carry ≈0 energy and cancel in the products. The
**conditioning-stress case** (§8.1) is deliberately tuned so `min-eig(D)` stays *above*
the lock (so neither engine locks and the §6.1 gate passes) — but note it is still
compared via §9 aggregate invariants + end-products, **not** element-wise `DI/PPlus/G`,
because there `DI` element error rides at the `1e-8` boundary and an element gate would be
flaky (§8.1/§9/§10, J1-B). Element-wise `DI/PPlus/G` diffs apply to the *well-conditioned*
structural cases only. The **true sub-`1e-12` singular hinge** is a **port-only** invariant
test (§8.2 #2), never a differential — assert it yields finite `udot` with a zero component
on the locked direction, with no Simbody reference.

**Priority targets:** BendStretch, SphericalCoords, FreeLine (q-dependent `H_FM`, native
Simbody mobilizers of the same name per §4.2, `src/JointKernels.cpp:168-199`,
`src/RobotEngine.cpp:648-660`) — the
joints §3 shows are covered only by kernel-FD and energy conservation. Test these first,
on `buildSingle` specs where an analytic check also applies.

## 7. Quantity contract

### 7.1 Inputs stored in the fixture (set identically on both sides)
`q` (quaternion blocks stored **already unit-normalized**, §6 stage 2), `u`, per-body
applied spatial force `bodyForceG` (`[torque; force]` about the body origin, Ground),
`mobilityForce`. `bodyMassScale` **must be 1.0** on both sides (it scales
`Mk_G` → `P/DI/logDetM/KE`, `RobotModel.hpp:98-113`).

**Applied-force injection into Simbody (S3 — the only path that exercises the force term
in `Z`/reactions).** The generator must feed `bodyForceG` as Simbody's
`Vector_<SpatialVec>& appliedBodyForces` (Ground frame, about each body origin — exactly
the port convention) directly through
`calcAccelerationIgnoringConstraints(state, appliedMobilityForces, appliedBodyForces,
udot, A_GB)` (or the equivalent explicit-force realize), **not** via a `Force::Custom`
element. Add a readback assert that the force Simbody reports for each body equals the
fixture `bodyForceG` before trusting any `Z`/reaction diff — a sign flip or wrong station
here would agree on the two zero-force states and only surface on the third.

### 7.2 Outputs — port array ↔ Simbody accessor (only the §3 "ADD" rows)

| Port (`RobotState`) | Simbody accessor (via `getMatterSubsystem()` / `State`) | Availability |
|---|---|---|
| `X_GB[b]`, `V_GB[b]`, `A_GB[b]` | `getMobilizedBody(b).getBodyTransform/getBodyVelocity/getBodyAcceleration(s)` | public |
| `X_FM[b]`, `V_FM[b]` (risky joints only) | `getMobilizedBody(b).getMobilizerTransform(s)`, `getMobilizerVelocity(s)` (`MobilizedBody.h:341,370`) | public |
| `P[b]` | `getArticulatedBodyInertia(s, b)` (singular; `SimbodyMatterSubsystem.h:2829`) | public |
| `PPlus[b]`, `DI[b]`, `G[b]` | clone-getter patch §4.3 | **requires patch** |
| `Z[b]`, `zPlus[b]`, `eps[b]` | acceleration-pass capture §4.3 | **requires patch/capture** |
| `udot`, `qdot` | `state.getUDot()`, `getQDot()` | public |
| dense `M`, `logDetM` | `calcM(s,M)`; `calcDetM(s,…,detM)`→`ln` | public |
| reactions at `Bo`/`Mo` in Ground | `getMobilizedBody(b).findMobilizerReactionOnBodyAtOriginInGround/AtMInGround(s)` (`MobilizedBody.h:860,868`) | public |

`Phi/H` stay internal on the Simbody side: validate transitively via the observable
outputs above. Spatial ordering is `[angular; linear]` on both sides
(`SpatialVec=Vec<2,Vec3>`; port `robot_math.hpp:606-639`; port force `[torque;force]`
confirmed `RobotEngine.cpp:1172`).

## 8. Test battery

Deterministic, serialized specs (no per-side RNG). Per spec, ≥3 states: a **near-zero
reference** state (`u=0`, `q` at the joint's natural rest — **not literally `q=0` for
SphericalCoords/BendStretch, which is a coordinate singularity, §8.2 #1**); a fixed
pseudo-random `(q,u)` (seeded once, values baked in); and a state with nonzero
`bodyForceG`. The edge-case states (§8.2) and the fuzz batch (§8.3) extend this.

**Scope A — synthetic (structural + stress):**
- One `buildSingle` per `JointType` (all 10).
- One mixed `buildChain` (`Free→Torsion→Ball→BendStretch→SphericalCoords…`).
- One branched `buildForest` (a body with ≥2 children — inward `Σ_children Φ P⁺ Φᵀ`).
- **Depth-stress:** a long linear chain (≥50 bodies) — targets kinematic product-chain
  drift (§8.1 mechanism 1).
- **Conditioning-stress:** a heavy rigid subtree hung off a light joint / high-branch hub
  giving `cond(D)~1e6–1e8` **with `min-eig(D)` bounded well above the lock** (`≥ 1e-9·scale`,
  so *neither* engine locks) — this tests the smooth `cond·eps` regime (§8.1 mechanism 2).
  The genuinely singular (`λ<1e-12`) hinge is **not** here — it is a port-only invariant
  test (§8.2 #2), because near the lock the divergence is unbounded and by-design (§2b/§6.1).

**Scope B — molecule (the three §4.4 classes, concrete molecules fixed):**
- **Rigid:** `examples/10ala.*` loaded with **all bond mobilities = Rigid** (no flexible
  internal DOF) — the connected 112-atom peptide (`ACE+10×ALA+NME`) collapses to a single
  rigid body; test at **Weld** root (0 total DOF, `udot=0`) and **Free** root (6-DOF
  whole-body, analytic free-rigid-body check).
- **Regular:** `examples/10ala.*` with the **default flexible** mobility — a pure
  linear/branched tree with **zero ring closures** (alanine has no ring side chain).
  *Deliberate control:* same molecule as the rigid case, so the only variable is the
  mobility selection, isolating the flexible-vs-welded build path.
- **Cyclic:** `examples/1APQ.*` — multiple ring closures (prolines + several disulfide
  `CYX` cross-links) → also exercises fused/nested rings (§8.2 #8, the tree-choice
  disagreement risk); compare pre-projection recursion only (§4.5).

**Generator must assert (fail-loud), or the cases do not exercise what they claim:** rigid
10ala → `numBodies==2` (Ground + one rigid body); regular 10ala → **0** active ring-closure
constraints; 1APQ → **≥1** active constraint (its 3 disulfides are inter-body regardless of
mobility, so this holds robustly) **and identical spanning tree (broken-bond set) on both
sides** (§4.4/§4.5 tree-agreement invariant). Report **both** the atom-graph
`#bondsRingClosing` and the active inter-body `ConstraintSet.numConstraints()` for 1APQ
(they may differ if a ring is intra-body-rigid, §5.1 T2 — do **not** assert them equal).

### 8.1 Why small hand-authored systems are not sufficient (scaling discussion)

The ABA is O(n), but *agreement between two independent implementations is not
size-independent* — a battery of only tiny systems can return a false "airtight" verdict.
Two mechanisms make the port↔clone difference grow, and both are latent (≈`1e-15`) in a
3-body toy:

1. **Kinematic product-chain drift — scales with tree DEPTH.** `X_GB[b]` is a compose
   chain from Ground down to `b`; round-off accumulates ~`O(depth)·eps`. This FP
   accumulation is present in **both** engines (neither re-orthonormalizes `Rotation`
   mid-chain — both use a plain `Mat33` multiply in the ABA compose path); the port↔clone
   *difference* is bounded by summation-order + root-quaternion normalization, **not** by
   the port's straight-copy `Rotation(Mat33)` (`robot_math.hpp:408`,
   `TestQuaternion.cpp:376-385`), which never sees a dirty `Mat33` in this path. Depth 50 ≈
   `1e-13`; a deep macromolecule can approach the stage-1 `1e-10` tol — survivable for the
   specified battery, but the reason to include the depth-stress chain rather than assume it.
2. **Inertia-reduction conditioning — scales with BRANCHING and mass dynamic range.**
   `D = HᵀP H`, `DI = D⁻¹`; when a heavy subtree hangs off a light joint or a hub has many
   children, `cond(D)` grows and (in the **smooth regime, `min-eig(D)` above the lock**) the
   port's Jacobi vs Simbody's factorization diverge ~`cond(D)·eps`; summation order over
   many children adds to it. A toy has `cond(D)~10`; a heavy-on-light configuration reaches
   `cond(D)~1e6–1e8` → `DI` error ~`1e-8`. **This is the "special case" that appears only
   with more rigid bodies / mass heterogeneity per molecule.** Two consequences: (i) the
   conditioning-stress case is sized to `cond(D)~1e6–1e8` but kept *above* the lock so this
   smooth model holds; (ii) because element `DI` error rides at ~`1e-8` (its own tolerance),
   the stage-3 **element-wise** `DI/PPlus/G` diff is *replaced* for this case by the §9
   aggregate invariants (`min-eig(D)`, `logDetM`) plus the well-conditioned end-products
   (`udot, A_GB, M`) — do not park an element gate on its own tolerance boundary (J1-B).
   The distinct sub-`1e-12` singular regime is out of the differential entirely (§6.1).

Implication for storage (§9): the structural cases are small enough to embed as source
constants; the two stress cases have large per-body outputs, resolved by storing
**aggregate/extremal invariants** rather than full per-body dumps.

### 8.2 Edge-case battery (ranked; a naive battery misses these)

Ranked by likelihood × severity × oracle-blindness. Cases marked **PORT-ONLY** are
un-buildable or ill-defined on the Simbody side and are tested by port invariants, never
diffed — they must be *documented as such*, not silently dropped.

1. **SphericalCoords/BendStretch at `q=0` is a coordinate singularity — the first mandated
   state hits it.** SphericalCoords: zenith `q1=0` collinearizes the azimuth (Fz) and
   radial (Mz) axes and radius `q2=0` puts the body at the F origin → `H` cols dependent,
   `D` near-singular (`JointKernels.cpp:43-67`). BendStretch: `q1=0` (zero stretch) gives
   `H₀=[Fz;0]`, rotational DOF loses inertial coupling. **Fix:** the near-zero reference
   state offsets `q1,q2` off the singularity for these joints (§8). *Highest — else a
   priority joint false-FAILs on state #1.*
2. **Massless / near-massless body carrying a DOF — PORT-ONLY.** `D=HᵀPH→0` → the lock
   discontinuity (§2b); Simbody rejects a non-SPD hinge inertia at realize or emits inf.
   Test port-only: finite `udot`, zero on the locked direction. Document as un-buildable
   on the oracle.
3. **FreeLine spin-suppression — do NOT diff the raw `qdot` quaternion block.** `nq=7,nu=5`;
   the suppressed DOF is spin about the body line, so the `u→qdot` N-map is rank-deficient
   and *which* spin component is zeroed is convention-dependent (`JointKernels.cpp:58-66,
   129-139`). Compare `R_FM` + `V_GB` only for FreeLine, never the raw quaternion `qdot`.
4. **Wide star hub (≥8 children) — summation-order determinism.** The inward
   `Σ_children Φ P⁺ Φᵀ` accrues children in a possibly different order than Simbody →
   `PPlus` at the hub varies by reduction order. The §8 branched forest covers this only
   weakly; add a high-fan-out hub. Stage 3.
5. **Extreme anisotropic inertia (needle/disk).** One principal moment ≪ others feeds
   `cond(D)`; a needle on a Ball joint has a near-null rotational DOF about its own axis.
   Pairs with the conditioning-stress case (kept above the lock).
6. **Zero-DOF Rigid (Weld) body mid-chain.** `nq=nu=0` interior node; inertia/force must
   still transmit through it via `Φ`. Cheap catcher for DOF-index / inward-accumulation bugs.
7. **Applied-force discriminators.** A generic mixed force agrees under a
   `[torque;force]`↔`[force;torque]` slot swap; a **pure-torque** state and a **pure-force**
   state isolate the top-3 vs bottom-3 slots (the §10 convention pitfall). **Force on
   Ground (body 0)** must be ignored by both — assert it does not leak into `udot`.
8. **Fused/nested rings sharing atoms (Scope B cyclic).** Multiple constraints; the
   spanning-tree choice is where the port and Molmodel could disagree on *which* bond to
   break — the §4.4/§4.5 tree-agreement invariant, backed by the §5.1 T1/T2 loop-count
   preflight, earns its keep here (T1/T2 catch a wrong loop count *even if both engines
   agree*, i.e. a shared-Python decomposition bug). 1APQ (multiple disulfides) exercises it.
9. **Duplicate molecules in a forest.** Identical topologies stress per-molecule
   `qIndex/uIndex` offsets and the §5 atom-set correspondence map (must not alias copies).
10. **Ground-only system (`nq=nu=0` total).** Empty `udot`, `0×0` `M`, `logDetM=0`; guards
    the comparator and the corrector (`den+1e-30`) against empty-array UB.
11. **Ball/Free double cover + near-180°** — LOW for this oracle (`q` is a stored *input*
    consumed identically both sides; Shepperd branch selection is off the compared path;
    §6 already mandates a double-cover-aware `R_FM` metric).
12. **Wrapped torsion `q>2π`** — LOW (`X_FM` uses cos/sin, periodic and safe).

### 8.3 Randomized fuzz batch (decision: seeded, frozen, generated by the clone)

Fuzzing is **worth integrating**, in the *only* form compatible with Architecture B's
frozen fixtures and the deterministic fail-loud gate (Rule 9): the clone (which owns the
oracle) draws the randoms, computes references, and bakes **both** inputs and outputs.
Rejected: live property-based fuzzing (needs the deferred in-process Architecture A —
scope creep) and drawing fresh randoms at gate time (no baked reference → flaky). The
existing randomization (`randomizeState`, the 200-trial `TestJointKernels`, the
`.hypothesis/` dir) is **all self-referential** (FD vs the port's own kinematics); a
random batch against the *external* Simbody oracle is new coverage, not redundant.

- **Reuse** `tests/RobotBuilders.hpp::randomizeState` (Shoemake `q` + Gaussian `u`) and
  `tests/TestHelpers.hpp::Rng` with a **pinned seed**.
- **Random states over fixed topology (first):** per structural spec, N≈32–64 seeded
  `(q,u)` states — catches the state-dependent `H_FM(q)`/Coriolis/velocity-coupled `udot`
  bugs the 3 hand states miss.
- **Small seeded random-topology batch (second — directly attacks §8.1 scaling-blindness):**
  random depth∈[2,60], random branching, random mass ratio∈[1e-3,1e3], seeded once —
  samples the `depth×cond` plane the two stress cases only sample at two points.
- **Singularity filter (MANDATORY):** the generator **rejects/resamples** any state with
  `min-eig(D) < 10·lockTol` on the Simbody side, else it would bake a lock-discontinuity
  false-FAIL (§2b) as a "reference." Filtered configs are resampled or recorded as
  port-only.
- **Storage (§9):** bake **aggregate invariants per state** (`logDetM, KE, ‖udot‖`,
  max-residual body index), not full per-body dumps; a failing aggregate escalates to a
  per-seed `.bin` for bisection.
- **Counterexample → regression:** a failure is already deterministic (fixed seed). The
  clone-side generator re-runs the failing seed and **minimizes** it (reduce body count,
  zero `u` components, halve `q`) until the residual drops below tol, then emits the
  reduced `BodySpec` as a **new named structural fixture** — the fuzz batch discovers, the
  reduced fixture guards, the gate stays flake-free.

## 9. Storage format and gate wiring

NOTE: an earlier revision of this section embedded the reference data as `constexpr`
arrays in a generated C++ header. That was replaced (Phase-1b retrospective): the header
reached 1.1 MB / ~22k lines for Scope A alone, compiled into the test TU on every build,
produced unreviewable diffs, and would not survive Scope B (1APQ → hundreds of bodies) or
the fuzz batch. The provenance of each datum is documented separately in
`docs/specs/robotics-oracle-data-provenance.md`.

Storage format (decided): **one NumPy `.npz` per case for the reference arrays, plus one
`manifest.json` per case** for the reviewable metadata, under
`tests/fixtures/robotics_oracle/`. See the schema in
`tests/fixtures/robotics_oracle/RoboticsOracleTypes.hpp` (the in-memory POD both sides
fill; it no longer carries the baked data, only the field layout).

- The reference arrays (`X_GB, V_GB, A_GB, X_FM, P, PPlus, DI, G, Z, zPlus, eps, udot,
  qdot, Mdense, logDetM, reactions, minEigD`, and the §8.1 aggregates) SHALL be written as
  named arrays in the case's `.npz`. Ragged sizes SHALL be stored at their true extent
  (no `kMax` zero-padding).
- `manifest.json` SHALL carry: `schema_version`, the model spec (bodies/joints/frames/mass
  — the §4.1 correspondence, small and human-reviewable), the list of arrays with shapes,
  the RNG seed for fuzz cases (§8.3), the generator and Simbody-patch git SHAs, and a
  SHA-256 over the case's arrays.
- The read/write path SHOULD use `cnpy` (MIT, single-file, numpy-compatible) vendored
  under `tests/`; it uses zlib, already present in the toolchain. NOTE: this is not a
  heavyweight dependency in the sense the port avoids (no SimTK); a reviewer can also open
  any fixture with `numpy.load`.
- **Stress cases** (§8.1) store only aggregate/extremal invariants (`logDetM`, `KE`,
  `‖udot‖`, worst `min-eig(D)`, and the report-body kinematics), not full per-body arrays.
  A stress or fuzz case that fails on an aggregate MAY escalate to a full per-body `.npz`
  for that one case, for bisection.
- Generator `Robosample/tools/gen_robotics_oracle.cpp` (built by the clone) writes the
  `.npz` + `manifest.json` per case, deterministically (no live RNG; §5 of the provenance
  doc). Regenerated only when the battery changes; the regen command is recorded in each
  manifest.
- Port test `tests/TestRoboticsOracle.cpp` (gtest) loads the fixtures at runtime and runs
  the §6 staged comparison. It reuses the §3 scaffolding and adds no SimTK dependency.
- Versioning: fixtures are committed to plain git for now (matching
  `tests/fixtures/loader_golden/*.pkl`, ~1.9 MB plain). NOTE: move to Git LFS once the fuzz
  batch pushes the total into tens of MB (git-lfs is not currently installed).
- No-silent-gap guard: the port test SHALL check the manifest's array list against an
  explicit expected-field list keyed to `schema_version`; a newly added quantity fails
  until the fixtures are regenerated.

## 10. Pitfalls (encode as explicit asserts)

- **Null-space lock is a step discontinuity, not smooth** → never element-diff `DI/PPlus/G`
  unless `min-eig(D) > lockTol` on both sides (§6.1); at/below the lock the divergence is
  unbounded and by-design → port-only invariant test (§8.2 #2).
- SphericalCoords/BendStretch `q=0` is a coordinate singularity → the reference state must
  offset off it (§8.2 #1); do not use literal `q=0` for these joints.
- FreeLine → compare `R_FM`/`V_GB`, never the raw `qdot` quaternion block (§8.2 #3).
- `D⁻¹` method divergence in the smooth regime → stage-3/4 tol `1e-8`; lean on invariant
  eig(`M`)/`logDetM`; for the conditioning-stress case use invariants, not element `DI` (J1-B).
- F/M frame ambiguity → Ground quantities first; gate frame-dependent diffs on frame
  equality.
- Quaternion double cover (`q ≡ −q`) → double-cover-aware metric on quaternion blocks
  (the port already does this in `checkReversibility`).
- Summation order → tolerances, never float `==`.
- `bodyMassScale ≠ 1.0` silently changes the inertia layer → assert 1.0 both sides.
- Applied-force convention → confirm both treat `bodyForceG` as `[torque; force]` about
  the body origin in Ground before comparing `Z`/reactions.
- Native custom-kinematics joints (BendStretch/SphericalCoords/FreeLine) — the native
  Simbody mobilizer chosen in the generator (§4.2, *not* a Custom mobilizer) must match
  the port's `X_FM` convention (`src/JointKernels.cpp:8-79`); if `X_FM` diverges but
  `X_GB` agrees, the frames differ and only Ground quantities are comparable for that joint.
- Cyclic molecules — compare the **pre-projection** `calcUDot` output and the
  *unconstrained* `calcLogDetM` only (§4.5); diffing the port's SHAKE/RATTLE-corrected
  udot or constraint-augmented Fixman log-det against the clone's unconstrained dynamics
  is a guaranteed false-FAIL, not a bug.

## 11. Deliverables and order

1. Build the clone once (`./Robosample`, no `build/` yet); confirm `getMatterSubsystem()`
   returns realized caches, the §7.2 public accessors link, and land the §4.3 clone-getter
   patch (verify whether `Z/zPlus/eps` persist in the realized `State` or need the capture
   path).
2. Scope A generator + source-embedded references (§9) for the §8 battery; author the
   deterministic `BodySpec` sets and the §4.2 native `JointType→MobilizedBody` map (state
   the SphericalCoords/FreeLine ctor args explicitly). Land the §4.2 frame-equality gate
   before any frame-dependent comparison.
3. `tests/TestRoboticsOracle.cpp`: correspondence (§5) → staged comparison (§6) incl. the
   §6.1 lock gate → fail-loud reporting; wire into `cuda-tests`. Start with BendStretch/
   SphericalCoords/FreeLine (using off-singularity reference states, §8.2 #1/#3).
4. Bring stages 1→5 green on the structural cases; add the §8.2 edge cases (mark the
   PORT-ONLY ones as invariant tests, not diffs); then the depth/conditioning stress cases
   (§8.1, kept above the lock).
5. Scope B (in-gate): the three §4.4 topology classes — rigid = 10ala all-Rigid (Weld +
   Free root), regular = 10ala flexible (assert 0 ring closures), cyclic = 1APQ
   (pre-projection recursion only, §4.5, with the tree-agreement invariant; report ring-
   closure count). Run the §5.1 topological preflight (T1 tree-acyclicity, T2 loop-count)
   on every case before numeric comparison; add the T3 CGK DOF-bookkeeping check to the
   port-only constraint tests for the cyclic case.
6. §8.3 fuzz batch (seeded random states, then random topologies) with the mandatory
   singularity filter and counterexample-minimization. (Deferred, opt-in) Architecture A
   preset — the only place live property-based fuzzing belongs.

## 12. Open questions for review

- Whether the conditioning-stress case (§8.1 mechanism 2) is best expressed synthetically
  (Scope A, tuned mass ratio) as specified, or additionally as a real large molecule in
  Scope B — the synthetic form is the controlled worst case; a real large system adds
  realism at higher fixture cost.

*Resolved in review / by user direction:* (a) storage is source-embedded generator-emitted
flat arrays, not NPZ (§9); (b) Scope B (all three topology classes) is in this spec's
gate; (c) the `PPlus/DI/G/Z/zPlus/eps` accessor question is closed — the clone-getter/
capture patch is binding (§4.3); (d) fuzzing is IN as a seeded frozen clone-generated batch
(§8.3), live property-based fuzzing deferred with Architecture A; (e) the null-space lock is
handled as a step discontinuity — smooth-regime element diffs gated on `min-eig(D)` (§6.1),
singular regime tested port-only (§2b review finding J1-A).
