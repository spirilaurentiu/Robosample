# Reaction-force monitoring (per-frame per-body forces to CSV, for analysis + VMD)

Status: implemented. Builds on already-shipped machinery; no new physics.

## 0. Purpose and scope

Reinstate, in the current SimTK-free engine, the spatial-force diagnostic that the pre-refactor
`World::calcSpatialForces()` provided (`Robosample/src/World.cpp:4355`), streaming it to a CSV at the
trajectory-write cadence for later analysis and VMD arrow rendering.

The recorded quantity, per interesting body, is a SINGLE spatial force `(force, torque)` about the
body origin `Bo`, in Ground: the SUM of whichever term(s) the reporter world was configured to include
(§1.2). Two terms are available:

- **`bodyForceG`** (Terms) — the OpenMM net applied force. Included by default.
- **the static (`u=0`) mobilizer reaction** (Terms) — the constraint force transmitted through the
  body's inboard joint. OPTIONAL, off by default.

Both terms are spatial forces about the SAME point (`Bo`, in Ground), so summing them is a valid
spatial-force addition — which term(s) are enabled is a per-reporter-world CHOICE, not a change to the
CSV's shape: the output is ALWAYS one `(force, torque)` pair per interesting body, per frame (§4).
Recorded for a sparse set of interesting bodies (conceptually: coarse domain-boundary bodies such as TM
helices), at the DCD write frequency, appended to a `.csv` and rendered by `vmd/arrows.tcl`.

Intended use: correlate the per-body force with geometric order parameters (e.g. TM-helix-to-membrane
angle) computed post-hoc from the paired DCD frame, and render force/torque arrows. See the GPCR
world-design specs (`docs/specs/gpcr-world-design/`) for which regions are worth reporting.

This spec covers only the **recording/export pipeline**. It does not add or alter any dynamics.

### Terms (defined once)

- **`bodyForceG`** — the net applied spatial force `(torque, force)` on a body, reduced from OpenMM's
  per-atom Cartesian forces about the body origin `Bo`, in Ground. Output of
  `ForceBridge::getForcesFromOpenMM` (`include/ForceBridge.hpp:56`): `force = Σ_a f_a`, `torque =
  Σ_a (r_a − Bo) × f_a`, summed over the body's atoms. This is the exact per-body reduction the
  articulated-body solver's dynamics consume every step (`docs/README` ForceBridge contract); the
  reporter reads it directly, with no forward-dynamics step. Velocity-free: a pure function of `q`.
- **Static (`u=0`) mobilizer reaction** — the spatial force `(torque, force)` a body's inboard
  mobilizer transmits to hold its outboard subtree in the field, reported at `Bo`, in Ground, with all
  generalized speeds set to zero. Output of `RobotEngine::calcMobilizerReactionForces`
  (`src/RobotEngine.cpp:1257`), ported from Simbody and oracle-validated (`docs/specs/
  robotics-oracle-reactions.md`, `tests/TestReactionForces.cpp`). Unlike `bodyForceG`, it is
  CUMULATIVE over the body's entire outboard subtree and depends on kinematic rooting — a different,
  complementary quantity, not a refinement of `bodyForceG`. Forced to `u=0` (hence also velocity-free)
  precisely so it can be summed with `bodyForceG` without mixing in a trajectory-endpoint-dependent
  signal (see "Why static", historically retained in `src/World.cpp`'s implementation comments).
- **Reporter world** — a world flagged to record forces. Its state supplies the `q` (and hence `X_GB`,
  atom positions, and both terms above) at which the snapshot is taken.
- **Interesting body** — a RobotModel body index whose force is recorded (Ground, body 0, is never
  included). SHOULD be sparse and at domain boundaries (§2.1).
- **Frame** — one trajectory-write event; a force CSV row and a DCD frame share a frame index and are
  joined on it.

## 1. What is recorded, and from which state

1. **Quantity — ONE `(force, torque)` pair per interesting body, per frame, in Ground.** The pair is
   the SUM of the enabled term(s) (§1.2): `bodyForceG` if `includeOpenmm`, the static `u=0` reaction if
   `includeReaction`, added together if both. No `u`/`uDot` output. No normalization at record time
   (§4).
2. **State — the accepted conformation.** Both terms SHALL be evaluated on the reporter world's
   accepted post-round `q`. `bodyForceG` is velocity-free by construction (a sum of per-atom OpenMM
   forces at fixed positions); the reaction term, when included, is forced to `u=0` for the SAME
   reason (§1.2) — so the recorded sum never reads a trajectory-endpoint velocity.
3. **State source SHALL be the reporter world, not `replicaCoords_`.** `Context::writeOutputs`
   (`src/Context.cpp:530`) writes the DCD frame from `replicaCoords_` (positions only). The recorded
   force needs the reporter world's articulated state (`X_GB`, hence atom positions in Ground), so it
   MUST be computed from that world on its accepted `q`. The accepted `q` used here SHALL be the same
   conformation written to the paired DCD frame, so force rows and DCD frames align 1:1.

### 1.1 Deltas from the old writer (deliberate)

- **Net applied per-body force is the default term, not a transmitted joint reaction.** The
  originally-shipped version of this reporter recorded ONLY the static (`u=0`) mobilizer reaction
  (`RobotEngine::calcMobilizerReactionForces`) — the spatial force an inboard mobilizer transmits to
  hold its outboard subtree in the field. That quantity folds in the *whole* outboard subtree's load
  through the kinematic chain. The structurally meaningful quantity for this diagnostic is the net
  EXTERNAL force the OpenMM field exerts directly on each rigid body — `bodyForceG` — so the reporter
  now includes that by default, with the static reaction available as an OPTIONAL additional term
  (§1.2) that SUMS into the same row rather than replacing it. With the default term selection
  (`includeOpenmm=true, includeReaction=false`) this also drops the per-snapshot forward-dynamics step
  (`calcUDot`) and the `u` save/zero/restore dance entirely: `bodyForceG` is read straight off
  `ForceBridge::evaluate`, a strict subset of what every dynamics step already does.
- **Per-body, not per-joint.** The old row was `(inboard_idx, outboard_idx, force, torque)` — an
  inboard/outboard atom PAIR bounding a joint. The new row is `(body_idx, atom_idx, force, torque)` —
  ONE body identity and ONE representative atom, since neither term has a notion of "outboard" atom in
  this reporter's row shape.
- **Streamed CSV, not accumulate-then-dump.** Rows are appended as frames are written, so a killed run
  keeps its data and large runs do not grow an unbounded in-memory history.
- **Raw values in the CSV.** The old `writeSpatialForces` normalized to global unit max before
  writing (lossy). The CSV SHALL hold raw values; normalization is a render-only transform in
  `arrows.tcl` (§5.2).
- **Geometry joined via the DCD frame, not duplicated.** The geometric CV is computed from the paired
  DCD frame in Python and joined on frame index; the CSV stores no poses.

### 1.2 Term selection: `includeOpenmm` / `includeReaction`

`enableReactionReporter`'s `includeOpenmm` (default `true`) and `includeReaction` (default `false`)
arguments select which term(s) are SUMMED into each interesting body's single `(force, torque)` row:

- `(true, false)` — DEFAULT. `bodyForceG` only. Reproduces the OpenMM-only output this reporter
  shipped with, byte-for-byte: no forward-dynamics step, `u` never touched.
- `(false, true)` — the static `u=0` reaction only.
- `(true, true)` — the SUM of both ("net" applied-plus-constraint load).
- `(false, false)` — INVALID: `captureReactionSnapshot` SHALL throw (nothing to report).

Both terms are spatial forces about the SAME point (`Bo`, in Ground, Terms), so their sum is a
well-defined spatial force — not an average, not a normalized combination, a plain vector sum. The CSV
format (§4) does NOT encode which term(s) went into a row: term selection is a property of the
reporter WORLD (fixed for that world's lifetime, set once via `enableReactionReporter`), not of the
row or the file.

When `includeReaction=true`, computing that term reintroduces, for that snapshot only, the static
`u=0` forward-dynamics evaluation the original writer used ("Why static": the accepted conformation is
a position sample with no intrinsic velocity, so the reaction is evaluated at `u=0` to get the
field-induced constraint load rather than a trajectory-endpoint-contaminated one): save `u`, zero it,
`realizeArticulatedBodyInertias`/`realizeVelocity`/`calcUDot` so `A_GB`/`zPlus` are consistent with
`u=0` and the field already evaluated for `bodyForceG`, then `RobotEngine::calcMobilizerReactionForces`
(one O(n) inward sweep), then restore `u` (and re-sync every `u`-derived cache) before returning — so
the no-perturbation property (§3) holds regardless of term selection.

## 2. Configuration surface

### 2.0 Consumer is `python/robosample/run.py`

`run.py` is the sole driver. A world in `run.py` is flagged as the reporter via a keyword on its
world-add call (the old name was `want_spatial_force_history`; keep it or rename, NOTE). The keyword
SHALL be restored on the world-add path `run.py` uses (`add_ncmc_world` / `add_robotic_world` /
`add_torsional_world` / `add_cartesian_world` in `python/robosample/context.py`, and the underlying
C++ binding in `src/PyBind11.cpp`, stub in `python/robosample/robo_bindings.pyi`). Legacy
`run_ffar1.py` is out of scope.

- The feature SHALL be **opt-in** and off by default; a run with no reporter world allocates nothing
  and does zero extra work (mirrors the ForceBridge "welded engine allocates nothing" discipline).
- Recording follows the DCD cadence: active only when a frame is written
  (`production && writeFreq > 0 && (round − equilRounds) % writeFreq == 0`, `src/Context.cpp:491`),
  one sample per such frame, per reporter world.

### 2.1 Interesting-body selection (conceptual guidance; NOTE)

Two consequences shape selection, though the engine does not enforce them:

- **Sparse, at domain boundaries.** A body at a domain's inboard base (plus its parent) reports that
  region's own loading; dense sampling gives redundant, locally near-identical readings for a rigid
  domain (all its welded-together bodies share ~the same net field load pattern up to how atoms are
  partitioned among them). Coarse is *more* representative.
- **Welded (rigid) bodies are the natural unit.** `bodyForceG` is a direct field reduction, meaningful
  for ANY body regardless of its joint's mobility. The static reaction, if included, is most
  interpretable at welded boundaries (a free mobilizer transmits ~zero reaction along its free axis).
  The reporter bodies SHOULD still bound rigid domains (e.g. welded TM-helix bodies) so each recorded
  value summarizes one coherent structural unit rather than a single atom's share of the load.

## 3. Snapshot computation

At the point the reporter world's accepted `q` is available (after its round resolves), for each
interesting body, the sample SHALL be produced without perturbing the sampler state:

1. Realize position on the reporter world's accepted state (`RobotEngine::realizePosition`) so
   `X_GB`/atom positions in Ground are current for this `q`.
2. Refresh the field: `bridge.evaluate(s)` — this pushes the accepted `q`'s atom positions to OpenMM,
   reads back per-atom forces, and reduces them into `bodyForceG` (`include/ForceBridge.hpp:112`,
   `getForcesFromOpenMM`). Velocity-independent. ALWAYS run, regardless of `includeOpenmm`:
   `bodyForceG` is also the `F_ext` the reaction term needs as a precondition (Terms).
3. OPTIONAL (`includeReaction=true` only, §1.2): save `u`, zero it, `realizeArticulatedBodyInertias` +
   `realizeVelocity` + `calcUDot` so `A_GB`/`zPlus` are consistent with `u=0` and the field already
   evaluated in step 2, then `RobotEngine::calcMobilizerReactionForces(m, s, reactionAtBo, nullptr)`
   (one O(n) inward sweep, at-Bo only), then restore `u` (and re-sync every `u`-derived cache) before
   proceeding. Skipped entirely when `includeReaction=false` — no forward dynamics, no `u` touched.
4. For each interesting body `b`, buffer the SUM of the enabled term(s) —
   `(includeOpenmm ? bodyForceG[b] : 0) + (includeReaction ? reactionAtBo[b] : 0)` — with `b`'s own
   body index and a representative atom identity (§4).
5. On the DCD-write rounds (§2.0), append the buffered rows to the `.csv` tagged with the frame index.

The procedure MUST NOT change the reporter world's `q`/`u` as seen by the next round. With
`includeReaction=false`, `bridge.evaluate` only writes arena scratch (`bodyForceG`/`atomForceG`/
`mobilityForce`), and every subsequent round's own `bridge.evaluate` call recomputes that scratch from
that round's `q` before it is read for dynamics — so no save/restore step is needed. With
`includeReaction=true`, step 3 DOES zero `u` for its forward-dynamics evaluation, but saves and
restores it (and every `u`-derived cache) before returning — so the no-perturbation property holds
regardless of term selection. This is the primary test (§6).

NOTE (integrator guard): a reporter world SHALL refuse (clear message) if its sampler is Cartesian-only
(`OpenMMVelocityVerlet`), where the reporter's body-index-based interesting-body set is not
meaningful (a Cartesian world has no articulated body indexing to select from).

NOTE (invalid selection): `captureReactionSnapshot` SHALL throw if BOTH `includeOpenmm` and
`includeReaction` are false (nothing to report) — the flags can still be left in that state via
`enableReactionReporter`, so the guard lives at snapshot time, not construction time.

## 4. CSV format (VMD-compatible)

- One CSV **per replica** (the filename carries the replica, matching the per-replica `.dcd`), so
  `vmd/arrows.tcl` (single-molecule) loads the file paired with that replica's trajectory. A
  `replica` column is still written (redundant with the filename, but keeps the row self-describing
  if files are ever concatenated for analysis).
- Comma-delimited, one row per (frame, interesting body). ALWAYS **10 columns**, regardless of term
  selection (§1.2):
  `frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz`.
- `body_idx` is the RobotModel body index. `atom_idx` SHALL be a **VMD/DCD atom index** (the prmtop
  atom order the DCD is written in, `systemTopology.atomsPrmtopIndex`) of a representative atom of
  that body — because `arrows.tcl` draws the force arrow at that atom's position.
- `fx fy fz` = raw SUMMED force (Ground); `tx ty tz` = raw SUMMED torque (Ground, about the body
  origin) — the sum of whichever term(s) the reporter world enabled (§1.2). No normalization,
  thresholding, or row-dropping at write time.
- Appended as frames are produced; a header comment line (`# frame,replica,body_idx,...`) MAY precede
  the rows (`arrows.tcl` skips `#` lines).

## 5. Export

### 5.1 CSV (primary)

The streamed per-replica `.csv` of §4 is the primary artifact; the GPCR analysis joins it to the DCD
frames by frame index.

### 5.2 VMD arrows (`vmd/arrows.tcl`)

- `arrows.tcl` renders red force arrows and blue torque arcs per frame from this file. It parses
  **comma**-delimited rows with the §4 columns (`frame,replica,body_idx,atom_idx,fx,fy,fz,tx,ty,tz`).
  The force arrow is drawn at `atom_idx`'s position; the torque arc is drawn about ITS OWN axis (the
  torque vector direction, `vecnorm(tx,ty,tz)`) centered at `atom_idx`'s position — there is no
  outboard atom to define a bond axis anymore (the quantity is per-body, not per-joint). Its
  global-max/manual-scale normalization is a render-only scaling and MUST NOT be pushed back into the
  CSV. `arrows.tcl` renders whatever SUM the CSV holds; it has no visibility into which term(s) went
  into that sum.

## 6. Test gate

- **No-perturbation (primary).** Recording a sample then advancing further rounds SHALL produce
  byte-identical sampler behavior to a run with recording disabled, for EVERY term selection. A
  minimal reproducer wired into the suite SHALL assert this. With `includeReaction=false`, `bodyForceG`
  is velocity-free and `bridge.evaluate` only overwrites arena scratch that every round's own
  `evaluate` call recomputes fresh (§3), so this reduces to checking that `captureReactionSnapshot`
  never touches `q`/`u`. With `includeReaction=true`, `u` IS zeroed for the reaction's forward-dynamics
  step (§1.2/§3 step 3), so the test SHALL additionally assert `u` is bit-identical to its
  pre-snapshot value after the call.
- **Cadence + format.** With a reporter world enabled, a short production run SHALL emit exactly one
  CSV row per interesting body per DCD frame, comma-delimited with the fixed 10 columns (§4) and raw
  values, and zero rows during equilibration; frame indices align 1:1 with DCD frames, for EVERY term
  selection.
- **Correctness (OpenMM term).** With `includeReaction=false`, the recorded `(fx,fy,fz,tx,ty,tz)` SHALL
  equal a direct `ForceBridge::evaluate` + `RobotState::bodyForceG()[b]` read on the same accepted `q`.
- **Correctness (reaction term).** With `includeOpenmm=false, includeReaction=true`, the recorded
  `(fx,fy,fz,tx,ty,tz)` SHALL equal a direct `u=0` `RobotEngine::calcMobilizerReactionForces` call on
  the same accepted `q` (reuses the oracle of `tests/TestReactionForces.cpp`).
- **Correctness (sum).** With `includeOpenmm=true, includeReaction=true`, the recorded
  `(fx,fy,fz,tx,ty,tz)` SHALL equal the elementwise sum of the two single-term reads above, on the
  same accepted `q`.
- **Invalid-selection refusal.** `includeOpenmm=false, includeReaction=false` SHALL cause
  `captureReactionSnapshot` to throw with a clear message.
- **Cartesian refusal.** A reporter flag on a Cartesian-only sampler SHALL be refused with a clear
  message.
- The full authoritative gate (`nox -s tests`) SHALL pass with the feature disabled (default) and
  enabled under every valid term selection.

## 7. Non-goals

- No per-substep or trajectory-interior sampling; the sample is on the accepted conformation at the
  DCD cadence.
- No 16-column / split-term CSV output; the row is ALWAYS the fixed 10-column SUM (§4). A caller
  wanting the two terms separately runs the reporter twice (once per term selection) into two
  differently-named runs/files -- this reporter does not support recording both a single-term AND a
  summed row in the same pass.
- No transmitted-joint-reaction quantity BY DEFAULT (`includeReaction=false`); `RobotEngine::
  calcMobilizerReactionForces` and its oracle (`docs/specs/robotics-oracle-reactions.md`,
  `tests/TestReactionForces.cpp`) remain in the engine for other uses, and are OPTIONALLY summed in by
  this reporter (§1.2) when `includeReaction=true` is requested explicitly.
- No `u`/`uDot` output; both terms are velocity-free by construction (the reaction term by being
  forced to `u=0`).
- No collective-variable computation in C++ (geometry lives in Python, joined by frame index).
- No change to dynamics, acceptance, or any sampled quantity.
