# Spec: Replica Exchange, Nonequilibrium WORK, and REBAS — port and cleanup

Status: **revision 2** — reviewed (4-agent hostile review, 2026-07-12) and revised to resolve the
confirmed Blocking findings. Ready for coder. Author: derived from a full read of the original
SimTK-molmodel/Simbody Robosample under `./Robosample/`. Target: port the subsystem into
the current engine (parent repo, `disasm` branch), removing dead/stubbed paths and fixing
the identified defects.

Requirement keywords follow RFC 2119 (`styles/spec.md`): SHALL / SHOULD / MAY / NOTE.

**Revision 2 changelog (what the review changed).** The two highest-risk physics claims held under
independent re-derivation (D6 corrected Jacobian; D3/INV-7 Fixman biconditional). Four Blocking findings
required revision:

- **D7 (new), INV-8 rewritten, RQ-1 withdrawn.** The "measure `U` at `x'`" fix was insufficient: the
  driven MD runs at the *source* temperature, so post-scale MD injects uncancelled heat. Driven worlds
  now run `mdSteps = 0` (exact deterministic-map limit; Nilmeier eq:20). Target-temperature NCMC is the
  deferred performance path.
- **D2 corrected, INV-9 (new), V10 (new).** The affine scaling map is an involution — hence
  `correctionTerm = 1` — only if forward and reverse share one frozen anchor mean; the original uses
  per-state means (`μ_C ≠ μ_H`), which breaks it. The scaling anchor moves to shared/global ownership.
- **Port-target baseline (new), citation-tree convention (new), Interface I1–I5 re-baselined.** The
  spec's "current engine" citations pointed at the *original* tree; the target already ships a
  coordinate-swap `runREX`. The port is a from-scratch reimplementation, guided by the original.
- **F4 upgraded High→Critical** (two bugs: stale-WORK REMC revert + non-atomic energy commit);
  **F5 downgraded High→Nit** (the `rewindReplica` throw is unreachable dead code); **F13 re-described.**
- **INV-10 (new)** ties RENE (BAT scaling, with Jacobian) vs RENEMC (volume-preserving velocity drive,
  no Jacobian) to the `distortOption` sign, resolving the INV-1 tension. **D5 selection rule** made
  concrete per `JointType`. **V8 corrected** (blind to the Jacobian; V5 is the sole Jacobian guard),
  **V5 strengthened** (≥3-body chain), **V9 (new)** fails loud on Fixman-off. B7 parity and B12
  accumulation fixed. Reproducers delivered: `tests/bat_jacobian_scaling_check.py` (V5),
  `tests/test_rex_swap_acceptance_algebra.py` (V8).

Primary citations use `file:line` into `Robosample/src` and `Robosample/include` (the original).
Literature keys reference `references/index.yaml`.

**Citation-tree convention (added in revision 2, 2026-07-12).** Unless a citation is explicitly tagged
`[CURRENT]`, every `file:line` in this spec points into the **original** tree `Robosample/src` — it
describes the *design being ported*, not code that exists in the target. The current engine (`src/`,
`include/`, `python/robosample/`) has none of `Replica`, `ThermodynamicState`, `attemptREXSwap`, the
two inverse maps, `WORK`, or `REBAS`; grep confirms zero occurrences of `RUN_TYPE`, `REBASONTOP`,
`RENEMC`, `zMatrixTable`, or the removed-symbol list. **The port is a from-scratch reimplementation of
the subsystem in the current engine, guided by the original's design — not an in-place uncommenting.**
See "Port-target baseline" below for what actually exists in the target.

---

## Motivation

### Problem restatement (user's phrasing)

The user asked, verbatim in substance:

- Document the run types `Default, REMC, RENEMC, RENE, REBASONTOP`; Python currently exposes
  only `DEFAULT/REMC/RENEMC/RENE` (`PyBind11.cpp:412-416`) — expose `REBASONTOP` too, with
  documentation.
- Document the Q-scale-factors and the distort options (introduce an enum instead of magic ints).
- Ask whether `RunReplicaWorldRange()` and the neighbouring functions can be simplified.
- Explain why there are special handles for equilibrium vs work, why `WORK` exists at all across
  `Context` (the run) and `Replica`, and what the "Q statistics" are.
- Document, in detail, how a replica swap is performed.
- Decide where coordinates and Z-matrices should live — replicas or worlds.
- State how many worlds/replicas we store: one per robot and one per thermodynamic state, or
  `W` worlds and `W*T` thermodynamic states.
- Document `prepareExchangePairs` and the replica-exchange apparatus in general.
- Wire `replicaMixingScheme` and `prepareExchangePairs` into real code (currently undocumented and
  disconnected).
- Document `attemptREXSwap` properly.
- State clearly, in code, that Fixman never enters acceptance.
- Confirm whether `REBASONTOP` is identical to `RENE`.
- Remove the hardcoded molecule and testing modes in REBAS, but preserve that machinery as a
  spec that later becomes a test.
- Propose fixes for everything dead/stubbed.

### The subsystem in the codebase's vocabulary

Robosample runs blocked Gibbs sampling: each **world** is one robot factorization (a Gibbs block)
that exposes a subset of generalized (BAT) coordinates and evolves them with Hamiltonian Monte
Carlo. Replica exchange wraps that inner sampler in an outer Markov chain over temperature.

The subsystem has three object families (`Context.cpp`):

- **World** — a physical Gibbs-block integration engine. There are `W` of them
  (`addWorld`, `Context.cpp:493-518`), held in one shared `std::vector<World> worlds`. Worlds are
  transient scratch: within a round each replica loads its coordinates into the worlds, integrates,
  and reads them back. Worlds hold no persistent per-replica state between rounds.
- **Replica** — a persistent molecular configuration plus its energies
  (`addReplica`, `Context.cpp:870-910`). There are `R` of them in `std::vector<Replica> replicas`.
- **ThermodynamicState** — a temperature and a per-world simulation schedule
  (`addThermodynamicState`, `Context.cpp:912-948`). There are `T` of them.

Two mutually-inverse maps connect replicas to states: `replica2ThermoIxs` and `thermo2ReplicaIxs`
(`Context.hpp:408,412`), initialised to identity by `loadReplica2ThermoIxs` (`Context.cpp:983-998`).

### The single most important fact

**The nonequilibrium half of this subsystem is checked in disabled.** As shipped:

- `RunREX` runs each replica over its equilibrium worlds only and performs **zero swaps**: every
  call to `prepareExchangePairs` is commented out (`Context.cpp:2692,2704,2745`), so
  `exchangePairList` is empty when `mixReplicas` runs (`Context.cpp:2801`).
- The two nonequilibrium propagation blocks ("RUN A / RUN B") are commented out
  (`Context.cpp:2701-2799`), so the `WORK_*` buffers are never populated at runtime.
- `replicaMixingScheme` is never read; `mixAllReplicas` is never called.
- The REBAS proposal-correction term is hardcoded to 1 (`Context.cpp:1266-1267`).
- Fixman is computed but never enters acceptance (`Context.cpp:1211-1234,1244,1248`).

The port therefore reconstructs an intended design that exists only as scaffolding. This spec
distinguishes, for every element, **observed live behavior**, **intended behavior inferable from
commented code**, and **proposed target behavior**.

NOTE (user directive, 2026-07-12): the nonequilibrium code "was commented for other reasons" and
**SHALL be reachable through Robosample**. The port therefore reinstates the commented driving loops,
swap wiring, and `REBASONTOP` interleaving (F1, F2, D4) and exposes all four exchange run types — it
does not treat "commented" as "intended to be absent". Where a choice arises, the port stays **as close
to the live original as possible** (see Resolved decisions D1–D6).

**Fidelity principle.** "As close to the original as possible" binds the code that **actually ran** in the
original — REMC, the object model, the label-swap, the two-endpoint acceptance structure, the WORK/REBAS
algorithms. The **dormant scaffolding** (the commented-out nonequilibrium driving, swap wiring, and
`REBASONTOP` interleaving) never executed, so reinstating it is a fresh, correct implementation, not a
change to working behaviour: where the researcher found the dormant code carries a bug that would bias the
sampled distribution (the BAT-Jacobian composition, D6/F7; the work-config timing, INV-8), the port
implements it **correctly** per theory and flags the deviation. The user did not write this code and asks
that bugs be flagged, not preserved. Faithfulness is to the original's *design intent*, not to bugs in code
that never ran.

### Port-target baseline (current engine) `[CURRENT]`

The target already ships a working replica exchange, which this spec REPLACES/EXTENDS. Confirmed by grep
against the current tree:

- `Context::runREX` (`src/Context.cpp:396-503` `[CURRENT]`) runs a **coordinate-swap** REMC: it holds
  `replicaCoords_`, one temperature per replica in `temperatures_` (`include/Context.hpp:124-126`
  `[CURRENT]`), evaluates each replica serially through the OpenMM singleton
  (`openmmPotential(replicaCoords_[r])`, `:494`), accepts on `(β_A−β_B)(E_a−E_b)`, and on accept does
  `std::swap(replicaCoords_[r], replicaCoords_[r+1])` (`:500`). There is no run-type argument
  (`include/Context.hpp:81` `[CURRENT]`), no `Replica`/`ThermodynamicState` abstraction, no label-swap.
- `DistortOption` is **already an `enum class`** (`include/World.hpp:59` `[CURRENT]`, current values
  `{NMA=0}`), carried as `std::optional<DistortOption>` on `SamplerConfig` (`:112`) and as a typed
  `add_sampler` parameter (`:267`). There is no `-6..-1 / 0 / 1..6` magic integer to migrate away from.
- `World` fixes timestep, MD steps, and accept/reject mode at `add_sampler` time; `SamplerConfig
  sampler_` is private (`include/World.hpp:515,578` `[CURRENT]`) with **no** `setTimeStep`/`setMdSteps`/
  `setAcceptRejectMode` and no mutable accessor. `setTemperature` **does** exist (`include/World.hpp:390`
  `[CURRENT]`) and is already driven per-replica each sweep (`src/Context.cpp:458`).
- The shared molecular connectivity Z-matrix lives on `SystemTopology`
  (`include/TopologyElements.hpp:48-49,137-138` `[CURRENT]`); there is no per-replica BAT-value storage
  and no per-state BAT-statistics object.
- Per-body joint identity and DOF counts are available via `RobotModel::bodyJoint` (`JointType`),
  `bodyNU`, reachable through `World::model()` (`include/RobotModel.hpp:68,78`; `include/World.hpp:303`
  `[CURRENT]`) — this is the API the topology-driven REBAS selection (D5) uses.

**Reconciliation (label-swap vs the existing coordinate-swap).** The port replaces the current
coordinate-swap `runREX` with the original's label-swap `Replica`/`ThermodynamicState` model (B6). The
two are provably equivalent for pure-temperature REMC (same permutation distribution over states); the
existing coordinate-swap `runREX` is retained as an **equivalence oracle** (INVARIANT-EQUIV: a
label-swap REMC run and the legacy coordinate-swap run SHALL produce the same per-state energy
statistics). Every Interface item (I1–I5) and object-model claim (B0/B1) below describes the *original*
design to be reimplemented unless tagged `[CURRENT]`; any handoff item resolving to `Robosample/` is a
from-scratch add, not an edit-in-place.

### Binding constraints

- Target scale: up to 1M atoms, 100k rigid bodies, 10k robots (`CLAUDE.md`). The storage model and
  the swap loop SHALL be `O(R)` per round in replica count and SHALL NOT duplicate per-replica
  Cartesian state across worlds.
- Correctness precedence (`CLAUDE.md` decision policy): published theory > project specs > tests >
  implementation > intuition. Where the original code contradicts theory (e.g. the missing
  proposal correction), theory wins.

---

## Behavior

### B0. Object model and counts (answers "how many worlds/replicas")

Live wiring, confirmed from `python/robosample/context.py:1061-1143`:

- Worlds are added once. Each world carries exactly one sampler with a `distortOption`, `flow`,
  `integratorType`, `timeStep`, `mdSteps`, `acceptRejectMode` (`context.py:1073-1088`).
- Every thermodynamic state is created with `worldIndexes = [0, 1, ..., W-1]` — the **same** full
  ordered world list (`context.py:1069,1086,1140`).
- One `addReplica()` and one `addThermodynamicState(temp, ...)` are created **per temperature**
  (`context.py:1130-1143`). Hence **`R == T == number of temperatures`**.

Therefore the model is **`W` shared worlds, `R = T` replicas/states — not `W*T`**. Each state
schedules the same `W` worlds at its own temperature. The per-world distort/flow/integrator vectors
are identical across states (extracted once from the worlds, `context.py:1062-1088`) and copied into
every `ThermodynamicState`.

NOTE (consequence): because the `W` worlds are shared scratch, replicas are integrated **serially**
within a round — a world cannot hold two replicas' states at once. Parallel replica execution is not
possible without either replicating the worlds per replica or making worlds re-entrant. The port keeps
the serial model (resolved in D1: parallel replicas are feasible only in multi-GPU deployments, out of
scope here).

The spec adopts these definitions:

- `W` = number of worlds = number of Gibbs blocks / robot factorizations in one sweep.
- `R` = number of replicas = number of persistent configurations.
- `T` = number of thermodynamic states = number of temperatures. Canonically `R = T`.

### B1. Where state lives (answers "replicas or worlds?")

Target ownership SHALL be:

| Data | Owner | Rationale |
| --- | --- | --- |
| Cartesian coordinates (committed) | **Replica** (`atomsLocations`) | Persistent across rounds; worlds are shared scratch and cannot retain it. |
| Cartesian coordinates (nonequilibrium trial) | **Replica** (`WORK_atomsLocations`) | A proposal held apart from the committed state until accept/reject. |
| Flattened x/y/z write buffers | **Replica** | Derived from `atomsLocations`; used for DCD output. |
| Z-matrix **table** (connectivity/topology) | **Shared / global** (`Context::zMatrixTable`) | Identical for every replica and world — it is molecular topology, not configuration. Passed by reference (`Context.cpp:876-880,923-929`). |
| Z-matrix **BAT values** (per configuration) | **Replica** (derived) | A function of the replica's Cartesian coordinates; computed on demand or cached on the replica. Never on a world. |
| BAT running statistics (means/vars) | **ThermodynamicState** | Accumulated over the ensemble at one temperature. |
| Transient integrator state (Simbody `State`, velocities) | **World** | Valid only during that world's integration. |

This matches where the original already puts most things; the port SHALL make it explicit and
SHALL remove the empty `Replica::restoreCoordinates`/`storeCoordinates` stubs (`Replica.cpp:115-120`)
in favour of the `Context` transfer helpers (`transferCoordsFromReplicaToWorld`,
`transferCoordsFromWorldToReplica`, `Context.cpp:13-36`).

### B2. Why `WORK` exists, and the equilibrium/work duality (answers "why special handles")

`WORK` implements **nonequilibrium candidate Monte Carlo** for the swap move. Plain
parallel-tempering (REMC) swaps two equilibrium configurations at different temperatures and accepts
on `-Δβ ΔU`; its acceptance collapses when the two configurational ensembles overlap poorly, which is
the regime of large systems. The nonequilibrium variants instead **drive** each replica with a short
work-generating perturbation that scales its BAT coordinates toward the neighbour's temperature, then
accept on the accumulated work (Crooks/Jarzynski / NCMC form). See Derivation sketch.

The driven trajectory is a **proposal**. It SHALL be kept separate from the committed state until the
swap decision, which is exactly why `Replica` carries two of everything (`Replica.hpp:166-181`):

| Committed (equilibrium) | Trial (`WORK_`, nonequilibrium) |
| --- | --- |
| `atomsLocations` | `WORK_atomsLocations` |
| `potential`, `referencePotential` | `WORK_potential`, `referenceWORK_potential` |
| `FixmanPotential` | `WORK_FixmanPotential` |
| — | `WORK`, `workJacobiansContributions` (accumulators) |

On **accept**, the trial is promoted to committed:
`updAtomsLocationsInGround_FromWORK()` (`Replica.cpp:50-74`) and `setPotentialEnergy_FromWORK()`
(`Replica.cpp:110-112`), invoked via `Context::set_WORK_CoordinatesAsFinal` /
`set_WORK_PotentialAsFinal` (`Context.cpp:1596-1604`). On **reject**, the trial is discarded and the
committed state stands.

`WORK` lives at both the `Replica` level (storage of the per-replica accumulated work) and the
`Context` level (the accumulation loop and the swap consumer) because accumulation happens across the
worlds of one replica's sweep (`Context`) while the result is owned by the replica.

### B3. What the Q statistics are (answers "what are those q statistics")

`ThermodynamicState::calcQStats` (`ThermodynamicState.cpp:402-592`) maintains single-pass running
means and variances, per world and per generalized coordinate, of three families:

- **Qs** — the generalized coordinates: `Qmeans`, `Qdiffs`, `Qvars`.
- **BMps** — bond-related per-coordinate terms (inferred: bond mid-point / boost-momentum scale):
  `BMps_means`, `BMps_diffs`.
- **PFrs** — pin-frame / angle-related per-coordinate terms: `PFrs_means`, `PFrs_diffs`.

They are updated only from equilibrium worlds (`Context.cpp:2477-2492`). Their **sole live consumer**
is the BAT-scaling perturbation (B4), which displaces each scaled coordinate by
`(scaleFactor - 1) * (current_BAT_value - running_mean)` (`HMCSampler.cpp:570,584-586,597,606`). The
statistics therefore provide the **reference point** the distortion scales around; without them the
drive has no anchor. `calcZMatrixBATStats` (`ThermodynamicState.cpp:265-318`) is the BAT-coordinate
analogue.

NOTE (confirmed by researcher against `World::getBMps`/`getPFrs`, `World.cpp:3110-3153`): `BMps` and
`PFrs` are per-body **generalized-coordinate reference values**. `BMps[mbx] = X_BM.p()[0]` is the
x-component of the outboard (Body→Mobilizer, "B-M") frame origin — the **bond-length coordinate** (used at
`HMCSampler.cpp:570` as `BONDLengths[row] − BMps_mean`). `PFrs[mbx] = acos(X_PF.R()(0)(0))` is the
inboard (Parent→Fixed, "P-F") frame angle — the **bond-angle coordinate** (used at `:597` as
`ANGLEBends[row] − (π − PFrs_mean)`). The earlier "bond mid-point / boost-momentum" reading for `BMps`
was incorrect. These running means are the reference the affine drive scales the deviation around.

### B4. The distort mechanism and Q-scale-factors (answers "document Q-scale-factors, distort options")

**Distort option** is a per-world integer (`HMCSampler.hpp:632`, set from the schedule at
`Context.cpp:1858-1859`) whose sign selects the machinery:

- `== 0`: equilibrium HMC (no drive).
- `> 0`: velocity/NMA distortion — `setVelocitiesToNMA` (`HMCSampler.cpp:1045-1520`).
- `< 0`: position (BAT bond/angle) distortion — `positionsPerturbMethod()` maps `-1..-6` to
  `BendStretch1..6` (`HMCSampler.cpp:3326-3350`), applied by `perturbPositions` before the HMC MD
  segment (`HMCSampler.cpp:3557-3567`).

The port SHALL replace the magic integer with an enum (Interface I2). Proposed semantics:

```
enum class DistortOption {
    None                = 0,   // equilibrium HMC
    ScaleBendStretch    = -6,  // live BAT bond/angle scaling (was BendStretch6)
    // NMA velocity variants (prototype; port only if velocity-drive is in scope)
    NMAAlternateSign    = 1,
    NMARandomSign       = 2,
    // ... 3..6 as needed
};
```

**Q-scale-factor** `s` is the amount by which a driven world scales BAT coordinates.
`PrepareNonEquilibriumParams_Q` (`Context.cpp:1694-1746`) computes, for neighbouring state pairs,

```
s_i = sqrt(T_j / T_i)
```

with even/odd variants (`qScaleFactorsEven/Odd`, `Context.cpp:1706-1736`). This is the
thermodynamic-length-inspired scaling: to nudge replica `i` toward neighbour `j`'s temperature,
stretch its bonds/angles by `sqrt(T_j/T_i)`. `perturbScalingFactor` (`Context.cpp:1893-1940`) may
further randomise `s` (`deterministic` / `Gauss` / `uniform` / `Bernoulli`). The scalar reaches the
sampler through `setBendStretchStdevScaleFactor` (`HMCSampler.cpp:4169-4170`) as `QScaleFactor`
(default 1.0, `HMCSampler.hpp:682`).

In the live path (`HMCSampler.cpp:516-663`, `testingMode == false`), only `BendStretch6` performs
real work; `BendStretch1..4` write precomputed deltas and `BendStretch5` only prints. For each
mobilized body selected for scaling, each Q component is displaced by
`(scaleFactor - 1) * (current - running_mean)` and the scaling log-Jacobian accumulates in `J_scale`.

### B5. Accumulating WORK and its Jacobian

Per driven world, work is `W_world = (U_curr - U_prev) + (F_fixman,curr - F_fixman,prev)`, counted
only when the world's distort option is negative (`World::getWork`, `World.cpp:1915-1935`). The
Context sums it into the replica when the world is nonequilibrium
(`transferCoordsFromWorldToReplica(..., intoWORK=true)`, `Context.cpp:18-35`):

```
destReplica.updWORK()          += srcWorld.getWork();
destReplica.upd_WORK_Jacobian()+= srcWorld.getSampler(0)->getDistortJacobianDetLog();
```

Accumulators are zeroed per driven range (`Context.cpp:2425-2426`). The original computes
`bendStretchJacobianDetLog = J_ini + J_scale - J_fin` (`HMCSampler.cpp:659`), where `J(x) = Σ_bodies
[2 ln r + ln sin θ]` (`calcMobodsBATJacobianDetLog_NEW`, `HMCSampler.cpp:5147-5226`) is the BAT→Cartesian
volume element.

**CRITICAL (researcher REFUTE, D6): this composition is wrong and would ship a biased Jacobian.** The
scaling map factors `x^0 →(Cart→BAT) q^0 →(scale) q^τ →(BAT→Cart) x'`, so the forward Cartesian
log-Jacobian is `lnJac = (J(x') − J(x^0)) + N_scaled·ln s`. The code inverts the `(J_fin − J_ini)` sign
and its `J_scale` (value-ratios `log(r'/r)+log(θ'/θ)` over every Q, bond double-counted at `:590`/`:624`)
is not the scaling Jacobian. Numeric check — a single isotropic bond scale `r→s·r` has Cartesian
`|∂x'/∂x| = s³` so `lnJac = 3 ln s`; the correct formula gives `2 ln s + 1·ln s = 3 ln s`, the code gives
`≈0` (double-counted) or `−ln s` (F8-fixed). See D6 for the corrected prescription. The `2 ln r + ln sin θ`
element itself is CONFIRMED correct (torsion contributes nothing to the radial/angular volume element).
The two other conventions are deleted: `4 ln r + 2 ln sin θ` (`Replica.cpp:565-586`; `Context.cpp:3634-3655`,
which is `J²`, dead, no-op NaN guard) and `internNdofs·log s` (`HMCSampler.cpp:4399-4494`).

### B6. The replica swap in detail (answers "how is replica swap performed", documents `attemptREXSwap`)

`attemptREXSwap(thermoState_C, thermoState_H)` (`Context.cpp:1133-1420`). Nomenclature: `C` = cold,
`H` = hot; `X = thermo2ReplicaIxs[C]`, `Y = thermo2ReplicaIxs[H]` (`Context.cpp:1141-1142`). Suffix
`set` = equilibrium endpoint (X⁰); suffix `tau` = end of the nonequilibrium drive (Xτ).

1. **Record the attempt** in the symmetric `nofAttemptedSwapsMatrix` (`Context.cpp:1145-1146`).
2. **Gather energies** for both replicas (`Context.cpp:1149-1209`): `beta_C/beta_H`; committed and
   reference potentials (`U_*set`, `refU_*set`); driven-endpoint potentials (`U_*tau`, `refU_*tau`);
   accumulated work `W_* = getWORK()`; log-Jacobians `lnJac_* = get_WORK_Jacobian()`; Fixman terms.
   Reduce each energy by the relevant β.
3. **Assemble the three log-terms** (`Context.cpp:1237-1252`):
   - `ETerm_equal    = -[(β_H-β_C)(refU_Xset - refU_Yset)]`  — standard PT exponent `-Δβ ΔU`.
   - `ETerm_nonequil` — same form on the driven-endpoint (`tau`) reference potentials.
   - `Work_X = (ref_uH_Xtau - ref_uC_Xset) - lnJac_X`, `Work_Y = (ref_uC_Ytau - ref_uH_Yset) - lnJac_Y`,
     `WTerm = -(Work_X + Work_Y)` — the Crooks/NCMC work term.
4. **Proposal correction** `correctionTerm` (`Context.cpp:1254-1267`) — intended as the Hastings ratio
   of scale-factor proposal densities; **hardcoded to 1** as shipped.
5. **Select the acceptance exponent by run type** (`Context.cpp:1322-1331`):

   | RUN_TYPE | `log_p_accept` |
   | --- | --- |
   | `REMC` | `ETerm_equal` |
   | `RENEMC` | `ETerm_nonequil + log(correctionTerm)` |
   | `RENE`, `REBASONTOP` | `WTerm + log(correctionTerm)` |

6. **Accept/reject**: accept iff `log_p_accept >= 0 || U(0,1) < exp(log_p_accept)`
   (`Context.cpp:1356`).
7. **On accept** (`Context.cpp:1357-1405`): increment sample counters; **commit** the driven state of
   both replicas (`set_WORK_CoordinatesAsFinal`, `set_WORK_PotentialAsFinal`); **swap labels** via
   `swapThermodynamicStates` (below). On **reject**: nothing (`Context.cpp:1406-1418`).

`swapThermodynamicStates(replica_i, replica_j)` (`Context.cpp:1078-1099`) is a **label swap, not a
coordinate swap**:

```
swap(replica2ThermoIxs[replica_i], replica2ThermoIxs[replica_j]);
swap(thermo2ReplicaIxs[thermoState_i], thermo2ReplicaIxs[thermoState_j]);
thermodynamicStates[thermoState_i].setZMatrixBATPointer(replicas[replica_j].getZMatrixBATPointer());
```

Configurations never move; each replica keeps its coordinates and is simulated at a different
temperature next round. This is the Chodera–Shirts Gibbs-sampling-over-permutations formulation
(`chodera_2011_gibbs_replica_exchange`, cited in-code at `Context.cpp:1126-1130`).

### B7. The exchange apparatus and mixing (documents `prepareExchangePairs`; wires `replicaMixingScheme`)

`ReplicaMixingScheme { All = 0, Neighboring = 1 }` (`bgeneral.hpp:81-84`).

- **Neighboring** — `prepareExchangePairs(rexRound, oddity)` (`Context.cpp:1458-1480`) builds
  nearest-neighbour pairs with an alternating offset:
  ```
  startIdx = (rexRound + oddity) % 2;
  for (thIx = startIdx; thIx + 1 < K; thIx += 2)
      exchangePairList.emplace_back(thIx, thIx + 1);
  ```
  i.e. even rounds pair `(0,1),(2,3),…`; odd rounds pair `(1,2),(3,4),…`. `mixReplicas`
  (`Context.cpp:1482-1497`, gated by `swapEvery`) iterates `exchangePairList` and calls
  `attemptREXSwap` on each pair.
- **All** — `mixAllReplicas(nSwapAttempts)` (`Context.cpp:1439-1456`) draws `nSwapAttempts` random
  state pairs and attempts each.

The port SHALL **wire the scheme selection** (currently `replicaMixingScheme` is never read). Target
`mixReplicas`:

```
if ((mixi % swapEvery) != 0) return;
if (runType == RUN_TYPE::Default || nofThermodynamicStates <= 1) return;
if (replicaMixingScheme == ReplicaMixingScheme::Neighboring) {
    prepareExchangePairs(exchangeRound, /*oddity=*/0);   // parity from the EXCHANGE-round counter, not mixi
    for (auto [i, j] : exchangePairList) attemptREXSwap(i, j);
    ++exchangeRound;                                     // one increment per EXECUTED mix
} else { // All
    mixAllReplicas(nSwapAttempts);                       // nSwapAttempts SHALL be configurable
}
```

**Parity SHALL derive from a dedicated exchange-round counter, not from `mixi % 2` (revision 2, reviewer
Should-fix R5).** `prepareExchangePairs` sets `startIdx = (round + oddity) % 2`, so even/odd rounds
alternate the neighbour pairing `(0,1),(2,3),…` ↔ `(1,2),(3,4),…`. If parity is taken from the raw
`mixi` *after* the `mixi % swapEvery == 0` gate and `swapEvery` is even, every surviving `mixi` shares
the same parity → `startIdx` is frozen → one pairing class (e.g. `(1,2),(3,4)`) is **never** attempted,
breaking exchange ergodicity. Incrementing a separate `exchangeRound` once per *executed* mix guarantees
both parities cycle regardless of `swapEvery`. NOTE: the original `mixReplicas` guards on `nofReplicas <=
1` (`Context.cpp:1489`) while the target above writes `nofThermodynamicStates <= 1`; these are equal only
under `R == T` (B0), which the port SHALL assert, not assume.

`RunREX` SHALL call `mixReplicas` (which calls `prepareExchangePairs`) **live** — the original's
commented-out call sites (`Context.cpp:2692,2704,2745`) are the root cause of "no swaps happen". The
legacy `setReplicaExchangePairs` (`Context.cpp:1658-1683`) has no place in the port.

### B8. Equilibrium vs nonequilibrium propagation and the Partitioning

`computeNonequilPartitioning` (`ThermodynamicState.cpp:129-159`) scans the per-world distort options
and splits the schedule: worlds `[0, N1_wCnt)` are equilibrium, `[N1_wCnt, W)` are nonequilibrium,
where `N1_wCnt` is the first world with a nonzero distort option, `N2_wCnt = N1_wCnt - 1` is the last
equilibrium world, and `nofEquilibriumWorlds`/`nofNonequilibriumWorlds` count each segment. The
intended round runs the equilibrium segment, then the nonequilibrium segment (which populates
`WORK_*`), then attempts swaps. The port SHALL implement the nonequilibrium segment
(`RunReplicaWorldRange(replicaIx, N1_wCnt, nofNonequilibriumWorlds, shouldPrint)`, using the
position-clean signature of B12 — the `isNonEquilibrium` parameter is removed), guarded by
`nofNonequilibriumWorlds > 0`. Per D7, every world in this segment SHALL be configured `mdSteps = 0`.

### B9. Run types (answers "document run types"; expose REBASONTOP)

| RUN_TYPE | Meaning | Acceptance | Drives? |
| --- | --- | --- | --- |
| `Default` | No exchange; independent replicas. | — | No |
| `REMC` | Replica Exchange Monte Carlo (parallel tempering). | `-Δβ ΔU` on committed reference potentials. | No |
| `RENEMC` | Replica Exchange Non-Equilibrium MC. | `-Δβ ΔU` on driven-endpoint potentials (no Jacobian). | Yes — **volume-preserving** (velocity/NMA, `distortOption > 0`) |
| `RENE` | Replica Exchange Non-Equilibrium (work-based). | `-(W_X + W_Y)` (Ballard–Jarzynski / NCMC), includes `lnJac`. | Yes — **volume-changing** (BAT scaling, `distortOption < 0`) |
| `REBASONTOP` | RENE with interleaved REMC neighbour swaps "on top". | Same as `RENE` for the work swap; plus a REMC sub-loop. | Yes — BAT scaling |

`REBASONTOP` SHALL be exposed to Python with documentation (Interface I1).

**RENEMC vs RENE — the drive-sign biconditional (revision 2, resolves the INV-1 tension).** The two
non-equilibrium types are distinguished by the **sign** of their world's `distortOption` (B4), which
fixes whether the configurational Jacobian is unity:

- **RENE / REBASONTOP** drive with **position (BAT) scaling** (`distortOption < 0`). The scaling map is
  volume-changing (`|∂x'/∂x| = s^{N_scaled·…} ≠ 1`), so the acceptance SHALL carry `lnJac` — this is
  `WTerm` (B6, D6). Omitting `lnJac` here is the F7 bias.
- **RENEMC** drives with **velocity / NMA distortion** (`distortOption > 0`). Provided the drive is a
  `π_source`-preserving move at the source temperature (velocity randomisation + a proper HMC Metropolis
  accept), the driven endpoint remains a source-temperature Boltzmann sample with **unit configurational
  Jacobian**, so the standard PT exponent `-Δβ ΔU` on the driven endpoints (`ETerm_nonequil`, no
  `lnJac`) is the correct acceptance (Nilmeier Eq 19 with `π`-preserving propagation ⇒ `w` reduces to
  the endpoint energy difference and no volume term). RENEMC is then REMC with a `π`-preserving
  decorrelation drive between swaps.

**INV-1 is satisfied by construction:** BAT scaling (volume-changing) always carries its Jacobian
(RENE); the velocity drive (volume-preserving) has no configurational Jacobian to carry (RENEMC).
**PRECONDITION (Critical):** RENEMC SHALL be paired only with a volume-preserving, `π_source`-preserving
drive. If RENEMC is ever configured with a position-scaling drive (`distortOption < 0`), its acceptance
omits the required Jacobian and is biased; the port SHALL assert against that pairing (INV-1 guard).

### B10. REBASONTOP vs RENE (answers "is this correct?")

**Confirmed: as shipped, `REBASONTOP` is behaviourally identical to `RENE`.** Every live branch
treats them the same (`Context.cpp:1329,1357`). The only intended distinction — temporarily setting
`runType = REMC`, running six neighbour-swap sub-rounds, then restoring `REBASONTOP` — is dead code
(`Context.cpp:2688-2699`). The name "on top" refers to layering REMC swaps on top of the
nonequilibrium exchange. **Resolved (D4): reinstate the interleaved-REMC sub-loop so `REBASONTOP` is
distinct and Python-exposed** — `RENE` work-swaps plus periodic REMC neighbour-swap sub-rounds.

### B11. REBAS body selection and testing modes (answers "remove hardcoded molecule + testing modes")

`REBAS_Scale_Mbx(molName, mbx)` (`HMCSampler.cpp:337-352`) decides whether a mobilized body is scaled,
from **hardcoded per-topology index tables** for `{ETHANE, ALA1, TRPCH}` (`HMCSampler.hpp:281-286`).
The molecule is hardcoded to `TRPCH` inside `perturbPositions` (`HMCSampler.cpp:387`), and
`testingMode` is hardcoded `false` (`HMCSampler.cpp:389`), gating an inert debug block
(`HMCSampler.cpp:399-514`) whose `scaleFactor` is chosen ad hoc (1.25/0.80 by parity, etc.) and which
prints `"SCALING IN TESTING MODE"` — a mode that does **not** preserve detailed balance.

The port SHALL:
- Remove the hardcoded `MOLECULE_NAME_Ix = TRPCH`, the `REBAS_Scale_Mbx` index tables, and the
  `testingMode` block from production source.
- Replace body selection with a **topology-driven** rule: the set of scaled mobilized bodies SHALL be
  derived from the molecular topology / world flexibility specification, not from literal indices.
  (Proposed: scale exactly the bodies whose parent world marks them flexible under the driven
  factorization; the researcher SHALL fix the precise rule.)
- Preserve the removed ad-hoc scaling schedules (1.25/0.80 alternation, Bernoulli 1.01, by-temperature)
  as **test fixtures** in a companion spec (see "REBAS test machinery" appendix), to become a
  regression test that checks the work/Jacobian bookkeeping against a known scale schedule.

### B12. Simplify `RunReplicaWorldRange` and neighbours (answers "can we simplify?")

`RunReplicaWorldRange(replicaIx, startWorldCnt, nofWorldsCounted, isNonEquilibrium, shouldPrint)`
(`Context.cpp:2410-2498`) has three defects that the port SHALL fix:

1. **The range is not enforced.** `nofWorldsCounted` is never used in the body; the loop iterates the
   entire `thermoWorldIxs` regardless of `startWorldCnt`/`nofWorldsCounted`. It is benign only because
   all live worlds are equilibrium; with a nonequilibrium segment it would run the wrong worlds. The
   loop SHALL respect `[startWorldCnt, startWorldCnt + nofWorldsCounted)`.
2. **Double indexing.** `for (worldScheduleIndex : thermoWorldIxs)` iterates the world-index *values*,
   then reads `worldIndex = thermoWorldIxs[worldScheduleIndex]` — indexing the list by a value
   (`Context.cpp:2430-2432`). Benign only because `thermoWorldIxs == [0..W-1]`; a genuine bug if the
   schedule is ever a permutation or subset. Iterate positions cleanly.
3. **Dead parameter.** `isNonEquilibrium` is never read; the equilibrium decision is per-world via
   `distortIx == 0` (`Context.cpp:2473`). Remove the parameter.

Proposed shape (position-clean, range-respecting):

```
void Context::RunReplicaWorldRange(int replicaIx, int startPos, int count, bool shouldPrint) {
    auto& replica = replicas[replicaIx];
    const auto& sched = thermodynamicStates[replica2ThermoIxs[replicaIx]].getWorldIndexes();
    const auto& distortOpts = ...getDistortOptions();
    replica.updWORK() = 0.0; replica.upd_WORK_Jacobian() = 0.0;
    for (int pos = startPos; pos < startPos + count; ++pos) {
        const int worldIx = sched[pos];
        if (pos > startPos) transferCoordsFromWorldToWorld(worlds[sched[pos-1]], worlds[worldIx]);
        RunWorld(worldIx, header, shouldPrint);
        const bool equilibrium = (distortOpts[pos] == 0);
        transferCoordsFromWorldToReplica(worlds[worldIx], replica, /*intoWORK=*/!equilibrium);
        if (equilibrium) thermoState.calcQStats(...);
        replica.incrementWorldsNofSamples(1);
    }
}
```

The many thin `store_WORK_*` / `storeReplica*` wrappers (`Context.cpp:1535-1594`) SHOULD be collapsed
into `transferCoordsFromWorldToReplica`. **Caution (revision 2, reviewer Should-fix R4): these wrappers
are DEAD (no callers), and they encode a *different, incompatible* aggregation scheme** — a single
front-world / back-world read (coordinates+potential from the front world, Jacobian from the back world)
— whereas the **live** scheme is per-world `+=` accumulation across the driven range (reset at
`Context.cpp:2425-2426`, accumulate at `:24-25`). The rewrite SHALL collapse toward the **live `+=`
accumulation** shown above (`updWORK() += ...; upd_WORK_Jacobian() += ...` per driven world) and delete
the front/back wrappers; it SHALL NOT resurrect their single-read semantics, which would drop (or
double-apply) the contributions of all-but-one driven world. With D7 (`mdSteps = 0`) a driven range is
typically a single scaling world, but the accumulation SHALL remain correct for a multi-world driven
segment. Guard: V5 on a **two-driven-world** schedule — `lnJac` SHALL be the sum over both, not the last
world's alone.

### Derivation sketch (nonequilibrium candidate exchange)

The paired drive + label swap is a deterministic **involution** `T` on the joint space `z = (x_X, x_Y)`:
scale X by `s`, scale Y by `1/s`, swap the temperature labels. With the sampler's configurational marginal
Cartesian-Boltzmann (INV-7), the target joint density is `π(z) ∝ exp[−β_C U(x_X) − β_H U(x_Y)]`. For a
Metropolis move with a deterministic involution, acceptance is `min(1, [π(Tz)/π(z)]·|det ∂T/∂z|)`
(`minh_2014_hastings` Eq 3). Because `T` is block-diagonal in the two replicas,

  `log α = −[β_H U(x_X^τ) − β_C U(x_X^0)] − [β_C U(x_Y^τ) − β_H U(x_Y^0)] + lnJac_X + lnJac_Y`.

Defining `Work_X = [β_H U(x_X^τ) − β_C U(x_X^0)] − lnJac_X` gives exactly `log α = −(Work_X + Work_Y) =
WTerm`. The **form, the β-assignment (`x^0` at source β, `x^τ` at target β), and `correctionTerm = 1` are
CONFIRMED correct**, and the no-drive limit (`x^τ = x^0`, `lnJac = 0`) collapses to `ETerm_equal` (V4).
This instantiates **Ballard & Jarzynski 2009** (replica exchange with nonequilibrium switches, `P_acc =
min{1, e^{−(w_A+w_B)}}`, reduced work `w = Δ(reduced potential) − ln J`), composed with the NCMC
deterministic-map Jacobian (**Nilmeier et al. 2011**). The in-code "Ballard-Jarzinski nomenclature"
(`Context.cpp:1148`) names it directly; the Chodera–Shirts citation (`Context.cpp:1126-1130`) covers only
the permutation/label-swap loop.

**Work-config timing (INV-8, revision 2 — REFUTES the "measure at x'" resolution).** In
`perturb_Q_QDot_QDotDot` a driven world runs `perturbPositions` (deterministic scaling `x^0 → x'`)
**then** `integrateTrajectory_BoundHMC` (MD, `x' → x^τ`), and `x^τ` is what is committed on accept
(`HMCSampler.cpp:3557-3567`; `Context.cpp:1389`). This is a **one-step NCMC move** (Nilmeier `T = 1`:
one perturbation followed by one propagation). Nilmeier's acceptance for such a move (eq:16, eq:19)
uses the **perturbation work** `w = u_target(x') − u_source(x^0)` plus the perturbation Jacobian, and
the propagation `x' → x^τ` cancels out of the acceptance (`ΔS = −q`, eq:18) **only if the propagation
kernel preserves the post-perturbation target** `π_target = e^{−u_target}` (eq:15).

The original violates this: the drive MD runs at the replica's **source** temperature
(`Context.cpp:1771-1779` sets each driven world's sampler to the replica's own `thermoStateTemperature`;
the scale `s = √(T_j/T_i)` is the only temperature-change mechanism, with **no** re-thermostatting to
the neighbour). So the propagation preserves `π_source`, not `π_target`; the `x' → x^τ` heat is
uncancelled and the two-endpoint `Δ(βU) − lnJac` — whether `U` is read at `x^τ` (original) or at `x'`
(the earlier RQ-1 patch) — is biased whenever `mdSteps > 0`. **Measuring `U` at `x'` does not fix this;
committing `x^τ` while scoring `x'` is internally inconsistent.**

Two theoretically exact regimes exist (Nilmeier eq:20 vs eq:19):

1. **`mdSteps = 0` for driven worlds (the SHALL for this port, D7).** With no propagation the move is a
   pure deterministic scaling map, `x^τ = x'`, and acceptance `= −(Work_X + Work_Y)` with `Work =
   β_target U(x') − β_source U(x^0) − lnJac` is **exact** (Nilmeier eq:20 symplectic/deterministic
   limit; Ballard–Jarzynski instantaneous-switch limit; the Nilmeier "instantaneous MC" move eq:27).
   V8 proves it. This retains RENE's advantage over REMC — the scaling pre-adapts the geometry toward
   the neighbour's ensemble, reducing the effective `ΔU` — while needing no new thermostat machinery.
2. **MD propagation re-thermostatted to `β_target` (deferred, Consequences).** Run the driven MD at the
   *post-perturbation* target temperature with a proper `π_target`-preserving HMC kernel and momentum
   inversion on reject (Nilmeier eq:12, eq:19); then the propagation relaxation improves acceptance and
   the work is still the single perturbation work at `x'`. This is the genuine finite-time
   nonequilibrium switch (Ballard–Jarzynski eq:13–14) and SHOULD be the follow-on optimisation.

The port SHALL satisfy INV-8 by regime 1 (`mdSteps = 0` on driven worlds, D7).

NOTE (revision 2): both references are now ingested — `ballard_2009_rens`
(`references/papers/ballard_2009_rens/`, DOI 10.1073/pnas.0900406106) and `nilmeier_2011_ncmc`
(`references/papers/nilmeier_2011_ncmc/`, DOI 10.1073/pnas.1106094108), both in `references/index.yaml`.
The acceptance form above cites `ballard_2009_rens` eq:2, eq:4.1/4.2 and `nilmeier_2011_ncmc`
eq:15–20, eq:27.

---

## Invariants

- **INV-1 (detailed balance).** For every run type, the outer exchange chain SHALL satisfy detailed
  balance over the joint (configuration, permutation) space. In particular the Jacobian of any BAT
  scaling SHALL appear in acceptance, and any stochastic scale-factor proposal SHALL carry its Hastings
  correction.
- **INV-2 (marginal preservation).** Each thermodynamic state's configurational marginal SHALL remain
  the Boltzmann distribution at its temperature; swaps and drives SHALL leave it invariant.
- **INV-3 (state/coordinate separation).** Committed coordinates live only on replicas; worlds hold no
  persistent per-replica state between rounds. A swap SHALL move labels, never coordinates.
- **INV-4 (commit atomicity).** On an accepted nonequilibrium swap, both the coordinates and all
  energies (potential, reference, Fixman) of the promoted state SHALL be committed together; on reject,
  none SHALL change. (The current asymmetric commit — coordinates + potential only — is a defect, F4.)
- **INV-5 (schedule integrity).** Propagating a world range SHALL touch exactly the scheduled worlds in
  order; `WORK`/`WORK_Jacobian` SHALL be reset once per driven range and accumulated only over
  nonequilibrium worlds.
- **INV-6 (energy consistency).** Reduced potentials in acceptance SHALL use each state's own β; a
  replica's stored `potential` SHALL equal the energy of its stored coordinates at swap time.
- **INV-7 (Fixman-in-sampler ⇒ Fixman-out-of-acceptance).** REX SHALL run with the Fixman potential
  enabled in every sampler; the swap acceptance SHALL then use the Fixman-free physical potential. Violating
  either side breaks INV-2 (D3).
- **INV-8 (work-config consistency, revision 2).** A driven world SHALL run zero post-scale MD steps
  (`mdSteps = 0`, D7), so the committed endpoint equals the scaled configuration `x^τ = x'` and the
  reduced potential entering `Work` is evaluated at exactly the configuration whose scaling Jacobian is
  in `lnJac`. Scoring `x'` while committing a post-MD `x^τ ≠ x'` — the earlier "measure at x'" patch —
  is a defect: the propagation heat is uncancelled unless the driven MD is re-thermostatted to the
  target temperature with a `π_target`-preserving kernel (deferred regime 2, Consequences).
  (Derivation sketch; `nilmeier_2011_ncmc` eq:15–20.)
- **INV-9 (scaling-anchor sharing).** The affine scaling anchor `μ` (the running BAT mean the drive
  scales the deviation around, B4) SHALL be a single value shared by both partners of a swap pair and
  frozen across the forward and reverse drive of that pair. Only then is the paired map an involution
  (`M_{1/s,μ} ∘ M_{s,μ} = id`) and `correctionTerm = 1` valid (D2). Using each thermodynamic state's own
  running mean (`μ_C ≠ μ_H`) breaks the involution and biases acceptance (D2, revision 2).
- **INV-10 (drive/run-type pairing).** RENE/REBASONTOP SHALL drive with volume-changing BAT scaling
  (`distortOption < 0`) and carry `lnJac`; RENEMC SHALL drive with a volume-preserving,
  `π_source`-preserving velocity/NMA move (`distortOption > 0`) and omit `lnJac`. Any other pairing
  biases acceptance (B9, INV-1).

---

## Interface

**NOTE (revision 2): none of the symbols in this section exist in the current tree.** Every `RUN_TYPE`
value, setter, and "removed" symbol below is **new work** (I1–I3) or an original-tree cleanup that does
**not apply** to the target (I5). Line numbers are original-tree unless tagged `[CURRENT]`.

### I1. Run types (PyBind) — NEW

The current engine has no `RUN_TYPE` enum and `runREX` takes no run-type argument
(`include/Context.hpp:81` `[CURRENT]`). The port SHALL add the enum, add a run-type parameter to the
exchange driver, and bind it. `src/PyBind11.cpp:412-416` `[CURRENT]` is currently `SystemTopology` CMAP
bindings — the new binding goes wherever the `Context` bindings live, not at that line. Target binding:

```cpp
py::enum_<RUN_TYPE>(m, "RunType")
    .value("DEFAULT",    RUN_TYPE::Default,    "No exchange; independent replicas.")
    .value("REMC",       RUN_TYPE::REMC,       "Replica Exchange MC (parallel tempering): accept on -Δβ·ΔU.")
    .value("RENEMC",     RUN_TYPE::RENEMC,     "Replica Exchange Non-Equilibrium MC (volume-preserving velocity drive): -Δβ·ΔU on driven endpoints, no Jacobian (INV-10).")
    .value("RENE",       RUN_TYPE::RENE,       "Replica Exchange Non-Equilibrium (BAT-scaling drive): accept on nonequilibrium work -(W_X+W_Y), includes lnJac; driven worlds run mdSteps=0 (INV-8).")
    .value("REBASONTOP", RUN_TYPE::REBASONTOP, "RENE work-swaps plus interleaved REMC neighbour swaps layered on top (D4).");
```

The Python enum docstrings SHALL be regenerated accordingly.

### I2. Distort option enum — EXTEND (not migrate)

The current engine **already** has `enum class DistortOption` (`include/World.hpp:59` `[CURRENT]`,
values `{NMA=0}`). There is no magic integer to replace. The port SHALL **extend** this existing enum
with the position-scaling variant used by RENE and the velocity variants used by RENEMC (B4), e.g. add
`ScaleBendStretch` for the BAT bond/angle scaling drive. The Python sampler config
(`python/robosample/context.py` `[CURRENT]`, the `add_sampler` call site — **not** `context.py:33-84`,
which is the `Context` class and dihedral helpers) SHALL accept the enum value directly; there is no
`-6..-1 / 0 / 1..6` wire convention to preserve. The RENE/RENEMC drive-sign semantics (INV-10) SHALL be
expressed as named enum members, not a raw sign test.

### I3. Mixing configuration — NEW

None of `setReplicaMixingScheme`, `nSwapAttempts`, `swapEvery`, `swapFixman` exist in the current
`include/Context.hpp` (131 lines `[CURRENT]`). All four SHALL be added. `swapEvery` SHALL gate exchange
frequency and `setReplicaMixingScheme(ReplicaMixingScheme)` SHALL select Neighboring/All; for `All`,
`nSwapAttempts` SHALL be configurable.

### I4. Observable behaviour that changes

- The label-swap `RunREX` (new) performs exchanges via the two inverse maps; the legacy coordinate-swap
  `runREX` (`src/Context.cpp:396-503` `[CURRENT]`) is retained as the INVARIANT-EQUIV oracle. A
  `PrintNofAcceptedSwapsMatrix`-style acceptance matrix becomes available.
- `RENE`/`RENEMC`/`REBASONTOP` are new functional run types (they never ran in the original).
- REBAS body selection is topology-driven (D5): no hardcoded per-molecule table.

### I5. Removed / deprecated symbols (original-tree cleanup — N/A to the target)

**NOTE: every symbol below has zero occurrences in the current tree** — `Replica::restoreCoordinates`,
`Replica::storeCoordinates`, `Replica::calcZMatrixBAT`, `calcZMatrixBAT_WORK`, `setReplicaExchangePairs`,
`swapPotentialEnergies`, `swapReferencePotentialEnergies`, `REBAS_Scale_Mbx`, `testingMode`,
`rexFlowOptions`, `rexWorkOptions`. There is nothing to remove from the target. This list is retained
only as a record of original-tree dead code the port SHALL NOT reintroduce when reimplementing the
subsystem.

---

## Validation strategy

Oracles tagged PRECONDITION / INVARIANT / LEMMA per `CLAUDE.md`.

- **V1 (PRECONDITION) — swap loop runs.** After the port, a 2-temperature `REMC` run over the
  alanine-dipeptide example SHALL record a non-zero attempted-swap count and a plausible acceptance
  ratio. Guards against the "no swaps happen" regression.
- **V2 (INVARIANT) — REMC marginal.** Two-replica `REMC` at temperatures `T` and `T` (degenerate)
  SHALL accept swaps with probability → 1 and leave each replica's energy histogram unchanged
  (self-consistency). At `T1 ≠ T2`, the measured potential-energy distribution of the replica *visiting*
  state `k` SHALL match a single long canonical run at `T_k` (KS test).
- **V3 (LEMMA) — swap acceptance symmetry.** `attemptREXSwap(i, j)` and `attemptREXSwap(j, i)` SHALL
  yield identical accept probability for `REMC` (symmetric proposal).
- **V4 (LEMMA) — work-acceptance reduces to REMC.** With all distort options zero (no drive), `RENE`
  and `RENEMC` acceptance SHALL equal `REMC` acceptance to within floating-point tolerance (`WTerm`
  degenerates to `ETerm_equal` when `x^τ = x^0` and `lnJac = 0`).
- **V5 (LEMMA) — Jacobian correctness, analytic targets.** Finite-difference of the *Cartesian*
  configuration under the scaling map SHALL match `getDistortJacobianDetLog()` to tolerance, against known
  values: single isotropic bond scale `s` → `lnJac = 3 ln s`; single angle scale → `ln(sin θ'/sin θ) + ln s`;
  general → `J(x') − J(x^0) + N_scaled·ln s`. A flipped `(J_fin − J_ini)` sign or a wrong `J_scale` fails
  this (D6, F7). **V5 is the *sole* guard on the Jacobian sign and magnitude** (see V8 NOTE below), so it
  SHALL NOT be dropped as redundant. V5 SHALL include at least one **≥3-body chain** that scales an
  *upstream* bond and an upstream angle, to confirm the geometric BAT tree in
  `calcMobodsBATJacobianDetLog_NEW` (which reaches a grandparent atom in an upstream body,
  `HMCSampler.cpp:5203-5207`) coincides with the scaled Simbody mobilizer coordinate — a single-DOF test
  cannot see a mismatch there. A spec-level standalone half (single-body, numpy-only) is delivered at
  `tests/bat_jacobian_scaling_check.py`; the coder extends it to the engine's real Cartesian map.
- **V6 (INVARIANT) — detailed balance under drive.** A driven `RENE` run on a small system (ethane or
  alanine dipeptide) SHALL reproduce the analytic canonical BAT-marginal (or a trusted long
  unconstrained run) within statistical error. This is the test that a plausible-but-biased
  implementation fails.
- **V7 (regression) — REBAS scale schedule.** The removed ad-hoc scaling schedules (B11) SHALL be
  reinstated as a fixture that drives a known `s(t)` and checks that accumulated `WORK` and
  `WORK_Jacobian` match a hand-computed reference.
- **V8 (LEMMA) — exact acceptance algebra.** On a one-body analytic system (single bond or angle,
  closed-form `U` and `J`), the computed `log α` SHALL satisfy `π(z)·P(z→Tz) = π(Tz)·P(Tz→z)` for the
  deterministic involution swap, compared **against the true `π`** — exercising the β-assignment and
  `correctionTerm` together, independent of MD sampling. **Correction (revision 2): V8 is BLIND to the
  Jacobian.** The symmetric paired map (scale X by `s`, Y by `1/s`) has joint `|det ∂T/∂z| = s^N·s^{−N}
  = 1`, so the Jacobian is identically absent from the *swap* acceptance; a flipped `lnJac` sign or a
  dropped `ln s` is invisible here. The earlier claim that V8 catches a Jacobian error was wrong. V8
  therefore guards β-assignment + `correctionTerm`; the Jacobian is guarded by **V5** and, additionally,
  by a **single-replica NCMC acceptance** where `|det| = s^{N} ≠ 1` (`nilmeier_2011_ncmc` eq:27/28) — a
  flipped sign or dropped `ln s` fails *that*. Both the V8 involution and the single-replica `|det|≠1`
  check are delivered in `tests/test_rex_swap_acceptance_algebra.py` (13 cases). The self-consistency
  form `log α(z) + log α(Tz) = 0` alone is insufficient — it is satisfied by any `T`-antisymmetric rule,
  including a wrong β-assignment — hence the comparison against true `π`.
- **V9 (PRECONDITION) — Fixman-on guard (INV-7).** A REX run configured with Fixman **disabled** in the
  sampler SHALL abort before sampling (a runtime assertion, not a comment). This is the discriminating
  guard for the INV-7 biconditional: with Fixman off, the configurational marginal is
  `|M_{N_f}|^{1/2} e^{−βU}` and the bare-`U` swap targets the wrong joint distribution (D3, Critical).
  Without V9 a biased Fixman-off run passes V5/V6/V8 silently.
- **V10 (LEMMA) — anchor round-trip (INV-9).** A coordinate-level assertion `M_{1/s,μ}(M_{s,μ}(q)) == q`
  SHALL hold to floating point for the shared frozen anchor `μ`, and SHALL **fail** if the two members
  of the pair use different state means (`μ_C ≠ μ_H`). Guards the D2 involution / `correctionTerm = 1`.

NOTE: V5 is the sharp Jacobian oracle; V8 the sharp acceptance-algebra oracle (β-assignment +
`correctionTerm`, blind to the Jacobian — see its correction); V10 the anchor-involution oracle; all
noise-free. V6 is the end-to-end correctness oracle but slow; V2 catches PT wiring; V4 checks the
no-drive limit; V9 fails-loud on the Fixman precondition.

---

## Consequences and trade-offs

- **Serial replicas.** Keeping `W` shared worlds preserves memory but forecloses parallel replica
  execution. Parallelism would require per-replica OpenMM Contexts (feasible — one shared System) plus a
  concurrency layer, and pays off only multi-GPU — see D1.
- **Redundant per-state option arrays.** Copying identical distort/flow/integrator vectors into every
  `ThermodynamicState` (`context.py:1062-1088`) is wasteful and error-prone. Moving per-world config to
  the `World`/sampler and leaving only temperature + `worldIndexes` on the state simplifies the model but
  changes the `addThermodynamicState` signature — an API break to weigh.
- **RENEMC vs RENE.** They are distinct estimators bound to distinct drives (INV-10): RENEMC =
  volume-preserving velocity/NMA drive with `-Δβ ΔU` on the decorrelated endpoints; RENE = volume-changing
  BAT scaling with the work term. Keeping both is cheap but doubles the validation surface. RENEMC's
  correctness depends on its drive being `π_source`-preserving; if the velocity/NMA drive lacks a
  Metropolis accept it only approximately preserves `π` and RENEMC becomes an approximate estimator —
  V2/V6 on the RENEMC path SHALL confirm the marginal.
- **`mdSteps = 0` on driven worlds (D7) forecloses in-drive relaxation.** The exact, minimal port makes
  the RENE drive a pure deterministic scaling map (INV-8), so acceptance falls off as the temperature
  gap widens with no relaxation to recover it. The follow-on **regime 2** — re-thermostatting the driven
  MD to the *post-perturbation target* temperature with a `π_target`-preserving HMC kernel and momentum
  inversion on reject (Nilmeier eq:12, eq:19; Ballard–Jarzynski eq:13–14) — recovers finite-time
  switching and higher acceptance. It requires: (1) per-drive target-temperature assignment distinct
  from the replica's own state temperature; (2) a proven `π_target`-preserving propagation with its
  Metropolis accept; (3) momentum reversal on reject. Deferred as a performance enhancement, not a
  correctness fix; recorded here and under a future `docs/decisions/` record.
- **Feasibility gap — `World` runtime setters.** `[CURRENT]` The current `World` fixes timestep, MD
  steps, and accept/reject mode at `add_sampler` time; `SamplerConfig sampler_` is private
  (`include/World.hpp:515,578` `[CURRENT]`) with no `setTimeStep`/`setMdSteps`/`setAcceptRejectMode` and
  no mutable accessor. `setTemperature` **does** exist (`include/World.hpp:390` `[CURRENT]`) and the
  existing coordinate-swap `runREX` already drives it per replica each sweep (`src/Context.cpp:458`
  `[CURRENT]`) — so the current engine performs REMC varying *only* temperature. The added setters become
  necessary once the port introduces the original's per-`ThermodynamicState` schedule (distinct
  timestep/mdSteps/mode per state); they are **not** a precondition for the existing REMC, which already
  swaps. NOTE the D7 constraint interacts here: driven worlds SHALL be configured `mdSteps = 0`, which
  the new `setMdSteps` (or the per-world schedule) SHALL enforce for `distortOption < 0` worlds.

---

## Resolved decisions

User directive (2026-07-12): keep the port **as close to the live original as possible**; the
nonequilibrium code "was commented for other reasons" and **SHALL be reachable through Robosample** —
all four exchange run types and the driven paths are reinstated, not dropped.

### D1 — Parallel replica execution: deferred; feasible only multi-GPU (was Q1)

**Verdict: not in scope for the port. Keep the serial, shared-engine model (which also matches the
original). It is architecturally feasible but the payoff is hardware-gated.**

Analysis of the current engine's OpenMM integration:

- Energy/force/integration is a **process-wide singleton**: `OpenMMContext::get()` returns one
  `static OpenMMContext` (`OpenMMContext.hpp:19-22`) holding one `OpenMM::System`, one
  `OpenMM::Context`, one `Integrator` (`OpenMMContext.hpp:261-263`); `initialize()` is defined once and
  called once (`OpenMMContext.cpp:417`, `Context.cpp:677`). Every consumer routes through
  `OpenMMContext::get()` — `ForceBridge` (every world's sampler), `NMA`, `World`, `Context`, `PyBind11`.
- Evaluation is **stateful and mutating**: `setPositions()` then compute then `getState()`
  (`OpenMMContext.hpp:116-126,195-199`). Two replicas evaluated on this one instance would race on the
  shared positions — hence today's replicas are necessarily serial.
- Additional process-global state compounds this: the CUDA kinematics device sub-object is a
  file-static "keyed to the singleton" (`OpenMMContext.cpp:220-221`), plus alchemy state, force-group
  labels, and `ForceBridge::posCache_`.

Is concurrent evaluation possible at all? **Yes, in principle.** OpenMM's
`Context(const System& system, Integrator& integrator)` takes the System by **const reference** and does
not own it (`openmm/openmmapi/include/openmm/Context.h:73`). Multiple `OpenMM::Context` objects can
therefore share one immutable `System` (topology is identical across replicas), each with its own
`Integrator` and its own device buffers — so `R` replicas could each own a Context and evaluate without
racing. Realising it requires: (1) de-singletonising `OpenMMContext` into a per-replica instance (or an
`R`-sized pool) and threading a handle through `ForceBridge`/`NMA`/`World` instead of
`OpenMMContext::get()`; (2) making the file-static CUDA kinematics device state per-context; (3) a
concurrency layer (threads or per-context CUDA streams).

**Why deferred:** the payoff is hardware-bound. For the target systems (10⁵–10⁶ atoms, explicit solvent)
a single energy evaluation already saturates one GPU, so `R` Contexts time-slicing one GPU give ~no
wall-clock gain. Real speedup needs one GPU per Context (`CudaDeviceIndex` pinning) or CUDA MPS; on the
CPU platform each Context is already multithreaded, so `R` of them oversubscribe cores. Net: a
medium-to-large refactor that pays off only in multi-GPU deployments. Record the path here; revisit under
`docs/decisions/` if multi-GPU replica parallelism becomes a goal. The serial model stands (B0 NOTE).

### D2 — Scale-factor proposal correction: `correctionTerm = 1` **conditional on INV-9** (was Q2; revised)

Keep the original deterministic `s = sqrt(T_j/T_i)` scaling with `correctionTerm = 1`, **subject to the
shared-anchor condition INV-9 below**. The scale factor itself is self-reverse — the even/odd pairing
gives `s_{C→H} = sqrt(T_H/T_C)`, `s_{H→C} = sqrt(T_C/T_H)`, whose product is 1 — but that argument alone
is **not sufficient**: the affine *shift* by the mean must also cancel, which it does only under INV-9.
The `perturbScalingFactor` randomisers (`Gauss`/`Bernoulli`/`uniform`, `Context.cpp:1893-1940`) are
preserved but default OFF; if any is enabled the Hastings `correctionTerm` SHALL be implemented (no longer
unity). Code SHALL `NOTE`/assert this precondition where `correctionTerm` is set.

Researcher CONFIRM: the affine map `M_s(q) = s·q − (s−1)·mean` satisfies `M_{1/s} ∘ M_s = id` exactly
**only when both compositions use the same `mean`** — `M_{1/s,μ'}(M_{s,μ}(q)) = q + (s−1)/s·(μ'−μ)`,
which equals `q` iff `μ' = μ`.

**REFINE (revision 2, reviewer REFUTE): the ported code violates this precondition.** The original feeds
each drive the *current thermodynamic state's* running mean (`Context.cpp:2026-2027`,
`getBMps_means`/`getPFrs_means`). After the forward drive of the cold config uses `μ_C` and the label
swap, the reverse drive of that same config (now at the hot state) uses `μ_H`; since `μ_C ≠ μ_H`
generically, `T ∘ T ≠ id`, the paired map is **not** an involution, and `correctionTerm = 1` with `lnJac`
only is **wrong**. The simple two-endpoint (involution) acceptance the port uses — no momentum reversal —
is valid only for a genuine involution.

**Resolution (INV-9): the scaling anchor `μ` SHALL be a single value shared by both partners of a swap
pair and frozen across the forward/reverse drive.** Concretely, the drive SHALL scale the deviation
around a **state-independent** running BAT mean per coordinate (accumulated across equilibrium samples of
all states, owned by `Context`, not by `ThermodynamicState`), or equivalently a per-pair anchor frozen at
round start. The anchor is only the pivot of the affine scaling; a temperature-independent pivot makes
the paired map an exact involution, so `correctionTerm = 1` holds with no proposal-density term beyond
`lnJac`. This revises the B1/B3 ownership: the *scaling-anchor* mean moves to shared/global ownership;
the per-state `Q`-statistics (Qmeans/Qvars, B3) remain per-`ThermodynamicState` as diagnostics. V10
(anchor round-trip) gates this. The `perturbScalingFactor` randomisers (`Gauss`/`Bernoulli`/`uniform`,
`Context.cpp:1893-1940`) are preserved but default OFF; if any is enabled the Hastings `correctionTerm`
SHALL be implemented (no longer unity). Code SHALL `NOTE`/assert INV-9 where `correctionTerm` is set.

### D3 — Fixman excluded from the driven work (was Q3)

Keep the original's live "variant 2" (`Context.cpp:1245,1249`): the nonequilibrium work uses the reduced
reference potential and the log-Jacobian only; **Fixman does not enter the work term**. Rationale from the
prior campaign: the sampler's `+½RT ln|M|` cancels the momentum Jacobian, so each state's configurational
marginal is flat-Cartesian Boltzmann and Fixman need not appear in the outer acceptance. `swapFixman` is
ported as an **off-by-default diagnostic** (computed, logged, never added to `log_p_accept`). NOTE: engine
`calcFixman` sign `+½RT ln|M_φ/M_3N|` is authoritative; the `EnergySnapshot.hpp` comment `−½kT ln det M`
is stale.

Researcher CONFIRM (extends to the driven RENE case), with a biconditional the spec now states explicitly.
From `spiridon_2017` (eq:2, eq:3): with Fixman `U_F` in the acceptance Hamiltonian, the sampled
configurational marginal is Cartesian-Boltzmann `e^{−βU(x)}` (because `|M_{3N}|^{1/2} = |M|^{1/2}·|det J|`
absorbs the internal-coordinate Jacobian). This holds under the drive too: `lnJac` is the Cartesian map
Jacobian and `x^τ` is drawn by Fixman-on HMC, so no `det M` term enters the work. **The exclusion is correct
*because* Fixman is enabled in the sampler (INV-7).** PRECONDITION (Critical): if Fixman is disabled in the
sampler, the marginal is `|M_{N_f}|^{1/2} e^{−βU}` and a bare-`U` swap targets the wrong joint distribution
— a Critical bias. REX SHALL assert Fixman-on, or fold `det M` back into both the sampler and the work
consistently.

### D4 — REBASONTOP: reinstate as distinct, exposed (was Q4)

Per "SHALL be reachable" and "closest to original intent", reinstate the commented interleaved-REMC block
(`Context.cpp:2688-2699`): `REBASONTOP` = `RENE` work-swaps **plus** periodic REMC neighbour-swap
sub-rounds on top. It SHALL be a first-class, Python-exposed run type (Interface I1), not a silent alias.
This is the one place we reconstruct commented behaviour rather than mirror the live (== RENE) behaviour,
because that live identity is exactly the "commented for other reasons" state the directive says to undo.

### D5 — REBAS body selection: topology-driven; original tables become a fixture (was Q5)

Deliberate, justified departure from the original. The original gates even its production `BendStretch6`
path by the hardcoded `REBAS_Scale_Mbx(TRPCH, mbx)` table (`HMCSampler.cpp:569,596`; molecule hardcoded at
`:387`), encoding literal atom indices for three specific test molecules — it cannot run any other system.
The **scaling algorithm** is kept verbatim (displace each scaled body's bond/angle by
`(s−1)·(current − running_mean)`); only the **selection** becomes topology-driven. The exact
`{ETHANE, ALA1, TRPCH}` index tables are preserved as a regression fixture (V7 / appendix) proving the
topology-driven rule reproduces the original selection on those three systems.

**Concrete selection rule (revision 2, resolves the "underspecified for torsional worlds" gap).** The BAT
scaling acts only on **bond-length** (`r`) and **bond-angle** (`θ`) generalized coordinates; torsions and
root DOFs are not scaled. The scaled DOFs are enumerated per flexible body from its `JointType`
(`include/RobotModel.hpp:27-38` `[CURRENT]`, via `RobotModel::bodyJoint`):

| `JointType` | DOFs | Scaled (bond `r` / angle `θ`) | `N_scaled` contribution |
| --- | --- | --- | --- |
| `Rigid` | 0 | — | 0 |
| `Torsion` | 1 (φ about bond) | none | 0 |
| `Slider` | 1 (translate along bond) | `r` | 1 |
| `Cylinder` | 2 (φ + translate) | `r` (φ excluded) | 1 |
| `BendStretch` | 2 (⊥ rotation θ + translate `r`) | `θ`, `r` | 2 |
| `SphericalCoords` | 3 (azimuth, zenith, radius) | zenith `θ`, radius `r` (azimuth excluded) | 2 |
| `Cartesian` / `Ball` / `FreeLine` / `Free` | 3–6 (root/orientation) | none (no bond/angle coordinate) | 0 |

`N_scaled` (D6/RQ-2) is the sum of these contributions over the driven world's flexible bodies, and the
Jacobian `J(x) = Σ_bodies [2 ln r + ln sin θ]` sums only the `r`/`θ` coordinates a body's joint actually
exposes. **A purely torsional Gibbs block** (all `Torsion` joints — the canonical dihedral factorization)
therefore has `N_scaled = 0`, `lnJac = 0`, and no scaled coordinate: the RENE drive is inert for that
block and its swap **correctly degenerates to REMC** (a pure-dihedral world has no bond/angle volume to
change). This is the intended behaviour, not a gap. The selection is well-defined for every joint type;
the only implementation freedom is which of the scalable joint types (`Slider`/`BendStretch`/`Cylinder`/
`SphericalCoords`) a given driven factorization actually instantiates — a property of the world's
mobility spec, read at Jacobian time.

### D6 — Canonical BAT-Jacobian: corrected composition (was Q6; researcher REFUTE)

The per-body volume element `J(x) = Σ_bodies [2 ln r + ln sin θ]` (`calcMobodsBATJacobianDetLog_NEW`,
`HMCSampler.cpp:5147-5226`) is CONFIRMED correct (torsion contributes nothing). But the **original's
composition `J_ini + J_scale − J_fin` (`HMCSampler.cpp:659`) is REFUTED — it ships a biased Jacobian.**
The correct forward Cartesian log-Jacobian of the deterministic scaling map `x^0 → x'` is

  `lnJac = ( J(x') − J(x^0) ) + N_scaled · ln s`,

where `s = QScaleFactor` and `N_scaled` is the number of scaled bond+angle DOFs (if per-coordinate scale
factors are ever used, `N_scaled·ln s` becomes `Σ_dof ln s_dof`). Numeric check: a single isotropic bond
scale gives Cartesian `|∂x'/∂x| = s³`, and this formula returns `2 ln s + 1·ln s = 3 ln s`; the original
returns `≈0` (double-counted) or `−ln s` (even after F8). The port SHALL: set `J_ini = J(x^0)`,
`J_fin = J(x')`, `lnJac = (J_fin − J_ini) + N_scaled·ln s`; delete the value-ratio `J_scale` block
(`HMCSampler.cpp:588-644`), the dead `internNdofs·log s` path (`:4399-4494`), and the `4 ln r + 2 ln sin θ`
variant (`Replica.cpp:565-586`; `Context.cpp:3634-3655`); fix the no-op `x != SimTK::NaN` guard →
`std::isnan`. V5 (with the analytic targets below) gates this. `N_scaled` requires the D5 selection rule to
enumerate the scaled DOFs at Jacobian time (RQ-2).

### Resolved by the researcher theory pass (2026-07-12)

- BMps/PFrs semantics CONFIRMED and corrected (B3 NOTE): `BMps` = bond-length coordinate,
  `PFrs` = bond-angle coordinate.
- Work acceptance CONFIRMED as the Ballard–Jarzynski involution rule (Derivation sketch); two references
  to ingest (below). Two residual questions remain (RQ-1, RQ-2).

### RQ-1, RQ-2 — finalised (2026-07-12)

Both are resolved under the **Fidelity principle** (Motivation): the nonequilibrium/WORK/REBAS path never
executed in the original (it is commented out, F1/F2), so reinstating it is a fresh, correct implementation
of dormant scaffolding — not a change to working original behaviour. "As close to the original as possible"
binds the parts that actually ran (REMC, the data model, the label-swap, the two-endpoint acceptance
*structure*); the dormant path is implemented **correctly** per the Ballard–Jarzynski derivation, with every
deviation from the scaffolding flagged.

- **RQ-1 — SUPERSEDED by D7 (revision 2).** The earlier decision ("keep two-endpoint, measure `U` at
  post-scale/pre-MD `x'`") is **withdrawn**: it does not fix the bias. With `mdSteps > 0` the driven MD
  (`integrateTrajectory_BoundHMC`, at the *source* temperature, `Context.cpp:1771-1779`) moves `x' →
  x^τ`, and `x^τ` is committed on accept while the work is scored at `x'` — the propagation heat is
  uncancelled because the kernel preserves `π_source`, not the post-perturbation `π_target` (Nilmeier
  eq:15, eq:18–19). See the Derivation-sketch work-config-timing analysis. The correct resolution is D7.
- **RQ-2 — DECIDED (unchanged).** `N_scaled` is the count of scaled bond+angle DOFs enumerated from the D5
  topology-driven selection (the driven world's flexible bodies) at Jacobian time; with per-coordinate
  factors, `N_scaled·ln s` becomes `Σ_dof ln s_dof`. Mechanical consequence of the D5 rule table.

### D7 — Driven worlds run zero post-scale MD (`mdSteps = 0`) (revision 2, reviewer REFUTE of RQ-1)

The RENE/REBASONTOP two-endpoint work acceptance `−(Work_X + Work_Y)` with `Work = β_target U(x') −
β_source U(x^0) − lnJac` is exact **iff the driven move is the pure deterministic scaling map** `x^τ =
x'` — i.e. the driven world runs **no** post-scale MD (Nilmeier eq:20/eq:27 deterministic limit; the
Ballard–Jarzynski instantaneous-switch limit). The original's post-scale MD at the source temperature
breaks the `ΔS = −q` cancellation (Nilmeier eq:15, eq:18) and biases acceptance whenever `mdSteps > 0`.

**Decision:** for every driven world (`distortOption < 0`), the port SHALL set `mdSteps = 0`. The RENE
drive is then a single deterministic BAT-scaling perturbation whose Jacobian `lnJac` (D6) fully accounts
for the volume change, and V8 proves the resulting swap satisfies detailed balance. This keeps RENE's
advantage over REMC — the scaling pre-adapts geometry toward the neighbour's ensemble, shrinking the
effective `ΔU` — with no new machinery. The finite-time alternative (post-scale MD re-thermostatted to
the *target* temperature with a `π_target`-preserving kernel and momentum inversion on reject) is a
correct, higher-acceptance follow-on (Consequences, regime 2) and is **deferred**, not adopted now.
INV-8 is a hard SHALL; V6 SHALL exercise the driven path to confirm the canonical marginal, and the
implementation SHALL assert `mdSteps == 0` for `distortOption < 0` worlds.

### References (ingested — revision 2)

Both are in `references/index.yaml` with extracted equations under `references/papers/`:

1. **Ballard, A. J.; Jarzynski, C.** "Replica exchange with nonequilibrium switches." *PNAS*
   **106**(30):12224-12229 (2009). DOI `10.1073/pnas.0900406106`. Key `ballard_2009_rens`. The result the
   code's "Ballard-Jarzinski nomenclature" (`Context.cpp:1148`) instantiates: `P_acc = min{1,
   e^{−(w_A+w_B)}}` (eq:2); reduced work `w = h_target(x_τ) − h_source(x_0) − ln J` (eq:4.1/4.2); the
   finite-time switch dynamics and Jacobian (eq:13–14); quasi-static limit `w → Δf`, `P_acc → 1` (eq:15).
2. **Nilmeier, J. P.; Crooks, G. E.; Minh, D. D. L.; Chodera, J. D.** "Nonequilibrium candidate Monte Carlo
   is an efficient tool for equilibrium simulation." *PNAS* **108**(45):E1009-E1018 (2011). DOI
   `10.1073/pnas.1106094108`. Key `nilmeier_2011_ncmc`. Perturbation work `w = Σ_t[u_t(x_t^*) −
   u_{t-1}(x_{t-1})]` (eq:16); heat and `ΔS = −q` (eq:17–18); reversible-MCMC acceptance `∝ e^{−w}`
   (eq:19); symplectic/deterministic limit `∝ e^{−Δu}·|det|` (eq:20, eq:27/28); momentum reversal on
   reject (eq:12). These are the basis of D7/INV-8.

---

## Appendix: fixes for dead/stubbed code

| ID | Location | Defect | Proposed fix | Severity |
| --- | --- | --- | --- | --- |
| F1 | `Context.cpp:2692,2704,2745` | `prepareExchangePairs` never called → no swaps | Call it live (or inside `mixReplicas`, B7). | Critical (no exchange) |
| F2 | `Context.cpp:2701-2799` | Nonequilibrium propagation commented out → `WORK_*` never populated | Reinstate the driven segment guarded by `nofNonequilibriumWorlds>0` (B8). | Critical (RENE/RENEMC inert) |
| F3 | `Context.cpp:1266-1267` | `correctionTerm` hardcoded to 1 | Keep = 1 **only under INV-9** (shared frozen anchor, D2 revision 2); with per-state means the deterministic map is NOT an involution and 1 is wrong. Implement the Hastings ratio only if a randomiser is enabled. | Correct as-is **iff** INV-9 holds; Critical otherwise (D2 revision 2) |
| F4 | `Context.cpp:1386-1395` | Two bugs: (a) WORK commit is **unconditional** across all run types → an accepted REMC swap calls `set_WORK_CoordinatesAsFinal` on a never-populated `WORK_atomsLocations`, reverting coordinates; (b) INV-4 asymmetry — only coords + `potential` committed, **not** `referencePotential`/`FixmanPotential` from their WORK counterparts | Gate the WORK commit on driven run types only; add `Replica` helpers to commit coords + `potential` + `referencePotential` + `FixmanPotential` **atomically** (INV-4). Assert INV-6 (stored reference potential == energy of committed coords) at swap time. | **Critical (wrong distribution + biased subsequent swaps)** — upgraded from High |
| F5 | `Context.cpp:1116-1124,1406-1418` | `rewindReplica()` throws "Not implemented", but its **only** call site (`Context.cpp:1407`) is commented out → the throw is **unreachable dead code**. The reject branch returns `false` and mutates nothing (the trial touched only `WORK_*`, never committed `atomsLocations`); the next round reseeds worlds from committed coords and resets `WORK` — so "no restore needed" is already the live behaviour. | Delete dead `rewindReplica`; no restore logic to implement. | **Nit** — downgraded from High (reviewer REFUTE of the "throws at runtime" description) |
| F6 | `Context.cpp:1211-1234,1244,1248` | Fixman computed but never enters acceptance; silent | Correct by design given INV-7 (Fixman on in sampler); state the biconditional in code + assert Fixman-on (D3). | Medium |
| F7 | `HMCSampler.cpp:659,588-644,5147-5226,4399-4494`; `Replica.cpp:565-586`; `Context.cpp:3634-3655` | Jacobian **composition** wrong: `(J_fin−J_ini)` sign inverted and `J_scale` bogus → biased `lnJac` (single-bond gives ~0 vs correct `3 ln s`); plus two dead variants and a no-op NaN guard | Set `lnJac = (J_fin−J_ini) + N_scaled·ln s` (D6); delete `J_scale` block and both dead variants; `std::isnan` guard. | **Critical (biased acceptance)** |
| F8 | `HMCSampler.cpp:621-634` | `perturbPositions` double-counts `J_scale` | Subsumed by F7 / D6: the entire `J_scale` value-ratio block is deleted in favour of `N_scaled·ln s`. | Subsumed by F7 |
| F9 | `Context.cpp:2410-2498` | `RunReplicaWorldRange` ignores its range, double-indexes, dead param | Rewrite per B12. | High (latent) |
| F10 | `HMCSampler.cpp:337-352,387,389,399-514` | Hardcoded molecule/tables + inert testing block | Remove; topology-driven selection; schedules → test fixture (B11, V7). | High |
| F11 | `bgeneral`/`Context.hpp:420` | `replicaMixingScheme` never read; `mixAllReplicas` never called | Wire scheme selection (B7). | Medium |
| F12 | `Replica.cpp:115-120,189-190,289-560` | Empty/assert-false stubs (`restore/storeCoordinates`, `calcZMatrixBAT*`, printers) | Remove; rely on Context transfer helpers and derived BAT. | Medium (maintainability) |
| F13 | `ThermodynamicState` | `rexFlowOptions`/`rexWorkOptions` set but never read. (Correction: the `nonequilibrium` flag **is** set nonzero — `setThermostatesNonequilibrium`, `Context.cpp:1007`, called at `:2591`; the real defect is that its reader `hasNonequilibriumMoves()` (`ThermodynamicState.hpp:81`) has **no caller**.) | Remove or define a consumer for all three. | Low |
| F14 | `Context.cpp:1658-1683,1101-1113` | `setReplicaExchangePairs`, `swap{Reference,}PotentialEnergies` legacy/unused | Remove. | Low |

### REBAS test machinery (companion spec — implement later)

Extract the removed scaling schedules into a regression fixture:

- **Constant** `s = 1` (identity; `WORK` and `WORK_Jacobian` SHALL be 0).
- **Alternating** `s ∈ {1.25, 0.80}` by sample parity (self-reverse over two samples).
- **By-temperature** `s = QScaleFactor` vs `1/QScaleFactor`.
- **Bernoulli** `s ∈ {1.01, 1/1.01}`.

The fixture drives each schedule on ethane (5 heavy bodies) and alanine dipeptide, and asserts the
accumulated `WORK` and `WORK_Jacobian` equal a hand-computed reference (V7). This becomes the test that
guards the work/Jacobian bookkeeping once the production REBAS path is topology-driven.
