# World selection for class-A GPCRs (receptor-agnostic)

Status: research spec (conceptual, portable). Read-only design guidance for
Robosample world / `Context` configuration (Python run scripts + CV analysis).
No engine (C++/CUDA) change is implied. General companion to the
FFAR1-specific spec `docs/specs/gpcr-world-design/00-ffar1-activation-worlds.md`
(read that first; this file reuses its general parts and drops all FFAR1 residue
numbers — the designer maps Ballesteros–Weinstein positions to their own receptor).

Anchoring references:
- `spiridon_2020_robosample` (primary Robosample paper: a *world* = rigid bodies
  + joint types mapped onto one molecular graph; blocked Gibbs over worlds;
  Fixman correction eq:7; guidance to define rigid bodies by secondary structure
  and to include at least one fully-flexible world per cycle).
- `spiridon_2017_cdhmc_gibbs` (CDHMC as a Gibbs move; few-DOF constrained moves
  mixed with unconstrained HMC shorten mean-first-passage time to the rarest
  state, §3.4, Tables 3–4).
- Codebase API: `python/robosample/context.py`
  (`build_flexibilities`, `add_torsional_world`/`add_robotic_world`,
  `add_cartesian_world`), `python/robosample/robo_bindings.pyi`
  (`JointType` enum: `Rigid`/`Torsion`/`Cylinder`/`Ball`/`SphericalCoords`/`Free`),
  `python/robosample/molecule_prototype.py:94-114` (root selection).
- Companion diagnostic spec: `docs/specs/reaction-force-monitoring.md`.
- External structural biology cited with URLs in §9 (same general sources as the
  FFAR1 spec).

---

## 1. Problem restatement

**User phrasing.** "How do I choose which mobilizers/joints to make flexible in
Robosample worlds so a generic class-A GPCR actually samples its
activation-relevant conformational transitions?"

**Restatement in codebase vocabulary.** A GPCR loads as one *robot*: a single
BFS-rooted kinematic tree (`atoms_root_index`; the molecule's compound index 0)
whose rigid bodies are joined by *mobilizers*. A *world* assigns a `JointType`
to a chosen set of bonds (`build_flexibilities`); every unselected non-terminal,
non-ring bond is `Rigid` (welded). Blocked Gibbs sampling cycles worlds; each
world draws momenta and integrates only *its* mobile DOFs, then
Metropolis-Hastings-accepts against the full OpenMM energy plus the Fixman
correction for the reduced-coordinate marginal (`spiridon_2020_robosample`
eq:7–8). The design lever is therefore **which bonds get which mobility in each
world, and how worlds compose across a Gibbs cycle**, chosen so the *slow*
collective coordinates that gate activation are proposed as few-DOF,
large-amplitude moves, while a fully-flexible world preserves ergodicity.

The underlying problems, in order:

1. **Identify the receptor's slow activation coordinates as concrete
   bonds/atoms** — using the conserved class-A architecture (§4) only as a
   *template*, then verifying each element exists and moves in the target
   receptor (§4.2).
2. **Map each coordinate to a mobility class** — backbone φ/ψ hinge, sidechain
   χ rotamer, or welded-helix rigid-body split — respecting Claim C1 (§3, §5).
3. **Compose a minimal, ergodic set of worlds** (§6).
4. **Define frame-invariant CVs** that confirm the intended state moved (§7),
   and use the force-monitoring diagnostic to refine the joint choice (§8).

**Reading chosen.** "Activation-relevant transitions" is read as *the
apo/inactive ↔ agonist/active transition of the 7TM bundle* (the hallmark TM6
cytoplasmic outward swing, connector/microswitch rearrangements, ligand-pocket
opening), not ligand (un)binding kinetics and not lipid dynamics. Rigid-body
coarse-graining does not reproduce kinetics (`spiridon_2020_robosample`
limitation); the goal is basin-to-basin *configurational* sampling. If the user
instead means binding kinetics or allosteric-modulator pathways, this spec's CVs
and worlds are not the right instrument and the framing must be revisited.

---

## 2. Terms (defined once)

- **Robot / kinematic tree.** One molecule loaded as a spanning tree of rigid
  bodies joined by mobilizers, BFS-rooted at compound index 0.
- **World.** A Gibbs block: a `JointType` assignment over the tree's bonds. A
  *torsional world* assigns `Torsion` (1-DOF pin about the bond axis) to a
  selected set and welds the rest; a *Cartesian world* frees all atoms in flat
  Cartesian OpenMM MD.
- **Mobility class.** The kind of DOF a motion is represented by: backbone φ/ψ
  torsion, sidechain χ torsion, or a rigid-body helix split (welded core +
  freed flanking loops).
- **Ballesteros–Weinstein (BW) `x.50` numbering.** Receptor-agnostic helix
  coordinates; `x.50` is the most conserved residue of helix TM`x`. Used
  throughout so the spec is portable; the designer maps `x.50` anchors to native
  residue indices of their construct (PRECONDITION P1).
- **Microswitch.** A small, conserved packing/rotamer element whose toggle is
  coupled to the global inactive↔active equilibrium (DRY, NPxxY, CWxP, PIF, Na⁺
  pocket; §4.1).
- **Activation CV.** A frame-invariant order parameter that increases (or
  toggles) on activation; primary one is the TM3–TM6 intracellular distance
  (§7).
- **Tree split.** The two-component partition a freed backbone joint induces
  (Claim C1).

---

## 3. The kinematic frame (why world design is constrained) — Claim C1

**Claim C1 (tree split).** A GPCR is one spanning tree rooted at a single fixed
atom (BFS compound index 0). Freeing a set of backbone joints and welding the
rest partitions the receptor, *per freed edge*, into exactly two maximal rigid
components: everything *inboard* (root side) of the edge stays fixed relative to
the base body; everything *outboard* (tip side) moves rigidly. There is **no
tree move that displaces an interior helix segment while holding both of its
covalent flanks fixed** — that would require a cycle, and the bundle's mutual
packing (TM6 against TM3/TM5/TM7) is non-covalent, hence absent from the tree
and enforced only energetically by OpenMM.

Consequences that constrain every design below:

- A conceptual "TM6 outward swing" world SHALL be read as *"free a backbone hinge
  and let the whole outboard component ride"*, never as "move TM6 alone."
  Whichever hinge is freed, the entire outboard component moves rigidly. The
  design job is to place the hinge so the two components are a *meaningful*
  collective mode and the moving component is small enough to accept (§6.4).
- **The root is not user-relocatable.** The BFS root is auto-chosen as the
  heaviest real (massive) *leaf* of the spanning tree
  (`molecule_prototype.py:94-114`; a massive degree-1 atom so the Z-matrix root
  triplet forms, ties broken by mass), and is always compound 0. The
  `atoms_root_index` setter stores a global index into an ordering already
  computed from compound 0; it does not re-run BFS from a chosen atom. So one
  cannot re-root to make the cytoplasmic helix ends the moving tips without an
  engine/API change (resolved Q1 in the FFAR1 spec). Which side of a given hinge
  is "fixed" vs "moving" is determined by the auto-selected leaf and the
  TM1→…→TM7 chain threading, not chosen by the designer.
- **Welding a helix interior is how it becomes a rigid body.** Not selecting a
  helix's bonds in a torsional world welds them; this realizes
  `spiridon_2020_robosample`'s "define rigid bodies by secondary structure."
  Rigid-body helix worlds (§6) rely on this.

NOTE. Because Claim C1 holds regardless of root, the §6 baseline worlds are
correct for any auto-chosen root. Root choice affects *which component is inboard*
(a labelling of the split), not whether a split exists. Exposing a root override
is a possible future engine task if a hinge world's acceptance proves
root-limited (see LEMMA L1 drift check).

---

## 4. The canonical activation architecture as a TEMPLATE (with a loud caveat)

### 4.1 Default candidate regions (conserved class-A microswitches)

Present the following as the DEFAULT candidate set — the regions a designer
should first consider for CVs and worlds, expressed in BW anchors so they port
across receptors:

- **DRY motif / ionic lock** — TM3 `3.49-3.50-3.51` (Asp-Arg-Tyr), with the
  R`3.50`–E`6.30` intracellular salt bridge ("ionic lock") that stabilizes the
  inactive state; breaking it accompanies TM6 opening.
- **NPxxY motif** — TM7 `7.49-…-7.53`, with the conserved Y`7.53` whose rotamer
  and TM7 position rearrange toward the bundle core on activation.
- **CWxP toggle** — TM6 `6.47`(Cys)-`6.48`(Trp toggle)-`6.50`(Pro); the W`6.48`
  rotamer ("toggle switch") repacks with agonist binding.
- **PIF / connector (transmission switch)** — P`5.50`, I`3.40`, F`6.44`; a
  hydrophobic connector that repacks between the ligand pocket and the
  cytoplasmic switches.
- **Na⁺ pocket** — centered on D`2.50`; a conserved allosteric microswitch that
  collapses on activation.
- **Hallmark TM6 cytoplasmic outward swing** — the largest activation motion;
  the cytoplasmic end of TM6 swings away from TM3 to open the G-protein cavity,
  canonically pivoting at the P`6.50` (CWxP) proline kink, which splits TM6 into
  two segments and amplifies the swing.

### 4.2 The caveat (SHALL)

The designer SHALL verify that each templated element *exists and moves* in
their receptor before building a world or CV around it. Specific class-A
receptors delete or immobilize some of these switches; a world built on an
absent or static element samples DOFs that do not gate activation, and a CV
built on it reports nothing.

Cautionary case (cite): FFAR1/GPR40 lacks the DRY ionic lock (has K, not E, at
`6.30`), lacks the conserved NPxxY tyrosine (Y`7.53` absent, TM7 disordered),
keeps the `6.48` toggle position static (Val, not Trp), has a static PIF/
connector, and pivots TM6 at a receptor-specific diglycine *cytoplasmic to* the
canonical proline kink — with a smaller opening amplitude. See
`docs/specs/gpcr-world-design/00-ffar1-activation-worlds.md` §3.2 and the PNAS
2023 source (§9). The lesson generalizes: **treat §4.1 as a hypothesis to test
per receptor, not a fact to transplant.**

NOTE (how to verify cheaply). If two endpoint structures (inactive and active,
or apo and agonist-bound) are available, verify a switch "moves" by measuring its
CV (§7) between the two structures; if only one structure exists, verify each
switch residue *exists* at its BW position (PRECONDITION P1) and treat mobility
as an assumption to be confirmed by the per-world CV-response check (§10).

---

## 5. Mapping mechanism → mobility class (the reusable rule)

The general, receptor-independent heuristic: classify a candidate motion by what
its natural reaction coordinate *is*, then pick the matching mobility class.

| Motion character | Reaction coordinate | Mobility class | Realization |
|---|---|---|---|
| Rigid-segment pivot / helix swing about a kink (TM6 outward swing) | backbone dihedral(s) at the pivot residue | **backbone φ/ψ hinge** | free φ/ψ at the pivot + adjacent loop; weld the helix cores (Claim C1 two-component split) |
| Conserved packing/rotamer toggle (W`6.48`, Y`7.53`, R`3.50` ionic-lock make/break, D`2.50` Na⁺ coordination, PIF repack) | a sidechain χ angle | **sidechain χ rotamer** | free χ of the switch residues only |
| Whole-helix repositioning as a body (e.g. TM5/TM7 shifting toward TM3) | a rigid displacement of a secondary-structure element | **rigid-body helix** | weld the helix core, free only the flanking loop φ/ψ so the helix reorients as a body |
| Local backbone breathing (ligand-pocket / extracellular helix ends gating ligand access) | φ/ψ of a short backbone stretch (+ pocket χ) | backbone φ/ψ (local) + χ | free the local loop/helix-end φ/ψ and the pocket χ |

Rules of thumb:

- **A rotamer toggle is a χ move, not a backbone move.** Freeing φ/ψ to sample a
  microswitch that is actually a sidechain flip wastes DOFs and injects backbone
  strain. SHOULD prefer the narrowest mobility class that carries the motion.
- **A large-amplitude helix swing is a φ/ψ hinge, not a χ move and not a
  Cartesian free-for-all.** Concentrate the freedom at the pivot; weld the helix
  core so the move is a genuine collective mode (few-DOF, large-amplitude).
- **A rigid-body helix reposition needs its flanks free and its core welded.**
  Freeing the core turns a collective mode into many small local moves and loses
  the large-amplitude advantage.
- **Backbone moves need companion χ repacking.** Any backbone or rigid-body
  helix move breaks and re-forms non-covalent bundle contacts that live only in
  the OpenMM energy (Claim C1); a repack world (§6, W1) must be in the cycle so
  those contacts can re-satisfy, or the coarse move rejects.

---

## 6. A composable world cycle (generic)

Notation: `{selection → JointType}`. `add_torsional_world(selection)` /
`add_robotic_world(selection)` build a torsional world (selected bonds →
`Torsion`, all else `Rigid`); `build_flexibilities` maps a bond/atom selection to
that mobile set; `add_cartesian_world()` builds the fully-flexible Cartesian
world. `Ball`/`Cylinder`/`SphericalCoords` are available for richer hinges
(§6.4).

The receptor SHALL be sampled by a Gibbs cycle drawn from the worlds below, and
the cycle SHALL include the ergodicity world (W0). Worlds are listed by scope,
not execution order.

### W0 — Ergodicity / fully-flexible (required)

`{all eligible bonds → Cartesian}` via `add_cartesian_world()`, or a
fully-flexible torsional world. Purpose: guarantee every DOF is reachable and let
the bundle relax the junction strain the coarse worlds inject. Short
trajectories, small `dt`.

### W1 — All-sidechain repack

`{χ of all receptor residues (+ lipid/tail χ) → Torsion}`. Purpose: rotamer
repacking that must accompany any backbone move for the bundle to re-form
contacts; cheap and high-acceptance. Moves no activation CV by itself, but
unblocks the coarse backbone worlds (§5 rule).

### W2 — Microswitch-χ world

`{χ of the conserved microswitch residues verified present in this receptor
(subset of DRY R`3.50`, NPxxY Y`7.53`, CWxP W`6.48`, PIF I`3.40`/F`6.44`, Na⁺
D`2.50`) + ligand-pocket χ → Torsion}`. Purpose: sample the toggles that couple
to the global equilibrium. Per §4.2, include only switches confirmed to exist
and move; a switch that is deleted or static in this receptor SHALL be omitted.

### W3 — TM6-hinge backbone world (the activation swing)

`{φ/ψ at the TM6 cytoplasmic pivot (default the P`6.50` kink; the
receptor-specific bend point where it differs) + the ICL3 / TM6 cytoplasmic
base → Torsion}`, with the TM6 core and all other TM cores welded. Purpose:
propose the TM6 cytoplasmic opening as a few-DOF, large-amplitude move about the
mechanistically correct hinge. Per Claim C1, this opens the cytoplasmic face by
moving the outboard component; it does not isolate TM6.

### W4 — Rigid-helix / flexible-loop world (global bundle repacking)

`{φ/ψ of all non-helix (loop) residues → Torsion}`, all TM cores welded. Purpose:
let the seven helices reorient as rigid bodies connected by flexible loops — the
coarse, large-step complement to W0, for collective helix repositioning
(TM5/TM7 toward TM3).

### 6.1 Composition rationale (from theory, not intuition)

`spiridon_2017_cdhmc_gibbs` §3.4 (Tables 3–4) shows the efficiency gain comes
from *few mobile DOFs taking large steps*: mixing constrained few-DOF moves with
unconstrained HMC shortens the mean-first-passage time to the rarest state
(isolated α_L) per integrator step. `spiridon_2020_robosample` requires a
fully-flexible world to restore ergodicity (the reduced-coordinate marginal
eq:6 differs from the flexible marginal eq:4 by the mass-matrix determinant; only
a world that reaches every DOF closes the gap). The cycle therefore alternates:
each coarse proposal (W3 hinge, W4 rigid-helix, W2 microswitch) is followed by a
repack (W1) and a flexible relax (W0) so the junction strain that Claim-C1 splits
inject at the weld boundaries is annealed before the next coarse proposal.

NOTE (ordering to anneal strain). A workable order is
microswitch/repack → hinge → rigid-helix → repack → flexible relax
(W2 → W1 → W3 → W4 → W1 → W0), mirroring a causal chain
switch → hinge → global reposition → repack → ergodic mix. The order is a
tuning choice, not a correctness requirement; the correctness requirement is
that W0 is present (INVARIANT I1) and that a repack world separates coarse
backbone proposals from the next coarse proposal.

### 6.2 Acceptance tuning

Per-world `dt` / trajectory length SHOULD be tuned toward the HMC target
acceptance (`spiridon_2020_robosample` acceptance-tuning workflow; the standard
HMC optimum is ≈ 0.65). A coarse world stuck near zero acceptance signals a
moving component that is too large (Claim C1 — shrink the freed hinge or add a
repack before it) or a welded-wall clash against the environment (INVARIANT I2).

### 6.3 Environment coupling (portable)

Lipid/solvent molecules SHALL keep mobile (`Free` roots + tail/sidechain χ in
W1). A coarse helix move that swings into *welded* lipid is a hard wall and
rejects — the welded-environment failure mode diagnosed for NCMC
(`docs/specs/ncmc-explicit-solvent/`). Any molecule given a `Weld` root
(nanodisc scaffold, a bound transducer) SHALL be an intended fixed reference, not
mobile lipid or the receptor (PRECONDITION P2).

### 6.4 NOTE on richer hinge joints

At a hinge, a `Ball` (3-rotational) or `Cylinder` (rotation+slide) mobilizer
samples a bend more completely than a single `Torsion`, at the cost of a larger
Fixman correction and a smaller stable `dt`. Optional and orthogonal to
correctness; start with `Torsion` (the paper default) and escalate only if the
hinge CV is under-sampled.

---

## 7. Frame-invariant collective variables (generic)

CVs SHALL be **frame-invariant** (internal distances/angles/dihedrals), because
the receptor root is `Free` (6-DOF) and the whole robot translates and rotates in
the box; a lab-frame RMSD would confound rigid-body drift with activation. This
is not a stylistic preference: under Claim C1, a freed hinge moves the outboard
component *rigidly through the lab frame*, so a lab-frame RMSD registers "change"
for a pure drift that did no activation work.

Portable CV set:

1. **TM3–TM6 intracellular distance — primary activation CV.** R`3.50` Cα to a
   cytoplasmic TM6 Cα (canonically the `6.30`/`6.34` region, or the receptor's
   nearest ordered cytoplasmic TM6 residue). Increases on activation. Use the
   *relative* increase from the run's own inactive starting frame unless an
   endpoint value is known.
2. **Hinge-dihedral-actuation check.** The backbone φ/ψ at the TM6 pivot that W3
   directly actuates. The opening CV (1) SHALL be checked to move *with* this
   hinge dihedral; if the opening CV moves with no hinge-dihedral change, the
   displacement is pure rigid-body drift of the outboard component (a Claim-C1
   artifact), not a hinge swing (LEMMA L1).
3. **Interhelical centroid distances** TM3–TM6, TM5–TM6, TM3–TM5 (Cα-centroid,
   frame-invariant): TM3–TM6 increases on opening.
4. **Microswitch CVs** for the switches confirmed present (§4.2): ionic-lock
   distance R`3.50`–E`6.30`; NPxxY Y`7.53` χ / TM7 position; W`6.48` χ; Na⁺-pocket
   D`2.50` coordination. Each SHALL be omitted where its switch is absent/static
   in this receptor (LEMMA L2).
5. **Helix tilt/rotation** vs a bundle-fixed axis (e.g. a PCA axis of the TM Cα
   set) or, in explicit membrane, the lipid-phosphate normal. In implicit
   solvent the membrane normal is undefined; a bundle-fixed axis SHALL be used
   instead.

---

## 8. The force-monitoring feedback loop (tie-in, DIAGNOSTIC)

`docs/specs/reaction-force-monitoring.md` streams, per coarse rigid-domain joint,
the **static (`u=0`) mobilizer reaction** (spatial force + torque at the body
origin, in Ground) evaluated on the accepted conformation, at the DCD write
cadence. Its role in world design:

- **Where the field concentrates load flags candidate joints to free.** A welded
  domain boundary carries a full 6-DOF reaction; its reaction torque about a
  candidate torsion axis equals `−∂U/∂φ` for that torsion
  (`reaction-force-monitoring.md` §2.1). A large weld torque about a plausible
  hinge axis flags a joint where added flexibility would relieve real strain — a
  design hint for where to place a §6 hinge.
- **A free domain stuck at ~0 motion flags poor sampling.** A free mobilizer
  transmits ~zero reaction along its free axis by construction, so a fully-free
  joint's reaction is uninformative about *suppressed* motion; conversely, a
  freed domain that shows negligible CV displacement over many rounds flags an
  ineffective world (too-large moving component, or a welded-wall clash), to be
  read from the CV trace, not the force.

This is a DIAGNOSTIC that *informs* world design, not an automatic efficiency
driver. Two properties SHALL be kept in mind:

- **Load ≠ slow mode.** A high reaction magnitude marks where the field pushes,
  not where the free-energy barrier gating a rare event lies. Rare-event gating
  is established by displacement/CV analysis (§7, §10), not force magnitude;
  strain magnitude is a heuristic for *candidate* joints only.
- **The static reaction is velocity-free**, hence a clean deterministic function
  of the accepted conformation `q` (`reaction-force-monitoring.md` §0). A single
  frame's static reaction *is* a legitimate conformational signal, so it can be
  correlated directly against a geometric CV computed on the same frame — no
  thermal averaging needed (this is the reason the design records the static
  `u=0` reaction rather than the instantaneous one).

---

## 9. Sources (external, portable)

- Common class-A activation pathway (CWxP–PIF–Na⁺–NPxxY–DRY chain; TM6 outward /
  TM3-contact switch): Zhou et al., *eLife* 2019 —
  https://elifesciences.org/articles/50279 ;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6954041/
- P`6.50` proline kink as TM6 hinge splitting the helix into two segments:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7484155/ ;
  https://www.sciencedirect.com/science/article/pii/S0021925819722133
- Na⁺ pocket / D`2.50` microswitch coupling:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5810373/
- TM3–TM6 (R`3.50`–TM6) distance as GPCR activation CV; NPxxY-RMSD and ionic-lock
  as MSM order parameters:
  https://pubs.acs.org/doi/10.1021/acs.jctc.5c00600 ;
  https://www.nature.com/articles/s41594-021-00674-7
- Cautionary case (receptor that deletes/immobilizes canonical switches;
  receptor-specific TM6 hinge): FFAR1/GPR40, Ho et al., *PNAS* 2023 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10235965/ ;
  https://www.pnas.org/doi/10.1073/pnas.2219569120 ; and the FFAR1 spec
  `docs/specs/gpcr-world-design/00-ffar1-activation-worlds.md`.

Robosample keys: `spiridon_2020_robosample` (eq:4/6/7/8; worlds; ergodicity),
`spiridon_2017_cdhmc_gibbs` (§3.4, Tables 3–4; MFPT reduction from few-DOF moves).

---

## 10. Correctness conditions

- **PRECONDITION P1 (BW ↔ native numbering).** The `resid`/`resSeq` in the
  prmtop equals the native numbering into which the BW `x.50` anchors are mapped
  for this construct. Guard: spot-check landmark identities (e.g. `x.50`
  positions are the expected conserved residue types). A fusion/truncation that
  shifts numbering shifts every §4–§7 index with it. Runtime guard on the run
  script, not a test.
- **PRECONDITION P2 (welded molecules).** Molecules given `Weld` roots are the
  intended fixed scaffold (nanodisc MSP / a bound transducer), not the receptor
  or mobile lipid. Guard before `initialize`.
- **INVARIANT I1 (ergodicity world present).** Every Gibbs cycle contains ≥1
  fully-flexible world (W0). A cycle without it can leave welded DOFs permanently
  unsampled — a silent sampling bug. Test: over a cycle, every receptor DOF is
  mobile in ≥1 world.
- **INVARIANT I2 (no welded environment adjacent to moving helices).** Lipid/
  solvent retain `Free` roots; coarse helix worlds (W3/W4) do not weld the
  surrounding lipid. Breaking this collapses acceptance to ~0 by the welded-wall
  mechanism.
- **INVARIANT I3 (frame-invariant CVs).** Activation CVs are internal
  distances/angles/dihedrals (§7); no lab-frame RMSD declares a state change
  while the receptor root is `Free`.
- **LEMMA L1 (hinge actuation).** A W3 move that increases the TM3–TM6 opening CV
  must show a correlated change in the TM6 pivot backbone dihedral it actuates.
  Opening CV moving with *no* hinge-dihedral change ⇒ rigid-body drift of the
  outboard component (Claim-C1 artifact), not a hinge swing ⇒ reject the world as
  mis-designed. This is the discriminator that a hinge world does hinge work,
  not drift; it must exercise the opening CV and the hinge dihedral *together*.
- **LEMMA L2 (verified-switch gate).** Any world or CV targeting a conserved
  microswitch (DRY ionic lock, NPxxY Y`7.53`, CWxP W`6.48`, PIF, Na⁺ D`2.50`)
  SHALL be present only if that switch is confirmed to exist and move in this
  receptor (§4.2). Presence of a switch world/CV for a deleted or static switch
  indicates a canonical template was transplanted without per-receptor
  verification. Discriminator: a world built on a dead switch produces no CV
  motion.

---

## 11. Touch list (configuration, not engine)

- Run script (`python/robosample/run.py`, the sole driver): world definitions
  W0–W4 (§6); microswitch list in W2 filtered by §4.2 verification; lipid/solvent
  kept mobile (I2); scaffold welds checked (P2). No engine change.
- CV/analysis: the §7 frame-invariant CV set; drop CVs for switches absent in
  this receptor (L2); use a bundle-fixed axis when no membrane normal exists
  (implicit solvent).
- Conventions at risk: BW ↔ native numbering (P1); inboard/outboard labelling of
  the tree split (Claim C1) determining which component a hinge world moves;
  `JointType` at the hinge (`Torsion` vs `Ball`/`Cylinder`, §6.4);
  membrane-normal vs bundle-axis reference for tilt CVs (§7).
- Optional diagnostic: enable the reporter world (`reaction-force-monitoring.md`)
  on coarse domain boundaries to inform hinge placement (§8).

---

## 12. Verification plan

Discriminating checks (correct design vs plausible-but-non-gating):

1. **Per-world CV response (LEMMA L1 / L2).** Run each coarse world in isolation
   (short single-world trials, small `dt`) and record §7 CVs. Correct: W3/W4
   move the TM3–TM6 opening CV and interhelical distances *with* a correlated
   TM6 pivot dihedral change; W2 moves its microswitch χ / CV. Biased: opening CV
   moves without the hinge dihedral (drift), or a world built on a §4.2 dead
   switch produces no CV motion.
2. **Acceptance tuning (INVARIANT, from paper).** Each world's `dt` / trajectory
   length tuned toward the HMC optimum (≈0.65); a coarse world stuck near 0
   signals a welded-wall clash (I2) or too-large a moving component (Claim C1 —
   shrink the freed hinge / add a repack before it).
3. **Ergodicity spot-check (INVARIANT I1).** Over a full cycle, every receptor
   DOF is mobile in ≥1 world (W0 covers all). A DOF welded in every world is
   unsampled.
4. **Transition oracle (integration test).** From the inactive structure, the
   composed cycle SHOULD, over long runs, increase the primary TM3–TM6 opening CV
   toward the active basin more often than a fully-flexible-only control in equal
   wall-clock — the receptor analogue of the ala-dipeptide α_L-basin MFPT
   reduction (`spiridon_2020_robosample`; `spiridon_2017_cdhmc_gibbs` §3.4). This
   is the end-to-end discriminator that the world composition accelerates the
   activation-relevant rare event, not just local fluctuation. The expected
   relation is *shorter MFPT to the open basin than the flexible-only control*,
   not an absolute rate.

---

## 13. Open questions / blocking unknowns (portable)

- **Q1 (root override) — not user-selectable today.** The BFS root is the
  auto-chosen heaviest massive leaf. The §6 baseline worlds are correct
  regardless, but which side of a hinge is fixed vs moving is not a config knob.
  If a hinge world's acceptance proves root-limited (L1 drift dominates),
  exposing a root override is a future engine task, not a configuration change.
- **Q2 (endpoint CV values).** The active/inactive numeric values of the primary
  TM3–TM6 CV are receptor-specific and set the basin boundary in check 4. Obtain
  from the receptor's own inactive and active structures; until then use the
  *relative* increase from the run's starting frame.
- **Q3 (which switches are live).** The W2 residue list and the §7 microswitch
  CVs depend on which canonical switches this receptor retains and moves (§4.2).
  Resolve from two endpoint structures or a sequence/structure survey before
  building W2; not blocking for W0/W1/W3/W4.
- **Q4 (loaded state).** World emphasis differs for apo vs orthosteric-agonist vs
  allosteric-modulator runs; a pocket/allosteric world is meaningful only if that
  site is occupied. Needed to finalize the pocket-χ portion of W2. Not blocking
  for the backbone/rigid worlds.
