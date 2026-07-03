# FFAR1 (class-A GPCR) world design for sampling activation-relevant transitions

Status: research spec. Read-only design guidance for Robosample `Context` /
world configuration (Python run scripts + CV analysis). No engine (C++/CUDA)
change is implied.

Anchoring references:
- `spiridon_2020_robosample` (primary Robosample paper: worlds = rigid bodies +
  joint types over one molecular graph; Gibbs sampling; Fixman correction;
  "define rigid bodies by secondary structure or domains; sample torsions more
  than bond lengths/angles; include at least one fully flexible world per cycle").
- `spiridon_2017_cdhmc_gibbs` (CDHMC-as-Gibbs-move; Ramachandran dynamics with
  only φ/ψ mobile substantially shortens the mean-first-passage time to the
  rarest transition — §3.4, Tables 3–4; qualitative MFPT reduction, not a
  specific factor).
- Codebase API: `python/robosample/context.py`, `python/robosample/robo_bindings.pyi`
  (`JointType` enum, `build_flexibilities`, `add_torsional_world`,
  `add_cartesian_world`, `set_root_mobility`), and the existing
  `python/robosample/run_ffar1.py` / `analyse_ffar1.py` scaffold.
- External structural biology cited inline with URLs in §7.

---

## 1. Problem restatement

**User phrasing.** "How do I design a representative Robosample `world` (Gibbs
block = a chosen set of flexible mobilizers over the kinematic tree) to sample
FFAR1 so that sampling actually explores the functionally relevant conformational
transitions?"

**Restatement in codebase vocabulary.** FFAR1 is loaded as one robot: a single
BFS-rooted kinematic tree (`atoms_root_index` = molecule compound index 0;
`context.py` §"Per-molecule root") whose bodies are joined by mobilizers. A
*world* selects a `JointType` for a chosen set of bonds (`build_flexibilities`);
every unselected non-terminal, non-ring bond is `Rigid` (welded). Blocked Gibbs
sampling cycles worlds; each world's HMC move draws momenta and integrates only
that world's mobile DOFs, then Metropolis-Hastings-accepts against the full
OpenMM energy (plus the Fixman correction for the reduced-coordinate marginal).

The design lever is therefore: **which bonds get which mobility in each world,
and how worlds compose across a Gibbs cycle**, chosen so that the *slow*
collective coordinates that gate FFAR1 activation are proposed in few-DOF,
large-amplitude moves (the mechanism by which Ramachandran dynamics beat
fully-flexible sampling in `spiridon_2017_cdhmc_gibbs`), while a fully-flexible
world preserves ergodicity.

The underlying problems to solve, in order:

1. **Identify FFAR1's slow activation coordinates** as concrete bonds/atoms
   (§3), *not* by analogy to canonical class-A switches — FFAR1 lacks several of
   them (§3.2).
2. **Map each coordinate to a mobility class** — backbone torsion (φ/ψ),
   sidechain torsion (χ), or a welded-helix rigid-body split — respecting the
   fact that a tree move splits the receptor into exactly two rigid components
   (§2, Claim C1).
3. **Compose a minimal set of worlds** that together make the activation path
   accessible while remaining ergodic and giving usable acceptance (§4).
4. **Define frame-invariant CVs** to confirm a world moved the intended state
   (§5), and **environment constraints** for the nanodisc/implicit case (§6).

**Reading chosen.** "Functionally relevant transitions" is read as *the
apo/inactive ↔ agonist/active transition of the receptor 7TM bundle* (TM6
cytoplasmic movement, TM4–TM5 initiation, ligand-pocket opening), not
ligand (un)binding kinetics and not lipid dynamics. Rigid-body coarse-graining
cannot reproduce kinetics anyway (`spiridon_2020_robosample` §4 limitation); the
goal is basin-to-basin *configurational* sampling.

---

## 2. Load-bearing kinematic fact (why "clean" helix swings are not free)

**Claim C1 (tree split).** FFAR1 is one spanning tree rooted at a single fixed
atom (BFS compound index 0). Freeing a set of backbone joints and welding the
rest partitions the receptor into exactly two maximal rigid components per
freed edge: everything *inboard* (root side) of the edge stays fixed relative to
the base body; everything *outboard* (tip side) moves rigidly. There is **no
tree move that displaces an interior helix segment while holding both of its
covalent flanks fixed** — that would require a cycle, and the bundle's mutual
packing (TM6 against TM3/TM5/TM7) is non-covalent, hence absent from the tree
and enforced only energetically by OpenMM.

Consequences that constrain every design below:

- A "TM6 outward swing" world cannot isolate TM6: whichever backbone hinge is
  freed, the entire outboard component rides along. The design job is to choose
  the hinge so the two components are a *meaningful* collective mode and the
  moving component is small enough to accept (§4.2).
- The root atom is fixed by BFS at compound 0 (≈ chain start), and is **not**
  currently a per-residue user choice (see Open Question Q1). So one cannot
  re-root to make the cytoplasmic helix ends the moving tips; the moving
  component is determined by the chain threading TM1→…→TM7 and by which bonds
  are welded.
- Welding a helix interior (not selecting its bonds in a torsional world) is how
  a helix becomes a rigid body; this is `spiridon_2020_robosample`'s "define
  rigid bodies by secondary structure." The existing `run_ffar1.py` World 4
  already does this (welds all seven TM cores, frees all loop φ/ψ).

---

## 3. FFAR1 activation mechanism mapped to atoms

Ballesteros–Weinstein (BW) `x.50` anchors: TM1 N1.50, TM2 D2.50, TM3 R3.50,
TM4 W4.50, TM5 P5.50, TM6 P6.50, TM7 P7.50. Residue numbers below are FFAR1
(human) native numbering as used in the PNAS 2023 structure paper (§7); the
existing `analyse_ffar1.py` `HELIX_RANGES` (TM6 = 223–248, etc.) already use the
same scale — see PRECONDITION P1 for the required verification.

### 3.1 What moves on activation (canonical class-A backbone)

- **TM6 cytoplasmic movement** is the largest activation motion; its cytoplasmic
  end swings away from TM3, opening the intracellular G-protein cavity. In
  canonical receptors the pivot is the P6.50 (CWxP) proline kink, which splits
  TM6 into two segments and amplifies the swing (§7). TM3↔TM6 contact is
  replaced by TM3↔TM7; TM5 and TM7 rearrange toward TM3.

### 3.2 FFAR1-specific deviations (do **not** transplant canonical switches)

FFAR1 lacks or immobilizes several canonical microswitches — designing worlds
around them samples DOFs that do **not** gate FFAR1 activation:

- **No DRY / no ionic lock.** FFAR1 has G3.49-R3.50-Y3.51 and **K6.30** (not
  Glu) at 6.30, so the R3.50–E6.30 salt bridge does not exist. The receptor is
  energetically pre-sensitized toward the active state. → The standard "ionic
  lock" CV and any world built to break it are **inapplicable**.
- **Altered NPxxY.** TM7 carries **NPLVT (272–276)**; the conserved Y7.53 is
  absent and TM7 is disordered / not part of the Gq interface. → No NPxxY-Tyr
  rotamer world; no NPxxY-RMSD CV.
- **Immobile toggle.** **V237 (6.48)** occupies the W6.48 toggle position and
  stays put between inactive and active. → No CWxP toggle world.
- **Static connector/PIF.** P194(5.50)/G94(3.40)/L233(6.44) do not repack much.
  → PIF is a poor lever for FFAR1.
- **Reduced-amplitude TM6 pivot at a diglycine, not (only) the proline kink.**
  Active-state R104(3.50)Cα–A222(6.33)Cα ≈ **10.6 Å** (vs 14.6 Å for β2AR-Gs),
  a *smaller* opening, facilitated by a **bend at the diglycine G227(6.38)-
  G228(6.39)**, which is *cytoplasmic to* the canonical P6.50 kink. → The FFAR1
  TM6 hinge world should free backbone φ/ψ at/around **226–229 (the diglycine)**
  plus the cytoplasmic ICL3/TM6 base, not the mid-helix toggle region.

### 3.3 What FFAR1 actually uses (the real levers)

- **TM4–TM5 initiation microswitch** (the trigger): fatty-acid binding between
  TM3–TM4 reorients L144(4.63)↔S178(5.34), V141(4.60)↔A182(5.38),
  H137(4.56)↔L186(5.42); TM5 then shifts toward TM6.
- **Core polar network**: R183(5.39), R258(7.35), N244(6.55), E172(ECL2),
  F87(3.33), W174(ECL2) coordinate the ligand carboxylate.
- **Ligand entry**: fatty acids enter the TM3–TM4 orthosteric site *through the
  membrane*, gated by "fluctuations in the position of the extracellular half of
  TM3."
- **Two agonist sites**: orthosteric (TM3–TM4 lipid interface; fatty acids,
  partial agonist TAK-875) and an inner-leaflet allosteric AgoPAM site at ICL2 +
  lipid faces of TM3/4/5 (full agonists AP8/compound 1).

### 3.4 Mechanism → mobility class

| Motion | Best mobility representation | Rationale |
|---|---|---|
| TM6 cytoplasmic swing | welded TM6 core + freed backbone φ/ψ at diglycine 226–229 &amp; ICL3 (Claim C1: two-component split) | it is a rigid-segment pivot about a backbone hinge |
| TM5/TM7 rearrangement toward TM3 | freed loop φ/ψ flanking welded TM5, TM7 cores | collective helix repositioning |
| TM4–TM5 initiation repack | sidechain χ of the interface residues (+ optional local φ/ψ) | rotamer/packing change, not backbone |
| Core polar network reorientation | sidechain χ of R183, R258, N244, E172, F87 | rotamer toggles |
| Ligand-pocket / TM3-EC fluctuation | φ/ψ of the extracellular half of TM3 + pocket χ | local backbone breathing + rotamers |

---

## 4. World designs (composable Gibbs cycle)

Notation for a world: `{selection → JointType}`. A torsional world assigns
`Torsion` (1-DOF pin) to the selected bonds; all other eligible bonds are
`Rigid`. `add_torsional_world(selection)` builds it; `build_flexibilities` maps a
bond selection (filter `standard_dihedral_bonds` by `dihedral_type` and `resid`)
to that selection. `Ball`/`Cylinder`/`SphericalCoords` are available if a
richer joint is wanted at a hinge (see §4.2 note).

The receptor SHALL be sampled by a **Gibbs cycle of the worlds below**, and the
cycle SHALL include one ergodicity world (W0). Worlds are ordered here by scope,
not by execution order.

### W0 — Ergodicity / fully-flexible (required)

`{all eligible bonds → Cartesian}` via `add_cartesian_world()`, or a
fully-flexible torsional world. Purpose: guarantee every DOF is reachable
(`spiridon_2020_robosample` §2.1.1 condition 2) and let the bundle relax the
junction strain that the coarse worlds inject. Short trajectories, small `dt`.

### W1 — All-sidechain repack

`{χ1..χ5 of all receptor residues (+ lipid χ) → Torsion}`. Purpose: rotamer
repacking that must accompany any backbone move for the bundle to re-form
contacts; cheap and high-acceptance. This is `run_ffar1.py` World 2b, kept.
Moves no CV by itself but unblocks the backbone worlds.

### W2 — Ligand-pocket / core-network world

`{χ of {F87, E172(ECL2), R183, N244, R258, W174} + χ of TAK-875/AgoPAM pocket
residues + φ/ψ of the extracellular half of TM3 → Torsion}`. Purpose: sample
the orthosteric carboxylate network and the TM3-EC breathing that gates ligand
accommodation. CV moved: TM3–TM4 pocket geometry, core-network contacts.

### W3 — TM4–TM5 initiation-microswitch world

`{χ (+ local φ/ψ) of {H137, V141, L144, S178, A182, L186} → Torsion}`. Purpose:
the FFAR1 activation *trigger* (§3.3). CV moved: TM4–TM5 interface distances
(H137–L186, V141–A182), TM5→TM6 shift.

### W4 — TM6 diglycine-hinge world (the activation swing)

`{φ/ψ of residues ≈ 209–229 (ICL3 + TM6 cytoplasmic base + diglycine
G227-G228) → Torsion}`, with the TM6 core (≈230–248) and all other TM cores
**welded**. Purpose: propose the TM6 cytoplasmic opening in a few-DOF,
large-amplitude move about the mechanistically correct FFAR1 hinge.

Design notes tied to Claim C1:
- Because the chain threads TM5→ICL3→TM6(cyto→EC)→ECL3→TM7, freeing the diglycine
  hinge moves the *outboard* component (everything tipward of the hinge). This
  world SHOULD be read as "open the cytoplasmic face," not "move TM6 alone."
- The existing `run_ffar1.py` World 3c frees 235–242 (around V237/6.48, the
  *immobile* toggle region) — that samples a non-gating DOF for FFAR1 and SHOULD
  be replaced by (or supplemented with) the 226–229 diglycine hinge. World 3b
  (209–223) already covers the ICL3/TM6 base and SHOULD be kept.

### W5 — Rigid-helix / flexible-loop world (global bundle repacking)

`{φ/ψ of all non-helix (loop) residues → Torsion}`, TM cores welded. This is
`run_ffar1.py` World 4, kept: it lets all seven helices reorient as rigid bodies
connected by flexible loops — the coarse, large-step complement to W0.

### 4.1 Composition rationale (from theory, not intuition)

`spiridon_2017_cdhmc_gibbs`/`spiridon_2020_robosample` show the efficiency gain
comes from *few mobile DOFs taking large steps* (Ramachandran dynamics markedly
shorten the mean-first-passage time to the rarest transition; §3.4, Tables 3–4),
while a fully-flexible world restores ergodicity.
W3→W4→W1→W5→W0 mirrors the mechanistic causal chain (trigger → hinge → repack →
global relax → ergodic mix): each coarse proposal (W3/W4/W5) is followed by
repacking (W1) and a flexible relax (W0) so junction strain from Claim-C1 splits
is annealed before the next coarse proposal. Per-world `dt`/trajectory length
SHALL be tuned to acceptance ≈ 0.651 (`spiridon_2020_robosample` §2.1.2).

### 4.2 NOTE on richer hinge joints

At the diglycine hinge a `Ball` (3-rotational) or `Cylinder` (rotation+slide)
mobilizer would sample the bend more completely than a single `Torsion`, at the
cost of a larger Fixman correction and smaller stable `dt`. This is optional and
orthogonal to correctness; start with `Torsion` (matches the paper's default and
the existing scaffold) and only escalate if the hinge CV is under-sampled.

---

## 5. Collective variables / order parameters (world validation)

CVs SHALL be **frame-invariant** (internal distances/angles), because the
receptor root is `Free` (6-DOF) and the whole robot translates/rotates in the
box; a lab-frame RMSD would confound rigid-body drift with activation.
`analyse_ffar1.py` already computes Cα-centroid interhelical distances and
per-helix tilt/rotation vs a lipid-phosphate membrane normal — reuse that path.

FFAR1 CV set (with the canonical ones that must be dropped, per §3.2):

1. **TM6 opening — primary.** R104(3.50)Cα – A222(6.33)Cα distance. Active FFAR1
   ≈ 10.6 Å; a productive W4/W5 move increases this toward the active value
   relative to the inactive (TAK-875) reference. Use *this*, **not** the
   R3.50–E6.30 ionic-lock distance (no E6.30 in FFAR1).
2. **TM6 hinge dihedral.** Backbone φ/ψ at G227/G228 — the coordinate W4 directly
   actuates; confirms the hinge, not just downstream drag, moved.
3. **Interhelical centroid distances** TM3–TM6, TM5–TM6, TM3–TM5 (already in
   `analyse_ffar1.py PAIRS`): TM3–TM6 increases and TM5–TM6 changes on opening.
4. **TM4–TM5 initiation** distances H137–L186, V141–A182 (W3 target).
5. **TM6 rotation/tilt** vs membrane normal (already computed) — the cytoplasmic
   register change.

Drop for FFAR1: ionic-lock (R3.50–E6.30), NPxxY-Tyr rotamer/RMSD, CWxP-toggle
χ — the residues are absent or static.

---

## 6. Environment / membrane caveats (nanodisc, implicit solvent)

- **Absolute reference frame.** In `run_ffar1.py`, molecules 1–2 are `Weld`ed to
  ground and molecule 0 (FFAR1) is `Free`. The welded molecules provide a fixed
  lab frame; this is fine because §5 CVs are internal (frame-invariant) anyway.
  If those welded molecules are the nanodisc scaffold (MSP belts), welding them
  is appropriate (they are not the sampling target). PRECONDITION P2 checks their
  identity.
- **Lipids must stay mobile.** Lipids/solvent SHALL keep `Free` roots and be
  given sidechain/tail torsional freedom (W1's lipid χ). A helix rigid-body move
  that swings into *welded* lipid would be a hard wall and reject — the same
  failure mode diagnosed for welded-environment NCMC (see
  `docs/specs/ncmc-explicit-solvent/`). Welding lipid around a moving helix is a
  correctness-relevant mistake, not a tuning choice.
- **Implicit-solvent variant.** With no box, GBSA-OBC2 and no explicit membrane;
  then the lipid-phosphate membrane-normal CV in `analyse_ffar1.py` is undefined
  and helix tilt/rotation must be referenced to a bundle-fixed axis (e.g. TM
  Cα-PCA of the whole bundle) instead of the phosphate leaflets.
- **Ground-frame stability.** Helix motions are only physical relative to a
  stable bundle core; W4/W5 rely on the TM cores being welded and on W0/W1
  annealing junction strain, per §4.1.

---

## 7. Derivation sketch and external sources

The efficiency argument is `spiridon_2017_cdhmc_gibbs` / `spiridon_2020_robosample`:
the reduced-coordinate marginal is corrected by the Fixman potential
U' = kT·ln(|M_tot|/|M|)^{1/2} (eq. 7), and few-DOF worlds take large steps that
cross barriers a fully-flexible world crosses far more slowly (Ramachandran
result: a qualitative MFPT reduction, §3.4 Tables 3–4 — not a specific factor).
Claim C1 is elementary graph theory on the BFS spanning tree
(`context.py`, `atoms_root_index`/`root_mobilities`). The structural-biology
mapping (§3) is external:

- Common class-A activation pathway (CWxP–PIF–Na⁺–NPxxY–DRY chain; TM6 outward /
  TM3 contact switch): Zhou et al., *eLife* 2019 —
  https://elifesciences.org/articles/50279 ;
  https://pmc.ncbi.nlm.nih.gov/articles/PMC6954041/
- P6.50 proline kink as TM6 hinge splitting the helix into two segments:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC7484155/ (rhodopsin proline hinge);
  https://www.sciencedirect.com/science/article/pii/S0021925819722133 (β2AR
  proline-kink rotamer toggle).
- Na⁺ pocket / D2.50 microswitch coupling:
  https://pmc.ncbi.nlm.nih.gov/articles/PMC5810373/
- FFAR1/GPR40 molecular mechanism, diglycine hinge, absent DRY/NPxxY, V6.48
  toggle static, TM4–TM5 initiation, dual agonist sites, R104(3.50)–A222(6.33)
  10.6 Å: Ho et al., *PNAS* 2023 —
  https://pmc.ncbi.nlm.nih.gov/articles/PMC10235965/ ;
  https://www.pnas.org/doi/10.1073/pnas.2219569120
- FFAR1 lacks conserved DRY/NPXXY; two allosterically coupled sites; TAK-875
  lipid-facing pocket between TM3/TM4:
  https://pubmed.ncbi.nlm.nih.gov/32851580/ ;
  https://www.pnas.org/doi/pdf/10.1073/pnas.2219569120
- TM3–TM6 (R3.50–6.x) distance as GPCR activation CV; NPxxY-RMSD and ionic-lock
  as MSM order parameters (canonical, contrasted here for FFAR1):
  https://pubs.acs.org/doi/10.1021/acs.jctc.5c00600 ;
  https://www.nature.com/articles/s41594-021-00674-7

---

## 8. Correctness conditions

- **PRECONDITION P1 (numbering).** The `resid` column of
  `standard_dihedral_bonds` and `resSeq` in the prmtop equal the FFAR1 native
  numbering in which BW positions are expressed. Guard: spot-check landmarks —
  `resid 104` is Arg (R3.50), `resid 222` is Ala (A6.33), `resid 227` and `228`
  are Gly-Gly (the diglycine hinge), `resid 237` is Val (V6.48). If a construct
  fusion/truncation shifts numbering, every §3–§5 residue index shifts with it.
  This is a runtime guard on the run script, not a test.
- **PRECONDITION P2 (welded molecules).** Molecules given `Weld` roots
  (currently indices 1–2) are the intended fixed scaffold (nanodisc MSP / a bound
  transducer), not the receptor or mobile lipid. Guard before `initialize`.
- **INVARIANT I1 (ergodicity world present).** Every Gibbs cycle contains ≥1
  fully-flexible world (W0). A cycle without it can leave welded DOFs
  permanently unsampled — a silent sampling bug.
- **INVARIANT I2 (no welded lipid adjacent to moving helices).** Lipid/solvent
  molecules retain `Free` roots; the coarse helix worlds (W4/W5) do not weld the
  surrounding lipid. Breaking this collapses acceptance to ~0 by the
  welded-wall mechanism.
- **INVARIANT I3 (frame-invariant CVs).** Activation CVs are internal
  distances/angles (§5); no lab-frame RMSD is used to declare a state change
  while the receptor root is `Free`.
- **LEMMA L1 (hinge actuation).** A W4 move that increases the TM6-opening CV
  (R104Cα–A222Cα) must show a correlated change in the G227/G228 backbone
  dihedral it actuates. If the opening CV moves with *no* hinge-dihedral change,
  the displacement is pure rigid-body drift of the outboard component (Claim C1
  artifact), not a hinge swing — reject the world as mis-designed.
- **LEMMA L2 (canonical-switch null).** Any CV or world targeting R3.50–E6.30,
  NPxxY-Y7.53, or the W6.48 toggle SHALL be absent from the FFAR1 configuration;
  their presence indicates a canonical-GPCR template was transplanted without the
  FFAR1 corrections of §3.2.

---

## 9. Touch list

- `python/robosample/run_ffar1.py` — world definitions: replace/augment the
  mid-TM6 world (currently 235–242) with the 226–229 diglycine hinge (W4); add
  the TM4–TM5 initiation world (W3) and the core-network world (W2); keep W0/W1/
  W5. No engine change.
- `python/robosample/analyse_ffar1.py` — CV set: add R104Cα–A222Cα and the
  G227/G228 hinge dihedral; the ionic-lock/NPxxY interpretation strings in
  `_SECTION_TEXT` (e.g. "TM6-TM7 = NPxxY-motif region", "TM6 rotation …
  DRY-motif exposure") are canonical-GPCR boilerplate that is **false for FFAR1**
  and should be corrected (Lemma L2).
- Conventions at risk: BW↔native residue numbering (P1); inboard/outboard
  orientation of the tree split (Claim C1) determining which component the
  diglycine world actually moves; `JointType` choice at the hinge (Torsion vs
  Ball, §4.2); membrane-normal reference validity in the implicit-solvent variant
  (§6).

---

## 10. Verification plan

Discriminating checks (a correct design vs a plausible-but-non-gating one):

1. **CV response, per world (LEMMA L1 / L2).** Run each coarse world in isolation
   (short single-world trials, `dt`=1 fs, as in `spiridon_2020_robosample`
   §2.1.2 workflow) and record the §5 CVs. Correct: W4/W5 move
   R104Cα–A222Cα and TM3–TM6 with a correlated G227/G228 dihedral change; W3
   moves H137–L186 / V141–A182. Biased/mis-designed: the opening CV moves without
   the hinge dihedral (drift), or a world built on a §3.2 dead switch produces no
   CV motion at all.
2. **Acceptance tuning (INVARIANT, from paper).** Each world's `dt`/trajectory
   length is tuned to acceptance ≈ 0.651; a coarse world stuck near 0 acceptance
   signals a welded-wall clash (check I2) or too-large a moving component
   (Claim C1 — shrink the freed hinge / add W1 repack before it).
3. **Ergodicity spot-check (INVARIANT I1).** Over a full cycle, confirm every
   receptor DOF is mobile in ≥1 world (W0 covers all). A DOF welded in every
   world is unsampled.
4. **Transition oracle (integration test).** Starting from the inactive
   (TAK-875-bound) structure, the composed cycle should, over long runs, increase
   the TM6-opening CV toward the active ~10.6 Å basin more often than a
   fully-flexible-only control in equal wall-clock — the FFAR1 analogue of the
   ala-dipeptide αL-basin MFPT reduction (`spiridon_2020_robosample` §3.1).
   This is the end-to-end discriminator that the world composition actually
   accelerates the activation-relevant rare event, not just local fluctuation.

---

## 11. Open questions / blocking unknowns

- **Q1 (root atom selection) — RESOLVED: not user-selectable today.** The BFS
  root is chosen automatically as the **heaviest real (massive) leaf** of the
  molecule's spanning tree (`molecule_prototype.py:94-114`: must be a massive,
  degree-1 atom so the Z-matrix root triplet forms; ties broken by mass).
  `MoleculePrototype.__init__` takes only `(molecule, dihedral_classifier)` — no
  root argument — and the whole compound/BFS ordering, Z-matrix, and
  `atoms_root_index` derive from that auto-chosen leaf. The
  `system_topology.atoms_root_index` setter (`robo_bindings.pyi:470`) only stores
  a global index into an ordering already computed from compound-0; it does not
  re-run BFS from a different atom. **Consequence:** the root cannot be placed at
  the ECL2 disulfide (or any chosen anchor) without an engine/API change, so the
  "optimal" hinge variant (cytoplasmic tips as the moving component) is not a
  config knob. The §4 baseline worlds are correct regardless of root — Claim C1's
  two-component split still holds — but which side of the diglycine hinge is
  "fixed" vs "moving" is determined by the auto-selected leaf + the TM1→…→TM7
  chain threading, not chooseable. Exposing a root override is a possible future
  engine task if W4 acceptance proves root-limited (see Lemma L1 drift check).
- **Q2 (FFAR1 inactive reference CV value).** The active R104Cα–A222Cα ≈ 10.6 Å
  is cited; the inactive-structure value (TAK-875 complex) is needed to set the
  numeric basin boundary in check 4. Obtain from the inactive FFAR1 structure
  (e.g. PDB 4PHU) before hard-coding a threshold; until then, use the *relative*
  increase from the run's own starting frame.
- **Q3 (which agonist/state is loaded).** The world emphasis differs for an
  apo run vs a TAK-875 (partial, orthosteric) vs AP8 (full, allosteric) run:
  the AgoPAM world (W2 allosteric variant) is only meaningful if the inner-leaflet
  site is occupied. Needed to finalize W2's residue list. Not blocking for
  W0/W1/W3/W4/W5.
