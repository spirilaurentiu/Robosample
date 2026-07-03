# Monitoring selections for β-adrenergic receptors (reaction-force / velocity reporter)

Status: selections/design spec for the three `examples/febs` β-adrenergic models. Read-only
guidance for configuring the Robosample reaction-force reporter world (via `python/robosample/run.py`)
and CV analysis. **No engine (C++/CUDA/Python-library) change is implied or designed here.**

Applies the method of the conceptual class-A spec; does not re-derive it. Read first:
- `docs/specs/gpcr-world-design/10-classA-world-selection-conceptual.md` (method, Claim C1 tree-split,
  BW→native discipline P1, frame-invariant CVs, force-monitoring tie-in §8).
- `docs/specs/gpcr-world-design/00-ffar1-activation-worlds.md` (cautionary counter-example: a receptor
  that *deletes* canonical switches — the opposite of the three receptors here).
- `docs/specs/reaction-force-monitoring.md` (what the reporter records: static `u=0` mobilizer
  reaction, force+torque in Ground, per interesting joint, at DCD cadence; sparse domain-boundary
  joints; self-derived from a flagged world's flexed bodies).

**All residue numbers in this spec are in the MODEL (prmtop) numbering of the `examples/febs` files**
— NOT deposited-PDB/UniProt numbering. §2 gives the exact mapping and how it was verified.

---

## 1. Purpose

Give ready-to-use monitoring selections for the three active-state β-adrenergic receptor models in
`examples/febs/`, so the reaction-force reporter can print the per-domain spatial load of each
receptor. Per the user's choice, **each model gets two reporter worlds** (§5): a **coarse**
boundary-monitor world (inter-helix loads) and a **focused** switch-monitor world (the activation
machinery). "Spatial velocity" monitoring uses the *same* selections but is not yet emitted by the
engine (§8, Q2).

---

## 2. The `examples/febs` models — ground truth (verified against the prmtops)

Each `<PDB>.{lig,noLig}.nanodisc.prmtop` is **receptor-only in a nanodisc** — there is no G protein,
no T4-lysozyme fusion, no nanobody. Molecule/segment layout (verified with ParmEd):

| molecule_idx | contents |
|---|---|
| **0** | the receptor (ACE cap + chain + NME cap) — the ONLY monitored molecule |
| 1, 2 | two MSP nanodisc-belt copies (apolipoprotein scaffold) |
| 3 … | membrane lipids |
| last | the agonist (only in the `.lig` files; absent in `.noLig`) |

Receptor chain is **contiguous** (loops modeled; no internal gaps), renumbered from 1:

| Model | Receptor | Model resid range | Offset (model = auth − OFF) | Agonist |
|---|---|---|---|---|
| **3SN6** | human β2AR (ADRB2, P07550) | 1–312 | **29** | BI-167107 |
| **7JJO** | **turkey** β1AR (ADRB1, P07700) | 1–318 | **39** | isoproterenol |
| **7DH5** | **dog** β3AR (ADRB3, O02662) | 1–325 | **35** | mirabegron |

**How the mapping was verified.** For each model, the conserved motifs DRY (`DRY`), CWxP (`CW.P` →
W6.48/P6.50), and NPxxY (`NP..Y` → Y7.53) were located directly in the prmtop sequence; the offset
they imply is identical across TM3 (pre-ICL3) and TM6/TM7 (post-ICL3), and every BW anchor
(D2.50, R3.50, E6.30, W6.48, P6.50, Y7.53) was confirmed to carry the expected residue identity at
`auth − OFF`. Because the offset is constant across the ICL3 region and the chain is contiguous, the
**ICL3/ECL2 gaps present in the deposited PDBs are NOT present in these models** — the loops were
built. NOTE: built loops may be lower-confidence geometry than the resolved TM core; interpret forces
on loop/ICL3 joints with that in mind.

NOTE (numbering used downstream): the model resid here is the ParmEd residue number (1..N over the
receptor). It is for human readability and for the §4 helper; the actual selection is by atom index
(§3), so it is robust to any off-by-one in the resid convention.

---

## 3. Robosample selection mechanics (read before using the tables)

The current selection surface constrains HOW these selections are expressed:

- **`standard_dihedral_bonds` (df_bonds) has NO usable residue column.** `residue_idx` is a hardcoded
  `-1` placeholder and `residue_name` is `"UNK"` (`context.py:393-394`). Filtering bonds by residue
  through the DataFrame does not work. (This is also why the legacy `run_ffar1.py` `["resid"]` filter
  is dead — that script is out of scope.)
- **Two working selection paths, both through `build_flexibilities`, neither needing a source change:**
  1. **By secondary structure (`dss` column)** — populated by real MDTraj DSSP
     (`molecule_prototype.py:242`). Use it for the **coarse world**: keep helix bonds rigid, free the
     rest. The reporter then self-derives the inter-helix boundary joints.
  2. **By explicit global atom-index `(i,j)` pairs** — `build_flexibilities` accepts an iterable of
     global atom-index bond pairs directly (`context.py:79,96`). Use it for the **focused world**:
     map a model-resid list → that residue's backbone/χ bonds → `(i,j)` pairs. §4 gives a helper that
     does this from the prmtop.
- **Receptor = molecule 0**, whose atoms are the first block (global atom index == prmtop atom index),
  so the atom-index mapping is direct.
- The reporter reads the reaction at whichever joints the flagged world **flexes**; a flexed joint
  transmits the load in its five constrained components (≈0 only about its own free axis), which is
  the informative per-domain reading (`reaction-force-monitoring.md` §2.1).

---

## 4. Selection helper (put under `examples/febs/`, not in the library)

A standalone helper that turns a model-resid list into the atom-index bonds `build_flexibilities`
wants, and the coarse helix-weld set from `dss`. It reads only the prmtop; it is NOT a Robosample
source file.

```python
# examples/febs/febs_select.py  — selection helper, not a library module.
import parmed as pmd

BACKBONE = {"phi": ("N", "CA"), "psi": ("CA", "C")}

def receptor_residues(prmtop_path):
    """Molecule-0 (receptor) residues as (model_resid, ParmEd residue). Skips leading caps."""
    s = pmd.load_file(prmtop_path)
    out, started = [], False
    for r in s.residues:
        aa = len(r.name.strip()) == 3 and r.name.strip() not in ("ACE", "NME", "NHE")
        if aa:
            started = True
            out.append((r.number, r))
        elif started:
            break  # first non-AA after the chain => receptor ended
    return out

def backbone_bond_pairs(prmtop_path, model_resids):
    """Global atom-index (i,j) pairs for the N-CA (phi) and CA-C (psi) bonds of the given residues."""
    pairs = []
    for resid, r in receptor_residues(prmtop_path):
        if resid not in model_resids:
            continue
        idx = {a.name: a.idx for a in r.atoms}   # a.idx is the GLOBAL prmtop atom index
        for names in BACKBONE.values():
            if names[0] in idx and names[1] in idx:
                pairs.append((idx[names[0]], idx[names[1]]))
    return pairs

def sidechain_bond_pairs(prmtop_path, model_resids):
    """Global (i,j) pairs for rotatable sidechain bonds (all intra-residue bonds beyond CA-CB)."""
    pairs = []
    for resid, r in receptor_residues(prmtop_path):
        if resid not in model_resids:
            continue
        ridx = {a.idx for a in r.atoms}
        for b in r.atoms[0].residue.bonds if hasattr(r, "bonds") else []:
            pass  # see NOTE below
    # NOTE: sidechain-chi enumeration is model-dependent; simplest is to pass these residues
    # through the df_bonds dihedral_type=='chi' rows filtered by atom membership in `ridx`
    # (atom1_idx is molecule-local == global for molecule 0). Left explicit rather than guessed.
    return pairs
```

For the coarse world, filter df_bonds on `dss` instead of atom indices (keep helix bonds rigid).
NOTE: the `sidechain_bond_pairs` body is intentionally left as a documented approach rather than a
guessed χ-atom enumeration — resolve χ bonds via the `dihedral_type == "chi"` rows of
`standard_dihedral_bonds` filtered by atom membership in the target residues, since those rows already
encode the rotatable χ bonds.

---

## 5. Per-structure selections (MODEL numbering)

For each model: (A) identity, (B) seven TM-helix cores to weld, (C) functional anchors, (D) the two
reporter worlds, (E) frame-invariant CVs. Ranges/residues are model resid; `auth` is shown once in §2.

### 5A. 3SN6 — human β2AR · agonist BI-167107 · active (offset 29)

**(B) TM cores (weld; model resid):**
TM1 1–27 · TM2 42–67 · TM3 75–100 · TM4 121–142 · TM5 169–191 · TM6 238–269 · TM7 276–297.

**(C) Functional anchors (BW → model resid):**

| BW | element | model resid |
|---|---|---|
| D2.50 | Na⁺ pocket | **50** |
| D3.32 | ligand amine anchor | **84** |
| I3.40 | PIF | **92** |
| D3.49/R3.50/Y3.51 | DRY | **101 / 102 / 103** |
| S5.42/S5.43/S5.46 | catechol serines | **174 / 175 / 178** |
| P5.50 | PIF | **182** |
| E6.30 | ionic-lock partner | **239** |
| F6.44 | PIF | **253** |
| C6.47/W6.48/P6.50 | CWxP toggle | **256 / 257 / 259** |
| N6.55 | pocket | **264** |
| N7.49/P7.50/Y7.53 | NPxxY | **293 / 294 / 297** |

**(D) Reporter worlds.**
- **Coarse (inter-helix loads):** weld the seven TM cores above (or keep all helix bonds rigid via
  `dss`), free the loops; flag as reporter. Self-derived monitors = the six cytoplasmic boundary
  joints: TM6 base **236–239**, TM5 base **202–210**, TM3/DRY base **103–107**, TM7/NPxxY base
  **297–300**, TM2 base **32–37** (reference), TM4 base **109–113** (reference).
- **Focused (activation switches):** free only — TM6 cytoplasmic hinge φ/ψ **236–239** (primary),
  and the χ of W6.48 **257**, Y7.53 **297**, R3.50 **102**, D2.50 **50**, F6.44 **253**. Monitors =
  those joints, reporting each switch's field load at the accepted `q`.

**(E) CVs (model resid):** primary TM3–TM6 opening **R102 Cα – E239 Cα**; ionic lock 102–239;
NPxxY Y297 χ; CWxP W257 χ; Na⁺ D50; PIF 182-92-253; interhelical Cα-centroid TM3–TM6/TM5–TM6/TM3–TM5;
hinge dihedral φ/ψ at 236–239 (LEMMA L1).

### 5B. 7JJO — turkey β1AR · agonist isoproterenol · active (offset 39)

**(B) TM cores (weld; model resid):**
TM1 2–25 · TM2 40–65 · TM3 73–98 · TM4 121–143 · TM5 166–189 · TM6 251–276 · TM7 285–304.

**(C) Functional anchors (BW → model resid):**

| BW | element | model resid |
|---|---|---|
| D2.50 | Na⁺ pocket | **48** |
| D3.32 | ligand amine anchor | **82** |
| I3.40 | PIF | **90** |
| D3.49/R3.50/Y3.51 | DRY | **99 / 100 / 101** |
| S5.42/S5.43/S5.46 | catechol serines | **172 / 173 / 176** |
| P5.50 | PIF | **180** |
| E6.30 | ionic-lock partner | **246** |
| F6.44 | PIF | **260** |
| C6.47/W6.48/P6.50 | CWxP toggle | **263 / 264 / 266** |
| N6.55 | pocket | **271** |
| N7.49/P7.50/Y7.53 | NPxxY | **300 / 301 / 304** |

**(D) Reporter worlds.**
- **Coarse:** weld TM cores; monitors = TM6 base **246–251**, TM5 base **190–196**, TM3/DRY base
  **101–105**, TM7/NPxxY base **304–308**, TM2 base **30–35** (ref), TM4 base **107–111** (ref).
- **Focused:** free TM6 hinge φ/ψ **246–251**; χ of W6.48 **264**, Y7.53 **304**, R3.50 **100**,
  D2.50 **48**, F6.44 **260**.

**(E) CVs:** primary **R100 Cα – E246 Cα**; ionic lock 100–246; NPxxY Y304; CWxP W264; Na⁺ D48;
PIF 180-90-260; interhelical centroids; hinge φ/ψ at 246–250.

### 5C. 7DH5 — dog β3AR · agonist mirabegron · active (offset 35)

**(B) TM cores (weld; model resid):**
TM1 5–25 · TM2 39–65 · TM3 74–98 · TM4 119–139 · TM5 168–188 · TM6 260–283 · TM7 287–311.

**(C) Functional anchors (BW → model resid):**

| BW | element | model resid |
|---|---|---|
| D2.50 | Na⁺ pocket | **48** |
| D3.32 | ligand amine anchor | **82** |
| I3.40 | PIF | **90** |
| D3.49/R3.50/Y3.51 | DRY | **99 / 100 / 101** |
| S5.42/S5.43/S5.46 | catechol serines | **173 / 174 / 177** |
| P5.50 | PIF | **181** |
| E6.30 | ionic-lock partner | **252** |
| F6.44 | PIF | **266** |
| C6.47/W6.48/P6.50 | CWxP toggle | **269 / 270 / 272** |
| N6.55 | pocket | **277** |
| N7.49/P7.50/Y7.53 | NPxxY | **307 / 308 / 311** |

**(D) Reporter worlds.**
- **Coarse:** weld TM cores; monitors = TM6 base **252–256**, TM5 base **205–211**, TM3/DRY base
  **101–105**, TM7/NPxxY base **311–315**, TM2 base **29–34** (ref), TM4 base **105–110** (ref).
- **Focused:** free TM6 hinge φ/ψ **252–256**; χ of W6.48 **270**, Y7.53 **311**, R3.50 **100**,
  D2.50 **48**, F6.44 **266**.

**(E) CVs:** primary **R100 Cα – E252 Cα**; ionic lock 100–252; NPxxY Y311; CWxP W270; Na⁺ D48;
PIF 181-90-266; interhelical centroids; hinge φ/ψ at 252–256.
NOTE: mirabegron's β3 selectivity involves an ECL2 exosite; those residues are in ECL2 (≈ model
141–149) and are present in this model (loops built), but their identities are from Nagiri et al.
(§7) — select them only if the extracellular pocket, not the activation core, is the target.

---

## 6. Cross-structure comparison (β1/β2/β3 are homologs)

**Directly comparable across the three runs:**
- All three retain the full microswitch set (verified per receptor, §2), so the same coarse-world
  boundary joints and the same focused-world switches have the same mechanistic meaning in each.
- The **primary CV is homologous**: R3.50 Cα – E6.30 Cα, i.e. **102–239 (β2) / 100–246 (β1) /
  100–252 (β3)** in model numbering. This is the one CV to compare head-to-head. All three are
  active-state, so the comparison target is the open-TM6 load field, not a transition.

**Do NOT transplant numbers — per-receptor and per-species:**
- Offsets differ (29 / 39 / 35); **7JJO is turkey**, **7DH5 is dog**. Always use the §5 model
  numbers, verified per file (P1).
- Because loops are built in all three (contiguous), ICL3-interior joints are selectable in every
  model — unlike the deposited PDBs where ICL3 is partly unresolved.

---

## 7. Sources

Model files: `examples/febs/{3SN6,7DH5,7JJO}.{lig,noLig}.nanodisc.{prmtop,rst7}` (ground-truth layout
and numbering verified with ParmEd; §2). BW→native anchors: GPCRdb (https://gpcrdb.org/) and UniProt
**P07550** (β2 human), **P07700** (β1 turkey), **O02662** (β3 dog); model offsets verified by conserved
motif (DRY/CWxP/NPxxY) against each prmtop sequence. Structures: RCSB
[3SN6](https://www.rcsb.org/structure/3SN6) (Rasmussen 2011, https://doi.org/10.1038/nature10361),
[7JJO](https://www.rcsb.org/structure/7JJO) (Su 2020, https://doi.org/10.1016/j.molcel.2020.08.001),
[7DH5](https://www.rcsb.org/structure/7DH5) (Nagiri 2021,
https://www.cell.com/molecular-cell/fulltext/S1097-2765(21)00504-9). Canonical class-A switches: as in
spec 10 §9. Robosample: `spiridon_2020_robosample`, `spiridon_2017_cdhmc_gibbs`;
`docs/specs/reaction-force-monitoring.md`; `docs/specs/gpcr-world-design/10-classA-world-selection-conceptual.md`.

---

## 8. Correctness conditions

- **P1 (anchor identity per model).** The §5 model residues carry the expected identities (already
  verified: e.g. 3SN6 257=Trp/102=Arg/239=Glu; 7JJO 264=Trp/100=Arg; 7DH5 270=Trp/100=Arg). Re-check
  after any re-build of the models. Runtime guard on the run script.
- **P2 (no residue outside the receptor chain).** Selections name only molecule-0 receptor residues
  (1..312 / 1..318 / 1..325); never MSP/lipid/agonist. The models are contiguous (no receptor gaps),
  so no gap-exclusion is needed — but a selection SHALL still stay within the receptor's model range.
- **P3 (reporter integrator).** The reporter world's sampler is not Cartesian-only
  (`reaction-force-monitoring.md` §3 NOTE).
- **I1 (welded-domain monitors).** Coarse-world monitored joints bound welded TM cores; a monitored
  joint whose outboard subtree is itself free carries ≈0 (mis-selection).
- **I2 (environment mobile).** MSP belt and lipids keep their default `Free` roots (load_amber gives
  every molecule a Free root for a boxed/nanodisc system, `context.py:125`); do not weld them around
  the receptor (spec 10 I2 welded-wall failure).
- **I3 (frame-invariant CVs).** §5E CVs are internal distances/dihedrals; no lab-frame RMSD.
- **LEMMA L1 (TM6-base joint tracks opening).** On frames where R3.50Cα–E6.30Cα changes, the TM6-base
  monitored joint's reaction and the TM6 hinge dihedral change together; opening-CV change with no
  hinge-dihedral change is rigid-body drift (Claim-C1 artifact), not a hinge event.
- **NOTE (built loops).** ICL3/ECL2 joints sit on modeled loops; treat their absolute force values as
  lower-confidence than TM-core-boundary joints.

---

## 9. Open questions

- **Q1 (active-state only).** All three files are agonist-bound active complexes (G protein stripped);
  they characterize the active-state load/velocity field, not an inactive→active transition. A
  transition study needs a paired inactive structure (e.g. β2AR 2RH1) — out of scope.
- **Q2 (velocity `V_GB`).** The reporter emits only the static `u=0` reaction force today; per-body
  spatial velocity `V_GB` of the welded coarse bodies is the natural companion the user asked for and
  is a small future engine extension (companion to `reaction-force-monitoring.md`). Its selection is
  the SAME as here. Not designed in this spec.
- **Q3 (`.lig` vs `.noLig`).** Selections are identical for both variants (same receptor chain); the
  agonist is a separate molecule (present only in `.lig`). Comparing the two force fields (with/without
  bound agonist) is a natural use of these selections.
