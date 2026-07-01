# FlexWizard — step-by-step tutorial

FlexWizard is a PyMOL plugin for **authoring Robosample `.flex` files**: you
decide which bonds are joints (mobilizers) and of what type, and it writes the
`.flex` Robosample reads. This walkthrough uses `2ala.2ala` (two identical
alanine dipeptides) so you can see every feature, including the
prototype-copy tools.

A bond is stored by its **0-based prmtop atom indices** (`PyMOL atom id − 1`),
the same indices Robosample uses — so the wizard's picture and the engine agree.

> **Prerequisite for dihedral detection.** FlexWizard detects AMBER dihedrals
> *on the fly* by importing Robosample inside PyMOL. Run PyMOL **from the
> Robosample conda env** (e.g. `conda install -c conda-forge pymol-open-source`
> into the same env, then launch `pymol`). If Robosample is not importable, the
> wizard still works in fully manual mode — you just lose the named selections,
> ring-closing locks, and copy propagation.

---

## 1. Launch PyMOL and load the structure

```
# in the PyMOL command line
load examples/2ala.2ala.prmtop, mol
load examples/2ala.2ala.rst7, mol         # coordinates onto the same object
```

Make bonds visible — FlexWizard communicates through **bond color**, so you must
be in a representation that draws bonds:

```
hide everything
show sticks
set valence, 0            # optional: flat bonds, easier to read colors
bg_color white            # optional
orient
```

Enable exactly **one** object (FlexWizard maps the single enabled object).

---

## 2. Start the wizard (with the topology → auto-detect)

```
run tools/FlexWizard.py
start_flex_wizard(prmtop="examples/2ala.2ala.prmtop")
```

Passing `prmtop` makes the wizard build the bond map **and immediately detect
dihedrals** in-process (no file is written). The coordinate file is found
automatically next to the prmtop; pass it explicitly if needed:
`start_flex_wizard(prmtop="...prmtop", inpcrd="...rst7")`.

If you started with a bare `start_flex_wizard()`, click **Detect Dihedrals
(AMBER)** at any time (it will ask for the prmtop), or set it first:
`cmd.get_wizard().set_topology("examples/2ala.2ala.prmtop")`.

A control panel appears on the right. On startup every bond is `Rigid`. Buttons:

| Button | What it does |
|---|---|
| **Joint type: …** | Pick the active joint type (Rigid / Torsion / BallF / BallM / Cartesian). |
| **Manual Selection On/Off** | Toggle click-two-atoms picking. |
| **Rebuild Bond Map** | Recompute the graph (after changing the enabled object). |
| **Save Flexibility File** | Write the `.flex`. |
| **Load Flexibility File** | Read a `.flex` back in to edit. |
| **Detect Dihedrals (AMBER)** | Generate bond metadata on the fly (no file). |
| **Apply to all copies On/Off** | Mirror each assignment onto every prototype copy (step 4). |
| **phi / psi / Backbone / Sidechain (chi) → Torsion** | One-click named selections (step 3). |
| **Set Selection** | Set every bond inside the current PyMOL selection to the active type. |
| **Open Joint Control** | Live torsion-rocking preview (step 6). |
| **Help Me!** | Print a short reminder to the log. |
| **Done** | Close the wizard and clear coloring. |

After detection:
- The wizard checks the detected atom count against the object and **warns** if
  they differ (that means indices won't line up — the enabled object must be the
  same topology).
- Bonds are **pre-colored**: chemically rigid / ring-closing bonds go grey,
  rotatable-but-unassigned candidates go **white**.
- The log reports how many multi-copy groups were found (31 for `2ala.2ala`).

You now have a color-coded map of *where it is even worth putting a torsion*.

---

## 3. Assign torsions the fast way (named selections)

With dihedrals detected, use the named buttons:

- **Backbone → Torsion** sets every φ and ψ bond to Torsion.
- **phi → Torsion** / **psi → Torsion** do one class at a time.
- **Sidechain (chi) → Torsion** sets χ₁…χ₅.

For `2ala.2ala`, **Backbone → Torsion** turns 8 bonds green (4 per copy). The
log prints the count. Ring-closing bonds are silently refused (kept rigid),
because Robosample forces them rigid anyway.

---

## 4. Prototype copies — do it once, apply everywhere

`2ala.2ala` is two identical molecules. You rarely want to hand-edit each copy
(imagine 602 waters). Turn **Apply to all copies → On**, then assign a bond on
*one* copy — the wizard mirrors it to the equivalent bond on every other copy of
that prototype.

Example: with the toggle on, set the φ bond `16-18` (copy 0) to Torsion → its
partner `48-50` (copy 1) becomes Torsion too. (The named selections in step 3
already cover all copies, because detection lists every copy — the toggle is for
**manual picks and Set Selection**.)

### Manual picking (any bond)

1. Set **Joint type** to the type you want.
2. **Manual Selection → On**.
3. Click the two atoms of a bond. The bond recolors to the type's color.
4. Repeat. **Manual Selection → Off** when done.

### Set Selection (a region at once)

1. Select atoms in PyMOL (e.g. `select loop, resi 5-9`).
2. Set the **Joint type**.
3. Click **Set Selection** — every bond fully inside the selection gets that
   type.

---

## 5. Save

Click **Save Flexibility File** and name it (e.g. `2ala.flex`). The file lists
`atom1  atom2  jointType  # comment` for every bond, sorted by first atom. You
can reopen it later with **Load Flexibility File**.

---

## 6. How do we visualize?

FlexWizard has **two** visual channels.

### A. Static — the color-coded kinematic tree

Every bond is tinted by its assignment, so the picture *is* the rigid-body
decomposition. You must be showing bonds (`show sticks` or `show lines`) for the
colors to appear. Legend:

| Color | Meaning |
|---|---|
| **Dim grey** `0x696969` | Rigid (welded — part of a rigid body) |
| **White** `0xFFFFFF` | Rotatable *candidate*, not yet assigned (from detection) |
| **Dark green** `0x006400` | Torsion (1-DOF dihedral) |
| **Deep-sky blue** `0x00BFFF` | Ball (mobile body) |
| **Spring green** `0x00FFBF` | Ball (fixed body) |
| **Light yellow** `0xFFFF66` | Cartesian (3 translations) |

Reading it: contiguous runs of grey bonds are single rigid bodies; the colored
bonds between them are the joints. This is the fastest way to sanity-check that
your rigid bodies are what you intended before spending GPU time.

Tips:
```
show sticks
set stick_radius, 0.15
set valence, 0
```

### B. Dynamic — the torsion-rocking preview

To *see the motion a torsion unlocks*, click **Open Joint Control**. Every bond
you set to Torsion becomes a row with a slider (set the angle), a **Speed**, and
an **Amp**:

1. Set a **Speed** (°/step) and **Amp** (°) — or use "Set all speeds/amplitudes".
2. Click **Start All**. The selected torsions oscillate live in the viewer,
   sinusoidally around their current value.
3. **Commit Changes** keeps the current pose as the new rest angle; **Close**
   restores the original angles.

Use this to pre-screen joints: a torsion that swings widely without clashing is
a good large-amplitude (basin-hopping) degree of freedom; one that immediately
clashes is not worth the sampling budget.

**One caveat on the preview.** PyMOL rotates the atoms on *one* side of the bond;
Robosample instead holds the side toward the molecule's root atom fixed and
moves the outboard subtree. A torsion is a *relative* coordinate, so the sampled
conformations (and their energies) are identical either way — but the *absolute*
motion you see in the preview may hold the opposite side still from what a
Robosample trajectory frame would show. Judge shapes and clashes, not absolute
placement.

---

## Recap

```
# PyMOL launched from the Robosample conda env
load examples/2ala.2ala.prmtop, mol
load examples/2ala.2ala.rst7, mol
show sticks
run tools/FlexWizard.py
start_flex_wizard(prmtop="examples/2ala.2ala.prmtop")   # auto-detects dihedrals
#   Apply to all copies On   (for manual picks)
#   Backbone -> Torsion
#   (optional) Open Joint Control -> Start All
#   Save Flexibility File    -> 2ala.flex
```

### Headless / remote PyMOL (Robosample not importable in the viewer)

If PyMOL cannot import Robosample, precompute a sidecar in the Robosample env and
skip on-the-fly detection:

```python
from robosample import AmberDihedralClassifier, flex_export
ctx = robosample.Context("2ala", 42, AmberDihedralClassifier())
ctx.load_amber("examples/2ala.2ala.prmtop", "examples/2ala.2ala.rst7")
flex_export.export_flex_meta(ctx, "examples/2ala.2ala.prmtop", "2ala.flexmeta")
```

The in-viewer detection and this file writer share the same code
(`flex_export.compute_flex_meta`), so they produce identical metadata.
```
