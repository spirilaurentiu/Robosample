# Spec: Fast, parmed-free AMBER/CHAMBER topology loader

Status: ready for implementation
Owner: coder (after review)
Gate: differential vs current loader + per-force-group energy ≤ 1e-6 rel, all example systems

## 1. Problem

`Context.load_amber` (`python/robosample/context.py`) takes ~30 s to load a
>80k-atom system. The cost is **not** in reading the file — the existing
`prmtop_reader.parse_prmtop` is already numpy-vectorized and fast. It is in the
full-system use of parmed and pandas:

1. `pmd.load_file(prmtop, xyz=inpcrd)` (context.py:153) — builds a full Python
   object graph, one object per atom/bond/angle/dihedral, over the whole system.
2. `parm.split()` (context.py:166) — molecule-identity grouping over the full
   system.
3. The per-total-atom loop (context.py:242-327) — per-atom parmed attribute
   access (`parm.atoms[i].xx`), dict inserts, and f-string formatting, all in
   Python, over every atom.
4. `pd.concat` inside the per-bond loop — O(N²). **Already fixed** (one-shot
   `pd.DataFrame` build); this spec covers the remaining rewrite.

Goal: reload in ~1-2 s for 80k atoms; drop the parmed dependency; keep CHAMBER
(CMAP/Urey-Bradley/harmonic impropers) working; reproduce energies exactly.

**Scope directive (binding):** the loader is an *exact replica of OpenMM's AMBER
readers* — `openmm/app/amberprmtopfile.py` +
`openmm/app/internal/amber_file_parser.py` — reproducing their **entire**
functionality, not merely the subset the current parmed path exercises. Every
prmtop feature OpenMM's reader supports (all GB models, NBFIX/CustomNonbonded,
PME/Ewald/CutoffPeriodic/CutoffNonPeriodic/NoCutoff, implicit-solvent variants,
CHAMBER CMAP/Urey-Bradley/harmonic impropers/1-4 LJ tables, virtual sites,
constraints/rigid water flags, polarizable/other sections OpenMM handles) must be
parsed and turned into the corresponding `SystemTopology` data. Where a feature
has **no consumer yet in `OpenMMContext.cpp`**, still implement the full Python
parse/representation but **raise a clear exception** at the point of use (see §4a)
rather than silently dropping or approximating it. Nothing OpenMM supports may be
silently ignored.

## 2. Correctness gate (this is the whole safety story)

The current parmed path already fills `SystemTopology` correctly — every bonded,
nonbonded, and CMAP per-force-group energy matches native OpenMM today (the only
open failure is a **pre-existing, C++-side GBSA-OBC2 discrepancy of ~15 kJ/mol,
out of scope here** — see §8). Therefore the new loader is a *refactor with a
mechanical check*, not a blind reimplementation:

**Primary gate — differential.** For every example system, run the current
(parmed) loader and the new loader and assert that **every `SystemTopology`
array is identical**: exact equality for integer/index/bool arrays, `≤1e-9`
relative for float arrays; and that `df_bonds` is equal. Ship this as a test
(§7). If the arrays are identical, energies are identical by construction.

**Backstop gate — energy.** `tests/test_openmm_potential_energy.py`
(`compare_by_force_group`, ≤1e-6 relative per force class) must still pass for
every system it covers, plus a **new CHAMBER system** added to the suite (§7),
since CHAMBER paths (CMAP, Urey-Bradley, harmonic impropers, `LENNARD_JONES_14_*`)
may be uncovered by the differential test if no CHAMBER example is wired in.

Energy is **invariant to how molecules are grouped into prototypes** and to atom
ordering within the global arrays (it is a sum over bonded terms and pairwise
nonbonded terms indexed by atom). This is the key de-risker; see §5.

## 3. Architecture

Four components. Keep the public surface of `Context.load_amber` and
`standard_dihedral_bonds` / `build_flexibilities` unchanged.

```
inpcrd/rst7 ─┐
prmtop ──────┤ (A) numpy parse ──► raw arrays (prmtop order)
             │                       │
             │        (B) scipy connected-components + fingerprint dedup
             │                       │  ► prototypes + instance groups
             │                       ▼
             │        (C) PrototypeTopology shim  (per UNIQUE prototype, cheap)
             │                       │   quacks like the parmed Structure subset
             │                       ▼
             │        MoleculePrototype / acyclic_graph / z_matrix /
             │        amber_dihedral_classifier  (UNCHANGED — run on the shim)
             ▼                       │
     (D) context.load_amber: vectorized per-instance flatten ─► SystemTopology
```

### (A) Full-system numpy parse — `amber_loader.py` (new) / extend `prmtop_reader`
- Reuse/extend `parse_prmtop` to return every section the pipeline needs (atoms,
  bonds inc/without H, angles, dihedrals, impropers, exclusions, LJ, GB radii/
  screen, CHAMBER `CHARMM_*` + `CMAP_*` sections, POINTERS).
- Read coordinates from inpcrd/rst7 **directly** into an `(N,3)` float64 numpy
  array (Å→nm). Do **not** go through parmed. Handle both formatted rst7/inpcrd
  and the box line. This replaces every `parm.atoms[i].xx/.xy/.xz`.
- All arrays are in **prmtop order**.

### (B) Vectorized dedup — replaces `parm.split()`
- Build the bond graph edge list from `BONDS_INC_HYDROGEN` +
  `BONDS_WITHOUT_HYDROGEN` (remember the **÷3 atom-index encoding**; 1-based
  section values). Include ring-closing bonds — connectivity only.
- Connected components via `scipy.sparse.csgraph.connected_components` (C-fast).
  Each component = one molecule instance. Instances are numbered by first
  appearance in prmtop atom order (matches parmed; see §5).
- Canonical fingerprint per component and group identical ones into prototypes.
  Fingerprint MUST reproduce parmed's identity criterion **or be strictly
  finer** (never merge non-identical): tuple of, in atom order within the
  component, `(residue_name, atom_name, round(charge,6), round(rmin,6),
  round(epsilon,6))`, plus — to be safe and unambiguous — the intra-component
  bond set remapped to component-local indices. (parmed keys single-residue
  molecules on the first five; adding the bond set only ever *splits* more, which
  is safe. See §5.)
- Output: `prototypes: list[component_atom_indices]`, and
  `molecules: list[(instance_idx, prototype_idx)]` sorted by instance_idx,
  identical in meaning to today's `self.molecules`.

### (C) `PrototypeTopology` shim — replaces the parmed `Structure` per prototype
A lightweight, numpy-backed object exposing exactly the parmed attribute surface
consumed downstream, so `MoleculePrototype`, `acyclic_graph`, `z_matrix`, and
`amber_dihedral_classifier` run **unchanged**. Built once per **unique**
prototype (few of them → cheap even in pure Python).

The coder MUST first enumerate the exact surface by grepping the consuming
modules; the known surface (verify and extend) is:
- `molecule.atoms` → list of atom shims, each with:
  `idx, name, mass, charge, sigma, epsilon, solvent_radius, screen, nb_idx,
  atomic_number, element_name, residue (.idx, .name), bond_partners, xx, xy, xz`.
- `molecule.bonds` → each with `.atom1`, `.atom2`, `.type (.k, .req)`.
- `molecule.angles` → `.atom1/2/3`, `.type (.k, .theteq)`.
- `molecule.dihedrals` → `.atom1/2/3/4`, `.improper` (bool), `.ignore_end`,
  `.type` which may be a `DihedralType` or a `DihedralTypeList`
  (multi-term — must iterate); `.type` exposes `.phi_k, .per, .phase, .scee,
  .scnb`.
- `molecule.impropers` (CHAMBER) → `.atom1/2/3/4`, `.type (.psi_k, .psi_eq)`.
- `molecule.urey_bradleys` (CHAMBER) → `.atom1`, `.atom2`, `.type (.k, .req)`.
- residues: `molecule.residues` with `.idx, .name`.

Values come straight from prmtop sections (units: raw prmtop; conversions stay
where they already are in `MoleculePrototype`). `bond_partners` is derivable from
the component-local bond list.

The shim exposes atoms in **prmtop-local** order (as parmed's split fragment
does); `MoleculePrototype` builds its own BFS/compound order on top. **Do not
change** the BFS root selection, compound ordering, z-matrix walk, or dihedral
classification — the differential gate depends on byte-identical output.

`load_nonbonded_exceptions` (`prmtop_reader`) currently takes
`(molecule.parm_data, molecule)`; re-plumb it to read from the shim / raw arrays
(the 1-4 pair + exclusion logic itself is unchanged and already matches OpenMM).

### (D) `context.load_amber` rewrite
- Replace `pmd.load_file` and `parm.split()` with (A) and (B).
- **Vectorize** the per-total-atom loop (context.py:242-327):
  - coordinates: gather with a single fancy-index from the `(N,3)` array using
    the per-instance prmtop index ranges → `atoms_x/y/z` in compound order.
  - `atoms_unique_name`: build with vectorized string ops (numpy `char` /
    list comprehension over arrays), not per-atom f-strings inside the graph
    walk. Format unchanged: `"{resname}{res}_{atomname}_{atom}"`, global 1-based
    numbers (§ field contract).
  - `atoms_prmtop_index` / `prmtop_to_global_index`: build from the per-instance
    `compound_to_local` maps via array math, not a Python dict-per-atom (a dict
    or an int array is fine as long as CMAP/virtual-site lookups still work).
- Bulk-flush `SystemTopology` fields exactly as today (the `_RANGE_SPECS` /
  `_FIELD_SPECS` machinery in `topology.py` is unchanged).
- Remove `import parmed as pmd` from `context.py` and `molecule_prototype.py`
  and `prmtop_reader.py` once the shim is in. Grep for any remaining
  `pmd.`/`parmed` in the runtime path (note: `roborun.py`, `flex_export.py`,
  `autoblock.py` also import parmed for trajectory/export — those are **out of
  scope**; leave them, but the *load* path must be parmed-free).

## 4. Field contract (must be reproduced byte-for-byte)

The complete `SystemTopology` contract is documented at length in the research
appendix (§9). The differential gate enforces it mechanically; the highlights the
coder MUST NOT get wrong:

- **Index spaces:** all flat atom-indexed arrays are in **compound/BFS order**
  within a molecule, concatenated globally, with `atom_offset` applied per
  `_FIELD_SPECS`. `atoms_prmtop_index[global_bfs] = global_prmtop` is the
  permutation. `df_bonds` `atom1_idx/atom2_idx` are **molecule-local compound**
  indices (offset via `atoms_begin[molecule_idx]`).
- **Z-matrix sentinels:** `z_matrix_j/k/l` use `-1` for absent references and
  **must not receive `atom_offset`**. Preserve the existing sentinel-aware
  offset handling (see `topology.py` / the `_FIELD_SPECS` offset logic) — a
  blind `x + atom_off` corrupts sentinels.
- **Coordinates are per-instance** (read at the instance's prmtop indices), never
  from the shared prototype.
- **Units:** lengths nm, energies kJ/mol, angles rad, charge e (÷18.2223).
- Per-molecule arrays (`atoms_begin/end`, `*_begin/end`, `atoms_root_index`,
  `root_mobilities`), scalars (`num_*`, `num_nb_types`, dielectrics), box/PME
  fields, virtual sites, and CMAP arrays: unchanged semantics.

## 4a. OpenMM feature parity + unsupported-feature policy

The reader must mirror OpenMM's `AmberPrmtopFile.createSystem` feature set. For
each feature, one of two outcomes:

1. **Consumed by `OpenMMContext.cpp` today** → parse it and populate the matching
   `SystemTopology` fields (differential + energy gates apply).
2. **Not yet consumed by `OpenMMContext.cpp`** → still parse it fully and build
   the Python-side representation, but raise a precise, greppable exception
   **only when that feature is actually present in the input** (so unaffected
   systems are unimpaired). Use a dedicated exception type, e.g.
   `UnsupportedTopologyFeature`, with a message naming the feature and the
   OpenMM behavior being deferred. Do **not** raise at import or unconditionally.

Before writing code, produce a **feature-coverage table** (as a module docstring
or short `docs/specs/loader-feature-matrix.md`) enumerating every feature in
OpenMM's Amber reader and marking each: *consumed* / *parsed-but-raises* /
*n-a*. This table is part of the deliverable and the reviewer will check it
against the OpenMM source for completeness — "exact replica" is graded against
it. Candidates likely to be *parsed-but-raises* today: NBFIX → CustomNonbonded,
non-default GB models (OBC1/GBn/GBn2), polarizable terms, certain constraint/
rigid-water requests, PME parameter overrides not currently wired — verify each
against `OpenMMContext.cpp`.

## 5. Dedup correctness (the one subtle risk)

Energy does not depend on prototype grouping — only on each atom getting its
correct parameters and coordinates, and each bonded/nonbonded term referencing
the correct atoms. Consequences:
- **Safe:** over-splitting (treating two truly-identical molecules as two
  prototypes). Costs a little speed, never correctness.
- **Unsafe:** merging two molecules that differ in topology or any parameter →
  one instance would inherit the wrong prototype topology → wrong energy.
- Therefore the fingerprint must be **conservative**: identical fingerprint ⇒
  provably identical molecules. Include the intra-component bond set (local
  indices) in the fingerprint in addition to parmed's per-atom key. The
  differential gate will catch any accidental over-merge (arrays would differ)
  and any instance-ordering drift.
- Instance numbering/order: components numbered by first-appearing atom in prmtop
  order; `self.molecules` sorted by instance_idx — reproduces today's ordering.

## 6. Implementation order (land speed first, keep the gate green each step)

1. **(done)** one-shot `df_bonds`.
2. **Coordinates + parse off parmed:** add direct rst7/inpcrd reader and the
   extended numpy parse; feed coordinates and per-atom scalar arrays from numpy
   while still using parmed for `split()`/prototype. Vectorize the context.py
   atom loop. → differential gate green. Big share of the 30s gone here.
3. **Dedup off parmed:** replace `parm.split()` with scipy CC + fingerprint. →
   differential gate green.
4. **Shim + remove parmed:** introduce `PrototypeTopology`, point
   `MoleculePrototype`/`acyclic_graph`/`z_matrix`/classifier/
   `load_nonbonded_exceptions` at it, delete parmed from the load path. →
   differential + energy + CHAMBER gates green.

Each step is independently shippable and independently verified by the
differential gate.

## 7. Tests to add (part of the deliverable)

- `tests/test_loader_differential.py`: for each example system (incl. explicit-
  solvent/PME and implicit/GBSA), load with the pre-rewrite path (git-pinned
  reference arrays, or a parmed-based reference builder kept in the test) and the
  new path; assert every `SystemTopology` array + `df_bonds` match (exact ints,
  1e-9 floats). This is the authoritative gate for the refactor.
- A **CHAMBER** example (CMAP + Urey-Bradley + harmonic impropers +
  `LENNARD_JONES_14_*`) wired into `test_openmm_potential_energy.py` so the
  CHAMBER numeric path is covered by the ≤1e-6 per-force-group check.
- A load-time smoke assertion on the largest available system (informational; not
  a hard perf gate, but log the wall-clock so regressions are visible).

## 8. Out of scope (flag, do not fix here)

- **GBSA-OBC2 ~15 kJ/mol discrepancy** in the C++ `GBSAOBCForce`
  (`test_openmm_potential_energy.py[implicit]` is red on `disasm` today, on the
  clean tree, independent of this work). Explicitly deferred by the user — to be
  addressed when the OpenMM interface is rewritten. The loader must keep feeding
  the same `atoms_radius/screen` as today (the differential gate ensures this).
- `roborun.py` is being removed (unmaintained). The supported driver is
  `python/robosample/run.py`; use it for smoke/manual checks, not `roborun.py`.

## 9. Research appendix (reference; sources to verify against)

Two authoritative extractions back this spec:
- **OpenMM numeric contract** (from
  `openmm/app/internal/amber_file_parser.py`, `customgbforces.py`,
  `amberprmtopfile.py`): exclusions from `NUMBER_EXCLUDED_ATOMS` /
  `EXCLUDED_ATOMS_LIST`; 1-4 pairs from DIHEDRALS with both 3rd & 4th indices
  positive, atoms via ÷3; SCEE/SCNB per-dihedral (defaults AMBER 1.2/2.0,
  CHAMBER 1.0/1.0), applied as divisors; LJ via `NONBONDED_PARM_INDEX` with
  `rMin=(2A/B)^(1/6)`, `ε=B²/(4A)`, `σ=rVdw·2^(-1/6)`; NBFIX ⇒ CustomNonbonded;
  CHAMBER `LENNARD_JONES_14_*`, harmonic impropers via CustomTorsion, CMAP grid
  circular phase-shift by `ngrid//2` with i-outer/j-inner transpose; charge
  factor 18.2223, 4.184, 0.1.
- **parmed `split()` identity + `SystemTopology` field contract** (from
  `parmed/structure.py:1359-1428`, `parmed/utils`, and the Robosample
  `topology.py` / `molecule_prototype.py` / `z_matrix.py`): connectivity-only
  components; single-residue fast-path key `(resname, len, atom_names, charges,
  rmins, epsilons)` @6dp; full field-by-field contract with index spaces, units,
  and sentinel rules.

Both are captured in the conversation that produced this spec; the coder should
re-verify any line-number-specific claim against the installed sources, since the
differential gate — not the appendix — is the binding correctness criterion.
