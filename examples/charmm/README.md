# CHARMM36 CHAMBER test asset

`charmm36_ala5.prmtop` / `charmm36_ala5.rst7` is a small, clean CHAMBER-format
AMBER topology/coordinate pair for a capped pentapeptide, ACE-(ALA)5-NME,
parameterized with the CHARMM36 protein force field. It exists to exercise the
CHAMBER-specific sections of the AMBER prmtop format end to end (against both
native OpenMM and Robosample's own loader) without the confounders present in
`examples/GfcDstrippedMin.prmtop` (NBFIX, a degenerate box, 45595 atoms):

* `%FLAG CTITLE` (CHAMBER format marker)
* `%FLAG CHARMM_CMAP_COUNT` / `CHARMM_CMAP_RESOLUTION` / `CHARMM_CMAP_PARAMETER_nn`
  / `CHARMM_CMAP_INDEX` (backbone CMAP correction maps)
* `%FLAG CHARMM_UREY_BRADLEY*` (1-3 Urey-Bradley terms)
* `%FLAG CHARMM_IMPROPERS` / `CHARMM_IMPROPER_*` (CHARMM harmonic impropers)
* `%FLAG LENNARD_JONES_14_ACOEF/BCOEF` with `SCEE_SCALE_FACTOR` =
  `SCNB_SCALE_FACTOR` = 1.0 (CHARMM's 1-4 nonbonded convention: scaling is
  baked into the separate 1-4 LJ table rather than a global AMBER-style
  1.2/2.0 scalar)

It deliberately has **no NBFIX** (`prmtop_reader.has_nbfix_fast` is `False`)
and **no box** (`IFBOX == 0` in `POINTERS`, and the `.rst7` carries no box
line at all -- genuinely absent, not the degenerate all-zero box seen in
`GfcDstrippedMin.rst7`), so both `openmm.app.AmberPrmtopFile` /
`AmberInpcrdFile` and Robosample's `Context.load_amber` can build a System and
compute a single-point energy without any workarounds. See
`tests/test_openmm_potential_energy.py` (the `charmm36_ala5` case) for the
full per-force-group validation against native OpenMM.

## Regenerating

```
python3 examples/charmm/generate_charmm36_ala5.py
```

Requires (all present in the `robo_cuda13.0` conda env):

* `tleap` (AmberTools, on `PATH`) -- builds a clash-free extended-conformation
  seed geometry with ff14SB (used only for coordinates, not parameters).
* `vmd` with the `psfgen` Tcl plugin (on `PATH`) -- builds the CHARMM-typed
  PSF from `top_all36_prot.rtf`.
* `parmed` (Python) -- loads the PSF + `par_all36_prot.prm` and writes the
  final CHAMBER-format `prmtop`/`rst7`.
* CHARMM36 protein parameters at `$CONDA_PREFIX/dat/chamber/{top,par}_all36_prot.*`.

See the module docstring in `generate_charmm36_ala5.py` for the full 5-step
pipeline (tleap seed -> AMBER->CHARMM atom-name remap -> psfgen PSF -> a small
deterministic coordinate jitter -> ParmEd CHAMBER prmtop/rst7 export) and why
each step is needed -- in particular, why the raw tleap-built geometry (exactly
planar along the backbone) is jittered before export: an exactly-planar
structure makes every CHARMM improper and CMAP term evaluate to exactly zero,
which would make the corresponding energy comparison pass trivially even if
the improper/CMAP implementation were broken.

## Tools that were evaluated and NOT used

* `psfgen` as a **Python** module: not installed in this environment
  (`ModuleNotFoundError`). VMD 1.9's Tcl `psfgen1.4` plugin (`vmd -dispdev
  text -e script.tcl`) was used instead -- same underlying tool, different
  entry point.
* AmberTools' `chamber` executable (the tool the CHAMBER format is named
  after, which converts a CHARMM PSF + parameters directly to a CHAMBER
  prmtop): not present on `PATH` in this environment. ParmEd's own Amber
  writer (`Structure.save(..., format='amber')` on a `CharmmPsfFile`-loaded,
  CHARMM-parameterized `Structure`) produces byte-for-byte the same section
  layout (verified against `ParmEd`'s `ChamberParm._cmap_prefix ==
  "CHARMM_"`) and was used instead.
* `pdbfixer`: importable, but not needed -- this asset is built de novo
  (via `tleap`), not repaired from an existing PDB.
