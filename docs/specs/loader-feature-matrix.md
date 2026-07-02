# Loader feature-coverage matrix

Companion to `docs/specs/fast-amber-loader.md` §4a. Enumerates every feature
in OpenMM's AMBER reader (`openmm/app/amberprmtopfile.py::AmberPrmtopFile` +
`openmm/app/internal/amber_file_parser.py::PrmtopLoader`/`readAmberSystem`,
installed at
`$CONDA_PREFIX/lib/python3.12/site-packages/openmm/app/{amberprmtopfile.py,internal/amber_file_parser.py}`,
OpenMM 8.5.1) and marks whether Robosample's loader (`python/robosample/
context.py`, `prmtop_reader.py`, `molecule_prototype.py`) + `OpenMMContext.cpp`
consume it, parse-but-raise, or don't apply. Verified against
`src/OpenMMContext.cpp` (function names cited below).

Legend: **consumed** = parsed and fed into `SystemTopology` and a live
`OpenMM::Force` in `OpenMMContext.cpp`. **parsed-but-raises** = the Python
side must still parse/represent it fully, but `OpenMMContext.cpp` raises a
clear exception when the feature is actually present (§4a policy).
**n/a** = not applicable to Robosample's architecture (e.g. constraints,
which Robosample replaces with rigid-body mobilizers) or not part of the
prmtop/inpcrd numeric contract (removeCMMotion, platform selection).
**GAP** = a real, pre-existing discrepancy from the policy, flagged for
follow-up, not fixed in this pass (Rule 3: surgical, don't fix adjacent bugs
unless asked).

## Topology / bonded

| Feature (OpenMM) | Robosample status | Notes |
|---|---|---|
| Atoms, masses, charges, `ATOM_NAME`, `RESIDUE_LABEL`/`RESIDUE_POINTER` | consumed | `atoms_mass/charge/...`, `atoms_unique_name` |
| `ATOMIC_NUMBER` (element) / name-based element guess fallback | consumed | `atoms_atomic_number` always from `ATOMIC_NUMBER`; Robosample does not implement OpenMM's name-based guess fallback for prmtops lacking `ATOMIC_NUMBER` (rare; all AMBER `tleap` prmtops carry it) | n/a in practice — GAP if ever hit (no example exercises it) |
| Bonds (`BONDS_INC_HYDROGEN` / `BONDS_WITHOUT_HYDROGEN`) | consumed | `bonds_i/j`, `bonds_stiffness/equilibrium`; special-cased water H-H "bond" is skipped by OpenMM's `Topology`-building step only (cosmetic, doesn't affect the `HarmonicBondForce`) — n/a for Robosample (no `Topology` object) |
| Angles | consumed | `angles_i/j/k`, stiffness/equilibrium |
| Proper dihedrals (multi-term `DihedralTypeList`) | consumed | `periodic_torsions_*`, one record per term |
| Impropers (plain AMBER periodic-torsion-style, 4th pointer negative) | consumed | folded into `periodic_torsions_improper` flag, same `PeriodicTorsionForce` path as OpenMM |
| CHAMBER harmonic impropers (`CHARMM_IMPROPERS`) | consumed | `harmonic_torsions_*` -> `createImproperHarmonicTorsionForce` (`CustomTorsionForce`, matches OpenMM's `k*min(dtheta,2pi-dtheta)^2` form) |
| CHAMBER Urey-Bradley (`CHARMM_UREY_BRADLEY`) | consumed | `urey_bradley_*` -> `createUreyBradleyForce` (`HarmonicBondForce`, matches OpenMM's separate `UreyBradleyForce`) |
| CHAMBER CMAP (`CMAP_RESOLUTION`/`CMAP_PARAMETER_*`/`CMAP_INDEX`) | consumed | `cmap_*` -> `createCMAPTorsionForce`; same `ngrid//2` phase shift + i/j transpose as OpenMM |
| CHAMBER `LENNARD_JONES_14_A/BCOEF` (separate 1-4 LJ table) | consumed | `prmtop_reader.load_nonbonded_exceptions` prefers `LENNARD_JONES_14_*` over the plain table when present |
| SCEE/SCNB per-dihedral scaling (defaults AMBER 1.2/2.0, CHAMBER 1.0/1.0) | consumed | `scaling14_charge_product`/`epsilon` divide by `SCEE_SCALE_FACTOR`/`SCNB_SCALE_FACTOR` from the prmtop; Robosample has no `scee=`/`scnb=` *override* (OpenMM's `createSystem(scee=..., scnb=...)`) — n/a, `load_amber` doesn't expose that knob |
| 1-4 exceptions from dihedral pointers (3rd/4th index sign convention) | consumed | `prmtop_reader.load_nonbonded_exceptions` |
| Explicit exclusions (`NUMBER_EXCLUDED_ATOMS`/`EXCLUDED_ATOMS_LIST`) | consumed | same function, `exclusion_i/j` |
| 10-12 H-bond terms (`HBOND_ACOEF`/`HBOND_BCOEF`, legacy pre-1994 force fields) | n/a | OpenMM itself raises `Exception('10-12 interactions are not supported')` if any nonzero coefficient is present; Robosample never reads `HBOND_ACOEF`/`HBOND_BCOEF` at all -- GAP: a legacy 10-12 prmtop would be silently accepted by Robosample instead of raising. No example system carries nonzero 10-12 terms. |
| `IFPERT` (perturbation) / `IFCAP` (CAP) | n/a | OpenMM raises `Exception` unconditionally if set; Robosample never reads these `POINTERS`. GAP, same rationale as 10-12 (untested/unused legacy AMBER features) |

## Nonbonded

| Feature (OpenMM) | Robosample status | Notes |
|---|---|---|
| Standard combining-rule LJ (`sigma`/`epsilon` per atom) | consumed | `atoms_sigma/epsilon` via `load_lj_coefs` diagonal terms |
| NBFIX / off-diagonal LJ -> `CustomNonbondedForce` | **parsed-but-raises (fixed in Step 4b)** | `prmtop_reader.has_nbfix_fast` implements the exact same combining-rule-deviation check as OpenMM's `getNonbondTerms`; `context.py::load_amber` now calls it (on the raw, unconverted `NONBONDED_PARM_INDEX`/`LENNARD_JONES_A/BCOEF` tables) and sets `system_topology.has_nb_fix`, which makes `OpenMMContext.cpp::createCustomNonbondedForce`'s existing unconditional `"not implemented yet"` guard (gated on `systemTopology.hasNBfix`, `initialize()` line ~218) actually fire for NBFIX prmtops instead of silently falling back to combining-rule LJ. Verified on `examples/GfcDstrippedMin.prmtop` (CHAMBER, genuinely NBFIX): loads with `has_nb_fix == True`, then `Context.initialize()` raises `RuntimeError: createCustomNonbondedForce (NBFIX) not implemented yet`. The differential gate's 16 non-NBFIX cases all still set `has_nb_fix == False` (no false positives). This was the pre-existing GAP flagged below Step 2; now closed. |
| `LENNARD_JONES_CCOEF` (1-6-4 potential, ion parameters) -> `CustomNonbondedForce` with `c/r^4` term | n/a (not parsed) | No example carries this section; Robosample doesn't read `LENNARD_JONES_CCOEF` and `OpenMMContext.cpp` has no matching custom force. Same class of gap as NBFIX -- not exercised by any current example, not fixed here. |
| `NoCutoff` | consumed | `NonbondedMethod::NoCutoff`, default for implicit/vacuum |
| `CutoffNonPeriodic` | consumed | `createNonbondedForce`, `createGBSAOBCForce` |
| `CutoffPeriodic` | consumed | |
| `Ewald` | consumed | |
| `PME` | consumed | default for explicit solvent |
| `LJPME` | parsed-but-raises | not a `NonbondedMethod` enumerator in `robo_bindings`; passing it raises at the Python/pybind boundary (`ValueError`/`TypeError`), which satisfies §4a ("raise ... only when actually requested") even though there's no dedicated message. No example uses it. |
| `ewaldErrorTolerance` | consumed | `system_topology.ewald_error_tolerance` |
| `switchDistance` (LJ switching function) | n/a (not parsed) | Robosample's `NonbondedForce` never calls `setUseSwitchingFunction`; no `load_amber` knob. Not exercised by any example (all use full LJ cutoff or NoCutoff). |
| Dispersion correction (`setUseDispersionCorrection`) | consumed | `createNonbondedForce`, method-dependent (mirrors OpenMM defaults per method) |
| `hydrogenMass` (HMR) | n/a (not parsed) | No `load_amber` knob; Robosample's timestep is chosen per mobilizer, not via HMR |

## Implicit solvent (GB)

| Feature (OpenMM) | Robosample status | Notes |
|---|---|---|
| `RADII`/`SCREEN` per-atom GB parameters | consumed | `atoms_radius`/`atoms_screen` (unconditionally populated regardless of `use_gbsa_obc2`, per the field contract) |
| `OBC2` (kappa==0, the built-in `GBSAOBCForce` path) | consumed | `createGBSAOBCForce` -- the ONLY GB model wired into `OpenMMContext.cpp`. This is the path with the known ~15 kJ/mol discrepancy (out of scope, §8 of the loader spec / `test_openmm_potential_energy.py[implicit]`) |
| `HCT` | parsed-but-raises | `use_gbsa_obc2` is a bool (on/off), not a model selector; no code path builds `GBSAHCTForce`. Selecting HCT isn't expressible via `load_amber`'s current API at all -- effectively n/a until a model-selection knob is added. Not attempted this pass (§4a: implement the parse+raise only where the input format could plausibly select it; `use_gbsa_obc2` cannot). |
| `OBC1` | parsed-but-raises | same as HCT |
| `GBn` | parsed-but-raises | same as HCT |
| `GBn2` (needs `GBn2` per-atom extra columns: `S_alpha, S_beta, S_gamma`) | parsed-but-raises | same as HCT; additionally Robosample never reads the extra GBn2 columns |
| Implicit solvent salt / Debye kappa (`implicitSolventSaltConc`/`implicitSolventKappa`) | n/a (not parsed) | No `load_amber` knob; GBSA is always run at kappa=0 (no ionic screening), matching the one GB model that IS wired (`OBC2` kappa==0 branch) |
| `soluteDielectric`/`solventDielectric` | consumed | `system_topology.gbsa_solute_dielectric`/`gbsa_solvent_dielectric`, hardcoded to AMBER defaults (1.0/78.5) in `load_amber`; no override knob (n/a: not exposed, not a parsing gap) |
| `sasaMethod` ('ACE'/'LCPO'/None) | n/a (not parsed) | `createGBSAOBCForce` always uses the built-in `GBSAOBCForce`'s default (ACE-equivalent) SASA term; no LCPO path. Matches the one GB model wired in. |

## Constraints / rigid water / virtual sites

| Feature (OpenMM) | Robosample status | Notes |
|---|---|---|
| `constraints=HBonds/AllBonds/HAngles` (SHAKE) | n/a | Robosample's articulated-body solver represents rigid units via mobilizer topology (`decomposeRigidUnits`), not OpenMM `System::addConstraint`. This is an architectural substitution, not a missing feature -- see `docs/specs/fast-amber-loader.md` §3, `bonds_ring_closing` |
| `rigidWater` | n/a | same substitution; water rigidity comes from the mobilizer/rigid-body decomposition, not OpenMM constraints |
| 3-particle-average virtual sites (4-point water: OPC/TIP4P/TIP4P-Ew) | consumed | `vs_site/vs_atom1-3/vs_weight1-3` -> `ThreeParticleAverageSite`; as of Step 4b, built by `amber_loader.extract_virtual_sites` directly from raw `BONDS_*_HYDROGEN`/`ANGLES_*_HYDROGEN` arrays + instance coordinates (replaces ParmEd's `ExtraPoint.frame_type`/`ThreeParticleExtraPointFrame`, reproduced bit-for-bit -- verified against the `tip4pew` differential golden). An atom is a virtual site iff `ATOMIC_NUMBER == 0` (ParmEd's own criterion). |
| Out-of-plane virtual sites (5-point water: TIP5P) / 2-particle frames / any parent bond count other than 3 | parsed-but-raises | `amber_loader.extract_virtual_sites` raises `amber_loader.UnsupportedTopologyFeature` naming the parent's bond count when it is not exactly 3 (the only frame shape it implements), and again if a 3-bond frame's geometry can't be resolved (no `a1-parent-a2` angle and no direct `a1-a2` bond) -- satisfies §4a. Superseded the pre-Step-4b `context.py` `NotImplementedError` (same policy, now parmed-free). |
| `flexibleConstraints` | n/a | only meaningful together with SHAKE constraints (n/a above) |

## Box / periodicity

| Feature (OpenMM) | Robosample status | Notes |
|---|---|---|
| `IFBOX` / box-vectors-from-restart auto-detection | consumed | As of Step 4b: `has_box = POINTERS[27] (IFBOX) > 0`, matching OpenMM's own `PrmtopLoader.getIfBox()`-driven periodicity decision exactly -- **not** "a box line is present in the coordinate file" (the pre-Step-4b rule, `parm.box_vectors is not None`, which was actually a latent bug: ParmEd's `parm.box = f.box` setter re-derives `IFBOX` from the mere presence of an rst7 box line and overwrites the prmtop's own `IFBOX`, so a non-periodic (`IFBOX==0`) system with a degenerate placeholder box line in its rst7 -- e.g. `examples/GfcDstrippedMin.rst7`'s `0 0 0 90 90 90` -- was silently treated as periodic). Verified: `GfcDstrippedMin` now loads non-periodic. |
| Reduced (lower-triangular) box vectors, incl. truncated-octahedral (`IFBOX==2`) | consumed | As of Step 4b: `amber_loader.box_vectors_from_lengths_angles` reduces the rst7's `(a,b,c,alpha,beta,gamma)` exactly as ParmEd's `geometry.box_lengths_and_angles_to_vectors` does (verified bit-for-bit for orthorhombic differential cases, and against a synthetic truncated-octahedral case -- no octahedral example exists in the repo). If the rst7/inpcrd carries no box line but `IFBOX > 0` (e.g. `examples/b1-1n.rst7`), falls back to the prmtop's own legacy `BOX_DIMENSIONS` section (`[beta, a, b, c]`, all three angles == beta), matching both ParmEd's `readparm.LoadParm` and OpenMM's `PrmtopLoader.getBoxBetaAndDimensions` fallback; `OpenMMContext.cpp::initialize` consumes the 9 floats directly, unchanged. |
| `periodicBoxVectors=`/`unitCellDimensions=` explicit override (`AmberPrmtopFile.__init__`) | n/a | `load_amber` always derives the box from the topology/coordinate files; no override knob (deliberate simplification, not a parsing gap) |

## Top-level `load_amber` policy (not in OpenMM's reader, Robosample-specific)

| Feature | Status | Notes |
|---|---|---|
| `use_gbsa_obc2` auto-detect from box presence | consumed | see `load_amber` docstring |
| CHAMBER-vs-plain-AMBER auto-detect (`CTITLE` flag) | consumed | `prmtop_reader.parse_prmtop["chamber"]`; branches only affect which optional sections exist, not a separate code path |

## Step 2 (that pass) scope note

Step 2 touched only the per-instance coordinate/name/index flattening
(`context.py` lines ~278-304 pre-rewrite) and added a parmed-free coordinate
reader (`amber_loader.py`). It did not change which features are parsed or
which raise -- this matrix reflected the state of the loader both before and
after Step 2 (no entries moved between consumed/parsed-but-raises/n-a/GAP as
a result of that pass).

## Step 4b (this pass) scope note

Step 4b removed the LAST parmed usages from the load path: `pmd.load_file`
in `context.py`, and the two features that were still reading a live ParmEd
`Structure` (box-vector construction and virtual-site/`ExtraPoint` frame
detection). Three entries above moved as a direct result:

* **NBFIX GAP closed**: `context.py` now calls `prmtop_reader.has_nbfix_fast`
  and sets `system_topology.has_nb_fix`, so the pre-existing
  `OpenMMContext.cpp::createCustomNonbondedForce` guard actually fires for
  NBFIX prmtops (previously stuck at its C++ default `false` for every
  system -- see the Nonbonded table above).
* **Box/periodicity** moved from "box line present in the coordinate file"
  (a ParmEd quirk that could silently misclassify a non-periodic system as
  periodic) to "`IFBOX` (POINTERS) > 0", matching OpenMM's own reader -- see
  the Box / periodicity table above.
* **Virtual sites**: frame extraction moved from ParmEd's
  `ExtraPoint.frame_type` to `amber_loader.extract_virtual_sites` (raw
  bond/angle arrays); the parsed-but-raises policy for non-3-particle frames
  is unchanged in substance, only in implementation.

The remaining GAPs (`LENNARD_JONES_CCOEF`, 10-12 H-bond, `IFPERT`/`IFCAP`)
are untouched -- still flagged for a separate, dedicated pass (Rule 3).
