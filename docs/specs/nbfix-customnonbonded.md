# Spec: NBFIX support via `CustomNonbondedForce` (AMBER/CHAMBER off-diagonal LJ)

Status: ready for implementation (after review)
Owner: coder
Gate: per-force-group energy ≤1e-6 rel (double precision) on `GfcDstrippedMin`, no regression on non-NBFIX systems

## Problem restatement

Make Robosample run AMBER/CHAMBER systems containing NBFIX (off-diagonal
Lennard-Jones pairs violating Lorentz-Berthelot), reproducing OpenMM energies.
`src/OpenMMContext.cpp::createCustomNonbondedForce` is a stub that throws; the
concrete target is `examples/GfcDstrippedMin.prmtop` (CHAMBER, NBFIX, NoCutoff,
GBSA-OBC2).

When `SystemTopology.hasNBfix` is true, the standard `OpenMM::NonbondedForce`
built by `createNonbondedForce` no longer carries van-der-Waals: every particle
is added with `(charge, sigma=1.0, epsilon=0.0)` (`src/OpenMMContext.cpp:441-442`),
so the `NonbondedForce` contributes only electrostatics plus the 1-4/exclusion
exceptions. The full pairwise LJ must instead be evaluated by an
`OpenMM::CustomNonbondedForce` that looks A/B coefficients up in a per-LJ-type-pair
table (`SystemTopology.aCoef`/`bCoef`, dimension `numNBTypes^2`), keyed by each
atom's 0-based LJ type index (`SystemTopology.atomsNonbondedIndex`). `initialize()`
already calls `createCustomNonbondedForce` and adds the result as a slow force
group named `"CustomNonbondedForce"` (`src/OpenMMContext.cpp:218-220`); only the
factory body is missing.

## Claims

1. **C1.** All data needed is already in `SystemTopology` and crosses PyBind11:
   `numNBTypes`, `aCoef`/`bCoef`, `atomsNonbondedIndex` (0-based), the
   exclusion/1-4 pair lists. **No new `SystemTopology` field, no Python-loader
   plumbing.** (`src/PyBind11.cpp:208,262-270,282-285`.)
2. **C2.** `aCoef`/`bCoef` are already in OpenMM's `Discrete2DFunction` form:
   `aCoef[k]=sqrt(A_raw)·sqrt(kcal→kJ)·(Å→nm)^6`, `bCoef[k]=B_raw·(kcal→kJ)·(Å→nm)^6`
   (`prmtop_reader.py:407-411` ≡ `amber_file_parser.py:958-965`). Feed verbatim; no
   re-scaling in C++.
3. **C3.** `NONBONDED_PARM_INDEX` and the A/B tables are symmetric, so Robosample's
   row-major `aCoef[i·NT+j]` and OpenMM's column-major `acoef[i+NT·j]` index-
   identically under `Discrete2DFunction(NT,NT,·)`. Use the Robosample arrays
   directly (derivation below).
4. **C4.** 1-4 LJ and all 1-2/1-3/1-4 exclusions are already handled for the NBFIX
   case by `createNonbondedForce` (`src/OpenMMContext.cpp:449-462`): real 1-4
   sigma/epsilon (CHAMBER `LENNARD_JONES_14_*`) as exceptions + zeroed exclusion
   exceptions. The new `CustomNonbondedForce` must **exclude every one of those
   pairs**; `addStandardExclusions` (`src/OpenMMContext.cpp:547-558`) produces
   exactly that set.
5. **C5.** `GfcDstrippedMin` is `NoCutoff` with no `LENNARD_JONES_CCOEF` (no 12-6-4),
   so the plain `(a/r6)^2-b/r6` expression applies; no `ccoef` table.
6. **C6.** The only genuinely new C++ code is the body of
   `createCustomNonbondedForce`.

## PART 1 — Exact OpenMM construction

Citations: `.../site-packages/openmm/app/`.

**Trigger.** `getNonbondTerms()` raises `NbfixPresent` (`amber_file_parser.py:305-307`)
when any off-diagonal `(i,j)` deviates from Lorentz-Berthelot by >1e-6, or a
one-sided zero (`:341-360`); `readAmberSystem` catches it → `nbfix=True`
(`:944-947`). Robosample reproduces this in `prmtop_reader.has_nbfix_fast` →
`system_topology.has_nb_fix` (`context.py:198-205`).

**(a) LJ → `CustomNonbondedForce`** (`amber_file_parser.py:948-991`):
- `numTypes=NTYPES`; `ene_conv=kcal→kJ`, `length_conv=Å→nm`,
  `afac=sqrt(ene_conv)·length_conv^6`, `bfac=ene_conv·length_conv^6` (`:956-959`).
- `acoef[i+NT·j]=sqrt(parm_acoef[nbidx[NT·i+j]−1])·afac`,
  `bcoef[i+NT·j]=parm_bcoef[nbidx[NT·i+j]−1]·bfac`, skip `nbidx−1<0` (`:960-965`).
- No-1264 energy expr (`:979-981`):
  `"(a/r6)^2-b/r6; r6=r^6;a=acoef(type1, type2);b=bcoef(type1, type2);"`.
- `addTabulatedFunction('acoef', Discrete2DFunction(NT,NT,acoef))`, same for bcoef
  (`:982-985`); `addPerParticleParameter('type')`; per atom
  `addParticle((ATOM_TYPE_INDEX-1,))` — **0-based** (`:989-991`).

**Energy identity.** With `a=sqrt(A)·afac`: `(a/r6)^2=A/r^12` (SI), `b/r6=B/r^6` (SI).
**Trap:** the table holds `sqrt(A)`; the expression squares it. Writing `a/r^12`
is wrong by `sqrt(A)`.

**(b) `NonbondedForce` keeps charges only** (`:948-949`): `addParticle(charge,1.0,0.0)`;
`getCharges` = raw/18.2223 (`:235-241`). Robosample mirrors this at
`src/OpenMMContext.cpp:441-442`.

**(c) Exclusions / 1-4 split (load-bearing).**
- **1-4 exceptions live on `NonbondedForce`, not the custom force**
  (`:1016-1027`): CHAMBER reads `LENNARD_JONES_14_ACOEF/BCOEF` (`:589-591`);
  `epsilon=b^2/(4a)`, `rMin=(2a/b)^(1/6)` (`:610-612`); per-dihedral SCEE/SCNB
  (CHAMBER default 1.0; `:616-624`); `chargeProd/=scee`, `epsilon/=scnb`,
  `sigma=rMin·2^(-1/6)`; `force.addException(i,l,chargeProd,sigma,epsilon)`.
- **1-2/1-3 exclusions** → zeroed exceptions `addException(i,j,0.0,0.1,0.0)`,
  skipping pairs already added as 1-4 (`:1029-1035`).
- **Custom force excludes every standard-force exception** (`:1041-1044`):
  for each exception `(ii,jj,...)`, `cforce.addExclusion(ii,jj)`. So 1-4 LJ comes
  solely from the standard-force exception (CHAMBER 1-4 tables), never from the
  NBFIX off-diagonal table; 1-2/1-3 LJ is fully off. Reproduce exactly.

**(d) Method/cutoff/periodicity/switching** (`:1046-1056`):
periodic → `CutoffPeriodic`, `setCutoffDistance(cutoff)`,
`setUseLongRangeCorrection(true)`; `CutoffNonPeriodic` → cutoff, no LRC; `NoCutoff`
→ nothing. Switching applied to both forces only when `switchDistance>0` and
method≠NoCutoff (`amberprmtopfile.py:325-339`); Robosample sets no switching, so
build the reference with default `switchDistance=0`. gfcd is NoCutoff → no
switching.

**Units.** charge raw/18.2223→e; energy kcal→kJ; length Å→nm; A-table carries
`sqrt(ene_conv)·length_conv^6`, B-table `ene_conv·length_conv^6`. All already in
`aCoef`/`bCoef`.

**C3 layout derivation.** `Discrete2DFunction(NT,NT,v)` evaluates `v[x+NT·y]`. OpenMM
`acoef(type1,type2)=g(A_{type1,type2})`. Robosample `aCoef[k]=g(A_{i,j})` for
`k=i·NT+j`; read through the same function: `aCoef[type1+NT·type2]=g(A_{type2,type1})
=g(A_{type1,type2})` by symmetry. Directly usable.

## PART 2 — Robosample integration

**Existing (`src/OpenMMContext.cpp`).** `createNonbondedForce` (`:400-464`): NBFIX
branch `(charge,1,0)` (`:441-442`); 1-4 from `scaling14*` (`:449-455`); exclusions
`(0,0.1,0)` (`:456-462`). `createGBSAOBCForce` (`:466-500`) independent of LJ.
`createCustomNonbondedForce` (`:502-505`) is the stub. `addStandardExclusions`
(`:547-558`) adds `scaling14I/L` then `exclusionI/J` — exactly the standard force's
exception set. Force-group wiring + name `"CustomNonbondedForce"` (`:218-220`).

**Data in `SystemTopology` (`include/TopologyElements.hpp`).** `hasNBfix` (`:188`),
`numNBTypes` (`:189`), `aCoef`/`bCoef` (`:190-191`, size `numNBTypes^2`),
`atomsNonbondedIndex` (`:65`, 0-based, per-atom BFS order aligned with
`atomsCharge`; `molecule_prototype.py:178`, `topology.py:117`),
`scaling14I/L/ChargeProduct/Sigma/Epsilon` (`:156-161`), `exclusionI/J` (`:167-169`,
disjoint from scaling14). All bound in PyBind11.

**C++ change — body of `createCustomNonbondedForce(sys)`, in order:**
1. Energy string `"(a/r6)^2-b/r6; r6=r^6; a=acoef(type1, type2); b=bcoef(type1, type2);"`.
2. `addTabulatedFunction("acoef", new OpenMM::Discrete2DFunction(sys.numNBTypes, sys.numNBTypes, sys.aCoef))`; same for `"bcoef"`/`sys.bCoef`.
3. `addPerParticleParameter("type")`.
4. per atom `addParticle({ double(sys.atomsNonbondedIndex[i]) })`.
5. method/cutoff mirroring `createNonbondedForce`: periodic → `CutoffPeriodic` +
   `setCutoffDistance(sys.nonbondedCutoff)` + `setUseLongRangeCorrection(true)`;
   `CutoffNonPeriodic` → cutoff, no LRC; `NoCutoff` → nothing (gfcd path).
6. `addStandardExclusions(cforce, sys)` (reuse verbatim).
7. return (ownership passes via `addForce`).

**Guards.** P1: `numNBTypes>0`, `aCoef.size()==bCoef.size()==numNBTypes^2`,
`atomsNonbondedIndex.size()==numAtoms`, each index ∈ `[0,numNBTypes)`. P2: reject
12-6-4 (`LENNARD_JONES_CCOEF`) — no `cCoef` field exists (does not fire for gfcd).

## PART 3 — Validation

**Reference-build obstacle.** `openmm_validation._build_reference_system`
(`openmm_validation.py:95-96`) uses `app.AmberInpcrdFile`, which raises
`ZeroDivisionError` on gfcd's degenerate rst7 box. Build the gfcd reference by:
- `prmtop = app.AmberPrmtopFile("examples/GfcDstrippedMin.prmtop")` (NOT AmberInpcrdFile);
- coords via `amber_loader.read_amber_coordinates("examples/GfcDstrippedMin.rst7")` (ignore box);
- `system = prmtop.createSystem(nonbondedMethod=app.NoCutoff, constraints=None, implicitSolvent=app.OBC2, rigidWater=False, removeCMMotion=False)`;
- set the same positions (nm) on both reference and Robosample contexts.
Reuse per-force-group aggregation by `type(force).__name__`; Robosample side uses
`set_separate_force_groups(True)` + `calc_openmm_potential_energy_by_group()`.

**Discriminating groups for gfcd:** `NonbondedForce` (charges + exceptions) and
`CustomNonbondedForce` (off-diagonal LJ). NBFIX-independent groups
(`GBSAOBCForce`, `CMAPTorsionForce`, `HarmonicBondForce`(+UB), `HarmonicAngleForce`,
`PeriodicTorsionForce`, improper `CustomTorsionForce`) validate immediately.

**Acceptance.** `|E_cpp − E_ref| ≤ 1e-6·max(|E_ref|,1)` kJ/mol per group and total,
**on a double-precision platform** (Reference/CPU or CUDA `Precision=double`). Under
CUDA mixed precision (current default, `src/OpenMMContext.cpp:279`) expect ~1e-4;
state the platform when reporting.

## Correctness conditions

- **P1/P2**: runtime guards above.
- **I1 (sqrt(A))**: isolated off-diagonal pair at fixed r ⇒ `E=A_ij/r^12−B_ij/r^6`.
  Fails if implementer writes `a/r^12`.
- **I2 (no double-count)**: custom-force exclusion set == standard-force exception
  set. Test: a 1-4 pair that is also an off-diagonal type pair; total LJ must match
  reference; dropping exclusions double-counts.
- **I3 (end-to-end)**: gfcd `NonbondedForce` + `CustomNonbondedForce` + total match
  reference within tolerance on a double-precision platform.
- **L1**: gfcd triggers `has_nb_fix==True`; exactly one `CustomNonbondedForce` group,
  finite energy, on both sides.
- **L2**: NBFIX-independent groups match before AND after the change (no perturbation).

## Touch list

- **Change:** `src/OpenMMContext.cpp::createCustomNonbondedForce` body only (may call
  `addStandardExclusions`). Rebuild `cuda-tests`/`cuda-release`.
- **No change:** `TopologyElements.hpp`, `PyBind11.cpp`, `context.py`,
  `prmtop_reader.py`, `initialize()` wiring, `createNonbondedForce`.
- **New (validation, Python):** gfcd reference builder bypassing `AmberInpcrdFile`,
  injecting `read_amber_coordinates` positions.

## Verification plan

1. Pre-impl (L2): gfcd per-force-group with current build — Nonbonded/Custom
   throw/fail (stub) but bonded/CMAP/GBSA groups already match.
2. Unit (I1): two-particle off-diagonal pair at fixed r.
3. Unit (I2): 1-4 pair that is an off-diagonal type pair; double-count check.
4. End-to-end (I3/L1): gfcd reference (NoCutoff, OBC2/ACE, injected coords), ≤1e-6
   double precision.
5. Regression: `ala-dipeptide` smoke + `nox -s tests` — factory is dead code for
   non-NBFIX systems.

## Open questions (non-blocking for gfcd)

- **OQ1**: 12-6-4 (`LENNARD_JONES_CCOEF`) unsupported (no `cCoef` field); `has_1264`
  not surfaced to C++. gfcd unaffected. Decide later whether to reject in loader or
  add a `cCoef` + `-c/r^4` term (matters only if a target uses 12-6-4, e.g. divalent
  ions).
- **OQ2**: `load_lj_coefs` raises on negative `NONBONDED_PARM_INDEX` (10-12 HBOND);
  OpenMM tolerates (cell 0). gfcd has none; OpenMM rejects 10-12 anyway. No action.
