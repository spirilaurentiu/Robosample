"""Validate the C++ OpenMM potential energy against a native OpenMM reference.

These tests are the pytest form of the ``--validate`` path in
``python/robosample/run.py``: for each solvent model they build a
``robosample.Context``, let ``load_amber`` auto-detect the model from the
periodic box (OpenMM-style), and then diff the C++ single-point energy against a
**native OpenMM** reference built from the same ``prmtop``/``rst7`` with
``openmm.app.AmberPrmtopFile`` + ``AmberInpcrdFile`` (NOT ParmEd). For implicit
solvent this means the reference is the built-in ``GBSAOBCForce`` produced by
``createSystem(implicitSolvent=OBC2)`` with the default ``sasaMethod='ACE'`` --
the exact force the C++ engine replicates.

``openmm_validation.compare_by_force_group`` is authoritative: it compares the
total PE (and every force group) and *raises* ``ValueError`` on any mismatch, so
a passing call is the assertion. The reference runs on the OpenMM ``Reference``
platform for a deterministic, full-precision comparison -- exactly as
``run.py`` invokes it.

Four models are covered:

* **implicit** -- ``2ala.implicit`` has no box, so ``load_amber`` selects
  GBSA-OBC2 (implicit solvent, non-periodic).
* **explicit** -- ``2ala.tip3p`` carries a periodic box, so ``load_amber``
  selects PME under PBC and turns GBSA off.
* **vacuum**   -- the same box-less ``2ala.implicit`` files loaded with
  ``use_gbsa_obc2=False``: no implicit solvent, no periodicity (gas phase).
* **charmm36_ala5** -- a clean CHARMM36 CHAMBER-format prmtop/rst7
  (``examples/charmm/``, see ``generate_charmm36_ala5.py``), gas phase. Unlike
  the three AMBER-native cases above, this one exercises the CHAMBER-specific
  prmtop sections: CMAP correction maps, Urey-Bradley terms, and CHARMM
  (harmonic) impropers.

On a bare dev run the tests skip cleanly when the compiled extension, OpenMM, or
the example inputs are unavailable, rather than failing collection. Under the
AUTHORITATIVE gate (``nox -s tests``, which exports ``ROBOSAMPLE_REQUIRE_OPENMM=1``)
these same conditions are HARD FAILURES instead: a missing OpenMM/compiled
extension/2ala input must fail the gate loudly rather than silently reduce it to
fewer oracles.

Two more CHAMBER-specific tests live at the bottom of this module (not part of
the ``CASES`` table above, since they need bespoke setup):

* ``test_charmm36_ala5_per_force_class_tight_tolerance`` -- the same
  ``charmm36_ala5`` case, but asserting a tight per-force-class tolerance
  (see ``_LOCKED_REL_TOL`` / ``_LOCKED_ABS_FLOOR_KJ_PER_MOL`` below) rather
  than the looser gate above.
* ``test_charmm_chamber_gfcdstrippedmin_matches_reference`` -- full
  per-force-class validation on the real, large ``examples/GfcDstrippedMin.prmtop``
  (CHAMBER + NBFIX + a degenerate all-zero box in its ``.rst7``), including the
  two NBFIX-carrying classes: ``NonbondedForce`` (charges + 1-4/exclusion
  exceptions) and ``CustomNonbondedForce`` (the off-diagonal LJ table built by
  ``OpenMMContext::createCustomNonbondedForce``). GBSA stays off on both sides
  (``use_gbsa_obc2=False`` / ``implicitSolvent=None``): validating GBSA-OBC2 on
  this system is a separate, out-of-scope concern, deliberately not exercised
  here so this test isolates the NBFIX gap it targets.

Both of the above assert against a **locked mixed-precision tolerance**
(``_LOCKED_REL_TOL`` / ``_LOCKED_ABS_FLOOR_KJ_PER_MOL``) instead of an
``xfail``: the C++ engine's OpenMM ``Context`` is built on the CUDA platform
with the compile-time default ``Precision='mixed'`` (``OpenMMContext.cpp``,
``registerPlatform``) -- effectively single-precision force/energy kernels --
which caps agreement with the double-precision ``Reference`` platform used
for the comparison at roughly float32 precision for an individual force
class, regardless of AMBER vs CHAMBER. This is a deliberate, kept choice (not
a bug to fix): CUDA mixed precision is faster and the accuracy loss is
confined to individual per-force-class bookkeeping, not the total energy used
for actual sampling. Measured worst-case per-class relative diff across two
runs of every case in this module (implicit/explicit_pme/vacuum/
charmm36_ala5/GfcDstrippedMin; single-point energy evaluation on a fixed
CUDA-mixed kernel path is deterministic run-to-run here, but sampling runs on
different GPUs/drivers could shift the floor slightly) was ``1.398e-05``
(``HarmonicBondForce`` on the small ``2ala`` systems, where the force class's
own magnitude is tiny -- ~0.09 kJ/mol -- so a ~1e-6 kJ/mol rounding wobble
reads as a large *relative* diff despite being physically negligible); the
next-worst was ``7.130e-06`` (``CustomTorsionForce`` on ``GfcDstrippedMin``, a
genuinely large-magnitude force class, i.e. real float32-vs-float64 noise).
``_LOCKED_REL_TOL = 5e-5`` sits ~3.6x above the true observed floor (and ~7x
above the large-magnitude-class floor), with ``_LOCKED_ABS_FLOOR_KJ_PER_MOL =
1e-5`` kJ/mol as a backstop for any future near-zero-magnitude force class
where the relative metric alone would be meaningless. A class passes if
EITHER bound is met, matching the abs-OR-rel convention already used by
``compare_by_force_group``.
"""

from __future__ import annotations

import math
import os
import pathlib

import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))

# The comparison needs native OpenMM (the reference is built with
# openmm.app.AmberPrmtopFile, not ParmEd) and the compiled robosample extension
# (the .so built by `cmake --build --preset cuda-release`). Bare dev run: skip --
# don't error -- when either is missing. Under the authoritative gate: import for
# real, so a missing dependency is a hard collection-time failure (not a silent skip).
if _REQUIRE_OPENMM:
    import openmm  # noqa: F401
    import robosample
else:
    pytest.importorskip("openmm")
    robosample = pytest.importorskip("robosample")
from robosample import openmm_validation  # noqa: E402  (after importorskip)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "examples" / "2ala"
EXAMPLES_DIR = REPO_ROOT / "examples"
CHARMM_DIR = REPO_ROOT / "examples" / "charmm"
GFC_PRMTOP = REPO_ROOT / "examples" / "GfcDstrippedMin.prmtop"
GFC_RST7 = REPO_ROOT / "examples" / "GfcDstrippedMin.rst7"

# (test id, data dir, prmtop, rst7, load_amber kwargs). Mirrors the run.py
# invocations:
#   run.py 2ala.implicit examples/2ala/2ala.implicit.prmtop ...  -> implicit
#   run.py 2ala.tip3p    examples/2ala/2ala.tip3p.prmtop    ...  -> explicit PME
# plus a gas-phase case (implicit inputs, GBSA forced off) and a CHAMBER case
# (CHARMM36 CMAP/Urey-Bradley/impropers, gas phase -- see the module docstring).
CASES = [
    ("implicit", DATA_DIR, "2ala.implicit.prmtop", "2ala.implicit.rst7", {}),
    ("explicit_pme", DATA_DIR, "2ala.tip3p.prmtop", "2ala.tip3p.rst7", {}),
    (
        "vacuum",
        DATA_DIR,
        "2ala.implicit.prmtop",
        "2ala.implicit.rst7",
        {"use_gbsa_obc2": False},
    ),
    (
        "charmm36_ala5",
        CHARMM_DIR,
        "charmm36_ala5.prmtop",
        "charmm36_ala5.rst7",
        {"use_gbsa_obc2": False},
    ),
    # Tier-2 torsion-conformational fixtures (docs/specs/ensemble-validation/
    # 30-tier2-*.md touch list): single-point-PE-validate them here before the
    # sampling tests use them. Loaded exactly as test_torsion_conformational.py
    # does -- plain load_amber (default kwargs, GBSA-OBC2 on) -- so this checks
    # the same load path those tests depend on.
    ("ethane", EXAMPLES_DIR, "ethane.prmtop", "ethane.rst7", {}),
    ("butane", EXAMPLES_DIR, "butane.prmtop", "butane.rst7", {}),
    ("2butanol", EXAMPLES_DIR, "2butanol.prmtop", "2butanol.rst7", {}),
]


@pytest.mark.parametrize(
    "name, data_dir, prmtop_name, rst7_name, load_kwargs",
    CASES,
    ids=[c[0] for c in CASES],
)
def test_openmm_potential_energy_matches_reference(
    name: str,
    data_dir: pathlib.Path,
    prmtop_name: str,
    rst7_name: str,
    load_kwargs: dict,
) -> None:
    """C++ single-point PE must match the native OpenMM reference for each model."""
    prmtop = data_dir / prmtop_name
    rst7 = data_dir / rst7_name
    if _REQUIRE_OPENMM:
        assert prmtop.exists() and rst7.exists(), (
            f"example inputs not found: {prmtop} / {rst7} -- the authoritative "
            f"gate (ROBOSAMPLE_REQUIRE_OPENMM=1) requires {data_dir} to be present"
        )
    elif not prmtop.exists() or not rst7.exists():
        pytest.skip(f"example inputs not found: {prmtop} / {rst7}")

    context = robosample.Context(name, 0, robosample.AmberDihedralClassifier())
    # Solvent model is auto-detected from the box; `use_gbsa_obc2` (when passed)
    # overrides the implicit-solvent force, matching run.py's load_amber call.
    context.load_amber(str(prmtop), str(rst7), **load_kwargs)
    # Keep whole molecules across worlds (never wrap); harmless for the box-less
    # cases. This mirrors run.py and does not affect energies.
    context.set_enforce_periodic_box(False)

    # Raises ValueError on any total- or per-group mismatch beyond tolerance.
    ok, cpp_pe, ref_pe, _per_class = openmm_validation.compare_by_force_group(
        context, str(prmtop), str(rst7), platform_name="Reference"
    )

    assert ok, (
        f"{name}: C++ OpenMM PE = {cpp_pe:.6f} kJ/mol, "
        f"reference = {ref_pe:.6f} kJ/mol"
    )


def _print_and_max_relative_diff(
    cpp_by_class: dict[str, float], ref_by_class: dict[str, float], header: str
) -> float:
    """Print a per-class relative-diff breakdown; return the worst relative diff.

    Shared by the two locked-tolerance tests below. Every class in
    *cpp_by_class* is required to also be in *ref_by_class* (and vice versa) --
    a class present on only one side is itself a hard failure, distinct from
    (and more severe than) a numeric mismatch.
    """
    names = sorted(set(cpp_by_class) | set(ref_by_class))
    print(header)
    print(f"  {'force class':<20} {'C++':>16} {'reference':>16} {'rel diff':>12}")
    worst_rel = 0.0
    for name in names:
        c = cpp_by_class.get(name, math.nan)
        r = ref_by_class.get(name, math.nan)
        assert not math.isnan(c) and not math.isnan(r), (
            f"force class {name!r} present on only one side "
            f"(cpp={cpp_by_class.get(name)!r}, ref={ref_by_class.get(name)!r})"
        )
        rel = abs(c - r) / max(abs(r), 1e-12)
        worst_rel = max(worst_rel, rel)
        print(f"  {name:<20} {c:>16.6f} {r:>16.6f} {rel:>12.3e}")
    print(f"  worst per-class relative diff: {worst_rel:.3e}")
    return worst_rel


# Locked CUDA-mixed-precision per-force-class tolerance (see the module
# docstring for the measured floor and the margin behind these numbers). A
# force class passes if EITHER bound holds -- same abs-OR-rel convention as
# ``compare_by_force_group``'s own (looser) default gate.
_LOCKED_REL_TOL = 5e-5
_LOCKED_ABS_FLOOR_KJ_PER_MOL = 1e-5


def _assert_locked_tolerance(
    cpp_by_class: dict[str, float], ref_by_class: dict[str, float], header: str
) -> None:
    """Print the per-class breakdown, then assert every class is within the
    locked mixed-precision tolerance (``_LOCKED_REL_TOL`` OR
    ``_LOCKED_ABS_FLOOR_KJ_PER_MOL``).

    Replaces the previous literal ``<=1e-6``-relative bar (which was never
    actually met by the CUDA 'mixed'-precision build and was marked
    ``xfail(strict=False)``): that bar was tighter than the platform can
    deliver and gave no real regression signal. This one is calibrated to the
    measured noise floor with margin, so it PASSES today and FAILS loudly on
    a genuine regression (e.g. GBSA/NBFIX physics breaking, not just precision
    wobble).
    """
    worst_rel = _print_and_max_relative_diff(cpp_by_class, ref_by_class, header)
    failures = []
    for name in sorted(set(cpp_by_class) | set(ref_by_class)):
        c, r = cpp_by_class[name], ref_by_class[name]
        abs_diff = abs(c - r)
        rel_diff = abs_diff / max(abs(r), 1e-12)
        if abs_diff > _LOCKED_ABS_FLOOR_KJ_PER_MOL and rel_diff > _LOCKED_REL_TOL:
            failures.append(f"{name} (abs={abs_diff:.3e}, rel={rel_diff:.3e})")
    assert not failures, (
        f"per-class mixed-precision tolerance exceeded (need rel<={_LOCKED_REL_TOL:.1e} "
        f"or abs<={_LOCKED_ABS_FLOOR_KJ_PER_MOL:.1e} kJ/mol; worst overall rel "
        f"{worst_rel:.3e}): " + "; ".join(failures)
    )


def test_charmm36_ala5_per_force_class_tight_tolerance() -> None:
    """Same case as ``charmm36_ala5`` above, at the locked mixed-precision bar.

    See the module docstring for ``_LOCKED_REL_TOL`` /
    ``_LOCKED_ABS_FLOOR_KJ_PER_MOL`` and the measured floor behind them.
    """
    prmtop = CHARMM_DIR / "charmm36_ala5.prmtop"
    rst7 = CHARMM_DIR / "charmm36_ala5.rst7"
    if _REQUIRE_OPENMM:
        assert prmtop.exists() and rst7.exists()
    elif not prmtop.exists() or not rst7.exists():
        pytest.skip(f"example inputs not found: {prmtop} / {rst7}")

    context = robosample.Context(
        "charmm36_ala5_tight", 0, robosample.AmberDihedralClassifier()
    )
    context.load_amber(str(prmtop), str(rst7), use_gbsa_obc2=False)
    context.set_enforce_periodic_box(False)

    # ok=True here only means the DEFAULT (abs 1e-3 / rel 1e-4) tolerance was
    # met -- compare_by_force_group is reused purely to build the per_class
    # breakdown; the assertion of interest is the locked one below.
    _ok, _cpp_pe, _ref_pe, per_class = openmm_validation.compare_by_force_group(
        context, str(prmtop), str(rst7), platform_name="Reference"
    )
    cpp_by_class = {name: c for name, (c, _r) in per_class.items()}
    ref_by_class = {name: r for name, (_c, r) in per_class.items()}

    _assert_locked_tolerance(
        cpp_by_class,
        ref_by_class,
        "charmm36_ala5 per-force-class tight-tolerance breakdown [kJ/mol]:",
    )


def test_charmm_chamber_gfcdstrippedmin_matches_reference() -> None:
    """Full CHAMBER + NBFIX validation on the real GfcDstrippedMin system.

    ``examples/GfcDstrippedMin.prmtop`` is a real CHAMBER (CHARMM36) system
    with off-diagonal NBFIX LJ pairs and a degenerate (all-zero) box line in
    its ``.rst7``. Two facts about the tooling shape this test:

    1. ``AmberInpcrdFile`` cannot parse this ``.rst7``: the box line is
       ``0 0 0 90 90 90`` (the prmtop's own ``IFBOX == 0``, i.e. there is no
       real box -- this is degenerate data, not periodicity), and OpenMM's
       box-vector construction divides by the (zero) box volume, raising
       ``ZeroDivisionError``. Coordinates are read directly instead, with
       ``amber_loader.read_amber_coordinates`` (already unit-tested
       elsewhere).
    2. NBFIX (off-diagonal LJ) is handled by
       ``OpenMMContext::createCustomNonbondedForce``, which looks up
       ``sqrt(A)``/``B`` coefficients in a per-LJ-type-pair table
       (``SystemTopology.aCoef``/``bCoef``) keyed by each atom's 0-based LJ
       type (``SystemTopology.atomsNonbondedIndex``), and excludes every pair
       the main ``NonbondedForce`` already accounts for (1-4 exceptions +
       1-2/1-3 exclusions) so LJ is never double-counted. This test compares
       that force, and the charge-only ``NonbondedForce`` it complements,
       directly against native OpenMM -- no workaround needed anymore.

    GBSA stays off on both sides (``use_gbsa_obc2=False`` / ``implicitSolvent=
    None``): validating GBSA-OBC2 on this system is a separate, out-of-scope
    concern (see the module docstring), deliberately not exercised here.
    """
    if _REQUIRE_OPENMM:
        assert GFC_PRMTOP.exists() and GFC_RST7.exists()
    elif not GFC_PRMTOP.exists() or not GFC_RST7.exists():
        pytest.skip(f"example inputs not found: {GFC_PRMTOP} / {GFC_RST7}")

    import openmm.app as app
    import openmm as mm
    import openmm.unit as unit

    from robosample.amber_loader import read_amber_coordinates

    # ---- C++ side ----------------------------------------------------------
    context = robosample.Context(
        "gfc_chamber", 0, robosample.AmberDihedralClassifier()
    )
    context.load_amber(str(GFC_PRMTOP), str(GFC_RST7), use_gbsa_obc2=False)
    assert context.system_topology.has_nb_fix, (
        "examples/GfcDstrippedMin.prmtop is expected to carry NBFIX pairs "
        "(prmtop_reader.has_nbfix_fast); if this ever fails the system "
        "changed and the NBFIX path in this test is no longer exercised."
    )
    context.set_enforce_periodic_box(False)
    context.set_separate_force_groups(True)
    assert context.initialize_openmm()
    _cpp_total, cpp_groups_raw = context.calc_openmm_potential_energy_by_group()
    cpp_by_class: dict[str, float] = {}
    for entry in cpp_groups_raw:
        name = str(entry.name)
        cpp_by_class[name] = cpp_by_class.get(name, 0.0) + float(entry.energy)

    # ---- Native OpenMM reference (coordinates read directly, no AmberInpcrdFile) --
    prmtop = app.AmberPrmtopFile(str(GFC_PRMTOP))
    coords_nm, _box = read_amber_coordinates(str(GFC_RST7))
    system = prmtop.createSystem(
        nonbondedMethod=app.NoCutoff,
        constraints=None,
        implicitSolvent=None,
        rigidWater=False,
        removeCMMotion=False,
    )
    all_forces = []
    next_group = 0
    for force in system.getForces():
        force.setForceGroup(next_group)
        all_forces.append((next_group, force))
        next_group += 1

    integrator = mm.VerletIntegrator(0.001)
    platform = mm.Platform.getPlatformByName("Reference")
    mm_context = mm.Context(system, integrator, platform)
    mm_context.setPositions(unit.Quantity(coords_nm, unit.nanometer))
    ref_by_class: dict[str, float] = {}
    for group, force in all_forces:
        e = (
            mm_context.getState(getEnergy=True, groups={group})
            .getPotentialEnergy()
            .value_in_unit(unit.kilojoule_per_mole)
        )
        name = type(force).__name__
        ref_by_class[name] = ref_by_class.get(name, 0.0) + e

    # ---- Compare -------------------------------------------------------------
    names = sorted(set(cpp_by_class) | set(ref_by_class))
    assert names, "no force classes found on either side"
    for name in names:
        assert name in cpp_by_class and name in ref_by_class, (
            f"force class {name!r} present on only one side "
            f"(cpp={cpp_by_class.get(name)!r}, ref={ref_by_class.get(name)!r})"
        )

    # Locked mixed-precision per-class tolerance (see the module docstring for
    # the measured floor and margin behind _LOCKED_REL_TOL /
    # _LOCKED_ABS_FLOOR_KJ_PER_MOL). Covers the two NBFIX-carrying classes
    # (NonbondedForce, CustomNonbondedForce) as well as the CHAMBER-specific
    # ones (CMAPTorsionForce, CustomTorsionForce).
    _assert_locked_tolerance(
        cpp_by_class,
        ref_by_class,
        "GfcDstrippedMin per-force-class breakdown [kJ/mol] (CUDA mixed vs. "
        "Reference double precision):",
    )
