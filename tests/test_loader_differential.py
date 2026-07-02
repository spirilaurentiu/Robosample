"""Differential gate for the AMBER/CHAMBER loader rewrite.

Spec: ``docs/specs/fast-amber-loader.md``. For each example system, this
compares a *live* ``Context.load_amber`` load against a golden
``SystemTopology`` snapshot pickled under ``tests/fixtures/loader_golden/``
(regenerated via ``tests/fixtures/loader_golden/_generate.py`` -- see that
file's docstring for when to run it).

Right now -- before the loader is rewritten -- both sides of the comparison
run through the SAME ``Context.load_amber`` code path, so this is a
self-consistency check that must pass trivially: it locks the
``SystemTopology`` field contract (docs/specs/fast-amber-loader.md §4) before
any refactor. Once the loader is rewritten in place, the golden fixtures stay
pinned to the pre-rewrite output and this becomes the real differential gate:
every array must still match exactly (int/bool/str/enum) or to <=1e-9
relative (float), and ``df_bonds`` must match exactly.

Per the spec (§2): identical ``SystemTopology`` arrays imply identical
energies by construction (a sum over bonded/nonbonded terms indexed by atom),
so this array-level comparison -- not a numeric energy comparison -- is the
primary correctness oracle for the rewrite. ``tests/test_openmm_potential_
energy.py`` remains the energy-level backstop gate.
"""

from __future__ import annotations

import os
import pickle
import time

import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))
# Reuses the same authoritative-gate flag as test_openmm_potential_energy.py
# (the only such flag defined in this suite): under the authoritative gate
# (`nox -s tests`) a missing compiled extension or example input is a hard
# collection-time failure rather than a silent skip.
if _REQUIRE_OPENMM:
    import robosample  # noqa: F401
else:
    pytest.importorskip("robosample")

from loader_differential_lib import (  # noqa: E402  (after importorskip)
    CASES,
    FIXTURE_DIR,
    LARGEST_CASE,
    Case,
    capture_topology,
    compare_topologies,
    load_context,
)


def _require_inputs(case: Case) -> None:
    ok = case.prmtop_path.exists() and case.rst7_path.exists()
    if _REQUIRE_OPENMM:
        assert ok, (
            f"{case.name}: example inputs not found ({case.prmtop_path} / "
            f"{case.rst7_path}) -- the authoritative gate requires examples/ "
            "to be present"
        )
    elif not ok:
        pytest.skip(f"{case.name}: example inputs not found")


@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_loader_matches_golden(case: Case) -> None:
    """Every ``SystemTopology`` array + ``df_bonds`` must match the pinned golden output."""
    _require_inputs(case)

    golden_path = FIXTURE_DIR / f"{case.name}.pkl"
    if _REQUIRE_OPENMM:
        assert golden_path.exists(), (
            f"missing golden fixture {golden_path}; regenerate via "
            "tests/fixtures/loader_golden/_generate.py"
        )
    elif not golden_path.exists():
        pytest.skip(f"{case.name}: golden fixture not found ({golden_path})")

    with open(golden_path, "rb") as fh:
        golden = pickle.load(fh)

    context = load_context(case)
    fresh = capture_topology(context)

    compare_topologies(golden, fresh)


def test_load_time_smoke(capsys) -> None:
    """Informational load-time check on the largest available example.

    Not a hard performance gate (per docs/specs/fast-amber-loader.md §7) --
    only asserts the load succeeds and logs the wall-clock time so
    regressions (or improvements) are visible in the test output.
    """
    _require_inputs(LARGEST_CASE)

    t0 = time.time()
    context = load_context(LARGEST_CASE)
    elapsed = time.time() - t0

    with capsys.disabled():
        print(
            f"\n[load-time smoke] {LARGEST_CASE.name}: {elapsed:.3f}s for "
            f"{context.system_topology.num_atoms} atoms, "
            f"{context.system_topology.num_molecules} molecules"
        )

    assert context.system_topology.num_atoms > 0
