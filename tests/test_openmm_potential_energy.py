"""Validate the C++ OpenMM potential energy against a Python/OpenMM reference.

These tests are the pytest form of the ``--validate`` path in
``python/robosample/run.py``: for each solvent model they build a
``robosample.Context``, let ``load_amber`` auto-detect the model from the
periodic box (OpenMM-style), and then diff the C++ single-point energy against a
ParmEd + OpenMM reference built from the same ``prmtop``/``rst7``.

``openmm_validation.compare_by_force_group`` is authoritative: it compares the
total PE (and every force group) and *raises* ``ValueError`` on any mismatch, so
a passing call is the assertion. The reference runs on the OpenMM ``Reference``
platform for a deterministic, full-precision comparison -- exactly as
``run.py`` invokes it.

Three models are covered:

* **implicit** -- ``2ala.implicit`` has no box, so ``load_amber`` selects
  GBSA-OBC2 (implicit solvent, non-periodic).
* **explicit** -- ``2ala.tip3p`` carries a periodic box, so ``load_amber``
  selects PME under PBC and turns GBSA off.
* **vacuum**   -- the same box-less ``2ala.implicit`` files loaded with
  ``use_gbsa_obc2=False``: no implicit solvent, no periodicity (gas phase).

On a bare dev run the tests skip cleanly when the compiled extension,
OpenMM/ParmEd, or the example inputs are unavailable, rather than failing
collection. Under the AUTHORITATIVE gate (``nox -s tests``, which exports
``ROBOSAMPLE_REQUIRE_OPENMM=1``) these same conditions are HARD FAILURES
instead: a missing OpenMM/ParmEd/compiled extension/2ala input must fail the
gate loudly rather than silently reduce it to fewer oracles.
"""

from __future__ import annotations

import os
import pathlib

import pytest

_REQUIRE_OPENMM = bool(os.environ.get("ROBOSAMPLE_REQUIRE_OPENMM"))

# The comparison needs the OpenMM/ParmEd reference stack and the compiled
# robosample extension (the .so built by `cmake --build --preset cuda-release`).
# Bare dev run: skip -- don't error -- when any of them is missing. Under the
# authoritative gate: import for real, so a missing dependency is a hard
# collection-time failure (not a silent skip).
if _REQUIRE_OPENMM:
    import openmm  # noqa: F401
    import parmed  # noqa: F401
    import robosample
else:
    pytest.importorskip("openmm")
    pytest.importorskip("parmed")
    robosample = pytest.importorskip("robosample")
from robosample import openmm_validation  # noqa: E402  (after importorskip)

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA_DIR = REPO_ROOT / "examples" / "2ala"

# (test id, prmtop, rst7, load_amber kwargs). Mirrors the run.py invocations:
#   run.py 2ala.implicit examples/2ala/2ala.implicit.prmtop ...  -> implicit
#   run.py 2ala.tip3p    examples/2ala/2ala.tip3p.prmtop    ...  -> explicit PME
# plus a gas-phase case (implicit inputs, GBSA forced off).
CASES = [
    ("implicit", "2ala.implicit.prmtop", "2ala.implicit.rst7", {}),
    ("explicit_pme", "2ala.tip3p.prmtop", "2ala.tip3p.rst7", {}),
    ("vacuum", "2ala.implicit.prmtop", "2ala.implicit.rst7", {"use_gbsa_obc2": False}),
]


@pytest.mark.parametrize(
    "name, prmtop_name, rst7_name, load_kwargs",
    CASES,
    ids=[c[0] for c in CASES],
)
def test_openmm_potential_energy_matches_reference(
    name: str,
    prmtop_name: str,
    rst7_name: str,
    load_kwargs: dict,
) -> None:
    """C++ single-point PE must match the ParmEd+OpenMM reference for each model."""
    prmtop = DATA_DIR / prmtop_name
    rst7 = DATA_DIR / rst7_name
    if _REQUIRE_OPENMM:
        assert prmtop.exists() and rst7.exists(), (
            f"example inputs not found: {prmtop} / {rst7} -- the authoritative "
            "gate (ROBOSAMPLE_REQUIRE_OPENMM=1) requires examples/2ala/ to be present"
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
