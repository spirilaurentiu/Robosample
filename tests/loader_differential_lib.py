"""Shared helpers for the AMBER/CHAMBER loader differential gate.

Used by both ``tests/test_loader_differential.py`` (the pytest gate) and
``tests/fixtures/loader_golden/_generate.py`` (the one-off golden-fixture
generator). This module has no ``test_`` prefix, so pytest never collects it
as a test file on its own.

See ``docs/specs/fast-amber-loader.md`` §2/§7 for why the comparison is
array-level rather than energy-level: identical ``SystemTopology`` arrays
imply identical energies by construction (a sum over bonded/nonbonded terms
indexed by atom), so this is the primary correctness oracle for the loader
rewrite.
"""

from __future__ import annotations

import pathlib
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

REPO_ROOT = pathlib.Path(__file__).resolve().parents[1]
FIXTURE_DIR = pathlib.Path(__file__).resolve().parent / "fixtures" / "loader_golden"


@dataclass(frozen=True)
class Case:
    name: str
    prmtop: str  # path relative to REPO_ROOT
    rst7: str  # path relative to REPO_ROOT
    kwargs: dict = field(default_factory=dict)

    @property
    def prmtop_path(self) -> pathlib.Path:
        return REPO_ROOT / self.prmtop

    @property
    def rst7_path(self) -> pathlib.Path:
        return REPO_ROOT / self.rst7


# Small/medium example systems covering: no-box implicit (ala-dipeptide,
# 2ala.implicit), explicit-solvent PME with a periodic box (2ala.tip3p, which
# is also the existing test_openmm_potential_energy.py PME case), a
# multi-molecule / duplicate-molecule system that exercises parm.split()
# dedup (2ala.2ala, host-guest), nucleic acids (DNA/RNA), a branched
# carbohydrate (glycan), and several single-residue/small-molecule systems.
# All are fast to load (< 0.1s each with the pre-rewrite loader), so the full
# suite is cheap even though it is the authoritative correctness gate.
#
# NOT included: GfcDstrippedMin.prmtop (the one CHAMBER example under
# examples/) -- it fails with ZeroDivisionError on the CURRENT,
# pre-rewrite loader already (verified independently of this work), so it
# cannot serve as a golden reference; out of scope here, flagged for
# separate investigation. b1-1n.prmtop (the largest example, 88628 atoms) is
# deliberately excluded from the exact-array gate too -- its golden fixture
# would be several MB -- and is instead used only by the informational
# load-time smoke check (LARGEST_CASE below); every code path it exercises is
# already covered by the smaller CASES.
CASES: tuple[Case, ...] = (
    Case("ala-dipeptide", "examples/ala-dipeptide.prmtop", "examples/ala-dipeptide.rst7"),
    Case(
        "2ala.implicit",
        "examples/2ala/2ala.implicit.prmtop",
        "examples/2ala/2ala.implicit.rst7",
    ),
    Case("2ala.tip3p", "examples/2ala/2ala.tip3p.prmtop", "examples/2ala/2ala.tip3p.rst7"),
    Case("10ala", "examples/10ala.prmtop", "examples/10ala.rst7"),
    Case("host-guest", "examples/host-guest.prmtop", "examples/host-guest.rst7"),
    Case("DNA-AAGCTA", "examples/DNA-AAGCTA.prmtop", "examples/DNA-AAGCTA.rst7"),
    Case("RNA-UAGCUU", "examples/RNA-UAGCUU.prmtop", "examples/RNA-UAGCUU.rst7"),
    Case("glycan", "examples/glycan.prmtop", "examples/glycan.rst7"),
    Case("2ala.2ala", "examples/2ala.2ala.prmtop", "examples/2ala.2ala.rst7"),
    Case("2but", "examples/2but.prmtop", "examples/2but.rst7"),
    Case("phe", "examples/phe.prmtop", "examples/phe.rst7"),
    Case("pro", "examples/pro.prmtop", "examples/pro.rst7"),
    Case("2M6A", "examples/2M6A.prmtop", "examples/2M6A.rst7"),
    Case("1APQ", "examples/1APQ.prmtop", "examples/1APQ.rst7"),
    # 4-point TIP4P-Ew water (1 virtual site): the only case exercising the
    # ExtraPoint/virtual-site extraction path (docs/specs/fast-amber-loader.md
    # Step 4b). Non-periodic (IFBOX==0, single isolated water, no box line).
    Case("tip4pew", "water/tip4pew.prmtop", "water/tip4pew.rst7"),
)

# Largest available example (88628 atoms); used only for the informational
# load-time smoke check, not the exact-array gate (see CASES docstring above).
LARGEST_CASE = Case("b1-1n", "examples/b1-1n.prmtop", "examples/b1-1.equil.rst7")


def load_context(case: Case):
    """Build a Context and load ``case`` via the CURRENT ``Context.load_amber``."""
    import robosample as rb

    context = rb.Context(case.name, 0, rb.AmberDihedralClassifier())
    context.load_amber(str(case.prmtop_path), str(case.rst7_path), **case.kwargs)
    return context


def capture_topology(context) -> dict:
    """Snapshot every ``SystemTopology`` field plus ``df_bonds``.

    Introspects the live pybind11 object (``dir()``) rather than a hardcoded
    field list, so a field added on the C++ side is captured automatically
    instead of silently skipped (see CODING_RULES Rule 11, fail loud).
    """
    st = context.system_topology
    names = [n for n in dir(st) if not n.startswith("_")]
    data = {name: getattr(st, name) for name in names}
    data["df_bonds"] = context.df_bonds.copy()
    return data


def _compare_field(name: str, gold, fresh, rtol: float, atol: float) -> None:
    if isinstance(gold, list) or isinstance(fresh, list):
        if not (isinstance(gold, list) and isinstance(fresh, list)):
            raise AssertionError(f"{name}: type differs ({type(gold)} vs {type(fresh)})")
        if len(gold) != len(fresh):
            raise AssertionError(f"{name}: length differs ({len(gold)} vs {len(fresh)})")
        if not gold:
            return
        if isinstance(gold[0], float):
            g = np.asarray(gold, dtype=np.float64)
            f = np.asarray(fresh, dtype=np.float64)
            close = np.isclose(g, f, rtol=rtol, atol=atol)
            if not np.all(close):
                bad = np.flatnonzero(~close)
                raise AssertionError(
                    f"{name}: {bad.size} float mismatch(es) (of {len(g)}), "
                    f"first at index {bad[0]}: {g[bad[0]]!r} vs {f[bad[0]]!r}"
                )
        else:
            for i, (g, f) in enumerate(zip(gold, fresh)):
                if g != f:
                    raise AssertionError(f"{name}[{i}]: {g!r} != {f!r}")
    else:
        if isinstance(gold, float):
            if not np.isclose(gold, fresh, rtol=rtol, atol=atol):
                raise AssertionError(f"{name}: {gold!r} != {fresh!r}")
        elif gold != fresh:
            raise AssertionError(f"{name}: {gold!r} != {fresh!r}")


def compare_topologies(
    golden: dict, fresh: dict, *, rtol: float = 1e-9, atol: float = 1e-12
) -> None:
    """Assert every captured field matches between ``golden`` and ``fresh``.

    Exact equality for int/bool/str/enum (scalars and lists); ``<= rtol``
    relative for float lists/scalars; exact equality for ``df_bonds``. Raises
    ``AssertionError`` naming the first offending field/index on any mismatch,
    so a failure is immediately actionable.
    """
    gold_keys = set(golden) - {"df_bonds"}
    fresh_keys = set(fresh) - {"df_bonds"}
    if gold_keys != fresh_keys:
        raise AssertionError(
            "SystemTopology field set changed: "
            f"missing={sorted(gold_keys - fresh_keys)} extra={sorted(fresh_keys - gold_keys)}"
        )
    for name in sorted(gold_keys):
        _compare_field(name, golden[name], fresh[name], rtol, atol)

    pd.testing.assert_frame_equal(
        golden["df_bonds"].reset_index(drop=True),
        fresh["df_bonds"].reset_index(drop=True),
    )
