"""Regenerate the golden ``SystemTopology`` fixtures for the loader
differential gate (``tests/test_loader_differential.py``).

This is a manual dev tool, NOT part of the pytest suite (no ``test_`` prefix,
so pytest never collects it). Run it by hand:

    python3 tests/fixtures/loader_golden/_generate.py

ONLY run this after deliberately reviewing a change to ``Context.load_amber``'s
output (docs/specs/fast-amber-loader.md). The checked-in fixtures pin the
*pre-rewrite* (parmed-based, per-atom-Python-loop) loader's output, captured
before docs/specs/fast-amber-loader.md Step 2 landed; regenerating them
against a broken loader would silently defeat the differential gate.
"""

from __future__ import annotations

import pathlib
import pickle
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[2]))  # tests/

from loader_differential_lib import CASES, FIXTURE_DIR, capture_topology, load_context  # noqa: E402


def main() -> None:
    FIXTURE_DIR.mkdir(parents=True, exist_ok=True)
    for case in CASES:
        context = load_context(case)
        data = capture_topology(context)
        out_path = FIXTURE_DIR / f"{case.name}.pkl"
        with open(out_path, "wb") as fh:
            pickle.dump(data, fh)
        print(f"wrote {out_path} ({context.system_topology.num_atoms} atoms)")


if __name__ == "__main__":
    main()
