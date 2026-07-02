"""Dump the PORT's own ``SystemTopology`` (docs/specs/robotics-oracle-
differential.md Scope B, §4.2) for a molecule into a C++-loadable ``.npz``.

This is the disasm-side half of the Scope-B correspondence: it reuses the
SAME loading path ``test_loader_differential.py``/``loader_golden`` already
exercise (``Context.load_amber`` -> ``context.system_topology``), and dumps
only the SUBSET of ``SystemTopology`` fields ``World::buildModel``
(``src/World.cpp``) and ``Context::buildFlexibilities`` (``src/Context.cpp``)
actually read: atom positions/mass/bond-degree, bond endpoints/ring-closing
flag, molecule ranges, and the per-molecule root-mobility default. This is
NOT a general-purpose SystemTopology serializer (no force-field arrays) --
Scope B only needs the DECOMPOSITION (§4.2), never OpenMM energetics.

Manual dev tool (no ``test_`` prefix -> never pytest-collected). Run by hand:

    python3 tests/fixtures/robotics_oracle_molecules/_generate_port_topology.py

Requires the disasm ``robo_bindings`` extension to be built and importable
(``python/robosample/``, i.e. ``cuda-tests``/``cuda-release`` already built).
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

REPO_ROOT = pathlib.Path(__file__).resolve().parents[3]
FIXTURE_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(REPO_ROOT / "python"))

import robosample as rb  # noqa: E402


CASES = {
    "10ala": (REPO_ROOT / "examples" / "10ala.prmtop", REPO_ROOT / "examples" / "10ala.rst7"),
    "1APQ": (REPO_ROOT / "examples" / "1APQ.prmtop", REPO_ROOT / "examples" / "1APQ.rst7"),
}


def dump_topology(name: str, prmtop: pathlib.Path, rst7: pathlib.Path) -> None:
    context = rb.Context(name, 0, rb.AmberDihedralClassifier())
    context.load_amber(str(prmtop), str(rst7))
    st = context.system_topology

    root_mobilities = np.asarray([int(x) for x in st.root_mobilities], dtype=np.int32)
    bonds_ring_closing = np.asarray([bool(x) for x in st.bonds_ring_closing], dtype=np.bool_)

    out = FIXTURE_DIR / f"{name}.systopo.npz"
    np.savez(
        out,
        num_molecules=np.asarray([st.num_molecules], dtype=np.int32),
        atoms_begin=np.asarray(st.atoms_begin, dtype=np.int32),
        atoms_end=np.asarray(st.atoms_end, dtype=np.int32),
        root_mobilities=root_mobilities,
        num_atoms=np.asarray([st.num_atoms], dtype=np.int32),
        atoms_x=np.asarray(st.atoms_x, dtype=np.float64),
        atoms_y=np.asarray(st.atoms_y, dtype=np.float64),
        atoms_z=np.asarray(st.atoms_z, dtype=np.float64),
        atoms_mass=np.asarray(st.atoms_mass, dtype=np.float64),
        atoms_num_bonds_involved=np.asarray(st.atoms_num_bonds_involved, dtype=np.int32),
        atoms_prmtop_index=np.asarray(st.atoms_prmtop_index, dtype=np.int32),
        num_bonds=np.asarray([st.num_bonds], dtype=np.int32),
        bonds_i=np.asarray(st.bonds_i, dtype=np.int32),
        bonds_j=np.asarray(st.bonds_j, dtype=np.int32),
        bonds_ring_closing=bonds_ring_closing,
    )
    print(f"wrote {out} (numAtoms={st.num_atoms}, numBonds={st.num_bonds}, "
          f"numMolecules={st.num_molecules}, ringClosing={int(bonds_ring_closing.sum())})")


def main() -> None:
    for name, (prmtop, rst7) in CASES.items():
        dump_topology(name, prmtop, rst7)


if __name__ == "__main__":
    main()
