"""Export per-bond flexibility metadata for the FlexWizard PyMOL plugin.

This bridges Robosample's AMBER dihedral detection
(:mod:`robosample.amber_dihedral_types`, surfaced as ``context.df_bonds``) and
bond-rigidity classification (:mod:`robosample.bond_util`) to the standalone
``tools/FlexWizard.py`` PyMOL wizard.

The wizard runs inside PyMOL, whose Python interpreter generally has neither
``parmed`` nor Robosample available, so the two sides communicate through a
small text sidecar instead of a live call.  The sidecar is keyed by **0-based
prmtop atom indices** -- exactly what PyMOL exposes (``atom id - 1``) and what
the FlexWizard ``.flex`` writer already emits -- so the wizard can pre-color
rigid / ring-closing bonds and offer named dihedral selections (phi, psi, chi,
...) with no dependency on this package.

The index bridge
----------------
``df_bonds`` stores bond endpoints as *compound* (BFS) atom indices that are
molecule-local (``atom1_idx`` / ``atom2_idx``) together with the instance's
``molecule_idx``.  Two public topology arrays turn those into the prmtop order
PyMOL sees:

* ``system_topology.atoms_begin[molecule_idx]`` -- the molecule's start in the
  global compound array (the same offset :meth:`Context.build_flexibilities`
  uses), so ``global_compound = atoms_begin[molecule_idx] + atomN_idx``.
* ``system_topology.atoms_prmtop_index[global_compound]`` -- the
  compound -> prmtop permutation, i.e. the 0-based prmtop index of that atom.

Usage
-----
    from robosample import flex_export

    # `context` after context.load_amber(prmtop, inpcrd)
    flex_export.export_flex_meta(context, prmtop_path, "system.flexmeta")

Then, in PyMOL, load the object and use FlexWizard's "Load Bond Metadata".
"""

from __future__ import annotations

import os

import parmed as pmd

from . import bond_util
from .amber_dihedral_types import DihedralType

SIDECAR_VERSION = 2


def _type_name(value) -> str:
    """Resolve a stored ``dihedral_type`` value to a ``DihedralType`` name.

    ``df_bonds`` stores the type as its bare :class:`DihedralType` integer value
    (the enum identity is dropped by pandas), so map it back to the member name
    (e.g. ``2 -> "PROTEIN_PHI"``). User-defined codes (>= 100) that are not enum
    members are passed through as their integer string.
    """
    try:
        int_value = int(value)
    except (TypeError, ValueError):
        return str(value)
    try:
        return DihedralType(int_value).name
    except ValueError:
        return str(int_value)


def compute_flex_meta(context, prmtop_path):
    """Compute FlexWizard bond metadata for a built ``context`` (in memory).

    Parameters
    ----------
    context :
        A Robosample ``Context`` on which ``load_amber`` has already run, so
        that ``context.standard_dihedral_bonds`` and ``context.system_topology``
        are populated.
    prmtop_path :
        The AMBER topology the context was built from.  Reloaded here (typing
        only, no coordinates) to classify each bond's rigidity with
        :func:`robosample.bond_util.is_rigid_bond`.

    Returns
    -------
    tuple[list[tuple[int, int, str, int, int, int]], int]
        ``(rows, n_atoms)``.  Each row is
        ``(atom1, atom2, dihedral_type, is_ring_closing, is_rigid, group)`` with
        0-based prmtop atom indices, ordered by ``(atom1, atom2)``.

        ``dihedral_type`` is the :class:`DihedralType` member name (e.g.
        ``PROTEIN_PHI``, or ``UNKNOWN``).  ``is_ring_closing`` / ``is_rigid`` are
        ``0`` / ``1``.  ``group`` is shared by every copy of the same
        prototype-local bond (repeated chains, waters, ...), so equal-``group``
        bonds can be assigned together ("Apply to all copies").
    """
    topo = context.system_topology
    atoms_begin = list(topo.atoms_begin)
    prmtop_index = list(topo.atoms_prmtop_index)
    n_atoms = len(prmtop_index)

    # instance molecule_idx -> prototype index (from parm.split() grouping).
    instance_to_prototype = dict(context.molecules)

    # Reload for atom types only; coordinates are irrelevant to classification.
    parm = pmd.load_file(str(prmtop_path))

    rows: list[tuple[int, int, str, int, int, int]] = []
    seen: set[tuple[int, int]] = set()
    # (prototype, local_a, local_b) -> dense group id, so every copy of the same
    # prototype-local bond collapses onto one id.
    group_of: dict[tuple[int, int, int], int] = {}
    for row in context.standard_dihedral_bonds.itertuples(index=False):
        mol = int(row.molecule_idx)
        local1, local2 = int(row.atom1_idx), int(row.atom2_idx)
        global_compound_1 = atoms_begin[mol] + local1
        global_compound_2 = atoms_begin[mol] + local2
        p1 = prmtop_index[global_compound_1]
        p2 = prmtop_index[global_compound_2]
        atom1, atom2 = (p1, p2) if p1 < p2 else (p2, p1)
        if (atom1, atom2) in seen:
            continue
        seen.add((atom1, atom2))

        proto_key = (
            instance_to_prototype[mol],
            min(local1, local2),
            max(local1, local2),
        )
        group = group_of.setdefault(proto_key, len(group_of))

        type_name = _type_name(row.dihedral_type)
        ring = 1 if bool(row.is_ring_closing) else 0
        rigid = 1 if bond_util.is_rigid_bond(parm.atoms[atom1], parm.atoms[atom2]) else 0
        rows.append((atom1, atom2, type_name, ring, rigid, group))

    rows.sort()
    return rows, n_atoms


def compute_flex_meta_from_files(
    prmtop_path: str | os.PathLike[str],
    inpcrd_path: str | os.PathLike[str],
    seed: int = 0,
):
    """Build a throwaway ``Context`` from files and return its bond metadata.

    This is the entry point FlexWizard uses to detect dihedrals **on the fly**
    (no sidecar): it constructs a Robosample context in-process, so the caller
    (PyMOL) needs Robosample importable but writes no intermediate file.

    Returns the same ``(rows, n_atoms)`` as :func:`compute_flex_meta`.
    """
    # Imported here so importing this module stays cheap and side-effect free.
    from .amber_dihedral_classifier import AmberDihedralClassifier
    from .context import Context

    context = Context("flexwizard", seed, AmberDihedralClassifier())
    context.load_amber(str(prmtop_path), str(inpcrd_path))
    return compute_flex_meta(context, prmtop_path)


def export_flex_meta(
    context,
    prmtop_path: str | os.PathLike[str],
    out_path: str | os.PathLike[str],
) -> str:
    """Write a FlexWizard bond-metadata sidecar (headless / batch use).

    FlexWizard itself generates metadata in memory via
    :func:`compute_flex_meta_from_files`; this file writer is kept for headless
    or remote-PyMOL setups where Robosample is not importable in the viewer, and
    for precomputing very large systems.  See :func:`compute_flex_meta` for the
    column meaning.
    """
    rows, n_atoms = compute_flex_meta(context, prmtop_path)
    with open(out_path, "w") as f:
        f.write(f"# FlexWizard bond metadata v{SIDECAR_VERSION}\n")
        f.write(f"# n_atoms={n_atoms}\n")
        f.write(
            "# columns: atom1 atom2 dihedral_type is_ring_closing is_rigid group\n"
        )
        for atom1, atom2, type_name, ring, rigid, group in rows:
            f.write(
                f"{atom1}\t{atom2}\t{type_name}\t{ring}\t{rigid}\t{group}\n"
            )

    print(
        f"[flex_export] Wrote {len(rows)} bond records for {n_atoms} atoms "
        f"to {out_path}."
    )
    return str(out_path)
