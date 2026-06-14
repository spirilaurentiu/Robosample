"""z_matrix.py

Z-matrix (BAT-style internal-coordinate tree) construction.

The builder returns four equal-length lists (``z_i, z_j, z_k, z_l``) of *local*
atom indices, with the sentinel ``-1`` in reference positions that do not exist
for the first few rows.  Remapping to compound indices is the caller's job.

Note on ordering
----------------
Row ``r`` is the atom *placed* at step ``r``; ``j/k/l`` are back-references
toward already-placed atoms (i.e. toward the root).  This is the opposite of
the "closest-to-root-first" orientation applied to bonds/angles/torsions, and
it is intentional: the Z-matrix defines each atom relative to atoms placed
before it.  The traversal order here is its own mass-prioritised walk and need
not coincide with the BFS/compound atom order.
"""

from __future__ import annotations

from collections import deque

import networkx as nx
import parmed as pmd

# Sentinel stored in z_j/z_k/z_l when a reference atom does not exist.
ZM_SENTINEL: int = -1


def _sort_atoms_by_mass(atoms: list[pmd.Atom]) -> list[pmd.Atom]:
    """Sort by descending mass; ascending atom index breaks ties (reproducible)."""
    return sorted(atoms, key=lambda a: (-a.mass, a.idx))


def build_z_matrix(
    molecule: pmd.Structure,
    acyclic_graph: nx.Graph,
    root_local_index: int,
) -> tuple[list[int], list[int], list[int], list[int]]:
    """
    Build the Z-matrix rooted at *root_local_index*.

    Returns
    -------
    (z_i, z_j, z_k, z_l) : four lists of length ``num_atoms`` (local indices).

    Raises
    ------
    ValueError
        On an empty molecule, a root absent from the graph, an unroutable root
        triplet, a disconnected graph, or a stalled traversal.
    """
    atom_by_idx: dict[int, pmd.Atom] = {a.idx: a for a in molecule.atoms}
    n: int = len(molecule.atoms)

    if n == 0:
        raise ValueError("Molecule contains no atoms; cannot build Z-matrix.")
    if root_local_index not in acyclic_graph:
        raise ValueError(
            "Root atom (local idx=%d) is not a node in the acyclic graph."
            % root_local_index
        )

    root: pmd.Atom = atom_by_idx[root_local_index]

    z_i: list[int] = []
    z_j: list[int] = []
    z_k: list[int] = []
    z_l: list[int] = []

    tree_adj: dict[int, list[pmd.Atom]] = {
        a.idx: _sort_atoms_by_mass(
            [atom_by_idx[nb] for nb in acyclic_graph.neighbors(a.idx)]
        )
        for a in molecule.atoms
    }

    # -- n == 1 ----------------------------------------------------------
    if n == 1:
        z_i.append(root.idx)
        z_j.append(ZM_SENTINEL)
        z_k.append(ZM_SENTINEL)
        z_l.append(ZM_SENTINEL)
        return z_i, z_j, z_k, z_l

    root_nbrs: list[pmd.Atom] = tree_adj[root.idx]
    if not root_nbrs:
        raise ValueError(
            "Root atom (local idx=%d) has no tree-neighbors in a %d-atom "
            "molecule." % (root.idx, n)
        )

    # -- n == 2 ----------------------------------------------------------
    if n == 2:
        second: pmd.Atom = root_nbrs[0]
        z_i.extend([root.idx, second.idx])
        z_j.extend([ZM_SENTINEL, root.idx])
        z_k.extend([ZM_SENTINEL, ZM_SENTINEL])
        z_l.extend([ZM_SENTINEL, ZM_SENTINEL])
        return z_i, z_j, z_k, z_l

    # -- Root triplet for n >= 3 -----------------------------------------
    initial_atom: pmd.Atom = root

    non_term_root: list[pmd.Atom] = [a for a in root_nbrs if len(a.bond_partners) > 1]
    second_atom: pmd.Atom = (_sort_atoms_by_mass(non_term_root) or root_nbrs)[0]

    second_excl: list[pmd.Atom] = [
        a for a in tree_adj[second_atom.idx] if a.idx != initial_atom.idx
    ]
    if not second_excl:
        raise ValueError(
            "Cannot build Z-matrix root triplet: second atom (local idx=%d) has "
            "no tree-neighbors other than root (local idx=%d) in a %d-atom "
            "molecule. Choose a different root atom."
            % (second_atom.idx, initial_atom.idx, n)
        )

    if n > 3:
        non_term_excl: list[pmd.Atom] = [
            a for a in second_excl if len(a.bond_partners) > 1
        ]
        third_atom: pmd.Atom = (
            _sort_atoms_by_mass(non_term_excl) or _sort_atoms_by_mass(second_excl)
        )[0]
    else:
        third_atom = _sort_atoms_by_mass(second_excl)[0]

    z_i.append(initial_atom.idx)
    z_j.append(ZM_SENTINEL)
    z_k.append(ZM_SENTINEL)
    z_l.append(ZM_SENTINEL)

    z_i.append(second_atom.idx)
    z_j.append(initial_atom.idx)
    z_k.append(ZM_SENTINEL)
    z_l.append(ZM_SENTINEL)

    z_i.append(third_atom.idx)
    z_j.append(second_atom.idx)
    z_k.append(initial_atom.idx)
    z_l.append(ZM_SENTINEL)

    if n == 3:
        return z_i, z_j, z_k, z_l

    # -- General case: n >= 4 (work-queue / deferred-retry) --------------
    selected: set[int] = {initial_atom.idx, second_atom.idx, third_atom.idx}

    work: deque[tuple[pmd.Atom, pmd.Atom]] = deque()
    for seed in (initial_atom, second_atom, third_atom):
        for nb in tree_adj[seed.idx]:
            if nb.idx not in selected:
                work.append((nb, seed))

    deferred: list[tuple[pmd.Atom, pmd.Atom]] = []
    placed_at_last_retry: int = 3

    while len(z_i) < n:
        if not work:
            if not deferred:
                missing: list[int] = [
                    idx for idx in acyclic_graph.nodes if idx not in selected
                ]
                raise ValueError(
                    "Z-matrix construction terminated at %d/%d atoms with empty "
                    "queues. Disconnected acyclic graph? Unplaced local atom "
                    "indices: %s." % (len(z_i), n, missing)
                )
            if len(z_i) == placed_at_last_retry:
                raise ValueError(
                    "Z-matrix stalled at %d/%d atoms. %d deferred pair(s) could "
                    "not be resolved: %s. Likely cause: degenerate topology or an "
                    "unsuitable root atom."
                    % (
                        len(z_i),
                        n,
                        len(deferred),
                        [(a0.idx, a1.idx) for a0, a1 in deferred],
                    )
                )
            placed_at_last_retry = len(z_i)
            work.extend(deferred)
            deferred.clear()

        next_work: deque[tuple[pmd.Atom, pmd.Atom]] = deque()

        while work:
            a0, a1 = work.popleft()
            if a0.idx in selected:
                continue

            placed: bool = False
            for a2 in tree_adj[a1.idx]:
                if a2.idx == a0.idx or a2.idx not in selected:
                    continue
                if len(a2.bond_partners) <= 1:
                    continue  # a2 must be non-terminal in the full graph
                for a3 in tree_adj[a2.idx]:
                    if a3.idx == a1.idx or a3.idx not in selected:
                        continue
                    z_i.append(a0.idx)
                    z_j.append(a1.idx)
                    z_k.append(a2.idx)
                    z_l.append(a3.idx)
                    selected.add(a0.idx)
                    placed = True
                    for nb in tree_adj[a0.idx]:
                        if nb.idx not in selected:
                            next_work.append((nb, a0))
                    break
                if placed:
                    break

            if not placed:
                deferred.append((a0, a1))

        work = next_work

    return z_i, z_j, z_k, z_l
