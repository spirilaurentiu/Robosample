"""indexing.py

Mapping between the two atom-index spaces used in the package:

* **local**    -- ParmEd / AMBER prmtop atom index (``atom.idx``).
* **compound** -- index assigned by a breadth-first traversal of the spanning
  forest, rooted at the chosen root atom.  Compound index 0 is the root;
  every other atom receives the next index in BFS *discovery* order.

Why compound order encodes "closeness to root"
----------------------------------------------
BFS visits every atom at depth ``d`` before any atom at depth ``d + 1``.
Therefore the compound index is monotonic with tree depth: if
``compound(a) < compound(b)`` then ``a`` is at most as deep as ``b`` (and was
discovered first within its level).  Consequently, *orienting* a topology
element so its first atom is closest to the root and its last atom is farthest
reduces to a single comparison of compound indices -- no separate shortest-path
computation is required, and the result is fully deterministic.
"""

from __future__ import annotations

import networkx as nx


class CompoundIndex:
    """
    Bidirectional local <-> compound atom-index map built from a BFS traversal.

    Parameters
    ----------
    acyclic_graph : nx.Graph
        Spanning forest of the molecule (ring-closing bonds already removed).
        Nodes are *local* atom indices.
    root_local_index : int
        Local index of the BFS root.  Becomes compound index 0.

    Raises
    ------
    ValueError
        If the BFS from ``root_local_index`` does not reach every node (i.e.
        the graph is disconnected at the root).
    """

    compound_to_local: list[int]
    """``compound_to_local[c]`` is the local index of compound atom ``c``."""

    local_to_compound: dict[int, int]
    """Inverse map: local index -> compound index."""

    def __init__(self, acyclic_graph: nx.Graph, root_local_index: int) -> None:
        self.compound_to_local = [root_local_index] + [
            child
            for _parent, child in nx.bfs_edges(acyclic_graph, source=root_local_index)
        ]
        self.local_to_compound = {
            local: compound for compound, local in enumerate(self.compound_to_local)
        }

        if len(self.compound_to_local) != acyclic_graph.number_of_nodes():
            missing = [
                n for n in acyclic_graph.nodes if n not in self.local_to_compound
            ]
            raise ValueError(
                "BFS from root (local idx=%d) reached %d of %d atoms; the "
                "spanning forest is disconnected at the root. Unreached local "
                "atom indices: %s."
                % (
                    root_local_index,
                    len(self.compound_to_local),
                    acyclic_graph.number_of_nodes(),
                    missing,
                )
            )

    def to_compound(self, local_index: int) -> int:
        """Return the compound index of *local_index*."""
        return self.local_to_compound[local_index]

    def to_compound_tuple(self, local_indices: tuple[int, ...]) -> tuple[int, ...]:
        """Map a tuple of local indices to compound indices (order preserved)."""
        return tuple(self.local_to_compound[i] for i in local_indices)

    def orient_local(self, local_indices: tuple[int, ...]) -> tuple[int, ...]:
        """
        Return *local_indices* oriented so the first atom is closest to the root
        (smallest compound index) and the last atom is farthest.

        The element is treated as a path: it is either kept as-is or reversed
        in full, so any intermediate atom (e.g. the central atom of an angle)
        keeps its position.  Comparison uses compound indices, which are unique,
        so there is never a tie.
        """
        if (
            self.local_to_compound[local_indices[0]]
            > self.local_to_compound[local_indices[-1]]
        ):
            return tuple(reversed(local_indices))
        return tuple(local_indices)
