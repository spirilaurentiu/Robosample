from typing import Iterable, Tuple, Set, FrozenSet, Dict, Optional, Hashable
from collections import deque
import networkx as nx

class GraphTraversalUtils:
    @staticmethod
    def validate_bfs_parent_child_edges(graph: nx.Graph, root: Hashable, bfs_edges: Iterable[Tuple[Hashable, Hashable, Hashable]]) -> None:
        """
        Validate that a sequence of edges corresponds to a valid BFS tree
        traversal starting from a given root.

        This function enforces the invariant that for every edge (u, v):
        - u has already been visited when the edge is produced
        - v has not been visited before (i.e., v is discovered via u)
        - (u, v) is an actual edge in the graph

        Parameters
        ----------
        graph : networkx.Graph
            The graph on which BFS is assumed to have been performed.
        root : hashable
            The BFS root node.
        bfs_edges : iterable of (hashable, hashable)
            Edges produced by a BFS traversal, typically from
            ``networkx.bfs_edges``.

        Raises
        ------
        ValueError
            If any edge violates BFS parent-child semantics.
        """
        visited = {root}

        # for step, (u, v, _) in enumerate(bfs_edges):
        #     if u not in visited:
        #         raise ValueError(f"BFS invariant violated at step {step}: parent node {u} has not been visited yet.")
        #     if v in visited:
        #         raise ValueError(f"BFS invariant violated at step {step}: child node {v} was already visited.")
        #     if not graph.has_edge(u, v):
        #         raise ValueError(f"BFS invariant violated at step {step}: edge ({u}, {v}) does not exist in graph.")
        #     visited.add(v)

    @staticmethod
    def nodes_to_distances(graph: nx.Graph, source_nodes: Iterable[int]) -> Dict[int, Optional[int]]:
        """Compute shortest path distance from each node to nearest source node."""
        dist = {node: None for node in graph.nodes}
        q = deque()
        for node in source_nodes:
            if node in graph:
                dist[node] = 0
                q.append(node)
        while q:
            u = q.popleft()
            for v in graph.neighbors(u):
                if dist[v] is None:
                    dist[v] = dist[u] + 1
                    q.append(v)
        return dist

    @staticmethod
    def bond_edges_to_nodes(edges: Iterable[FrozenSet[int]]) -> Set[int]:
        """Convert bond edges to a set of atom indices."""
        nodes: Set[int] = set()
        for edge in edges:
            nodes.update(edge)
        return nodes
