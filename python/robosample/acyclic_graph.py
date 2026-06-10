"""acyclic_graph.py

Construction of the molecule's spanning forest by selecting ring-closing
(cotree) bonds.

The ring-closing heuristic is the *maximum-weight spanning tree* approach: the
ring-closing bonds are chosen as the complement of a spanning tree whose edge
weights encode how desirable each bond is to KEEP as a normal tree edge.  This
never dead-ends on heavily fused polycyclic systems (e.g. cucurbiturils) the
way greedy cycle-basis breaking does, and topological bridges (e.g. an
inter-chain disulfide that is the only link between two chains) are guaranteed
to remain in the tree, so connectivity is preserved automatically.

This module deals exclusively in *local* (ParmEd) atom indices; BFS / compound
ordering is layered on top afterwards (see ``indexing.CompoundIndex``).
"""

from __future__ import annotations

import logging

import networkx as nx
import parmed as pmd

from .amber_dihedral_types import DihedralType
from .bond_util import is_rigid_bond

logger: logging.Logger = logging.getLogger(__name__)


def atom_label(molecule: pmd.Structure, idx: int) -> str:
    """Return a human-readable label for atom *idx* (used in log messages)."""
    a = molecule.atoms[idx]
    return "%s%d_%s_%d" % (a.residue.name, a.residue.idx + 1, a.name, idx + 1)


def find_forbidden_bonds(
    molecule: pmd.Structure, atom_type_pairs: list[tuple[str, str]]
) -> set[frozenset[int]]:
    """
    Identify bonds whose two endpoint atom *names* match any forbidden
    atom-type pair (case-insensitive, order-independent).

    Forbidden bonds are preferentially KEPT in the spanning tree rather than
    turned into ring closures -- e.g. disulfides (``("SG", "SG")``), which carry
    no useful torsional degree of freedom and may be the only covalent link
    between two chains.
    """
    forbidden_pairs: set[frozenset[str]] = {
        frozenset((a.lower(), b.lower())) for a, b in atom_type_pairs
    }
    forbidden_edges: set[frozenset[int]] = set()
    for bond in molecule.bonds:
        a1, a2 = bond.atom1, bond.atom2
        if frozenset((a1.name.lower(), a2.name.lower())) in forbidden_pairs:
            forbidden_edges.add(frozenset((a1.idx, a2.idx)))
    return forbidden_edges


def build_acyclic_graph(
    molecule: pmd.Structure,
    dihedral_classifier,
) -> tuple[nx.Graph, dict[frozenset[int], DihedralType]]:
    """
    Build the spanning forest and classify every bond's dihedral type.

    Pipeline
    --------
    1. Build the full molecular bond graph and classify every bond's dihedral
       type, populating the returned cache.  Bonds whose best classification is
       ``PROTEIN_RING_DIHEDRAL`` are recorded as *candidate* ring-closing bonds.
    2. First pass -- cut each candidate ring-closing bond, but only where safe:
       a bond whose removal would disconnect the molecule (a bridge) is kept and
       a warning logged; at most one closure is anchored per atom (a soft
       preference -- deferred bonds are reconsidered in step 3).
    3. Reduce residual cycles to a spanning forest as the complement of a
       maximum-weight spanning tree (``keep_weight``: forbidden > one-per-atom >
       flexible-over-rigid > near-forbidden-atom > higher-degree).
    4-5. Validate the ring-closing bond set and run final acyclicity /
       connectivity sanity checks.

    Returns
    -------
    (nx.Graph, dict[frozenset[int], DihedralType])
        The spanning forest (nodes are local atom indices) and the dihedral
        type cache keyed by the frozenset of each bond's two atom indices.

    Raises
    ------
    ValueError
        If the result is still cyclic or disconnected, or if ring-closing bond
        validation fails.
    """
    num_atoms: int = len(molecule.atoms)
    atom_by_idx: dict[int, pmd.Atom] = {a.idx: a for a in molecule.atoms}
    dihedral_type_cache: dict[frozenset[int], DihedralType] = {}

    # ----------------------------------------------------------------
    # Step 1: full graph + classify + cache
    # ----------------------------------------------------------------
    full_graph: nx.Graph = nx.Graph()
    full_graph.add_nodes_from(a.idx for a in molecule.atoms)

    ring_dihedral_bonds: set[frozenset[int]] = set()

    for bond in molecule.bonds:
        parent: pmd.Atom = bond.atom2
        child: pmd.Atom = bond.atom1
        key: frozenset[int] = frozenset((parent.idx, child.idx))
        full_graph.add_edge(parent.idx, child.idx)

        best_type: DihedralType = DihedralType.UNKNOWN
        is_ring_dihedral: bool = False

        for grandparent in parent.bond_partners:
            if grandparent is child:
                continue
            for gchild in child.bond_partners:
                if gchild is parent:
                    continue
                candidate = pmd.Dihedral(grandparent, parent, child, gchild)
                dtype: DihedralType = dihedral_classifier.classify(candidate)
                if dtype == DihedralType.PROTEIN_RING_DIHEDRAL:
                    is_ring_dihedral = True
                    best_type = dtype
                    break
                if dtype != DihedralType.UNKNOWN:
                    best_type = dtype
            if is_ring_dihedral:
                break

        dihedral_type_cache[key] = best_type
        if is_ring_dihedral:
            ring_dihedral_bonds.add(key)

    # ----------------------------------------------------------------
    # Step 2: cut classified ring dihedrals where safe
    # ----------------------------------------------------------------
    g: nx.Graph = full_graph.copy()
    forbidden_edges: set[frozenset[int]] = find_forbidden_bonds(
        molecule, [("SG", "SG")]
    )

    ring_closing_set: set[frozenset[int]] = set()
    closures_per_atom: dict[int, int] = {}

    def _can_cut(u: int, v: int) -> bool:
        """Remove edge (u, v) iff the graph stays connected (it lies on a cycle)."""
        g.remove_edge(u, v)
        if nx.has_path(g, u, v):
            return True
        g.add_edge(u, v)  # bridge: restore
        return False

    # Iterate bonds in ParmEd order (not over the set) for reproducibility.
    for bond in molecule.bonds:
        u: int = bond.atom2.idx
        v: int = bond.atom1.idx
        key = frozenset((u, v))
        if key not in ring_dihedral_bonds:
            continue
        if closures_per_atom.get(u, 0) >= 1 or closures_per_atom.get(v, 0) >= 1:
            continue  # soft one-per-atom: defer to spanning step
        if not _can_cut(u, v):
            logger.warning(
                "Bond %d -- %d is classified as ring-closing but is a bridge "
                "(its removal would disconnect the molecule); keeping it in the "
                "spanning tree to preserve connectivity.",
                u,
                v,
            )
            continue
        ring_closing_set.add(key)
        for idx in (u, v):
            closures_per_atom[idx] = closures_per_atom.get(idx, 0) + 1

    # ----------------------------------------------------------------
    # Step 3: reduce residual cycles via maximum-weight spanning tree
    # ----------------------------------------------------------------
    forbidden_nodes: set[int] = {n for edge in forbidden_edges for n in tuple(edge)}
    forbidden_dist: dict[int, int] = (
        nx.multi_source_shortest_path_length(g, forbidden_nodes)
        if forbidden_nodes
        else {}
    )

    def _keep_weight(u: int, v: int) -> float:
        """Soft preference for KEEPING bond (u, v) as a tree edge (higher == keep)."""
        weight: float = 0.0
        if frozenset((u, v)) in forbidden_edges:
            weight += 1e9
        if closures_per_atom.get(u, 0) or closures_per_atom.get(v, 0):
            weight += 1e6
        if not is_rigid_bond(atom_by_idx[u], atom_by_idx[v]):
            weight += 1e3
        if forbidden_dist:
            dist = max(forbidden_dist.get(u, 0), forbidden_dist.get(v, 0))
            weight += -0.1 * float(dist)
        weight += 1e-3 * (g.degree[u] + g.degree[v])
        return weight

    for u, v in g.edges:
        g[u][v]["keep_weight"] = _keep_weight(u, v)

    spanning: nx.Graph = nx.maximum_spanning_tree(g, weight="keep_weight")
    kept_edges: set[frozenset[int]] = {frozenset(e) for e in spanning.edges}

    for u, v in list(g.edges):
        if frozenset((u, v)) in kept_edges:
            continue
        g.remove_edge(u, v)
        ring_closing_set.add(frozenset((u, v)))
        for idx in (u, v):
            closures_per_atom[idx] = closures_per_atom.get(idx, 0) + 1

    overloaded: list[int] = sorted(
        idx for idx, count in closures_per_atom.items() if count > 1
    )
    if overloaded:
        shown = ", ".join(atom_label(molecule, idx) for idx in overloaded[:10])
        logger.warning(
            "%d atom(s) participate in more than one ring-closing bond. This is "
            "expected for fused polycyclic systems (e.g. cucurbituril) and is "
            "treated as a soft preference rather than an error. Atoms: %s%s.",
            len(overloaded),
            shown,
            ", ..." if len(overloaded) > 10 else "",
        )

    # ----------------------------------------------------------------
    # Step 4-5: validation
    # ----------------------------------------------------------------
    if full_graph.number_of_nodes() != num_atoms:
        raise ValueError(
            "Acyclic graph node count %d != atom count %d."
            % (full_graph.number_of_nodes(), num_atoms)
        )
    validate_ring_closing_bonds(full_graph, g, ring_closing_set)
    check_residual_cycles(molecule, g)
    check_disconnected_graph(molecule, g)

    return g, dihedral_type_cache


def validate_ring_closing_bonds(
    full_graph: nx.Graph,
    spanning_tree: nx.Graph,
    ring_closing_set: set[frozenset[int]],
) -> None:
    """
    Validate the ring-closing bond set: (1) its size equals the cyclomatic
    number, (2) each ring-closing bond's endpoints are connected in the spanning
    tree, and (3) no ring-closing bond lies inside another's fundamental cycle.
    """
    n_comp: int = nx.number_connected_components(spanning_tree)
    cyclomatic: int = (
        full_graph.number_of_edges() - full_graph.number_of_nodes() + n_comp
    )
    if len(ring_closing_set) != cyclomatic:
        raise ValueError(
            "Ring-closing bond count (%d) != cyclomatic number (%d). Either a "
            "cycle has no ring-closing bond, or a non-cyclic (bridge) bond was "
            "incorrectly removed." % (len(ring_closing_set), cyclomatic)
        )

    fundamental_cycle_edges: dict[frozenset[int], frozenset[frozenset[int]]] = {}
    for rc_bond in ring_closing_set:
        u, v = tuple(rc_bond)
        if not nx.has_path(spanning_tree, u, v):
            raise ValueError(
                "Ring-closing bond %d -- %d has no path in the spanning tree. "
                "The bond may be a bridge, or the spanning tree is disconnected "
                "at this node." % (u, v)
            )
        path: list[int] = nx.shortest_path(spanning_tree, u, v)
        tree_edges: frozenset[frozenset[int]] = frozenset(
            frozenset((path[i], path[i + 1])) for i in range(len(path) - 1)
        )
        fundamental_cycle_edges[rc_bond] = tree_edges | {rc_bond}

    for rc_bond, cycle_edges in fundamental_cycle_edges.items():
        for other_rc in ring_closing_set:
            if other_rc == rc_bond:
                continue
            if other_rc in cycle_edges:
                u1, v1 = tuple(rc_bond)
                u2, v2 = tuple(other_rc)
                raise ValueError(
                    "Ring-closing bond %d -- %d appears inside the fundamental "
                    "cycle of ring-closing bond %d -- %d. Ring-closing bonds must "
                    "be independent. This indicates a bug in the classifier or "
                    "bond-removal logic." % (u2, v2, u1, v1)
                )

    logger.debug(
        "Ring-closing bond validation passed: %d bond(s), cyclomatic number %d.",
        len(ring_closing_set),
        cyclomatic,
    )


def check_residual_cycles(molecule: pmd.Structure, g: nx.Graph) -> None:
    """Raise ``ValueError`` if *g* contains cycles, logging each one."""
    if nx.is_forest(g):
        return
    cycles = nx.cycle_basis(g)
    logger.error("Graph contains %d residual cycle(s):", len(cycles))
    for ci, cycle in enumerate(cycles, start=1):
        logger.error("  Cycle %d:", ci)
        for j in range(len(cycle)):
            n1 = atom_label(molecule, cycle[j])
            n2 = atom_label(molecule, cycle[(j + 1) % len(cycle)])
            logger.error("    %s -- %s", n1, n2)
    raise ValueError("Acyclic graph still contains cycles after ring-closure removal.")


def check_disconnected_graph(molecule: pmd.Structure, g: nx.Graph) -> None:
    """Raise ``ValueError`` if *g* is disconnected, logging each component."""
    if nx.is_connected(g):
        return
    components = list(nx.connected_components(g))
    logger.error(
        "Graph is disconnected (%d components) after ring-closure removal.",
        len(components),
    )
    for ci, comp in enumerate(components, start=1):
        logger.error("  Component %d (%d atoms):", ci, len(comp))
        for idx in comp:
            logger.error("    - %s", atom_label(molecule, idx))
    raise ValueError("Acyclic graph is disconnected after ring-closure removal.")
