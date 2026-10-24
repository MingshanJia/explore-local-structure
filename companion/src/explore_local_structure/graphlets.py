"""Two- to four-node edge-typed graphlets with stable research coordinates.

Implements Jia et al., Encoding edge type information in graphlets (2022),
https://doi.org/10.1371/journal.pone.0273609. See docs/graphlets.md for the
distinction between edge histograms, type sets and type multisets.
"""

from itertools import combinations, combinations_with_replacement
from numbers import Integral

import networkx as nx

from ._common import _selected_nodes

__all__ = [
    "induced_graphlet_degree_vector",
    "induced_graphlet_degree_vector_ego",
    "graphlet_degree_vector_ego",
    "typed_edge_induced_graphlet_degree_vector",
    "typed_edge_induced_graphlet_degree_vector_ego",
    "typed_edge_graphlet_degree_vector_ego",
    "colored_graphlet_vector_for_typed_edge",
    "hetero_graphlet_vector_for_typed_edge",
    "colored_ego_graphlet_vector_for_typed_edge",
    "hetero_ego_graphlet_vector_for_typed_edge",
    "graphlet_feature_names",
    "get_type_from_index",
    "get_type_from_index_for_ego_network",
    "get_total_number_of_types_of_colored_graphlets_for_typed_edge",
    "get_total_number_of_types_of_hetero_graphlets_for_typed_edge",
    "get_total_number_of_types_of_colored_ego_graphlets_for_typed_edge",
    "get_total_number_of_types_of_hetero_ego_graphlets_for_typed_edge",
]

_EGO_ORBITS = (0, 2, 3, 7, 11, 13, 14)
_EDGE_COUNTS = (1, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 6)
# A connected graph on <=4 vertices is identified by its degree sequence;
# the focal degree then identifies its node orbit (standard 0..14 numbering).
_ORBITS = {
    (1, 1): {1: 0},
    (1, 1, 2): {1: 1, 2: 2},
    (2, 2, 2): {2: 3},
    (1, 1, 2, 2): {1: 4, 2: 5},
    (1, 1, 1, 3): {1: 6, 3: 7},
    (2, 2, 2, 2): {2: 8},
    (1, 2, 2, 3): {1: 9, 2: 10, 3: 11},
    (2, 2, 3, 3): {2: 12, 3: 13},
    (3, 3, 3, 3): {3: 14},
}


def _num_types(num_type):
    if isinstance(num_type, bool) or not isinstance(num_type, Integral) or num_type < 1:
        raise ValueError("num_type must be a positive integer")
    return int(num_type)


def _prepare(G, nodes, num_type=None, edge_type="edge_type"):
    if G.is_directed() or G.is_multigraph():
        raise nx.NetworkXNotImplemented("Graphlets require a simple undirected graph")
    selected, _ = _selected_nodes(G, nodes)
    adjacency = {u: set(G[u]) - {u} for u in G}
    types = {}
    if num_type is not None:
        num_type = _num_types(num_type)
        for u, v, data in G.edges(data=True):
            if u == v:
                continue
            value = data.get(edge_type)
            if isinstance(value, str):
                try:
                    value = int(value)
                except ValueError:
                    value = None
            if (
                isinstance(value, bool)
                or not isinstance(value, Integral)
                or not 0 <= value <= num_type
            ):
                raise ValueError(
                    f"Every non-loop edge needs {edge_type!r} in 0..{num_type} "
                    "(an integer or integer string)"
                )
            types[u, v] = types[v, u] = int(value)
    return selected, adjacency, types


def _orbit(vertices, edges):
    degree = dict.fromkeys(vertices, 0)
    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
    return _ORBITS[tuple(sorted(degree.values()))][degree[vertices[0]]]


def _connected_sets(adjacency, focal):
    """Enumerate each rooted connected vertex set once, with at most 4 nodes."""

    def extend(vertices, candidates, excluded):
        candidates = set(candidates)
        excluded = set(excluded)
        while candidates:
            v = candidates.pop()
            chosen = (*vertices, v)
            yield chosen
            if len(chosen) < 4:
                new = (candidates | adjacency[v]) - set(chosen) - excluded
                yield from extend(chosen, new, excluded)
            excluded.add(v)

    yield from extend((focal,), adjacency[focal], set())


def _occurrences(adjacency, focal, ego, induced=True):
    if not ego:
        for vertices in _connected_sets(adjacency, focal):
            edges = tuple(
                (u, v) for u, v in combinations(vertices, 2) if v in adjacency[u]
            )
            yield _orbit(vertices, edges), edges
        return
    for size in (1, 2, 3):
        for alters in combinations(adjacency[focal], size):
            spokes = tuple((focal, v) for v in alters)
            links = tuple(
                (u, v) for u, v in combinations(alters, 2) if v in adjacency[u]
            )
            # Non-induced graphlets retain each subset of present alter edges.
            subsets = (
                (links,)
                if induced
                else (
                    subset
                    for k in range(len(links) + 1)
                    for subset in combinations(links, k)
                )
            )
            for subset in subsets:
                edges = spokes + subset
                yield _orbit((focal, *alters), edges), edges


def graphlet_feature_names(method, num_type, *, ego=False):
    """Return flat vector coordinates as ``(orbit_id, type_tuple)`` pairs.

    ``method`` is 'typed', 'colored' or 'hetero'. Typed coordinates are row-major
    entries of the 15*T (or 7*T ego) edge-count matrix; types are 1..T.
    Colored tuples are sets sorted by cardinality then lexicographically;
    hetero tuples are lexicographic multisets. Ego colored/hetero coordinates
    include type 0 in the same blocks as the original research implementation.
    """
    num_type = _num_types(num_type)
    if method not in ("typed", "colored", "hetero"):
        raise ValueError("method must be 'typed', 'colored' or 'hetero'")
    orbits = _EGO_ORBITS if ego else range(15)
    features = []
    for orbit in orbits:
        if method == "typed":
            features.extend((orbit, (t,)) for t in range(1, num_type + 1))
            continue
        allow_zero = ego and orbit in (3, 11, 13, 14)
        types = range(0 if allow_zero else 1, num_type + 1)
        count = _EDGE_COUNTS[orbit]
        if method == "colored":
            signatures = (
                c
                for k in range(1, min(count, len(types)) + 1)
                for c in combinations(types, k)
            )
        else:
            signatures = combinations_with_replacement(types, count)
            if allow_zero:
                zero_limit = {3: 1, 11: 1, 13: 2, 14: 3}[orbit]
                signatures = (c for c in signatures if c.count(0) <= zero_limit)
        features.extend((orbit, signature) for signature in signatures)
    return features


def _vectors(G, nodes, num_type, edge_type, method, ego, induced=True):
    selected, adjacency, types = _prepare(G, nodes, num_type, edge_type)
    orbits = _EGO_ORBITS if ego else tuple(range(15))
    rows = {orbit: row for row, orbit in enumerate(orbits)}
    if method not in ("untyped", "typed"):
        coordinates = graphlet_feature_names(method, num_type, ego=ego)
        lookup = {coordinate: idx for idx, coordinate in enumerate(coordinates)}
        if not ego and any(value == 0 for value in types.values()):
            raise ValueError("Full colored/hetero vectors require types in 1..num_type")
        if ego and any(types[u, v] == 0 for u in selected for v in adjacency[u]):
            raise ValueError(
                "Ego colored/hetero vectors require positive focal-edge types"
            )
    result = {}
    for focal in selected:
        if method == "typed":
            vector = [[0] * num_type for _ in orbits]
        elif method == "untyped":
            vector = [0] * len(orbits)
        else:
            vector = [0] * len(coordinates)
        for orbit, edges in _occurrences(adjacency, focal, ego, induced):
            row = rows[orbit]
            if method == "untyped":
                vector[row] += 1
            elif method == "typed":
                for u, v in edges:
                    t = types[u, v]
                    if t:
                        vector[row][t - 1] += 1
            else:
                values = [types[u, v] for u, v in edges]
                signature = tuple(
                    sorted(set(values) if method == "colored" else values)
                )
                vector[lookup[orbit, signature]] += 1
        result[focal] = vector
    return result


def induced_graphlet_degree_vector(G, nodes=None):
    """Return dictionaries of 15 induced node-orbit counts (2..4 vertices)."""
    return _vectors(G, nodes, None, None, "untyped", False)


def induced_graphlet_degree_vector_ego(G, nodes=None):
    """Return seven induced ego counts in orbit order 0,2,3,7,11,13,14."""
    return _vectors(G, nodes, None, None, "untyped", True)


def graphlet_degree_vector_ego(G, nodes=None):
    """Return seven non-induced ego counts, counting retained edge subsets."""
    return _vectors(G, nodes, None, None, "untyped", True, False)


def typed_edge_induced_graphlet_degree_vector(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return 15*T TyE-GDV matrices: edges of each type in each induced orbit.

    Results are always node-keyed dictionaries of integer lists, including for
    one selected node. Types 1..num_type select columns; type 0 affects topology
    but is not accumulated. Integer strings from legacy CSV inputs are accepted.
    Missing/out-of-range attributes raise ValueError. Loops are ignored and
    directed/multigraph inputs rejected. Weights are ignored; input is unchanged.

    Enumerates rooted connected vertex sets up to size four, so dense graphs
    can be expensive. See docs/graphlets.md for complexity and vector semantics.
    """
    return _vectors(G, nodes, num_type, edge_type, "typed", False)


def typed_edge_induced_graphlet_degree_vector_ego(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return 7*T induced TyE-EGDV matrices; conventions match the full version.

    Rows correspond to orbits 0,2,3,7,11,13,14. Every selected occurrence has
    the focal node adjacent to all its other vertices. Type-0 edges are retained
    topologically and ignored in type histograms. No ego-network copy is needed.
    """
    return _vectors(G, nodes, num_type, edge_type, "typed", True)


def typed_edge_graphlet_degree_vector_ego(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return non-induced 7*T ego type counts, retaining each eligible edge subset."""
    return _vectors(G, nodes, num_type, edge_type, "typed", True, False)


def colored_graphlet_vector_for_typed_edge(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return flat induced orbit/type-set counts; full graph types must be 1..T."""
    return _vectors(G, nodes, num_type, edge_type, "colored", False)


def hetero_graphlet_vector_for_typed_edge(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return flat induced orbit/type-multiset counts; full types must be 1..T."""
    return _vectors(G, nodes, num_type, edge_type, "hetero", False)


def colored_ego_graphlet_vector_for_typed_edge(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return induced ego/type-set counts in legacy order; alter types may be 0."""
    return _vectors(G, nodes, num_type, edge_type, "colored", True)


def hetero_ego_graphlet_vector_for_typed_edge(
    G, num_type, nodes=None, *, edge_type="edge_type"
):
    """Return induced ego/type-multiset counts; positive focal-edge types required."""
    return _vectors(G, nodes, num_type, edge_type, "hetero", True)


def _indices(method, num_type, idx_list, ego):
    names = graphlet_feature_names(method, num_type, ego=ego)
    indices = list(idx_list)
    if any(
        isinstance(i, bool) or not isinstance(i, Integral) or not 0 <= i < len(names)
        for i in indices
    ):
        raise ValueError("Every feature index must be an integer within the vector")
    return [names[i] for i in indices]


def get_type_from_index(method, num_type, idx_list):
    """Decode full colored/hetero flat indices using the original coordinate order."""
    if method not in ("colored", "hetero"):
        raise ValueError("method must be 'colored' or 'hetero'")
    return _indices(method, num_type, idx_list, False)


def get_type_from_index_for_ego_network(method, num_type, idx_list):
    """Decode ego indices; method is 'colored_ego' or 'hetero_ego'."""
    if method not in ("colored_ego", "hetero_ego"):
        raise ValueError("method must be 'colored_ego' or 'hetero_ego'")
    return _indices(method.removesuffix("_ego"), num_type, idx_list, True)


def get_total_number_of_types_of_colored_graphlets_for_typed_edge(num_type):
    """Return the full colored vector length, including all orbit blocks."""
    return len(graphlet_feature_names("colored", num_type))


def get_total_number_of_types_of_hetero_graphlets_for_typed_edge(num_type):
    """Return the full heterogeneous vector length."""
    return len(graphlet_feature_names("hetero", num_type))


def get_total_number_of_types_of_colored_ego_graphlets_for_typed_edge(num_type):
    """Return the colored ego vector length, including legacy zero-type blocks."""
    return len(graphlet_feature_names("colored", num_type, ego=True))


def get_total_number_of_types_of_hetero_ego_graphlets_for_typed_edge(num_type):
    """Return the heterogeneous ego vector length with legacy zero-type limits."""
    return len(graphlet_feature_names("hetero", num_type, ego=True))
