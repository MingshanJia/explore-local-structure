"""Inner- and outer-node quadrangle formation, including chorded cycles."""

from networkx.utils import not_implemented_for

from ._common import _mean, _normalized_adjacency, _result, _selected_nodes, _sum

__all__ = [
    "inner_quadrangle_coefficient",
    "outer_quadrangle_coefficient",
    "quadrangle_coefficients",
    "average_inner_quadrangle_coefficient",
    "average_outer_quadrangle_coefficient",
    "global_quadrangle_coefficient",
    "number_of_quadrangles",
]


def _quadrangle_counts(G, selected, weight):
    adjacency = _normalized_adjacency(G, weight)
    neighbors = {v: set(nbrs) for v, nbrs in adjacency.items()}
    for i in selected:
        ni = neighbors[i]
        closed_terms, inner_terms, outer_terms = [], [], []
        for j in ni:
            for k in neighbors[j] - {i}:
                nk = neighbors[k]
                if weight is None:
                    closed_terms.append(len((ni & nk) - {j}))
                    inner_terms.append(len(ni - {j, k}))
                    outer_terms.append(len(nk - {i, j}))
                else:
                    path = adjacency[i][j] * adjacency[j][k]
                    closed_terms.append(
                        path
                        * _sum(
                            adjacency[i][ell] * adjacency[k][ell]
                            for ell in (ni & nk) - {j}
                        )
                    )
                    inner_terms.append(
                        path * _sum(adjacency[i][ell] for ell in ni - {j, k})
                    )
                    outer_terms.append(
                        path * _sum(adjacency[k][ell] for ell in nk - {i, j})
                    )
        if weight is None:
            yield i, sum(closed_terms), sum(inner_terms), sum(outer_terms)
        else:
            yield (
                i,
                _sum(closed_terms),
                _sum(inner_terms),
                _sum(outer_terms),
            )


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def inner_quadrangle_coefficient(G, nodes=None, weight=None):
    r"""Return the inner quadrangle (i-quad) coefficient.

    It is :math:`2Q_i/OQI_i`, the fraction of distinct-vertex three-edge
    paths with the focal node in an inner position that close to four-cycles.
    Chorded cycles are included.

    Parameters
    ----------
    G : Graph
        A simple undirected graph.
    nodes : node or iterable, optional
        Nodes to evaluate; a single existing node returns a float.
    weight : string or None, optional
        Finite nonnegative edge weights, normalized by the maximum non-loop
        weight in the full graph. Missing attributes are one. A closed cycle
        contributes the product of four normalized weights; an inner path
        contributes the product of its three weights.

    Returns
    -------
    float or dict
        Local i-quad coefficients; zero if the denominator is zero.

    Notes
    -----
    Self-loops are ignored; weights are not mutated. For simple, unweighted,
    loopless undirected graphs, this is exactly ``2*s/(1+s)`` where
    ``s = nx.square_clustering(G, i)`` in NetworkX 3.7. This equivalent
    transformation does not cover the weighted definition. I-quad is exposed
    alongside o-quad for a consistent path-closure interpretation and API.
    Work scans two-edge paths and neighbor-set intersections; weighted counting
    additionally sums weights over eligible third edges. Preparation uses
    O(n + m) space.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.cycle_graph(4)
    >>> G.add_edge(1, 4)
    >>> inner_quadrangle_coefficient(G, 0)
    0.6666666666666666
    >>> nx.square_clustering(G, 0)
    0.5

    References
    ----------
    .. [1] Jia, Gabrys and Musial, Measuring Quadrangle Formation in Complex
       Networks, IEEE TNSE 9(2), 538-551, 2022.
       https://doi.org/10.1109/TNSE.2021.3123735
       Accessible manuscript: https://arxiv.org/abs/2011.10763
    """
    selected, single = _selected_nodes(G, nodes)
    values = {
        v: float(q / inner) if inner else 0.0
        for v, q, inner, _ in _quadrangle_counts(G, selected, weight)
    }
    return _result(values, selected, single)


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def outer_quadrangle_coefficient(G, nodes=None, weight=None):
    r"""Return the outer quadrangle (o-quad) coefficient.

    It is :math:`2Q_i/OQO_i`, the fraction of distinct-vertex three-edge paths
    starting at the focal node that close to four-cycles, including chords.

    Parameters
    ----------
    G : Graph
        A simple undirected graph.
    nodes : node or iterable, optional
        Nodes to evaluate; None evaluates all nodes.
    weight : string or None, optional
        Nonnegative finite edge-weight attribute. Weights are normalized by
        the full graph's maximum non-loop weight; missing values are one.
        The numerator sums four-edge products, the denominator three-edge
        products along paths starting at the focal node.

    Returns
    -------
    float or dict
        Coefficients, zero when there are no weighted/unweighted opportunities.

    Notes
    -----
    Self-loops are ignored and input weights are not changed. Unlike i-quad's
    unweighted transformation from square clustering, o-quad requires the
    outer-node three-edge-path denominator. It can differ from i-quad even
    on a bipartite graph. Complexity is the same as the inner coefficient.

    See Also
    --------
    inner_quadrangle_coefficient

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.cycle_graph(4)
    >>> G.add_edge(1, 4)
    >>> outer_quadrangle_coefficient(G, 0)
    1.0

    References
    ----------
    .. [1] Jia, Gabrys and Musial, Measuring Quadrangle Formation in Complex
       Networks, IEEE TNSE 9(2), 538-551, 2022.
       https://doi.org/10.1109/TNSE.2021.3123735
    """
    selected, single = _selected_nodes(G, nodes)
    values = {
        v: float(q / outer) if outer else 0.0
        for v, q, _, outer in _quadrangle_counts(G, selected, weight)
    }
    return _result(values, selected, single)


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def quadrangle_coefficients(G, nodes=None, weight=None):
    """Compute (i-quad, o-quad) pairs using a shared counting pass.

    Parameters follow :func:`inner_quadrangle_coefficient`. Returns a pair for
    a single node, otherwise a dictionary of node-to-pair entries.
    """
    selected, single = _selected_nodes(G, nodes)
    values = {
        v: (float(q / inner) if inner else 0.0, float(q / outer) if outer else 0.0)
        for v, q, inner, outer in _quadrangle_counts(G, selected, weight)
    }
    return _result(values, selected, single)


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def average_inner_quadrangle_coefficient(G, nodes=None, weight=None, count_zeros=True):
    """Return the arithmetic mean of local i-quad coefficients.

    Parameters
    ----------
    G : Graph
        A simple undirected graph.
    nodes : node or iterable, optional
        Nodes to average; the full graph determines each coefficient.
    weight : string or None, optional
        Nonnegative edge-weight attribute, globally normalized as in the local
        coefficient. Missing weights default to one.
    count_zeros : bool, optional
        Include zero local coefficients, default True.

    Returns
    -------
    float
        Arithmetic mean, or zero when no values are retained.

    See Also
    --------
    inner_quadrangle_coefficient
    average_outer_quadrangle_coefficient
    """
    return _mean(inner_quadrangle_coefficient(G, nodes, weight), count_zeros)


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def average_outer_quadrangle_coefficient(G, nodes=None, weight=None, count_zeros=True):
    """Return the arithmetic mean of local o-quad coefficients.

    Parameters
    ----------
    G : Graph
        A simple undirected graph.
    nodes : node or iterable, optional
        Nodes to average; the full graph determines each coefficient.
    weight : string or None, optional
        Nonnegative edge-weight attribute; see the local coefficient.
    count_zeros : bool, optional
        Include zero local coefficients, default True.

    Returns
    -------
    float
        Arithmetic mean, or zero when no values are retained.

    See Also
    --------
    outer_quadrangle_coefficient
    average_inner_quadrangle_coefficient
    """
    return _mean(outer_quadrangle_coefficient(G, nodes, weight), count_zeros)


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def global_quadrangle_coefficient(G, weight=None):
    """Ratio of summed closed quadrangle counts to summed path opportunities.

    This is distinct from averaging local values. Inner and outer global
    denominators agree by reindexing paths. For unweighted bipartite graphs
    it equals ``nx.bipartite.robins_alexander_clustering``, including an explicit
    zero convention when no three-edge paths exist. ``weight`` uses the same
    normalization as the local functions. The weighted ratio-of-sums is a
    derived extension implemented by the legacy code, not a separate weighted
    global definition stated in the quadrangle paper.
    """
    counts = list(_quadrangle_counts(G, list(G), weight))
    numerator = _sum(row[1] for row in counts)
    denominator = _sum(row[2] for row in counts)
    return float(numerator / denominator) if denominator else 0.0


@not_implemented_for("directed")
@not_implemented_for("multigraph")
def number_of_quadrangles(G):
    """Count undirected four-cycles, including cycles with chordal edges."""
    return sum(q for _, q, _, _ in _quadrangle_counts(G, list(G), None)) // 8
