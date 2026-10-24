"""Triangle closure measured from the endpoint of a two-edge path."""

import networkx as nx
from networkx.utils import not_implemented_for

from ._common import _mean, _normalized_adjacency, _result, _selected_nodes, _sum

__all__ = ["closure", "average_closure", "directed_closure", "average_directed_closure"]


def _total_closure_counts(G, selected, weight):
    adjacency = _normalized_adjacency(G, weight)
    if G.is_directed():
        symmetric = {v: {} for v in G}
        for u, nbrs in adjacency.items():
            for v, value in nbrs.items():
                symmetric[u][v] = symmetric[u].get(v, 0) + value
                symmetric[v][u] = symmetric[v].get(u, 0) + value
        adjacency = symmetric
    factor = 2 if G.is_directed() else 1
    for i in selected:
        nbrs = adjacency[i]
        opportunities = _sum(
            wij * _sum(wjk for k, wjk in adjacency[j].items() if k != i)
            for j, wij in nbrs.items()
        )
        closed = _sum(
            wij * _sum(wjk * nbrs[k] for k, wjk in adjacency[j].items() if k in nbrs)
            for j, wij in nbrs.items()
        )
        yield i, closed, factor * opportunities


@not_implemented_for("multigraph")
def closure(G, nodes=None, weight=None):
    r"""Return local endpoint-based triangle closure coefficients.

    For an undirected graph, the coefficient of node :math:`i` is
    :math:`2T_i / \sum_{j\in N(i)}(d_j-1)`. For a directed graph,
    write :math:`M=A+A^T`; it is
    :math:`(M^3)_{ii}/(2\sum_j M_{ij}\sum_{k\ne i}M_{jk})`.
    Reciprocal arcs count separately. These opportunities start at the focal
    node, unlike the centered opportunities used by ``clustering``.

    Parameters
    ----------
    G : graph
        A simple undirected or directed NetworkX graph.
    nodes : node or iterable of nodes, optional
        Nodes to evaluate. A single existing node returns a float; an iterable
        or None returns a dictionary. A subset only selects output nodes.
    weight : string or None, optional
        Edge attribute for finite nonnegative weights, normalized by the
        maximum non-loop edge weight in the full graph. Missing values are one.
        Products of triangle weights are divided by products of path weights.

    Returns
    -------
    float or dict
        Coefficients of the requested nodes; zero when opportunities vanish.

    Notes
    -----
    Self-loops are ignored. Empty graphs return an empty dictionary. All-zero
    weighted graphs give zero. Input weights are never modified. Signed weights
    are outside the supported convention. The directed coefficient reduces to
    its undirected counterpart on reciprocal graphs with equal reverse weights.
    Work is proportional to the neighbor-of-neighbor scans of requested nodes;
    graph preparation takes O(n + m) space and time.

    Examples
    --------
    >>> G = nx.DiGraph([(0, 1), (0, 2), (1, 2), (1, 3)])
    >>> closure(G, 0)
    0.3333333333333333
    >>> nx.clustering(G, 0)
    0.5

    References
    ----------
    .. [1] Yin et al., The Local Closure Coefficient: A New Perspective On
       Network Clustering, WSDM 2019. https://doi.org/10.1145/3289600.3290991
    .. [2] Jia, Gabrys and Musial, Directed closure coefficient and its patterns,
       PLOS ONE 16(6), e0253822, 2021.
       https://doi.org/10.1371/journal.pone.0253822
    """
    selected, single = _selected_nodes(G, nodes)
    values = {
        v: float(numerator / denominator) if denominator else 0.0
        for v, numerator, denominator in _total_closure_counts(G, selected, weight)
    }
    return _result(values, selected, single)


@not_implemented_for("multigraph")
def average_closure(G, nodes=None, weight=None, count_zeros=True):
    """Return the arithmetic mean of local closure coefficients.

    Parameters
    ----------
    G : graph
        A simple undirected or directed graph.
    nodes : node or iterable, optional
        Nodes to average; the full graph determines each coefficient.
    weight : string or None, optional
        Weight attribute; see :func:`closure`.
    count_zeros : bool, optional
        Include zero coefficients in the mean (default True).

    Returns
    -------
    float
        Mean coefficient, or zero for an empty selected/retained set.

    See Also
    --------
    closure

    Examples
    --------
    >>> average_closure(nx.complete_graph(3))
    1.0
    """
    return _mean(closure(G, nodes, weight), count_zeros)


def _directed_pattern_counts(G, selected):
    successors = {v: set(G.successors(v)) - {v} for v in G}
    predecessors = {v: set(G.predecessors(v)) - {v} for v in G}
    for i in selected:
        si, pi = successors[i], predecessors[i]
        head = sum(len(si & predecessors[j]) + len(si & successors[j]) for j in si)
        end = sum(len(pi & predecessors[j]) + len(pi & successors[j]) for j in pi)
        mid = 2 * sum(len(pi & predecessors[j]) for j in si)
        cyclic = 2 * sum(len(pi & successors[j]) for j in si)
        head_den = sum(
            len(successors[j]) + len(predecessors[j]) - 1 - (j in pi) for j in si
        )
        end_den = sum(
            len(successors[j]) + len(predecessors[j]) - 1 - (j in si) for j in pi
        )
        mid_den = sum(len(predecessors[j]) - 1 for j in si) + sum(
            len(successors[j]) - 1 for j in pi
        )
        cyc_den = sum(len(successors[j]) - (j in pi) for j in si) + sum(
            len(predecessors[j]) - (j in si) for j in pi
        )
        opportunities = head_den + end_den
        mixed = (mid + cyclic) // 2
        yield (
            i,
            {
                "source": (head + mixed, opportunities),
                "target": (end + mixed, opportunities),
                "head": (head, head_den),
                "mid": (mid, mid_den),
                "end": (end, end_den),
                "cyclic": (cyclic, cyc_den),
            },
        )


@not_implemented_for("undirected")
@not_implemented_for("multigraph")
def directed_closure(G, nodes=None, weight=None, *, pattern="total"):
    """Return directed closure for an overall, closing-edge, or role pattern.

    Parameters
    ----------
    G : DiGraph
        A simple directed graph; self-loops are ignored.
    nodes : node or iterable, optional
        Output nodes, following the scalar/dictionary convention of ``closure``.
    weight : string or None, optional
        Only the total pattern supports nonnegative weighted input.
    pattern : {'total', 'source', 'target', 'head', 'mid', 'end', 'cyclic'}
        Total uses :func:`closure`. Source/target distinguish the direction of
        the closing edge. Head/mid/end distinguish the focal node's role in a
        transitive triangle; cyclic selects directed cycles. Directed triangle
        orientations embedded in reciprocal arcs are counted separately.

    Returns
    -------
    float or dict
        Requested coefficients, with zero for undefined opportunities.

    Raises
    ------
    NetworkXNotImplemented
        For undirected/multigraph inputs or weighted non-total patterns.
    ValueError
        For an unknown pattern or unsupported weight values.

    Notes
    -----
    Patterns use their own endpoint-opportunity denominators. They are not
    conventional directed clustering or induced triad census values. For unweighted inputs at every
    node, total equals (source + target) / 2. Reversing all arcs swaps source with
    target and head with end, and preserves total, mid and cyclic. Pattern
    counting scans intersections of predecessor/successor sets along incident
    arcs; graph preparation takes O(n + m) space and time.

    Examples
    --------
    >>> G = nx.DiGraph([(0, 1), (0, 2), (1, 2), (1, 3)])
    >>> directed_closure(G, 0, pattern="source")
    0.6666666666666666
    >>> directed_closure(G, 0, pattern="target")
    0.0

    References
    ----------
    .. [1] Jia, Gabrys and Musial, Directed closure coefficient and its patterns,
       PLOS ONE 16(6), e0253822, 2021.
       https://doi.org/10.1371/journal.pone.0253822
    """
    valid = {"total", "source", "target", "head", "mid", "end", "cyclic"}
    if pattern not in valid:
        raise ValueError(f"pattern must be one of {sorted(valid)}")
    if pattern == "total":
        return closure(G, nodes, weight)
    if weight is not None:
        raise nx.NetworkXNotImplemented(
            "Weighted non-total closure patterns are not defined here."
        )
    selected, single = _selected_nodes(G, nodes)
    values = {}
    for v, counts in _directed_pattern_counts(G, selected):
        numerator, denominator = counts[pattern]
        values[v] = float(numerator / denominator) if denominator else 0.0
    return _result(values, selected, single)


@not_implemented_for("undirected")
@not_implemented_for("multigraph")
def average_directed_closure(
    G, nodes=None, weight=None, count_zeros=True, *, pattern="total"
):
    """Return the mean directed closure for a selected pattern.

    Parameters
    ----------
    G : DiGraph
        A simple directed graph.
    nodes : node or iterable, optional
        Nodes to average.
    weight : string or None, optional
        Supported only for the total pattern.
    count_zeros : bool, optional
        Include zero values, default True.
    pattern : string, optional
        A pattern accepted by :func:`directed_closure`, default 'total'.

    Returns
    -------
    float
        Mean coefficient, or zero if no values are retained.

    See Also
    --------
    directed_closure
    """
    return _mean(directed_closure(G, nodes, weight, pattern=pattern), count_zeros)
