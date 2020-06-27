"""Algorithms to characterize the number of triangles in a graph."""

from itertools import chain
from itertools import combinations
from collections import Counter
import networkx as nx

from networkx.utils import not_implemented_for

__all__ = ['triangles', 'average_clustering', 'clustering', 'transitivity',
           'square_clustering', 'quadrangle_coefficient', 'quadrangle_coefficient_2', 'generalized_degree', 'average_closure', 'closure',
           'src_closure', 'tgt_closure', 'head_closure', 'mid_closure', 'end_closure', 'cyc_closure']


@not_implemented_for('directed')
def triangles(G, nodes=None):
    """Compute the number of triangles.

    Finds the number of triangles that include a node as one vertex.

    Parameters
    ----------
    G : graph
       A networkx graph
    nodes : container of nodes, optional (default= all nodes in G)
       Compute triangles for nodes in this container.

    Returns
    -------
    out : dictionary
       Number of triangles keyed by node label.

    Examples
    --------
    >>> G=nx.complete_graph(5)
    >>> print(nx.triangles(G,0))
    6
    >>> print(nx.triangles(G))
    {0: 6, 1: 6, 2: 6, 3: 6, 4: 6}
    >>> print(list(nx.triangles(G,(0,1)).values()))
    [6, 6]

    Notes
    -----
    When computing triangles for the entire graph each triangle is counted
    three times, once at each node.  Self loops are ignored.

    """
    # If `nodes` represents a single node in the graph, return only its number
    # of triangles.
    if nodes in G:
        return next(_triangles_and_degree_iter(G, nodes))[2] // 2
    # Otherwise, `nodes` represents an iterable of nodes, so return a
    # dictionary mapping node to number of triangles.
    return {v: t // 2 for v, d, t, _ in _triangles_and_degree_iter(G, nodes)}


@not_implemented_for('multigraph')
def _triangles_and_degree_iter(G, nodes=None):
    """ Return an iterator of (node, degree, triangles, generalized degree).

    This double counts triangles so you may want to divide by 2.
    See degree(), triangles() and generalized_degree() for definitions
    and details.

    """
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for v, v_nbrs in nodes_nbrs:
        vs = set(v_nbrs) - {v}
        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        ntriangles = sum(k * val for k, val in gen_degree.items())
        yield (v, len(vs), ntriangles, gen_degree)

# for clo-co
@not_implemented_for('multigraph')
def _triangles_and_opentriads_iter(G, nodes=None):

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for v, v_nbrs in nodes_nbrs:
        vs = set(v_nbrs) - {v}
        ot = 0
        for w in vs:
            ns = set(G[w]) - {w}
            ot += len(ns) - 1

        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        ntriangles = sum(k * val for k, val in gen_degree.items())
        yield (v, ntriangles, ot)


# @not_implemented_for('multigraph')
# def _weighted_triangles_and_degree_iter(G, nodes=None, weight='weight'):
#     """ Return an iterator of (node, degree, weighted_triangles).
#
#     Used for weighted clustering.
#
#     """
#     if weight is None or G.number_of_edges() == 0:
#         max_weight = 1
#     else:
#         max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))
#     if nodes is None:
#         nodes_nbrs = G.adj.items()
#     else:
#         nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))
#
#     def wt(u, v):
#         return G[u][v].get(weight, 1) / max_weight
#
#     for i, nbrs in nodes_nbrs:
#         inbrs = set(nbrs) - {i}
#         weighted_triangles = 0
#         seen = set()
#         for j in inbrs:
#             seen.add(j)
#             # This prevents double counting.    Note: prevent double here, but double later... funny..
#             jnbrs = set(G[j]) - seen
#             # Only compute the edge weight once, before the inner inner
#             # loop.
#             wij = wt(i, j)
#             weighted_triangles += sum((wij * wt(j, k) * wt(k, i)) ** (1 / 3)
#                                       for k in inbrs & jnbrs)
#         yield (i, len(inbrs), 2 * weighted_triangles)


# another way to calculate weighted clustering-co; to replace _weighted_triangles_and_degree_iter
# paper:A general framework for weighted gene co-expression  network  analysis
@not_implemented_for('multigraph')
def _weighted_triangles_and_otc_iter(G, nodes=None, weight='weight'):
    """ Return an iterator of (node, weighted_triangles, weighted center-opentriads).

    Used for weighted clustering.

    """
    if weight is None or G.number_of_edges() == 0:
        max_weight = 1
    else:
        max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    def wt(u, v):
        return G[u][v].get(weight, 1) / max_weight

    for i, nbrs in nodes_nbrs:
        inbrs = set(nbrs) - {i}
        weighted_triangles = 0
        weighted_center_opentriads = 0
        seen = set()
        for j in inbrs:
            seen.add(j)
            # This prevents double counting.    Note: prevent double here, but double later... funny..
            jnbrs = set(G[j]) - seen
            # Only compute the edge weight once, before the inner inner
            # loop.
            wij = wt(i, j)
            weighted_triangles += sum((wij * wt(j, k) * wt(k, i))
                                      for k in inbrs & jnbrs)

            weighted_center_opentriads += sum(wij * wt(i, k) for k in (inbrs - {j}))

        yield (i, 2 * weighted_triangles, weighted_center_opentriads)




@not_implemented_for('multigraph')
def _directed_triangles_and_degree_iter(G, nodes=None):
    """ Return an iterator of
    (node, total_degree, reciprocal_degree, directed_triangles).

    Used for directed clustering.

    """
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        directed_triangles = 0
        for j in chain(ipreds, isuccs):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            directed_triangles += sum(1 for k in
                                       chain((ipreds & jpreds),
                                             (ipreds & jsuccs),
                                             (isuccs & jpreds),
                                             (isuccs & jsuccs)))
        dtotal = len(ipreds) + len(isuccs)
        dbidirectional = len(ipreds & isuccs)
        yield (i, dtotal, dbidirectional, directed_triangles)

# for closure-co
def _directed_triangles_and_opentriads_iter(G, nodes=None):
    """ Return an iterator of
    (node, directed_triangles, open_triads).

    Used for directed clustering and closure.

    """
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        ts = 0
        tt = 0
        directed_triangles = 0
        for j in chain(ipreds, isuccs):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            ts += sum(1 for k in
                      chain((isuccs & jpreds),
                            (isuccs & jsuccs)))

            tt += sum(1 for k in
                      chain((ipreds & jpreds),
                            (ipreds & jsuccs)))

        directed_triangles += (ts + tt)

        open_triads = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 4 * (dj - 2)
        for j in ((ipreds | isuccs) - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 2 * (dj - 1)

        yield (i, ts, tt,  open_triads)


#for src-clo
def _directed_src_triangles_and_opentriads_iter(G, nodes=None):

    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        ts = 0
        for j in chain(ipreds, isuccs):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            ts += sum(1 for k in
                      chain((isuccs & jpreds),
                            (isuccs & jsuccs)))

        open_triads = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 4 * (dj - 2)
        for j in ((ipreds | isuccs) - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 2 * (dj - 1)

        yield (i, ts, open_triads)

# for tgt-clo
def _directed_tgt_triangles_and_opentriads_iter(G, nodes=None):

    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        tt = 0
        for j in chain(ipreds, isuccs):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            tt += sum(1 for k in
                      chain((ipreds & jpreds),
                            (ipreds & jsuccs)))

        open_triads = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 4 * (dj - 2)
        for j in ((ipreds | isuccs) - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            open_triads += 2 * (dj - 1)

        yield (i, tt, open_triads)


# for head-of-path pattern
def _head_pattern_iter(G, nodes=None):
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        t_head = 0
        for j in isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            t_head += sum(1 for k in
                      chain((isuccs & jpreds),
                            (isuccs & jsuccs)))

        ot_head = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            ot_head += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ot_head += dj - 1
        yield (i, t_head, ot_head)


# for mid-of-path pattern
def _mid_pattern_iter(G, nodes=None):
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        t_mid = 0
        for j in isuccs:
            jpreds = set(G._pred[j]) - {j}

            t_mid += sum(1 for k in ipreds & jpreds)
        for j in ipreds:
            jsuccs = set(G._succ[j]) - {j}

            t_mid += sum(1 for k in isuccs & jsuccs)

        ot_mid = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ot_mid += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj_in = len(jpreds)

            ot_mid += dj_in - 1

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj_out = len(jsuccs)

            ot_mid += dj_out - 1

        yield (i, t_mid, ot_mid)


# for end-of-path pattern
def _end_pattern_iter(G, nodes=None):
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        t_end = 0
        for j in ipreds:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            t_end += sum(1 for k in
                         chain((ipreds & jpreds),
                               (ipreds & jsuccs)))

        ot_end = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            ot_end += dj - 2

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ot_end += dj - 1

        yield (i, t_end, ot_end)


#for cyclic pattern
def _cyc_pattern_iter(G, nodes=None):
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        t_cyc = 0
        for j in isuccs:
            jsuccs = set(G._succ[j]) - {j}

            t_cyc += sum(1 for k in ipreds & jsuccs)

        for j in ipreds:
            jpreds = set(G._pred[j]) - {j}

            t_cyc += sum(1 for k in isuccs & jpreds)

        ot_cyc = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ot_cyc += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jsuccs = set(G._succ[j]) - {j}
            dj_out = len(jsuccs)

            ot_cyc += dj_out

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            dj_in = len(jpreds)

            ot_cyc += dj_in

        yield (i, t_cyc, ot_cyc)


# @not_implemented_for('multigraph')
# def _directed_weighted_triangles_and_degree_iter(G, nodes=None, weight='weight'):
#     """ Return an iterator of
#     (node, total_degree, reciprocal_degree, directed_weighted_triangles).
#
#     Used for directed weighted clustering.
#
#     """
#     if weight is None or G.number_of_edges() == 0:
#         max_weight = 1
#     else:
#         max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))
#
#     nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))
#
#     def wt(u, v):
#         return G[u][v].get(weight, 1) / max_weight
#
#     for i, preds, succs in nodes_nbrs:
#         ipreds = set(preds) - {i}
#         isuccs = set(succs) - {i}
#
#         directed_triangles = 0
#         for j in ipreds:
#             jpreds = set(G._pred[j]) - {j}
#             jsuccs = set(G._succ[j]) - {j}
#             directed_triangles += sum((wt(j, i) * wt(k, i) * wt(k, j))**(1 / 3)
#                                       for k in ipreds & jpreds)
#             directed_triangles += sum((wt(j, i) * wt(k, i) * wt(j, k))**(1 / 3)
#                                       for k in ipreds & jsuccs)
#             directed_triangles += sum((wt(j, i) * wt(i, k) * wt(k, j))**(1 / 3)
#                                       for k in isuccs & jpreds)
#             directed_triangles += sum((wt(j, i) * wt(i, k) * wt(j, k))**(1 / 3)
#                                       for k in isuccs & jsuccs)
#
#         for j in isuccs:
#             jpreds = set(G._pred[j]) - {j}
#             jsuccs = set(G._succ[j]) - {j}
#             directed_triangles += sum((wt(i, j) * wt(k, i) * wt(k, j))**(1 / 3)
#                                       for k in ipreds & jpreds)
#             directed_triangles += sum((wt(i, j) * wt(k, i) * wt(j, k))**(1 / 3)
#                                       for k in ipreds & jsuccs)
#             directed_triangles += sum((wt(i, j) * wt(i, k) * wt(k, j))**(1 / 3)
#                                       for k in isuccs & jpreds)
#             directed_triangles += sum((wt(i, j) * wt(i, k) * wt(j, k))**(1 / 3)
#                                       for k in isuccs & jsuccs)
#
#         dtotal = len(ipreds) + len(isuccs)
#         dbidirectional = len(ipreds & isuccs)
#         yield (i, dtotal, dbidirectional, directed_triangles)


# another way to calculate weighted clustering co; replace _directed_weighted_triangles_and_degree_iter
# paper:A general framework for weighted gene co-expression  network  analysis
@not_implemented_for('multigraph')
def _directed_weighted_triangles_and_otc_iter(G, nodes=None, weight='weight'):

    if weight is None or G.number_of_edges() == 0:
        max_weight = 1
    else:
        max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))

    def wt(u, v):
        return G[u][v].get(weight, 1) / max_weight

    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        directed_triangles = 0
        directed_center_opentriads = 0
        for j in ipreds:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            directed_triangles += sum((wt(j, i) * wt(k, i) * wt(k, j))
                                      for k in ipreds & jpreds)
            directed_triangles += sum((wt(j, i) * wt(k, i) * wt(j, k))
                                      for k in ipreds & jsuccs)
            directed_triangles += sum((wt(j, i) * wt(i, k) * wt(k, j))
                                      for k in isuccs & jpreds)
            directed_triangles += sum((wt(j, i) * wt(i, k) * wt(j, k))
                                      for k in isuccs & jsuccs)

            directed_center_opentriads += sum(abs(wt(j, i) * wt(k, i)) for k in (ipreds - {j}))
            directed_center_opentriads += sum(abs(wt(j, i) * wt(i, k)) for k in (isuccs - {j}))


        for j in isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            directed_triangles += sum((wt(i, j) * wt(k, i) * wt(k, j))
                                      for k in ipreds & jpreds)
            directed_triangles += sum((wt(i, j) * wt(k, i) * wt(j, k))
                                      for k in ipreds & jsuccs)
            directed_triangles += sum((wt(i, j) * wt(i, k) * wt(k, j))
                                      for k in isuccs & jpreds)
            directed_triangles += sum((wt(i, j) * wt(i, k) * wt(j, k))
                                      for k in isuccs & jsuccs)

            directed_center_opentriads += sum(abs(wt(i, j) * wt(k, i)) for k in (ipreds - {j}))
            directed_center_opentriads += sum(abs(wt(i, j) * wt(i, k)) for k in (isuccs - {j}))

        yield (i, directed_triangles, directed_center_opentriads)



# for clo-co
@not_implemented_for('multigraph')
def _directed_weighted_triangles_and_opentriads_iter(G, nodes=None, weight='weight'):
    """ Return an iterator of
    (node, total_degree, reciprocal_degree, directed_weighted_triangles).

    Used for directed weighted clustering.

    """
    if weight is None or G.number_of_edges() == 0:
        max_weight = 1
    else:
        max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))

    # normalize weight in [0,1]
    for u, v, data in G.edges(data=True):
        data[weight] /= max_weight

    # for signed weighted, get a copy of G with weight in absolute value
    G_abs = nx.DiGraph()
    e_list = list(G.edges(data=True))
    G_abs.add_edges_from(e_list)
    for u, v, data in G_abs.edges(data=True):
        data[weight] = abs(data[weight])

    def wt(u, v):
        return G[u][v].get(weight, 1)

    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}

        directed_triangles = 0
        for j in ipreds:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            directed_triangles += sum((wt(j, i) * wt(k, i) * wt(k, j))
                                      for k in ipreds & jpreds)
            directed_triangles += sum((wt(j, i) * wt(k, i) * wt(j, k))
                                      for k in ipreds & jsuccs)
            directed_triangles += sum((wt(j, i) * wt(i, k) * wt(k, j))
                                      for k in isuccs & jpreds)
            directed_triangles += sum((wt(j, i) * wt(i, k) * wt(j, k))
                                      for k in isuccs & jsuccs)

        for j in isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            directed_triangles += sum((wt(i, j) * wt(k, i) * wt(k, j))
                                      for k in ipreds & jpreds)
            directed_triangles += sum((wt(i, j) * wt(k, i) * wt(j, k))
                                      for k in ipreds & jsuccs)
            directed_triangles += sum((wt(i, j) * wt(i, k) * wt(k, j))
                                      for k in isuccs & jpreds)
            directed_triangles += sum((wt(i, j) * wt(i, k) * wt(j, k))
                                      for k in isuccs & jsuccs)
        # three conditions
        ot = 0
        for j in ipreds & isuccs:
            ot += 2 * (abs(wt(i, j)) + abs(wt(j, i))) * (G_abs.degree(j, weight) - (abs(wt(i, j)) + abs(wt(j, i))))
        for j in isuccs - (ipreds & isuccs):
            ot += 2 * abs(wt(i, j)) * (G_abs.degree(j, weight) - abs(wt(i, j)))
        for j in ipreds - (ipreds & isuccs):
            ot += 2 * abs(wt(j, i)) * (G_abs.degree(j, weight) - abs(wt(j, i)))

        yield (i, directed_triangles, ot)



def average_clustering(G, nodes=None, weight=None, count_zeros=True):

    c = clustering(G, nodes, weight=weight).values()
    if not count_zeros:
        c = [v for v in c if v > 0]
    return sum(c) / len(c)

# for average closure-co
def average_closure(G, nodes=None, weight=None, count_zeros=True):

    ce = closure(G, nodes, weight=weight)
    list_ce = []
    for k, v in ce.items():
        list_ce.append(v[0])

    if not count_zeros:
        list_ce = [v for v in list_ce if v > 0]
    return sum(list_ce) / len(list_ce)


def clustering(G, nodes=None, weight=None):

    if G.is_directed():
        # change to another way
        if weight is not None:
            td_iter = _directed_weighted_triangles_and_otc_iter(
                G, nodes, weight)
            clusterc = {v: 0 if t == 0 else t / (2 * otc)
                        for v, t, otc in td_iter}
        else:
            td_iter = _directed_triangles_and_degree_iter(G, nodes)
            clusterc = {v: 0 if t == 0 else t / ((dt * (dt - 1) - 2 * db) * 2)
                        for v, dt, db, t in td_iter}
    else:
        # change to another way
        if weight is not None:
            td_iter = _weighted_triangles_and_otc_iter(G, nodes, weight)
            clusterc = {v: 0 if t == 0 else t / otc for
                        v, t, otc in td_iter}
        else:
            td_iter = _triangles_and_degree_iter(G, nodes)
            clusterc = {v: 0 if t == 0 else t / (d * (d - 1)) for
                        v, d, t, _ in td_iter}
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return clusterc[nodes]
    return clusterc


# KEYFUNC: for closure-co
def closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            td_iter = _directed_weighted_triangles_and_opentriads_iter(G, nodes, weight)
            closurec = {v: [0] if t == 0 else [t / ot]
                        for v, t, ot in td_iter}
        else:
            td_iter = _directed_triangles_and_opentriads_iter(G, nodes)

            closurec = {v: [0, 0, 0] if (ts == 0 and tt == 0) else [(ts + tt) / ot, ts / ot, tt / ot]
                        for v, ts, tt, ot in td_iter}

    # undirected:
    else:
        if weight is not None:
            pass
            # td_iter = _weighted_triangles_and_degree_iter(G, nodes, weight)
            # clusterc = {v: 0 if t == 0 else t / (d * (d - 1)) for
            #             v, d, t in td_iter}
        else:
            td_iter = _triangles_and_opentriads_iter(G, nodes)
            closurec = {v: [0] if t == 0 else [t / ot] for
                        v, t, ot in td_iter}
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return closurec[nodes]
    return closurec

# not used
def src_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            td_iter = _directed_src_triangles_and_opentriads_iter(G, nodes)

            closurec = {v: 0 if ts == 0 else ts / ot
                        for v, ts, ot in td_iter}

    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return closurec[nodes]
    return closurec


# not used
def tgt_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            td_iter = _directed_tgt_triangles_and_opentriads_iter(G, nodes)

            closurec = {v: 0 if tt == 0 else tt / ot
                        for v, tt, ot in td_iter}

    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return closurec[nodes]
    return closurec


def head_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _head_pattern_iter(G, nodes)

            closurec = {v: 0 if th == 0 else th / oth
                        for v, th, oth in pattern_iter}

    if nodes in G:
        return closurec[nodes]
    return closurec


def mid_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _mid_pattern_iter(G, nodes)

            closurec = {v: 0 if tm == 0 else tm / otm
                        for v, tm, otm in pattern_iter}

    if nodes in G:
        return closurec[nodes]
    return closurec


def end_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _end_pattern_iter(G, nodes)

            closurec = {v: 0 if te == 0 else te / ote
                        for v, te, ote in pattern_iter}

    if nodes in G:
        return closurec[nodes]
    return closurec


def cyc_closure(G, nodes=None, weight=None):

    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _cyc_pattern_iter(G, nodes)

            closurec = {v: 0 if tc == 0 else tc / otc
                        for v, tc, otc in pattern_iter}

    if nodes in G:
        return closurec[nodes]
    return closurec



def transitivity(G):
    r"""Compute graph transitivity, the fraction of all possible triangles
    present in G.

    Possible triangles are identified by the number of "triads"
    (two edges with a shared vertex).

    The transitivity is

    .. math::

        T = 3\frac{\#triangles}{\#triads}.

    Parameters
    ----------
    G : graph

    Returns
    -------
    out : float
       Transitivity

    Examples
    --------
    >>> G = nx.complete_graph(5)
    >>> print(nx.transitivity(G))
    1.0
    """
    triangles = sum(t for v, d, t, _ in _triangles_and_degree_iter(G))
    contri = sum(d * (d - 1) for v, d, t, _ in _triangles_and_degree_iter(G))
    return 0 if triangles == 0 else triangles / contri



# proposed algorithm
# should be faster or same?
def quadrangle_coefficient_2(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))
    quad_co = {}

    for v, v_nbrs in nodes_nbrs:
        quad_co[v] = 0
        quad = 0
        potential = 0
        vs = set(v_nbrs) - {v}
        for u in vs:
            u_nbrs = set(G[u]) - {v}
            for w in u_nbrs:
                if w in G[v]:
                    potential += len(G[v]) - 2
                potential += len(G[v]) - 1
                quad += len((set(G[w]) & set(G[v])) - {v}) - 1     # numerator: 2 times number of quadrangles
        if potential > 0:
            quad_co[v] = quad / potential

    if nodes in G:
        return quad_co[nodes]
    return quad_co


# based on square_co below
def quadrangle_coefficient(G, nodes=None):
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    clustering = {}
    for v in node_iter:
        clustering[v] = 0
        potential = 0
        for u, w in combinations(G[v], 2):
            squares = len((set(G[u]) & set(G[w])) - {v}) * 2
            clustering[v] += squares
            degm = 1
            if w in G[u]:
                degm += 1
            potential += (len(G[u]) - degm) + (len(G[w]) - degm)
        if potential > 0:
            clustering[v] /= potential
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return clustering[nodes]
    return clustering


def square_clustering(G, nodes=None):
    r""" Compute the squares clustering coefficient for nodes.


    Parameters
    ----------
    G : graph

    nodes : container of nodes, optional (default=all nodes in G)
       Compute clustering for nodes in this container.

    Returns
    -------
    c4 : dictionary
       A dictionary keyed by node with the square clustering coefficient value.

    Examples
    --------
    >>> G=nx.complete_graph(5)
    >>> print(nx.square_clustering(G,0))
    1.0
    >>> print(nx.square_clustering(G))
    {0: 1.0, 1: 1.0, 2: 1.0, 3: 1.0, 4: 1.0}

    Notes
    -----
    While :math:`C_3(v)` (triangle clustering) gives the probability that
    two neighbors of node v are connected with each other, :math:`C_4(v)` is
    the probability that two neighbors of node v share a common
    neighbor different from v. This algorithm can be applied to both
    bipartite and unipartite networks.

    References
    ----------
    .. [1] Pedro G. Lind, Marta C. González, and Hans J. Herrmann. 2005
        Cycles and clustering in bipartite networks.
        Physical Review E (72) 056127.
       [2] Peng Zhang 2008
        Clustering coefficient and community structure of bipartite networks

    """
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    clustering = {}
    for v in node_iter:
        clustering[v] = 0
        potential = 0
        for u, w in combinations(G[v], 2):
            squares = len((set(G[u]) & set(G[w])) - {v})
            clustering[v] += squares
            degm = squares + 1
            if w in G[u]:
                degm += 1
            potential += (len(G[u]) - degm) + (len(G[w]) - degm) + squares  # changed multiply to addition
        if potential > 0:
            clustering[v] /= potential
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return clustering[nodes]
    return clustering


@not_implemented_for('directed')
def generalized_degree(G, nodes=None):
    r""" Compute the generalized degree for nodes.

    For each node, the generalized degree shows how many edges of given
    triangle multiplicity the node is connected to. The triangle multiplicity
    of an edge is the number of triangles an edge participates in. The
    generalized degree of node :math:`i` can be written as a vector
    :math:`\mathbf{k}_i=(k_i^{(0)}, \dotsc, k_i^{(N-2)})` where
    :math:`k_i^{(j)}` is the number of edges attached to node :math:`i` that
    participate in :math:`j` triangles.

    Parameters
    ----------
    G : graph

    nodes : container of nodes, optional (default=all nodes in G)
       Compute the generalized degree for nodes in this container.

    Returns
    -------
    out : Counter, or dictionary of Counters
       Generalized degree of specified nodes. The Counter is keyed by edge
       triangle multiplicity.

    Examples
    --------
    >>> G=nx.complete_graph(5)
    >>> print(nx.generalized_degree(G,0))
    Counter({3: 4})
    >>> print(nx.generalized_degree(G))
    {0: Counter({3: 4}), 1: Counter({3: 4}), 2: Counter({3: 4}), 3: Counter({3: 4}), 4: Counter({3: 4})}

    To recover the number of triangles attached to a node:

    >>> k1 = nx.generalized_degree(G,0)
    >>> sum([k*v for k,v in k1.items()])/2 == nx.triangles(G,0)
    True

    Notes
    -----
    In a network of N nodes, the highest triangle multiplicty an edge can have
    is N-2.

    The return value does not include a `zero` entry if no edges of a
    particular triangle multiplicity are present.

    The number of triangles node :math:`i` is attached to can be recovered from
    the generalized degree :math:`\mathbf{k}_i=(k_i^{(0)}, \dotsc,
    k_i^{(N-2)})` by :math:`(k_i^{(1)}+2k_i^{(2)}+\dotsc +(N-2)k_i^{(N-2)})/2`.

    References
    ----------
    .. [1] Networks with arbitrary edge multiplicities by V. Zlatić,
        D. Garlaschelli and G. Caldarelli, EPL (Europhysics Letters),
        Volume 97, Number 2 (2012).
        https://iopscience.iop.org/article/10.1209/0295-5075/97/28005
    """
    if nodes in G:
        return next(_triangles_and_degree_iter(G, nodes))[3]
    return {v: gd for v, d, t, gd in _triangles_and_degree_iter(G, nodes)}
