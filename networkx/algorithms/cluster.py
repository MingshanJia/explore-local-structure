"""Algorithms to characterize the number of triangles in a graph."""

from itertools import chain
from itertools import combinations
from collections import Counter
import networkx as nx

from networkx.utils import not_implemented_for

__all__ = ['triangles', 'number_of_triangles', 'average_clustering', 'clustering', 'transitivity', 'triangles_and_otc',
           'triangles_and_ote', 'global_clustering', 'generalized_degree', 'average_closure', 'closure', 'clustering_closure_coefs',
           'src_closure', 'tgt_closure', 'head_closure', 'mid_closure', 'end_closure', 'cyc_closure',
           'four_clustering_patterns', 'average_eight_patterns', 'four_closure_patterns', 'eight_patterns']



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
        tri = 0
        for w in vs:
            ns = set(G[w]) - {w}
            ot += len(ns) - 1
            tri += len(vs & ns)

        yield (v, tri, ot)


@not_implemented_for('multigraph')
def _weighted_triangles_and_opentriads_iter(G, nodes=None, weight='weight'):

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

    for i, inbrs in nodes_nbrs:
        inbrs = set(inbrs) - {i}
        weighted_triangles = 0
        weighted_end_opentriads = 0
        for j in inbrs:
            jnbrs = set(G[j]) - {j}
            weighted_triangles += sum((wt(i, j) * wt(j, k) * wt(k, i))
                                      for k in inbrs & jnbrs)
            weighted_end_opentriads += sum(wt(i, j) * wt(j, k) for k in (jnbrs- {i}))
        yield (i, weighted_triangles, weighted_end_opentriads)


@not_implemented_for('multigraph')
def _triangles_and_otc_ote_iter(G, nodes=None):

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for v, v_nbrs in nodes_nbrs:
        vs = set(v_nbrs) - {v}
        ote = 0
        tri = 0
        for w in vs:
            ns = set(G[w]) - {w}
            ote += len(ns) - 1
            tri += len(vs & ns)

        yield (v, tri, len(vs), ote)


# calculate clustering and closure coefs at the same time
def clustering_closure_coefs(G, nodes=None, weight=None):
    if G.is_directed():
        pass
    else:
        if weight is not None:
            pass
        else:
            td_iter = _triangles_and_otc_ote_iter(G, nodes)
            cluclo = {v: [0, 0] if t == 0 else [t / (d * (d-1)), t / ote] for
                        v, t, d, ote in td_iter}
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return cluclo[nodes]
    return cluclo


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


##------------------------------------------ directed triangle formation patterns ------------------------------------------------------------
# Key Function
def _directed_triangle_patterns_iter(G, nodes=None):
    nodes_nbrs = ((n, G._pred[n], G._succ[n]) for n in G.nbunch_iter(nodes))

    for i, preds, succs in nodes_nbrs:
        ipreds = set(preds) - {i}
        isuccs = set(succs) - {i}
        di_out = len(isuccs)
        di_in = len(ipreds)
        di_in_out = len(ipreds & isuccs)

        t_head = 0
        t_end = 0
        t_mid = 0
        t_cyc = 0

        ote_head = 0
        ote_end = 0
        ote_mid = 0
        ote_cyc = 0

        # otc_head and otc_end are actually two times the number of open triads, because the closing edge can take two directions
        otc_head = di_out * (di_out - 1)
        otc_end = di_in * (di_in - 1)
        otc_mid_cyc = di_in * di_out - di_in_out

        # calculating triangles in four patterns: t_head, t_end, t_mid and t_cyc here are actually two times the number of triangles of each pattern.
        # for the four clustering patterns, we need to divide them by two.
        # TODO: debug t_mid, t_cyc
        for j in isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            t_head += sum(1 for k in
                          chain((isuccs & jpreds),
                                (isuccs & jsuccs)))
            t_mid += sum(1 for k in ipreds & jpreds)
            t_cyc += sum(1 for k in ipreds & jsuccs)

        for j in ipreds:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}

            t_end += sum(1 for k in
                         chain((ipreds & jpreds),
                               (ipreds & jsuccs)))
            t_mid += sum(1 for k in isuccs & jsuccs)
            t_cyc += sum(1 for k in isuccs & jpreds)

        # calculating open triads in four patterns
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ote_head += dj - 2
            ote_end += dj - 2
            ote_mid += dj -2
            ote_cyc += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jsuccs = set(G._succ[j]) - {j}
            jpreds = set(G._pred[j]) - {j}
            dj_out = len(jsuccs)
            dj_in = len(jpreds)
            dj = dj_out + dj_in

            ote_head += dj - 1
            ote_mid += dj_in - 1
            ote_cyc += dj_out

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj_in = len(jpreds)
            dj_out = len(jsuccs)
            dj = dj_out + dj_in

            ote_end += dj - 1
            ote_mid += dj_out - 1
            ote_cyc += dj_in
        # for testing
        # print("t_head:{};  t_end:{};  t_mid:{};  t_cyc:{};   ote_head:{};   ote_end:{};   ote_mid:{};   ote_cyc:{};  \
        #  otc_head:{};   otc_end:{};   otc_mid_cyc:{}\n".format(t_head, t_end, t_mid, t_cyc, ote_head, ote_end, ote_mid, ote_cyc, otc_head, otc_end, otc_mid_cyc))
        yield (i, t_head, t_end, t_mid, t_cyc, ote_head, ote_end, ote_mid, ote_cyc, otc_head, otc_end, otc_mid_cyc)


# not used, for head-of-path pattern
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

        ote_head = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            ote_head += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ote_head += dj - 1
        yield (i, t_head, ote_head)


# not used, for end-of-path pattern
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

        ote_end = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)
            ote_end += dj - 2

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ote_end += dj - 1

        yield (i, t_end, ote_end)


# not used, for mid-of-path pattern
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

        ote_mid = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ote_mid += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj_in = len(jpreds)

            ote_mid += dj_in - 1

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj_out = len(jsuccs)

            ote_mid += dj_out - 1

        yield (i, t_mid, ote_mid)


# not used, for cyclic pattern
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

        ote_cyc = 0
        for j in ipreds & isuccs:
            jpreds = set(G._pred[j]) - {j}
            jsuccs = set(G._succ[j]) - {j}
            dj = len(jpreds) + len(jsuccs)

            ote_cyc += dj - 2

        for j in (isuccs - (ipreds & isuccs)):
            jsuccs = set(G._succ[j]) - {j}
            dj_out = len(jsuccs)

            ote_cyc += dj_out

        for j in (ipreds - (ipreds & isuccs)):
            jpreds = set(G._pred[j]) - {j}
            dj_in = len(jpreds)

            ote_cyc += dj_in

        yield (i, t_cyc, ote_cyc)


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


# closure coefficient
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
            td_iter = _weighted_triangles_and_opentriads_iter(G, nodes, weight)
            closurec = {v: [0] if t == 0 else [t / ot] for
                         v, t, ot in td_iter}
        else:
            td_iter = _triangles_and_opentriads_iter(G, nodes)
            closurec = {v: [0] if t == 0 else [t / ot] for
                        v, t, ot in td_iter}
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return closurec[nodes]
    return closurec


# to get number of triangels and centre-node-based open triads
def triangles_and_otc(G, nodes=None):
    td_iter =  _triangles_and_degree_iter(G, nodes)
    tri_otc = {v : [t, d*(d-1)] for v, d, t, _ in td_iter}
    if nodes in G:
        return tri_otc[nodes]
    return tri_otc

# to get number of triangels and end-node-based open triads
def triangles_and_ote(G, nodes=None):
    td_iter = _triangles_and_opentriads_iter(G, nodes)
    tri_ote = {v : [t, ote] for v, t, ote in td_iter}
    if nodes in G:
        return tri_ote[nodes]
    return tri_ote


# get global clustering, same for global closure
def global_clustering(G, nodes=None, weight=None):
    if G.is_directed():
        if weight is not None:
            pass
        else:
            pass
    else:
        if weight is not None:
            pass
        else:
            td_iter = _triangles_and_degree_iter(G, nodes)
            t_ot_list = [[t, d*(d-1)] for _, d, t, _ in td_iter]
            t_sum = sum(row[0] for row in t_ot_list)
            ot_sum = sum(row[1] for row in t_ot_list)

    return t_sum / ot_sum


@not_implemented_for('directed')
def number_of_triangles(G):
    td_iter = _triangles_and_degree_iter(G)
    t_ot_list = [[t, d * (d - 1)] for _, d, t, _ in td_iter]
    t_sum = sum(row[0] for row in t_ot_list)
    return t_sum / 6


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


# obsolete, Calculate four closure patterns at the same time
def four_closure_patterns(G, nodes=None, weight=None):
    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _directed_triangle_patterns_iter(G, nodes)
            res = {}
            for v, t_h, t_e, t_m, t_c, ote_h, ote_e, ote_m, ote_c, _, _, _ in pattern_iter:
                head = 0 if t_h == 0 else t_h / ote_h
                end = 0 if t_e == 0 else t_e / ote_e
                mid = 0 if t_m == 0 else t_m / ote_m
                cyc = 0 if t_c == 0 else t_c / ote_c
                res[v] = [head, end, mid, cyc]
    if nodes in G:
        return res[nodes]
    return res


# obsolete, Calculate four clustering patterns at the same time
def four_clustering_patterns(G, nodes=None, weight=None):
    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _directed_triangle_patterns_iter(G, nodes)
            res = {}
            for v, t_h, t_e, t_m, t_c, _, _, _, _, otc_h, otc_e, otc_mc in pattern_iter:
                head = 0 if t_h == 0 else t_h / otc_h
                end = 0 if t_e == 0 else t_e / otc_e
                mid = 0 if t_m == 0 else t_m / otc_mc
                cyc = 0 if t_c == 0 else t_c / otc_mc
                res[v] = [head, end, mid, cyc]
    if nodes in G:
        return res[nodes]
    return res


# Calculate four closure patterns and four clustering patterns at the same time
def eight_patterns(G, nodes=None, weight=None):
    if G.is_directed():
        if weight is not None:
            pass
        else:
            pattern_iter = _directed_triangle_patterns_iter(G, nodes)
            res = {}
            for v, t_h, t_e, t_m, t_c, ote_h, ote_e, ote_m, ote_c, otc_h, otc_e, otc_mc in pattern_iter:
                clo_head = 0 if t_h == 0 else t_h / ote_h
                clo_end = 0 if t_e == 0 else t_e / ote_e
                clo_mid = 0 if t_m == 0 else t_m / ote_m
                clo_cyc = 0 if t_c == 0 else t_c / ote_c

                # divide clustering patterns by two
                clu_head = 0 if t_h == 0 else t_h / (2 * otc_h)
                clu_end = 0 if t_e == 0 else t_e / (2 * otc_e)
                clu_mid = 0 if t_m == 0 else t_m / (2 * otc_mc)
                clu_cyc = 0 if t_c == 0 else t_c / (2 * otc_mc)
                res[v] = [clo_head, clo_end, clo_mid, clo_cyc, clu_head, clu_end, clu_mid, clu_cyc]
    if nodes in G:
        return res[nodes]
    return res


def average_eight_patterns(G, nodes=None, weight=None, count_zeros=True):
    patterns = eight_patterns(G, nodes, weight=weight).values()
    clo_head = []
    clo_end = []
    clo_mid = []
    clo_cyc = []
    clu_head = []
    clu_end = []
    clu_mid = []
    clu_cyc = []
    for v in patterns:
        clo_head.append(v[0])
        clo_end.append(v[1])
        clo_mid.append(v[2])
        clo_cyc.append(v[3])
        clu_head.append(v[4])
        clu_end.append(v[5])
        clu_mid.append(v[6])
        clu_cyc.append(v[7])
    if not count_zeros:
        clo_head = [v for v in clo_head if v > 0]
        clo_end = [v for v in clo_end if v > 0]
        clo_mid = [v for v in clo_mid if v > 0]
        clo_cyc = [v for v in clo_cyc if v > 0]
        clu_head = [v for v in clu_head if v > 0]
        clu_end = [v for v in clu_end if v > 0]
        clu_mid = [v for v in clu_mid if v > 0]
        clu_cyc = [v for v in clu_cyc if v > 0]
    return [sum(clo_head) / len(clo_head), sum(clo_end) / len(clo_end), sum(clo_mid) / len(clo_mid), sum(clo_cyc) / len(clo_cyc), \
           sum(clu_head) / len(clu_head), sum(clu_end) / len(clu_end), sum(clu_mid) / len(clu_mid), sum(clu_cyc) / len(clu_cyc)]



# not used
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

# not used
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

# not used
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

# not used
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
