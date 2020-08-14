from itertools import chain
from itertools import combinations
from collections import Counter
import networkx as nx

from networkx.utils import not_implemented_for

__all__ = ['quadrangle_coefficient_iter', 'global_quadrangle',
           'square_clustering', 'quadrangle_coefficient', 'inner_quadrangle_coefficient','outer_quadrangle_coefficient',
           'quad_iquad_oquad', 'average_inner_quad_co', 'average_outer_quad_co']


#TODO: caculate the two simultaneously
# **************************************************************************** Quadrangle Coefficient ************************************************************
# for calculating inner-quad-co and outer-quad-co
def quadrangle_coefficient_iter(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for i, nbrs in nodes_nbrs:
        quad = 0
        inner_quad = 0
        outer_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                quad += len((knbrs & inbrs) - {j})  # numerator: 2 times number of quadrangles
                inner_quad += len(inbrs - {j} - {k})
                outer_quad += len(knbrs - {j} - {i})
            # for k in jnbrs:
            #     if k in G[i]:
            #         inner_quad += len(G[i]) - 2
            #         outer_quad += len(set(G[k])) - 2
            #     else:
            #         inner_quad += len(G[i]) - 1
            #         outer_quad += len(set(G[k])) - 1
            #     quad += len((set(G[k]) & set(G[i])) - {i} - {k}) - 1  # numerator: 2 times number of quadrangles
        yield (i, quad, inner_quad, outer_quad)


# print number of quad, iquad and oquad
def quad_iquad_oquad(G, nodes=None):

    qc_iter = quadrangle_coefficient_iter(G, nodes)
    res = {v: [q, iq, oq]  for v, q, iq, oq in qc_iter}
    if nodes in G:
        return res[nodes]
    return res


def i_quad_coef_iter(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for i, nbrs in nodes_nbrs:
        quad = 0
        inner_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                quad += len((knbrs & inbrs) - {j})  # numerator: 2 times number of quadrangles
                inner_quad += len(inbrs - {j} - {k})
        yield (i, quad, inner_quad)


def o_quad_coef_iter(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for i, nbrs in nodes_nbrs:
        quad = 0
        outer_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                quad += len((knbrs & inbrs) - {j})  # numerator: 2 times number of quadrangles
                outer_quad += len(knbrs - {j} - {i})
        yield (i, quad, outer_quad)


# for calculating weighted inner-quad-co and outer-quad-co
# def weighted_quadrangle_coefficient_iter(G, nodes=None, weight='weight'):
#
#     if weight is None or G.number_of_edges() == 0:
#         max_weight = 1
#     else:
#         max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))
#
#     def wt(u, v):
#         return G[u][v].get(weight, 1) / max_weight
#
#     if nodes is None:
#         nodes_nbrs = G.adj.items()
#     else:
#         nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))
#
#     for i, nbrs in nodes_nbrs:
#         weighted_quad = 0
#         weighted_inner_quad = 0
#         weighted_outer_quad = 0
#         inbrs = set(nbrs) - {i}
#         for j in inbrs:
#             jnbrs = set(G[j]) - {i} - {j}
#             for k in jnbrs:
#                 knbrs = set(G[k]) - {k}
#                 weighted_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) * wt(k, l) for l in
#                                      (inbrs & knbrs - {j}))  # numerator: 2 times number of quadrangles
#                 weighted_inner_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) for l in (inbrs - {j} - {k}))
#                 weighted_outer_quad += sum(wt(i, j) * wt(j, k) * wt(k, l) for l in (knbrs - {j} - {i}))
#         yield (i, weighted_quad, weighted_inner_quad, weighted_outer_quad)


def i_weighted_quad_coef_iter(G, nodes=None, weight='weight'):

    if weight is None or G.number_of_edges() == 0:
        max_weight = 1
    else:
        max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))

    def wt(u, v):
        return G[u][v].get(weight, 1) / max_weight

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for i, nbrs in nodes_nbrs:
        weighted_quad = 0
        weighted_inner_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                weighted_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) * wt(k, l) for l in
                                     (inbrs & knbrs - {j}))  # numerator: 2 times number of quadrangles
                weighted_inner_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) for l in (inbrs - {j} - {k}))
        yield (i, weighted_quad, weighted_inner_quad)


def o_weighted_quad_coef_iter(G, nodes=None, weight='weight'):

    if weight is None or G.number_of_edges() == 0:
        max_weight = 1
    else:
        max_weight = max(d.get(weight, 1) for u, v, d in G.edges(data=True))

    def wt(u, v):
        return G[u][v].get(weight, 1) / max_weight

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    for i, nbrs in nodes_nbrs:
        weighted_quad = 0
        weighted_outer_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                weighted_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) * wt(k, l) for l in
                                     (inbrs & knbrs - {j}))  # numerator: 2 times number of quadrangles
                weighted_outer_quad += sum(wt(i, j) * wt(j, k) * wt(k, l) for l in (knbrs - {j} - {i}))
        yield (i, weighted_quad, weighted_outer_quad)


# inner quadrangle coefficient
def inner_quadrangle_coefficient(G, nodes=None, weight=None):
    if weight is not None:
        qc_iter = i_weighted_quad_coef_iter(G, nodes, weight)
        inner_quad_co = {v: 0 if q == 0 else q / in_q for
                         v, q, in_q in qc_iter}
    else:
        qc_iter = i_quad_coef_iter(G, nodes)
        inner_quad_co = {v: 0 if q == 0 else q / in_q for
                    v, q, in_q, in qc_iter}
    if nodes in G:
        return inner_quad_co[nodes]
    return inner_quad_co


# outer quadrangle coefficient
def outer_quadrangle_coefficient(G, nodes=None, weight=None):
    if weight is not None:
        qc_iter = o_weighted_quad_coef_iter(G, nodes, weight)
        outer_quad_co = {v: 0 if q == 0 else q / out_q for
                         v, q, out_q in qc_iter}
    else:
        qc_iter = o_quad_coef_iter(G, nodes)
        outer_quad_co = {v: 0 if q == 0 else q / out_q for
                    v, q, out_q in qc_iter}
    if nodes in G:
        return outer_quad_co[nodes]
    return outer_quad_co


# average inner quadrangle coefficient
def average_inner_quad_co(G, nodes=None, weight=None, count_zeros=True):

    c = inner_quadrangle_coefficient(G, nodes, weight=weight).values()
    if not count_zeros:
        c = [v for v in c if v > 0]
    return sum(c) / len(c)


# average outer quadrangle coefficient
def average_outer_quad_co(G, nodes=None, weight=None, count_zeros=True):

    c = outer_quadrangle_coefficient(G, nodes, weight=weight).values()
    if not count_zeros:
        c = [v for v in c if v > 0]
    return sum(c) / len(c)


# get global quadrangle coef., same for i-quad and o-quad
def global_quadrangle(G, nodes=None, weight=None):
    if weight is not None:
        wqd_iter = i_weighted_quad_coef_iter(G, nodes, weight)
        q_iq_list = [[q, iq] for _, q, iq in wqd_iter]
        q_sum = sum(row[0] for row in q_iq_list)
        iq_sum = sum(row[1] for row in q_iq_list)
    else:
        qd_iter = i_quad_coef_iter(G, nodes)
        q_iq_list = [[q, iq] for _, q, iq in qd_iter]
        q_sum = sum(row[0] for row in q_iq_list)
        iq_sum = sum(row[1] for row in q_iq_list)
    return q_sum / iq_sum


# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# obsolete
# based on square_co below
# same value as inner_quad_co
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