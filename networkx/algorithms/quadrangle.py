from itertools import chain
from itertools import combinations
from collections import Counter
import networkx as nx

from networkx.utils import not_implemented_for

__all__ = ['quadrangle_coefficient_iter', 'global_quadrangle', 'number_of_quadrangles', 'square_clustering', 'order_two_clustering',
           'square_clustering_2', 'quadrangle_coefficient', 'inner_quadrangle_coefficient','outer_quadrangle_coefficient',
           'iquad_oquad_coefs','quad_iquad_oquad', 'average_inner_quad_co', 'average_outer_quad_co', 'average_inner_and_outer_quad_co',
           'order_three_clustering_coef','primary_grid_coef']


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


# calculate iquad and oquad at the same time
def iquad_oquad_coefs(G, nodes=None, weight=None):

    if weight is not None:
        wqc_iter = weighted_iquad_oquad_coef_iter(G, nodes, weight)
        res = {v: [0, 0] if q == 0 else [q / iq, q / oq] for v, q, iq, oq in wqc_iter}
    else:
        qc_iter = quadrangle_coefficient_iter(G, nodes)
        res = {v: [0, 0] if q == 0 else [q / iq, q / oq] for v, q, iq, oq in qc_iter}
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


# calculate iquad and oquad at the same time
def weighted_iquad_oquad_coef_iter(G, nodes=None, weight='weight'):

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
        weighted_outer_quad = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs:
                knbrs = set(G[k]) - {k}
                weighted_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) * wt(k, l) for l in
                                     (inbrs & knbrs - {j}))  # numerator: 2 times number of quadrangles
                weighted_inner_quad += sum(wt(i, j) * wt(j, k) * wt(i, l) for l in (inbrs - {j} - {k}))
                weighted_outer_quad += sum(wt(i, j) * wt(j, k) * wt(k, l) for l in (knbrs - {j} - {i}))
        yield (i, weighted_quad, weighted_inner_quad, weighted_outer_quad)


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


# average inner and outer quadrangle coefficient
def average_inner_and_outer_quad_co(G, nodes=None, weight=None, count_zeros=True):

    iq_oq = iquad_oquad_coefs(G, nodes, weight=weight)
    list_iquad = []
    list_oquad = []
    for k, v in iq_oq.items():
        list_iquad.append(v[0])
        list_oquad.append(v[1])
    if not count_zeros:
        list_iquad = [v for v in list_iquad if v > 0]
        list_oquad = [v for v in list_oquad if v > 0]
    return sum(list_iquad) / len(list_iquad), sum(list_oquad) / len(list_oquad)


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


@not_implemented_for('directed')
def number_of_quadrangles(G):
    qd_iter = i_quad_coef_iter(G)
    q_iq_list = [[q, iq] for _, q, iq in qd_iter]
    q_sum = sum(row[0] for row in q_iq_list)
    return q_sum / 8


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


#Higher order clustering: another square clustering 2002 (problem : could larger than 1)
def order_two_clustering(G, nodes=None):
    """Paper: Higher-order clustering coefficients in barabasi–albert networks"""

    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    clustering = {}
    for v in node_iter:
        clustering[v] = 0
        for u, w in combinations(G[v], 2):
            squares = len((set(G[u]) & set(G[w])) - {v})
            clustering[v] += squares
        degree = len(set(G[v]) - {v})
        denominator = degree * (degree - 1) / 2
        if denominator > 0:
            clustering[v] /= denominator
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return clustering[nodes]
    return clustering


# Pedro G. Lind, Marta C. González, and Hans J. Herrmann. 2005 Cycles and clustering in bipartite networks. Physical Review E (72) 056127.
def square_clustering(G, nodes=None):

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
            potential += (len(G[u]) - degm) * (len(G[w]) - degm) + squares
        if potential > 0:
            clustering[v] /= potential
    if nodes in G:
        # Return the value of the sole entry in the dictionary.
        return clustering[nodes]
    return clustering


# Peng Zhang 2008 Clustering coefficient and community structure of bipartite networks
def square_clustering_2(G, nodes=None):
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


def primary_grid_coef(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))
    res = {}
    for i, i_nbrs in nodes_nbrs:
        res[i] = 0
        four_cycle_plus = 0
        vs = set(i_nbrs) - {i}
        k = len(vs)

        for u, v, w in combinations(vs, 3):
            if (w in (set(G[u]) - {u}) & (set(G[v]) - {v})):
                four_cycle_plus += 1
            if (u in (set(G[w]) - {w}) & (set(G[v]) - {v})):
                four_cycle_plus += 1
            if (v in (set(G[u]) - {u}) & (set(G[w]) - {w})):
                four_cycle_plus += 1
        if k > 2:
            val = 2 * four_cycle_plus / (k * (k - 1) * (k - 2))
            res[i] = val
    return res


# Defined in paper: Higher-order clustering in networks 2018
def order_three_clustering_coef(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))
    res = {}
    for i, i_nbrs in nodes_nbrs:
        res[i] = 0
        four_clique = 0
        vs = set(i_nbrs) - {i}
        k = len(vs)
        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        T = sum(k * val for k, val in gen_degree.items()) // 2
        denominator = T * (k - 2)
        for u, v, w in combinations(vs, 3):
            if (w in (set(G[u]) - {u}) & (set(G[v]) - {v})) and (u in (set(G[v]) - {v})):
                four_clique += 1
    
        if denominator > 0:
            val = 3 * four_clique / denominator
            res[i] = val
    
    if nodes in G:
        return res[nodes]
    return res


