from itertools import combinations
from collections import Counter

__all__ = ['graphlet_vector_ego', 'three_wedge', 'four_clique', 'four_cycle_plus', 'four_cycle_plus_2']


# GV = [2-clique, 2-wedge, 3-clique, 3-star, 3-wedge, 4-cycle+, 4-clique]
def graphlet_vector_ego(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        k = G.degree(i)
        vs = set(i_nbrs) - {i}
        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        T = sum(k * val for k, val in gen_degree.items()) // 2

        four_cycle_plus = 0
        four_clique = 0
        for u, v, w in combinations(G[i], 3):
            if (w in (set(G[u]) - {u}) & (set(G[v]) - {v})):
                four_cycle_plus += 1
            if (u in (set(G[w]) - {w}) & (set(G[v]) - {v})):
                four_cycle_plus += 1
            if (v in (set(G[u]) - {u}) & (set(G[w]) - {w})):
                four_cycle_plus += 1

            if (w in (set(G[u]) - {u}) & (set(G[v]) - {v})) and (u in (set(G[v]) - {v})):
                four_clique += 1

        vec = [k, k * (k - 1) // 2, T, k * (k - 1) * (k - 2) // 6, T * (k - 2), four_cycle_plus, four_clique]
        res[i] = vec
    if nodes in G:
        return res[nodes]
    return res


# 3-wedge is 3-clique + 1 link
def three_wedge(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for v, v_nbrs in nodes_nbrs:
        k = G.degree(v)
        vs = set(v_nbrs) - {v}
        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        # Counter({2:3, 0:1}) means that 3 node form 2 triangles, 1 node forms 0 triangle
        ntriangles = sum(k * val for k, val in gen_degree.items())
        res[v] = (ntriangles // 2) * (k - 2)
    if nodes in G:
        return res[nodes]
    return res


# four cycle plus (4-cycle+) is 4-cycle with one diagonal link
def four_cycle_plus(G, nodes=None):
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    res = {}
    for v in node_iter:
        res[v] = 0

        for u, w in combinations(G[v], 2):
            squares = len(set(G[u]) & set(G[w]) & set(G[v]))
            res[v] += squares

    if nodes in G:
        return res[nodes]
    return res


def four_clique(G, nodes=None):
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    res = {}
    for i in node_iter:
        res[i] = 0

        for u, v, w in combinations(G[i], 3):
            if (w in set(G[u]) & set(G[v])) and (u in set(G[v])):
                res[i] += 1
    if nodes in G:
        return res[nodes]
    return res


# another way of calculating 4-cycle+
def four_cycle_plus_2(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, nbrs in nodes_nbrs:
        res[i] = 0
        inbrs = set(nbrs) - {i}
        for j in inbrs:
            jnbrs = set(G[j]) - {i} - {j}
            for k in jnbrs & inbrs:
                knbrs = set(G[k]) - {k}
                res[i] += len((knbrs & inbrs) - {j})  # numerator: 2 times number of quadrangles
        res[i] = res[i] // 2
    if nodes in G:
        return res[nodes]
    return res


