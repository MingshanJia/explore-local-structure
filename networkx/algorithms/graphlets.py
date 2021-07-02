from itertools import combinations

__all__ = ['four_clique', 'four_cycle_plus', 'four_cycle_plus_2']


# four cycle plus (4-cycle+) is 4-cycle with one diagonal link
def four_cycle_plus(G, nodes=None):
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    res = {}
    for v in node_iter:
        # k = G.degree(v)
        res[v] = 0

        for u, w in combinations(G[v], 2):
            squares = len(set(G[u]) & set(G[w]) & set(G[v]))
            res[v] += squares

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


def four_clique(G, nodes=None):
    if nodes is None:
        node_iter = G
    else:
        node_iter = G.nbunch_iter(nodes)
    res = {}
    for i in node_iter:
        # k = G.degree(v)
        res[i] = 0

        for u, v, w in combinations(G[i], 3):
            if (w in set(G[u]) & set(G[v])) and (u in set(G[v])):
                res[i] += 1
    if nodes in G:
        return res[nodes]
    return res