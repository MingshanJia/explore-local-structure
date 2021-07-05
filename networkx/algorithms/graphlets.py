from itertools import combinations
from collections import Counter

__all__ = ['typed_edge_graphlet_degree_vector_ego', 'graphlet_degree_vector_ego', 'three_wedge', 'four_clique',
           'four_cycle_plus', 'four_cycle_plus_2']


# TyE-GDV : GDV with edge type information, returning a 7*num_type 2D list
# note that the type of edges not incident to ego node is not considered
# because there is no such information in the targeted dataset.
# TODO: take all edges of graphlets into consideration
def typed_edge_graphlet_degree_vector_ego(G, num_type, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        vec = [[0] * num_type for _ in range(7)]  # 7 ego-centric graphlets by number of types of edges
        inbrs = set(i_nbrs) - {i}

        # 2-clique, this is also the classic
        for u in inbrs:
            iu_type = G.get_edge_data(i, u)['edge_type']
            edit_vec(vec, 0, iu_type, None, None)

        # 3-node graphlets
        for u, v in combinations(inbrs, 2):
            u_nbrs = set(G[u]) - {u}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']

            # 2-path
            edit_vec(vec, 1, iu_type, iv_type, None)

            # 3-clique
            if v in u_nbrs:
                edit_vec(vec, 2, iu_type, iv_type, None)

                # 4-node graphlets
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']
            iw_type = G.get_edge_data(i, w)['edge_type']

            # 3-star
            edit_vec(vec, 3, iu_type, iv_type, iw_type)

            # tailed-tri
            if w in u_nbrs:
                edit_vec(vec, 4, iu_type, iv_type, iw_type)
            if v in w_nbrs:
                edit_vec(vec, 4, iu_type, iv_type, iw_type)
            if v in u_nbrs:
                edit_vec(vec, 4, iu_type, iv_type, iw_type)

            # 4-cycle+
            if w in u_nbrs & v_nbrs:
                edit_vec(vec, 5, iu_type, iv_type, iw_type)
            if u in w_nbrs & v_nbrs:
                edit_vec(vec, 5, iu_type, iv_type, iw_type)
            if v in u_nbrs & w_nbrs:
                edit_vec(vec, 5, iu_type, iv_type, iw_type)

            # 4-clique
            if (w in u_nbrs & v_nbrs) and (u in v_nbrs):
                edit_vec(vec, 6, iu_type, iv_type, iw_type)

        res[i] = vec
    if nodes in G:
        return res[nodes]
    return res


# edge type is coded from 1, edges without type information is coded as 0.
def edit_vec(vec, idx, link_1, link_2, link_3):
    if int(link_1) > 0:
        vec[idx][int(link_1) - 1] += 1
    if link_2 != None and int(link_2) > 0:
        vec[idx][int(link_2) - 1] += 1
    if link_3 != None and int(link_3) > 0:
        vec[idx][int(link_3) - 1] += 1


# Ego-GDV = [2-clique, 2-wedge, 3-clique, 3-star, 3-wedge, 4-cycle+, 4-clique]
def graphlet_degree_vector_ego(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        vs = set(i_nbrs) - {i}
        k = len(vs)
        gen_degree = Counter(len(vs & (set(G[w]) - {w})) for w in vs)
        T = sum(k * val for k, val in gen_degree.items()) // 2

        four_cycle_plus = 0
        four_clique = 0
        for u, v, w in combinations(vs, 3):
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
        for u, w in combinations(set(G[v]) - {v}, 2):
            tmp = len((set(G[u]) - {u}) & (set(G[w]) - {w}) & (set(G[v]) - {v}))
            res[v] += tmp
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
        for u, v, w in combinations(set(G[i]) - {i}, 3):
            if (w in (set(G[u]) -{u}) & (set(G[v]) - {v})) and (u in (set(G[v]) - {v})):
                res[i] += 1
    if nodes in G:
        return res[nodes]
    return res


# obselete: another way of calculating 4-cycle+
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


