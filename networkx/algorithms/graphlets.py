from itertools import combinations
from collections import Counter

__all__ = ['typed_edge_induced_graphlet_degree_vector_ego', 'typed_edge_graphlet_degree_vector_ego',
           'induced_graphlet_degree_vector', 'typed_edge_induced_graphlet_degree_vector',
           'induced_graphlet_degree_vector_ego', 'graphlet_degree_vector_ego', 'three_wedge', 'four_clique',
           'four_cycle_plus', 'four_cycle_plus_2', 'induced_graphlet_degree_vector_v2', 'colored_graphlet_vector_for_typed_edge']

# v_2 is implemented without using combination, and therefore include a lot of repetitions in calculation
def induced_graphlet_degree_vector_v2(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        inbrs = set(i_nbrs) - {i}

        orbit_0 = len(inbrs)
        # 3-node graphlets
        orbit_1 = 0
        orbit_2 = 0
        orbit_3 = 0
        # i-quad, o-quad, 4-cycle
        orbit_4 = 0
        orbit_5 = 0
        orbit_8 = 0
        # based on orbit 6
        orbit_6 = 0
        orbit_9 = 0
        orbit_10 = 0
        orbit_12 = 0
        orbit_13 = 0
        orbit_14 = 0
        # based on orbit-7
        orbit_7 = 0
        orbit_11 = 0

        for j in inbrs:
            jnbrs = set(G[j]) - {j}

            # orbit-2, orbit-3 are two times of the actual number
            for k in inbrs - {j}:
                knbrs = set(G[k]) - {k}
                if k not in jnbrs:
                    orbit_2 += 1
                else:
                    orbit_3 += 1

                # orbit 7 is six times the actual number
                for l in inbrs - {j} - {k}:
                    lnbrs = set(G[l]) - {l}
                    if (j not in knbrs) and (j not in lnbrs) and (k not in lnbrs):
                        orbit_7 += 1
                    # orbit 11 is six times the actual number
                    if (j in knbrs) and (l not in jnbrs) and (l not in knbrs):
                        orbit_11 += 1
                    if (l in knbrs) and (j not in lnbrs) and (j not in knbrs):
                        orbit_11 += 1
                    if (j in lnbrs) and (k not in jnbrs) and (k not in lnbrs):
                        orbit_11 += 1

            for k in (jnbrs - {i}):
                knbrs = set(G[k]) - {k}

                # orbit-1
                if i not in knbrs:
                    orbit_1 += 1

                # orbit-4, orbit-8
                for l in (knbrs - {i} - {j}):
                    if l not in inbrs and l not in jnbrs and k not in inbrs:
                        orbit_4 += 1
                    # orbit-8 are two times of the actual number)
                    if l in inbrs and l not in jnbrs and k not in inbrs:
                        orbit_8 += 1

                # orbit-5
                for l in (inbrs - {j}):
                    if l not in jnbrs and l not in knbrs and k not in inbrs:
                        orbit_5 += 1

                # based on orbit 6 (uninduced)
                # orbit 6, 9, 10, 13 are counted twice
                for l in jnbrs - {i} - {k}:
                    lnbrs = set(G[l]) - {l}
                    # orbit 6
                    if k not in inbrs and l not in inbrs and k not in lnbrs:
                        orbit_6 += 1
                    # orbit 9
                    if k not in inbrs and l not in inbrs and k in lnbrs:
                        orbit_9 += 1
                    # orbit 10
                    if (k in inbrs and k not in lnbrs and l not in inbrs) or (
                            l in inbrs and l not in knbrs and k not in inbrs):
                        orbit_10 += 1
                    # orbit 12 is counted 4 times!
                    if (k in lnbrs and l in inbrs and k not in inbrs) or (k in lnbrs and l not in inbrs and k in inbrs):
                        orbit_12 += 1
                    # orbit #13
                    if k in inbrs and l in inbrs and k not in lnbrs:
                        orbit_13 += 1
                    # orbit-14 is six times of the actual number)
                    if k in inbrs and l in inbrs and k in lnbrs:
                        orbit_14 += 1

        vec = [orbit_0, orbit_1, orbit_2//2, orbit_3//2, orbit_4, orbit_5, orbit_6//2, orbit_7//6, orbit_8//2,
               orbit_9//2, orbit_10//2, orbit_11//6, orbit_12//4, orbit_13//2, orbit_14//6]
        res[i] = vec
    return res


# taking into account all 15 orbits
# compare this with v_2, discuss the time complexity
# use combination to avoid repetition
def induced_graphlet_degree_vector(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        inbrs = set(i_nbrs) - {i}

        orbit_0 = len(inbrs)
        # 3-node graphlets
        orbit_1 = 0
        orbit_2 = 0
        orbit_3 = 0
        # i-quad, o-quad, 4-cycle
        orbit_4 = 0
        orbit_5 = 0
        orbit_8 = 0
        # based on orbit 6
        orbit_6 = 0
        orbit_9 = 0
        orbit_10 = 0
        orbit_12 = 0
        orbit_13 = 0
        orbit_14 = 0
        # based on orbit-7
        orbit_7 = 0
        orbit_11 = 0

        for j in inbrs:
            jnbrs = set(G[j]) - {j}

            # orbit-2, orbit-3 are two times of the actual number
            for k in inbrs - {j}:
                if k not in jnbrs:
                    orbit_2 += 1
                else:
                    orbit_3 += 1

            for k in (jnbrs - {i}):
                knbrs = set(G[k]) - {k}

                # orbit-1
                if i not in knbrs:
                    orbit_1 += 1

                # orbit-4, orbit-8
                for l in (knbrs - {i} - {j}):
                    if l not in inbrs and l not in jnbrs and k not in inbrs:
                        orbit_4 += 1
                    # orbit-8 are two times of the actual number)
                    if l in inbrs and l not in jnbrs and k not in inbrs:
                        orbit_8 += 1

                # orbit-5
                for l in (inbrs - {j}):
                    if l not in jnbrs and l not in knbrs and k not in inbrs:
                        orbit_5 += 1

            # # based on orbit 6 (uninduced)
            for k, l in combinations((jnbrs - {i}), 2):
                knbrs = set(G[k]) - {k}
                lnbrs = set(G[l]) - {l}

                if k not in inbrs and l not in inbrs and k not in lnbrs:
                    orbit_6 += 1
                if k not in inbrs and l not in inbrs and k in lnbrs:
                    orbit_9 += 1
                if (k in inbrs and k not in lnbrs and l not in inbrs) or (l in inbrs and l not in knbrs and k not in inbrs):
                    orbit_10 += 1
                # orbit-12 are two times of the actual number)
                if (k in lnbrs and l in inbrs and k not in inbrs) or (k in lnbrs and l not in inbrs and k in inbrs):
                    orbit_12 += 1
                # orbit #13
                if k in inbrs and l in inbrs and k not in lnbrs:
                    orbit_13 += 1
                # orbit-14 are three times of the actual number)
                if k in inbrs and l in inbrs and k in lnbrs:
                    orbit_14 += 1

        # orbit-7, orbit-11
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}

            if (u not in v_nbrs) and (u not in w_nbrs) and (v not in w_nbrs):
                orbit_7 += 1

            if (w in u_nbrs) and (v not in w_nbrs) and (v not in u_nbrs):
                orbit_11 += 1
            if (v in u_nbrs) and (w not in v_nbrs) and (w not in u_nbrs):
                orbit_11 += 1
            if (w in v_nbrs) and (u not in w_nbrs) and (u not in v_nbrs):
                orbit_11 += 1

        vec = [orbit_0, orbit_1, orbit_2//2, orbit_3//2, orbit_4, orbit_5, orbit_6, orbit_7, orbit_8//2,
               orbit_9, orbit_10, orbit_11, orbit_12//2, orbit_13, orbit_14//3]
        res[i] = vec
    return res

# extended version: taking into account all 15 orbits
def typed_edge_induced_graphlet_degree_vector(G, num_type, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        inbrs = set(i_nbrs) - {i}
        vec = [[0] * num_type for _ in range(15)]  # 15 graphlets by number of types of edges

        for j in inbrs:
            jnbrs = set(G[j]) - {j}
            ij_type = G.get_edge_data(i, j)['edge_type']
            edit_vec(vec, 0, ij_type)

            # orbit-2, orbit-3 are two times of the actual number
            for k in inbrs - {j}:
                ik_type = G.get_edge_data(i, k)['edge_type']
                if k not in jnbrs:
                    edit_vec(vec, 2, ij_type, ik_type)
                else:
                    jk_type = int(G.get_edge_data(j, k)['edge_type'])
                    edit_vec(vec, 3, ij_type, ik_type, jk_type)

            for k in (jnbrs - {i}):
                knbrs = set(G[k]) - {k}
                jk_type = G.get_edge_data(j, k)['edge_type']

                # orbit-1
                if i not in knbrs:
                    edit_vec(vec, 1, ij_type, jk_type)

                # orbit-4, orbit-8
                for l in (knbrs - {i} - {j}):
                    kl_type = G.get_edge_data(l, k)['edge_type']
                    if l not in inbrs and l not in jnbrs and k not in inbrs:
                        edit_vec(vec, 4, ij_type, jk_type, kl_type)
                    # orbit-8 are two times of the actual number)
                    if l in inbrs and l not in jnbrs and k not in inbrs:
                        il_type = int(G.get_edge_data(i, l)['edge_type'])
                        edit_vec(vec, 8, ij_type, jk_type, kl_type, il_type)

                # orbit-5
                for l in (inbrs - {j}):
                    il_type = G.get_edge_data(l, i)['edge_type']
                    if l not in jnbrs and l not in knbrs and k not in inbrs:
                        edit_vec(vec, 5, ij_type, jk_type, il_type)

            # # based on orbit 6 (uninduced)
            for k, l in combinations((jnbrs - {i}), 2):
                knbrs = set(G[k]) - {k}
                lnbrs = set(G[l]) - {l}
                jk_type = G.get_edge_data(j, k)['edge_type']
                jl_type = G.get_edge_data(j, l)['edge_type']

                if k not in inbrs and l not in inbrs and k not in lnbrs:
                    edit_vec(vec, 6, ij_type, jk_type, jl_type)

                if k not in inbrs and l not in inbrs and k in lnbrs:
                    kl_type = G.get_edge_data(l, k)['edge_type']
                    edit_vec(vec, 9, ij_type, jk_type, jl_type, kl_type)

                if k in inbrs and k not in lnbrs and l not in inbrs:
                    ik_type = G.get_edge_data(i, k)['edge_type']
                    edit_vec(vec, 10, ij_type, jk_type, jl_type, ik_type)
                if l in inbrs and l not in knbrs and k not in inbrs:
                    il_type = G.get_edge_data(i, l)['edge_type']
                    edit_vec(vec, 10, ij_type, jk_type, jl_type, il_type)

                # orbit-12 are two times of the actual number)
                if k in lnbrs:
                    kl_type = G.get_edge_data(l, k)['edge_type']
                    if l in inbrs and k not in inbrs:
                        il_type = G.get_edge_data(i, l)['edge_type']
                        edit_vec(vec, 12, ij_type, jk_type, jl_type, kl_type, il_type)
                    if l not in inbrs and k in inbrs:
                        ik_type = G.get_edge_data(i, k)['edge_type']
                        edit_vec(vec, 12, ij_type, jk_type, jl_type, kl_type, ik_type)

                if k in inbrs and l in inbrs and k not in lnbrs:
                    il_type = G.get_edge_data(i, l)['edge_type']
                    ik_type = G.get_edge_data(i, k)['edge_type']
                    edit_vec(vec, 13, ij_type, jk_type, jl_type, il_type, ik_type)
                # orbit-14 are three times of the actual number)
                # if k in inbrs and l in inbrs and k in lnbrs:
                #     il_type = G.get_edge_data(i, l)['edge_type']
                #     ik_type = G.get_edge_data(i, k)['edge_type']
                #     kl_type = G.get_edge_data(l, k)['edge_type']
                #     edit_vec(vec, 14, ij_type, jk_type, jl_type, kl_type, il_type, ik_type)

        # orbit-7, orbit-11
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']
            iw_type = G.get_edge_data(i, w)['edge_type']

            if (u not in v_nbrs) and (u not in w_nbrs) and (v not in w_nbrs):
                edit_vec(vec, 7, iu_type, iv_type, iw_type)

            if (w in u_nbrs) and (v not in w_nbrs) and (v not in u_nbrs):
                wu_type = G.get_edge_data(w, u)['edge_type']
                edit_vec(vec, 11, iu_type, iv_type, iw_type, wu_type)
            if (v in u_nbrs) and (w not in v_nbrs) and (w not in u_nbrs):
                vu_type = G.get_edge_data(v, u)['edge_type']
                edit_vec(vec, 11, iu_type, iv_type, iw_type, vu_type)
            if (w in v_nbrs) and (u not in w_nbrs) and (u not in v_nbrs):
                wv_type = G.get_edge_data(w, v)['edge_type']
                edit_vec(vec, 11, iu_type, iv_type, iw_type, wv_type)

            if (w in v_nbrs) and (u in w_nbrs) and (u in v_nbrs):
                wv_type = int(G.get_edge_data(w, v)['edge_type'])
                wu_type = int(G.get_edge_data(w, u)['edge_type'])
                vu_type = int(G.get_edge_data(v, u)['edge_type'])
                edit_vec(vec, 14, iu_type, iv_type, iw_type, wv_type, wu_type, vu_type)

        # deal with duplicated count:
        for x in range(len(vec)):
            if x == 2 or x == 3 or x == 8 or x == 12:
                for y in range(len(vec[x])):
                    vec[x][y] //= 2

        res[i] = vec
    return res



# TyE-GDV : GDV with edge type information, returning a 7*num_type 2D list;
# The first element graphlet 2-clique is the classic degree vector with edge type information;
def typed_edge_graphlet_degree_vector_ego(G, num_type, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        vec = [[0] * num_type for _ in range(7)]  # 7 ego-centric graphlets by number of types of edges
        inbrs = set(i_nbrs) - {i}

        # 2-clique, this is also the classic node degree
        for u in inbrs:
            iu_type = G.get_edge_data(i, u)['edge_type']
            edit_vec(vec, 0, iu_type)

        # 3-node graphlets
        for u, v in combinations(inbrs, 2):
            u_nbrs = set(G[u]) - {u}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']

            # 2-path
            edit_vec(vec, 1, iu_type, iv_type)

            # 3-clique
            if v in u_nbrs:
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 2, iu_type, iv_type, uv_type)

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
                uw_type = G.get_edge_data(u, w)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, uw_type)
            if v in w_nbrs:
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, vw_type)
            if v in u_nbrs:
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, uv_type)

            # 4-cycle+
            if w in u_nbrs & v_nbrs:
                uw_type = G.get_edge_data(u, w)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uw_type, vw_type)
            if u in w_nbrs & v_nbrs:
                uw_type = G.get_edge_data(u, w)['edge_type']
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uw_type, uv_type)
            if v in u_nbrs & w_nbrs:
                uv_type = G.get_edge_data(u, v)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uv_type, vw_type)

            # 4-clique
            if (w in u_nbrs & v_nbrs) and (u in v_nbrs):
                uw_type = G.get_edge_data(u, w)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 6, iu_type, iv_type, iw_type, uw_type, vw_type, uv_type)

        res[i] = vec
    return res


# modify the above function with induced graphlet (instead of partial graphlets)
def typed_edge_induced_graphlet_degree_vector_ego(G, num_type, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        vec = [[0] * num_type for _ in range(7)]  # 7 ego-centric graphlets by number of types of edges
        inbrs = set(i_nbrs) - {i}

        # 2-clique (node degree)
        for u in inbrs:
            iu_type = G.get_edge_data(i, u)['edge_type']
            edit_vec(vec, 0, iu_type)

        # 3-node graphlets
        for u, v in combinations(inbrs, 2):
            u_nbrs = set(G[u]) - {u}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']

            # 2-path
            if v not in u_nbrs:
                edit_vec(vec, 1, iu_type, iv_type)
            # 3-clique
            else:
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 2, iu_type, iv_type, uv_type)

        # 4-node graphlets
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}
            iu_type = G.get_edge_data(i, u)['edge_type']
            iv_type = G.get_edge_data(i, v)['edge_type']
            iw_type = G.get_edge_data(i, w)['edge_type']

            # 3-star
            if (u not in v_nbrs) and (u not in w_nbrs) and (v not in w_nbrs):
                edit_vec(vec, 3, iu_type, iv_type, iw_type)

            # tailed-tri
            if (w in u_nbrs) and (v not in w_nbrs) and (v not in u_nbrs):
                uw_type = G.get_edge_data(u, w)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, uw_type)
            if (w in v_nbrs) and (u not in w_nbrs) and (u not in v_nbrs):
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, vw_type)
            if (v in u_nbrs) and (w not in v_nbrs) and (w not in u_nbrs):
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 4, iu_type, iv_type, iw_type, uv_type)

            # 4-cycle+
            if (w in u_nbrs & v_nbrs) and (u not in v_nbrs):
                uw_type = G.get_edge_data(u, w)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uw_type, vw_type)
            if (u in w_nbrs & v_nbrs) and (w not in v_nbrs):
                uw_type = G.get_edge_data(u, w)['edge_type']
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uw_type, uv_type)
            if (v in u_nbrs & w_nbrs) and (u not in w_nbrs):
                uv_type = G.get_edge_data(u, v)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                edit_vec(vec, 5, iu_type, iv_type, iw_type, uv_type, vw_type)

            # 4-clique
            if (w in u_nbrs & v_nbrs) and (u in v_nbrs):
                uw_type = G.get_edge_data(u, w)['edge_type']
                vw_type = G.get_edge_data(v, w)['edge_type']
                uv_type = G.get_edge_data(u, v)['edge_type']
                edit_vec(vec, 6, iu_type, iv_type, iw_type, uw_type, vw_type, uv_type)

        res[i] = vec
    return res


# edge type is coded from 1, edges without type information is coded as 0.
def edit_vec(vec, idx, link_1, link_2=None, link_3=None, link_4=None, link_5=None, link_6=None):
    if int(link_1) > 0:
        vec[idx][int(link_1) - 1] += 1
    if link_2 != None and int(link_2) > 0:
        vec[idx][int(link_2) - 1] += 1
    if link_3 != None and int(link_3) > 0:
        vec[idx][int(link_3) - 1] += 1
    if link_4 != None and int(link_4) > 0:
        vec[idx][int(link_4) - 1] += 1
    if link_5 != None and int(link_5) > 0:
        vec[idx][int(link_5) - 1] += 1
    if link_6 != None and int(link_6) > 0:
        vec[idx][int(link_6) - 1] += 1

def edit_1dvec(vec, link_1, link_2=None, link_3=None, link_4=None, link_5=None, link_6=None):
    if int(link_1) > 0:
        vec[int(link_1) - 1] += 1
    if link_2 != None and int(link_2) > 0:
        vec[int(link_2) - 1] += 1
    if link_3 != None and int(link_3) > 0:
        vec[int(link_3) - 1] += 1
    if link_4 != None and int(link_4) > 0:
        vec[int(link_4) - 1] += 1
    if link_5 != None and int(link_5) > 0:
        vec[int(link_5) - 1] += 1
    if link_6 != None and int(link_6) > 0:
        vec[int(link_6) - 1] += 1


# Ego-GDV = [2-clique, 2-path, 3-clique, 3-star, tailed-tri, 4-chordal-cycle, 4-clique] (partial graphlet degree)
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
    return res


# The original definition of graphlet is induced graph: all edges between nodes of the subgraph considered
def induced_graphlet_degree_vector_ego(G, nodes=None):
    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        inbrs = set(i_nbrs) - {i}
        k = len(inbrs)

        # 3-node graphlets
        two_path = 0
        triangle = 0
        for u, v in combinations(inbrs, 2):
            u_nbrs = set(G[u]) - {u}
            if v not in u_nbrs:
                two_path += 1
            else:
                triangle += 1

        # 4-node graphlets
        three_star = 0
        tailed_tri = 0
        four_cycle_plus = 0
        four_clique = 0
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}

            if (u not in v_nbrs) and (u not in w_nbrs) and (v not in w_nbrs):
                three_star += 1

            if (w in u_nbrs) and (v not in w_nbrs) and (v not in u_nbrs):
                tailed_tri += 1
            if (v in u_nbrs) and (w not in v_nbrs) and (w not in u_nbrs):
                tailed_tri += 1
            if (w in v_nbrs) and (u not in w_nbrs) and (u not in v_nbrs):
                tailed_tri += 1

            if (w in u_nbrs & v_nbrs) and (u not in v_nbrs):
                four_cycle_plus += 1
            if (u in w_nbrs & v_nbrs) and (w not in v_nbrs):
                four_cycle_plus += 1
            if (v in u_nbrs & w_nbrs) and (u not in w_nbrs):
                four_cycle_plus += 1

            if (w in u_nbrs & v_nbrs) and (u in v_nbrs):
                four_clique += 1

        vec = [k, two_path, triangle, three_star, tailed_tri, four_cycle_plus, four_clique]
        res[i] = vec
    return res


# 3-wedge is 3-clique + 1 link (tailed-tri)
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


# 4-cycle+ is 4-cycle with one diagonal link (4-chordal-cycle)
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


# helper method for generating all possible combinations
def get_all_comb_from_list(type_list, num_edge):
    res = []
    for i in range(1, min(num_edge, len(type_list)) + 1):
        res = res + list(combinations(type_list, i))
    return res

# Baseline implementation: colored graphlet approach for typed edge (2-4 node graphlets, 1 - 6 edges)
def colored_graphlet_vector_for_typed_edge(G, num_type, nodes=None):
    type_list = list(range(1,num_type+1))
    comb_1_edge = get_all_comb_from_list(type_list, 1)
    comb_2_edge = get_all_comb_from_list(type_list, 2)
    comb_3_edge = get_all_comb_from_list(type_list, 3)
    comb_4_edge = get_all_comb_from_list(type_list, 4)
    comb_5_edge = get_all_comb_from_list(type_list, 5)
    comb_6_edge = get_all_comb_from_list(type_list, 6)

    if nodes is None:
        nodes_nbrs = G.adj.items()
    else:
        nodes_nbrs = ((n, G[n]) for n in G.nbunch_iter(nodes))

    res = {}
    for i, i_nbrs in nodes_nbrs:
        inbrs = set(i_nbrs) - {i}
        # initialise vec0 to vec14, each representing an orbit
        vec_0 = [0] * len(get_all_comb_from_list(type_list, 1))
        vec_1 = [0] * len(get_all_comb_from_list(type_list, 2))
        vec_2 = [0] * len(get_all_comb_from_list(type_list, 2))
        vec_3 = [0] * len(get_all_comb_from_list(type_list, 3))
        vec_4 = [0] * len(get_all_comb_from_list(type_list, 3))
        vec_5 = [0] * len(get_all_comb_from_list(type_list, 3))
        vec_6 = [0] * len(get_all_comb_from_list(type_list, 3))
        vec_7 = [0] * len(get_all_comb_from_list(type_list, 3))
        vec_8 = [0] * len(get_all_comb_from_list(type_list, 4))
        vec_9 = [0] * len(get_all_comb_from_list(type_list, 4))
        vec_10 = [0] * len(get_all_comb_from_list(type_list, 4))
        vec_11 = [0] * len(get_all_comb_from_list(type_list, 4))
        vec_12 = [0] * len(get_all_comb_from_list(type_list, 5))
        vec_13 = [0] * len(get_all_comb_from_list(type_list, 5))
        vec_14 = [0] * len(get_all_comb_from_list(type_list, 6))

        for j in inbrs:
            jnbrs = set(G[j]) - {j}
            ij_type = int(G.get_edge_data(i, j)['edge_type'])
            t0 = (ij_type, )
            vec_0[comb_1_edge.index(t0)] += 1

            # orbit-2, orbit-3 are two times of the actual number
            for k in inbrs - {j}:
                ik_type = int(G.get_edge_data(i, k)['edge_type'])
                if k not in jnbrs:
                    t2 = tuple(sorted(set((ij_type, ik_type))))
                    vec_2[comb_2_edge.index(t2)] += 1
                else:
                    jk_type = int(G.get_edge_data(j, k)['edge_type'])
                    t3 = tuple(sorted(set((ij_type, ik_type, jk_type))))
                    vec_3[comb_3_edge.index(t3)] += 1

            for k in (jnbrs - {i}):
                knbrs = set(G[k]) - {k}
                jk_type = int(G.get_edge_data(j, k)['edge_type'])
                t2 = tuple(sorted(set((ij_type, jk_type))))
                # orbit-1
                if i not in knbrs:
                    vec_1[comb_2_edge.index(t2)] += 1

                # orbit-4, orbit-8
                for l in (knbrs - {i} - {j}):
                    kl_type = int(G.get_edge_data(l, k)['edge_type'])
                    if l not in inbrs and l not in jnbrs and k not in inbrs:
                        t3 = tuple(sorted(set((ij_type, jk_type, kl_type))))
                        vec_4[comb_3_edge.index(t3)] += 1
                    # orbit-8 are two times of the actual number)
                    if l in inbrs and l not in jnbrs and k not in inbrs:
                        il_type = int(G.get_edge_data(i, l)['edge_type'])
                        t4 = tuple(sorted(set((ij_type, jk_type, kl_type, il_type))))
                        vec_8[comb_4_edge.index(t4)] += 1

                # orbit-5
                for l in (inbrs - {j}):
                    il_type = int(G.get_edge_data(l, i)['edge_type'])
                    t3 = tuple(sorted(set((ij_type, jk_type, il_type))))
                    if l not in jnbrs and l not in knbrs and k not in inbrs:
                        vec_5[comb_3_edge.index(t3)] += 1

            # # based on orbit 6 (uninduced)
            for k, l in combinations((jnbrs - {i}), 2):
                knbrs = set(G[k]) - {k}
                lnbrs = set(G[l]) - {l}
                jk_type = int(G.get_edge_data(j, k)['edge_type'])
                jl_type = int(G.get_edge_data(j, l)['edge_type'])

                if k not in inbrs and l not in inbrs and k not in lnbrs:
                    t3 = tuple(sorted(set((ij_type, jk_type, jl_type))))
                    vec_6[comb_3_edge.index(t3)] += 1

                if k not in inbrs and l not in inbrs and k in lnbrs:
                    kl_type = int(G.get_edge_data(l, k)['edge_type'])
                    t4 = tuple(sorted(set((ij_type, jk_type, jl_type, kl_type))))
                    vec_9[comb_4_edge.index(t4)] += 1

                if k in inbrs and k not in lnbrs and l not in inbrs:
                    ik_type = int(G.get_edge_data(i, k)['edge_type'])
                    t4 = tuple(sorted(set((ij_type, jk_type, jl_type, ik_type))))
                    vec_10[comb_4_edge.index(t4)] += 1
                if l in inbrs and l not in knbrs and k not in inbrs:
                    il_type = int(G.get_edge_data(i, l)['edge_type'])
                    t4 = tuple(sorted(set((ij_type, jk_type, jl_type, il_type))))
                    vec_10[comb_4_edge.index(t4)] += 1

                # orbit-12 are two times of the actual number)
                if k in lnbrs:
                    kl_type = int(G.get_edge_data(l, k)['edge_type'])
                    if l in inbrs and k not in inbrs:
                        il_type = int(G.get_edge_data(i, l)['edge_type'])
                        t5 = tuple(sorted(set((ij_type, jk_type, jl_type, kl_type, il_type))))
                        vec_12[comb_5_edge.index(t5)] += 1
                    if l not in inbrs and k in inbrs:
                        ik_type = int(G.get_edge_data(i, k)['edge_type'])
                        t5 = tuple(sorted(set((ij_type, jk_type, jl_type, kl_type, ik_type))))
                        vec_12[comb_5_edge.index(t5)] += 1

                if k in inbrs and l in inbrs and k not in lnbrs:
                    il_type = int(G.get_edge_data(i, l)['edge_type'])
                    ik_type = int(G.get_edge_data(i, k)['edge_type'])
                    t5 = tuple(sorted(set((ij_type, jk_type, jl_type, il_type, ik_type))))
                    vec_13[comb_5_edge.index(t5)] += 1

        # orbit-7, orbit-11, orbit-14
        for u, v, w in combinations(inbrs, 3):
            u_nbrs = set(G[u]) - {u}
            v_nbrs = set(G[v]) - {v}
            w_nbrs = set(G[w]) - {w}
            iu_type = int(G.get_edge_data(i, u)['edge_type'])
            iv_type = int(G.get_edge_data(i, v)['edge_type'])
            iw_type = int(G.get_edge_data(i, w)['edge_type'])

            if (u not in v_nbrs) and (u not in w_nbrs) and (v not in w_nbrs):
                t3 = tuple(sorted(set((iu_type, iv_type, iw_type))))
                vec_7[comb_3_edge.index(t3)] += 1

            if (w in u_nbrs) and (v not in w_nbrs) and (v not in u_nbrs):
                wu_type = int(G.get_edge_data(w, u)['edge_type'])
                t4 = tuple(sorted(set((iu_type, iv_type, iw_type, wu_type))))
                vec_11[comb_4_edge.index(t4)] += 1
            if (v in u_nbrs) and (w not in v_nbrs) and (w not in u_nbrs):
                vu_type = int(G.get_edge_data(v, u)['edge_type'])
                t4 = tuple(sorted(set((iu_type, iv_type, iw_type, vu_type))))
                vec_11[comb_4_edge.index(t4)] += 1
            if (w in v_nbrs) and (u not in w_nbrs) and (u not in v_nbrs):
                wv_type = int(G.get_edge_data(w, v)['edge_type'])
                t4 = tuple(sorted(set((iu_type, iv_type, iw_type, wv_type))))
                vec_11[comb_4_edge.index(t4)] += 1

            if (w in v_nbrs) and (u in w_nbrs) and (u in v_nbrs):
                wv_type = int(G.get_edge_data(w, v)['edge_type'])
                wu_type = int(G.get_edge_data(w, u)['edge_type'])
                vu_type = int(G.get_edge_data(v, u)['edge_type'])
                t6 = tuple(sorted(set((iu_type, iv_type, iw_type, wv_type, wu_type, vu_type))))
                vec_14[comb_6_edge.index(t6)] += 1

        # deal with duplicated count:
        vec_2 = [i // 2 for i in vec_2]
        vec_3 = [i // 2 for i in vec_3]
        vec_8 = [i // 2 for i in vec_8]
        vec_12 = [i // 2 for i in vec_12]
        vec = vec_0 + vec_1 + vec_2 + vec_3 + vec_4 + vec_5 + vec_6 + vec_7 + vec_8 + vec_9 + vec_10 + vec_11 + vec_12 + vec_13 + vec_14
        res[i] = vec
    return res



