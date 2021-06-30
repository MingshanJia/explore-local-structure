

__all__ = ['edge_clustering_coef', 'edge_clustering_vector', 'average_edge_clustering']

def _edge_clustering_iter(G, edges=None):
    if edges is None:
        edges = G.edges()
    for e in edges:
        tri = len(set(G[e[0]]) & set(G[e[1]]))
        poss_tri = min(G.degree(e[0]), G.degree(e[1])) - 1
        yield (e, tri, poss_tri)


def edge_clustering_coef(G, edges=None):
    '''Wang, Jianxin, et al. "Identification of essential proteins based on edge clustering coefficient.'''
    ec_iter = _edge_clustering_iter(G, edges)
    ecc = {e: 0 if t == 0 else t / pt for e, t, pt in ec_iter}
    if edges in G:
        return ecc[edges]
    return ecc


def average_edge_clustering(G, edges=None, count_zeros=True):
    ecc = edge_clustering_coef(G, edges).values()
    if not count_zeros:
        ecc = [v for v in ecc if v > 0]
    return sum(ecc) / len(ecc)


# ECV contains 13 coordinates: each coordinate indicate one type of relationship
# None in return value means that there is no such relationship at all.
def edge_clustering_vector(G, edges=None):
    vec = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    type_list = [[],[],[],[],[],[],[],[],[],[],[],[],[]]
    for u,v,data in G.edges(data=True):
        if data.get('edge_type') == '1':
            type_list[0].append((u,v))
        if data.get('edge_type') == '2':
            type_list[1].append((u,v))
        if data.get('edge_type') == '3':
            type_list[2].append((u,v))
        if data.get('edge_type') == '4':
            type_list[3].append((u,v))
        if data.get('edge_type') == '5':
            type_list[4].append((u,v))
        if data.get('edge_type') == '6':
            type_list[5].append((u,v))
        if data.get('edge_type') == '7':
            type_list[6].append((u,v))
        if data.get('edge_type') == '8':
            type_list[7].append((u,v))
        if data.get('edge_type') == '9':
            type_list[8].append((u,v))
        if data.get('edge_type') == '10':
            type_list[9].append((u,v))
        if data.get('edge_type') == '11':
            type_list[10].append((u,v))
        if data.get('edge_type') == '12':
            type_list[11].append((u,v))
        if data.get('edge_type') == '13':
            type_list[12].append((u,v))
    print(type_list)
    for i in range(13):
        if len(type_list[i]) != 0:
            vec[i] = average_edge_clustering(G, edges=type_list[i])
        else:
            vec[i] = None
    return vec