import networkx as nx

__all__ = ['degree_coefs_corr']

# if the last chunk is less than bin width, add it to the second to last bin
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        if len(lst[i : i + n]) == n:
            if len(lst[i + n : i + 2 * n]) < n:
                yield lst[i : i + 2 * n]
            else:
                yield lst[i : i + n]


def degree_coefs_corr(G):
    degree_list = list(G.degree())
    sorted_degree_list = [(k,v) for k, v in sorted(degree_list, key=lambda item: item[1])]
    bin_width = len(sorted_degree_list) // 10
    bin_it = chunks(sorted_degree_list, bin_width)
    x = []
    clu = []
    clo = []
    iquad = []
    oquad = []
    for one_bin in bin_it:
        node_list = [n for n, _ in one_bin]
        avg_degree = sum(d for _, d in one_bin) / len(one_bin)
        avg_clustering = nx.average_clustering(G, node_list)
        avg_closure = nx.average_closure(G, node_list)
        avg_iquad = nx.average_inner_quad_co(G, node_list)
        avg_oquad = nx.average_outer_quad_co(G, node_list)

        if avg_degree > 1:
            x.append(avg_degree)
        if avg_clustering > 0:
            clu.append(avg_clustering)
        if avg_closure > 0:
            clo.append(avg_closure)
        if avg_iquad > 0:
            iquad.append(avg_iquad)
        if avg_oquad > 0:
            oquad.append(avg_oquad)
    return x, clu, clo, iquad, oquad