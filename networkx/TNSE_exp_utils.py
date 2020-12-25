import networkx as nx
import pandas as pd

__all__ = ['degree_coefs_corr', 'get_four_coefs', 'closeness_coefs_corr', 'get_feature_score', 'get_iquad_wiquad_df', 'get_oquad_woquad_df']


# if the last chunk is less than bin width, add it to the second to last bin
def chunks(lst, n):
    """Yield successive n-sized chunks from lst."""
    for i in range(0, len(lst), n):
        if len(lst[i : i + n]) == n:
            if len(lst[i + n : i + 2 * n]) < n:
                yield lst[i : i + 2 * n]
            else:
                yield lst[i : i + n]


# for correlation of node degree and four coefs.
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
        avg_iquad, avg_oquad = nx.average_inner_and_outer_quad_co(G, node_list)

        x.append(avg_degree)
        clu.append(avg_clustering)
        clo.append(avg_closure)
        iquad.append(avg_iquad)
        oquad.append(avg_oquad)
    return x, clu, clo, iquad, oquad


def closeness_coefs_corr(G):
    closeness_list = list(nx.closeness_centrality(G).items())
    sorted_list = [(k,v) for k, v in sorted(closeness_list, key=lambda item: item[1])]
    bin_width = len(sorted_list) // 10
    bin_it = chunks(sorted_list, bin_width)
    x = []
    clu = []
    clo = []
    iquad = []
    oquad = []
    for one_bin in bin_it:
        node_list = [n for n, _ in one_bin]
        avg_closeness = sum(d for _, d in one_bin) / len(one_bin)
        avg_clustering = nx.average_clustering(G, node_list)
        avg_closure = nx.average_closure(G, node_list)
        avg_iquad, avg_oquad = nx.average_inner_and_outer_quad_co(G, node_list)

        x.append(avg_closeness)
        clu.append(avg_clustering)
        clo.append(avg_closure)
        iquad.append(avg_iquad)
        oquad.append(avg_oquad)
    return x, clu, clo, iquad, oquad


def get_four_coefs(G):
    clu = nx.clustering(G).values()
    clo = nx.closure(G)
    iquad_oquad = nx.iquad_oquad_coefs(G)
    clo_list = []
    iquad_list = []
    oquad_list = []
    for v in clo.values():
        clo_list.append(v[0])
    for v in iquad_oquad.values():
        iquad_list.append(v[0])
        oquad_list.append(v[1])
    return list(clu), clo_list, iquad_list, oquad_list


def get_feature_score(filename):
    with open(filename,"r") as fi:
        res_string = []
        for ln in fi:
            if ln.startswith("0"):
                fields = ln.split(" ")
                for item in fields:
                    res_string.append(item)
    res = []
    for i in res_string[:7]:
        res.append(float(i))
    return res


def get_iquad_wiquad_df(G):
    iquad = nx.inner_quadrangle_coefficient(G, weight = None)
    wiquad = nx.inner_quadrangle_coefficient(G, weight = 'weight')
    node_iquad_wiquad = []
    for k, v1, v2 in common_entries(iquad, wiquad):
        node_iquad_wiquad.append([k, v1, v2])
    node_iquad_wiquad = sorted(node_iquad_wiquad, key=lambda t:t[1])
    node_iquad_wiquad = pd.DataFrame(node_iquad_wiquad, columns=['node-id','i-quad', 'weighted-i-quad'])
    return node_iquad_wiquad


def get_oquad_woquad_df(G):
    oquad = nx.outer_quadrangle_coefficient(G, weight = None)
    woquad = nx.outer_quadrangle_coefficient(G, weight = 'weight')
    node_oquad_woquad = []
    for k, v1, v2 in common_entries(oquad, woquad):
        node_oquad_woquad.append([k, v1, v2])
    node_oquad_woquad = sorted(node_oquad_woquad, key=lambda t:t[1])
    node_oquad_woquad = pd.DataFrame(node_oquad_woquad, columns=['node-id','o-quad', 'weighted-o-quad'])
    return node_oquad_woquad


def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)