import networkx as nx
import pandas as pd

__all__ = ['average_normalized_patterns_app', 'get_key_info', 'get_network_info', 'get_cc_ce_df']


# get four coefs.
def get_network_info(G, filename="", weight=None):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = 2 * m / n
    clusering = nx.average_clustering(G, weight=weight)
    closure = nx.average_closure(G, weight=weight)
    iquad, oquad = nx.average_inner_and_outer_quad_co(G, weight=weight)
    if filename:
        with open(filename, 'w') as f:
            f.write("|V|:      %d\n" % n)
            f.write("|E|:      %d\n" % m)
            f.write("k:        %.2f\n" % k)
            f.write("clustering: %.3f\n" % clusering)
            f.write("closure:    %.3f\n" % closure)
            f.write("i-quad:     %.3f\n" % iquad)
            f.write("o-quad:     %.3f\n" % oquad)
            if closure > 0:
                f.write("clustering/closure: %.3f\n" % (clusering/closure))
            if oquad > 0:
                f.write("i-quad/o-quad:      %.3f\n" % (iquad/oquad))
            if clusering > 0:
                f.write("i-quad/clustering:  %.3f\n" % (iquad/clusering))
            if closure > 0:
                f.write("o-quad/closure:     %.3f\n" % (oquad/closure))
    res = [n, m, k, clusering, closure, iquad, oquad]
    signatures = []
    signatures.append(clusering / closure if closure > 0 else 0)
    signatures.append(iquad / oquad if oquad > 0 else 0)
    signatures.append(iquad / clusering if clusering > 0 else 0)
    signatures.append(oquad / closure if closure > 0 else 0)
    res.append(signatures)
    return res


def get_key_info(G, weight = None):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = m / n
    r = nx.overall_reciprocity(G)
    cc = nx.average_clustering(G, weight=weight)
    ce = nx.average_closure(G, weight=weight)

    res = [n, m, k, r, cc, ce]
    return res


def average_normalized_patterns_app(G, nodes = None):
    res = []
    dict_head = nx.head_closure(G, nodes)
    dict_mid = nx.mid_closure(G, nodes)
    dict_end = nx.end_closure(G, nodes)
    dict_cyc = nx.cyc_closure(G, nodes)

    dict_all_pattern = dict()
    for (key, h, m, e, c) in zip(dict_head.keys(), dict_head.values(), dict_mid.values(), dict_end.values(),
                                 dict_cyc.values()):
        dict_all_pattern[key] = (h + m + e + c)

    normalized_head = dict()
    normalized_mid = dict()
    normalized_end = dict()
    normalized_cyc = dict()
    for (k, head, all_pattern) in zip(dict_head.keys(), dict_head.values(), dict_all_pattern.values()):
        if all_pattern > 0:
            if head == 0:
                normalized_head[k] = 0
            else:
                normalized_head[k] = head / all_pattern

    for (k, mid, all_pattern) in zip(dict_mid.keys(), dict_mid.values(), dict_all_pattern.values()):
        if all_pattern > 0:
            if mid == 0:
                normalized_mid[k] = 0
            else:
                normalized_mid[k] = mid / all_pattern

    for (k, end, all_pattern) in zip(dict_end.keys(), dict_end.values(), dict_all_pattern.values()):
        if all_pattern > 0:
            if end == 0:
                normalized_end[k] = 0
            else:
                normalized_end[k] = end / all_pattern

    for (k, cyc, all_pattern) in zip(dict_cyc.keys(), dict_cyc.values(), dict_all_pattern.values()):
        if all_pattern > 0:
            if cyc == 0:
                normalized_cyc[k] = 0
            else:
                normalized_cyc[k] = cyc / all_pattern

    length = len(normalized_head)
    res.append(sum(normalized_head.values()) / length)
    res.append(sum(normalized_mid.values()) / length)
    res.append(sum(normalized_end.values()) / length)
    res.append(sum(normalized_cyc.values()) / length)

    return res

def get_cc_ce_df(G, weight = None):
    cc = nx.clustering(G, weight = weight)
    ce = nx.closure(G, weight = weight)
    node_cc_ce = []
    for k, v1, v2 in common_entries(cc, ce):
        node_cc_ce.append([k, v1, v2[0]])
    node_cc_ce = sorted(node_cc_ce, key=lambda t:t[1])
    df_cc_ce = pd.DataFrame(node_cc_ce, columns=['node-id','cc', 'ce'])
    return df_cc_ce


def common_entries(*dcts):
    for i in set(dcts[0]).intersection(*dcts[1:]):
        yield (i,) + tuple(d[i] for d in dcts)