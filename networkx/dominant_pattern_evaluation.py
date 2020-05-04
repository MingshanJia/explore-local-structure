import networkx as nx
__all__ = ['average_normalized_patterns_app', 'get_key_info']


def get_key_info(G, weight = None):
    n = G.number_of_nodes()
    m = G.number_of_edges()
    k = m / n
    r = nx.overall_reciprocity(G)
    cc = nx.average_clustering(G, weight)
    ce = nx.average_closure(G, weight)

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