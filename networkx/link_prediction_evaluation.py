import random
import networkx as nx

__all__ = ['link_pred_app',
           'link_pred_sample_app',
           'link_pred_sample_app_2']


# for not very large networks, nodes < 10K.
# input graph G, repeat = 1 for dataset with timestamp, set repeat > 1 for dataset without timestamp
# old_pct is the the percentage of old edges, based on which we predict new edges
def link_pred_app(G, repeat=1, old_pct=0.5):
    print('old edges percentage: %.1f' % old_pct)
    print('repeat time: %d' % repeat)
    res = []
    rg = 0 #random guess
    cn = 0
    aa = 0
    ra = 0
    clo1 = 0
    clo2 = 0
    #dgr = 0
    e_all = list(G.edges(data=True))
    k = round(len(e_all) * old_pct)

    for n in range(0, repeat):

        print('n = %d' % n)

        if repeat == 1:
            e_all = sorted(e_all, key=lambda t: t[2].get('sec'))  # for dataset with timestamp
        else:
            random.shuffle(e_all)

        e_old = e_all[: k]
        e_new = e_all[k:]
        G_old = nx.DiGraph()
        G_new = nx.DiGraph()
        G_old.add_edges_from(e_old)
        G_new.add_edges_from(e_new)

        dict_ce = nx.closure(G_old)

        rg += nx.random_guess(G_old, G_new)
        print('rg: %.3f' % (rg))
        cn += nx.perform_link_prediction(G_old, G_new, 'cn', dict_ce)
        print('cn: %.3f' % (cn))
        aa += nx.perform_link_prediction(G_old, G_new, 'aa', dict_ce)
        print('aa: %.3f' % (aa))
        ra += nx.perform_link_prediction(G_old, G_new, 'ra', dict_ce)
        print('ra: %.3f' % (ra))
        clo1 += nx.perform_link_prediction(G_old, G_new, 'clo1', dict_ce)
        print('clo1: %.3f' % (clo1))
        clo2 += nx.perform_link_prediction(G_old, G_new, 'clo2', dict_ce)
        print('clo2: %.3f' % (clo2))
        #dgr += nx.perform_link_prediction(G_old, G_new, 'dgr', dict_ce)
        #print('dgr: %.3f' % (dgr))

    rg /= (repeat)
    cn /= (repeat)
    aa /= (repeat)
    ra /= (repeat)
    clo1 /= (repeat)
    clo2 /= (repeat)
    #dgr /= repeat

    res.append(rg)
    res.append(cn)
    res.append(aa)
    res.append(ra)
    res.append(clo1)
    res.append(clo2)
    #res.append(dgr)
    return res


# sample 5K as G, then split into G_old and G_new
def link_pred_sample_app_2(G, repeat=1, sample_time=10, sample_size=5000, old_pct=0.5):
    print('sample time : %d' % sample_time)     
    print('repeat split time: %d' % repeat)
    print('old edges percentage: %.1f' % old_pct)
    res = []
    rg = 0 #random guess
    cn = 0
    aa = 0
    ra = 0
    clo1 = 0
    clo2 = 0

    e_all = list(G.edges(data=True))
           
    for i in range(0, sample_time):
        print('sample: %d' % i)       
        random_edge = random.choice(e_all)
        random_node = random_edge[0]
        sample_nodes = get_sample_nodes(G, random_node, sample_size)        
        print("sampling done")
        G_sampled = G.subgraph(sample_nodes)
        e_sampled = list(G_sampled.edges(data=True))
        k = round(len(e_sampled) * old_pct)  
        print('k = %d' % k)    

        for n in range(0, repeat):
            print('repeat = %d' % n)

            if repeat == 1:
                e_sampled = sorted(e_sampled, key=lambda t: t[2].get('sec'))  # for dataset with timestamp
            else:
                random.shuffle(e_sampled)

            e_old = e_sampled[: k]
            e_new = e_sampled[k:]
            G_old = nx.DiGraph()
            G_new = nx.DiGraph()
            G_old.add_edges_from(e_old)
            G_new.add_edges_from(e_new)

            dict_ce = nx.closure(G_old)

            rg += nx.random_guess(G_old, G_new)
            print('rg: %.3f' % (rg))
            cn += nx.perform_link_prediction(G_old, G_new, 'cn', dict_ce)
            print('cn: %.3f' % (cn))
            aa += nx.perform_link_prediction(G_old, G_new, 'aa', dict_ce)
            print('aa: %.3f' % (aa))
            ra += nx.perform_link_prediction(G_old, G_new, 'ra', dict_ce)
            print('ra: %.3f' % (ra))
            clo1 += nx.perform_link_prediction(G_old, G_new, 'clo1', dict_ce)
            print('clo1: %.3f' % (clo1))
            clo2 += nx.perform_link_prediction(G_old, G_new, 'clo2', dict_ce)
            print('clo2: %.3f' % (clo2))

    rg /= (repeat * sample_time)
    cn /= (repeat * sample_time)
    aa /= (repeat * sample_time)
    ra /= (repeat * sample_time)
    clo1 /= (repeat * sample_time)
    clo2 /= (repeat * sample_time)

    res.append(rg)
    res.append(cn)
    res.append(aa)
    res.append(ra)
    res.append(clo1)
    res.append(clo2)
    return res


# for very large networks, nodes > 10K. input graph G, repeat = 1 for dataset with timestamp, set repeat > 1 for dataset without timestamp
# old_pct is the the percentage of old edges, based on which we predict new edges
def link_pred_sample_app(G, repeat=1, sample_time=5, sample_size=1000, old_pct=0.5):
    print('old edges percentage: %.1f' % old_pct)
    print('repeat time: %d' % repeat)
    res = []
    rg = 0 #random guess
    cn = 0
    aa = 0
    ra = 0
    clo1 = 0
    clo2 = 0
    #dgr = 0

    e_all = list(G.edges(data=True))
    k = round(len(e_all) * old_pct)

    for n in range(0, repeat):

        print('n = %d' % n)

        if repeat == 1:
            e_all = sorted(e_all, key=lambda t: t[2].get('sec'))  # for dataset with timestamp
            print("sort done\n")
        else:
            random.shuffle(e_all)

        e_old = e_all[: k]
        e_new = e_all[k:]
        G_old = nx.DiGraph()
        G_new = nx.DiGraph()
        G_old.add_edges_from(e_old)
        G_new.add_edges_from(e_new)

        for i in range(0, sample_time):
            print('sample %d' % i)

            random_edge = random.choice(e_old)
            random_node = random_edge[0]

            sample_nodes = get_sample_nodes(G_old, random_node, sample_size)
            print("sampling done")
            G_old_sampled = G_old.subgraph(sample_nodes)

            dict_ce = nx.closure(G_old_sampled)

            rg += nx.random_guess(G_old_sampled, G_new)
            print('rg: %.3f' % (rg))
            cn += nx.perform_link_prediction(G_old_sampled, G_new, 'cn', dict_ce)
            print('cn: %.3f' % (cn))
            aa += nx.perform_link_prediction(G_old_sampled, G_new, 'aa', dict_ce)
            print('aa: %.3f' % (aa))
            ra += nx.perform_link_prediction(G_old_sampled, G_new, 'ra', dict_ce)
            print('ra: %.3f' % (ra))
            clo1 += nx.perform_link_prediction(G_old_sampled, G_new, 'clo1', dict_ce)
            print('clo1: %.3f' % (clo1))
            clo2 += nx.perform_link_prediction(G_old_sampled, G_new, 'clo2', dict_ce)
            print('clo2: %.3f' % (clo2))
            #dgr += nx.perform_link_prediction(G_old, G_new, 'dgr', dict_ce)
            #print('dgr: %.3f' % (dgr))

    rg /= (repeat * sample_time)
    cn /= (repeat * sample_time)
    aa /= (repeat * sample_time)
    ra /= (repeat * sample_time)
    clo1 /= (repeat * sample_time)
    clo2 /= (repeat * sample_time)
    #dgr /= (repeat * sample_time)

    res.append(rg)
    res.append(cn)
    res.append(aa)
    res.append(ra)
    res.append(clo1)
    res.append(clo2)
    #res.append(dgr)
    return res


# get sample: connected nodes
def get_sample_nodes(G, random_node, sample_size):
    res_list = []  # result list of connected nodes
    work_list = []
    visited_list = []
    res_list.append(random_node)
    work_list.append(random_node)

    while (len(res_list) < sample_size) and (len(work_list) > 0):

        for n in work_list:
            for nbr in nx.all_neighbors(G, n):
                if nbr not in res_list:
                    res_list.append(nbr)
                    # return faster
                    if len(res_list) == sample_size:
                        return res_list          
            visited_list.append(n)

        work_list = [n for n in res_list if n not in visited_list]

    return res_list[: sample_size]


#for networks with timestamps
# method = 'closure' or 'degree'
# worse than pure degree! not meaningful to use..
# def link_direction_prediction_app(G, method = 'closure', repeat = 10, old_pct=0.5):
#
#     e_all = list(G.edges(data=True))
#     e_all = sorted(e_all, key=lambda t: t[2].get('sec'))    #for networks with timestamps
#     k = round(len(e_all) * old_pct)
#
#     e_old = e_all[: k]
#     e_new = e_all[k:]
#     G_old = nx.DiGraph()
#     G_old.add_edges_from(e_old)
#
#     avg_precision = 0
#     for i in range(0, repeat):
#         random.shuffle(e_new)
#         e_new_sample = e_new[: 1000]
#         e_target = get_e_possible_to_predict(e_new_sample, G_old)  # ground truth
#         target_num_links = len(e_target)
#         print("target number of links: %d" % target_num_links)
#
#         G_new = nx.DiGraph()
#         G_new.add_edges_from(e_target)  # ground truth
#         G_new_undirected = G_new.to_undirected()
#
#         if (method == 'closure'):
#             dict_e_with_di_score = get_direction_score(G_new_undirected, G_old)
#         if (method == 'degree'):
#             dict_e_with_di_score = get_direction_score_two(G_new_undirected, G_old)
#         if (method == 'mixed'):
#             dict_e_with_di_score = get_direction_score_mixed(G_new_undirected, G_old)
#
#         dict_e_with_zero = {k:v for k,v in dict_e_with_di_score.items() if v[0]==0}   # dict containg edges with zero direction score
#         ordered_list_of_zero_dict = [(k, v) for k, v in sorted(dict_e_with_zero.items(), key=lambda item: item[1][1], reverse=True)]  #sort according to L2R + R2L in descending order
#         len_zero = len(ordered_list_of_zero_dict)
#
#         res_G_directed = nx.DiGraph()  # what to return
#
#         for k, v in dict_e_with_di_score.items():
#             if v[0] > 0:
#                 res_G_directed.add_edge(k[0], k[1])
#             if v[0] < 0:
#                 res_G_directed.add_edge(k[1], k[0])
#         current_num_link = res_G_directed.number_of_edges()
#
#         index = 0
#         while current_num_link < target_num_links and index < len_zero:
#             e_key = ordered_list_of_zero_dict[index][0]
#             res_G_directed.add_edge(e_key[0], e_key[1])
#             current_num_link += 1
#             if current_num_link == target_num_links:
#                 break
#             res_G_directed.add_edge(e_key[1], e_key[0])
#             current_num_link += 1
#             index += 1
#
#         if current_num_link < target_num_links:
#             dict_abs_score = dict()  # used to sort, only contain edges with non_zero score
#             for key in dict_e_with_di_score:
#                 if dict_e_with_di_score[key][0] != 0:
#                     dict_abs_score[key] = abs(dict_e_with_di_score[key][0])
#
#             ordered_list_of_abs_dict = [(k, v) for k, v in
#                                         sorted(dict_abs_score.items(), key=lambda item: item[1])]  # in ascending order
#
#             length = len(ordered_list_of_abs_dict)
#             index = 0  # initialised at the begining of list
#
#             while current_num_link < target_num_links and index < length:
#                 e_key = ordered_list_of_abs_dict[index][0]
#                 e_value = dict_e_with_di_score[e_key][0]
#
#                 if e_value > 0:
#                     res_G_directed.add_edge(e_key[1], e_key[0])
#                     current_num_link += 1
#                 else:
#                     res_G_directed.add_edge(e_key[0], e_key[1])
#                     current_num_link += 1
#                 index += 1
#
#         numerator = 0
#         for e in res_G_directed.edges():
#
#             if res_G_directed.has_edge(e[0], e[1]) and G_new.has_edge(e[0], e[1]):
#                 numerator += 1
#         avg_precision += numerator / target_num_links
#     return avg_precision / repeat



# def get_e_possible_to_predict(e_new, G_old):
#     e_possible_to_predict = []
#     old_node_list = list(G_old.nodes())
#
#     for e in e_new:
#         if (e[0] in old_node_list) and (e[1] in old_node_list):
#             e_possible_to_predict.append(e)
#     return e_possible_to_predict


# get direction score for interested edges, using src-clo and tgt-clo
# def get_direction_score(G_new_sample, G_old):
#     dict_e_with_di_score = dict()
#
#     dict_ce = nx.closure(G_old, G_new_sample.nodes())   #only for nodes exist in G_new
#
#     for e in G_new_sample.edges():
#         left_src = dict_ce[e[0]][1]
#         left_tgt = dict_ce[e[0]][2]
#         right_src = dict_ce[e[1]][1]
#         right_tgt = dict_ce[e[1]][2]
#         L2R = left_src + right_tgt
#         R2L = right_src + left_tgt
#
#         dict_e_with_di_score[e] = [L2R - R2L, L2R + R2L]
#
#     return dict_e_with_di_score



# get direction score for interested edges, using out_degree and in_degree
# def get_direction_score_two(G_new_sample, G_old):
#     dict_e_with_di_score = dict()
#
#     for e in G_new_sample.edges():
#         left_src = G_old.out_degree(e[0])
#         left_tgt = G_old.in_degree(e[0])
#         right_src = G_old.out_degree(e[1])
#         right_tgt = G_old.in_degree(e[1])
#         L2R = left_src + right_tgt
#         R2L = right_src + left_tgt
#
#         dict_e_with_di_score[e] = [L2R - R2L, L2R + R2L]
#
#     return dict_e_with_di_score


# get direction score for interested edges, using clo and degree
# def get_direction_score_mixed(G_new_sample, G_old):
#     dict_e_with_di_score = dict()
#
#     dict_ce = nx.closure(G_old, G_new_sample.nodes())   #only for nodes exist in G_new
#
#     for e in G_new_sample.edges():
#         left_src = dict_ce[e[0]][1] * G_old.out_degree(e[0])
#         left_tgt = dict_ce[e[0]][2] * G_old.in_degree(e[0])
#         right_src = dict_ce[e[1]][1] * G_old.out_degree(e[1])
#         right_tgt = dict_ce[e[1]][2] * G_old.in_degree(e[1])
#         L2R = left_src + right_tgt
#         R2L = right_src + left_tgt
#
#         dict_e_with_di_score[e] = [L2R - R2L, L2R + R2L]
#
#     return dict_e_with_di_score



