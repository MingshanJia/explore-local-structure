import random
import networkx as nx
import numpy as np
from operator import itemgetter
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve

__all__ = ['link_predict_new_metrics', 'link_pred_app', 'undir_link_pred_app']

# metrics: ROC-AUC, PR-AUC, Average Precision
def link_predict_new_metrics(G, repeat = 10, old_pct = 0.8, method = 'cn'):
    roc_auc = 0
    pr_auc = 0
    ave_precision = 0
    all_edges = list(G.edges())
    k = round(len(all_edges) * old_pct)
    print("using {:f} % of edges to predict {:f} % edges".format(100 * old_pct, (100 * (1 - old_pct))))

    if G.number_of_nodes() > 10000:
        sample = True
        sample_time = 10
        print("sample 5K nodes for 10 times")
    else:
        sample = False
        sample_time = 1
        print("no sampling")

    for i in range(0, sample_time):
        if sample:
            print('sample: %d' % i)
            random_edge = random.choice(all_edges)
            random_node = random_edge[0]
            sample_nodes = get_sample_nodes(G, random_node, 5000)
            print("sampling done")
            G_sampled = G.subgraph(sample_nodes)
            all_edges = list(G_sampled.edges(data=True))

        for n in range(0, repeat):
            G_old = nx.Graph()
            G_new = nx.Graph()
            print('repeat time = %d' % n)
            random.shuffle(all_edges)
            old_edges = all_edges[:k]
            new_edges = all_edges[k:]
            G_old.add_edges_from(old_edges)
            G_new.add_edges_from(new_edges)
            G_new = G_new.subgraph(G_old.nodes())

            label_all, score_all = get_predicts_and_labels(G_old, G_new, method)
            precision, recall, _ = precision_recall_curve(label_all, score_all)
            roc_auc += roc_auc_score(label_all, score_all)
            pr_auc += auc(recall, precision)
            ave_precision += average_precision_score(label_all, score_all)

    roc_auc /= (repeat * sample_time)
    pr_auc /= (repeat * sample_time)
    ave_precision /= (repeat * sample_time)
    print("{} :\nROC-AUC: {};\nPR-AUC: {};\nAve_Precision: {}.".format(method, roc_auc, pr_auc, ave_precision))


def get_predicts_and_labels(G_old, G_new, method):
    predicts = get_predict_score(G_old, method)
    pos_label = []
    pos_score = []
    neg_label = []
    neg_score = []
    for p in predicts:
        if ((p[0], p[1]) in G_new.edges()):
            pos_label.append(1)
            pos_score.append(p[2])
        else:
            neg_label.append(0)
            neg_score.append(p[2])

    label_all = np.hstack([pos_label, neg_label])
    score_all = np.hstack([pos_score, neg_score])
    return label_all, score_all


def get_predict_score(G_old, method):
    if method == 'cn':
        predicts = nx.common_neighbor_index(G_old)
    if method == 'ra':
        predicts = nx.resource_allocation_index(G_old)
    if method == 'cn+clu':
        predicts = nx.common_neighbor_plus_clustering(G_old)
    if method == 'cn-l3':
        predicts = nx.common_neighbor_l3(G_old)
    if method == 'cn-l3-norm':
        predicts = nx.common_neighbor_l3_degree_normalized(G_old)

    max_score = max(predicts, key=itemgetter(2))[2]
    normalized_predicts = [(i[0], i[1], i[2] / max_score) for i in predicts]
    return normalized_predicts


def undir_link_pred_app(G, repeat=1, sample_time=10, sample_size=5000, old_pct=0.5):
    res = []
    rg = 0  # random guess
    cn = 0
    ra = 0
    cn_clu = 0
    cn_l3 = 0
    cn_l3_norm = 0

    e_all = list(G.edges(data=True))

    if G.number_of_nodes() > 10000:
        sample = True
        sample_time = sample_time
    else:
        sample = False
        sample_time = 1

    print('sample time : %d' % sample_time)
    print('repeat split time: %d' % repeat)
    print('old edges percentage: %.1f' % old_pct)

    for i in range(0, sample_time):
        print('sample: %d' % i)

        if sample:
            random_edge = random.choice(e_all)
            random_node = random_edge[0]
            sample_nodes = get_sample_nodes(G, random_node, sample_size)
            print("sampling done")
            G_sampled = G.subgraph(sample_nodes)
            e_sampled = list(G_sampled.edges(data=True))
        else:
            e_sampled = e_all

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
            G_old = nx.Graph()
            G_new = nx.Graph()
            G_old.add_edges_from(e_old)
            G_new.add_edges_from(e_new)

            rg += nx.random_guess(G_old, G_new)
            print('rg: %.3f' % rg)
            cn += nx.perform_link_prediction_undir(G_old, G_new, 'cn')
            print('cn: %.3f' % cn)
            ra += nx.perform_link_prediction_undir(G_old, G_new, 'ra')
            print('ra: %.3f' % ra)
            cn_clu += nx.perform_link_prediction_undir(G_old, G_new, 'cn+clu')
            print('cn+clu: %.3f' % cn_clu)
            cn_l3 += nx.perform_link_prediction_undir(G_old, G_new, 'cn-l3')
            print('cn-l3: %.3f' % cn_l3)
            cn_l3_norm += nx.perform_link_prediction_undir(G_old, G_new, 'cn-l3-norm')
            print('cn-l3-norm: %.3f' % cn_l3_norm)

    rg /= (repeat * sample_time)
    cn /= (repeat * sample_time)
    ra /= (repeat * sample_time)
    cn_clu /= (repeat * sample_time)
    cn_l3 /= (repeat * sample_time)
    cn_l3_norm /= (repeat * sample_time)
    res.append(rg)
    res.append(cn)
    res.append(ra)
    res.append(cn_clu)
    res.append(cn_l3)
    res.append(cn_l3_norm)
    return res



# when nodes > 10K, sampling.
# when graph has timestamps, repeat = 1
def link_pred_app(G, repeat=1, sample_time=10, sample_size=5000, old_pct=0.5):
    res = []
    rg = 0 #random guess
    cn = 0
    aa = 0
    ra = 0
    clo1 = 0
    clo2 = 0

    e_all = list(G.edges(data=True))

    if G.number_of_nodes() > 10000:
        sample = True
        sample_time = sample_time
    else:
        sample = False
        sample_time = 1

    print('sample time : %d' % sample_time)
    print('repeat split time: %d' % repeat)
    print('old edges percentage: %.1f' % old_pct)
           
    for i in range(0, sample_time):
        print('sample: %d' % i)

        if sample:
            random_edge = random.choice(e_all)
            random_node = random_edge[0]
            sample_nodes = get_sample_nodes(G, random_node, sample_size)
            print("sampling done")
            G_sampled = G.subgraph(sample_nodes)
            e_sampled = list(G_sampled.edges(data=True))
        else:
            e_sampled = e_all

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