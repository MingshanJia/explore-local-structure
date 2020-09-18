import random
import networkx as nx
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from operator import itemgetter
from sklearn.preprocessing import normalize
from sklearn.metrics import roc_auc_score
from sklearn.metrics import auc
from sklearn.metrics import average_precision_score
from sklearn.metrics import precision_recall_curve
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from xgboost import XGBClassifier

__all__ = ['link_predict_supervised_learning', 'get_dataset',
           'link_predict_similarity_based', 'undir_link_pred_app', 'link_pred_app', ]
# usage: 1. train_set = nx.get_train_set(G)    2. nx.link_predict_supervised_learning(train_set)


#APP 1
def link_predict_supervised_learning(train_set, method='log-reg', number_of_features=2):
    positive_ratio = 0
    roc_auc = 0
    pr_auc = 0
    ave_precision = 0
    feature_importance = np.zeros(round(number_of_features))
    n = len(train_set)
    for g in tqdm(train_set):
        positive_ratio += g[1].number_of_edges() / len(list(nx.non_edges(g[0])))
        label_all, score_all, feature_importance_of_g = get_predicts_labels_and_feature_importance(g[0], g[1],  method, number_of_features)
        precision, recall, _ = precision_recall_curve(label_all, score_all)
        pr_auc += auc(recall, precision)
        roc_auc += roc_auc_score(label_all, score_all)
        ave_precision += average_precision_score(label_all, score_all)
        feature_importance += feature_importance_of_g
    positive_ratio /= n
    roc_auc /= n
    pr_auc /= n
    ave_precision /= n
    feature_importance /= n
    print("{} :\nPositive_Ratio: {}\nROC-AUC: {};\nPR-AUC: {};\nAve_Precision: {}.".format(method, positive_ratio, roc_auc, pr_auc, ave_precision))

    # if number_of_features == 2:
    #     features = ['cn', 'l3']
    # if number_of_features == 4:
    #     features = ['cn', 'ra', 'l3', 'l3-norm']
    # if number_of_features == 6:
    #     features = ['cn', 'ra', 'l3', 'l3-norm', 'clu', 'clo']
    # if number_of_features == 8:
    #     features = ['cn', 'ra', 'l3', 'l3-norm', 'clu', 'clo', 'i-quad', 'o-quad']
    # plt.bar(features, feature_importance)
    # plt.show()
    # for feature, score in zip(features, feature_importance):
    #     print(feature, score)
    return roc_auc, pr_auc, ave_precision, feature_importance



def get_predicts_labels_and_feature_importance(G_old, G_new, method, number_of_features):
    X, y = get_features_and_labels(G_old, G_new, number_of_features)
    X = normalize(X)
    if method =='log-reg':
        model = LogisticRegression()
        model.fit(X, y)
        importances = model.coef_[0]
    if method == 'tree':
        model = DecisionTreeClassifier()
        model.fit(X, y)
        importances = model.feature_importances_
    if method == 'forest':
        model = RandomForestClassifier()
        model.fit(X, y)
        importances = model.feature_importances_
    #TODO: add importances for svm
    if method == 'svm':
        model = svm.SVC(kernel='linear', probability=True)
        model.fit(X, y)
        importances = model.coef_[0]
    if method == 'xgboost':
        model = XGBClassifier()
        X = np.asarray(X)
        model.fit(X, y)
        importances = model.feature_importances_

    score = model.predict_proba(X)
    return y, score[:, 1], importances


# compare with 4 and 6
def get_features_and_labels(G_old, G_new, number_of_features):
    possible_links = list(nx.non_edges(G_old))
    X = []
    y = []
    if number_of_features == 2:
        for u, v in possible_links:
            cn_score = len(list(nx.common_neighbors(G_old, u, v)))
            cn_l3_score = nx.common_neighbors_l3(G_old, u, v)
            X.append([cn_score, cn_l3_score])
            if (u, v) in G_new.edges():
                y.append(1)
            else:
                y.append(0)

    # if number_of_features == 4:
    #     for u, v in possible_links:
    #         cn_score = len(list(nx.common_neighbors(G_old, u, v)))
    #         ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
    #         cn_l3_score = nx.common_neighbors_l3(G_old, u, v)
    #         cn_l3_norm_score = nx.common_neighbors_l3_degree_normalized(G_old, u, v)
    #         X.append([cn_score, ra_score, cn_l3_score, cn_l3_norm_score])
    #         if (u, v) in G_new.edges():
    #             y.append(1)
    #         else:
    #             y.append(0)

    if number_of_features == 4:
        clu_clo_dict = nx.clustering_closure_coefs(G_old)
        for u, v in possible_links:
            cn_score = len(list(nx.common_neighbors(G_old, u, v)))
            ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
            clu_score = clu_clo_dict[u][0] + clu_clo_dict[v][0]
            clo_score = clu_clo_dict[u][1] + clu_clo_dict[v][1]
            X.append([cn_score, ra_score, clu_score, clo_score])
            if (u, v) in G_new.edges():
                y.append(1)
            else:
                y.append(0)

    # if number_of_features == 6:
    #     clu_clo_dict = nx.clustering_closure_coefs(G_old)
    #     for u, v in possible_links:
    #         cn_score = len(list(nx.common_neighbors(G_old, u, v)))
    #         ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
    #         cn_l3_score = nx.common_neighbors_l3(G_old, u, v)
    #         cn_l3_norm_score = nx.common_neighbors_l3_degree_normalized(G_old, u, v)
    #         clu_score = clu_clo_dict[u][0] + clu_clo_dict[v][0]
    #         clo_score = clu_clo_dict[u][1] + clu_clo_dict[v][1]
    #         X.append([cn_score, ra_score, cn_l3_score, cn_l3_norm_score, clu_score, clo_score])
    #         if (u, v) in G_new.edges():
    #             y.append(1)
    #         else:
    #             y.append(0)

    if number_of_features == 6:
        clu_clo_dict = nx.clustering_closure_coefs(G_old)
        iquad_oquad_dict = nx.iquad_oquad_coefs(G_old)
        for u, v in possible_links:
            cn_score = len(list(nx.common_neighbors(G_old, u, v)))
            ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
            clu_score = clu_clo_dict[u][0] + clu_clo_dict[v][0]
            clo_score = clu_clo_dict[u][1] + clu_clo_dict[v][1]
            iquad_score = iquad_oquad_dict[u][0] + iquad_oquad_dict[v][0]
            oquad_score = iquad_oquad_dict[u][1] + iquad_oquad_dict[v][1]
            X.append([cn_score, ra_score, clu_score, clo_score, iquad_score, oquad_score])
            if (u, v) in G_new.edges():
                y.append(1)
            else:
                y.append(0)

    # if number_of_features == 8:
    #     clu_clo_dict = nx.clustering_closure_coefs(G_old)
    #     iquad_oquad_dict = nx.iquad_oquad_coefs(G_old)
    #     for u, v in possible_links:
    #         cn_score = len(list(nx.common_neighbors(G_old, u, v)))
    #         ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
    #         cn_l3_score = nx.common_neighbors_l3(G_old, u, v)
    #         cn_l3_norm_score = nx.common_neighbors_l3_degree_normalized(G_old, u, v)
    #         clu_score = clu_clo_dict[u][0] + clu_clo_dict[v][0]
    #         clo_score = clu_clo_dict[u][1] + clu_clo_dict[v][1]
    #         iquad_score = iquad_oquad_dict[u][0] + iquad_oquad_dict[v][0]
    #         oquad_score = iquad_oquad_dict[u][1] + iquad_oquad_dict[v][1]
    #         X.append([cn_score, ra_score, cn_l3_score, cn_l3_norm_score, clu_score, clo_score, iquad_score, oquad_score])
    #         if (u, v) in G_new.edges():
    #             y.append(1)
    #         else:
    #             y.append(0)
    return X, y


def get_dataset(G, repeat = 10, old_pct = 0.7):
    dataset = []
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
            random_edge = random.choice(all_edges)
            random_node = random_edge[0]
            sample_nodes = get_sample_nodes(G, random_node, 5000)
            G_sampled = G.subgraph(sample_nodes)
            all_edges = list(G_sampled.edges(data=True))

        for n in range(0, repeat):
            G_old = nx.Graph()
            G_new = nx.Graph()
            random.shuffle(all_edges)
            old_edges = all_edges[:k]
            new_edges = all_edges[k:]
            G_old.add_edges_from(old_edges)
            G_new.add_edges_from(new_edges)
            G_new = G_new.subgraph(G_old.nodes())
            dataset.append([G_old, G_new])
    print("Number of dataset: {}".format(sample_time * repeat))
    return dataset


# APP2: metrics: ROC-AUC, PR-AUC, Average Precision
def link_predict_similarity_based(dataset, method = 'cn'):
    roc_auc = 0
    pr_auc = 0
    ave_precision = 0
    n = len(dataset)
    for g in tqdm(dataset):
        label_all, score_all = get_predicts_and_labels(g[0], g[1], method)
        precision, recall, _ = precision_recall_curve(label_all, score_all)
        roc_auc += roc_auc_score(label_all, score_all)
        pr_auc += auc(recall, precision)
        ave_precision += average_precision_score(label_all, score_all)

    roc_auc /= n
    pr_auc /= n
    ave_precision /= n
    print("{} :\nROC-AUC: {};\nPR-AUC: {};\nAve_Precision: {}.".format(method, roc_auc, pr_auc, ave_precision))


def get_predicts_and_labels(G_old, G_new, method):
    predicts = get_predict_score(G_old, method)
    label_all = []
    score_all = []
    for p in predicts:
        score_all.append(p[2])
        if (p[0], p[1]) in G_new.edges():
            label_all.append(1)
        else:
            label_all.append(0)
    return label_all, score_all


def get_predict_score(G_old, method):
    if method == 'cn':
        predicts = nx.common_neighbor_index(G_old)
    if method == 'ra':
        predicts = nx.resource_allocation_index(G_old)
    if method == 'cn+clu':
        predicts = nx.common_neighbor_plus_clustering(G_old)
    if method == 'cn-l3':
        predicts = nx.common_neighbor_l3_index(G_old)
    if method == 'cn-l3-norm':
        predicts = nx.common_neighbor_l3_degree_normalized_index(G_old)

    max_score = max(predicts, key=itemgetter(2))[2]
    normalized_predicts = [(i[0], i[1], i[2] / max_score) for i in predicts]
    return normalized_predicts


#obselete APP3: using precision
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