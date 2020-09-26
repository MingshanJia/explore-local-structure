import random
import networkx as nx
from math import log
from queue import Queue
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

__all__ = ['link_pred_supervised_learning', 'get_dataset', 'BFS_sampling',
           'link_pred_similarity_based', 'link_pred_in_precision', 'link_pred_directed_network',]


# compare with 4 feature sets:
# set 1: 5 baseline features
# set 2: baseline features + iquad
# set 3: baseline features + oquad
# set 4: baseline features + iquad + oquad
def link_pred_supervised_learning(G, method="xgboost", sample_size=3000, sample_time=5, repeat=10, train_pct=0.7, train_old_pct=0.7):
    dataset = get_dataset(G, sample_time, sample_size, repeat, train_pct, train_old_pct, supervised=True, directed=False)
    #positive_ratio = 0
    train_roc_auc_1 = 0
    train_roc_auc_2 = 0
    train_roc_auc_3 = 0
    train_roc_auc_4 = 0
    test_roc_auc_1 = 0
    test_roc_auc_2 = 0
    test_roc_auc_3 = 0
    test_roc_auc_4 = 0
    feature_importance = np.zeros(7)
    n = len(dataset)
    for g in tqdm(dataset):
        #positive_ratio += g[1].number_of_edges() / len(list(nx.non_edges(g[0])))
        y_train, train_score_1, train_score_2, train_score_3, train_score_4, \
        y_test, test_score_1, test_score_2, test_score_3, test_score_4, feature_importance_g \
            = get_predicts_labels_and_feature_importance(method, g[0], g[1], g[2], g[3])
        train_roc_auc_1 += roc_auc_score(y_train, train_score_1)
        train_roc_auc_2 += roc_auc_score(y_train, train_score_2)
        train_roc_auc_3 += roc_auc_score(y_train, train_score_3)
        train_roc_auc_4 += roc_auc_score(y_train, train_score_4)
        test_roc_auc_1 += roc_auc_score(y_test, test_score_1)
        test_roc_auc_2 += roc_auc_score(y_test, test_score_2)
        test_roc_auc_3 += roc_auc_score(y_test, test_score_3)
        test_roc_auc_4 += roc_auc_score(y_test, test_score_4)
        feature_importance += feature_importance_g
    #positive_ratio /= n
    train_roc_auc_1 /= n
    train_roc_auc_2 /= n
    train_roc_auc_3 /= n
    train_roc_auc_4 /= n
    test_roc_auc_1 /= n
    test_roc_auc_2 /= n
    test_roc_auc_3 /= n
    test_roc_auc_4 /= n
    feature_importance /= n
    train_compare_2_to_1 =  (train_roc_auc_2 - train_roc_auc_1) * 100 / train_roc_auc_1
    train_compare_3_to_1 = (train_roc_auc_3 - train_roc_auc_1) * 100 / train_roc_auc_1
    train_compare_4_to_1 = (train_roc_auc_4 - train_roc_auc_1) * 100 / train_roc_auc_1
    test_compare_2_to_1 =  (test_roc_auc_2 - test_roc_auc_1) * 100 / test_roc_auc_1
    test_compare_3_to_1 = (test_roc_auc_3 - test_roc_auc_1) * 100 / test_roc_auc_1
    test_compare_4_to_1 = (test_roc_auc_4 - test_roc_auc_1) * 100 / test_roc_auc_1
    print("Model: {}, result in ROC-AUC".format(method))
    print("Train result:")
    print("  with baseline features:                   {:.4f}".format(train_roc_auc_1))
    print("  with baseline features + i-quad:          {:.4f}, {:.3f}%".format(train_roc_auc_2, train_compare_2_to_1))
    print("  with baseline features + o-quad:          {:.4f}, {:.3f}%".format(train_roc_auc_3, train_compare_3_to_1))
    print("  with baseline features + i-quad + o-quad: {:.4f}, {:.3f}%".format(train_roc_auc_4, train_compare_4_to_1))
    print("Test result:")
    print("  with baseline features:                   {:.4f}".format(test_roc_auc_1))
    print("  with baseline features + i-quad:          {:.4f}, {:.3f}%".format(test_roc_auc_2, test_compare_2_to_1))
    print("  with baseline features + o-quad:          {:.4f}, {:.3f}%".format(test_roc_auc_3, test_compare_3_to_1))
    print("  with baseline features + i-quad + o-quad: {:.4f}, {:.3f}%".format(test_roc_auc_4, test_compare_4_to_1))
    features = ['cn', 'aa', 'ra', 'clu', 'clo', 'i-quad', 'o-quad']
    plt.bar(features, feature_importance)
    plt.show()
    for feature, score in zip(features, feature_importance):
        print(feature, score)
    return [train_roc_auc_1, [train_roc_auc_2, train_compare_2_to_1], [train_roc_auc_3, train_compare_3_to_1], [train_roc_auc_4, train_compare_4_to_1]],\
           [test_roc_auc_1, [test_roc_auc_2, test_compare_2_to_1], [test_roc_auc_3, test_compare_3_to_1], [test_roc_auc_4, test_compare_4_to_1]], feature_importance


def get_predicts_labels_and_feature_importance(method, G_train, G_test, G_train_old, G_train_new):
    X_train, y_train = get_data_targets(G_train_old, G_train_new)  # train is on G_train_old + G_train_new
    X_test, y_test = get_data_targets(G_train, G_test)   # test in on G_train + G_test
    X_train = normalize(X_train)
    X_test = normalize(X_test)
    X_train = np.array(X_train)
    X_test = np.array(X_test)

    if method == "xgboost":
        model = XGBClassifier()
    if method == "log-reg":
        model = LogisticRegression()
    # baseline features
    model.fit(X_train[:, :5], y_train)
    train_score_1 = model.predict_proba(X_train[:, :5])
    test_score_1 = model.predict_proba(X_test[:, :5])
    # baseline features + iquad
    model.fit(X_train[:, :6], y_train)
    train_score_2 = model.predict_proba(X_train[:, :6])
    test_score_2 = model.predict_proba(X_test[:, :6])
    # baseline features + oquad
    model.fit(X_train[:, [0, 1, 2, 3, 4, 6]], y_train)
    train_score_3 = model.predict_proba(X_train[:, [0, 1, 2, 3, 4, 6]])
    test_score_3 = model.predict_proba(X_test[:, [0, 1, 2, 3, 4, 6]])
    # baseline features + iquad + oquad
    model.fit(X_train, y_train)
    train_score_4 = model.predict_proba(X_train)
    test_score_4 = model.predict_proba(X_test)
    if method == "log-reg":
        importances = model.coef_[0]
    else:
        importances = model.feature_importances_
    return y_train, train_score_1[:, 1], train_score_2[:, 1], train_score_3[:, 1], train_score_4[:, 1], \
           y_test,  test_score_1[:, 1],  test_score_2[:, 1],  test_score_3[:, 1],  test_score_4[:, 1],  importances


# get data (with all 7 features) and targets
def get_data_targets(G_old, G_new):
    possible_links = list(nx.non_edges(G_old))
    X = []
    y = []
    clu_clo_dict = nx.clustering_closure_coefs(G_old)
    iquad_oquad_dict = nx.iquad_oquad_coefs(G_old)
    for u, v in possible_links:
        cn_score = len(list(nx.common_neighbors(G_old, u, v)))
        aa_score = sum(1 / log(G_old.degree(w)) for w in nx.common_neighbors(G_old, u, v))
        ra_score = sum(1 / G_old.degree(w) for w in nx.common_neighbors(G_old, u, v))
        clu_score = clu_clo_dict[u][0] + clu_clo_dict[v][0]
        clo_score = clu_clo_dict[u][1] + clu_clo_dict[v][1]
        iquad_score = iquad_oquad_dict[u][0] + iquad_oquad_dict[v][0]
        oquad_score = iquad_oquad_dict[u][1] + iquad_oquad_dict[v][1]
        X.append([cn_score, aa_score, ra_score, clu_score, clo_score, iquad_score, oquad_score])
        if (u, v) in G_new.edges():
            y.append(1)
        else:
            y.append(0)
    return X, y


# dataset with timestamp: set repeat = 1
# plit into train and test dataset, then split train set into train_old and train_new in order to fit into the model,
# and evaluate the model in test set.
def get_dataset(G, sample_time=5, sample_size=3000, repeat=10, train_pct=0.7, train_old_pct=0.7, supervised=True, directed=False):
    dataset = []
    if G.number_of_nodes() > 8000:
        sample = True
        print("sample {} nodes for {} times".format(sample_size, sample_time))
    else:
        sample = False
        sample_time = 1
        print("no sampling")

    print("train test split time: {}: using {:.2f} % of edges to predict {:.2f} % edges".
          format(repeat, 100 * train_pct, (100 * (1 - train_pct))))
    for i in range(0, sample_time):
        if sample:
            G_sampled = BFS_sampling(G, sample_size)
            all_edges = list(G_sampled.edges(data=True))
        else:
            all_edges = list(G.edges(data=True))
        k = round(len(all_edges) * train_pct)

        for n in range(0, repeat):
            if repeat == 1:
                all_edges = sorted(all_edges, key=lambda t: t[2].get('sec'))  # for dataset with timestamp
            else:
                random.shuffle(all_edges)

            if not directed:
                G_train = nx.Graph()
                G_test = nx.Graph()
            else:
                G_train = nx.DiGraph()
                G_test = nx.DiGraph()
            train_edges = all_edges[:k]
            test_edges = all_edges[k:]
            G_train.add_edges_from(train_edges)
            G_test.add_edges_from(test_edges)
            G_test = G_test.subgraph(G_train.nodes())

            if supervised:
                l = round(k * train_old_pct)
                if not directed:
                    G_train_old = nx.Graph()
                    G_train_new = nx.Graph()
                else:
                    G_train_old = nx.DiGraph()
                    G_train_new = nx.DiGraph()
                train_old_edges = train_edges[:l]
                train_new_edges = train_edges[l:]
                G_train_old.add_edges_from(train_old_edges)
                G_train_new.add_edges_from(train_new_edges)
                G_train_new = G_train_new.subgraph(G_train_old.nodes())
                dataset.append([G_train, G_test, G_train_old, G_train_new])
            else:
                dataset.append([G_train, G_test])
    print("Number of dataset: {}".format(sample_time * repeat))
    return dataset


# APP2: metrics: ROC-AUC, PR-AUC, Average Precision
def link_pred_similarity_based(G, method='cn', sample_size=5000, sample_time=5, repeat=10, old_pct=0.5):
    dataset = get_dataset(G, sample_size, sample_time, repeat, old_pct,supervised=False, directed=False)
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
    return roc_auc, pr_auc, ave_precision


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


# APP3: similarity based methods evaluated in precision
def link_pred_in_precision(G, sample_size=5000, sample_time=10, repeat=10, old_pct=0.5):
    dataset = get_dataset(G, sample_size, sample_time, repeat, old_pct, supervised=False, directed=False)
    rg = 0  # random guess
    cn = 0
    ra = 0
    cn_clu = 0
    cn_l3 = 0
    cn_l3_norm = 0
    n = len(dataset)
    for g in tqdm(dataset):
        rg += nx.random_guess(g[0], g[1])
        cn += nx.perform_link_prediction_undir(g[0], g[1], 'cn')
        ra += nx.perform_link_prediction_undir(g[0], g[1], 'ra')
        cn_clu += nx.perform_link_prediction_undir(g[0], g[1], 'cn+clu')
        cn_l3 += nx.perform_link_prediction_undir(g[0], g[1], 'cn-l3')
        cn_l3_norm += nx.perform_link_prediction_undir(g[0], g[1], 'cn-l3-norm')
    rg /= n
    cn /= n
    ra /= n
    cn_clu /= n
    cn_l3 /= n
    cn_l3_norm /= n
    print('rg: %.3f' % rg)
    print('cn: %.3f' % cn)
    print('ra: %.3f' % ra)
    print('cn+clu: %.3f' % cn_clu)
    print('cn-l3: %.3f' % cn_l3)
    print('cn-l3-norm: %.3f' % cn_l3_norm)
    return rg, cn, ra, cn_clu, cn_l3, cn_l3_norm


# APP4: link prediction for directed network
def link_pred_directed_network(G, sample_size=5000, sample_time=10, repeat=10, old_pct=0.5):
    dataset = get_dataset(G, sample_size, sample_time, repeat, old_pct, supervised=False, directed=True)
    rg = 0
    cn = 0
    aa = 0
    ra = 0
    clo1 = 0
    clo2 = 0
    n = len(dataset)
    for g in tqdm(dataset):
        dict_ce = nx.closure(g[0])
        rg += nx.random_guess(g[0], g[1])
        cn += nx.perform_link_prediction(g[0], g[1], 'cn', dict_ce)
        aa += nx.perform_link_prediction(g[0], g[1], 'aa', dict_ce)
        ra += nx.perform_link_prediction(g[0], g[1], 'ra', dict_ce)
        clo1 += nx.perform_link_prediction(g[0], g[1], 'clo1', dict_ce)
        clo2 += nx.perform_link_prediction(g[0], g[1], 'clo2', dict_ce)
    rg /= n
    cn /= n
    aa /= n
    ra /= n
    clo1 /= n
    clo2 /= n
    print('rg: %.3f' % rg)
    print('cn: %.3f' % cn)
    print('aa: %.3f' % aa)
    print('ra: %.3f' % ra)
    print('clo1: %.3f' % clo1)
    print('clo2: %.3f' % clo2)
    return rg, cn, aa, ra, clo1, clo2


# get sample graph
def BFS_sampling(G, sample_size=2000):
    all_nodes = list(G.nodes)
    sampled_nodes = []  # result list of connected nodes
    queue = Queue()
    while len(sampled_nodes) < sample_size:
        random_node = random.choice(all_nodes)
        sampled_nodes.append(random_node)
        queue.put(random_node)
        while (len(sampled_nodes) < sample_size) and (not queue.empty()):
            n = queue.get()
            neigbhors = list(nx.all_neighbors(G, n))
            random.shuffle(neigbhors)
            for nbr in neigbhors:
                if nbr not in sampled_nodes:
                    sampled_nodes.append(nbr)
                    queue.put(nbr)
                    if len(sampled_nodes) == sample_size:
                        return G.subgraph(sampled_nodes)
    return G.subgraph(sampled_nodes[: sample_size])
