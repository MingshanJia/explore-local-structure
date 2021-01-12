from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import homogeneity_score, completeness_score, v_measure_score, accuracy_score
from sklearn.tree import DecisionTreeClassifier
from tqdm import tqdm

__all__ = ['classify_networks_unsupervised', 'classify_networks_supervised_loo']

def classify_networks_unsupervised(data, labels, repeat=1000):
    data = scale(data)
    n_cluster = len(np.unique(labels))
    homo = 0
    compl = 0
    v_measure = 0
    for i in range(repeat):
        estimator = KMeans(init='random', n_clusters=n_cluster, n_init=1)
        estimator.fit(data)
        h = homogeneity_score(labels, estimator.labels_)
        c = completeness_score(labels, estimator.labels_)
        v = v_measure_score(labels, estimator.labels_)
        if h > homo:
            homo = h
        if c > compl:
            compl = c
        if v > v_measure:
            v_measure = v
            model = estimator
    print("homo:{}  compl:{}   v-measure:{}\n".format(homo, compl, v_measure))
    return model, homo, compl, v_measure


def classify_networks_supervised_loo(X, y_true, model, repeat=1000):
    l = len(y_true)
    assert(np.shape(X)[0] == l)
    y_pred = np.zeros((l,), dtype=int)
    acc_all = 0
    acc_best = 0
    y_pred_best = np.zeros((l,), dtype=int)
    FI_all = np.zeros(np.shape(X)[1])
    loo = LeaveOneOut()
    for n in tqdm(range(repeat)):
        for train_index, test_index in loo.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y_true[train_index], y_true[test_index]
            model_i = model
            model_i.fit(X_train,y_train)
            y_pred_i = model_i.predict(X_test)
            y_pred[test_index] = y_pred_i
        acc = accuracy_score(y_true, y_pred)
        acc_all += acc
        if acc > acc_best:
            acc_best = acc
            y_pred_best = y_pred

        # calculating average feature importance
        model_n = model
        model_n.fit(X, y_true)
        FI_n = model_n.feature_importances_
        FI_all = FI_all + FI_n

    acc_avg = acc_all / repeat
    FI_avg = FI_all / repeat
    print("Average Accuracy: {}\n".format(acc_avg))
    print("Best Accuracy: {}\n".format(acc_best))
    print("Best Prediction: {}\n".format(y_pred_best))
    print("Average FI: {}\n".format(FI_avg))

    return acc_avg, acc_best, y_pred_best, FI_avg

