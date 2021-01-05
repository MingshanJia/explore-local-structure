from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import homogeneity_score, completeness_score, v_measure_score, accuracy_score
from sklearn.tree import DecisionTreeClassifier

__all__ = ['classify_networks']

def classify_networks(data, labels, repeat=1000):
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


def classify_networks_supervised_loo(X, y_true, repeat, model='tree'):
    l = len(y_true)
    y_pred = np.zeros((l,), dtype=int)
    acc = 0
    loo = LeaveOneOut()
    for n in range(0, repeat):
        for train_index, test_index in loo.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y_true[train_index], y_true[test_index]
            clf = DecisionTreeClassifier()
            clf.fit(X_train,y_train)
            y_pred_i = clf.predict(X_test)
            y_pred[test_index] = y_pred_i
        acc += accuracy_score(y_true, y_pred)
    return acc / repeat