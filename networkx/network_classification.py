from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
import numpy as np
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import homogeneity_score, completeness_score, v_measure_score, accuracy_score
from tqdm import tqdm
from sklearn.base import clone
from sklearn.inspection import permutation_importance

__all__ = ['classify_networks_unsupervised', 'classify_networks_supervised_loo', 'impurity_decrease_importances', 'dropcols_importances']

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
    r, c = np.shape(X)
    assert(r == l)
    y_pred = np.zeros((l,), dtype=int)
    acc_all = 0
    acc_best = 0
    y_pred_best = np.zeros((l,), dtype=int)
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

    acc_avg = acc_all / repeat
    print("Average Accuracy: {}\n".format(acc_avg))
    print("Best Accuracy: {}\n".format(acc_best))
    print("Best Prediction: {}\n".format(y_pred_best))
    return acc_avg, acc_best, y_pred_best


def impurity_decrease_importances(X, y_true, model, repeat=1000):
    l = len(y_true)
    r, c = np.shape(X)
    assert(r == l)
    FI_all = np.zeros((repeat, c))
    for n in tqdm(range(repeat)):
        model_n = model
        model_n.fit(X, y_true)
        FI_n = model_n.feature_importances_
        FI_all[n] = FI_n
    FI_avg = np.mean(FI_all, axis=0)
    FI_std = np.std(FI_all, axis=0)
    print("FI Average: {}\n".format(FI_avg))
    print("FI Std: {}\n".format(FI_std))
    return FI_all, FI_avg, FI_std


def permutation_importances(X, y_true, model, scoring='accuracy', model_repeat=100, permutation_repeat=5):
    perm_imp_all = []
    for n in tqdm(range(model_repeat)):
        mod = model.fit(X, y_true)
        perm_imp = permutation_importance(mod, X, y_true, scoring, permutation_repeat)
        perm_imp_all.append(perm_imp['importances'])
    FI_all = np.hstack(tuple(perm_imp_all))
    FI_avg = np.mean(FI_all, axis=1)
    FI_std = np.std(FI_all, axis=1)
    print("FI Average: {}\n".format(FI_avg))
    print("FI Std: {}\n".format(FI_std))
    return FI_all, FI_avg, FI_std


# baseline accuracy score is caculated using all features.
# segs: features are grouped into several segments. For example, first segments contains 4 closure patterns, second segment contains 4 clustering features, etc
# segs = 2, when 8 patterns are used, segs = 3 when 4 random features are added
# then the features in one segment are excluded, fit into a model, and using LOO to caculate the accuracy score
# result means accuracy decrease after removing each segment of features
def dropcols_importances(X, y_true, model, segs=2, repeat=10):
    loo = LeaveOneOut()
    l = len(y_true)
    r, c = np.shape(X)
    y_pred = np.zeros((l,), dtype=int)
    y_pred_dropcols = np.zeros((segs, r), dtype=int)
    acc_baseline = 0
    acc_dropcols = np.zeros((segs,), dtype=float)
    for n in tqdm(range(repeat)):
        for train_index, test_index in loo.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y_true[train_index], y_true[test_index]
            model_i = clone(model)
            model_i.fit(X_train, y_train)
            y_pred[test_index] = model_i.predict(X_test)
            i = 0
            for cols in np.split(np.arange(c), segs):
                X_train_dropcols = np.delete(X_train, cols, 1)
                X_test_dropcols = np.delete(X_test, cols, 1)
                model_dropcols = clone(model)
                model_dropcols.fit(X_train_dropcols, y_train)
                y_pred_dropcols[i][test_index] = model_dropcols.predict(X_test_dropcols)
                i = i + 1
        acc_baseline += accuracy_score(y_true, y_pred)
        for seg in range(segs):
            acc_dropcols[seg] += accuracy_score(y_true, y_pred_dropcols[seg])
    acc_baseline = acc_baseline / repeat
    acc_dropcols = acc_dropcols / repeat
    res = acc_baseline - acc_dropcols
    return res