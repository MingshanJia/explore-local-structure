from sklearn import metrics
from sklearn.cluster import KMeans
from sklearn.preprocessing import scale
import numpy as np

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
        h = metrics.homogeneity_score(labels, estimator.labels_)
        c = metrics.completeness_score(labels, estimator.labels_)
        v = metrics.v_measure_score(labels, estimator.labels_)
        if h > homo:
            homo = h
        if c > compl:
            compl = c
        if v > v_measure:
            v_measure = v
            model = estimator
    print("homo:{}  compl:{}   v-measure:{}\n".format(homo, compl, v_measure))
    return model, homo, compl, v_measure