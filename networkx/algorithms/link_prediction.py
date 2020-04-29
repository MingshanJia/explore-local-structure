"""
Link prediction algorithms.
"""


from math import log

import networkx as nx
from networkx.utils import not_implemented_for


__all__ = ['random_guess',
           'common_neighbor_index',
           'closure_similarity_index',
           'closure_similarity_index_two',
           'closure_similarity_index_three',
           'resource_allocation_index',
           'jaccard_coefficient',
           'adamic_adar_index',
           'preferential_attachment',
           'cn_soundarajan_hopcroft',
           'ra_index_soundarajan_hopcroft',
           'within_inter_cluster',
           'perform_link_prediction']


# KeyFunc: return prediction precision
def perform_link_prediction(G_old, G_new, method, dict_ce, dict_src, dict_tgt):
    G_new = G_new.subgraph(G_old.nodes())
    k = G_new.number_of_edges()    # number of links chosen from prediction, also number of links in ground truth

    if method == 'cn':
        pred_links = common_neighbor_index(G_old)[0 : k]
    if method == 'ja':
        pred_links = jaccard_coefficient(G_old)[0 : k]
    if method == 'aa':
        pred_links = adamic_adar_index(G_old)[0 : k]
    if method == 'ra':
        pred_links = resource_allocation_index(G_old)[0 : k]
    if method == 'clo1':
        pred_links = closure_similarity_index(G_old, dict_ce, dict_src, dict_tgt)[0 : k]
    if method == 'clo2':
        pred_links = closure_similarity_index_two(G_old, dict_src, dict_tgt)[0 : k]
    if method == 'dgr':
        pred_links = degree_similarity_index(G_old)[0 : k]
    # clo3 is only for directed network
    # if method == 'clo3':
    #     pred_links = closure_similarity_index_three(G_old, dict_ce)[0 : k]

    correct = 0
    for e in pred_links:
        if (e[0], e[1]) in G_new.edges():
            correct += 1
    return 100 * correct / k


# random guess precision
def random_guess(G_old, G_new):
    G_new = G_new.subgraph(G_old.nodes())
    k = G_new.number_of_edges()  # number of links chosen from prediction, also number of links in ground truth
    possible_num_edges = len(list(nx.non_edges(G_old)))
    return k * 100 / possible_num_edges


# ChangeNote: return sorted 3-tuple list, according to function score
def _apply_prediction(G, func, ebunch=None):
    """Applies the given function to each edge in the specified iterable
    of edges.

    `G` is an instance of :class:`networkx.Graph`.

    `ebunch` is an iterable of pairs of nodes. If not specified, all
    non-edges in the graph `G` will be used.

    """
    if ebunch is None:
        ebunch = nx.non_edges(G)
    return sorted([(u, v, func(G, u, v)) for u, v in ebunch], key = lambda t:t[2], reverse = True)


# ChangeNote: newly added
def common_neighbor_index(G, ebunch=None):

    def predict(G, u, v):
        if G.is_directed():
            return len(list(nx.directed_common_neighbors(G, u, v)))
        else:
            return len(list(nx.common_neighbors(G, u, v)))

    return _apply_prediction(G, predict, ebunch)


# KeyFunc: newly introduced
def closure_similarity_index(G, dict_ce, dict_src, dict_tgt, ebunch=None):

    def predict(G, u, v):
        if G.is_directed():
            return len(list(nx.directed_common_neighbors(G, u, v))) * (dict_src[u] + dict_tgt[v])
        else:
            return len(list(nx.common_neighbors(G, u, v))) * (dict_ce[u] + dict_ce[v])

    return _apply_prediction(G, predict, ebunch)


# KeyFunc: only for directed network
def closure_similarity_index_two(G, dict_src, dict_tgt, ebunch=None):
# dict_Ce: {v: [clo, src_clo, tgt_clo]}
    def predict(G, u, v):
        return len(list(nx.directed_common_neighbors_two(G, u, v))) * (dict_src[u] + dict_tgt[v])

    return _apply_prediction(G, predict, ebunch)


# using out_degree(s), in_degree(t) and common nbrs info
def degree_similarity_index(G, ebunch=None):
# dict_Ce: {v: [clo, src_clo, tgt_clo]}
    def predict(G, u, v):
        return len(list(nx.directed_common_neighbors_two(G, u, v))) * (G.out_degree[u] + G.in_degree[v])

    return _apply_prediction(G, predict, ebunch)


# Not used.
def closure_similarity_index_three(G, dict_ce, ebunch=None):
# dict_Ce: {v: [clo, src_clo, tgt_clo]}
    def predict(G, u, v):
        if G.is_directed():
            return sum(dict_ce[w][0] for w in nx.directed_common_neighbors(G, u, v))
        else:
            return sum(dict_ce[w][0] for w in nx.common_neighbors(G, u, v))

    return _apply_prediction(G, predict, ebunch)




# ChangeNote: include directed
#@not_implemented_for('directed')
@not_implemented_for('multigraph')
def resource_allocation_index(G, ebunch=None):
    r"""Compute the resource allocation index of all node pairs in ebunch.

    Resource allocation index of `u` and `v` is defined as

    .. math::

        \sum_{w \in \Gamma(u) \cap \Gamma(v)} \frac{1}{|\Gamma(w)|}

    where $\Gamma(u)$ denotes the set of neighbors of $u$.

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        Resource allocation index will be computed for each pair of
        nodes given in the iterable. The pairs must be given as
        2-tuples (u, v) where u and v are nodes in the graph. If ebunch
        is None then all non-existent edges in the graph will be used.
        Default value: None.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their resource allocation index.

    References
    ----------
    .. [1] T. Zhou, L. Lu, Y.-C. Zhang.
       Predicting missing links via local information.
       Eur. Phys. J. B 71 (2009) 623.
       https://arxiv.org/pdf/0901.0553.pdf
    """
    def predict(G, u, v):
        if G.is_directed():
            return sum(1 / G.degree(w) for w in nx.directed_common_neighbors(G, u, v))
        else:
            return sum(1 / G.degree(w) for w in nx.common_neighbors(G, u, v))

    return _apply_prediction(G, predict, ebunch)


# ChangeNote: include directed
#@not_implemented_for('directed')
@not_implemented_for('multigraph')
def jaccard_coefficient(G, ebunch=None):
    r"""Compute the Jaccard coefficient of all node pairs in ebunch.

    Jaccard coefficient of nodes `u` and `v` is defined as

    .. math::

        \frac{|\Gamma(u) \cap \Gamma(v)|}{|\Gamma(u) \cup \Gamma(v)|}

    where $\Gamma(u)$ denotes the set of neighbors of $u$.

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        Jaccard coefficient will be computed for each pair of nodes
        given in the iterable. The pairs must be given as 2-tuples
        (u, v) where u and v are nodes in the graph. If ebunch is None
        then all non-existent edges in the graph will be used.
        Default value: None.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their Jaccard coefficient.

    References
    ----------
    .. [1] D. Liben-Nowell, J. Kleinberg.
           The Link Prediction Problem for Social Networks (2004).
           http://www.cs.cornell.edu/home/kleinber/link-pred.pdf
    """
    def predict(G, u, v):
        if G.is_directed():
            union_size = len(set(G._succ[u]) | set(G._pred[v]))
            if union_size == 0:
                return 0
            return len(list(nx.directed_common_neighbors(G, u, v))) / union_size
        else:
            union_size = len(set(G[u]) | set(G[v]))
            if union_size == 0:
                return 0
            return len(list(nx.common_neighbors(G, u, v))) / union_size

    return _apply_prediction(G, predict, ebunch)


# ChangeNote: include directed
#@not_implemented_for('directed')
@not_implemented_for('multigraph')
def adamic_adar_index(G, ebunch=None):
    r"""Compute the Adamic-Adar index of all node pairs in ebunch.

    Adamic-Adar index of `u` and `v` is defined as

    .. math::

        \sum_{w \in \Gamma(u) \cap \Gamma(v)} \frac{1}{\log |\Gamma(w)|}

    where $\Gamma(u)$ denotes the set of neighbors of $u$.
    This index leads to zero-division for nodes only connected via self-loops.
    It is intended to be used when no self-loops are present.

    Parameters
    ----------
    G : graph
        NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        Adamic-Adar index will be computed for each pair of nodes given
        in the iterable. The pairs must be given as 2-tuples (u, v)
        where u and v are nodes in the graph. If ebunch is None then all
        non-existent edges in the graph will be used.
        Default value: None.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their Adamic-Adar index.

    References
    ----------
    .. [1] D. Liben-Nowell, J. Kleinberg.
           The Link Prediction Problem for Social Networks (2004).
           http://www.cs.cornell.edu/home/kleinber/link-pred.pdf
    """
    def predict(G, u, v):
        if G.is_directed():
            return sum(1 / log(G.degree(w)) for w in nx.directed_common_neighbors(G, u, v))
        else:
            return sum(1 / log(G.degree(w)) for w in nx.common_neighbors(G, u, v))

    return _apply_prediction(G, predict, ebunch)


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def preferential_attachment(G, ebunch=None):
    r"""Compute the preferential attachment score of all node pairs in ebunch.

    Preferential attachment score of `u` and `v` is defined as

    .. math::

        |\Gamma(u)| |\Gamma(v)|

    where $\Gamma(u)$ denotes the set of neighbors of $u$.

    Parameters
    ----------
    G : graph
        NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        Preferential attachment score will be computed for each pair of
        nodes given in the iterable. The pairs must be given as
        2-tuples (u, v) where u and v are nodes in the graph. If ebunch
        is None then all non-existent edges in the graph will be used.
        Default value: None.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their preferential attachment score.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.complete_graph(5)
    >>> preds = nx.preferential_attachment(G, [(0, 1), (2, 3)])
    >>> for u, v, p in preds:
    ...     print(f'({u}, {v}) -> {p}')
    (0, 1) -> 16
    (2, 3) -> 16

    References
    ----------
    .. [1] D. Liben-Nowell, J. Kleinberg.
           The Link Prediction Problem for Social Networks (2004).
           http://www.cs.cornell.edu/home/kleinber/link-pred.pdf
    """
    def predict(u, v):
        return G.degree(u) * G.degree(v)
    return _apply_prediction(G, predict, ebunch)


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def cn_soundarajan_hopcroft(G, ebunch=None, community='community'):
    r"""Count the number of common neighbors of all node pairs in ebunch
        using community information.

    For two nodes $u$ and $v$, this function computes the number of
    common neighbors and bonus one for each common neighbor belonging to
    the same community as $u$ and $v$. Mathematically,

    .. math::

        |\Gamma(u) \cap \Gamma(v)| + \sum_{w \in \Gamma(u) \cap \Gamma(v)} f(w)

    where $f(w)$ equals 1 if $w$ belongs to the same community as $u$
    and $v$ or 0 otherwise and $\Gamma(u)$ denotes the set of
    neighbors of $u$.

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        The score will be computed for each pair of nodes given in the
        iterable. The pairs must be given as 2-tuples (u, v) where u
        and v are nodes in the graph. If ebunch is None then all
        non-existent edges in the graph will be used.
        Default value: None.

    community : string, optional (default = 'community')
        Nodes attribute name containing the community information.
        G[u][community] identifies which community u belongs to. Each
        node belongs to at most one community. Default value: 'community'.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their score.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.path_graph(3)
    >>> G.nodes[0]['community'] = 0
    >>> G.nodes[1]['community'] = 0
    >>> G.nodes[2]['community'] = 0
    >>> preds = nx.cn_soundarajan_hopcroft(G, [(0, 2)])
    >>> for u, v, p in preds:
    ...     print(f'({u}, {v}) -> {p}')
    (0, 2) -> 2

    References
    ----------
    .. [1] Sucheta Soundarajan and John Hopcroft.
       Using community information to improve the precision of link
       prediction methods.
       In Proceedings of the 21st international conference companion on
       World Wide Web (WWW '12 Companion). ACM, New York, NY, USA, 607-608.
       http://doi.acm.org/10.1145/2187980.2188150
    """
    def predict(u, v):
        Cu = _community(G, u, community)
        Cv = _community(G, v, community)
        cnbors = list(nx.common_neighbors(G, u, v))
        neighbors = (sum(_community(G, w, community) == Cu for w in cnbors)
                     if Cu == Cv else 0)
        return len(cnbors) + neighbors
    return _apply_prediction(G, predict, ebunch)


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def ra_index_soundarajan_hopcroft(G, ebunch=None, community='community'):
    r"""Compute the resource allocation index of all node pairs in
    ebunch using community information.

    For two nodes $u$ and $v$, this function computes the resource
    allocation index considering only common neighbors belonging to the
    same community as $u$ and $v$. Mathematically,

    .. math::

        \sum_{w \in \Gamma(u) \cap \Gamma(v)} \frac{f(w)}{|\Gamma(w)|}

    where $f(w)$ equals 1 if $w$ belongs to the same community as $u$
    and $v$ or 0 otherwise and $\Gamma(u)$ denotes the set of
    neighbors of $u$.

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        The score will be computed for each pair of nodes given in the
        iterable. The pairs must be given as 2-tuples (u, v) where u
        and v are nodes in the graph. If ebunch is None then all
        non-existent edges in the graph will be used.
        Default value: None.

    community : string, optional (default = 'community')
        Nodes attribute name containing the community information.
        G[u][community] identifies which community u belongs to. Each
        node belongs to at most one community. Default value: 'community'.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their score.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edges_from([(0, 1), (0, 2), (1, 3), (2, 3)])
    >>> G.nodes[0]['community'] = 0
    >>> G.nodes[1]['community'] = 0
    >>> G.nodes[2]['community'] = 1
    >>> G.nodes[3]['community'] = 0
    >>> preds = nx.ra_index_soundarajan_hopcroft(G, [(0, 3)])
    >>> for u, v, p in preds:
    ...     print(f'({u}, {v}) -> {p:.8f}')
    (0, 3) -> 0.50000000

    References
    ----------
    .. [1] Sucheta Soundarajan and John Hopcroft.
       Using community information to improve the precision of link
       prediction methods.
       In Proceedings of the 21st international conference companion on
       World Wide Web (WWW '12 Companion). ACM, New York, NY, USA, 607-608.
       http://doi.acm.org/10.1145/2187980.2188150
    """
    def predict(u, v):
        Cu = _community(G, u, community)
        Cv = _community(G, v, community)
        if Cu != Cv:
            return 0
        cnbors = nx.common_neighbors(G, u, v)
        return sum(1 / G.degree(w) for w in cnbors
                   if _community(G, w, community) == Cu)
    return _apply_prediction(G, predict, ebunch)


@not_implemented_for('directed')
@not_implemented_for('multigraph')
def within_inter_cluster(G, ebunch=None, delta=0.001, community='community'):
    """Compute the ratio of within- and inter-cluster common neighbors
    of all node pairs in ebunch.

    For two nodes `u` and `v`, if a common neighbor `w` belongs to the
    same community as them, `w` is considered as within-cluster common
    neighbor of `u` and `v`. Otherwise, it is considered as
    inter-cluster common neighbor of `u` and `v`. The ratio between the
    size of the set of within- and inter-cluster common neighbors is
    defined as the WIC measure. [1]_

    Parameters
    ----------
    G : graph
        A NetworkX undirected graph.

    ebunch : iterable of node pairs, optional (default = None)
        The WIC measure will be computed for each pair of nodes given in
        the iterable. The pairs must be given as 2-tuples (u, v) where
        u and v are nodes in the graph. If ebunch is None then all
        non-existent edges in the graph will be used.
        Default value: None.

    delta : float, optional (default = 0.001)
        Value to prevent division by zero in case there is no
        inter-cluster common neighbor between two nodes. See [1]_ for
        details. Default value: 0.001.

    community : string, optional (default = 'community')
        Nodes attribute name containing the community information.
        G[u][community] identifies which community u belongs to. Each
        node belongs to at most one community. Default value: 'community'.

    Returns
    -------
    piter : iterator
        An iterator of 3-tuples in the form (u, v, p) where (u, v) is a
        pair of nodes and p is their WIC measure.

    Examples
    --------
    >>> import networkx as nx
    >>> G = nx.Graph()
    >>> G.add_edges_from([(0, 1), (0, 2), (0, 3), (1, 4), (2, 4), (3, 4)])
    >>> G.nodes[0]['community'] = 0
    >>> G.nodes[1]['community'] = 1
    >>> G.nodes[2]['community'] = 0
    >>> G.nodes[3]['community'] = 0
    >>> G.nodes[4]['community'] = 0
    >>> preds = nx.within_inter_cluster(G, [(0, 4)])
    >>> for u, v, p in preds:
    ...     print(f'({u}, {v}) -> {p:.8f}')
    (0, 4) -> 1.99800200
    >>> preds = nx.within_inter_cluster(G, [(0, 4)], delta=0.5)
    >>> for u, v, p in preds:
    ...     print(f'({u}, {v}) -> {p:.8f}')
    (0, 4) -> 1.33333333

    References
    ----------
    .. [1] Jorge Carlos Valverde-Rebaza and Alneu de Andrade Lopes.
       Link prediction in complex networks based on cluster information.
       In Proceedings of the 21st Brazilian conference on Advances in
       Artificial Intelligence (SBIA'12)
       https://doi.org/10.1007/978-3-642-34459-6_10
    """
    if delta <= 0:
        raise nx.NetworkXAlgorithmError('Delta must be greater than zero')

    def predict(u, v):
        Cu = _community(G, u, community)
        Cv = _community(G, v, community)
        if Cu != Cv:
            return 0
        cnbors = set(nx.common_neighbors(G, u, v))
        within = {w for w in cnbors
                     if _community(G, w, community) == Cu}
        inter = cnbors - within
        return len(within) / (len(inter) + delta)

    return _apply_prediction(G, predict, ebunch)


def _community(G, u, community):
    """Get the community of the given node."""
    node_u = G.nodes[u]
    try:
        return node_u[community]
    except KeyError:
        raise nx.NetworkXAlgorithmError('No community information')
