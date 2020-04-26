import random
import networkx as nx

__all__ = ['link_pred_app']

# input graph G, repeat = 1 for dataset with timestamp, old_pct is the the percentage of old edges, based on which we predict new edges
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

    rg /= repeat
    cn /= repeat
    aa /= repeat
    ra /= repeat
    clo1 /= repeat
    clo2 /= repeat

    res.append(rg)
    res.append(cn)
    res.append(aa)
    res.append(ra)
    res.append(clo1)
    res.append(clo2)
    return res