from copy import deepcopy

import networkx as nx
import pytest
from oracles import quad_oracle

import explore_local_structure as els


def test_atlas_against_distinct_vertex_paths_and_square_transform():
    for G in nx.graph_atlas_g():
        pairs = els.quadrangle_coefficients(G)
        square = nx.square_clustering(G)
        for v in G:
            q, inner, outer = quad_oracle(G, v)
            assert pairs[v] == pytest.approx(
                (q / inner if inner else 0, q / outer if outer else 0)
            )
            assert pairs[v][0] == pytest.approx(2 * square[v] / (1 + square[v]))
        if nx.is_bipartite(G):
            assert els.global_quadrangle_coefficient(G) == pytest.approx(
                nx.bipartite.robins_alexander_clustering(G)
            )
        if len(G) == 4 and len(G.edges()) == 6:
            assert els.number_of_quadrangles(G) == 3


def test_weighted_rational_cycle_fixture():
    G = nx.Graph()
    G.add_weighted_edges_from([(0, 1, 1), (1, 2, 2), (2, 3, 3), (3, 0, 4)])
    assert els.inner_quadrangle_coefficient(G, weight="weight") == pytest.approx(
        {0: 3 / 5, 1: 6 / 7, 2: 2 / 5, 3: 1 / 3}
    )
    assert els.outer_quadrangle_coefficient(G, weight="weight") == pytest.approx(
        {0: 2 / 5, 1: 1 / 3, 2: 3 / 5, 3: 6 / 7}
    )
    assert els.average_inner_quadrangle_coefficient(
        G, weight="weight"
    ) == pytest.approx(23 / 42)
    assert els.average_outer_quadrangle_coefficient(
        G, weight="weight"
    ) == pytest.approx(23 / 42)
    assert els.global_quadrangle_coefficient(G, weight="weight") == pytest.approx(
        12 / 25
    )


def test_weighted_random_graphs_oracle_scale_invariance_and_no_mutation():
    for seed in range(20):
        G = nx.gnp_random_graph(6, 0.6, seed=seed)
        for index, (u, v) in enumerate(G.edges()):
            G[u][v]["w"] = index % 9
        original = deepcopy(G)
        pairs = els.quadrangle_coefficients(G, weight="w")
        for v in G:
            q, inner, outer = quad_oracle(G, v, "w")
            assert pairs[v] == pytest.approx(
                (q / inner if inner else 0, q / outer if outer else 0)
            )
        assert nx.utils.graphs_equal(G, original)
        nx.set_edge_attributes(
            G, {(u, v): d["w"] * 17 for u, v, d in G.edges(data=True)}, "w"
        )
        for v, pair in pairs.items():
            assert els.quadrangle_coefficients(G, v, weight="w") == pytest.approx(pair)
        nx.set_edge_attributes(G, 1, "w")
        for v in G:
            assert els.quadrangle_coefficients(G, v, weight="w") == pytest.approx(
                els.quadrangle_coefficients(G, v)
            )


def test_subsets_labels_loops_chords_and_degeneracy():
    G = nx.cycle_graph(4)
    G.add_edge(1, 4)
    assert els.quadrangle_coefficients(G, 0) == pytest.approx((2 / 3, 1))
    G = nx.relabel_nodes(G, {0: ("ego", 0), 1: "one"})
    before = els.quadrangle_coefficients(G)
    G.add_edge("one", "one", w=1e100)
    assert els.quadrangle_coefficients(G) == before
    assert els.quadrangle_coefficients(G, [("ego", 0), "absent"]) == {
        ("ego", 0): before[("ego", 0)]
    }
    assert els.quadrangle_coefficients(G, ("ego", 0), weight="w") == pytest.approx(
        before[("ego", 0)]
    )
    for function in (
        els.average_inner_quadrangle_coefficient,
        els.average_outer_quadrangle_coefficient,
    ):
        assert function(G, ("ego", 0)) in before[("ego", 0)]
        assert function(nx.Graph()) == 0
        assert function(nx.path_graph(3), count_zeros=False) == 0
    assert els.global_quadrangle_coefficient(nx.Graph()) == 0
    assert els.global_quadrangle_coefficient(nx.path_graph(3)) == 0
    assert els.inner_quadrangle_coefficient(nx.complete_graph(4), 0) == 1
    assert els.outer_quadrangle_coefficient(nx.complete_graph(4), 0) == 1
    nx.set_edge_attributes(G, 0, "w")
    assert set(els.inner_quadrangle_coefficient(G, weight="w").values()) == {0.0}


@pytest.mark.parametrize("G", [nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()])
def test_unsupported_graphs(G):
    for fn in (
        els.inner_quadrangle_coefficient,
        els.outer_quadrangle_coefficient,
        els.quadrangle_coefficients,
        els.global_quadrangle_coefficient,
        els.number_of_quadrangles,
    ):
        with pytest.raises(nx.NetworkXNotImplemented):
            fn(G)


@pytest.mark.parametrize("value", [-1, float("nan"), float("inf"), "one"])
def test_weight_validation(value):
    G = nx.Graph([(0, 1, {"w": value})])
    with pytest.raises(ValueError, match="finite nonnegative"):
        els.inner_quadrangle_coefficient(G, weight="w")
