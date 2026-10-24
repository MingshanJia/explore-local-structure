from copy import deepcopy

import networkx as nx
import pytest
from oracles import quad_oracle, weighted_closure_oracle

import explore_local_structure as els


def test_full_graph_normalization_even_for_output_subset():
    G = nx.cycle_graph(4)
    nx.set_edge_attributes(G, 1, "w")
    assert els.inner_quadrangle_coefficient(G, 0, "w") == 1
    # This separate component sets the paper's graph-wide normalization maximum.
    G.add_edge("outside", "component", w=10)
    assert els.inner_quadrangle_coefficient(G, 0, "w") == pytest.approx(0.1)
    assert els.outer_quadrangle_coefficient(G, 0, "w") == pytest.approx(0.1)
    H = nx.complete_graph(3)
    nx.set_edge_attributes(H, 1, "w")
    H.add_edge("outside", "component", w=10)
    assert els.closure(H, 0, "w") == pytest.approx(0.1)


def test_frozen_and_read_only_views_are_not_mutated():
    G = nx.complete_graph(5)
    nx.set_edge_attributes(G, 2, "w")
    before = deepcopy(G)
    for graph in (nx.freeze(G), G.subgraph(range(4))):
        assert els.closure(graph, weight="w")
        assert els.inner_quadrangle_coefficient(graph, weight="w")
        assert els.outer_quadrangle_coefficient(graph, weight="w")
    assert nx.utils.graphs_equal(G, before)


def test_duplicate_nodes_do_not_reweight_an_average():
    G = nx.cycle_graph(4)
    G.add_edge(1, 4)
    assert els.average_inner_quadrangle_coefficient(G, [0, 0, 1]) == (
        els.average_inner_quadrangle_coefficient(G, [0, 1])
    )


@pytest.mark.parametrize("maximum", [1.0, 1e308])
def test_weighted_products_and_normalization_do_not_underflow_before_division(maximum):
    G = nx.Graph()
    G.add_weighted_edges_from(
        [(0, 1, 1e-200), (1, 2, 1e-200), (2, 3, maximum), (3, 0, maximum)]
    )
    numerator, inner, outer = quad_oracle(G, 1, "weight", exact=True)
    pair = els.quadrangle_coefficients(G, 1, "weight")
    assert pair[0] == 1.0
    assert pair == pytest.approx(
        (float(numerator / inner), float(numerator / outer)), abs=0, rel=1e-12
    )
    assert all(isinstance(value, float) for value in pair)


@pytest.mark.parametrize("directed", [False, True])
def test_tiny_weighted_closure_against_exact_oracle(directed):
    G = nx.DiGraph() if directed else nx.Graph()
    G.add_weighted_edges_from([(0, 1, 1e-200), (1, 2, 1e-200), (2, 0, 1)])
    G.add_edge(1, 4, weight=1e-200)
    expected = weighted_closure_oracle(G, 1, "weight", exact=True)
    assert expected > 0
    actual = els.closure(G, 1, "weight")
    assert actual == pytest.approx(float(expected), abs=0, rel=1e-12)


def test_weighted_global_exact_fallback_with_leading_isolate():
    G = nx.Graph()
    G.add_node("isolate")
    G.add_weighted_edges_from([(0, 1, 1e-200), (1, 2, 1e-200), (2, 3, 1), (3, 0, 1)])
    counts = [quad_oracle(G, v, "weight", exact=True) for v in G]
    expected = sum(row[0] for row in counts) / sum(row[1] for row in counts)
    assert expected > 0
    assert els.global_quadrangle_coefficient(G, "weight") == pytest.approx(
        float(expected), abs=0, rel=1e-12
    )


def test_weighted_exact_local_counts_with_open_path_zero_terms():
    G = nx.Graph()
    G.add_weighted_edges_from(
        [
            (0, 1, 1e-200),
            (1, 2, 1e-200),
            (2, 3, 1),
            (3, 0, 1),
            (1, 4, 1e-200),
            (4, 5, 1),
        ]
    )
    for v in G:
        q, inner, outer = quad_oracle(G, v, "weight", exact=True)
        expected = (float(q / inner) if inner else 0, float(q / outer) if outer else 0)
        assert els.quadrangle_coefficients(G, v, "weight") == pytest.approx(
            expected, abs=0, rel=1e-12
        )
