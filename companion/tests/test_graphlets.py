import json
from copy import deepcopy
from itertools import combinations, product
from pathlib import Path

import networkx as nx
import pytest
from graphlet_oracle import template_counts

import explore_local_structure as els

EGO_ORBITS = (0, 2, 3, 7, 11, 13, 14)


def test_all_binary_typed_four_vertex_graphs_against_template_embeddings():
    pairs = list(combinations(range(4), 2))
    for state in product(range(3), repeat=6):
        G = nx.Graph()
        G.add_nodes_from(range(4))
        G.add_edges_from(
            (u, v, {"edge_type": t}) for (u, v), t in zip(pairs, state) if t
        )
        ordinary = els.induced_graphlet_degree_vector(G)
        typed = els.typed_edge_induced_graphlet_degree_vector(G, 2)
        ego = els.typed_edge_induced_graphlet_degree_vector_ego(G, 2)
        for focal in G:
            counts, histogram = template_counts(G, focal, 2)
            assert ordinary[focal] == counts
            assert typed[focal] == histogram
            assert ego[focal] == [histogram[o] for o in EGO_ORBITS]


def test_original_research_golden_vectors_and_coordinates():
    fixture = json.loads(
        (Path(__file__).parent / "data/graphlets_legacy.json").read_text()
    )
    for case in fixture["cases"]:
        G = nx.Graph()
        G.add_nodes_from(case["nodes"])
        G.add_edges_from((u, v, {"edge_type": t}) for u, v, t in case["edges"])
        for method, expected in case["expected"].items():
            kwargs = {"num_type": 3} if "typed_edge" in method else {}
            if "colored" in method or "hetero" in method:
                kwargs["num_type"] = 3
            actual = getattr(els, method)(G, nodes=case["focal"], **kwargs)
            assert actual == {case["focal"]: expected}, (case["name"], method)
    for entry in fixture["schemas"]:
        names = els.graphlet_feature_names(
            entry["method"], entry["types"], ego=entry["ego"]
        )
        assert [[orbit, list(types)] for orbit, types in names] == entry["names"]


def test_histogram_is_edge_count_not_graphlet_count_and_noninduced_multiplicity():
    G = nx.complete_graph(4)
    nx.set_edge_attributes(G, 1, "edge_type")
    assert els.induced_graphlet_degree_vector_ego(G, 0)[0] == [3, 0, 3, 0, 0, 0, 1]
    assert els.typed_edge_induced_graphlet_degree_vector_ego(G, 1, 0)[0] == [
        [3],
        [0],
        [9],
        [0],
        [0],
        [0],
        [6],
    ]
    assert els.graphlet_degree_vector_ego(G, 0)[0] == [3, 3, 3, 1, 3, 3, 1]
    assert els.typed_edge_graphlet_degree_vector_ego(G, 1, 0)[0] == [
        [3],
        [6],
        [9],
        [3],
        [12],
        [15],
        [6],
    ]


def test_zero_types_relabeling_views_subsets_and_immutability():
    G = nx.complete_graph(4)
    nx.set_edge_attributes(G, "2", "edge_type")
    G[1][2]["edge_type"] = "0"
    before = deepcopy(G)
    typed = els.typed_edge_induced_graphlet_degree_vector(G, 2)
    assert typed[0][14] == [0, 5]
    assert els.typed_edge_induced_graphlet_degree_vector_ego(G, 2)[0][-1] == [0, 5]
    names = els.graphlet_feature_names("hetero", 2, ego=True)
    hetero = els.hetero_ego_graphlet_vector_for_typed_edge(G, 2, 0)[0]
    assert hetero[names.index((14, (0, 2, 2, 2, 2, 2)))] == 1
    with pytest.raises(ValueError, match="positive focal-edge"):
        els.colored_ego_graphlet_vector_for_typed_edge(G, 2)
    with pytest.raises(ValueError, match="Full colored/hetero"):
        els.hetero_graphlet_vector_for_typed_edge(G, 2)
    H = nx.relabel_nodes(G, {0: ("ego", 0), 1: "alter"})
    H.add_edge("alter", "alter", edge_type="ignored invalid loop")
    assert els.typed_edge_induced_graphlet_degree_vector(
        nx.freeze(H), 2, ("ego", 0)
    ) == {("ego", 0): typed[0]}
    assert els.typed_edge_induced_graphlet_degree_vector(G, 2, [0, 0, "absent"]) == {
        0: typed[0]
    }
    assert els.typed_edge_induced_graphlet_degree_vector(G.subgraph(range(3)), 2)
    assert nx.utils.graphs_equal(G, before)


@pytest.mark.parametrize("value", [-1, 3, 1.5, float("nan"), True, "bad", None])
def test_invalid_type_values(value):
    G = nx.Graph([(0, 1, {"edge_type": value})])
    with pytest.raises(ValueError, match="Every non-loop edge"):
        els.typed_edge_induced_graphlet_degree_vector(G, 2)


@pytest.mark.parametrize("types", [0, -1, True, 1.5, "2"])
def test_invalid_num_types(types):
    with pytest.raises(ValueError, match="positive integer"):
        els.typed_edge_induced_graphlet_degree_vector(nx.Graph(), types)


@pytest.mark.parametrize("G", [nx.DiGraph(), nx.MultiGraph(), nx.MultiDiGraph()])
def test_unsupported_graphs(G):
    with pytest.raises(nx.NetworkXNotImplemented):
        els.typed_edge_induced_graphlet_degree_vector(G, 2)


def test_empty_inputs_attributes_and_feature_indices():
    assert els.typed_edge_induced_graphlet_degree_vector(nx.Graph(), 2) == {}
    G = nx.empty_graph(1)
    assert els.typed_edge_induced_graphlet_degree_vector(G, 2) == {
        0: [[0, 0] for _ in range(15)]
    }
    G.add_edge(0, 1, relation="1")
    assert els.typed_edge_induced_graphlet_degree_vector(G, 1, 0, edge_type="relation")[
        0
    ][0] == [1]
    with pytest.raises(ValueError, match="Every non-loop edge"):
        els.typed_edge_induced_graphlet_degree_vector(G, 1)
    for ego in (False, True):
        for method in ("colored", "hetero"):
            names = els.graphlet_feature_names(method, 2, ego=ego)
            if ego:
                decoded = els.get_type_from_index_for_ego_network(
                    method + "_ego", 2, range(len(names))
                )
                count = getattr(
                    els,
                    f"get_total_number_of_types_of_{method}_ego_graphlets_for_typed_edge",
                )(2)
            else:
                decoded = els.get_type_from_index(method, 2, range(len(names)))
                count = getattr(
                    els,
                    f"get_total_number_of_types_of_{method}_graphlets_for_typed_edge",
                )(2)
            assert decoded == names
            assert count == len(names)
    with pytest.raises(ValueError):
        els.get_type_from_index("colored", 2, [-1])
    with pytest.raises(ValueError):
        els.get_type_from_index_for_ego_network("colored_ego", 2, [10000])
    with pytest.raises(ValueError):
        els.graphlet_feature_names("unknown", 2)
