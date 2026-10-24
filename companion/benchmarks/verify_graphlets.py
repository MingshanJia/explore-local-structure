"""Compare modern graphlets with the pinned original implementation.

Run from companion/: python benchmarks/verify_graphlets.py --output results.json
Uses git to read the original source and public Chesapeake food-web dataset.
Types on that dataset are deterministic synthetic test labels, not observations.
The unavailable patient dataset/classification experiment is not reproduced.
"""

import argparse
import hashlib
import importlib.util
import json
import subprocess
from io import StringIO
from itertools import combinations, product
from pathlib import Path
from tempfile import TemporaryDirectory

import networkx as nx

import explore_local_structure as els

PIN = "9454ce0da9fd6e9399e17dc552797add44e25933"
METHODS = (
    "induced_graphlet_degree_vector",
    "induced_graphlet_degree_vector_ego",
    "graphlet_degree_vector_ego",
    "typed_edge_induced_graphlet_degree_vector",
    "typed_edge_induced_graphlet_degree_vector_ego",
    "typed_edge_graphlet_degree_vector_ego",
    "colored_graphlet_vector_for_typed_edge",
    "hetero_graphlet_vector_for_typed_edge",
    "colored_ego_graphlet_vector_for_typed_edge",
    "hetero_ego_graphlet_vector_for_typed_edge",
)


def read_original(repo, path):
    return subprocess.check_output(["git", "show", f"{PIN}:{path}"], cwd=repo)


def evaluate(module, method, G, focal=None, types=3):
    args = () if method in METHODS[:3] else (types,)
    return getattr(module, method)(G, *args, nodes=focal)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--write-fixtures", action="store_true")
    args = parser.parse_args()
    companion = Path(__file__).resolve().parents[1]
    repo = companion.parent
    source = read_original(repo, "networkx/algorithms/graphlets.py")
    # Loading this file defines functions only; no vendored NetworkX is imported.
    with TemporaryDirectory() as temporary:
        path = Path(temporary) / "original_graphlets.py"
        path.write_bytes(source)
        spec = importlib.util.spec_from_file_location("original_graphlets", path)
        legacy = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(legacy)
    comparisons = 0
    cases = []

    def compare(G, name, methods=METHODS, focal=None, types=3, record=False):
        nonlocal comparisons
        expected = {}
        for method in methods:
            old = evaluate(legacy, method, G, focal, types)
            new = evaluate(els, method, G, focal, types)
            if new != old:
                raise AssertionError(f"Mismatch: {name}, {method}")
            comparisons += 1
            if record:
                expected[method] = old[focal]
        if record:
            cases.append(
                {
                    "name": name,
                    "nodes": list(G),
                    "focal": focal,
                    "edges": [
                        (u, v, data["edge_type"]) for u, v, data in G.edges(data=True)
                    ],
                    "expected": expected,
                }
            )

    pairs = list(combinations(range(4), 2))
    for index, state in enumerate(product(range(3), repeat=6)):
        G = nx.Graph()
        G.add_nodes_from(range(4))
        G.add_edges_from(
            (u, v, {"edge_type": str(t)}) for (u, v), t in zip(pairs, state) if t
        )
        compare(G, f"exhaustive-{index}", types=2)
    # Overlapping graphlets on larger graphs, isolated nodes, and diverse topology.
    for seed in range(32):
        G = nx.gnp_random_graph(5 + seed % 4, (seed % 7 + 1) / 8, seed=seed)
        for index, (u, v) in enumerate(G.edges()):
            G[u][v]["edge_type"] = str(index % 3 + 1)
        compare(G, f"random-{seed}", focal=0, record=seed < 12)
    # All alter configurations: absent, unknown type 0, or positive type 1/2.
    ego_methods = tuple(m for m in METHODS if m.endswith("_ego") or "ego_graphlet" in m)
    for index, state in enumerate(product(range(4), repeat=3)):
        G = nx.star_graph(3)
        nx.set_edge_attributes(G, "1", "edge_type")
        for (u, v), value in zip(combinations(range(1, 4), 2), state):
            if value:
                G.add_edge(u, v, edge_type=str(value - 1))
        compare(G, f"zero-alters-{index}", ego_methods, 0, record=index in (1, 21, 63))

    schemas = []
    for types in range(1, 5):
        for method, ego in product(("colored", "hetero"), (False, True)):
            names = els.graphlet_feature_names(method, types, ego=ego)
            decoder = (
                legacy.get_type_from_index_for_ego_network
                if ego
                else legacy.get_type_from_index
            )
            assert names == decoder(
                method + ("_ego" if ego else ""), types, range(len(names))
            )
            schemas.append(
                {"method": method, "ego": ego, "types": types, "names": names}
            )

    data_path = "experiments_for_quad/datasets/FW-Chesa"
    raw_data = read_original(repo, data_path)
    G = nx.read_edgelist(StringIO(raw_data.decode()), data=False, create_using=nx.Graph)
    G.remove_edges_from(list(nx.selfloop_edges(G)))
    for index, (u, v) in enumerate(sorted(tuple(sorted(edge)) for edge in G.edges())):
        G[u][v]["edge_type"] = str(index % 3 + 1)
    compare(G, "FW-Chesa: undirected, synthetic three-type labels")
    report = {
        "legacy_commit": PIN,
        "legacy_source_sha256": hashlib.sha256(source).hexdigest(),
        "graph_method_comparisons": comparisons,
        "coordinate_schemas_compared": len(schemas),
        "matches": "exact integer equality for every returned vector",
        "max_absolute_difference": 0,
        "dataset": {
            "path": data_path,
            "sha256": hashlib.sha256(raw_data).hexdigest(),
            "nodes": len(G),
            "edges": G.number_of_edges(),
            "transformation": "undirected simple graph, loops removed",
            "types": "synthetic 1..3 cyclically over lexicographically sorted undirected edges",
            "note": "This is count validation, not a reproduction of the patient classification study.",
            "example_node": "1",
            "example_typed_ego": els.typed_edge_induced_graphlet_degree_vector_ego(
                G, 3, "1"
            )["1"],
        },
    }
    if args.write_fixtures:
        destination = companion / "tests/data/graphlets_legacy.json"
        destination.parent.mkdir(exist_ok=True)
        destination.write_text(
            json.dumps(
                {
                    "legacy_commit": PIN,
                    "legacy_source_sha256": report["legacy_source_sha256"],
                    "cases": cases,
                    "schemas": schemas,
                },
                separators=(",", ":"),
            )
            + "\n",
            encoding="utf-8",
        )
    if args.output:
        args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
