"""Reproducible wall-time/peak-memory samples; no promised speedup target."""

import json
import tracemalloc
from time import perf_counter

import networkx as nx

import explore_local_structure as els


def main():
    graphs = {
        "sparse_er": nx.gnp_random_graph(200, 0.03, seed=41),
        "preferential_attachment": nx.barabasi_albert_graph(200, 3, seed=41),
        "dense_neighborhoods": nx.gnp_random_graph(60, 0.4, seed=41),
    }
    results = []
    for name, graph in graphs.items():
        nx.set_edge_attributes(
            graph, {(u, v): 1 + k % 7 for k, (u, v) in enumerate(graph.edges())}, "w"
        )
        directed = graph.to_directed()
        for method, fn, G, kwargs in [
            ("total_closure", els.closure, directed, {}),
            ("head_closure", els.directed_closure, directed, {"pattern": "head"}),
            ("quadrangle_pair", els.quadrangle_coefficients, graph, {}),
            (
                "weighted_quadrangle_pair",
                els.quadrangle_coefficients,
                graph,
                {"weight": "w"},
            ),
        ]:
            tracemalloc.start()
            start = perf_counter()
            fn(G, **kwargs)
            elapsed = perf_counter() - start
            _, peak = tracemalloc.get_traced_memory()
            tracemalloc.stop()
            results.append(
                {
                    "graph": name,
                    "nodes": len(G),
                    "edges": G.number_of_edges(),
                    "method": method,
                    "seconds": elapsed,
                    "peak_bytes": peak,
                }
            )
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
