"""Independent rooted-template embedding oracle (no degree classification)."""

from itertools import combinations, permutations

# Root is vertex 0. Each entry is (number of vertices, template edges).
TEMPLATES = (
    (2, ((0, 1),)),
    (3, ((0, 1), (1, 2))),
    (3, ((0, 1), (0, 2))),
    (3, ((0, 1), (0, 2), (1, 2))),
    (4, ((0, 1), (1, 2), (2, 3))),
    (4, ((0, 1), (0, 2), (2, 3))),
    (4, ((0, 1), (1, 2), (1, 3))),
    (4, ((0, 1), (0, 2), (0, 3))),
    (4, ((0, 1), (1, 2), (2, 3), (0, 3))),
    (4, ((0, 1), (1, 2), (1, 3), (2, 3))),
    (4, ((0, 1), (0, 2), (1, 2), (2, 3))),
    (4, ((0, 1), (0, 2), (0, 3), (1, 2))),
    (4, ((0, 1), (0, 2), (1, 2), (1, 3), (2, 3))),
    (4, ((0, 1), (0, 2), (0, 3), (1, 2), (1, 3))),
    (4, tuple(combinations(range(4), 2))),
)


def template_counts(G, focal, num_type):
    counts, typed = [], []
    for size, edges in TEMPLATES:
        edge_set = {frozenset(e) for e in edges}
        automorphisms = sum(
            {frozenset((mapping[u], mapping[v])) for u, v in edges} == edge_set
            for tail in permutations(range(1, size))
            for mapping in [(0, *tail)]
        )
        count = 0
        histogram = [0] * num_type
        for tail in permutations([v for v in G if v != focal], size - 1):
            mapping = (focal, *tail)
            if all(
                G.has_edge(mapping[u], mapping[v]) == (frozenset((u, v)) in edge_set)
                for u, v in combinations(range(size), 2)
            ):
                count += 1
                for u, v in edges:
                    t = int(G[mapping[u]][mapping[v]]["edge_type"])
                    if t:
                        histogram[t - 1] += 1
        assert count % automorphisms == 0
        assert all(value % automorphisms == 0 for value in histogram)
        counts.append(count // automorphisms)
        typed.append([value // automorphisms for value in histogram])
    return counts, typed
