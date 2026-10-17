# Explore Local Structure

A companion package for endpoint-based triangle closure and inner/outer
quadrangle coefficients. It uses standard NetworkX graphs and does not replace
or install a second package named `networkx`.

```shell
python -m pip install .
```

Python 3.12+ (excluding 3.14.1) and NetworkX 3.7+ are supported.

```python
import networkx as nx
import explore_local_structure as els

D = nx.DiGraph([(0, 1), (0, 2), (1, 2), (1, 3)])
assert els.directed_closure(D, 0) == 1 / 3
assert els.directed_closure(D, 0, pattern="source") == 2 / 3

G = nx.cycle_graph(4)
G.add_edge(1, 4)
i_quad, o_quad = els.quadrangle_coefficients(G, 0)
assert (i_quad, o_quad) == (2 / 3, 1.0)
```

## Methods and scope

| Method | Input | Weights |
|---|---|---|
| `closure`, `average_closure` | Graph / DiGraph | Unweighted or finite nonnegative weights |
| `directed_closure`, `average_directed_closure` | DiGraph | Total as above; six non-total patterns are unweighted |
| `inner_quadrangle_coefficient`, `outer_quadrangle_coefficient` | Graph | Unweighted or finite nonnegative weights |
| `quadrangle_coefficients` | Graph | Computes the pair in one counting pass |
| `average_inner_quadrangle_coefficient`, `average_outer_quadrangle_coefficient` | Graph | Arithmetic local means |
| `global_quadrangle_coefficient` | Graph | Ratio of summed counts; weighted form is a derived legacy extension |
| `number_of_quadrangles` | Graph | Unweighted four-cycle count, including chords |

The directed patterns are `total`, `source`, `target`, `head`, `mid`, `end`,
and `cyclic`. Weighted non-total patterns raise `NetworkXNotImplemented` because
the cited paper does not specify their weighted formulas. Signed inputs raise
`ValueError`; a signed normalization convention needs separate scientific review.

For simple loopless **unweighted** undirected graphs, i-quad is exactly
`2*s/(1+s)` where `s = nx.square_clustering(G, v)` in NetworkX 3.7. We expose
i-quad alongside o-quad to keep both inner- and outer-path perspectives explicit,
and to provide the separately defined weighted versions. O-quad uses a different
denominator and is not this transformation. See [definitions](docs/methods.md).

Both general and bipartite graphs are accepted. Quadrangles with chords count.
Self-loops are ignored in both counting and weight normalization. Parallel edges
are unsupported. Missing weight attributes default to one; the normalization
maximum uses the full graph, even when selecting output nodes. Inputs are never
mutated. Undefined local/global values and empty means return zero by explicit
software convention.

## Research provenance

- Hao Yin, Austin R. Benson and Jure Leskovec. The Local Closure Coefficient:
  A New Perspective On Network Clustering. WSDM 2019.
  [DOI](https://doi.org/10.1145/3289600.3290991).
- Mingshan Jia, Bogdan Gabrys and Katarzyna Musial. Directed closure coefficient
  and its patterns. PLOS ONE 16(6), 2021.
  [DOI](https://doi.org/10.1371/journal.pone.0253822).
- Mingshan Jia, Bogdan Gabrys and Katarzyna Musial. Measuring Quadrangle
  Formation in Complex Networks. IEEE TNSE 9(2), 538–551, 2022.
  [DOI](https://doi.org/10.1109/TNSE.2021.3123735).

The original research snapshot is preserved at commit
`9454ce0da9fd6e9399e17dc552797add44e25933` and the local
`legacy-research-snapshot` branch. Its datasets, vendored NetworkX and experiment
scripts remain accessible there. This modern branch intentionally replaces that
vendored tree. Ego/typed/colored/heterogeneous graphlet ports remain a separate
track; this release does not claim to port those experimental APIs.

See [migration](docs/migration.md) and [contribution preparation](docs/upstream.md).

## Development

```shell
python -m pip install -e ".[test,dev]"
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m build
```

The tests compare against independent orientation/path enumeration, exhaustive
three-/four-node digraphs, graph-atlas graphs, weighted rational fixtures,
isomorphism/reversal identities, input immutability and edge cases. The repository
includes no patient datasets in the installable package.
