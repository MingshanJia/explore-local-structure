# Explore Local Structure

Research code and a modern Python companion package for directed closure
coefficients, their patterns, and inner/outer quadrangle coefficients.

## Modern package

The maintained package is in [`companion/`](companion). It uses standard
NetworkX graphs and installs as `explore-local-structure`, with the import
name `explore_local_structure`. It does not replace NetworkX.

Install directly from this repository:

```shell
python -m pip install "git+https://github.com/MingshanJia/explore-local-structure.git#subdirectory=companion"
```

Or install from a local checkout:

```shell
python -m pip install ./companion
```

Python 3.12+ (excluding 3.14.1) and NetworkX 3.7+ are supported. Run modern
examples from `companion/` or outside the repository root: the root contains
the original research copy of `networkx`, which Python would otherwise import
instead of your installed NetworkX.

```python
import networkx as nx
import explore_local_structure as els

G = nx.cycle_graph(4)
G.add_edge(1, 4)
assert els.quadrangle_coefficients(G, 0) == (2 / 3, 1.0)

D = nx.DiGraph([(0, 1), (0, 2), (1, 2), (1, 3)])
assert els.directed_closure(D, 0, pattern="source") == 2 / 3
```

The package includes i-quad and o-quad, their averages, weighted versions,
combined pairs, global quadrangle ratios and four-cycle counts. Directed
closure includes total/source/target/head/mid/end/cyclic patterns. See the
[API overview](companion/README.md), [definitions](companion/docs/methods.md)
and [migration guide](companion/docs/migration.md).

For simple loopless unweighted undirected graphs, i-quad is equivalent to
`2*s/(1+s)`, where `s = nx.square_clustering(G, v)` in NetworkX 3.7.
I-quad is provided alongside o-quad with this equivalence explicit; the
weighted definition is separate. These companion APIs are distributed by this
repository and are not part of official NetworkX.

## Original research

The original vendored NetworkX, experiments, datasets and graphlet work remain
in their existing locations. See the [original README](README.rst).
The modern package ports closure and quadrangle methods; ego/typed/colored/
heterogeneous graphlet APIs remain in the original research code.

- [Directed closure coefficient and its patterns (PLOS ONE, 2021)](https://doi.org/10.1371/journal.pone.0253822)
- [Measuring Quadrangle Formation in Complex Networks (IEEE TNSE, 2022)](https://doi.org/10.1109/TNSE.2021.3123735)

## Development

```shell
cd companion
python -m pip install -e ".[test,dev]"
python -m pytest
python -m ruff check .
python -m ruff format --check .
python -m build
```
