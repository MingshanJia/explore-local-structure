# Edge-typed graphlets

Two- to four-node representations from
[Jia et al., Encoding edge type information in graphlets, PLOS ONE (2022)](https://doi.org/10.1371/journal.pone.0273609),
using standard NetworkX graphs.

## Three representations

| Representation | One graphlet contributes | Output per node |
|---|---|---|
| TyE-GDV | One increment for **each edge** of each positive type | 15 × T integer matrix; ego: 7 × T |
| ColoredE-GDV | One increment for its orbit and **set** of types | Flat integer list |
| HeteroE-GDV | One increment for its orbit and **multiset** of types | Flat integer list |

A triangle with types `(1, 1, 2)` contributes `[2, 1]` to its TyE triangle row,
one to colored coordinate `(3, (1, 2))`, and one to heterogeneous coordinate
`(3, (1, 1, 2))`. These representations do not distinguish arrangements of
types on individual edge positions within a given graphlet orbit.

## Example

Run from `companion/` or outside the repository root.

```python
import networkx as nx
import explore_local_structure as els

G = nx.Graph()
G.add_edges_from(
    [
        ("ego", "a", {"edge_type": 1}),
        ("ego", "b", {"edge_type": "1"}),
        ("a", "b", {"edge_type": 2}),
    ]
)
typed = els.typed_edge_induced_graphlet_degree_vector_ego(G, 2, nodes="ego")
assert typed["ego"][2] == [2, 1]
flat = [value for row in typed["ego"] for value in row]
assert len(flat) == 14

colored = els.colored_ego_graphlet_vector_for_typed_edge(G, 2, nodes="ego")["ego"]
names = els.graphlet_feature_names("colored", 2, ego=True)
assert colored[names.index((3, (1, 2)))] == 1

hetero = els.hetero_ego_graphlet_vector_for_typed_edge(G, 2, nodes="ego")["ego"]
names = els.graphlet_feature_names("hetero", 2, ego=True)
assert hetero[names.index((3, (1, 1, 2)))] == 1
```

Results always use node-keyed dictionaries, including a single `nodes` argument,
matching the original API. `nodes=None` selects all nodes; iterable selections
ignore unknown entries and deduplicate nodes. Mixed, non-orderable node labels
are supported. The input graph is never changed.

## Orbit and feature order

Full vectors use standard induced node orbits **0..14**. Ego rows correspond to
**0, 2, 3, 7, 11, 13, 14**, in that order:

| Ego row | Standard orbit | Graphlet / focal role |
|---|---|---|
| 0 | 0 | Edge / 2-clique |
| 1 | 2 | Center of a 2-path |
| 2 | 3 | Triangle |
| 3 | 7 | Center of a 3-star |
| 4 | 11 | Universal node of a tailed triangle |
| 5 | 13 | Degree-three node of a diamond / 4-chordal-cycle |
| 6 | 14 | 4-clique |

An ego occurrence has its focal node adjacent to all other vertices. The input
need not already be an extracted ego graph. Full functions also count graphlet
roles with vertices two or three hops from the focal node.

TyE columns represent types **1..T**; flatten matrices row by row, as in the
notebooks. `graphlet_feature_names(method, T, ego=...)` decodes flat coordinates
as `(standard_orbit_id, type_tuple)`. Orbit blocks are concatenated in order.
Colored type sets sort by cardinality, then lexicographically; heterogeneous
multisets follow lexicographic combinations-with-replacement order.

The original `get_type_from_index`, `get_type_from_index_for_ego_network` and
four edge-typed vector-length helpers are available. Invalid method names and
negative/out-of-range indices now raise `ValueError`.

## Induced and non-induced functions

`induced_graphlet_degree_vector` and `typed_edge_induced_graphlet_degree_vector`
count connected induced graphlets on two to four vertices using 15 orbits.
Their `_ego` versions restrict the focal role to the seven orbits above.

`graphlet_degree_vector_ego` and `typed_edge_graphlet_degree_vector_ego` are
non-induced: they retain all focal-to-alter edges and count every eligible
subset of present alter-to-alter edges. One four-clique therefore contributes
one 3-star, three tailed triangles, three diamonds and one four-clique at its
ego. Colored and heterogeneous functions are induced only.

## Input and type-0 conventions

- Undirected simple graphs only; directed and multigraph inputs raise
  `NetworkXNotImplemented`. Self-loops are ignored. Weights and node attributes
  do not enter the definitions. Isolates return correctly sized zero vectors.
- `num_type` must be a positive integer. Every non-loop edge must have an
  `edge_type` attribute; use `edge_type="relation"` to select another attribute.
- Integer labels and integer strings from the original CSV inputs are accepted.
  Missing attributes, out-of-range types, floats, booleans and malformed strings
  raise `ValueError`. No type labels are inferred.
- TyE functions accept **0..T**. Type 0 affects topology but does not increment
  a type column, preserving the original one-based implementation convention.
  It is not an extra counted type. Negative types, formerly ignored by parts
  of the legacy implementation, are now rejected explicitly.
- Full colored/heterogeneous functions require **1..T** everywhere. Their ego
  variants require positive types on edges incident to **each selected focal
  node**, but allow type 0 between alters. Select the intended ego explicitly
  when needed; selecting every node may violate this restriction.
- Ego flat lists retain the legacy zero-type blocks, including unused slots.
  Heterogeneous triangle/tailed-triangle/diamond/clique blocks allow at most
  1/1/2/3 zeros respectively, reflecting the possible untyped alter edges.

## Verification

The new engine enumerates rooted connected vertex sets, then classifies their
induced topology by degree sequence and focal degree. That classification is
sufficient for connected graphs on at most four vertices. Ego counting enumerates
one to three alters. An independent test oracle instead matches rooted templates
over vertex permutations and divides by root-preserving automorphism counts.
It checks all **729** four-node graphs whose possible edges are absent, type 1
or type 2. Golden fixtures preserve original vectors and coordinate decodings.

Repeat the legacy comparison from a repository checkout:

```shell
cd companion
python -m pip install -e ".[test,dev]"
python benchmarks/verify_graphlets.py --output graphlet-verification.json
```

The script reads original source at commit
`9454ce0da9fd6e9399e17dc552797add44e25933` via Git. **8,004 graph/method cases**
and **16 coordinate schemas** match exactly. They include exhaustive four-node
graphs, larger random graphs, type-0 alter configurations and the repository's
Chesapeake food web (39 nodes, 170 undirected edges).
The [recorded result](graphlet-verification.json) includes source/data checksums
and an example seven-row vector.

The food-web check uses deterministic **synthetic** three-type edge labels for
software verification. They are not observed social categories. The patient
ego-network files referenced by the notebooks are absent from the repo; the
paper directs readers to request data access from its data author. These
checks do not reproduce the patient-classification experiment.

## Cost and scope

Ego enumeration grows with combinations of one to three neighbors per focal
node. Full enumeration visits rooted connected sets up to size four, with
frontier-set operations; dense graphs can approach quartic work across all
nodes. This implementation prioritizes explicit, verifiable semantics and does
not claim the paper's optimized counting cost or a speed improvement over the
original neighbor-pattern code. Select focal nodes on large graphs.

TyE output stores `15*T` or `7*T` integers per selected node. Colored and
especially heterogeneous flat outputs, feature-name tables and lookup tables
grow combinatorially with T. Five-node graphlets, edge-orbit vectors, typed-node
APIs and the experimental `colored_ego..._v2` convention allowing type-0 focal
edges are outside this port.
