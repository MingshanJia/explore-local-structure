# Definitions and software conventions

These functions implement the cited binary definitions and nonnegative weighted
definitions. Chorded structures count; an “open” path means a closure opportunity,
including paths whose closing edge is already present, not an induced open motif.

## Triangle closure

For undirected adjacency A and degree d:

`E(i) = sum_jk A_ij A_jk A_ki / sum_j A_ij (d_j - 1)`.

For directed adjacency A, let M = A + Aᵀ. Reciprocal arcs contribute twice:

`E_D(i) = (M³)_ii / [2 sum_j M_ij sum_(k != i) M_jk]`.

Weighted total closure replaces A by W, where every eligible raw edge weight is
divided by the maximum eligible raw edge weight in the whole graph. Numerators
use products of three normalized weights; opportunity denominators use products
of two. This is not the geometric-mean weighted clustering coefficient.

## Directed patterns

Enumerate a two-edge path i–j–k with three distinct nodes, keeping actual arc
orientations separately, including reciprocal alternatives. For each path there
are two candidate directions for the closing i–k edge. Source/target distinguish
whether the candidate leaves/enters i. Head/mid/end/cyclic classify the resulting
directed triangle by i's role. Each coefficient is the number of present candidate
closures of its category divided by that category's candidate count. This oracle
is independent of the implementation's neighbor-intersection equations.

The invariants are `(source + target)/2 == total`; reversing arcs swaps source
and target and swaps head and end, while preserving mid/cyclic/total. Pattern
weights are unsupported rather than silently borrowing an unpublished formula.

The paper's source/target equations have an apparent focal/neighbor-degree
typographic inconsistency. The implementation follows the endpoint-opportunity
definition, surrounding prose, overall equation, and research source code;
exhaustive orientation enumeration verifies that interpretation.

## Quadrangles

For each focal i sum over pairwise distinct i,j,k,l, with absent edges assigned
weight zero and undirected symmetric normalized W:

```text
N_i = sum W_ij W_jk W_kl W_li
D_inner_i = sum W_ij W_jk W_il
D_outer_i = sum W_ij W_jk W_kl
I(i) = N_i / D_inner_i
O(i) = N_i / D_outer_i
```

Binary input uses edge indicators instead of weights. N_i counts each incident
four-cycle twice, including distinct cycles embedded in chorded subgraphs.

Let Q(i) be the number of incident four-cycles and D(i) the binary inner
opportunities. Current NetworkX 3.7 square clustering is `s = Q/(D-Q)` and
i-quad is `2Q/D`. Therefore **i-quad = 2s/(1+s)**. This equivalence is explicitly
part of the proposed API documentation. It is not an equality to the raw square
coefficient and does not extend to the weighted function.

Arithmetic averages of local coefficients differ from the global ratio of sums.
For a full graph, inner/outer opportunity sums agree by path reindexing.
The binary global ratio equals Robins–Alexander clustering on its documented
bipartite domain. Weighted global is an explicitly labelled derived extension
already present in the research code, rather than a separate published equation.

## Contracts

- One existing scalar node gives a native Python float (a pair for the combined
  quadrangle API); an iterable/None gives a dictionary. Unknown nodes in an
  iterable are ignored according to NetworkX nbunch conventions.
- Selecting nodes selects outputs only; all topology and the weight normalization
  maximum come from the full graph. Duplicate selected nodes count once.
- Loops are ignored, including their weight attributes when finding the maximum.
  Parallel edges are rejected. Quadrangles reject directed inputs.
- Finite nonnegative real weights are accepted; booleans act as 0/1. Missing
  requested attributes default to one. Negative, complex, nonnumeric and nonfinite
  values are rejected. All-zero weights yield zero.
- Undefined local/global coefficients and empty retained averaging sets give
  zero. Empty local dictionaries remain empty. These degenerate-case rules are
  explicit software choices, not additional scientific claims.
- Normalization and all calculations leave input topology and attributes intact,
  including frozen graphs and read-only views.

Ordinary weight ranges use floating-point products and accurate summation.
Extreme dynamic ranges use exact ratios of the supplied binary floating-point
weights until the final coefficient is converted to a Python float. This avoids
returning zero merely because both a numerator and denominator underflow before
division. That uncommon fallback trades speed for numerical stability.

## Signed scope

The directed paper extends to signed graphs but its printed weighted equation
and absolute-strength prose do not fully specify a raw signed normalization
convention. The research source normalizes in place using max(weight), which is
unsafe for all-negative inputs and mutates the original graph. This release
supports nonnegative total weights only. Any future signed extension must name
and test its max-absolute/magnitude-opportunity convention and reconcile it with
the published formula, rather than silently changing the method.
