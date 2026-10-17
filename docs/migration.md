# Migrating from the research snapshot

Install this package alongside standard NetworkX. Import NetworkX as `nx` and
this package as `els`. Do not uninstall or shadow the standard `networkx` package.

| Research API | Modern API |
|---|---|
| `nx.closure(G)` with list-valued entries | `els.closure(G)` with scalar total values |
| Directed closure entry `[total, source, target]` | Separate `els.directed_closure(..., pattern=...)` calls |
| `src_closure`, `tgt_closure` | Patterns `source`, `target` |
| `head_closure`, `mid_closure`, `end_closure`, `cyc_closure` | Patterns `head`, `mid`, `end`, `cyclic` |
| `average_closure` | `els.average_closure`; averaging empty selections now returns zero |
| `inner_quadrangle_coefficient`, `outer_quadrangle_coefficient` | Same names under `els` |
| `iquad_oquad_coefs` | `els.quadrangle_coefficients`, returning pairs instead of lists |
| `average_inner_quad_co`, `average_outer_quad_co` | `average_inner_quadrangle_coefficient`, `average_outer_quadrangle_coefficient` |
| `global_quadrangle` | `global_quadrangle_coefficient`, full graph only, zero denominator gives zero |
| `number_of_quadrangles` | Same name under `els`, returning an integer |

Experimental counting iterators, joint clustering/closure pattern routines,
alternative grid/square measures, and graphlet/colored/heterogeneous baselines
are not advertised as supported APIs in 0.1.0. Use the archived research snapshot
to reproduce the original studies; the modern branch does not overwrite it.

Weighted non-total closure used to reach incomplete branches and raise
`UnboundLocalError`; it now raises an explicit `NetworkXNotImplemented`.
Negative/nonfinite weights are rejected before normalization. Loop weights do
not change normalization. Unlike the legacy directed weighted function, no
function edits graph edge weights. Published nonnegative product formulas remain
the basis of the weighted total and quadrangle functions.

The old `networkx` API namespace and list-valued return shape are deliberately
not emulated. This is a new companion-package API, not a drop-in NetworkX fork.
