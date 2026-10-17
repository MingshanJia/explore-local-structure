# Upstream contribution preparation

The first target is official NetworkX main, audited at
`abf2feda2252ee4e0d5a66964addf179225dfe76`. Source is prepared as separate changes:

1. Overall endpoint closure and its arithmetic average.
2. Source/target and head/mid/end/cyclic directed patterns through an explicit
   pattern API and shared validated counts.
3. **Both i-quad and o-quad**, their nonnegative weighted versions, and averages.
   The i-quad docstring states its exact binary square-clustering transformation.

The companion-only combined/global/count helpers can evolve without increasing
the first upstream proposal's API burden. Naming and final API scope remain
maintainer decisions. These files are private engineering notes, not a public
issue/PR description to paste verbatim.

Before submission, verify no current issue/PR covers the proposal and establish
project fit. The maintainer's automated contribution policy requires personal
review, understanding and testing, contributor-written descriptions, and AI
disclosure. The repository's AGENTS.md additionally sets conditions on issue age
and labels for PRs targeting existing issues. Read them when selecting a target:

- [AGENTS.md](https://github.com/networkx/networkx/blob/abf2feda2252ee4e0d5a66964addf179225dfe76/AGENTS.md)
- [Contributor policy](https://networkx.org/documentation/latest/developer/contribute.html#automated-contributions-policy)

Scientific review should focus on the endpoint opportunity denominators,
reciprocal-arc multiplicity, zero conventions, published equation typo noted in
methods.md, and separation of signed/weighted-pattern work. Compatibility review
should focus on new public names, dispatch registration and the companion API's
documented departure from research list-valued returns.
