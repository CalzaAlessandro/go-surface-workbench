# Research roadmap

These are proposed research directions. None is implemented by this snapshot.
Citing and reusing prior work properly is more important than adding cosmetic
features.

1. **Inverse XPS recipes with uncertainty.** Infer ensembles of surface
   compositions compatible with measured XPS peak areas, propagate peak-fit and
   sensitivity-factor uncertainty, and report non-identifiability rather than
   returning a single over-precise recipe.
2. **Chemical graph and reaction templates with provenance.** Represent each
   functional group and graft as explicit atom/bond graph edits with valence,
   charge, atom mapping, and parameter-source records; validate reversibility
   and conservation constraints.
3. **Controllable clustering, two-face decoration, and curvature.** Add spatial
   controls for clustered or dispersed oxidation, independent control of the two
   sheet faces, and curved or rippled substrates while preserving a consistent
   graph geometry.
4. **Modular force-field validation from `GO_PARAMETERIZATION`.** Define
   pluggable parameter backends and validate candidate parameter sets against
   matched QM fragments, conformers, and interaction energies before enabling
   new graft chemistries in MD exports.
5. **Ensembles and matched benchmarks.** Generate reproducible structural
   ensembles and compare methods on matched composition, protonation, geometry,
   boundary conditions, and reference data with uncertainty-aware summaries.
6. **Scalable incremental graph and spatial index.** Support large sheets with
   local graph updates, periodic neighbor indexing, and incremental collision
   checks, benchmarked for scaling and checked against the current reference
   behavior.
