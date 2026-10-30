# Material presets, charge choices and external GO preparation

Design update from the user's stated distribution goals. This document does
not claim a new importer, charge-correction CLI or final parameter release.

## Three separate choices with compatibility checks

- **Material recipe:** four LayerOne examples/quick creations plus custom
  geometry and group counts. A preset records its XPS/source assumptions,
  reference dimensions, composition and protonation. It is not a universal
  unique atomistic reconstruction of an XPS spectrum. The existing provisional
  presets and the unresolved dense P935 placement remain explicit; definitive
  rebuilding is deferred until parameter selection.
- **Charge assignment:** the standard reusable block library is the default.
  An optional, separately recorded structure-specific step can redistribute
  those charges after construction. Current, OLD and new/published controls
  must retain their status and provenance; diagnostic pilot values are not
  production presets.
- **Force-field parameters:** an explicit bonded/LJ/1–4/radii choice with a
  matrix of supported charge combinations. An OPLS or GAFF2–AIM charge table
  does not automatically define a valid variant of the current GAFF2 model.

Charge assignment and structure-specific correction are distinct settings:
choosing the reference charge library does not imply a correction ran. Preserve
both vectors and their exact atom map for review and reproducibility.

## One preparation engine, two entry paths

The native generator and a separate user-facing external-preparation script
should converge on the same explicit chemical graph and assignment engine.
Existing `import_Amber.py` already imports Amber topology/coordinates; that
does not establish support for generic coordinate-only PDB files.

External sources include GOPY, HierGO and online structures. Prefer chemical
formats with explicit connectivity and bond orders or an existing topology.
A coordinate-only PDB may lack bonds, bond orders, protonation and periodic
links. It cannot be treated as an unambiguous chemical model. Any reconstruction
must be checked against valence, ring context and boundary conditions, with
ambiguities reported instead of silently applying bond-by-distance.

The preparation stage should:

1. Preserve the original file, identify atoms and record connectivity/PBC
   evidence. Keep stable source-to-output atom identifiers.
2. Verify hydrogen/protonation and recognize groups by their chemical graph,
   independently of source residue names. Distinguish EPX/EPH/EPT, basal and
   phenolic OH, carbonyl, COOH/COO and other supported motifs.
3. Assign the selected standard library and force-field terms. Detect
   overlapping blocks, unsupported groups and incomplete cross-environment
   bonded terms. Graph context must distinguish an actual carbon triangle
   from a non-ring C–C–C angle between epoxidized carbons.
4. Optionally apply the separately selected static correction method, without
   overwriting shared global library files. Record total and per-block charge
   changes; global conservation alone does not guarantee neutral blocks.
5. Verify complete topology/parameter coverage and consistent coordinates,
   then export Amber topology/restart/PDB and the atom/group map. Retain the
   existing GROMACS conversion stage for supported complete topologies.

Expected review outputs are the normalized graph, source map, baseline/final
per-atom charges, assignment provenance and parameter/charge checks. Raw inputs
and baseline values remain available. Unknown lactone/lactol/ether motifs do
not get relabeled as the nearest existing group to suppress an error.

MD and GBSA use the same fixed final charges; subset topologies inherit them.
Functionalization chemistry, protonation or charge response must not be
re-inferred independently after stripping a ligand for an energy analysis.

## Published correction candidate

The current first candidate for evaluation is the SQE+qp family implemented in
ChargeFW2, because it includes reference type charges and a bond-based response.
See [method selection](STATIC_CHARGE_METHOD_SELECTION.md) for sources and the
important NPA-versus-RESP and GO-domain limitations. This is a recommendation
for testing, not selection of a validated GO parameter set.
