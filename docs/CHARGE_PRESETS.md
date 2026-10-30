# Charge presets and the RESP-first study

The source catalog is `resources/charge_presets.json`. It records charge family,
coupled parameter expectations, provenance and readiness. It is a research
catalog; the general preparation CLI does not yet implement every assignment.

| ID | Charge model | Status |
|---|---|---|
| `resp-current` | Existing user RESP-derived fragment libraries | Preserved legacy reference; current assignment path exists |
| `resp-go-new` | New RESP-derived GO blocks/context response | Under derivation, no production coefficients selected |
| `gaff2-aim` | Published environment-averaged DDEC6 plus AIM LJ | Source preset archived unchanged; generic finite-GO mapping pending |
| `opls-hiergo` | Published HierGO/OPLS model | Original source preserved; generic assignment/Amber conversion needs validation |

GAFF2–AIM source topologies/coordinates at nominal C/O 1.8 and 3.0 are retained
with their associated force-field files and license under
`resources/reference/gaff2_aim/`. The source manifest verifies the original
download hashes. These fixed source sheets are parameter references, not
material-size presets: users always choose generated sheet dimensions.
The published basal OH/epoxide data do not establish finite-edge, carbonyl,
COOH/COO or graft-junction parameters. No missing group charge is invented and
no source coefficient is silently renormalized.

The new derivation remains RESP-based. First measure group and surrounding-carbon
responses on isolated/balanced and mixed fragments at controlled graph distances,
faces, protonation and sizes. Then compare a context table, a compact learned
rule, and a suitable published redistribution method recalibrated to RESP.
Existing NPA/AM1-BCC tools remain comparison controls rather than default new
charge derivations. Source presets are not refitted.

Current execution and scientific limits are recorded under
`research/parameterization/reports/CAMPAIGN_STATUS.json`. Experimental constraints
and the role of coupled LJ terms are discussed in
`research/parameterization/docs/EXPERIMENTAL_CHARGE_CONSTRAINTS.md` and
`EXPERIMENTAL_NONBONDED_BENCHMARKS.md` in that same directory. New bonded/frcmod
derivation is deferred under the present user scope.
