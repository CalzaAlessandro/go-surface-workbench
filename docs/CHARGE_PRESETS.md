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

Published presets are complete, indivisible force-field bundles. Their charges,
bonded terms, LJ parameters, mixing rules and 1–4 conventions must be preserved
together. Selecting `gaff2-aim` or `opls-hiergo` does not permit a charges-only
override or an SQE/Espaloma correction under the same published name. A new
derivative would require a separate identity, provenance and validation.
For HierGO/OPLS, COO− is not covered by the archived charge table. No missing
chemistry is silently borrowed from another force field. Amber export also
requires a validated conversion and an explicit compatible GB-radius policy.

The assignment-method catalog distinguishes three research paths:

| Method | Current readiness |
|---|---|
| Fixed library | Legacy user library available; recalibrated modular HF–RESP is a limited-context candidate |
| Graph context | Minimum graph-distance/context approach under calibration; no final table/rule selected |
| Static redistribution | First SQE role-coefficient fit to molecular HF–RESP completed on 16 contexts; no production approval. Published Espaloma remains a BCC-oriented control |

All final charges remain fixed during MD. Structure-specific response is a
preparation step and must preserve the chemically determined integer total
charge, which can be nonzero. The catalog describes readiness and compatibility;
it does not implement a generic production selector or all published mappings.

For the new Amber-oriented preset, the canonical reference target is gas-phase
HF/6-31G* MK ESP with Amber RESP fitting. Current DFT–RESP fits are separate
comparison controls. The mixed-level controller has produced 16 conventional
HF–RESP references on DFT geometries. It checks the accepted DFT reference,
runs three rigid HF MK orientations, and calls external Amber RESP. These
rotations are not three independent chemical environments.
External RHF-to-UHF instability is retained as a scientific flag and is not
automatic Amber incompatibility. The [Amber RESP acceptance policy](../research/parameterization/docs/AMBER_RESP_ACCEPTANCE.md)
records this distinction and the outstanding force-field validation.

Current execution and scientific limits are recorded under
`research/parameterization/reports/CAMPAIGN_STATUS.json`. Experimental constraints
and the role of coupled LJ terms are discussed in
`research/parameterization/docs/EXPERIMENTAL_CHARGE_CONSTRAINTS.md` and
`EXPERIMENTAL_NONBONDED_BENCHMARKS.md` in that same directory. New bonded/frcmod
derivation is deferred under the present user scope.
