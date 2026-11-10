<a id="audit-delle-cariche-e-riferimento-elettrostatico"></a>
# Charge audit and electrostatic reference

The audit provides a descriptive comparison of the sources. It does not select a winner:
GO CURRENT, the two OLD systems and the HierGO table do not describe the same
composition, protonation, connectivity or total charge. Differences in partial
charges, even large ones, do not by themselves measure the accuracy of a parameterization.

The provenance reported by the user for CURRENT is standard Amber parameterization
of separate GO fragments. The original QM method, basis set, fitting scheme, RESP
equivalences, geometries and constraints remain to be recovered. HF/6-31G*, RESP
or AM1-BCC are not retroactively attributed to the available tables.

<a id="sorgenti-e-risultati-dellaudit"></a>
## Audit sources and results

`CURRENT` includes only the nine `*_gaff2.lib` files in the
`resources/reference/current` copy. The nine `*_amber.lib` files in the same folder
are a separate family, `BUNDLED_LEGACY_AMBER`; `EHH.lib` is
`BUNDLED_AUXILIARY`. Sharing a folder does not imply that all files contain
parameters currently in use.

The OFF parser extracts the actual atom tables and internal connectivity without
depending on `boundbox` or the perturbation sections. This allows the audit
of the historical CBD and of EPH with unrelated EPX metadata, but does not verify
that these files are fully valid for LEaP. Atom identity uses the index within the template:
duplicate names are not merged.

OLD is read from the compiled topologies `gox.parm7` and `rgo.parm7` through ParmEd.
Chemical motifs are described by local connectivity, keeping residue-name
clues separate. Where reliable bond orders are unavailable, assignments
remain explicitly ambiguous, for example `carbonyl_or_deprotonated_oxygen`.
The neighbors in an incomplete template and those in a complete topology are not
automatically equivalent contexts.

HierGO is pinned to commit `60382a1a1807ef0003ffcc3903392e913b1e5e9a` of the
[official repository](https://github.com/IFM-molecular-simulation-group/HierGO).
The GPL-3.0 license, original files, URLs and SHA256 hashes are preserved in
`resources/reference/opls/PROVENANCE.json`; the audit checks these hashes.
The original `Oplsatom2` table contains 17 contextual rows, not 17 atoms of a
system. The same `opls_180` appears as epoxide oxygen EPX/OE with −0.4000 e
and ether oxygen ETH/OT with −0.1700 e: the OPLS type alone does not identify the charge.

Explicit selection of GO residues 1–494 in the available CURRENT complex,
excluding ligand 495, contains 1007 atoms (C613 O253 H141), with charge −16 e.
OLD gox contains 940 atoms (C661 O214 H65), with charge −9.0821659936 e; OLD rgo
contains 814 atoms (C674 O47 H93), with charge −3.0000199991 e. The noninteger charge
of OLD gox warrants recovery of the original procedure; it is not automatically
corrected and does not justify a conclusion about relative accuracy.

```python
from go_builder.charges import audit_charge_sources

report = audit_charge_sources(
    "results/charge_audit",
    current_topology="/path/to/compl.parm7",
    current_go_residue_count=494,
)
```

API: `audit_charge_sources(output_dir, *, current_dir=None, old_dir=None,
opls_dir=None, current_topology=None, current_go_residue_count=494)`.
The default source directories are the copies under `resources/reference`.
The compiled CURRENT topology is optional; when supplied, selection
of the first N residues must be checked against the provided system.

Output: `audit.json`, `atoms.csv`, `motifs.csv`, `sources.csv`. Statistics are
separated by source and atomic context; sums over templates or parameter tables
are not total surface charges. Each data file has a SHA256 hash and path.
An existing output directory raises `FileExistsError`.

<a id="preparazione-di-un-confronto-verificabile"></a>
## Preparing a verifiable comparison

1. Choose finite GO fragments, with edges explicitly capped with H, that
   represent graphene edges, OH, epoxide, carbonyl, COOH, COO− and any
   ring-opening product. Define protonation, integer charge and multiplicity;
   some basal-group counts produce species with unpaired electrons.
2. Apply each source only where a documented chemical correspondence
   exists. Charges for capping atoms and new junctions require a justified
   choice. If a source does not cover the fragment, record the case as
   unassignable; do not add zeros or neutralize charges to obtain a comparison.
3. Use the same geometry and charge state for all candidates
   for each fragment. Prepare a declared QM reference and repeat the
   comparison across multiple conformations; the reference does not reconstruct
   the original CURRENT method.
4. Compare ESP over the same set of points and dipoles relative to the same
   origin. Subsequently check geometry, stability in Amber and sensitivity
   of observables, including MM/GBSA with comparable GB radii and sampling.

```python
from go_builder.qm import prepare_qm

manifest = prepare_qm(
    capped_model, "results/qm_epoxide", total_charge=0,
    multiplicity=1, method="HF", basis="6-31G*", nproc=4, memory="4GB",
)
```

API: `prepare_qm(model, output_dir, *, total_charge=None, multiplicity=1,
method='HF', basis='6-31G*', nproc=4, memory='4GB')`. If `total_charge` is absent,
formal charges must be present for all atoms. Coordinates are in Å.
Periodic models, bonds with remaining image offsets, open edges,
atomic overlaps and incompatible charge/spin combinations are rejected.
Periodic fragments are not automatically cut or capped.

The files `reference.gjf`, `qm_manifest.json`, `resp-stage1.in`,
`prepare_resp.sh` and instructions are written. The Gaussian file specifies a single point,
without optimization, with `Pop=MK`, tight SCF and ESP point output.
The default HF/6-31G* reference is an explicit benchmark choice.
Gaussian and AmberTools remain external software: **no calculation is launched**.
The user can run the RESP script after checking the QM output;
it refuses to overwrite results. The supplied one-stage RESP fit is
exploratory, without equivalences: production fitting requires chemical constraints
and any second stage to be defined for the chosen protocol.

<a id="da-espgen-al-confronto-comune"></a>
## From espgen to a common comparison

```python
from go_builder.qm import espgen_to_reference, compare_esp

espgen_to_reference(
    "results/qm_epoxide/reference.esp", capped_model,
    "results/epoxide_reference.json", total_charge=0,
    # Optional: QM dipole in the same orientation and with the same origin.
    dipole=[0.0, 0.0, 1.2], dipole_unit="debye",
    dipole_origin=[0.0, 0.0, 0.0],
)
compare_esp(
    "results/epoxide_reference.json",
    {"CURRENT": current_matched_model, "OLD": old_matched_model,
     "HIERGO": hiergo_matched_model},
    "results/epoxide_esp",
)
```

`espgen_to_reference(esp_path, model, output_path, *, total_charge,
geometry_tolerance=1e-5, dipole=None, dipole_unit='debye', dipole_origin=None)`
reads standard Amber ESP for one molecule, with coordinates in bohr and potential
in atomic units. The geometry must match the model's atom order
within the tolerance in Å. `NoSymm` preserves the required orientation;
if atoms have been rotated or reordered, the converter rejects the data and requires
an explicit mapping. `dipole_origin`, when supplied to the converter, is in Å.

The minimal JSON structure for an example analytical reference is:

```json
{
  "atoms": [{"id": 7, "element": "H", "position": [0, 0, 0]}],
  "total_charge": 1,
  "position_unit": "bohr",
  "potential_unit": "hartree/e",
  "grid": [{"position": [2, 0, 0], "potential": 0.5}],
  "dipole": [0, 0, 0],
  "dipole_unit": "debye",
  "dipole_origin": [0, 0, 0]
}
```

A candidate file for the ESP CLI uses this form, with coordinates in Å by default:

```json
{
  "analytical_unit_charge": {
    "atoms": [{"id": 7, "element": "H", "position": [0, 0, 0], "charge": 1}]
  }
}
```

This example checks units and Coulomb's law; it does not represent a completed
QM benchmark or a GO model. Real data require the corresponding charges
for all atoms, including capping atoms.

API: `compare_esp(reference_path, candidates, output_dir, *,
geometry_tolerance=1e-6, charge_tolerance=1e-6)`. Identity,
elements, geometry, finiteness, nonperiodicity and equality of total charge are checked.
No charge renormalization or compensation is applied.
The reference requires explicit units; it supports bohr/Å coordinates and dipoles in
debye, e·bohr or e·Å. In the reference JSON, the origin uses `position_unit`.

Output: `esp_comparison.json` and `esp_residuals.csv`, with RMSE, MAE, maximum error,
relative RMSE, dipoles and their respective errors if the QM dipole is provided. Potentials are in
Hartree/e; conversion to kcal/mol/e gives the energy per unit test charge.
The dipole of a charged system depends on the origin. The RMSE values are unweighted:
the distribution of reference points determines the spatial weighting.
Hashes of the reference and each candidate charge/geometry set allow
the exact comparison to be reconstructed. No overall recommendation is produced.

Primary operational references: the espgen format, including bohr coordinates and
fixed-width fields, is verified in the
[official Amber source](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/espgen.c).
The use of espgen and respgen for two-stage fitting is described in the
[Amber manual](https://dev-archive.ambermd.org/200603/att-0122/amber9.pdf).
