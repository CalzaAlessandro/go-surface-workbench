<a id="nuovi-gruppi-go-riferimenti-di-ricerca-non-parametrizzati"></a>
# New GO groups: unparameterized research references

Seven neutral C/H/O references have been prepared, with no radicals in the graph and even electron counts. A closed-shell singlet is an initial assumption: no QM calculation, ground-state check, RESP fit, or MD test has been performed.

| Reference | Motif | Formula | SMILES source |
|---|---|---|---|
| `dibenzofuran_ether_edge` | non_epoxide_ether, 5-membered ring | C12H8O | [dibenzofuran](https://www.ebi.ac.uk/chebi/CHEBI:28145) |
| `diphenyl_ether_bridge_control` | non_epoxide_ether, ring not required | C12H10O | [diphenyl ether](https://www.ebi.ac.uk/chebi/CHEBI:39258) |
| `phthalide_lactone5_edge` | lactone, 5-membered ring | C8H6O2 | [phthalide](https://www.ebi.ac.uk/chebi/CHEBI:38085) |
| `coumarin_lactone6_edge` | lactone, 6-membered ring | C9H6O2 | [coumarin](https://www.ebi.ac.uk/chebi/CHEBI:28794) |
| `oxolan2ol_lactol_control` | lactol, 5-membered ring | C4H8O2 | [oxolan-2-ol](https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=200907060640520321) |
| `phenol_aromatic_OH_control` | phenolic_OH, ring not required | C6H6O | [phenol](https://www.ebi.ac.uk/chebi/CHEBI:15882) |
| `naphthoquinone14_edge` | quinone, 6-membered ring | C10H6O2 | [1,4-naphthoquinone](https://www.ebi.ac.uk/chebi/CHEBI:27418) |

Dibenzofuran is a small fused aromatic reference for an edge/defect ether environment; diphenyl ether is a diaryl bridge control. Phthalide, coumarin, and 1,4-naphthoquinone are fused aromatic references; phenol is a monocyclic aromatic control. Oxolan-2-ol is only an **aliphatic** control for lactol connectivity. None is an equivalent substitute for a complete GO sheet or evidence of the presence/abundance of these groups in our XPS samples.

The ChEBI sources were opened on 5 November 2026. For oxolan-2-ol, the SMILES and InChIKey were retrieved from the indexed J-GLOBAL/JST page; direct access returned HTTP 503. This access limitation is preserved in the source manifest. The graph identities of the seven inputs were also checked using the formula and RDKit InChIKey, before assigning a single geometric stereoisomer where the SMILES does not specify stereochemistry.

The gate distinguishes non-epoxide ether, five-/six-membered cyclic ester, cyclic hemiacetal without a carbonyl, neutral phenol, and conjugated para-quinone. In the models, each atom has a stable ID, element, integer formal charge, `extended_…` role, graph equivalence, and coordinates. **There are no partial charges.** The new roles are `UNPARAMETERIZED_RESEARCH_REFERENCE`; `modular.partition_model` and the builder are not declared compatible.

This second export restores the RDKit aromatic representation after MMFF94s, before assigning roles and writing JSON/SDF. The coumarin round-trip test checks that bond orders do not change between the two formats. The first export remains historical and must not be used as the current panel.

The seven `model.json`/`initial.sdf`/`initial.xyz` sets contain converged ETKDGv3 + MMFF94s geometries with seed 20261105. These are exclusively initial geometries. `geometry_opt/reference.gjf` prepares Opt/Freq/MK and Stable at wB97XD/6-311++G(d,p). `HF_RESP_template/reference.gjf` prepares an HF/6-31G* MK/Stable **template** with initial coordinates: it must be regenerated using an accepted DFT geometry before any scientific use. No orientations, submissions, or charges have been generated.

HierGO includes ETH for non-epoxide ether, but no explicit set of lactones or COO−. The lactones and lactol here are research extensions inspired by the literature, not established HierGO blocks. The GO rationale and attribution limits are in the local review `docs/MISSING_GROUPS_LITERATURE.md`; the identities of the small references do not establish the identities of sites in the material.

Transferable parameters require accepted QM, preservation of connectivity and atom order, fitting and validation on contexts excluded from calibration, definitions of nonoverlapping blocks/anchors, and consistent bonded/LJ/scaling conventions. Neutrality of subsets is not automatically imposed on future molecular RESP charges, and EPX parameters are not transferred to ethers, nor COOH parameters to lactones.

Transferability from approximately **4 nm to approximately 10 nm**, required for proteins/DNA, remains a requirement to be checked. The final method must preserve the same graph classes and molecular/group charge constraints across sizes. Checks held out from calibration are needed on larger sheets, varying edge/interior environments, low and high oxidation at a fixed size, different sizes at comparable composition, nearby groups, and protonation, with relevant electrostatic/MD controls. **New QM for the entire sheet is not required at every production size.** A new unsupported environment requires targeted references and a new coverage check.

The actual graph tests were written and observed to fail before implementation. Complete verification: 265 tests passed, including 20 in the extension. Charges, previous scientific inputs, and existing assignment code remain outside the extension.

To regenerate the panel, use a new directory from the PARAMETERIZATION directory:

```sh
PYTHONPATH=src uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python -m go_param.extended_fragments --output-dir reports/EXTENDED_GO_GROUPS_20261105_attempt_003 --cpus 8 --mem-gb 6
```

To rebuild only this README, the 2D panel, and the report:

```sh
uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python reports/EXTENDED_GO_GROUPS_20261105_attempt_002/build_report.py
```

See also `catalog.json`, `source_manifest.json`, `artifact_manifest.json`, and `docs/EXTENDED_GO_GROUPS_20261105.md`.
