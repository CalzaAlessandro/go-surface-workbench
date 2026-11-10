<a id="frammenti-di-riferimento-per-gruppi-go-aggiuntivi"></a>
# Reference fragments for additional GO groups

This extension prepares seven neutral C/H/O molecular graphs and reproducible initial geometries for future QM comparisons. It adds no types or charges to the builder, does not modify L0 and runs neither QM nor molecular dynamics. All roles start with `extended_` and have status `UNPARAMETERIZED_RESEARCH_REFERENCE`. The format is compatible with `go_param.resp._validate_model` and the Gaussian input helpers; this does not mean that `modular.partition_model` supports the new roles.

| Model | Motif | Interpretation of the reference |
|---|---|---|
| `dibenzofuran_ether_edge` | O in a five-membered ring, not an epoxide | Small fused aromatic motif for a possible edge/defect ether environment; heteroaromatic O is not equivalent to every aliphatic ether |
| `diphenyl_ether_bridge_control` | Noncyclic Ar–O–Ar | Separate diaryl-bridge control, not a contiguous sheet fragment |
| `phthalide_lactone5_edge` | Five-membered cyclic ester | Fused aromatic lactone reference; includes the carbonyl and ester O of the same ring |
| `coumarin_lactone6_edge` | Six-membered conjugated cyclic ester | Fused aromatic reference; charges are not automatically transferred from the five-membered lactone |
| `oxolan2ol_lactol_control` | sp³ C with OH and cyclic OR | Aliphatic cyclic-hemiacetal control; not a lactone and does not establish an aromatic edge environment |
| `phenol_aromatic_OH_control` | Aromatic C–OH | Neutral phenol control, distinct from basal sp³ alcohol and phenolate |
| `naphthoquinone14_edge` | Two conjugated para carbonyls | Fused aromatic reference for a quinone context broader than a single carbonyl |

These molecules are connectivity references, not equivalent replacements for complete sheets. No presence or abundance is assigned to the user's XPS samples. A neutral graph with zero RDKit radical electrons and an even electron count permits an initial closed-shell singlet hypothesis, but does not validate the ground state or determinant stability.

<a id="rapporto-con-hiergo-e-letteratura"></a>
## Relationship to HierGO and the literature

The verified local inventory and [review of missing groups](MISSING_GROUPS_LITERATURE.md) distinguish HierGO ETH, a nonepoxide ether, from EPX/EPH/EPT. HierGO does not provide an explicit lactone or COO− set: the lactones prepared here are a research extension motivated by GO literature, rather than an established HierGO implementation. The lactol also remains a separate candidate.

The local review links the GO motivation to Gao et al., DOI [10.1038/nchem.281](https://doi.org/10.1038/nchem.281), Bagri et al., DOI [10.1038/nchem.686](https://doi.org/10.1038/nchem.686), and Aliyev et al., DOI [10.3390/nano9081180](https://doi.org/10.3390/nano9081180), with their respective material conditions and assignment limitations. The identities of the small references do not establish experimental GO site identities. SMILES sources and their identifiers are recorded in the catalog and report README.

<a id="protocollo-e-gate"></a>
## Protocol and acceptance gates

1. RDKit sanitization and verification of the formula and InChIKey of the documented graph. Explicit addition of H atoms, total formal charge 0, zero graph radicals, an even electron count, contiguous atom IDs and graph-equivalence classes without atom maps.
2. Motif-specific gates: ether O belongs neither to a triangle nor to the ester linkage of a carbonyl; in a lactone, the carbonyl C and ester O belong to the required ring; in a lactol, a cyclic sp³ C bonds to OH and OR and the graph contains no carbonyls; phenolic OH is on aromatic C; a quinone has two para carbonyl sites and conjugation in the same six-membered ring.
3. ETKDGv3 with seed 20261105 and one thread, followed by MMFF94s to convergence (maximum 3000 iterations). These are initial geometries, not accepted QM geometries. Where the SMILES leaves stereochemistry undefined, one deterministic stereoisomer is recorded from the geometry, without assigning populations or racemates.
4. `geometry_opt/reference.gjf`: existing helper for optimization, frequencies and stability at wB97XD/6-311++G(d,p). `HF_RESP_template/reference.gjf`: an HF/6-31G* MK/Stable input template at fixed geometry only. **The HF template still contains MMFF coordinates and must be regenerated from the accepted DFT geometry before scientific use.** No orientation or submission is prepared automatically.
5. Before fitting, actual optimization/frequency/stability evidence and checks that connectivity, bond order, stereochemistry and atom identity are preserved are required. The subsequent RESP fit conserves molecular charge; it does not invent atomic charges or subset neutrality in advance.
6. Before any production assignment, disjoint block and anchor ownership, an explicit local environment, validation on contexts excluded from fitting and charges coupled to consistent bonded/LJ/scaling terms are required. Torsions and barriers, equilibrium geometries, nonbonded interactions and relevant MD properties must be checked. EPX parameters are not automatically reused for ETH, HYD for phenol or COOH for lactone.

<a id="trasferimento-da-circa-4-nm-a-circa-10-nm"></a>
## Transfer from approximately 4 nm to approximately 10 nm

Size transferability is a future requirement of the final charge method, including for protein/DNA simulations. These small references and their inputs do not demonstrate it. Once validated, assignment must reuse the same local graph classes and total/group charge constraints at both sizes, without introducing dependence on the sheet name or total atom count.

Validation excluded from calibration is required on larger sheets, including low and high oxidation, edge/interior variation, group combinations and separations, protonation and geometry. Size and oxidation level must be varied separately: the same counts on a larger sheet do not imply the same coverage. Exact charge sums, class coverage, absence of new unsupported environments and stable assignments for equivalent environments must be checked; relevant electrostatic and MD properties must be compared with available independent references. The protocol does not require a new whole-sheet QM calculation for every production size. A new environment or failed validation requires targeted reference expansion, rather than automatically treating the small fragment as evidence of transferability. See the [acceptance criteria](SCALABILITY_ACCEPTANCE_20261105.md).

<a id="implementazione-e-verifica"></a>
## Implementation and verification

Isolated API: `go_param.extended_fragments.build_candidate`, `validate_motif` and `export_candidates`. Export rejects an existing path and builds/validates the entire panel before writing files. Tests were written before implementation and cover real incorrect graphs, SDF/ID reconstruction, geometric determinism and template preparation without executing QM.

For a new export, from the `PARAMETERIZATION` directory:

```sh
PYTHONPATH=src uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python -m go_param.extended_fragments --output-dir reports/EXTENDED_GO_GROUPS_20261105_attempt_003 --cpus 8 --mem-gb 6
```

Use a new directory for each attempt. The current panel is `attempt_002`, with consistent aromatic representation between JSON and SDF after MMFF94s. `attempt_001` preserves the first historical export, superseded because of inconsistent coumarin representation; it contains no QM results.
