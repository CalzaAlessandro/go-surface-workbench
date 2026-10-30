# Fixed-geometry HF probes: execution versus stability

Seven closed-shell HF/6-31G* probes at independently accepted DFT geometries
completed Gaussian normally with exit code 0. All seven reported an external
RHF-to-UHF instability. The strict research worker therefore exited with code
4, while retaining the successful Gaussian logs and the scientific flag.
This is a deliberate additional research gate, not a Gaussian SCF execution
failure or automatic Amber incompatibility.

All 59 collected source files matched their recorded hashes. The HF probe
gates were repeated locally against the actual seven Gaussian logs and the
fixed DFT geometry records. Each had only hf_stability_not_confirmed; method,
charge, atom identity, geometry, finite energy and termination gates passed.
No HF optimization, frequency calculation, automatic UHF change or RESP fit
was performed by these probes.

HF_PANEL_SUMMARY.json and per_model_status.csv keep native Gaussian outcomes
separate from strict project acceptance. SOURCE_SNAPSHOT_MANIFEST.json binds
the source file identities. The archive SHA-256 is
2d1a2ed03a99f56ee941d1665908dda7f6199649f1a02d3192f4b470dd25d078.

See ../../docs/DFT_LEVEL_AND_RHF_STABILITY.md and
../../docs/AMBER_RESP_ACCEPTANCE.md. The conventional RHF charge definition
remains a separate comparison with its stability flag retained; a lower UHF
solution alone does not select a new charge/parameter model.
