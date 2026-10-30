# Seven completed DFT–RESP comparison controls

All seven accepted-parent controls completed their three-orientation MK and
Amber RESP calculations. The independent postanalysis verified the actual
SCF-method/energy evidence and original fit artifact hashes. The 62-file
metadata/context snapshot was hash-verified locally and its context tables
were regenerated from the seven original atomic charge files.

Method: gas-phase omegaB97X-D/6-311++G(d,p) optimization, frequencies and
stability, followed by MK potentials at the same DFT level and Amber RESP.
There are seven physical contexts and 21 numerical orientations, not 21
independent chemical contexts. See DFT_RESP_CONTEXT_COMPLETE.json and
SOURCE_SNAPSHOT_MANIFEST.json for the evidence chain.

training_ESP_comparisons.csv retains all 210 candidate/orientation rows,
including unsupported coverage and mismatched-charge exclusions. New fitted
charges were trained on those ESP grids; those residuals do not validate
transfer, adsorption or a complete force field. Whole-sheet OLD and periodic
GAFF2–AIM charges remain separate references without invented finite-fragment
mappings. No published source coefficients were changed.

These are DFT–RESP controls, not the canonical new Amber-oriented RESP preset.
HF/6-31G* electronic-state probes at independently accepted DFT geometries
are a separate calculation branch. A stable HF probe is not itself a RESP
fit or a validated force field. See ../../docs/AMBER_RESP_ACCEPTANCE.md.

The source archive SHA-256 is
84a8a6024b0e23dbec09bb34a4b5cff6aeb18b0e7b85200b8531e0f1c06c5b29.
Raw calculation archives remain outside Git. dataset.json in the local study
contains machine-specific paths; regenerate it after relocating source files.
