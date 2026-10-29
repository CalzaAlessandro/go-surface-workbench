# Interim results, 29 October 2026

A real system-Gaussian / AmberTools pilot has completed, independently from
fixture tests. Pyrene HF/6-31G* optimization has 72 positive harmonic modes,
minimum 106.4879 cm⁻¹, but an RHF→UHF orbital instability. Strict reference
acceptance remains false. Three rotated fixed-geometry MK jobs and actual
Amber RESP produced a total charge of 0 e and verified symmetry/identity and
21 artifact hashes. This is a restricted-HF diagnostic atomic reference, not
final transferable GO parameters.

A separate UHF singlet calculation at the same RHF geometry lowered the
energy by 17.6289 kcal/mol and passed its stability check, but its determinant
has ⟨S²⟩ = 1.8279 before annihilation (a pure singlet would be 0); after
annihilation the reported value is 7.1070. This spin-contaminated solution
must not be silently adopted as the physical singlet or a final charge model.
The UHF calculation did not optimize nuclear geometry.

ESP-grid residuals are diagnostic metrics on training data. They do not
establish transferable charges or adsorption accuracy. Charge-method
comparisons include RESP and newly requested AM1-BCC, alongside published
controls; full GAFF/GAFF2/OPLS/GO-specific model comparisons require their own
bonded/LJ/1–4 terms and matched external validation.

See the adjacent report JSON files for the numerical diagnostic evidence.
Raw logs/checkpoints remain in the isolated local/Slurm study archive.
