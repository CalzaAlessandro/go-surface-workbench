# Bonded QM reference branch

Eight finite fragment geometries are independently optimized with Gaussian
ωB97X-D/6-311++G(d,p), an ultrafine DFT grid and tight SCF/geometry convergence.
Analytic harmonic frequencies and a separate fixed-geometry wavefunction
stability check are retained. The harmonic checkpoint is preserved separately
from the stability checkpoint. Both formatted checkpoints are exported;
Cartesian force constants must be present before the job records technical
completion. No vibrational scaling or force-field fitting is applied by this
worker, and completion still needs topology, stereochemistry, frequencies,
spin/electronic-state and Hessian QA.

This is a finite-molecule bonded-reference strategy, not a claim that it is a
validated periodic GO model. The choice is supported as a molecular force-field
benchmark by the earlier QUBE workflow and [QM-to-MM protocol comparisons](https://pubs.rsc.org/en/content/articlepdf/2022/cp/d2cp02864f).
[Allen, Payne and Cole (2018)](https://doi.org/10.1021/acs.jctc.7b00785)
provide the modified-Seminario approach for bond/angle extraction. Any future
use must verify implementation, symmetry/context pooling, AMBER harmonic
normalization and transferability. Torsions need suitable scans and the
intended nonbonded model; LJ terms cannot be inferred from the Hessian alone.

These calculations are independent from the conventional HF/MK/RESP charge
controls and do not silently replace them. All candidates and active MD
parameters remain distinct until comparison is complete.
