# Individual benchmark cases

All declared cases are retained. The overview bars show medians for QM/charge comparisons and means for pattern separation. No confidence intervals are inferred.

## Force RMSE, kcal/(mol·Å)

| Combination | Pyrene | Diol cis | Diol trans | Diol + EPX | COOH | COO⁻ |
|---|---:|---:|---:|---:|---:|---:|
| GAFF 1.81 + AM1-BCC | 20.975 | 52.391 | 52.430 | 17.881 | 20.258 | 21.970 |
| GAFF 1.81 + HF-RESP | 21.055 | 52.402 | 52.436 | 18.788 | 20.459 | 21.833 |
| GAFF2 2.2.20 + AM1-BCC | 16.653 | 39.457 | 39.508 | 14.057 | 16.184 | 17.888 |
| GAFF2 2.2.20 + HF-RESP | 16.747 | 39.475 | 39.533 | 15.120 | 16.425 | 17.793 |
| OpenFF 2.0 + AM1-BCC | 15.741 | 12.230 | 12.545 | 16.796 | 16.400 | 22.346 |
| Sage 2.3 + AshGC | 15.653 | 10.498 | 10.885 | 12.621 | 17.111 | 19.023 |
| Grappa 1.4 + AM1-BCC | 5.538 | 18.211 | 17.815 | 6.240 | 7.931 | 8.112 |
| Espaloma 0.3.2 + AM1-BCC | 9.772 | 8.282 | 7.765 | 7.975 | 9.795 | 9.269 |
| Espaloma 0.3.2 + EspalomaCharge | 9.803 | 8.390 | 7.855 | 7.962 | 9.917 | 9.412 |
| Q-Force (Hessian-only, in-sample)* | 2.097 | 5.763 | 4.374 | 5.692 | 3.482 | 2.400 |

`*` Q-Force (Hessian-only, in-sample) is fitted directly to each molecule's QM Hessian without torsion scans, using GAFF2 LJ/AM1-BCC.

## Cartesian-Hessian error, %

| Combination | Pyrene | Diol cis | Diol trans | Diol + EPX | COOH | COO⁻ |
|---|---:|---:|---:|---:|---:|---:|
| GAFF 1.81 + AM1-BCC | 24.973 | 31.555 | 31.500 | 23.984 | 25.928 | 25.992 |
| GAFF 1.81 + HF-RESP | 24.970 | 31.560 | 31.499 | 23.999 | 25.932 | 26.002 |
| GAFF2 2.2.20 + AM1-BCC | 18.002 | 25.543 | 25.314 | 22.138 | 19.705 | 18.929 |
| GAFF2 2.2.20 + HF-RESP | 18.002 | 25.546 | 25.316 | 22.161 | 19.707 | 18.936 |
| OpenFF 2.0 + AM1-BCC | 18.515 | 21.326 | 20.732 | 23.579 | 20.370 | 21.106 |
| Sage 2.3 + AshGC | 18.424 | 20.301 | 20.137 | 18.244 | 18.098 | 18.685 |
| Grappa 1.4 + AM1-BCC | 17.356 | 20.652 | 21.263 | 16.131 | 16.825 | 17.262 |
| Espaloma 0.3.2 + AM1-BCC | 17.295 | 45.937 | 45.077 | 21.605 | 16.691 | 19.798 |
| Espaloma 0.3.2 + EspalomaCharge | 17.301 | 45.949 | 45.087 | 21.607 | 16.697 | 19.796 |
| Q-Force (Hessian-only, in-sample)* | 13.701 | 16.947 | 17.635 | 9.590 | 13.419 | 13.808 |

`*` Q-Force (Hessian-only, in-sample) is fitted directly to each molecule's QM Hessian without torsion scans, using GAFF2 LJ/AM1-BCC.

## C54 ESP relative RMSE, %

Median of three orientations within each molecule; molecular RESP was fitted on the same grids. HierGO OPLS covers the five neutral C54 molecules only (`COO⁻` absent); `HierGO_OPLS_COO_HF_ESP_extension` adds a 2-parameter `COO⁻` block fitted on 10 non-C54 `COO⁻` fragments.

| Assignment | EPX2 | EPX4 | OH2 | OH4 | OH6 | COO⁻ |
|---|---:|---:|---:|---:|---:|---:|
| L0_strict40 | 16.58 | 18.06 | 14.80 | 21.44 | 15.55 | 16.52 |
| Local_BCC_ESP_strict40 | 14.11 | 17.30 | 12.89 | 22.50 | 14.48 | 13.39 |
| SQE_ESP_strict40 | 16.06 | 17.87 | 15.57 | 23.44 | 15.98 | 12.19 |
| HierGO_OPLS_COO_HF_ESP_extension | 20.99 | 22.13 | 16.76 | 28.55 | 22.85 | 16.62 |
| hiergo_opls | 20.99 | 22.13 | 16.76 | 28.55 | 22.85 | not covered |
| current_gaff2 | 32.48 | 27.25 | 33.71 | 37.51 | 32.30 | 18.96 |
| AshGC | 84.03 | 59.13 | 75.29 | 84.49 | 72.91 | 25.76 |
| EspalomaCharge | 202.40 | 163.24 | 301.25 | 220.93 | 209.27 | 43.00 |
| Atomic_RESP_in_sample_control | 9.02 | 8.41 | 8.28 | 7.88 | 6.04 | 1.17 |

## Published charge-table coverage and charge-conservation audit on C54

| Table / extension | Covered C54 molecules | Matched formal charge | Missing roles / scope limitation |
|---|---:|---:|---|
| current_gaff2 | 6/6 | 6/6 | Complete C54 role coverage |
| hiergo_opls | 5/6 | 5/6 | Missing: COO_C, COO_O, COO_anchor |
| HierGO_OPLS_COO_HF_ESP_extension | 6/6 | 6/6 | Charge-only extension; OPLS-AA lacks aromatic `CA-C_3-O2` / `CA-CA-C_3-O2` bonded terms |
| gopy_raw | 5/6 | 0/6 | Nonzero neutral motif sums (-1.476 to -0.454 e); missing: COO_C, COO_O, COO_anchor |
| stauffer_esp | 0/6 | 0/6 | Missing: COO_C, COO_O, COO_anchor, edge_C, edge_H, scaffold_C |
| stauffer_sesp | 0/6 | 0/6 | Missing: COO_C, COO_O, COO_anchor, edge_C, edge_H, scaffold_C |
| pinto_gaff2_aim | 0/6 | 0/6 | Periodic DDEC6 + AIM LJ (`GO_18`/`GO_30`); no automatic finite-fragment transfer |

## Pattern separation across 16 seeds

Distance is the shortest C-scaffold path between group anchor sets, averaged over group pairs.

| Pattern | Mean across seeds | Minimum seed value | Maximum seed value |
|---|---:|---:|---:|
| random | 7.583 | 5.643 | 9.750 |
| island | 4.337 | 3.571 | 4.929 |
| dispersed | 9.444 | 8.071 | 10.286 |
