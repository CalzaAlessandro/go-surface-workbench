# Dione: native AM1-BCC fallback diagnostic

The actual `pyrene_dione/attempt_003` remains **excluded for incomplete BCC parameter coverage**. Its existing failed marker and all input/native-output artifacts remain unchanged. No production policy was relaxed, no new calculation was run, and no charge candidate was exported.

Native antechamber exited successfully and SQM AM1 geometry converged. Atom names, IDs, elements, connectivity and constrained resonance representations pass inspection. Of 29 native bonds, 28 have explicit parameter-table matches. The only gap is **C013–C014 (atoms 13–14), BCC types 16–12, double-bond code 2**. The native AC and MOL2 both retain this double bond.

The public [Amber `bccparm` implementation](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/am1bcc.c#L113-L154) returns zero after direct/reversed and optional CORR parameter lookups fail; the exact short source context is recorded in the JSON report. This documents a native fallback, not a newly fitted zero parameter. The installed binary's output independently agrees with that behavior. Public-source inspection does not establish byte-identical installed source.

I independently reconstructed every bond increment from the actual snapshotted `BCCPARM.DAT`, allowing **only the documented diagnostic zero for that single unmatched bond**. The resulting per-atom increments reproduce every PRE→POST charge change, and POST AC charges match the raw MOL2. Maximum per-atom proof residual is 5.55e-17 e. The full 29-bond and 26-atom proof tables are in [the JSON report](am1_bcc_dione_native_fallback_diagnostic.json).

| Quantity | Actual result |
| --- | --- |
| Explicit BCC coverage | 28/29 bonds |
| Native implicit-zero gaps | 1 |
| Nonzero explicit bond increments | 4 |
| Raw Mulliken / PRE / POST / MOL2 total | -0.001000000 e |
| Target formal charge | 0 e |
| Missing-parameter method status | Unsupported coverage; diagnostic only |

A uniform projection would change each atomic charge by approximately 0.00003846 e and bring the total to zero within floating-point precision. **This is hypothetical and was not applied**: repairing rounded total charge cannot supply the missing BCC parameter. No derived dione candidate is promoted.

The parameter snapshot SHA-256 is `05f034901d058ad9331325b6ee3fafa25c75b1ce4f52839fe85032eb92941276`. All 20 source artifact digests are recorded; hashes were checked before and after the read-only audit. Seven other successful campaign cases are unaffected. This report supplies no force-field ranking or final bonded/Lennard-Jones validation.
