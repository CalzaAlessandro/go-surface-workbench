# GO Builder Implementation Plan

> Execute the approved user request continuously in this new folder; independent domains may be delegated. No Git. The separately requested alignment maintenance is recorded with backups in ALIGNMENT.md.

**Goal:** reproducible GO construction and anti epoxide opening, Amber-first preparation, GROMACS export and charge-validation tooling.

**Architecture:** explicit JSON graph; independent geometry, topology, electrostatics and CLI modules. Heavy external programs remain external.

**Tech Stack:** Python 3.11+, NumPy, ParmEd, optional RDKit where needed; uv/pytest/ruff.

**Spec:** docs/DESIGN.md.

## Constraints / review focus

- Periodic graph winding must never be silently converted to a finite Amber/GBSA receptor.
- Infeasible compositions must fail, not silently undershoot.
- New covalent junctions must not receive fabricated validated charges.
- Fragment charge tables are not whole-surface electrostatic validation.
- No inferred bonds during Amber import or export.
- Original current/OLD source files remain unchanged; hashes track copies.

## Tasks

Aggiornamento richiesto dall'utente il 29 ottobre: concentrare il lavoro su
generazione/funzionalizzazione. Rendere esplicita la selezione casuale dei
siti per PEG e frammenti generici: numero e ID facoltativi, frazione calcolata
dall'output, seme riproducibile; registrare quanti innesti sono realizzati.
Anche i gruppi standard sono casuali per default, con count/id per fissarli.
Richieste non raggiunte danno warning strutturato, tentativi alternativi e
arresto senza salvare composizioni parziali. La selezione casuale non modifica la
reazione anti, non elimina altri gruppi e non certifica una distribuzione
all'equilibrio. Testare frazioni derivate, riproducibilità, siti espliciti
e percorso completo da CLI prima della consegna.

- [x] 1. Graph model, finite/periodic construction, holes and functional-group placement. Write tests before implementation; run targeted tests.
- [x] 2. Anti epoxide opening with explicit fragment/transfer-H and chemistry validation. Test stereochemistry and conservation first. Count with optional IDs, derived fractions and incomplete-placement warnings verified.
- [x] 3. Current-fragment assignment, explicit LEaP input and periodic gates; import/export parameterized Amber/GROMACS; GBSA readiness. Test missing terms and round trips.
- [x] 4. Current/OLD/OPLS source audit and Gaussian/RESP/common-reference electrostatic validation tools. Test units, matching and provenance. QM execution and charge-model selection remain scientific work.
- [x] 5. Ordered CLIs, same worker serial/Slurm, examples, portable dependency setup and Italian README. Exercise complete finite and periodic demonstrations.
- [x] 6. Independent scientific/software review, fix substantial findings, full suite and documentation of tested versus unvalidated results. 258 tests pass; see VALIDATION.md.
