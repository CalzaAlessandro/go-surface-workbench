# Evidenza e protocollo QM per i frcmod GO

Data: 5 novembre 2026. Stato: proposta scientifica, nessun fit e nessun nuovo job QM. Sono stati letti solo documenti/fonti esistenti. Questo dossier non sostituisce l'audit numerico indipendente dei file né un confronto su topologie risolte.

## Decisione

Servono verifiche fisiche mirate prima di promuovere un frcmod a produzione; non è ancora dimostrato che serva rifittare tutti i termini. Prima vengono attribuzione esatta dei termini, contesto chimico e riuso dei riferimenti QM accettabili già disponibili. La scelta finale delle cariche resta un prerequisito al fit torsionale e alla validazione congiunta q–LJ. Hessiani e scansioni QM sono riferimenti indipendenti dal modello MM: si possono esaminare ora, ma non giustificano un nuovo frcmod mentre q resta aperto.

## Fonti primarie e limiti di accesso

| ID | Fonte esatta | Verifica e localizzatore | Supporto; limite |
|---|---|---|---|
| E01 | Pinto AV, Ferreira P, Fernandes PA, Magalhães AL, Ramos MJ. *Development of Nanoscale Graphene Oxide Models for the Adsorption of Biological Molecules*. J Phys Chem B 127(2), 557–566 (2023; online 25 ottobre 2022). DOI https://doi.org/10.1021/acs.jpcb.2c06037 ; https://pubmed.ncbi.nlm.nih.gov/36282235/ | Abstract originale e metadati bibliografici letti online il 5 novembre 2026; ACS diretto non accessibile dal tool. | Il lavoro cambia i nonbonded, eredita i bonded GAFF2 e confronta acqua/adsorbimento di amminoacidi; non certifica qualsiasi variante locale GAFF2 né bordi nuovi. |
| E02 | Supporting Information di E01; https://api.figshare.com/v2/articles/21395006 | Copia locale `references/papers/Pinto2023_SI.txt`, righe 12–18, 49–89, 111–132, 142–155; PDF locale preservato e metadata JSON. Endpoint figshare non riaperto dal browser in questa sessione. | Frammenti ACPYPE → libreria GAFF2; PBE/DZVP-MOLOPT-SR-GTH; DDEC6/Chargemol e volumi atomici → LJ; classi di primi vicini. Nessun nuovo fit bonded contro Hessiani/scansioni descritto qui. |
| E03 | Repository autori GAFF2-AIM https://github.com/arvpinto/GO_sheets ; snapshot https://github.com/arvpinto/GO_sheets/tree/1e59f43bd6c4d8c19d0a2a208edf53ce11e00bc5 | README online letto; snapshot locale `references/gaff2_aim/` e `references/gaff2_aim_source_parameters.json`. | Topologie complete e coordinate periodiche; la cartella ff modificata incorpora parametri GO provenienti da ACPYPE. Non trattare il nome amber99sb-ildn della cartella come origine autonoma di un fit bonded. |
| E04 | Wang J, Wolf RM, Caldwell JW, Kollman PA, Case DA. *Development and testing of a general amber force field*. J Comput Chem 25, 1157–1174 (2004). DOI https://doi.org/10.1002/jcc.20035 | Abstract editore aperto online; metadati e ambito generale. | GAFF copre organici e usa anche modelli empirici/euristici; questo non è un benchmark GO denso né una fonte versionata di ogni coefficiente GAFF2. |
| E05 | Cerutti DS, Debiec KT, Case DA, Chong LT. *Links between the charge model and bonded parameter force constants in biomolecular force fields*. J Chem Phys 147, 161730 (2017). DOI https://doi.org/10.1063/1.4985866 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC5669220/ | Articolo completo letto online, Abstract e Introduzione; in particolare discussione di superfici gas/fase acquosa. | Evidenza metodologica di accoppiamento q–bonded e 1–4; sistema proteico, nessuna dimensione dell'effetto trasferibile automaticamente al GO. |
| E06 | Allen AEA, Payne MC, Cole DJ. *Harmonic Force Constants for Molecular Mechanics Force Fields via Hessian Matrix Projection*. J Chem Theory Comput 14(1), 274–281 (2018). DOI https://doi.org/10.1021/acs.jctc.7b00785 ; https://www.repository.cam.ac.uk/items/aa3d9bf5-0a37-4322-9906-929dadd11ab8 | Record istituzionale e abstract autori consultati via ricerca; ACS 403; PDF accepted manuscript indicizzato https://api.repository.cam.ac.uk/server/api/core/bitstreams/06646b1c-4830-47ae-8eba-85e7b8c3549a/content | Modified Seminario deriva bond/angle dal Hessiano QM; non deriva LJ o una superficie di reazione e non rende automaticamente trasferibile un frammento a un foglio. |
| E07 | Garcia NA, Awuah JB, Zhao C, Vuković F, Walsh TR. *Simulation-ready graphene oxide structures with hierarchical complexity: a modular tiling strategy*. 2D Mater 10, 025007 (2023). DOI https://doi.org/10.1088/2053-1583/acb0e1 | Metadati editore Crossmark verificati online: https://crossmark.crossref.org/dialog-content?doi=10.1088%2F2053-1583%2Facb0e1&domain=pdf ; full text IOP bloccato da robots. | Attribuzione bibliografica verificata; nessuna affermazione su lattone dedotta dal titolo/abstract o da citazioni secondarie. |
| E08 | Codice/autori HierGO https://github.com/IFM-molecular-simulation-group/HierGO ; commit 60382a1a1807ef0003ffcc3903392e913b1e5e9a | README online e snapshot locale `references/local_snapshot/opls/PROVENANCE.json`, `Oplsatom2`, `atom-typing.py`, `OPLS-top.f90`. | Inventario reale del sottoinsieme preservato: OH, epossido, etere non epossidico, carbonile, COOH, C–H terminale e C sp²; nessuna tabella lattone o COO⁻ separata. Non prova l'assenza in ogni possibile risorsa upstream. |
| E09 | Implementazione primaria Amber-MD parmchk2 https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/parmchk2.c | Sorgente e definizioni CORR/penalty consultati online. | Un punteggio di somiglianza e una copertura di termini non equivalgono a una misura di accuratezza su GO. |

## Fonti locali lette e affermazioni legate

- `docs/CHARGES_BONDED_LJ_COUPLING.md`: distinzione coppia q–LJ, bonded ereditati GAFF2-AIM, validazione finale non eseguita.
- `docs/BONDED_QM_PROTOCOL.md`: riferimento finito ωB97X-D/6-311++G(d,p), Opt/Freq/Stability separati, QA da completare. Non suggerire una nuova campagna indistinta.
- `docs/NEW_FRCMOD_CANDIDATES.md`: candidati unfitted e collisione `cx-cx-cx` fuori da triangolo C–C–C; il contesto ha precedenza sul semplice match di tipi.
- `docs/FRCMOD_AUDIT.md:9,15–21,25–39,55–57`: provenienza parmchk + analogie manuali riferita dall'utente; override ≠ FF completo; zeri contestuali, continuazioni e valori non convalidati fisicamente.
- `docs/FORCE_FIELD_CANDIDATES.md`: convenzioni complete e scope dei preset; GAFF2-AIM basale OH/EPX, non un modello già parametrizzato dei bordi COOH/COO⁻.
- `parameters/comparison_sources/HOSHIN_OLD_CURRENT_20261105/README.md` e `MANIFEST.json`: fonti byte-identiche e origine Hoshin/OLD/current, senza prova di un nuovo fit QM.
- `references/go_prior_art_20261104/AmberGO_Tutorial.txt:743–772`: scelta GAFF ridotta per gestione/copertura pratica; l'aggettivo “validated” del tutorial non sostituisce un confronto fisico indipendente.
- `reports/dft_bonded_reference_QA_20261030.json`: rapporto datato 30 ottobre, 7/8 preliminari pass e `final_parameter_validation: false`. In `pyrene_epoxide` distanza del legame atteso 15–16 = 2.110471977 Å e topology screen fallito, pur con frequenze positive e stabilità elettronica. Non sono stati recuperati/verificati qui i raw Gaussian/fchk; nessuna conclusione sulla reazione fisica o stato attuale dei job.

## Inventario HierGO direttamente verificato

`Oplsatom2` contiene 17 righe, non 17 gruppi: HYD 3, EPX 2, ETH 2, CAB 2, CARB 5, terminale GRAT 2 e GRA 1. Sono ruoli di atomo; i nomi non definiscono automaticamente tutti i grafi, la protonazione o la topologia finale.

`atom-typing.py:204–238` usa i vicini geometrici: O a due C, con i due C reciprocamente vicini → OE/epossido; altrimenti OT/etere. Non presenta una classe separata per etere basale vs bordo né per lattone. Un lattone contiene C(=O)–O–C chiuso in anello e richiede contesto estereo/estere specifico: la sola presenza di OT + OC non ne prova l'assegnazione corretta. CARB include H acido e non equivale a lattone o COO⁻.

Discrepanza da mantenere visibile: `Oplsatom2:6` riporta ETH O −0.1700 e, mentre i commenti di `atom-typing.py:174,344` riportano −0.2300 e. `OPLS-top.f90:40–68` legge Oplsatom2. La tabella operativa preservata è dunque la fonte da seguire per il confronto, senza correggere silenziosamente i commenti. Non importare queste cariche in GAFF2 come un preset HierGO completo.

## Riferimenti aggiuntivi: scelta metodologica, non conclusioni bibliografiche

Tutte le quantità di scansioni, distanze e soglie riportate in REPORT.html sono proposte di disegno sperimentale, non risultati pubblicati o nuovi calcoli. I criteri numerici finali vanno fissati rispetto all'osservabile desiderata e alla convergenza QM; non sono stati inventati errori o tempi CPU.

Questa proposta riusa il metodo QM già documentato per coerenza di confronto, senza dichiararlo livello definitivo per reazioni, fogli periodici o anioni in soluzione. Per energie di interazione si propone dispersione, controllo BSSE/convergenza e geometrie dichiarate; la decisione effettiva sul livello e solvente va presa nel protocollo separato prima dei calcoli.

Memoria usata solo per routing: `MEMORY.md:47–60`, sessione 01a07be4-18ac-7111-bcc4-ab20caec3b9d; le affermazioni operative sopra sono state rilette dalle fonti correnti. Nessuna modifica alle memorie.
