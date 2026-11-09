# Nuovi gruppi GO: riferimenti di ricerca non parametrizzati

Sono stati preparati sette riferimenti C/H/O neutri, senza radicali nel grafo e con elettroni pari. Il singoletto chiuso è un’ipotesi iniziale: nessun QM, verifica dello stato fondamentale, fit RESP o test MD è stato eseguito.

| Riferimento | Motivo | Formula | Fonte dello SMILES |
|---|---|---|---|
| `dibenzofuran_ether_edge` | non_epoxide_ether, anello 5 | C12H8O | [dibenzofuran](https://www.ebi.ac.uk/chebi/CHEBI:28145) |
| `diphenyl_ether_bridge_control` | non_epoxide_ether, anello non richiesto | C12H10O | [diphenyl ether](https://www.ebi.ac.uk/chebi/CHEBI:39258) |
| `phthalide_lactone5_edge` | lactone, anello 5 | C8H6O2 | [phthalide](https://www.ebi.ac.uk/chebi/CHEBI:38085) |
| `coumarin_lactone6_edge` | lactone, anello 6 | C9H6O2 | [coumarin](https://www.ebi.ac.uk/chebi/CHEBI:28794) |
| `oxolan2ol_lactol_control` | lactol, anello 5 | C4H8O2 | [oxolan-2-ol](https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=200907060640520321) |
| `phenol_aromatic_OH_control` | phenolic_OH, anello non richiesto | C6H6O | [phenol](https://www.ebi.ac.uk/chebi/CHEBI:15882) |
| `naphthoquinone14_edge` | quinone, anello 6 | C10H6O2 | [1,4-naphthoquinone](https://www.ebi.ac.uk/chebi/CHEBI:27418) |

Il dibenzofurano è un piccolo riferimento aromatico fuso per un ambiente etereo di bordo/difetto; il difeniletere è un controllo di ponte diarile. Ftalide, cumarina e 1,4-naftochinone sono riferimenti aromatici fusi; fenolo è un controllo aromatico monociclico. Oxolan-2-ol è soltanto un controllo **alifatico** della connettività lattolica. Nessuno è un sostituto equivalente di un foglio GO completo o una prova della presenza/abbondanza dei gruppi nei nostri campioni XPS.

Le fonti ChEBI sono state aperte il 5 novembre 2026. Per oxolan-2-ol lo SMILES e InChIKey sono stati recuperati dalla pagina J-GLOBAL/JST indicizzata; l’accesso diretto ha restituito HTTP 503. Questo limite di accesso è conservato nel manifest delle fonti. Le identità grafiche dei sette input sono state controllate anche mediante formula e InChIKey RDKit, prima di assegnare un singolo stereoisomero geometrico dove lo SMILES non specifica la stereochimica.

Il gate separa etere non epossidico, estere ciclico a cinque/sei termini, emiacetale ciclico senza carbonile, fenolo neutro e chinone para coniugato. Nei modelli ogni atomo ha ID stabile, elemento, carica formale intera, ruolo `extended_…`, equivalenza del grafo e coordinate. **Non ci sono cariche parziali.** I nuovi ruoli sono `UNPARAMETERIZED_RESEARCH_REFERENCE`; `modular.partition_model` e il builder non sono dichiarati compatibili.

Questa seconda esportazione ripristina la rappresentazione aromatica RDKit dopo MMFF94s, prima di assegnare ruoli e scrivere JSON/SDF. Il test di round trip della cumarina verifica che non cambino gli ordini di legame fra i due formati. La prima esportazione rimane storica e non va usata come pannello corrente.

I sette `model.json`/`initial.sdf`/`initial.xyz` contengono geometrie ETKDGv3 + MMFF94s convergenti con seme 20261105. Sono esclusivamente geometrie iniziali. `geometry_opt/reference.gjf` prepara Opt/Freq/MK e Stable a wB97XD/6-311++G(d,p). `HF_RESP_template/reference.gjf` prepara un **template** HF/6-31G* MK/Stable con coordinate iniziali: va rigenerato usando una geometria DFT accettata prima di qualunque uso scientifico. Nessuna orientazione, sottomissione o carica è stata generata.

HierGO comprende ETH per l’etere non epossidico, ma non un set esplicito di lattoni o COO−. I lattoni e il lattolo qui sono estensioni di ricerca ispirate alla letteratura, non blocchi già stabiliti di HierGO. Le motivazioni GO e i limiti di attribuzione sono nella revisione locale `docs/MISSING_GROUPS_LITERATURE.md`; l’identità dei piccoli riferimenti non prova l’identità dei siti nel materiale.

Per arrivare a parametri trasferibili occorrono QM accettato, conservazione della connettività e dell’ordine atomico, fit e validazione su contesti esclusi dalla calibrazione, definizione di blocchi/anchor senza sovrapposizioni e bonded/LJ/scaling coerenti. Non si impone automaticamente la neutralità di sottoinsiemi a cariche RESP molecolari future e non si trasferiscono parametri EPX agli eteri o COOH ai lattoni.

La trasferibilità da circa **4 nm a circa 10 nm**, necessaria per proteine/DNA, resta un requisito da verificare. Il metodo finale dovrà conservare le stesse classi di grafo e gli stessi vincoli di carica molecolare/gruppo tra dimensioni. Servono verifiche tenute fuori dalla calibrazione su fogli più grandi, variando bordo/interno, ossidazione bassa e alta a dimensione fissa, dimensioni diverse a composizione comparabile, gruppi vicini e protonazione, con controlli elettrostatici/MD pertinenti. **Non si richiede nuovo QM dell’intero foglio a ogni dimensione produttiva.** Un nuovo ambiente non supportato richiede riferimenti mirati e una nuova verifica della copertura.

I test reali del grafo sono stati scritti e visti fallire prima dell’implementazione. Verifica completa: 265 test passati, di cui 20 nell’estensione. Le cariche, gli input scientifici precedenti e il codice di assegnamento esistente restano fuori dall’estensione.

Per rigenerare il pannello, dalla directory PARAMETERIZATION, usare una nuova directory:

```sh
PYTHONPATH=src uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python -m go_param.extended_fragments --output-dir reports/EXTENDED_GO_GROUPS_20261105_attempt_003 --cpus 8 --mem-gb 6
```

Per ricostruire solo questo README, la tavola 2D e il rapporto:

```sh
uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python reports/EXTENDED_GO_GROUPS_20261105_attempt_002/build_report.py
```

Vedere anche `catalog.json`, `source_manifest.json`, `artifact_manifest.json` e `docs/EXTENDED_GO_GROUPS_20261105.md`.
