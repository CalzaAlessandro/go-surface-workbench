# GO Builder (`go-surface-workbench`)

[![Status: Under Active Development](https://img.shields.io/badge/status-active%20development-orange.svg)]()
[![Python: >=3.11](https://img.shields.io/badge/python-3.11+-blue.svg)]()

> [!IMPORTANT]
> **Work in Progress / Active Research Repository**
> This repository contains the source code, parameterization protocols, and validation toolchains for the GO Builder workbench.
> The project is currently under active scientific development. Force-field parameters, atomic charges, and model topologies are preliminary and subject to revision prior to peer-reviewed publication.
> 
> Molecular engine outputs (Gaussian/ORCA/Amber logs and trajectories) are tracked separately. External dependencies and molecular dynamics engines must be configured as described in [Engine setup](docs/ENGINE_SETUP.md).

Nome provvisorio della repository: **go-surface-workbench**.
Vedere [direzione scientifica](docs/ROADMAP.md) e [prestazioni misurate](docs/PERFORMANCE.md).

Lo studio dei metodi di carica e dei tre candidati di force field è in
[research/parameterization](research/parameterization/README.md); i nuovi
parametri restano separati dalle librerie usate nelle campagne esistenti.

Il [catalogo dei preset di carica](docs/CHARGE_PRESETS.md) include ora il
riferimento pubblicato GAFF2–AIM con cariche e Lennard-Jones originali, oltre
alle librerie correnti e al riferimento HierGO/OPLS. La nuova derivazione
resta RESP; l'assegnazione generale di tutti i preset è ancora da completare.

Generatore generale di grafene ossido con connettività esplicita. Questa
cartella contiene codice, esempi, test e copie dei parametri con provenienza;
non include Amber o Gaussian. I modelli GOX/REGOX corrected e i parametri OLD
restano riferimenti: non vengono ricostruiti o sostituiti automaticamente.
L'unica manutenzione applicata al workflow precedente riguarda gli input
di allineamento richiesti dall'utente, con backup; vedere
[ALIGNMENT.md](docs/ALIGNMENT.md).

Il percorso principale per **creare e funzionalizzare il GO** è
`01_build_GO.py` → `02_open_epoxides.py` (opzionale) → `03_prepare_amber.py`.
`prepare_fragment.py` prepara una volta il reagente da innestare; la stessa
ricetta si può riusare per diversi caricamenti e semi casuali.

## Installazione e trasferimento

Richiede Python 3.11 o successivo. Dalla cartella:

```bash
uv sync
uv run python go.py --help
uv run pytest -q
```

`pyproject.toml` e `uv.lock` definiscono NumPy, ParmEd e RDKit. Per un sistema
senza uv si può usare l'esportazione generata `requirements.txt` in un ambiente
Python separato. Amber/AmberTools (`tleap`, `sander`, `MMPBSA.py`, `espgen`,
`resp`) e Gaussian (`g16`) restano installazioni di sistema.

Questa snapshot contiene sorgenti, esempi, test, risorse di riferimento,
documentazione e file delle dipendenze. Non include `software/`, gli eseguibili
di Vina/Open Babel, MGLTools o i relativi archivi, né gli output di `validation`.
Ricreare `.venv` sul sistema di destinazione: un ambiente binario macOS non si
trasferisce su Linux. Installare separatamente gli engine e gli strumenti di
sistema richiesti; vedere [Engine setup](docs/ENGINE_SETUP.md). Per importare le
API direttamente da questa cartella usare
`PYTHONPATH=src uv run python ...`; gli script della radice gestiscono già il
percorso dei sorgenti.

## Script e ordine d'uso

| Script | Uso |
|---|---|
| `prepare_fragment.py` | Opzionale: prepara un frammento da SDF/MOL/MOL2/SMILES e sceglie esplicitamente il sito reattivo |
| `01_build_GO.py` | Crea foglio, difetti, gruppi e terminazioni H |
| `02_open_epoxides.py` | Opzionale: apertura anti degli epossidi con il frammento scelto |
| `03_prepare_amber.py` | Assegna i template GAFF2 e prepara LEaP; `--run` esegue AmberTools esterno |
| `04_prepare_docking_receptor.py` | Estrae un GO nativo o equilibrato e crea box 3D senza appiattirlo |
| `05_dock_boxes.py` | Docking Vina oppure GNINA con scoring Vina e CNN disattivata |
| `06_assemble_docked_complex.py` | Assembla il recettore completo e ligandi già parametrizzati nelle coordinate delle pose |
| `07_export_gromacs.py` | Esporta una topologia Amber completa in `.top`, `.itp`, `.gro` |
| `audit_charges.py` | Confronto descrittivo di sorgenti CURRENT, OLD e HierGO/OPLS |
| `prepare_QM.py` | Prepara un riferimento Gaussian/MK su un frammento finito esplicitamente saturato |
| `compare_ESP.py` | Confronta cariche sulla stessa geometria e griglia quantistica |
| `import_Amber.py` | Ingresso alternativo a 01: importa una topologia esistente senza ricostruire legami |
| `check_model.py` | Controlla grafo, prerequisiti GBSA e maschere disponibili |
| `update_alignment.py` | Inventario e migrazione con backup degli input CPPTRAJ; non ricalcola traiettorie |

Gli stessi comandi sono disponibili da `go.py`: `fragment`, `build`,
`functionalize`, `amber`, `gromacs`, `audit`, `qm`, `esp-import`, `esp`, `import`,
`check`, `snapshot`, `dock`, `complex`, `alignment`. `--help` elenca le opzioni di ogni script.

## Generare un GO

```bash
uv run python 01_build_GO.py --recipe examples/finite_GO.json \
  --output work/GO.json --pdb work/GO.pdb
uv run python 03_prepare_amber.py --model work/GO.json \
  --output-dir work/amber --radii mbondi2 --run
uv run python 07_export_gromacs.py --topology work/amber/go.prmtop \
  --coordinates work/amber/go.inpcrd --output-prefix work/gromacs/GO
```

Senza `--run`, 03 prepara gli input ma non certifica una parametrizzazione
completa. Con `--run` richiede `tleap` nel PATH e verifica il risultato.
I risultati esistenti non vengono sovrascritti: scegliere un nuovo percorso.

La ricetta JSON contiene:

| Campo | Significato |
|---|---|
| `nx`, `ny` | Numero di celle rettangolari, ciascuna inizialmente con quattro C |
| `bond_length` | Distanza C–C iniziale in Å, default 1,42 |
| `pbc` | `[false,false,false]` per foglio finito; `[true,true,false]` per il grafo periodico |
| `vacuum` | Dimensione z della cella di costruzione in Å; non aggiunge solvente |
| `vacancies` | Numero oppure lista di ID dei carboni da rimuovere |
| `holes` | Lista di `{center:[x,y], radius:r}` in Å oppure `{size:n, seed_atom:id}` |
| `groups` | Conteggi di `OH`, `epoxide`, `carbonyl`, `COOH`, `COO`; ogni gruppo può specificare anche `count` e `id` |
| `faces` | `both`, `up`, `down` per i gruppi basali |
| `seed` | Seme casuale riproducibile |

Le dimensioni della cella sono 3 × bond_length × nx e √3 × bond_length × ny.
La rimozione dei carboni pendenti modifica il conteggio finale. OH ed epossidi
richiedono siti basali; carbonili e carbossili richiedono bordi o difetti.
Il programma riporta i conteggi effettivi; una richiesta non realizzabile o una
ricerca di posizionamento fallita termina con errore, senza ridurre i conteggi.

XPS può guidare la scelta dei conteggi, ma non è obbligatorio. Questa versione
riceve conteggi, non interpreta automaticamente uno spettro o una percentuale.
I conteggi sono gruppi chimici: un epossido coinvolge due carboni.

### Come avviene la costruzione

Il reticolo nasce da celle esagonali analitiche di quattro C, con C–C iniziale
di 1,42 Å. I legami vengono registrati durante la generazione. Vacanze e buchi
rimuovono atomi dal grafo; i carboni pendenti vengono eliminati e il foglio deve
restare connesso. Si scelgono poi i siti: OH ed epossidi sul piano basale,
carbonile/COOH/COO sui bordi esterni o dei difetti. Un epossido occupa una coppia
C–C; gli altri gruppi usano un sito, senza riutilizzi.

**La scelta è casuale per tutti i gruppi standard, salvo ID espliciti.**
Ogni tentativo sceglie l'intera disposizione, poi aggiunge e controlla nell'ordine
COOH, COO, carbonili, epossidi e OH. Se un posizionamento fallisce, riparte da
una nuova disposizione, fino a 128 tentativi. Non enumera tutte le permutazioni.
`seed` riproduce la ricerca; `faces: both` sceglie casualmente la faccia dei
gruppi basali, senza imporre esattamente metà sopra/metà sotto. I bordi liberi
sono infine saturati con H. Sono coordinate iniziali, da rilassare dopo aver
assegnato i parametri.

Per fissare siti specifici si può sostituire un conteggio con un oggetto:

```json
"groups": {
  "OH": {"count": 1, "id": [58]},
  "epoxide": 2,
  "COOH": 1
}
```

Qui OH usa il carbonio 58, mentre gli altri gruppi restano casuali. Per gli
epossidi `id` è una lista di coppie di carboni adiacenti, per esempio
`{"count": 1, "id": [[58, 59]]}` in una ricetta distinta. Sono gli ID stabili
dei carboni nel JSON, prima della decorazione, non necessariamente i seriali
del PDB. Devono esistere dopo buchi/pulizia ed essere siti chimicamente adatti.
I siti fissati sono riservati prima della scelta casuale degli altri gruppi.
Se manca `count` viene ricavato dagli ID; se entrambi sono presenti devono
coincidere. Gli ID fissati non vengono sostituiti quando sono ingombrati.
Un esempio completo è `examples/GO_mixed_sites.json`.

## Amminoacidi, PEG, peptidi e ciclodestrine

La funzionalizzazione può essere **casuale e riproducibile**, come il
posizionamento dei gruppi standard. Il sito reattivo del frammento resta
esplicito; i siti sulla superficie possono essere scelti automaticamente.

```bash
uv run python prepare_fragment.py --smiles 'NCC(=O)O' --inspect --json
uv run python prepare_fragment.py --smiles 'NCC(=O)O' \
  --attachment-atom 1 --seed 17 --output work/open_glycine.json
uv run python 02_open_epoxides.py --model work/GO.json \
  --recipe work/open_glycine.json --output work/GO_glycine.json \
  --pdb work/GO_glycine.pdb
```

Per un peptide, un PEG amminato o una ciclodestrina ammino-funzionalizzata usare
`--input molecola.sdf` e scegliere l'atomo dalla mappa `--inspect --json`.
L'identificativo è un indice **da 1**, conservato prima dell'aggiunta di H.
Il sito non viene scelto dal nome della molecola: per la lisina, ad esempio,
α-N ed ε-N sono scelte distinte. MOL2 richiede tutti gli H espliciti.

L'apertura implementata consuma un legame C–O dell'epossido, forma C–nucleofilo
e trasferisce un H del donatore all'ossigeno rimasto. L'attacco iniziale è dal
lato opposto al legame C–O rotto; l'OH rimane sulla faccia opposta al sostituente.
Atomi, carica formale, origine degli atomi e stereochimica del frammento vengono
tracciati. La ricerca considera siti e torsioni evitando sovrapposizioni gravi.

L'esempio seguente genera un GO e innesta casualmente due copie del PEG
amminato incluso, usando metà dei quattro epossidi:

```bash
uv run python 01_build_GO.py --recipe examples/GO_for_random_grafting.json \
  --output work/GO_for_PEG.json --pdb work/GO_for_PEG.pdb
uv run python 02_open_epoxides.py --model work/GO_for_PEG.json \
  --recipe examples/open_epoxide_PEG2_amine.json \
  --count 2 --seed 17 \
  --output work/GO_PEG.json --pdb work/GO_PEG.pdb
```

Opzioni valide anche per amminoacidi e altri frammenti:

| Opzione | Selezione degli innesti |
|---|---|
| `--count 2` | Due innesti; i siti disponibili sono scelti casualmente per default |
| `--id 2 3` | Opzionale: usa quei precisi anelli epossidici, identificati in `metadata.groups` |
| `--seed 17` | Riproduce scelta dei siti e ricerca dell'orientamento, a parità di input/versione |

La frazione è **calcolata dall'output**, come innesti ottenuti diviso epossidi
intatti all'inizio del comando: due su quattro danno 0,5 (50%); due su tre
danno circa 0,667 (66,7%). Non è una flag di input. `--id` è facoltativo e può
accompagnare `--count`: in tal caso il numero di ID deve coincidere. Con il
solo `--id` il numero viene ricavato dalla lista. La CLI sostituisce la scelta
della ricetta; senza opzioni usa il `count` della ricetta, oppure uno come
default. Numero richiesto/ottenuto, frazione effettiva, siti e seme sono salvati
in `metadata.functionalization_history` e riportati a schermo.
Una nuova esecuzione con altro frammento usa gli epossidi rimasti, permettendo
superfici con modifiche miste. In questa fase `--id` indica ID di **gruppo
epossidico**, mentre gli ID nella ricetta di costruzione indicano carboni.

Il campionamento usa siti mescolati e rifiuto degli ingombri: non promette
equiprobabilità fra siti accettati o una distribuzione termodinamica. Ruota
il frammento attorno al nuovo legame; la sua conformazione interna rimane
quella preparata o fornita. La faccia dell'innesto è opposta all'ossigeno
dell'epossido scelto. Per il foglio nativo `faces` nella ricetta di costruzione
controlla la distribuzione iniziale degli ossigeni sulle due facce.

### Superficie insufficiente o ingombri irrisolvibili

Per gruppi standard e frammenti, una richiesta non raggiunta produce
`[WARNING]` e stato `INCOMPLETE`, con conteggi richiesti, miglior disposizione
provata quando disponibile, capacità dei siti e numero di tentativi. La capacità
residua dei gruppi standard considera già i carboni riservati da ID espliciti.
I suggerimenti indicano se aumentare il foglio, ridurre il caricamento, cambiare
seme/siti o provare un altro conformero del frammento.

Le quantità superiori alla capacità nota vengono segnalate prima della ricerca.
Il fallback geometrico prova altre disposizioni, siti e orientamenti nei limiti
descritti. A esaurimento il comando termina con codice 1 e non salva un modello
parziale. Il miglior tentativo non dimostra il massimo caricamento fisicamente
possibile. Con `--json` la stessa diagnostica è un unico oggetto strutturato.

Il prodotto è una geometria iniziale, non una traiettoria di reazione o una
struttura equilibrata. Sono implementati donatori neutri amminici, alcolici e
tiolici con un H trasferibile e valenze supportate. Ammidi e ammine già
protonate non vengono trattate automaticamente come ammine nucleofile.
La protonazione del prodotto alle condizioni della MD va definita e
parametrizzata esplicitamente.

**La nuova giunzione deve essere parametrizzata.** 03 rifiuta il prodotto finché
non riceve `--parameter-templates` con una mappa completa e tracciabile di
atomi/tipi/cariche, legami e un `frcmod` dedicato. Non riutilizza automaticamente
le cariche dell'epossido intatto sul prodotto. Dettagli e schema in
[FRAGMENTS.md](docs/FRAGMENTS.md) e [CHARGES.md](docs/CHARGES.md).

I due articoli indicati dall'utente sono discussi in
[FUNCTIONALIZATION_LITERATURE.md](docs/FUNCTIONALIZATION_LITERATURE.md).
La modalità anti è la scelta costruttiva richiesta dall'utente; non viene
presentata come stereochimica misurata in quei lavori.

## Amber, periodicità, docking e MM/GBSA

Il percorso Amber supportato è **foglio covalente finito**, anche in una scatola
di solvente periodica. I legami che attraversano le immagini di un reticolo
covalente infinito vengono conservati nel JSON, ma l'esportazione di quel
reticolo verso Amber o attraverso Amber verso GROMACS è bloccata finché non
esiste una verifica specifica del motore. Il PDB non conserva gli spostamenti
di immagine. `examples/periodic_GO.json` dimostra il grafo, non una MD Amber.

Il recettore per il docking/GBSA è GO **più tutte le modifiche covalenti**.
Mantenere la topologia Amber come riferimento anche se si prepara un PDBQT:
cariche, ordine degli atomi e legami usati nel docking non devono sostituire
automaticamente quelli del modello MD.

### Recettore nativo oppure equilibrato e box tridimensionali

```bash
uv run python 04_prepare_docking_receptor.py \
  --topology compl_solv.parm7 --coordinates md_10ns.rst \
  --receptor-mask ':1-494' --box-size 20 20 20 --overlap 5 --padding 5 \
  --output-dir work/GO_snapshot
```

Per REGOX usare `--receptor-mask ':1-628'`. Una topologia costruita da zero può
avere numerazione diversa: leggere le maschere prodotte da 03, senza riutilizzare
questi conteggi. Si può passare anche un NetCDF con `--frame N` (indice da 1).
Il recettore deve contenere interamente le molecole covalenti selezionate;
una funzionalizzazione covalente non va esclusa dalla maschera.

04 riunisce le immagini del foglio finito seguendo i legami e conserva la sua
forma tridimensionale. Scrive `receptor.parm7`, `receptor.rst7`, `receptor.pdb`,
il grafo e `boxes.json`. I box coprono il volume della struttura, comprese
ondulazioni e ripiegamenti; non assumono due piani z. Anche cavità non
accessibili possono rientrare nei box: l'accessibilità non viene certificata.

### Vina e GNINA empirico

Con GNINA si possono usare direttamente PDB del recettore e SDF del ligando:

```bash
uv run python 05_dock_boxes.py --boxes work/GO_snapshot/boxes.json \
  --receptor work/GO_snapshot/receptor.pdb --ligand ligand.sdf \
  --engine gnina --gnina /percorso/gnina --cpus 8 \
  --runs 3 --exhaustiveness 16 --num-modes 10 \
  --output-dir work/dock_gnina --run
```

Il worker imposta `--scoring vina --cnn_scoring none --no_gpu
--pose_sort_order Energy`. Con Vina usare `--engine vina`, recettore e ligando
PDBQT già preparati e un eseguibile Vina installato (o indicato con `--vina /percorso/vina`).
Entrambi ricevono `--cpu` esplicito; `--exhaustiveness` controlla la ricerca e
non sostituisce la richiesta di risorse Slurm. I run e i box sono sequenziali.
Senza `--run` viene scritto soltanto il piano; `--slurm --run` invia il worker
che eseguirà il docking. Log, semi, hash e pose migliori sono registrati.
Dettagli in [DOCKING_BOXES.md](docs/DOCKING_BOXES.md).

Il vecchio comando `sed` ricordato dall'utente sostituiva `0.000 OA` con
`-0.342 OA`: cambiava la carica, non il tipo OA. Non è una correzione generale
dei tipi PDBQT e non viene applicato. Lo scoring Vina standard ignora le
cariche parziali fornite; le cariche Amber restano quelle della topologia.
Vedere la [FAQ ufficiale di Vina](https://autodock-vina.readthedocs.io/en/latest/faq.html).

### Assemblaggio per MD e MM/GBSA

```bash
uv run python 06_assemble_docked_complex.py \
  --receptor-topology work/amber/go.prmtop \
  --receptor-coordinates work/amber/go.inpcrd \
  --ligand-topology ligand.parm7 --ligand-coordinates docked_pose.rst7 \
  --output-dir work/complex_pose1
```

Le coordinate devono essere complete e già nello stesso riferimento della
posa. Non si ricostruiscono gli H da PDBQT. Per più ligandi simultanei ripetere
le due opzioni ligando; pose alternative sovrapposte richiedono complessi
separati. L'assemblaggio conserva i legami interni e non crea legami fra
recettore e contaminante. Le maschere sono calcolate dalla topologia risultante.
Dettagli in [DOCKING.md](docs/DOCKING.md).

03 permette `--radii bondi|mbondi|mbondi2|mbondi3`; scegliere raggi coerenti con
il modello GB e mantenerli uguali fra recettore, ligando e complesso. La
compatibilità tecnica GBSA non dimostra l'accuratezza del modello del GO.
Il legame covalente GO–funzionalizzazione resta dentro il recettore e non va
tagliato per un normale calcolo di binding del contaminante.

`import_Amber.py` importa topologie esistenti con coordinate esplicite. Conserva la cella
reale, oppure `null` se assente; non inventa legami o cariche formali. Riconosce
come siti di apertura soltanto epossidi basali `EPX` con il grafo locale atteso.
Altri nomi/motivi richiedono un'annotazione esplicita. Amber non conserva gli
ordini chimici di legame: l'importazione numerica non li ricostruisce.

## Cariche: confronto e validazione

```bash
uv run python audit_charges.py --output-dir work/charge_audit
uv run python 01_build_GO.py --recipe examples/QM_sheet.json --output work/QM.json
uv run python prepare_QM.py --model work/QM.json --output-dir work/QM_reference
```

`prepare_QM.py` prepara un single point HF/6-31G* con ESP MK e un fit RESP esplorativo.
Gaussian e RESP non partono automaticamente. Il protocollo proposto non viene
attribuito retroattivamente alle cariche CURRENT, di cui l'utente riferisce una
parametrizzazione Amber dei singoli frammenti.

Dopo Gaussian, controllare convergenza e terminazione; lo script generato
`prepare_resp.sh` usa AmberTools esterno. Importare la griglia e confrontare:

```bash
uv run python go.py esp-import --esp-file reference.esp --model work/QM.json \
  --total-charge 0 --output work/reference.json
uv run python compare_ESP.py --reference work/reference.json \
  --candidates candidates.json --output-dir work/ESP_comparison
```

`candidates.json` associa nomi a modelli con stessi ID, elementi e coordinate,
più le cariche da verificare. Geometria, protonazione e carica totale devono
coincidere. Non vengono inventate corrispondenze o rinormalizzate cariche.
Servono frammenti e conformazioni comuni prima di scegliere fra CURRENT, OLD e
OPLS; un punteggio ESP isolato non seleziona automaticamente il force field
migliore per l'adsorbimento. Vedere [CHARGES.md](docs/CHARGES.md).

## Esecuzione locale e Slurm

La creazione del grafene con `01_build_GO.py` o `go.py build` è locale:
accetta ricette personalizzate o preset e non espone opzioni Slurm. `--dry-run`
mostra il comando locale senza creare il modello o file di invio.

Gli altri stadi mantengono lo stesso worker per esecuzione locale e Slurm.
Per esempio, per preparare gli input Amber sul cluster:

```bash
uv run python 03_prepare_amber.py --model work/GO.json \
  --output-dir work/amber_slurm --slurm --cpus 1 --memory 2G \
  --time 00:10:00 --partition ciamician_cpu
```

Opzioni comuni agli stadi successivi alla creazione: `--slurm`, `--dry-run`, `--submission-dir`, `--cpus`, `--memory`,
`--time`, `--partition`, `--account`, `--setup`, `--python`, `--json`.
`--setup /percorso/site_setup.sh` inizializza i moduli del proprio sistema;
`--python` indica l'interprete disponibile sul worker. Non sono scelti moduli o
partizioni automaticamente. `--slurm --dry-run` scrive un worker verificabile
senza inviarlo; la ricevuta registra comando, hash e ID realmente restituito.
Gli input relativi sono risolti nella directory da cui si lancia il comando.

## Stato e differenze dagli originali

[VALIDATION.md](VALIDATION.md) distingue prove completate e scienza ancora da
validare. [COMPARISON.md](docs/COMPARISON.md) confronta il metodo con HierGO e
GOPY. I file originali copiati sono registrati in
`resources/reference/SOURCES.json`; HierGO/OPLS mantiene commit, hash e licenza
in `resources/reference/opls/PROVENANCE.json`.
