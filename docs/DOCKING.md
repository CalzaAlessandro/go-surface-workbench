# Complessi dopo docking: Amber e MM/GBSA

Il recettore è l'intera superficie parametrizzata: GO, gruppi ossigenati e
qualsiasi linker o funzionalizzazione covalente. Il contaminante adsorbito è
un componente distinto. La separazione non dipende da nomi come GOX o da
maschere fissate a priori.

## Input e comando

Servono topologie Amber complete di recettore e ligando e coordinate delle
pose nello stesso riferimento cartesiano. Tutti i formati richiedono ogni atomo della topologia, inclusi gli H.
Il comando conserva cariche, tipi, masse e parametri Amber; le eventuali
cariche contenute nel file di docking non vengono usate.

| Coordinate | Corrispondenza con il template Amber |
|---|---|
| Restart Amber (`.rst7`, `.inpcrd`, `.ncrst`, `.rst`, `.restrt`, `.crd`) | Ordine completo della topologia, fornito esplicitamente dall'utente |
| PDB/ENT | Ordine completo; nomi, elementi espliciti e confini dei residui verificati |
| PDBQT | Primo MODEL; nomi atomici univoci, un solo residuo; ordine dell'albero Vina permutabile |
| SDF/MOL | Prima posa; grafo esplicito ed elementi corrispondenti al template, tutti gli H presenti |

```bash
uv run python 06_assemble_docked_complex.py \
  --receptor-topology GO/go.prmtop --receptor-coordinates GO/go.inpcrd \
  --ligand-topology contaminante.parm7 --ligand-coordinates posa.rst7 \
  --output-dir complesso_posa1
```

Oppure aggiungere `--slurm` e le opzioni di partizione/risorse: la logica
scientifica è identica all'esecuzione locale.

Un PDBQT è accettato solo se contiene **tutti** gli atomi del template.
I nomi devono essere univoci in una singola molecola/residuo; i tipi AutoDock
verificano l'elemento. L'ordine di ROOT/BRANCH non determina l'ordine Amber.
Pseudoatomi di chiusura e dummy di idratazione non sono supportati.
[Definizioni ufficiali dei tipi Vina](https://github.com/ccsb-scripps/AutoDock-Vina/blob/develop/src/lib/atom_constants.h).

Gli SDF/MOL vengono letti con RDKit senza aggiungere H. La connettività
completa, inclusi gli H legati a ciascun atomo, deve coincidere. Gli ordini di
legame non vengono ricavati da Amber: il confronto usa adiacenza ed elementi,
mentre il modello elettronico/force field resta quello del template Amber.
Non è una verifica indipendente di aromaticità, stati elettronici o stereo E/Z.

Il mapping automatico enumera solo le corrispondenze degli atomi pesanti
(max 4096); gli H vengono associati per atomo genitore. Una simmetria che
scambia atomi con cariche RESP, tipi o masse diversi viene rifiutata. Gli H
dello stesso genitore devono essere equivalenti nei parametri atomici. Per
risolvere una corrispondenza nota, inserire nello **SDF** la proprietà
`GO_BUILDER_ATOM_NAMES`: un array JSON di tutti i nomi Amber, nell'ordine
degli atomi SDF. Tutti i nomi del template devono essere univoci. Esempio:

```text
> <GO_BUILDER_ATOM_NAMES>
["H2", "C1", "H1", "H3", "H4"]

```

Anche con la mappa esplicita, elementi e grafo completo vengono verificati.
Nessun indice SDF viene assunto automaticamente come indice Amber.

Per gli SDF con centri tetraedrici C/Si potenzialmente stereogenici servono
coordinate di riferimento complete del template, precedenti al docking.
La geometria della posa viene confrontata con quella del riferimento senza
allineamenti: inversioni e centri quasi planari vengono rifiutati. Fornire
`--ligand-reference-coordinates contaminante.rst7` per ciascun ligando.
Nell'API Python l'opzione è `ligand_reference_coordinates=[...]`, con un
file o `None` per ogni ligando. Il riferimento deve essere un restart Amber
o un PDB nell'ordine completo del template.

Gli output GNINA ordinari non garantiscono gli H richiesti: la documentazione
nel sorgente di `stripH` specifica che gli H non polari vengono sempre rimossi.
Non esiste una flag qui introdotta per mantenerli. Uno SDF GNINA privo di H
va completato **esternamente**, mantenendo identità, protonazione e geometria
stereochimica del template, prima di questo assemblaggio.
[Sorgente ufficiale GNINA](https://github.com/gnina/gnina/blob/master/gninasrc/main/main.cpp).
Lo stesso vincolo vale per i PDBQT che hanno perso H.

Questo comando non effettua allineamento, assegnazione di cariche,
ricostruzione del ligando da distanze o generazione di H.

## Output

| File | Contenuto |
|---|---|
| `compl.parm7`, `compl.rst7`, `compl.pdb` | Complesso secco completo |
| `rec.parm7`, `rec.rst7`, `rec.pdb` | Recettore, comprese le modifiche covalenti |
| `lig.parm7`, `lig.rst7`, `lig.pdb` | Ligando, oppure insieme dei ligandi simultanei |
| `lig_001.*`, ecc. | Componenti ligando individuali |
| `masks.json` | Maschere di recettore, ligando/i e complesso |
| `manifest.json` | Sorgenti e hash, carica, controlli di conservazione dei parametri |

Le estensioni `.parm7` e `.prmtop` indicano lo stesso formato Amber. I PDB
includono CONECT; la topologia resta il riferimento per tutti i parametri.

L'assemblaggio conserva legami, angoli, torsioni, esclusioni, interazioni 1–4,
parametri LJ, cariche e raggi/screen GB; non aggiunge legami recettore–ligando.
Nomi atomici vicini geometricamente non determinano la connettività.

## Pose alternative o caricamento simultaneo

Una posa alternativa va in un complesso separato. Ripetere le opzioni ligando
significa invece inserire più molecole **contemporaneamente** nello stesso
sistema. Le coordinate devono essere compatibili; sovrapposizioni gravi
vengono rifiutate.

L'insieme `lig.parm7` di più molecole descrive un gruppo ligando congiunto. Una
singola energia di binding per quel gruppo non è l'affinità individuale di
ciascuna copia. Per quelle analisi occorrono partizioni e traiettorie coerenti.

## Compatibilità con la dinamica e il GBSA

I file prodotti sono secchi. Solvatazione, controioni ed eventuale rilassamento
sono passaggi successivi del workflow MD; non vengono avviati da questo
comando. Il soluto covalente deve essere finito e le coordinate devono
contenere molecole intere.

Le componenti devono usare un protocollo coerente per raggi GB, fattori 1–4 e
tipi LJ omonimi. Dati incompatibili sono rifiutati: non viene adattato un
force field al volo. Il controllo numerico di compatibilità non convalida le
cariche della superficie, la correttezza della posa o l'interpretazione
termodinamica del GBSA.

Il legame GO–PEG/amminoacido/peptide/ciclodestrina resta interno al recettore.
Tagliarlo per trattare la modifica come un ligando non covalente cambierebbe
il problema fisico e non è l'analisi ordinaria del contaminante adsorbito.
