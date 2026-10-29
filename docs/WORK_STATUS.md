> **Nota sulla snapshot:** questo documento conserva lo stato di lavoro registrato il 29 ottobre 2026. La cartella `validation/`, gli archivi di output e i binari citati nei resoconti sono stati esclusi dalla snapshot; i relativi risultati qui riassunti non sono quindi verificabili dai soli file inclusi.

# Stato di lavoro — 29 ottobre 2026

Il focus richiesto dall'utente è creare e funzionalizzare il GO. Il percorso
principale è 01_build_GO.py → 02_open_epoxides.py (opzionale) →
03_prepare_amber.py, con prepare_fragment.py per i reagenti.

## Implementato e verificato

- Grafo esplicito, reticolo finito/xy-periodico, buchi, vacanze, conteggi esatti
  OH/epoxide/carbonyl/COOH/COO, terminazioni H e semi riproducibili.
- Apertura anti con frammento N/O/S supportato, H trasferito, connettività,
  carica formale e identità conservate; normale locale e gestione delle immagini.
- Interfaccia semplificata: --count e --id facoltativo, --seed; random per
  default e frazione calcolata dall'output. Rimossa la selezione per frazione.
  Anche i gruppi classici accettano {count,id}: ID dei carboni riservati prima
  dei siti casuali; per epossidi liste di coppie C-C.
- PlacementError su richieste valide ma non raggiunte: WARNING, INCOMPLETE,
  capacità nota, miglior conteggio provato e tentativi, suggerimenti; exit1,
  nessun modello parziale. Classic fallback128disposizioni; innesti cercano
  siti/regioisomeri/torsioni alternative. Non si deduce un massimo fisico.
- Prove CLI reali: GO con gruppi standard, due PEG su quattro epossidi, lisina
  tramite epsilon-N e Gly-Gly tramite terminale N. Dati e log in
  validation/final/generation_and_random_functionalization. Nuove prove con
  count/id, siti classici, warning lisina, capacità e128tentativi esauriti in
  validation/final/count_id_and_warnings; vecchie ricevute sono storiche.
- Revisione indipendente del nuovo selettore: corretto un caso di input ID
  non hashabile; riproducibilità, selettori e contabilità senza altri difetti.
- 258 test passati; Ruff pulito (33file). 14 entry point con help/Slurm, worker Slurm
  eseguito localmente produce lo stesso modello del comando seriale.
- Amber LEaP + Sander GB su cluster: job1374053 COMPLETED0:0, 179atomi,
  232legami, q−1; graph/q/coordinate preservati. GROMACS roundtrip numerico
  verificato; nessun grompp o confronto energetico fra motori.
- Audit finale CURRENT GAFF2/OLD/OPLS, provenienza e 51copie originali
  verificate. Gaussian/RESP preparati, non eseguiti.

## Integrazioni già presenti

04 snapshot nonplanare e box3D; 05 Vina/GNINA-Vina senza CNN;
06 assemblaggio topologie complete postdocking; 07 export GROMACS.
Sono integrazioni successive al costruttore. Prove esterne isolate: GNINA1.3.2
su GO reale+etanolo; MMPBSA1frame+idecomp1+EO/DEO, job1378877 COMPLETED0:0,
CSV finiti verificati e scaricati. Gli indici EO0/DEO1 richiedono mappa frame.
L'assemblaggio richiede tutti gli H; gli output GNINA ordinari possono
richiedere ricostruzione esterna. Nessun nuovo sviluppo su questi moduli è
necessario per usare il generatore.

## Manutenzione allineamento autorizzata separatamente

11533file HPC (11520campagna+12template+1helper), 4template locali aggiornati
align→rms con backup. I11520input precedenti erano align first senza maschera.
CPPTRAJ6.29.13 provato numericamente: bare/maskedalign/rms muovono tutti
l'intero frame, erroreligando1.57e−6Å. Nessun motivo per rifare GBSA solo per
quel comando, nessuna traiettoria o analisi di campagna ricalcolata.
Backup HPC e helper scaricati in validation/final/alignment_migration.
README/SHA256SUMS del vecchio VS_POLL_WORKFLOW aggiornati e verificati (147file);
le versioni precedenti dei due documenti sono nel backup locale. autoimage
anchor valido prima di boxnobox; limiti per fogli grandi/ripiegati in ALIGNMENT.md.

## Limiti scientifici aperti

Parametrizzare le nuove giunzioni; eseguire benchmark QM/ESP comune prima
di scegliere cariche CURRENT/OLD/OPLS; rilassare/equilibrare le strutture.
Nessun vincitore delle cariche selezionato. Reticoli covalenti avvolgenti
conservati nel grafo ma export Amber ancora bloccato; percorso testato finito.
Non confondere test di formato o un singolo snapshot GBSA con validazione fisica.

## Ambiente

Cartella di lavoro isolata GO_BUILDER; nessun Git al momento dello stato registrato.
Python3.14.10, uv lock; dipendenze esportate senza editable in requirements.txt.
CLI inseriscono src nel path; API da questa cartella con PYTHONPATH=src.
Software copiato byte-identico da VS_POLL_WORKFLOW: Vina/vina_split1.2.7,
MGLTools1.5.7p1 archivio Linuxx86_64, wrapper e licenze. Amber/Gaussian/GNINA
restano esterni. Non estrarre MGLTools su filesystem macOS case-insensitive.
Documentazione per l'utente: README.md, VALIDATION.md, docs/FRAGMENTS.md,
docs/COMPARISON.md e docs/CHARGES.md.

## Confronto richiesto con GOPY/HierGO

Aggiornato docs/COMPARISON.md in7componenti funzionali:3blocchi geometrici
comuni a entrambi (reticolo, difetti/pulizia, decorazione), periodicità affine
al tema di HierGO, innesti/Amber-ParmEd affini al tema di GOPY. Algoritmi e
parametri diversi; nessuna percentuale di codice attribuita. GOPY tutorial
usa analogieAMBER99SB, HierGO OPLStabulato; qui9librerieGAFF2utente+frcmod.
CinqueoriginaliHierGO archiviati per riferimento; soloOplsatom2letto dall'audit,
nessun programmaHierGO/GOPYeseguito dal generatore. Nuove giunzioni ancora da
parametrizzare, QM/RESP non eseguiti, nessun vincitore delle cariche scelto.
