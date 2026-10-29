# Componenti comuni e differenze rispetto a HierGO e GOPY

Il confronto seguente divide il workflow in **sette componenti funzionali**.
È una classificazione per spiegare le scelte del progetto; non misura una
percentuale di codice condiviso o di accuratezza fisica.

| Componente | GOPY | HierGO | GO Builder |
|---|---|---|---|
| 1. Reticolo | Generazione di grafene 2D | Celle/piastrelle replicabili | Reticolo rettangolare analitico, quattro C per cella |
| 2. Difetti | Costruzione di buchi e pulizia | Vacanze, buchi e difetti topologici | Vacanze/buchi e pulizia dei pendenti; nessuno Stone–Wales o bordo di grano |
| 3. Gruppi ossigenati | COOH, epossidi, OH e varianti | Decorazione basale/bordi, modelli fresh/aged | Conteggi OH/epossidi/carbonili/COOH/COO; casuale o siti fissati |
| 4. Innesti molecolari | Procedura PEG dedicata | Il percorso confrontato riguarda la decorazione del GO | Frammenti espliciti; apertura anti e conservazione O come OH |
| 5. Periodicità | Nessun percorso equivalente verificato qui | Costruzione e stitching di piastrelle | Immagini dei legami nel grafo; export Amber dei reticoli infiniti bloccato |
| 6. Parametri per MD | Librerie/frcmod nel percorso Amber | Atom typing e assegnazione OPLS | Nove librerie GAFF2 dell'utente e GO_gaff2.frcmod; LEaP con legami espliciti |
| 7. GROMACS | Conversione Amber tramite ParmEd | Percorso OPLS diretto | Conversione della medesima topologia Amber e confronto numerico nel round trip |

Le prime tre componenti condividono lo scopo con entrambi gli strumenti.
Con HierGO si aggiunge il tema della periodicità; con GOPY quello degli
innesti PEG e del percorso Amber/ParmEd. I dettagli degli algoritmi e dei
parametri restano distinti.
Fonti: [GOPY](https://github.com/Iourarum/GOPY),
[HierGO](https://github.com/IFM-molecular-simulation-group/HierGO).

## Costruzione e selezione dei siti nel nuovo codice

Il motore genera gli atomi e i legami insieme. I gruppi standard usano siti
casuali salvo ID espliciti, con prenotazione dei carboni e ricerca a tentativi.
Ogni tentativo sceglie la disposizione completa e poi verifica i gruppi in
sequenza; fino a 128 tentativi, senza enumerare tutte le permutazioni. Per gli
innesti si sceglie count e optional id; la frazione è calcolata dal conteggio.
Errori di capacità o ingombro danno warning e stato INCOMPLETE, senza modelli
parziali. Il modello è una geometria iniziale, non un ensemble equilibrato.

GOPY sceglie casualmente gruppi da rimuovere e costruisce una catena PEG
specifica. Qui si conserva l'ossigeno dell'epossido come OH durante l'apertura
anti con un frammento fornito: è una trasformazione diversa.
[Sorgente GOPY](https://github.com/Iourarum/GOPY/blob/master/GOPY.py).

## Cariche e parametri: origine e destinazione

Nel tutorial GOPY i parametri dei gruppi sono ricavati per analogia da
AMBER99SB e da studi sul grafene; il percorso carica librerie/frcmod, usa
bondbydistance e converte Amber in GROMACS con ParmEd. Non è il nostro set GAFF2.
[Tutorial originale](https://github.com/Iourarum/GOPY/blob/master/Tutorial_GO_topology.md).

Nel percorso HierGO, atom typing e Oplsatom2 assegnano tipi/cariche tabulati;
la guida usa VMD/TopoTools e OPLS-top.f90 per l'ITP, con i parametri OPLS
bonded/nonbonded da fornire. Non è un fitting quantistico nuovo della superficie.
[Guida ufficiale](https://github.com/IFM-molecular-simulation-group/HierGO/blob/master/Generating%20an%20OPLS-AA%20topology%20file%20from%20a%20pdb%20file.docx),
[generatore OPLS](https://github.com/IFM-molecular-simulation-group/HierGO/blob/master/OPLS-top.f90).

Nel nuovo percorso MD le cariche e i tipi vengono **letti dalle nove librerie
GAFF2 dell'utente**, con GO_gaff2.frcmod. LEaP riceve gli atomi e tutti i legami:
nessun bondbydistance. Non vengono ricalcolate cariche nuove su ogni foglio,
né imposte neutralizzazioni per ottenere un totale desiderato.

Per una giunzione GO–amminoacido/PEG ecc. i parametri dell'epossido intatto non
sono riutilizzati automaticamente. Servono tipi, cariche e termini della
chimica modificata; la preparazione MD resta bloccata senza un template
completo con provenienza. Questo limite non è risolto dal solo fatto che la
geometria del prodotto sia stata costruita.

OPLS e OLD sono **riferimenti per l'audit**, non alternative già selezionabili
per generare la MD. Non si sostituiscono le sole cariche OPLS dentro GAFF2:
anche potenziali bonded, Lennard–Jones, regole di combinazione e fattori 1–4
fanno parte del modello. Gaussian/MK e RESP sono preparati, ma il benchmark
QM/ESP comune non è stato eseguito. Nessuna famiglia è stata dichiarata migliore.
Vedere [CHARGES.md](CHARGES.md).

## Che cosa è riutilizzato materialmente

Il generatore e il percorso Amber non eseguono i programmi GOPY/HierGO.
Le librerie correnti sono copie tracciate dei file dell'utente, conservate in
resources/reference/current. Cinque file originali HierGO sono archiviati
separatamente con commit/hash/licenza: LICENSE, README.md, atom-typing.py,
OPLS-top.f90 e Oplsatom2. L'audit legge la tabella Oplsatom2 e controlla gli hash;
la costruzione MD GAFF2 non usa quelle cariche né esegue quei programmi.

Non sono implementati tutti gli strumenti originali: in particolare
Stone–Wales, mini grain boundary, stitching multistrato, bunching dedicato,
modelli fresh/aged e N-doping richiedono estensioni proprie.
