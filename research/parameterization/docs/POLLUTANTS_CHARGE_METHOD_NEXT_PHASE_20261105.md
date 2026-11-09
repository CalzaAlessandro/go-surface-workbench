# Cariche per il successivo database dei pollutants

Nota del 5 novembre 2026, su indicazione dell'utente. Non modifica il database,
le librerie di produzione o i parametri dei calcoli già avviati.

RESP è il modello di carica di riferimento nello sviluppo sia di GAFF sia
di GAFF2. Non esiste una divisione esclusiva «RESP per GAFF, BCC per GAFF2».
AM1-BCC classico può essere usato con entrambi. ABCG2, abbreviazione di
AM1-BCC-GAFF2, è una riparametrizzazione delle correzioni BCC specifica per
GAFF2, mirata soprattutto alle energie libere sperimentali di solvatazione.

Per la fase successiva del database:

1. Valutare **GAFF2 + ABCG2** come candidato per piccole molecole organiche,
   con **GAFF2 + RESP** come confronto di riferimento.
2. Registrare separatamente il metodo di carica e quello dei tipi/parametri:
   `antechamber -c abcg2 -at gaff2` seleziona ABCG2; `-c bcc -at gaff2`
   seleziona il BCC classico con tipi GAFF2. La seconda combinazione non
   diventa automaticamente ABCG2.
3. Verificare protonazione, tautomeri, stereochimica, carica totale intera,
   convergenza AM1 e copertura effettiva delle correzioni. Le evidenze più
   ampie per ABCG2 riguardano soluti neutri; per ioni e specie non ordinarie
   occorre una verifica specifica.
4. Confrontare molecole rappresentative delle famiglie del database,
   includendo idratazione e interazioni con acqua/GO. I risultati di
   solvatazione non dimostrano da soli accuratezza di adsorbimento o GBSA.
5. Conservare versione e hash di eseguibili, tabelle BCC e base GAFF2,
   geometrie e protocollo RESP, modello d'acqua e convenzioni non bonded.
   Nessuna sostituzione retroattiva di cariche nei sistemi già studiati.

La scelta definitiva rimane aperta. Non sono stati eseguiti qui benchmark
MD dei pollutants con ABCG2. I test sui frammenti GO hanno scopo diverso:
confrontano l'elettrostatica con griglie QM, senza validare il force field
completo in soluzione.

Fonti primarie:

- [He et al., 2020, sviluppo di ABCG2](https://pmc.ncbi.nlm.nih.gov/articles/PMC7728379/).
- [He et al., 2025, estensione e confronto con GAFF2/RESP](https://pmc.ncbi.nlm.nih.gov/articles/PMC11948320/).

