# Scelta delle cariche GO: riproducibilità, trasferibilità e scala

Decisioni confermate dall'utente il 3 novembre 2026. Questo è il protocollo del
confronto successivo, non un nuovo metodo già implementato o validato.

Aggiornamento del 4 novembre: il primo fit SQE è completato come studio pilota;
QM sui due fogli interi ed espansione dei frammenti sono in esecuzione.
Il [protocollo SQE](SQE_HF_RESP_PILOT_20261104.md) e il
[report corrente](../reports/QM_SQE_UPDATE_20261104_attempt_001/REPORT.html)
riportano risultati e limiti. Questa pagina resta il piano di scelta tra A/B/C,
non l'approvazione di un metodo finale.

## Obiettivo e perimetro

Derivare un modello di carica riproducibile per GO di dimensione scelta
dall'utente, utilizzabile in Amber con small molecules, proteine e successive
funzionalizzazioni. HF/6-31G* MK seguito da RESP resta il riferimento principale.
DFT–RESP e AM1-BCC sono controlli separati. Il fit GAFF/GAFF2, Lennard-Jones
e frcmod resta rinviato: ora si decide come assegnare le cariche.

Restano i gruppi classici già coperti: GGG/GRT, HYD, EPX, EPH/EPT, CAB, CBX e
CBD. Eteri, lattoni/lattoli e altre chimiche presenti in modelli come HierGO,
oltre a PEG, amminoacidi e peptidi, sono estensioni da parametrizzare, non
tipologie cui assegnare silenziosamente le cariche di gruppi simili.

## Vincolo obbligatorio sulla carica totale

La somma delle cariche parziali deve uguagliare la carica formale intera del
prodotto finale. Non si impone neutralità al GO. La protonazione, i gruppi
carichi e i prodotti della reazione determinano il totale.

Esempio contabile: 16 carbossilati e 4 centri ammonio da +1 e, senza altri
gruppi ionici, danno Q = −12 e. Un PEG neutro non cambia il totale se anche
il giunto formato e il bilancio di protoni sono neutri. Per un nuovo innesto
si usa il grafo del prodotto, con H e cariche formali espliciti: non basta
sommare cariche di precursori che potrebbero perdere atomi o protoni.

GO covalentemente funzionalizzato, molecole adsorbite e controioni hanno
contabilità separata. Ogni componente chimica indipendente mantiene il
proprio stato di carica; l'eventuale neutralità della scatola solvata è un
requisito distinto. Il totale non dipende dalla scelta L0/intorno/redistribuzione.
Si impone durante il fit o per costruzione, non con un riscalamento arbitrario
di tutti gli atomi alla fine.

La carica formale COO⁻ = −1 e non obbliga la somma delle cariche parziali su
quel particolare sottoinsieme di atomi a essere −1 e in un fit molecolare.
La neutralità dei singoli blocchi e il fondo basale esattamente a zero sono
vincoli aggiuntivi del modello modulare; vanno testati separatamente dalla
conservazione obbligatoria del totale.

## Tre famiglie da confrontare

| Famiglia | Assegnazione durante la preparazione | Cosa verificare |
|---|---|---|
| A · Libreria fissa L0 | Una carica per ruolo, comune a ogni occorrenza; baseline già calcolata | Trasferibilità senza classi ambientali, effetto del fondo a zero e dei confini dei blocchi |
| B · Libreria per intorno | Poche classi chimiche con vettori RESP condivisi; tipo dei vicini, distanza di legame, bordo e giunto | Beneficio su frammenti interi esclusi dal fit; copertura e numero di parametri |
| C · Redistribuzione statica | Valori di base più correzioni conservative assegnate al GO completo | Calibrazione RESP, identificabilità, dipendenza da dimensione/geometria e costo di assegnazione |

Tutte producono cariche fisse durante MD e MM/GBSA. La famiglia B non è
una dinamica a cariche variabili. Una dipendenza dalla distanza euclidea
istantanea durante MD richiederebbe un altro modello energetico.

Per B, la distanza nel grafo è il primo descrittore da provare: non cambia
quando il foglio si piega. Tipo dei gruppi, posizione basale/bordo e faccia
locale sono informazioni distinte. L'eventuale distanza 3D entra soltanto
come ipotesi aggiuntiva verificata su geometrie rappresentative; non diventa
automaticamente una correzione universale 1/r.

Per C, correzioni +δ/−δ su legami conservano il totale per costruzione.
Una baseline deve già avere il totale corretto; componenti disconnesse non
si scambiano carica attraverso legami inesistenti. Una regola locale di
incrementi e un solver SQE sono candidati diversi, non sinonimi. SQE+qp è
un precedente pubblicato, ma i suoi parametri pubblicati non sono una
calibrazione RESP per GO.

## Primo confronto: isolare l'effetto dei vincoli

Prima di moltiplicare le classi per distanza, confrontare:

1. L0 attuale: fondo interno a zero e blocchi con carica imposta.
2. Variante con intorno basale esteso, mantenendo controllati i confini di
   carica. Provare intorni di una e due distanze di legame come ipotesi di
   benchmark, non come raggio fisico già stabilito.
3. Variante con carica vincolata su una regione funzionale più ampia anziché
   sul solo blocco originale. Se due regioni si sovrappongono, fondere i
   domini o usare incrementi conservativi: un atomo non appartiene a due
   somme indipendenti costruite come se fossero disgiunte.

Cambiare un'ipotesi alla volta. Le RESP atomiche di ogni frammento rimangono
un riferimento più flessibile, distinto da una libreria applicabile a ogni
nuovo foglio. I valori dei carboni interni di frammenti finiti non sono
automaticamente valori trasferibili alla superficie macroscopica.

### Dato diagnostico sul carbossilato

Nel modello pyrene_COO già calcolato, il fit molecolare HF–RESP dà
COO_C = +0,807216 e e ciascun COO_O = −0,786679 e. La libreria attuale usa
+0,760000 e e −0,790000 e; il L0 HF condiviso dà +1,731213 e e −1,049356 e.
Il carbonio di ancoraggio differisce fra i tre modelli.

Il fit molecolare e L0 non hanno gli stessi vincoli: il confronto non prova
quale vincolo causi da solo la differenza. Richiede il confronto controllato
sopra descritto. In particolare, non dimostra che il protocollo HF–RESP
richieda di assegnare +1,73 e al carbonio carbossilico, né che la libreria
originale sia sbagliata. COOH e COO⁻ hanno ciascuno un solo contesto fisico:
tre rotazioni non costituiscono tre esempi indipendenti.

## Dati, generalizzazione e criterio di scelta

Usare gli stessi riferimenti HF, geometrie, griglie, pesi e convenzioni RESP
per tutti i candidati. Riutilizzare i dati reali già disponibili prima di
lanciare nuovi QM. Valutare residuo del potenziale generato dalle cariche
RESP, sensibilità dei parametri, somme dei domini e risposta al cambio di
dimensione. Analizzare carboni interni e bordi separatamente.

Escludere dal fit interi frammenti, separazioni e dimensioni, non una loro
rotazione. Una classe chimica assente dal training resta non supportata.
La successiva espansione mirata deve includere COOH/COO⁻ su più dimensioni
e vicino/lontano da OH/EPX; quindi nuovi giunti covalenti rappresentativi.
Per i conformeri, distinguere variazioni fisiche da rotazioni numeriche.

La scelta preferisce il modello più semplice che dia risultati stabili
fuori dall'addestramento e tempi/memoria adatti ai fogli previsti. Una
regola per intorno locale non richiede un nuovo HF sull'intero foglio a ogni
generazione. Per un solver globale, il costo e la dipendenza dalla taglia
vanno misurati: non si presume una scalabilità lineare.

## Riproducibilità e scalabilità della costruzione

Separare ricetta chimica, realizzazione geometrica e preset delle cariche.
Registrare seed, dimensioni, conteggi effettivi dei gruppi, protonazione,
connettività e IDs, versione/hash della libreria e carica di ogni atomo.
La stessa ricetta con lo stesso seed deve riprodurre la stessa realizzazione;
seed diversi producono un insieme di superfici per verificare la robustezza.

I gruppi basali scalano con l'area; i siti di bordo con perimetro e bordi
dei difetti. Con soli carbossilati al bordo esterno, mantenere la stessa
frazione COO/C o densità di carica su tutte le dimensioni non è garantito.
La ricetta deve dichiarare quali densità conserva e rifiutare conteggi non
realizzabili, senza inventare siti né modificare la carica per compensare.
Questa è una specifica da verificare, non una nuova funzione già aggiunta.

Per PEG o amminoacidi, caratterizzare separatamente catena ripetuta, terminali
e giunto GO–innesto. Il legame nuovo cambia l'intorno chimico: un modello
libero di PEG/amminoacido non valida le cariche del prodotto covalente.
Le topologie del GO, dei complessi e dei sottoinsiemi MM/GBSA devono ereditare
la stessa assegnazione finale, senza ricalcolarla separatamente.

## Appoggio sperimentale

Le cariche atomiche parziali non hanno un unico valore sperimentale diretto.
Gli esperimenti vincolano composizione, protonazione e proprietà collettive.
Non è quindi disponibile una tabella sperimentale universale che scelga da
sola tra i tre schemi di assegnazione.

- XPS e analisi complementari: ricetta e proporzioni dei gruppi, con
  incertezza di assegnazione; non una carica parziale univoca per atomo.
- Titolazione potentiometrica: quantità di siti dissociabili e risposta
  al pH del campione, utile per definire protonazione e totale del modello.
  [Tao et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10339393/) combinano
  titolazioni e ipotesi strutturali DFT; queste ultime restano interpretazioni.
- Funzionalizzazioni: [Orth et al.](https://doi.org/10.1016/j.jcis.2016.01.013)
  usano titolazione e spettroscopia per caratterizzare gruppi introdotti.
  Non sono misure dirette dei coefficienti RESP del giunto.
- Adsorbimento e calorimetria: [Pandit e De](https://doi.org/10.1021/acs.jpcc.6b11571)
  misurano interazioni GO–amminoacidi; sono un successivo controllo del modello
  completo, con acqua, protonazione e parametri non elettrostatici coerenti.
  Non si adattano ora le sole cariche a un'unica energia di adsorbimento.

Per i materiali LayerOne servono dati del campione corrispondente: risultati
su GO di altra sintesi sono precedenti e controlli di trasferibilità, non
misure delle proprietà del nostro campione. Proteine e peptidi rappresentano
una verifica successiva più complessa; le piccole sonde permettono prima
di separare meglio gli effetti del modello.

## Precedenti metodologici

[R.E.D., Dupradeau et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC2918240/)
documentano fit RESP su più molecole/conformazioni con vincoli per frammenti.
È un precedente per costruire librerie coerenti; non valida da solo il GO.
Il SI [GAFF2–AIM, Pinto et al.](https://doi.org/10.1021/acs.jpcb.2c06037),
S5, usa classi di primi vicini: è un precedente di intorno chimico, distinto
dal nostro target RESP. [SQE+qp](https://link.springer.com/article/10.1186/s13321-021-00528-w)
fornisce una famiglia di redistribuzione conservativa con riferimento di
carica; una calibrazione GO/RESP rimane da eseguire e validare.
