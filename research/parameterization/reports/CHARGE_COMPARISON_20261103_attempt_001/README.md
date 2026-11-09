> Aggiornamento successivo: completato anche il terzo C54. Il pannello DFT ha ora 19 contesti. Vedi il [nuovo report con Espaloma e SQE](../CHARGE_CONTROLS_UPDATE_20261103_attempt_001/REPORT.html). Questa cartella conserva lo snapshot precedente a 18 contesti.

# Cariche RESP per Amber: confronto del 3 novembre 2026

Sono completati 16 riferimenti atomici RESP HF/6-31G* su geometrie DFT,
16 corrispondenti controlli DFT–RESP, quattro fit RESP a blocchi condivisi L0
e 12 controlli AM1-BCC con copertura completa dei parametri BCC.
Due ulteriori C₅₄ portano il pannello DFT–RESP a 18 strutture.
Nessuna nuova libreria di produzione è stata selezionata.

HF e DFT condividono le stesse prime 16 geometrie: non sono 34 contesti chimici
indipendenti. COOH e COO⁻ hanno ciascuno un solo contesto fisico nel pannello.

La versione da leggere è la [pagina HTML completa](RESP_AMBER_REPORT.html),
con cariche, metodi, grafici selezionabili e considerazioni sui prossimi passi.

## RESP ed ESP: cosa viene calcolato

**Le cariche nuove sono RESP, ottenute con l'eseguibile RESP di Amber.**
ESP indica il potenziale elettrostatico quantistico usato come dato di ingresso.
Il residuo ESP misura quanto bene le cariche RESP riproducono quel potenziale:
non indica che siano state adottate le cariche ESP non ristrette di Gaussian.
Vedi anche la [documentazione RESP](https://q4md-forcefieldtools.org/RED/resp/).

Il riferimento principale è il potenziale MK in fase gas a HF/6-31G*,
calcolato sulle geometrie ottimizzate ωB97X-D/6-311++G(d,p), seguito dal fit
iperbolico RESP. La geometria DFT e il livello HF del potenziale sono registrati
separatamente. I controlli DFT–RESP usano invece il potenziale ωB97X-D/6-311++G(d,p).
AM1-BCC è un metodo alternativo, esplicitamente separato dalle RESP.

Il fit atomico usa qwt = 0,0005, ihfree = 1, vincoli sulla carica intera e sulle
equivalenze dichiarate. Le tre orientazioni rigide hanno pesi registrati che
compensano la diversa numerosità delle griglie. Sono una sola geometria fisica,
non tre conformeri indipendenti. Non essendoci gruppi CH₂/CH₃ da rifinire,
lo stadio RESP 2 è dichiarato non operativo: non viene simulato un secondo fit.
I dettagli sono in [RESP_PROTOCOL.md](../../docs/RESP_PROTOCOL.md).

L'instabilità esterna RHF→UHF osservata nei precedenti controlli è distinta
da un errore SCF o da un'incompatibilità automatica con Amber. Questo pannello
mantiene il riferimento HF ristretto senza dichiararne validata la stabilità
elettronica. Le geometrie DFT hanno superato i rispettivi controlli.

## L0 e gruppi analizzati

L0 è un **fit RESP con una carica condivisa per ruolo chimico**, indipendente
dalla distanza dagli altri gruppi. Il carbonio interno non funzionalizzato è
fissato a 0 e; ciascun blocco neutro somma a 0 e e il carbossilato a −1 e.
I carboni di ancoraggio di COOH/COO⁻ appartengono ai rispettivi blocchi.
I carboni e idrogeni del bordo sono ruoli distinti dal fondo interno.
È una prova della trasferibilità dell'ipotesi modulare, non un preset approvato.

Il pannello copre:

| Chimica | Unità della libreria | Contesti disponibili |
|---|---|---|
| Scheletro interno e bordo C–H | GGG, GRT | C₁₆, C₂₄ e controlli DFT C₅₄ |
| Idrossili basali | HYD | Due OH cis/trans, vicini/lontani, miscele OH/EPX |
| Epossidi basali | EPX | Singolo, coppia e miscele; coppia C₅₄ DFT recuperata |
| Epossidi al bordo | EPH, EPT | Due controlli C₂₄ dedicati |
| Carbonili | CAB | Frammento con due carbonili |
| Carbossile protonato | CBX | COOH con carbonio di ancoraggio |
| Carbossilato | CBD | COO⁻ con carbonio di ancoraggio |

C₁₆/C₂₄/C₅₄ indicano la dimensione dello scheletro aromatico di partenza;
il carbonio del carbossile è aggiuntivo. I casi chiamati OH_EPX contengono
due OH e un epossido. Questo non è ancora un pannello esaustivo di tutte
le distanze, densità, facce, protonazioni e combinazioni di gruppi del GO.
Lattoni e altre funzionalità proposte non sono stati parametrizzati qui.

## Risultati descrittivi

Si calcola prima la media delle tre orientazioni di ogni struttura e poi la
mediana tra strutture. RRMS è l'errore quadratico relativo sul potenziale;
valori minori indicano una migliore riproduzione di quel riferimento QM.

| Cariche valutate sul medesimo ESP HF | Mediana RRMS, 16 strutture |
|---|---:|
| RESP HF atomiche, specifiche di ciascun frammento | 0,0797 |
| RESP DFT atomiche sugli stessi frammenti | 0,1123 |
| RESP L0 HF, fit condiviso sui 16 casi | 0,1883 |
| Libreria attuale completa | 0,3505 |

I fit atomici sono valutati anche sui dati usati per determinarli. Il vantaggio
di questi fit non dimostra una migliore dinamica o affinità di adsorbimento.
Anche il fit L0 su tutti i 16 casi è un risultato di addestramento.

Una verifica separata esclude dal fit tutte le strutture C₂₄: L0 viene stimato
sui sette contesti con scheletro C₁₆ e applicato ai sette C₂₄ con ruoli coperti.
EPH ed EPT restano non disponibili perché non presenti nell'addestramento.

| Confronto sui medesimi 7 C₂₄ | Mediana RRMS, ESP HF |
|---|---:|
| RESP L0 HF addestrato sui soli C₁₆ | 0,1792 |
| RESP L0 HF addestrato sui 16 casi, quindi includendo i C₂₄ | 0,1777 |
| Libreria attuale completa | 0,3501 |

Questo è un primo controllo fuori dall'addestramento, favorevole al candidato
L0 sul pannello disponibile. Non dimostra trasferibilità a fogli macroscopici,
nuove chimiche o tutte le condizioni di ossidazione.

La tabella delle cariche mostra però compensazioni forti per COO⁻: nel L0 HF
su 16 casi C carbossilico = +1,731213 e, ciascun O = −1,049356 e e carbonio
di ancoraggio = −0,632501 e. Prima dell'adozione occorre verificare sensibilità
ai vincoli, correlazione delle variabili e contesti carbossilato aggiuntivi.
Una decomposizione RESP non misura in modo univoco il trasferimento di carica.

## AM1-BCC e modelli pubblicati

Quattro casi hanno un termine BCC nativo non coperto e sono esclusi dal
confronto principale: pyrene_dione, C24_single_basal_EPX_control,
C24_EPH_edge_control e C24_OH_OH_near_trans. Le tipologie mancanti sono
16–12–2 o 12–16–2. Non sono sostituite con zero né con parametri inventati.
Le cariche native sono conservate; il candidato indicato AM1-BCC* nei grafici
contiene una piccola correzione vincolata dell'arrotondamento della carica totale.
La proiezione L0 delle AM1-BCC è ancora un candidato diverso dal metodo nativo.

Il CSV completo include anche controlli da HierGO/OPLS, GOPY e Stauffer e
la copertura effettiva di ciascun candidato. I confronti con le librerie OLD
di foglio intero e GAFF2-AIM non hanno qui un mapping atomico completo sui
frammenti: **non hanno ricevuto uno score fittizio**. I controlli ibridi che
combinano fonti diverse sono etichettati come tali. La classifica dei soli
residui sul potenziale non seleziona una coppia cariche/Lennard-Jones.

## File da leggere

- [Grafico selezionabile HF/DFT](comparison.html), [PNG HF](HF_ESP.png), [PNG DFT](DFT_ESP.png); gli SVG sono esportabili.
- [Cariche RESP L0 per ruolo](L0_shared_role_charges.csv), con quattro fit separati HF/DFT e addestramento completo/C₁₆.
- [Risultati numerici sintetici](SUMMARY_RESULTS.json).
- [Copertura di tutti i candidati](CANDIDATE_COVERAGE.csv) e [confronti su sottoinsiemi identici](MATCHED_SUBSET_SCORES.csv).
- [Residui per struttura](RRMS_by_physical_context.csv), compresi i casi non disponibili e le rispettive ragioni.
- [Differenze RESP HF–DFT per ruolo](matched_HF_DFT_role_differences.csv).
- Cariche individuali, somme per gruppo, coppie e distanze dai gruppi: [HF, 16 contesti](context_HF_16/context_analysis.json) e [DFT, 18 contesti](context_DFT_18/context_analysis.json), con CSV nelle stesse cartelle.

## Recupero e provenienza

Il job 1452975 ha completato la sintesi complementare. Sono stati recuperati
fit reali già terminati, verificando gli hash degli artefatti e i log MK;
non sono stati rilanciati Gaussian, SQM o Antechamber. Il nuovo fit condiviso
usa RESP da sorgente AmberClassic identificato, con capacità aumentata per
le copie atomiche delle diverse molecole e orientazioni. Il solver è invariato;
il confronto con il precedente pilot del RESP installato passa entro 10⁻⁶ e.

Il job 1453077 ha recuperato i fit DFT–RESP C54_pristine_parent_reference e
C54_EPX_EPX_near_trans. L'header espgen espande il numero di punti oltre
99.999: sia il parser Python sia il lettore Fortran richiedevano una lettura
che non tronchi la cifra aggiuntiva. La modifica dell'eseguibile separato
riguarda la lettura dell'header e il reporting, non le equazioni del fit.
La regressione sul pilot coincide alla precisione di output.
Sono stati mantenuti tutti i punti: 120.100/120.005/120.107 nel primo caso e
122.272/122.477/122.602 nel secondo. Gli input e gli output originari restano
conservati e nessuna installazione di sistema è stata modificata.

Snapshot locali dei risultati in `PARAMETERIZATION/calculations/`:

| Archivio | SHA-256 |
|---|---|
| complementary_results_20261103_snapshot_001.tar.gz | 32ea75cec824cd06e740687215ec4d59267461f462ddca88e1509927c1d46552 |
| c54_recovered_results_snapshot_20261103_attempt_001.tar.gz | 471798f22cd9d6ae5306d9b301c878a466cc11ac03ba32a3b1d04a786cca4759 |

Le directory remote sono isolate sotto
`/export/scratch/alessandro.calza3/GO_PARAMETERIZATION_20261029_fqOFQv/calculations/`.
Il codice della sintesi e del report è conservato nelle rispettive directory
di tentativo, con i manifest degli input e i controlli sull'eseguibile.

All'ultimo controllo del 3 novembre il C54_single_basal_EPX_control
(1406363_3) era ancora in esecuzione, dopo circa 3 giorni e 11 ore, sul limite
di cinque giorni. Il job 1453111 ne attende la conclusione e tenterà il
postprocessing solo dopo aver verificato gli artefatti QM. La dipendenza è
afterany: anche se il padre fallisce, i controlli scientifici restano obbligatori.

Rimangono la scelta delle cariche, i controlli di sensibilità/trasferibilità,
il confronto coerente cariche–Lennard-Jones e le verifiche in solvente e di
adsorbimento. Il nuovo fit bonded/frcmod e la ricostruzione definitiva dei
quattro materiali restano rinviati.

## Decisioni successive confermate il 3 novembre

HF–RESP principale; DFT–RESP e AM1-BCC controlli. La carica totale deriva
dalla chimica e protonazione del prodotto ed è intera, anche non nulla.
La scelta fra libreria fissa, classi locali e redistribuzione statica resta
aperta. GAFF/frcmod rimandati. Il [protocollo del prossimo confronto](../../docs/CHARGE_STRATEGY_BENCHMARK_20261103.md)
include trasferibilità, scala, nuovi giunti e riferimenti sperimentali.

Il [diagnostico COO](COO_charge_diagnostic.csv) confronta la libreria attuale,
RESP molecolare HF e L0: il valore elevato di C compare in L0, mentre il
fit molecolare dà C = +0,807216 e e O = −0,786679 e. Servono confronti
controllati per attribuire la differenza ai vari vincoli. I [carboni interni](basal_C_summary_by_context.csv)
sono riassunti per contesto e distanza; non sono misure di trasferimento
fisico di carica. Non sono stati eseguiti nuovi fit in questo aggiornamento.
