# Direzione scientifica del progetto

Questa è una proposta di sviluppo, non una dichiarazione di novità assoluta o
un elenco di funzioni già implementate. Le fonti GOPY/HierGO restano citate in
[COMPARISON.md](COMPARISON.md); provenienza e licenze non cambiano aggiungendo
funzioni o riscrivendo nomi. La verifica della letteratura specifica andrà
ripetuta quando si definirà un contributo da pubblicare.

## Base già presente

Il modello principale è un grafo esplicito con ID atomici stabili e immagini
dei legami periodici. La scelta dei siti è casuale per default, riproducibile,
a conteggi esatti e con warning sui caricamenti irrealizzabili. Sono presenti
apertura anti degli epossidi con frammenti ammessi, Amber senza bondbydistance,
conversione GROMACS, import di snapshot equilibrati e docking su superfici 3D.
Il supporto geometrico non equivale a parametri MD validati per qualsiasi
innesto. L'export Amber dei fogli covalentemente infiniti resta bloccato.

## Estensioni che proporrei, in ordine

| Priorità | Funzione proposta | Risultato verificabile |
|---|---|---|
| 1 | Ricette inverse da XPS con incertezza | Insieme di conteggi interi compatibili con C/O e componenti spettrali, anziché una composizione univoca non giustificata. Rapporti/frazioni derivati, non duplicati degli input count. |
| 2 | Distribuzione spaziale controllabile | Gruppi dispersi o aggregati, regioni protette, vincoli indipendenti sulle due facce, statistiche di correlazione confrontabili tra campioni. |
| 3 | Catalogo di reazioni covalenti | Template con mappatura atomica, bilancio di atomi/carica e stereochimica: apertura epossidica, poi altre reazioni solo con chimica e parametri documentati. Non un innesto universale senza controlli. |
| 4 | Blocchi parametrizzati e giunzioni | Cariche e parametri con dominio chimico, protocollo QM, ambiente locale e test di trasferibilità; blocco dell'export se manca una giunzione. |
| 5 | Ensemble superficie–adsorbimento | Medesimo grafo su foglio piano e snapshot equilibrati; ID preservati, accessibilità delle due facce e box locali per curvatura, confronto di più semi/conformazioni. |
| 6 | Campagne riprendibili | Ricette e semi in un manifesto, stesso worker locale/Slurm, ripartenza solo da artefatti verificati e riepilogo dei fallimenti. |

La combinazione più utile per il tuo progetto sarebbe ricetta XPS → ensemble
spaziale → funzionalizzazione reattiva → parametrizzazione verificata → docking
su snapshot equilibrati. Serve dimostrare la trasferibilità dei parametri e
quantificare quanto cambia l'adsorbimento tra modelli, oltre a produrre i file.

## Velocità, memoria e semplicità

Il primo intervento misurato è già implementato: controlli di collisione
vettoriali e riuso delle viste del grafo durante la terminazione dei bordi.
[PERFORMANCE.md](PERFORMANCE.md) riporta velocità, memoria e limiti; gli output
restano identici nei casi di compatibilità verificati.

Il passo successivo, se giustificato da fogli più grandi, è un indice spaziale
locale con verifica esatta delle distanze candidate. Eviterei matrici dense
N×N e parallelismo su singoli atomi. Conviene parallelizzare i semi o i modelli
indipendenti, evitando oversubscription con i motori QM/MD.

Il nucleo di costruzione usa già NumPy e importa i backend specialistici solo
quando servono. Valuterei dipendenze opzionali installabili separatamente per
frammenti RDKit, topologie ParmEd e traiettorie NetCDF. Non è ancora una
modifica del pacchetto: va provata con installazioni minime prima di cambiare
le dipendenze. Un solo schema dati e un worker per locale/Slurm devono restare
la base; meno duplicazioni conta più di ridurre le righe a ogni costo.
