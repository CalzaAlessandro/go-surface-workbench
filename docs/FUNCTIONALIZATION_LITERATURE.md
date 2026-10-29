# Funzionalizzazione: riferimenti e scelte del costruttore

Moro et al. (2024) descrivono GO modificato con L-lisina, L-arginina e
L-metionina tramite apertura degli epossidi, in condizioni di sintesi a pH 9.
La composizione XPS guida il modello; la funzionalizzazione comporta anche
parziale riduzione del materiale. Il lavoro impiega GAFF, acqua TIP3P e MD in
Amber 16; per il glifosato riporta cariche RESP. Questi dati sostengono la
classe di trasformazione, senza determinare automaticamente il sito reattivo,
la protonazione o le cariche di ogni nuovo frammento del nostro generatore.
Fonte: [articolo e record ufficiale UNIBO](https://cris.unibo.it/handle/11585/983516),
DOI [10.1016/j.mtchem.2024.101936](https://doi.org/10.1016/j.mtchem.2024.101936).
È stato consultato il testo principale; il supplemento DOCX non è stato
recuperato durante questa verifica.

Tunioli et al. (2023) descrivono β-ciclodestrine ammino-funzionalizzate con
linker C0/C2/C6/C12 e il loro innesto sul GO tramite apertura degli epossidi.
Il controllo con βCD non modificata distingue la reattività dei precursori
amminati nelle condizioni studiate. I modelli sono guidati da XPS; vengono
riportati GAFF, cariche ottenute a livello AM1, TIP3P, Amber 16 e analisi
MM/GBSA. Lunghezza del linker e composizione del materiale sono variabili
dello studio: non forniscono parametri trasferibili a ogni ciclodestrina o
frammento arbitrario.
Fonti: [articolo](https://chemistry-europe.onlinelibrary.wiley.com/doi/full/10.1002/chem.202301854),
[record ufficiale e supplementi](https://cris.unibo.it/handle/11585/948766).

## Regole implementate

Il nostro costruttore richiede un sito N/O/S esplicito, neutro e con un H
trasferibile. Rompe un legame C–O dell'epossido, forma il legame C–nucleofilo e
trasferisce quell'H all'ossigeno conservato. L'attacco anti e un solo ancoraggio
per apertura sono regole costruttive richieste per questo progetto; non sono
presentate come una traiettoria di reazione osservata nei due articoli.

Per lisina, PEG, peptidi o βCD si sceglie l'atomo tramite la mappa del frammento.
α-N ed ε-N della lisina restano scelte distinte. Protonazione finale e nuova
giunzione richiedono parametri propri prima della MD; il programma non
trasferisce automaticamente le cariche dell'epossido intatto al prodotto.
