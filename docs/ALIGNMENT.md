# Immagini periodiche e allineamento delle traiettorie

Per il GO finito del workflow classico:

```text
parm compl_solv.parm7
trajin md_10ns.nc
strip :WAT
strip :Cl-
strip :Na+
autoimage anchor :1-494
box nobox
rms first :1-494
trajout md_10ns_imaged.nc
run
```

Per REGOX sostituire la maschera con `:1-628`. I conteggi valgono per i sistemi
corrected esistenti; nuovi modelli richiedono le proprie maschere. Gli strip
mostrati non rinumerano il GO perché nei sistemi osservati acqua e ioni lo
seguono. Se si cambia l'ordine dei residui verificare la maschera dopo strip.

`autoimage :1-494` e `autoimage anchor :1-494` sono sintassi equivalenti.
L'ancora centra il GO; le altre molecole vengono traslate per immagini
periodiche. Il comportamento predefinito tratta i ligandi come molecole fixed,
scegliendo l'immagine più vicina all'ancora quando li avvicina. L'istruzione
richiede la cella e quindi precede `box nobox`. `rms first :1-494` usa il GO
per il fit e applica rotazione/traslazione all'intero sistema; non usare
`nofit` per questa finalità.
Fonti: [autoimage](https://amberhub.chpc.utah.edu/autoimage/),
[RMSD](https://amberhub.chpc.utah.edu/rmsd/).

Autoimage non appiattisce il foglio, né garantisce da solo la ricostruzione di
un foglio covalente già spezzato fra immagini. Per un GO molto esteso o
ripiegato verificare continuità dei legami e immagini dei ligandi: il criterio
rispetto al centro dell'ancora non equivale sempre alla minima distanza dalla
superficie. Non certifica automaticamente i contatti di ogni frame. Conservare
traiettoria e topologia con cella per analisi periodiche; il file senza box
non permette di recuperare a posteriori informazioni periodiche eliminate.
Il nuovo estrattore `04_prepare_docking_receptor.py` usa il grafo covalente per
riunire il recettore prima di preparare i box.

## Verifica del vecchio align

CPPTRAJ 6.29.13 del cluster è stato provato con due frame collegati da una
rototraslazione nota, su GOX + ligando, 1029 atomi. `align first`,
`align first :1-494` e `rms first :1-494` hanno tutti spostato l'intero frame.
L'errore residuo sul ligando è circa 1,57 × 10⁻⁶ Å; gli ultimi due output sono
identici numericamente. Il log del comando senza maschera riporta `*`,
1029 atomi. Evidenza: `validation/final/alignment_behavior/verification.json`,
job Slurm 1378874 COMPLETED 0:0.

La [pagina align](https://amberhub.chpc.utah.edu/align/) descrive un default
`move` limitato alla selezione, che non coincide con il risultato di questa
versione installata. Non si può quindi attribuire una geometria GO–ligando
errata al solo uso della parola align. In particolare i vecchi input di
campagna usavano `align first` senza maschera: non è necessario ricalcolare
GBSA soltanto per quel comando. Il fit GO-specifico ora è esplicito per i
prossimi calcoli; un RMSD risultante può differire da un vecchio fit su tutti
gli atomi senza che cambino le distanze interne.

## Manutenzione applicata il 29 ottobre 2026

- Locale: quattro template attivi `align_{1,10,50,100}ns.in` di
  `VS_POLL_WORKFLOW/resources/LIB`, ora con `rms first MASK`.
- HPC: 11520 input di campagna (5760 GOX e 5760 REGOX), 12 template di runtime
  e un helper per MD del ligando libero; totale 11533 file. L'helper usa
  `rms first *`, perché non contiene GO.
- Tutti i 11520 input di campagna osservati usavano precedentemente il fit
  senza maschera; nessun caso masked-only trovato in questo inventario.
- Non sono state ricalcolate traiettorie o analisi di campagna.

Backup locali: `validation/alignment_local_20261029`. Backup e manifest HPC
scaricati: `validation/final/alignment_migration`; l'archivio originale è
`validation/final/alignment_backup_20261029.tar.gz`. Il backup dell'helper,
aggiunto dopo l'archivio, è scaricato separatamente nella cartella migration.
Origine remota: area di validazione HPC; il percorso assoluto è omesso dalla snapshot.

`update_alignment.py` è un'utilità di manutenzione: senza `--apply` produce
l'inventario; con `--apply` richiede una cartella di backup nuova e conserva
gli originali prima delle modifiche. Non è un launcher di rianalisi.
