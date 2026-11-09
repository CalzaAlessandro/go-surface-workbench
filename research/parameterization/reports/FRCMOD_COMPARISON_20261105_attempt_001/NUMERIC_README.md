# Audit numerico dei frcmod GO

Aprire `numeric_report.html`. Il report è separato dalla relazione integrata e
dall'analisi delle cariche. Tutti gli output numerici sono riproducibili con:

```sh
uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache \
  --python VS_POLL_HIGH_COVERAGE/.venv/bin/python \
  PARAMETERIZATION/reports/FRCMOD_COMPARISON_20261105_attempt_001/build_report.py
```

Lo script usa soltanto la libreria standard Python. Scrive esclusivamente nella
propria cartella; legge i sette frcmod archiviati, i relativi originali snapshot,
i sei DAT locali e le due topologie GAFF2-AIM. Verifica gli hash dei sette file
contro MANIFEST.json, la corrispondenza con gli originali e gli hash dei sei DAT
contro i manifest precedenti. Il lettore tratta `NONB` come `NONBON`, conserva
le righe incomplete senza completarle, e ricompone le serie DIHE tramite il
segno di continuazione. Non importa il parser di produzione e non lo modifica.

- `numeric_inventory.csv`: sezioni, hash, duplicati, scaling DIHE esplicito.
- `numeric_terms.csv`: ogni riga valida con commento, unità, riga originale.
- `numeric_invalid_rows.csv`: ogni riga incompleta, separata dai coefficienti.
- `numeric_duplicates.csv`: ripetizioni della stessa firma, confrontate come
  serie intere per DIHE. Reverse di BOND/ANGLE/DIHE canonico; IMPROPER ordinato.
- `numeric_improper_permutation_overlaps.csv`: screening aggiuntivo con tipo
  centrale in terza posizione e permutazioni dei tre esterni; nessuna deduzione
  automatica dell'effetto nel loader.
- `numeric_base_comparison.csv` e `numeric_base_summary.csv`: confronti contro
  tutte le quattro basi GAFF/GAFF2 conservate, oppure le due basi Amber per il
  namespace uppercase. Il confronto non implica caricamento simultaneo.
- `numeric_pairwise_terms.csv` e `numeric_pairwise_summary.csv`: firme native
  condivise o presenti in un solo override. L'assenza di override non è un gap
  di topologia dimostrato. AMBER uppercase non è associato a GAFF lowercase.
- `numeric_current_earlier_changes.csv`: dettaglio delle differenze locali.
- `numeric_aim_lj.csv`: σ/ε dalle topologie originali e conversione a Rmin/2/ε
  in unità Amber, mantenendo tutte le classi ambientali originali.
- `numeric_aim_bonds.csv`: legami espliciti delle due topologie AIM con numero
  di occorrenze e conversione della costante armonica, incluso il fattore 1/2.
- `numeric_highlights.csv`, `numeric_flags.csv`, `numeric_qm_priorities.csv`:
  risultati essenziali e priorità chimiche, senza classifica di accuratezza.
- `numeric_summary.json`: inventario strutturato, risultati e limiti.

`same_coefficients` indica uguaglianza numerica dopo normalizzazione DIHE
PK/IDIVF e modulo del numero d'onda. Si confrontano anche il termine costante
e i coefficienti seno/coseno dell'intera serie. `same_to_reported_precision`
indica compatibilità con metà dell'ultima cifra stampata più grossolana; non è
una tolleranza fisica. La massa e il campo opzionale di polarizzabilità sono
distinti. Il secondo campo MASS non è una nuova massa atomica.

Il report non risolve wildcard, ordine di caricamento storico, grafo/typing GO,
permutazioni improper effettive o topologia LEaP finale. Per questo non può
certificare copertura chimica, accuratezza o assenza di termini mancanti.
SCEE/SCNB espliciti sono riportati; quelli ereditati devono essere verificati
nel prmtop. Le cariche, LJ, regole di mixing e fattori 1–4 di un preset vanno
considerati insieme. Nessun nuovo QM, MD, fit o aggiornamento di produzione è
stato eseguito.
