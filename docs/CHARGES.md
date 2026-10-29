# Audit delle cariche e riferimento elettrostatico

L'audit confronta le sorgenti in modo descrittivo. Non assegna un vincitore:
GO CURRENT, i due sistemi OLD e la tabella HierGO non descrivono la stessa
composizione, protonazione, connettività e carica totale. Differenze nelle cariche
parziali, anche grandi, non misurano da sole l'accuratezza di una parametrizzazione.

La provenienza comunicata dall'utente per CURRENT è la parametrizzazione Amber
standard di frammenti GO separati. Metodo QM, base, schema di fitting, equivalenze
RESP, geometrie e vincoli originali restano da recuperare. Non si attribuiscono
retroattivamente HF/6-31G*, RESP o AM1-BCC alle tabelle disponibili.

## Sorgenti e risultati dell'audit

`CURRENT` comprende esclusivamente i nove file `*_gaff2.lib` nella copia
`resources/reference/current`. I nove `*_amber.lib` presenti nella stessa cartella
sono una famiglia separata, `BUNDLED_LEGACY_AMBER`; `EHH.lib` è
`BUNDLED_AUXILIARY`. La cartella condivisa non implica che tutti i file siano
parametri attualmente usati.

Il parser OFF estrae le tabelle atomiche effettive e la connettività interna senza
dipendere da `boundbox` o dalle sezioni di perturbazione. Questo permette l'audit
del CBD storico e di EPH con metadati EPX estranei, ma non certifica la validità
completa di questi file per LEaP. L'identità atomica usa l'indice nel template:
nomi duplicati non vengono fusi.

OLD viene letto dalle topologie compilate `gox.parm7` e `rgo.parm7` tramite ParmEd.
I motivi chimici sono descritti dalla connettività locale, separando gli indizi
dei nomi di residuo. Dove mancano ordini di legame affidabili, le assegnazioni
restano esplicitamente ambigue, per esempio `carbonyl_or_deprotonated_oxygen`.
I vicini di un template incompleto e quelli di una topologia completa non sono
contesti automaticamente equivalenti.

HierGO è fissato al commit `60382a1a1807ef0003ffcc3903392e913b1e5e9a` del
[repository ufficiale](https://github.com/IFM-molecular-simulation-group/HierGO).
Licenza GPL-3.0, file originali, URL e SHA256 sono conservati in
`resources/reference/opls/PROVENANCE.json`; l'audit verifica questi hash.
La tabella originale `Oplsatom2` contiene 17 righe contestuali, non 17 atomi di un
sistema. Lo stesso `opls_180` compare come ossigeno epossidico EPX/OE con −0.4000 e
ed etereo ETH/OT con −0.1700 e: il solo tipo OPLS non identifica la carica.

La selezione esplicita dei residui GO 1–494 del complesso CURRENT disponibile,
escludendo il ligando 495, contiene 1007 atomi (C613 O253 H141), carica −16 e.
OLD gox contiene 940 atomi (C661 O214 H65), carica −9.0821659936 e; OLD rgo
contiene 814 atomi (C674 O47 H93), carica −3.0000199991 e. La carica non intera
di OLD gox merita recupero della procedura originale; non viene corretta
automaticamente e non autorizza una conclusione sull'accuratezza relativa.

```python
from go_builder.charges import audit_charge_sources

report = audit_charge_sources(
    "results/charge_audit",
    current_topology="/path/to/compl.parm7",
    current_go_residue_count=494,
)
```

API: `audit_charge_sources(output_dir, *, current_dir=None, old_dir=None,
opls_dir=None, current_topology=None, current_go_residue_count=494)`.
Le directory di sorgente predefinite sono le copie sotto `resources/reference`.
La topologia CURRENT compilata è opzionale; quando specificata, la selezione
dei primi N residui deve essere verificata sul sistema fornito.

Output: `audit.json`, `atoms.csv`, `motifs.csv`, `sources.csv`. Le statistiche sono
separate per sorgente e contesto atomico; somme di template o tabelle di parametri
non sono cariche totali di una superficie. Ogni file dati ha SHA256 e percorso.
Una directory di output già presente genera `FileExistsError`.

## Preparazione di un confronto verificabile

1. Scegliere frammenti GO finiti, con bordi saturati esplicitamente con H, che
   rappresentino bordo grafenico, OH, epossido, carbonile, COOH, COO− e l'eventuale
   prodotto di apertura. Definire protonazione, carica intera e molteplicità;
   alcuni conteggi di gruppi basali generano specie con elettroni spaiati.
2. Applicare ciascuna sorgente solo dove esiste una corrispondenza chimica
   documentata. Cariche dei caps e dei nuovi giunti richiedono una scelta
   motivata. Se una sorgente non copre il frammento, registrare il caso come
   non assegnabile; non aggiungere zeri né neutralizzare per ottenere un confronto.
3. Usare la stessa geometria e lo stesso stato di carica per tutti i candidati
   di ciascun frammento. Preparare un riferimento QM dichiarato e ripetere il
   confronto su più conformazioni; il riferimento non ricostruisce il metodo
   originario di CURRENT.
4. Confrontare ESP sullo stesso insieme di punti e dipoli rispetto alla stessa
   origine. Successivamente verificare geometria, stabilità Amber e sensibilità
   delle osservabili, compreso MM/GBSA con raggi GB e campionamento confrontabili.

```python
from go_builder.qm import prepare_qm

manifest = prepare_qm(
    capped_model, "results/qm_epoxide", total_charge=0,
    multiplicity=1, method="HF", basis="6-31G*", nproc=4, memory="4GB",
)
```

API: `prepare_qm(model, output_dir, *, total_charge=None, multiplicity=1,
method='HF', basis='6-31G*', nproc=4, memory='4GB')`. Se manca `total_charge`,
devono essere presenti le cariche formali di tutti gli atomi. Coordinate in Å.
Sono rifiutati modelli periodici, legami con immagini residue, bordi aperti,
sovrapposizioni atomiche e combinazioni carica/spin incompatibili.
Non vengono tagliati o saturati automaticamente frammenti periodici.

Vengono scritti `reference.gjf`, `qm_manifest.json`, `resp-stage1.in`,
`prepare_resp.sh` e istruzioni. Il file Gaussian è un single point, senza
ottimizzazione, con `Pop=MK`, SCF stretto e output dei punti ESP.
Il riferimento predefinito HF/6-31G* è una scelta esplicita di benchmark.
Gaussian e AmberTools restano software esterni: **nessun calcolo viene avviato**.
Lo script RESP può essere lanciato dall'utente dopo la verifica dell'output QM e
rifiuta di sovrascrivere risultati. Il fitting RESP a uno stadio fornito è
esplorativo, senza equivalenze: il fitting produttivo richiede vincoli chimici
e l'eventuale secondo stadio definiti per il protocollo scelto.

## Da espgen al confronto comune

```python
from go_builder.qm import espgen_to_reference, compare_esp

espgen_to_reference(
    "results/qm_epoxide/reference.esp", capped_model,
    "results/epoxide_reference.json", total_charge=0,
    # Facoltativo: dipolo QM nella stessa orientazione e origine.
    dipole=[0.0, 0.0, 1.2], dipole_unit="debye",
    dipole_origin=[0.0, 0.0, 0.0],
)
compare_esp(
    "results/epoxide_reference.json",
    {"CURRENT": current_matched_model, "OLD": old_matched_model,
     "HIERGO": hiergo_matched_model},
    "results/epoxide_esp",
)
```

`espgen_to_reference(esp_path, model, output_path, *, total_charge,
geometry_tolerance=1e-5, dipole=None, dipole_unit='debye', dipole_origin=None)`
legge l'ESP standard Amber per una molecola, con coordinate in bohr e potenziale
in unità atomiche. La geometria deve corrispondere all'ordine degli atomi del
modello entro la tolleranza in Å. Con `NoSymm` si conserva l'orientazione richiesta;
in caso di rotazioni o riordini il convertitore rifiuta il dato e richiede una
mappatura esplicita. `dipole_origin`, se fornito al convertitore, è in Å.

La struttura JSON minima per un riferimento analitico di esempio è:

```json
{
  "atoms": [{"id": 7, "element": "H", "position": [0, 0, 0]}],
  "total_charge": 1,
  "position_unit": "bohr",
  "potential_unit": "hartree/e",
  "grid": [{"position": [2, 0, 0], "potential": 0.5}],
  "dipole": [0, 0, 0],
  "dipole_unit": "debye",
  "dipole_origin": [0, 0, 0]
}
```

Un file candidati per la CLI ESP usa questa forma, con coordinate predefinite in Å:

```json
{
  "analytical_unit_charge": {
    "atoms": [{"id": 7, "element": "H", "position": [0, 0, 0], "charge": 1}]
  }
}
```

Questo esempio verifica unità e legge di Coulomb; non rappresenta un benchmark
QM completato né un modello GO. Per dati reali occorre fornire le cariche
corrispondenti di tutti gli atomi, caps compresi.

API: `compare_esp(reference_path, candidates, output_dir, *,
geometry_tolerance=1e-6, charge_tolerance=1e-6)`. Vengono verificati identità,
elementi, geometria, finitezza, non periodicità e uguaglianza della carica totale.
Nessuna rinormalizzazione o compensazione delle cariche viene applicata.
Il riferimento richiede unità esplicite; supporta coordinate bohr/Å e dipolo
debye, e·bohr o e·Å. Nel JSON di riferimento l'origine usa `position_unit`.

Output: `esp_comparison.json` e `esp_residuals.csv`, con RMSE, MAE, errore massimo,
RMSE relativa, dipoli e rispettivi errori se è dato il dipolo QM. Potenziali in
Hartree/e; la conversione kcal/mol/e è l'energia per carica di prova unitaria.
Il dipolo di un sistema carico dipende dall'origine. Le RMSE sono non pesate:
la distribuzione dei punti del riferimento determina il peso spaziale.
Hash del riferimento e di ogni insieme candidato carica/geometria permettono
di ricostruire l'esatto confronto. Non viene prodotta una raccomandazione globale.

Riferimenti operativi primari: il formato espgen, comprese coordinate in bohr e
campi a larghezza fissa, è verificato nel
[sorgente Amber ufficiale](https://github.com/Amber-MD/AmberClassic/blob/main/src/antechamber/espgen.c).
L'uso di espgen e respgen per fitting a due stadi è descritto nel
[manuale Amber](https://dev-archive.ambermd.org/200603/att-0122/amber9.pdf).
