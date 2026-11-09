# Frammenti di riferimento per gruppi GO aggiuntivi

L’estensione prepara sette grafi molecolari C/H/O neutri e geometrie iniziali riproducibili per confronti QM futuri. Non aggiunge tipi o cariche al builder, non modifica L0 e non esegue QM o dinamica molecolare. Tutti i ruoli iniziano con `extended_` e hanno stato `UNPARAMETERIZED_RESEARCH_REFERENCE`. Il formato è compatibile con `go_param.resp._validate_model` e con gli helper degli input Gaussian; questo non significa che `modular.partition_model` supporti i nuovi ruoli.

| Modello | Motivo | Interpretazione del riferimento |
|---|---|---|
| `dibenzofuran_ether_edge` | O in anello a cinque termini, non epossido | Piccolo motivo aromatico fuso per un possibile ambiente etereo di bordo/difetto; l’O eteroaromatico non è equivalente a ogni etere alifatico |
| `diphenyl_ether_bridge_control` | Ar–O–Ar non ciclico | Controllo separato di ponte diarile, non un frammento continuo di foglio |
| `phthalide_lactone5_edge` | Estere ciclico a cinque termini | Riferimento aromatico fuso per lattone; comprende il carbonile e l’O estereo dello stesso anello |
| `coumarin_lactone6_edge` | Estere ciclico coniugato a sei termini | Riferimento aromatico fuso; non si trasferiscono automaticamente cariche dal lattone a cinque termini |
| `oxolan2ol_lactol_control` | C sp³ con OH e OR ciclico | Controllo alifatico di emiacetale ciclico; non è un lattone e non prova un ambiente aromatico di bordo |
| `phenol_aromatic_OH_control` | C aromatico–OH | Controllo del fenolo neutro, distinto dall’alcol basale sp³ e dal fenolato |
| `naphthoquinone14_edge` | Due carbonili para coniugati | Riferimento aromatico fuso per un contesto chinonico più ampio del singolo carbonile |

Le molecole sono riferimenti di connettività, non sostituti equivalenti di fogli completi. Nessuna presenza o abbondanza nei campioni XPS dell’utente è assegnata. Un grafo neutro con zero elettroni radicalici RDKit ed elettroni pari permette un’ipotesi iniziale di singoletto chiuso, ma non valida lo stato fondamentale o la stabilità del determinante.

## Rapporto con HierGO e letteratura

L’inventario locale verificato e la [revisione dei gruppi mancanti](MISSING_GROUPS_LITERATURE.md) distinguono ETH di HierGO, etere non epossidico, da EPX/EPH/EPT. HierGO non fornisce un set esplicito di lattoni o COO−: i lattoni preparati qui sono un’estensione di ricerca ispirata alla letteratura GO, non un’implementazione stabilita di HierGO. Anche il lattolo rimane un candidato distinto.

La revisione locale collega la motivazione GO a Gao et al., DOI [10.1038/nchem.281](https://doi.org/10.1038/nchem.281), Bagri et al., DOI [10.1038/nchem.686](https://doi.org/10.1038/nchem.686), e Aliyev et al., DOI [10.3390/nano9081180](https://doi.org/10.3390/nano9081180), con le rispettive condizioni di materiale e limiti delle assegnazioni. Le identità dei piccoli riferimenti non sono identità sperimentali dei siti GO. Le fonti SMILES e i rispettivi identificatori sono registrati nel catalogo e nel README del rapporto.

## Protocollo e gate

1. Sanitizzazione RDKit e verifica della formula e dell’InChIKey del grafo documentato. Aggiunta esplicita degli H, carica formale totale 0, zero radicali nel grafo, elettroni pari, ID atomici contigui e classi di equivalenza del grafo senza atom map.
2. Gate specifici: l’O etereo non appartiene a un triangolo né è l’O estereo di un carbonile; nel lattone C carbonilico e O estereo appartengono all’anello richiesto; nel lattolo un C sp³ ciclico lega OH e OR e il grafo non contiene carbonili; OH fenolico su C aromatico; chinone con due siti carbonilici para e coniugazione nello stesso anello a sei termini.
3. ETKDGv3 con seme 20261105 e un thread, quindi MMFF94s fino a convergenza (massimo 3000 iterazioni). Sono geometrie iniziali, non geometrie QM accettate. Dove lo SMILES lascia indefinita la stereochimica, si registra un solo stereoisomero deterministico dalla geometria, senza assegnare popolazioni o racemati.
4. `geometry_opt/reference.gjf`: helper esistente per ottimizzazione, frequenze e stabilità a wB97XD/6-311++G(d,p). `HF_RESP_template/reference.gjf`: solo modello di input HF/6-31G* MK/Stable a geometria fissata. **Il template HF contiene ancora le coordinate MMFF e va rigenerato dalla geometria DFT accettata prima del suo uso scientifico.** Nessuna orientazione o sottomissione è preparata automaticamente.
5. Prima del fit servono evidenze reali di ottimizzazione/frequenze/stabilità e verifica che connettività, ordine di legame, stereochimica e identità atomica siano conservati. Il fit RESP successivo conserva la carica molecolare; non inventa a priori cariche atomiche o neutralità dei sottoinsiemi.
6. Prima di qualsiasi assegnamento produttivo servono ownership disgiunta dei blocchi e degli anchor, ambiente locale esplicito, validazione su contesti esclusi dal fit e cariche accoppiate a bonded/LJ/scaling coerenti. Vanno controllati torsioni e barriere, geometrie d’equilibrio, interazioni non legate e proprietà MD pertinenti. Non si riusano automaticamente i parametri EPX per ETH, HYD per fenolo o COOH per lattone.

## Trasferimento da circa 4 nm a circa 10 nm

La trasferibilità di dimensione è un requisito futuro del metodo finale di carica, anche per simulazioni con proteine/DNA. Questi piccoli riferimenti e i loro input non la dimostrano. Una volta validata, l’assegnazione deve riutilizzare le stesse classi di grafo locali e gli stessi vincoli di carica totale/gruppo a entrambe le dimensioni, senza introdurre dipendenze dal nome o dal numero totale di atomi del foglio.

Occorre una validazione tenuta fuori dalla calibrazione su fogli più grandi, includendo ossidazione bassa e alta, variazioni di bordo/interno, combinazioni e distanze dei gruppi, protonazione e geometria. Dimensione e grado di ossidazione vanno variati separatamente: stessi conteggi su un foglio più grande non significano stessa copertura. Si devono verificare somma di carica esatta, copertura delle classi, assenza di nuovi ambienti non supportati e stabilità delle assegnazioni per ambienti equivalenti; le proprietà elettrostatiche e MD pertinenti vanno confrontate con riferimenti indipendenti disponibili. Il protocollo non richiede un nuovo QM del foglio intero per ogni dimensione produttiva. Un ambiente nuovo o una mancata validazione richiede un’estensione mirata dei riferimenti, non una promozione automatica del piccolo frammento a prova di trasferibilità. Vedere [i criteri di accettazione](SCALABILITY_ACCEPTANCE_20261105.md).

## Implementazione e verifica

API isolata: `go_param.extended_fragments.build_candidate`, `validate_motif` e `export_candidates`. L’esportazione rifiuta un percorso esistente e costruisce/valida l’intero pannello prima di scrivere i file. I test sono stati scritti prima dell’implementazione e comprendono grafi reali errati, ricostruzione SDF/ID, determinismo della geometria e preparazione dei template senza eseguire QM.

Per una nuova esportazione, dalla directory `PARAMETERIZATION`:

```sh
PYTHONPATH=src uv run --no-project --cache-dir /private/tmp/go-builder-uv-cache --python ../VS_POLL_HIGH_COVERAGE/.venv/bin/python -m go_param.extended_fragments --output-dir reports/EXTENDED_GO_GROUPS_20261105_attempt_003 --cpus 8 --mem-gb 6
```

Usare una nuova directory per ogni tentativo. Il pannello attuale è `attempt_002`, con rappresentazione aromatica coerente fra JSON e SDF dopo MMFF94s. `attempt_001` conserva la prima esportazione storica, sostituita per l'incoerenza di rappresentazione della cumarina; non contiene risultati QM.
