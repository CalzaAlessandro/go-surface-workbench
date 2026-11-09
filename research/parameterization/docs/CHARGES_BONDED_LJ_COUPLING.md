# Cariche, Lennard-Jones e bonded: accoppiamento e protocollo

Aggiornamento del 3 novembre 2026. Questa nota descrive la scelta del modello;
non introduce nuovi coefficienti GAFF, LJ o frcmod.

## Perché sono accoppiati

L'energia MM contiene termini bonded, Coulomb e Lennard-Jones. Le osservabili
dipendono dalla loro somma. Cariche e LJ contribuiscono insieme a idratazione,
orientazione dell'acqua, distanze di contatto e adsorbimento. Il cambiamento
di un termine può essere parzialmente compensato da un altro: un singolo
valore di adsorbimento non identifica univocamente cariche e LJ.

Anche i contributi intramolecolari elettrostatici e LJ, comprese le interazioni
1–4 e le rispettive convenzioni, entrano nel profilo energetico di angoli e
torsioni. I coefficienti bonded derivati nel contesto di un modello di carica
possono compensarne parte della superficie energetica. Questo non implica
che ogni modifica di carica richieda di rifittare tutti i bonded; richiede
di verificarne la coerenza, soprattutto per nuovi giunti covalenti.

Il lavoro [Links between the charge model and bonded parameter force constants
in biomolecular force fields](https://pmc.ncbi.nlm.nih.gov/articles/PMC5669220/)
esamina questo legame nei force field proteici. È un precedente metodologico,
non una misura dell'entità dell'effetto nel nostro GO.

## Cosa hanno fatto in GAFF2–AIM

Il paper e il Supporting Information di
[Pinto et al.](https://doi.org/10.1021/acs.jpcb.2c06037) descrivono:

| Parte | Origine dichiarata |
|---|---|
| Bonded | Parametrizzazione GAFF2 di frammenti GO tramite ACPYPE, riutilizzata per i fogli |
| Cariche | Densità elettronica DFT dei fogli, partizionata con DDEC6/Chargemol |
| LJ | Volumi atomici dalla stessa densità, con la procedura Tkatchenko–Scheffler descritta nel SI |
| Classi di intorno | Parametri mediati per classi definite dai primi vicini |

Il SI, S1 e S4–S5, non descrive un nuovo fit bonded AIM su scansioni QM.
Le etichette ambientali dei tipi atomici sono distinte dall'origine dei
coefficienti bonded. Nei file GROMACS pubblicati i termini bonded sono
espliciti; la cartella amber99sb-ildn.ff-modified è descritta dagli autori
come contenente parametri GO ottenuti da ACPYPE.

La derivazione coordinata di cariche e LJ è diversa da un fit simultaneo di
tutti i coefficienti contro esperimenti. Il paper verifica il modello completo
mediante proprietà dell'acqua e adsorbimento di amminoacidi. Non stabilisce
che DDEC6 abbinato a LJ arbitrari conservi quelle prestazioni.

Fonti locali preservate: references/papers/Pinto2023_SI.txt, S1 e S4–S5;
references/gaff2_aim/README.md e GO_30/GO_sheet_3_0.top.
L'ispezione non costituisce un confronto numerico di ogni coefficiente
pubblicato contro ogni versione del file gaff2.dat.

## Conseguenza per il nostro procedimento

Studiare prima le cariche ha senso per definire il riferimento HF–RESP,
verificare l'effetto dei vincoli e scegliere un'assegnazione trasferibile.
Confronti sul potenziale e sulle cariche sono test di elettrostatica, non una
selezione definitiva del force field per MD.

Il successivo confronto delle osservabili deve mantenere identificata la
coppia cariche–LJ, il modello d'acqua e le convenzioni 1–4. Un candidato
RESP con GAFF2 e un preset GAFF2–AIM devono essere valutati nel proprio
contesto; le combinazioni fra fonti sono controlli ibridi dichiarati.

Sequenza proposta:

1. Derivare e verificare le cariche con HF–RESP, conservando i controlli.
2. Conservare bonded/LJ di riferimento dichiarati durante i confronti
   diagnostici, senza attribuire loro validazione finale.
3. Verificare e, se necessario, calibrare la coppia cariche–LJ con più sonde,
   geometrie di interazione e osservabili pertinenti.
4. Verificare successivamente geometria e profili energetici bonded, con
   attenzione a epossidi, bordi e nuovi giunti GO–innesto.

La scelta del modello di carica e la scelta definitiva dei parametri restano
distinte fasi di ricerca, con validazione congiunta prima della produzione.
