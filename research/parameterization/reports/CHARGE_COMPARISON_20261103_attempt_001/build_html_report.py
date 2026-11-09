"""Build a portable, offline scientific report from completed charge artifacts."""
import argparse
import base64
import csv
import hashlib
import html
import json
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('--report', type=Path, required=True)
parser.add_argument('--source', type=Path, required=True)
parser.add_argument('--replace', action='store_true', help='Replace only the generated HTML report')
args = parser.parse_args()
out = args.report.resolve()
source = args.source.resolve()
summary = json.loads((out/'SUMMARY_RESULTS.json').read_text())
fits = json.loads((source/'SYNTHESIS_COMPLETE.json').read_text())['L0_fits']
current = json.loads((source/'l0_hf_all16/current_library.json').read_text())['shared_role_charges_e']
hf = fits['l0_hf_all16']['shared_roles']
dft = fits['l0_dft_all16']['shared_roles']
assert len(current) == len(hf) == len(dft) == 25

def q(value):
    return f'{value:+.4f}'.replace('-', '−').replace('.', ',') if value else '0,0000'

def image_uri(name):
    return 'data:image/png;base64,' + base64.b64encode((out/name).read_bytes()).decode()

names = {
    'scaffold_C': 'GGG · C interno', 'edge_C': 'GRT · C al bordo', 'edge_H': 'GRT · H al bordo',
    'hydroxyl_C': 'HYD · C', 'hydroxyl_O': 'HYD · O', 'hydroxyl_H': 'HYD · H',
    'epoxide_C': 'EPX · ciascun C', 'epoxide_O': 'EPX · O',
    'carbonyl_C': 'CAB · C', 'carbonyl_O': 'CAB · O',
    'COOH_anchor': 'COOH · C di ancoraggio', 'COOH_C': 'COOH · C carbossilico',
    'COOH_O': 'COOH · O carbonilico', 'COOH_OH': 'COOH · O dell’OH', 'COOH_H': 'COOH · H',
    'COO_anchor': 'COO⁻ · C di ancoraggio', 'COO_C': 'COO⁻ · C carbossilico',
    'COO_O': 'COO⁻ · ciascun O', 'EPH_CH': 'EPH · C con H', 'EPH_C': 'EPH · C senza H',
    'EPH_O': 'EPH · O', 'EPH_H': 'EPH · H', 'EPT_C': 'EPT · ciascun C',
    'EPT_O': 'EPT · O', 'EPT_H': 'EPT · ciascun H',
}

def charge_table(roles):
    rows = ''.join('<tr'+(' class="sensitive"' if role.startswith('COO_') else '')+'><th scope="row">'
                   +html.escape(names[role])+'</th><td>'+q(current[role])+'</td><td>'+q(hf[role])
                   +'</td><td>'+q(dft[role])+'</td></tr>' for role in roles)
    return '<div class="scroll"><table><thead><tr><th>Ruolo chimico</th><th>Tua libreria</th><th>RESP L0 HF</th><th>RESP L0 DFT</th></tr></thead><tbody>'+rows+'</tbody></table></div>'

coverage = list(csv.DictReader((out/'CANDIDATE_COVERAGE.csv').open()))
coverage_rows = ''.join('<tr><td>'+html.escape(r['candidate'])+'</td><td>'+r['covered_contexts']
                        +' / '+r['total_contexts']+'</td></tr>' for r in coverage)
hashes = {
    name: hashlib.sha256((out/name).read_bytes()).hexdigest()
    for name in ['SUMMARY_RESULTS.json', 'RRMS_by_physical_context.csv', 'L0_shared_role_charges.csv',
                 'COO_charge_diagnostic.csv', 'basal_C_summary_by_context.csv']
}
important = ['hydroxyl_O', 'epoxide_O', 'carbonyl_O', 'COOH_C', 'COO_anchor', 'COO_C', 'COO_O']
page = '''<!doctype html>
<html lang="it"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>GO · Cariche RESP per Amber · 3 novembre 2026</title>
<style>
@media(prefers-color-scheme:light){:root{--foreground:#182532;--muted-foreground:#506170;--secondary:#f2f5f8;--border:#ccd5de;--ring:#285f9f;--warning:#925a00;--warning-surface:#fff4d6;--radius:8px;--font-sans:system-ui,sans-serif}}
@media(prefers-color-scheme:dark){:root{--foreground:#e4eaf1;--muted-foreground:#afbfce;--secondary:#202c38;--border:#405262;--ring:#8dbfff;--warning:#f4c56a;--warning-surface:#352c1e;--radius:8px;--font-sans:system-ui,sans-serif}}
html,body{margin:0;color:var(--foreground);font-family:var(--font-sans);font-size:15px;line-height:1.55}
main{max-width:1000px;margin:auto}p{margin:.7em 0}h2{font-size:1.22rem;line-height:1.3;margin:1.3em 0 .55em}h3{font-size:1rem;margin:1.1em 0 .4em}a{color:var(--ring)}
.muted,small{color:var(--muted-foreground)}.counts{display:flex;gap:8px 20px;flex-wrap:wrap;margin:12px 0}.counts span{border-bottom:2px solid var(--border);padding:2px 0}.counts b{font-size:1.25rem}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:.88rem}th,td{padding:7px 9px;border-bottom:1px solid var(--border);text-align:right;white-space:nowrap}th:first-child,td:first-child{text-align:left}th[scope=row]{font-weight:500}thead th{font-weight:600;background:var(--secondary)}.scroll{overflow:auto}.sensitive{background:var(--warning-surface)}.notice{border-left:3px solid var(--warning);padding:0 14px;margin:16px 0}.sensitive-note{color:var(--warning)}
details{margin:14px 0}summary{cursor:pointer;font-weight:600;padding:5px 0}details[open] summary{margin-bottom:8px}button{font:inherit;color:var(--foreground);background:var(--secondary);border:1px solid var(--border);border-radius:var(--radius);padding:6px 12px;cursor:pointer}button[aria-pressed=true]{outline:2px solid var(--ring)}.tabs{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}.plot{display:block;width:100%;max-width:850px;height:auto}.note{font-size:.8rem;color:var(--muted-foreground)}li{margin:9px 0}code{font-size:.85em;overflow-wrap:anywhere}pre{overflow:auto;background:var(--secondary);padding:12px;border-radius:var(--radius)}
@media print{details>*{display:block!important}summary{list-style:none}.tabs{display:none}main{max-width:none}table{font-size:9pt}a{color:inherit}.plot{max-width:700px}}
</style></head><body><main>
<p class="muted">GO parameterization · 3 novembre 2026 · risultati preliminari verificati</p>
<p><strong>Stiamo derivando cariche RESP per Amber.</strong> Il candidato modulare L0 migliora la riproduzione del potenziale QM rispetto alla libreria attuale, anche nei primi casi esclusi dall’addestramento. I valori del carbossilato richiedono ulteriori controlli prima di scegliere una libreria.</p>
<div class="counts"><span><b>16</b> RESP HF//DFT</span><span><b>18</b> DFT–RESP</span><span><b>12</b> AM1-BCC coperti</span><span><b>4</b> fit RESP L0</span></div>
<p class="note">HF e DFT condividono le stesse prime 16 geometrie; il pannello DFT aggiunge due C₅₄. Non sono 34 contesti chimici indipendenti.</p>

<h2>Decisioni confermate: riferimento e carica del prodotto</h2>
<p>HF–RESP resta il riferimento principale; DFT–RESP e AM1-BCC sono controlli. La scelta tra libreria fissa, classi di intorno e redistribuzione statica è ancora aperta. La somma delle cariche parziali deve uguagliare la carica formale intera del GO funzionalizzato, determinata da composizione e protonazione, anche se non nulla. L’eventuale neutralità della scatola con controioni è una verifica separata.</p>
<p>Per ora si lavora sulle cariche; GAFF/GAFF2, Lennard-Jones e frcmod restano rinviati. Tutte le opzioni assegneranno un vettore di cariche fisso durante MD. Il <a href="../../docs/CHARGE_STRATEGY_BENCHMARK_20261103.md">protocollo del confronto</a> distingue i risultati già ottenuti dalle nuove varianti ancora da implementare e validare.</p>

<h2>Le cariche RESP che abbiamo ottenuto</h2>
<p>Valori in e. Le colonne L0 sono fit RESP con cariche condivise per ruolo: HF indica il potenziale HF/6-31G*; DFT indica il controllo ωB97X-D/6-311++G(d,p). Tutti e due i fit usano le medesime 16 strutture.</p>
__KEY_CHARGES__
<p class="note">Le righe evidenziate riguardano COO⁻. Questi numeri sono candidati di ricerca; le librerie originali sono conservate. Il CSV contiene sei decimali, la tabella ne mostra quattro.</p>
<details><summary>Tutti i 25 ruoli, compresi EPH/EPT e bordo</summary>__ALL_CHARGES__</details>

<h2>Carbossilato: il fit molecolare aiuta a interpretare L0</h2>
<p>Nel fit HF–RESP atomico di pyrene_COO, il carbonio carbossilico è +0,807216 e e ciascun ossigeno −0,786679 e, vicini ai valori della libreria attuale. Il valore +1,731213 e deriva dal fit modulare condiviso con altri vincoli. Questo richiede un confronto controllato dei vincoli; non dimostra che HF–RESP imponga di assegnare +1,73 e a quel carbonio.</p>
<img class="plot" src="__COO_IMAGE__" alt="Cariche di ancoraggio, carbonio e ossigeno COO nei tre candidati: libreria attuale, RESP HF atomico e RESP HF L0.">
<p class="note">Nel fit atomico il sottoinsieme ancoraggio + COO somma a −0,680093 e, mentre la carica dell’intera molecola resta vincolata a −1 e. Nel L0 il blocco stesso deve sommare a −1 e. La conservazione del totale e la neutralità/carica imposta a ogni blocco sono vincoli distinti. Un solo contesto fisico non stabilisce una libreria trasferibile.</p>

<h2>RESP, ESP e L0: tre cose distinte</h2>
<p><strong>ESP</strong> è il potenziale elettrostatico quantistico. <strong>RESP</strong> è il fit con restrizioni che determina le cariche atomiche a partire da quel potenziale. L’errore ESP riportato nei grafici valuta le cariche RESP: non significa che usiamo le cariche ESP non ristrette di Gaussian. Questo è il rapporto descritto nella <a href="https://q4md-forcefieldtools.org/RED/resp/">documentazione RESP</a> e nel <a href="https://doi.org/10.1021/j100142a004">lavoro originale di Bayly e colleghi</a>.</p>
<p><strong>L0</strong> è il nostro modello RESP modulare di base: un valore per ruolo chimico, C interno non funzionalizzato a 0 e, blocchi neutri a 0 e e COO⁻ a −1 e. COOH e COO⁻ includono il carbonio di ancoraggio. I vincoli mantengono i blocchi assemblabili; la loro accuratezza va verificata.</p>
<details><summary>Protocollo effettivamente eseguito</summary>
<ol><li>Geometrie ottimizzate in fase gas a ωB97X-D/6-311++G(d,p), con frequenze e controlli elettronici registrati.</li>
<li>Tre orientazioni rigide di ogni geometria. Il riferimento convenzionale usa il potenziale MK HF/6-31G*; il ramo DFT–RESP usa un potenziale DFT separato.</li>
<li>Fit con il programma RESP di Amber: qwt = 0,0005, ihfree = 1, equivalenze dichiarate, carica totale intera e pesi delle orientazioni registrati. I blocchi L0 aggiungono i vincoli modulari.</li>
<li>Nessun gruppo CH₂/CH₃ da rifinire in questi frammenti: lo stadio 2 è dichiarato non operativo, senza simulare un secondo fit.</li></ol>
<p>Le tre rotazioni sono una geometria fisica, non tre conformeri. Il ramo HF mantiene la densità ristretta; non ne dichiara validata la stabilità elettronica. L’instabilità RHF→UHF dei precedenti probe è distinta da un errore SCF o dall’incompatibilità automatica con Amber. Gli eseguibili separati del fit sono verificati contro il pilot del RESP installato; Gaussian e Amber di sistema sono preservati.</p>
</details>

<h2>Quali gruppi sono coperti</h2>
<p>Idrossili basali (HYD), epossidi basali (EPX), epossidi al bordo (EPH/EPT), carbonili (CAB), COOH (CBX), COO⁻ (CBD), scheletro interno e bordo C–H (GGG/GRT). Ci sono coppie OH cis/trans, OH ed epossidi vicini/lontani e miscele. Il pannello non esaurisce tutte le distanze, densità o protonazioni.</p>
<p class="note">C₁₆/C₂₄/C₅₄ indicano lo scheletro aromatico iniziale; il carbonio del carbossile è aggiuntivo. I casi OH_EPX hanno due OH e un epossido. Lattoni e altre chimiche proposte non sono stati parametrizzati qui.</p>

<h2>La prima verifica di trasferibilità</h2>
<p>Un fit L0 usa solo i sette contesti con scheletro C₁₆. Sui sette C₂₄ con ruoli coperti, esclusi da quel fit, l’errore relativo mediano sul potenziale HF è <strong>0,1792</strong>, contro <strong>0,3501</strong> della libreria attuale. Il fit L0 che include anche i C₂₄ dà 0,1777 sugli stessi sette casi. EPH/EPT restano non disponibili nel fit C₁₆.</p>
<p>È un risultato favorevole all’idea dei blocchi RESP condivisi su questo pannello. Non basta per dimostrare che ogni C interno possa avere carica zero o che le stesse cariche funzionino su qualunque foglio.</p>
<details><summary>Confronto del potenziale per struttura: griglia HF o DFT</summary>
<div class="tabs"><button id="bh" aria-pressed="true" onclick="showGrid('HF')">Potenziale HF</button><button id="bd" aria-pressed="false" onclick="showGrid('DFT')">Potenziale DFT</button></div>
<img class="plot" id="hf" src="__HF_IMAGE__" alt="Errore ESP su sedici strutture per cariche RESP atomiche, RESP L0, libreria attuale e AM1-BCC; valori numerici e casi N/D.">
<img class="plot" id="dft" style="display:none" src="__DFT_IMAGE__" alt="Identico confronto sul potenziale DFT, sulla medesima scala.">
<p class="note">Media delle tre orientazioni di ogni struttura; N/D indica copertura assente, non zero. Scala comune HF/DFT. Le RESP atomiche e il fit L0 sui 16 casi sono valutati anche sul proprio training: questi residui non sono una classifica dei force field. AM1-BCC* include una correzione di arrotondamento dichiarata, distinta dai valori nativi conservati.</p></details>

<h2>Come procederei</h2>
<ol><li><strong>Terrei RESP HF/6-31G* come riferimento principale.</strong> DFT–RESP e AM1-BCC restano controlli di sensibilità. L0 è il candidato di base da valutare; non sostituirei ancora la libreria attuale.</li>
<li><strong>Verificherei prima la robustezza di COO⁻ e dei vincoli.</strong> Il C carbossilico passa da +0,7600 a +1,7312 e, compensato da O e ancoraggio. Non dimostra da solo un errore né una redistribuzione fisica corretta. COOH e COO⁻ hanno ciascuno un solo contesto fisico nel pannello attuale. Confronterei vincoli di blocco alternativi dichiarati, sensibilità del fit e contesti carbossilato aggiuntivi.</li>
<li><strong>Metterei alla prova il fondo a zero.</strong> Confronterei L0 con una variante che permetta la compensazione su un intorno limitato di carboni, mantenendo carica totale e regole di assemblaggio. La scelta dipenderebbe dalle prestazioni su interi contesti esclusi dal fit, non sulle sole rotazioni.</li>
<li><strong>Completerei C₅₄ e i contesti mirati prima di una regola ambientale.</strong> Priorità a singolo/coppia EPX su dimensioni confrontabili e a COO⁻ vicino/lontano da gruppi basali. Aggiungerei classi locali solo se il miglioramento fuori dall’addestramento fosse riproducibile. Non ricaverei adesso una legge universale dalla sola distanza.</li>
<li><strong>Sceglierei infine una coppia cariche–parametri.</strong> Dopo il controllo delle RESP, confronterei cariche e Lennard-Jones insieme su acqua e sonde di adsorbimento, con protonazione e condizioni dichiarate. I modelli pubblicati restano nella propria definizione; OPLS e GAFF2-AIM non diventano RESP cambiando il nome. Solo dopo congelerei le nuove LIB e rifarei i quattro materiali.</li></ol>
<p class="notice"><strong>Decisione attuale:</strong> nessun preset finale selezionato, nessun nuovo frcmod approvato e nessun modello originale rimpiazzato. I calcoli disponibili sostengono ulteriori verifiche del modello modulare, non una parametrizzazione già conclusa.</p>

<h2>Quali esperimenti possono aiutarci</h2>
<p>Le cariche parziali non sono numeri univocamente misurati per ciascun atomo. XPS e spettroscopie vincolano la composizione; titolazioni vincolano siti ionizzabili e protonazione del campione. <a href="https://pmc.ncbi.nlm.nih.gov/articles/PMC10339393/">Tao et al.</a> combinano titolazioni e modelli DFT per interpretare i gruppi del GO. <a href="https://doi.org/10.1016/j.jcis.2016.01.013">Orth et al.</a> applicano titolazioni e spettroscopia anche a materiali funzionalizzati.</p>
<p>La <a href="https://doi.org/10.1021/acs.jpcc.6b11571">calorimetria GO–amminoacidi di Pandit e De</a> è un possibile controllo successivo per l’interazione con piccole sonde; verifica il modello completo, con acqua e parametri coerenti. Non determina da sola una carica atomica né sostituisce i dati del campione LayerOne specifico.</p>

<h2>Riproducibilità e scala: cosa deve rimanere controllato</h2>
<p>Registrare ricetta, seed, connettività e IDs, protonazione, versione delle cariche e conteggi effettivi. La stessa realizzazione deve essere riproducibile; seed diversi consentono di verificare che il risultato non dipenda da un solo foglio casuale. L’assegnazione ordinaria dovrà usare una libreria o una regola calibrata su frammenti, senza richiedere un nuovo HF sull’intero foglio a ogni dimensione.</p>
<p>I siti basali crescono con l’area, quelli dei bordi esterni con il perimetro. Con carbossilati confinati ai bordi non è automatico conservare la stessa frazione COO/C o densità di carica aumentando il foglio. Densità, bordi dei difetti e limiti di capacità devono essere dichiarati. PEG, amminoacidi e peptidi richiederanno cariche coerenti anche per il nuovo giunto covalente, oltre a catena e terminali.</p>

<details><summary>AM1-BCC, OPLS, GOPY, OLD e GAFF2-AIM: copertura reale</summary>
<p>I quattro casi AM1-BCC esclusi sono pyrene_dione, C24_single_basal_EPX_control, C24_EPH_edge_control e C24_OH_OH_near_trans: manca un parametro nativo 16–12–2 o 12–16–2. Non viene imputato uno zero.</p>
<p>La tabella distingue la copertura delle cariche sui frammenti dalla validazione del modello completo. La proiezione AM1-BCC L0 è un candidato derivato, diverso dalle AM1-BCC native. I controlli Stauffer ibridi combinano fonti dichiarate. OLD di foglio intero e GAFF2-AIM non hanno qui un mapping completo: non ricevono uno score inventato.</p>
<div class="scroll"><table><thead><tr><th>Candidato</th><th>Contesti coperti</th></tr></thead><tbody>__COVERAGE__</tbody></table></div>
<p>La distribuzione originale GAFF2-AIM è mantenuta come riferimento separato: <a href="https://github.com/arvpinto/GO_sheets">repository degli autori</a>. Non è stata dimostrata qui la superiorità di RESP, OPLS o GAFF2-AIM per l’adsorbimento.</p>
</details>

<details><summary>Stato dei job, recuperi e limiti</summary>
<p>Sintesi 1452975 e recupero C₅₄ 1453077 completati, con artefatti verificati. Due C₅₄ DFT–RESP recuperati: parent pristine e coppia EPX vicina trans. Sono stati mantenuti tutti i punti ESP, circa 120.000–123.000 per orientazione, correggendo la lettura dell’header oltre cinque cifre in Python e in un eseguibile RESP separato. La regressione sul pilot coincide alla precisione di output; le equazioni del fit non sono cambiate.</p>
<p>Al controllo del 3 novembre il job 1406363_3, C54_single_basal_EPX_control, era ancora in esecuzione da circa 3 giorni e 11 ore, con limite di 5 giorni. Il postprocessing 1453111 attende la conclusione del padre. La dipendenza afterany non sostituisce i controlli sui risultati QM.</p>
<p>Le somme delle cariche dei gruppi nei riferimenti RESP atomici non sono vincolate come in L0. Non sono una misura univoca di trasferimento elettronico. I dipoli ricostruiti dalle cariche non sono errori rispetto a un dipolo QM indipendente. Mancano ancora verifica in solvente, adsorbimento e selezione della coppia cariche–LJ.</p>
</details>

<details><summary>Dati, riproducibilità e provenienza</summary>
<p>Questo HTML contiene immagini e stile: si apre offline e può essere condiviso da solo. I CSV restano nella cartella del report.</p>
<ul><li><a href="L0_shared_role_charges.csv">Cariche dei quattro fit RESP L0</a> e <a href="matched_HF_DFT_role_differences.csv">differenze RESP atomiche HF–DFT</a>.</li>
<li><a href="COO_charge_diagnostic.csv">Diagnostico COO</a>, <a href="COO_charge_diagnostic.svg">figura vettoriale COO</a> e <a href="basal_C_summary_by_context.csv">carboni basali per contesto e distanza</a>. Sono descrizioni dei fit esistenti, non nuove parametrizzazioni.</li>
<li><a href="RRMS_by_physical_context.csv">Residui per struttura</a>, <a href="CANDIDATE_COVERAGE.csv">copertura</a>, <a href="MATCHED_SUBSET_SCORES.csv">confronti sui medesimi sottoinsiemi</a> e <a href="SUMMARY_RESULTS.json">sintesi JSON</a>.</li>
<li><a href="context_HF_16/context_analysis.json">Contesti HF</a> e <a href="context_DFT_18/context_analysis.json">contesti DFT</a>, con cariche per atomo/gruppo e distanze nei CSV adiacenti.</li>
<li><a href="README.md">Dettaglio tecnico</a>, <a href="HF_ESP.svg">grafico HF vettoriale</a>, <a href="DFT_ESP.svg">grafico DFT vettoriale</a> e <a href="FIGURE_MANIFEST.json">manifest delle figure</a>.</li></ul>
<p>Archivi locali dei dati verificati:</p><pre>complementary_results_20261103_snapshot_001.tar.gz
SHA-256: 32ea75cec824cd06e740687215ec4d59267461f462ddca88e1509927c1d46552

c54_recovered_results_snapshot_20261103_attempt_001.tar.gz
SHA-256: 471798f22cd9d6ae5306d9b301c878a466cc11ac03ba32a3b1d04a786cca4759</pre>
<p class="note">Fonti: Gaussian e RESP reali, manifest e controlli di copertura; nessun nuovo Gaussian/SQM/Antechamber lanciato per questa sintesi. Root remoto: /export/scratch/alessandro.calza3/GO_PARAMETERIZATION_20261029_fqOFQv/calculations/. Le directory dei tentativi sono separate dagli originali.</p>
</details>
</main><script>
function showGrid(which){document.getElementById('hf').style.display=which==='HF'?'block':'none';document.getElementById('dft').style.display=which==='DFT'?'block':'none';document.getElementById('bh').setAttribute('aria-pressed',which==='HF');document.getElementById('bd').setAttribute('aria-pressed',which==='DFT')}
</script></body></html>'''
for token, replacement in {
    '__KEY_CHARGES__': charge_table(important), '__ALL_CHARGES__': charge_table(hf),
    '__HF_IMAGE__': image_uri('HF_ESP.png'), '__DFT_IMAGE__': image_uri('DFT_ESP.png'),
    '__COO_IMAGE__': image_uri('COO_charge_diagnostic.png'),
    '__COVERAGE__': coverage_rows,
}.items():
    page = page.replace(token, replacement)
assert '__KEY_' not in page and '__HF_IMAGE__' not in page
target = out/'RESP_AMBER_REPORT.html'
with target.open('w' if args.replace else 'x') as handle:
    handle.write(page)
manifest = {'report': target.name, 'self_contained_images_and_styles': True,
            'data_files_sha256': hashes, 'html_sha256': hashlib.sha256(page.encode()).hexdigest(),
            'final_parameter_validation': False}
(out/'HTML_REPORT_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(f'Generated {target} ({target.stat().st_size} bytes).')
