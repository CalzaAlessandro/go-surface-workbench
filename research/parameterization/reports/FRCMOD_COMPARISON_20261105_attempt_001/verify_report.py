#!/usr/bin/env python3
"""Independent section-count, link, source-hash and deterministic rebuild checks."""
import csv
import hashlib
import json
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[1]
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROOT/"src"))
from go_param.frcmod_audit import parse_frcmod


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


files=sorted(OUT.glob("numeric_*.csv"))+[OUT/"numeric_report.html",OUT/"numeric_summary.json"]
before={p.name:digest(p) for p in files}
subprocess.run([sys.executable,"-B",str(OUT/"build_report.py")],check=True,capture_output=True,text=True)
after={p.name:digest(p) for p in files}
assert before==after,"Rebuild is not byte deterministic"
summary=json.loads((OUT/"numeric_summary.json").read_text())
for source in summary["sources"]:
    old=parse_frcmod(ROOT/source["file"])
    for section in ("MASS","BOND","ANGLE","DIHE","IMPROPER","NONBON"):
        if source["source"]=="Hoshin_original_GAFF" and section=="NONBON":
            assert not old["sections"][section] and source[section]==5
        else:
            assert len(old["sections"][section])==source[section],(source["source"],section)
    if source["source"]!="Hoshin_original_GAFF":
        assert len(old["issues"])==source["invalid_or_incomplete_rows"]
    else:
        assert len(old["issues"])==6 # NONB plus 5 LJ lines incorrectly retained in IMPROPER by old parser.


class Page(HTMLParser):
    def __init__(self):
        super().__init__();self.links=[];self.tables=0
    def handle_starttag(self,tag,attrs):
        if tag=="a":self.links.extend(v for k,v in attrs if k=="href")
        if tag=="table":self.tables+=1


page=Page();page.feed((OUT/"numeric_report.html").read_text())
assert all((OUT/link).is_file() for link in page.links)
assert len(list(csv.DictReader((OUT/"numeric_terms.csv").open())))==1205
assert len(summary["invalid_rows"])==1 and summary["invalid_rows"][0]["line"]==139
assert digest(OUT/"build_report.py")==summary["script_sha256"]
result={
    "deterministic_rebuild":True,
    "all_generator_source_hash_and_Fourier_assertions_passed":True,
    "existing_production_parser_section_counts_agree_except_documented_NONB_fix":True,
    "production_parser_modified":False,
    "Hoshin_NONB_alias_recovery_rows":5,
    "valid_numeric_parameter_rows":1205,
    "incomplete_rows_preserved":1,
    "html_tables":page.tables,
    "html_valid_local_artifact_links":len(page.links),
    "HTML_QA_scope":"structure and local link validation; no browser render",
    "generator_sha256":digest(OUT/"build_report.py"),
    "verifier_sha256":digest(Path(__file__)),
    "outputs_sha256":after,
}
(OUT/"numeric_verification.json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(result,ensure_ascii=False,indent=2))
