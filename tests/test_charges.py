import json
import shutil
from pathlib import Path

import pytest

from go_builder.charges import (
    audit_charge_sources,
    classify_connectivity,
    read_hiergo_charges,
    read_off_atoms,
)

ROOT = Path(__file__).resolve().parents[1]


def test_off_reads_atoms_despite_unrelated_malformed_boundbox():
    rows = read_off_atoms(ROOT / "resources/reference/current/CBD_amber.lib")
    assert len(rows) == 4
    assert sum(a["charge"] for a in rows) == pytest.approx(-0.5548)
    assert rows[1]["neighbors"] == ["C", "O", "O"]


def test_off_stray_epx_metadata_does_not_duplicate_eph_atoms():
    rows = read_off_atoms(ROOT / "resources/reference/current/EPH_amber.lib")
    assert len(rows) == 3
    assert {a["resname"] for a in rows} == {"EPH"}
    assert rows[2]["charge"] == pytest.approx(-0.3691)


def test_hiergo_same_type_can_have_different_contextual_charges():
    rows = read_hiergo_charges(ROOT / "resources/reference/opls/Oplsatom2")
    shared = [r for r in rows if r["atom_type"] == "opls_180"]
    assert {(r["resname"], r["charge"]) for r in shared} == {("EPX", -0.4), ("ETH", -0.17)}


def test_audit_preserves_source_context_and_refuses_unmatched_ranking(tmp_path):
    out = tmp_path / "audit"
    result = audit_charge_sources(out)
    assert result["ranking"] is None
    assert len([key for key in result["datasets"] if key.startswith("CURRENT:")]) == 9
    assert "BUNDLED_LEGACY_AMBER:EPX_amber.lib" in result["datasets"]
    assert "CURRENT:EPX_amber.lib" not in result["datasets"]
    assert result["current_charge_method"]["qm_method"] is None
    assert result["datasets"]["OLD:gox.parm7"]["total_charge_e"] != pytest.approx(
        result["datasets"]["OLD:rgo.parm7"]["total_charge_e"]
    )
    assert result["datasets"]["HIERGO:Oplsatom2"]["scope"] == "contextual_parameter_table"
    assert (out / "atoms.csv").is_file()
    assert (out / "motifs.csv").is_file()
    assert len(json.loads((out / "audit.json").read_text())["provenance"]) >= 20
    with pytest.raises(FileExistsError):
        audit_charge_sources(out)


def test_connectivity_motifs_do_not_confuse_epoxide_and_ether():
    atoms = [
        {"id": 1, "element": "C"},
        {"id": 2, "element": "C"},
        {"id": 3, "element": "O"},
        {"id": 4, "element": "C"},
        {"id": 5, "element": "O"},
    ]
    motifs = classify_connectivity(atoms, [(1, 2), (1, 3), (2, 3), (2, 5), (4, 5)])
    assert motifs[3] == "epoxide"
    assert motifs[5] == "ether"
    assert motifs[2] == "epoxide+ether"


def test_audit_detects_modified_official_hiergo_snapshot(tmp_path):
    copied = tmp_path / "opls"
    shutil.copytree(ROOT / "resources/reference/opls", copied)
    table = copied / "Oplsatom2"
    table.write_text(table.read_text().replace("-0.4000", "-0.9000"))
    with pytest.raises(ValueError, match="hash"):
        audit_charge_sources(tmp_path / "bad_audit", opls_dir=copied)
    assert not (tmp_path / "bad_audit").exists()
