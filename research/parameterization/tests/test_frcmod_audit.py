"""Parameter text audit tests: distinguish Fourier continuation from conflicts."""

import hashlib
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def test_multiterm_dihedrals_and_zero_ho_lj_are_not_invalid(tmp_path):
    from go_param.frcmod_audit import audit_file

    path = tmp_path / "sample.frcmod"
    path.write_text("""Synthetic audit fixture
MASS
c3 12.01
ho 1.008
BOND
c3-c3 300.0 1.5
ANGLE
c3-c3-c3 60.0 109.5
DIHE
X -c3-c3-X 1 0.3 180.0 -3.0
X -c3-c3-X 1 0.0 0.0 2.0
IMPROPER
X -X -c3-ho 1.1 180.0 2.0
NONBON
c3 1.908 0.1094
ho 0.0 0.0
""")
    result = audit_file(path)
    assert not [i for i in result["issues"] if i["severity"] == "error"]
    assert result["dihedral_series"][0]["component_count"] == 2
    assert result["sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert any(
        i["code"] == "contextual_zero_hydroxyl_hydrogen_lj" for i in result["issues"]
    )
    assert not result["duplicate_groups"]


def test_reversed_bond_conflicts_and_repeated_full_dihedral_definitions(tmp_path):
    from go_param.frcmod_audit import audit_file

    path = tmp_path / "conflict.frcmod"
    path.write_text("""Fixture
BOND
ca-cx 300.0 1.5
cx-ca 310.0 1.5
DIHE
X -ca-cx-X 1 1.0 0.0 2.0
X -ca-cx-X 1 2.0 0.0 2.0
""")
    result = audit_file(path)
    assert len(result["duplicate_groups"]) == 2
    assert all(g["conflicting"] for g in result["duplicate_groups"])
    assert any(
        i["code"] == "conflicting_definition" and i["section"] == "BOND"
        for i in result["issues"]
    )


def test_nonfinite_negative_lj_and_placeholder_are_explicit(tmp_path):
    from go_param.frcmod_audit import audit_file

    path = tmp_path / "bad.frcmod"
    path.write_text("""Fixture
MASS
ca nan
ANGLE
ca-ca-ca 0.0 120.0 ATTN needs revision
DIHE
X -ca-ca-X 0 1.0 0.0 -2.0
NONBON
ca -1.9 0.1
""")
    result = audit_file(path)
    codes = {i["code"] for i in result["issues"]}
    assert {
        "nonfinite_numeric_value",
        "negative_lj_parameter",
        "placeholder_comment",
        "invalid_dihedral_divisor",
        "unterminated_dihedral_continuation",
    } <= codes


def test_actual_current_epoxide_angles_and_library_namespace_coverage():
    from go_param.frcmod_audit import audit_file, audit_library_types

    path = ROOT / "references/local_snapshot/current/GO_gaff2.frcmod"
    result = audit_file(path)
    entries = result["sections"]["ANGLE"]
    row = next(e for e in entries if e["types"] == ["cx", "op", "cx"])
    assert row["values"] == [119.01, 61.76]
    libraries = audit_library_types(
        ROOT / "references/local_snapshot/current", result, suffix="_gaff2.lib"
    )
    assert any(
        i["file"].endswith("EPH_gaff2.lib") and i["code"] == "mixed_type_namespace"
        for i in libraries["issues"]
    )
    assert "CA" in libraries["override_coverage"]["types_without_mass_override"]
    assert (
        libraries["override_coverage"]["base_force_field_completeness"]
        == "not_assessed"
    )


def test_report_generation_is_read_only_for_sources(tmp_path):
    from go_param.frcmod_audit import generate_report

    source = ROOT / "references/local_snapshot/current/GO_gaff2.frcmod"
    original = source.read_bytes()
    result = generate_report(ROOT, output_dir=tmp_path)
    assert source.read_bytes() == original
    assert len(result["files"]) >= 5
    assert (tmp_path / "frcmod_audit.json").exists()
    assert (tmp_path / "frcmod_audit_issues.csv").exists()
    assert result["scientific_status"] == "text_audit_not_force_field_validation"
    assert json.loads((tmp_path / "frcmod_audit.json").read_text())["base_comparison"][
        "status"
    ] == (
        "available_versioned_direct_type_comparison"
        if (ROOT / "references/system_forcefield_sources/gaff2.dat").is_file()
        else "not_available"
    )


@pytest.mark.skipif(
    not all(
        (ROOT / "references/system_forcefield_sources" / name).is_file()
        for name in ("gaff.dat", "gaff2.dat", "parm99.dat")
    ),
    reason="External installed Amber DAT files are intentionally not vendored",
)
def test_installed_base_comparison_and_orphan_index_are_explicit(tmp_path):
    from go_param.frcmod_audit import generate_report

    report = generate_report(ROOT, output_dir=tmp_path)
    assert (
        report["base_comparison"]["status"]
        == "available_versioned_direct_type_comparison"
    )
    current = report["base_comparison"]["comparisons"]["GO_gaff2.frcmod"]
    row = next(
        r
        for r in current["terms"]
        if r["section"] == "BOND" and r["types"] == ["ca", "ca"]
    )
    assert row["base_candidates"][0]["values"] == [354.25, 1.3986]


def test_orphan_index_and_user_provenance_without_installed_data(tmp_path):
    from go_param.frcmod_audit import generate_report

    report = generate_report(ROOT, output_dir=tmp_path)
    lib = next(
        l
        for l in report["library_audits"]["_gaff2.lib"]["libraries"]
        if l["file"].endswith("EPH_gaff2.lib")
    )
    assert lib["indexed_units"] == ["EPH"]
    assert lib["orphan_unit_blocks"] == ["EPX"]
    assert "user_reported" in report["historical_parameter_provenance"]["status"]
