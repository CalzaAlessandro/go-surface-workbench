import json

import pytest

from go_builder.alignment import patch_inputs, rewrite_alignment


def test_missing_fit_mask_uses_explicit_go_anchor_and_records_previous_all_atom_fit():
    result = rewrite_alignment("autoimage :1-494\nalign first \ntrajout out.nc\n")
    assert result["text"] == "autoimage :1-494\nrms first :1-494\ntrajout out.nc\n"
    assert result["changes"][0]["previous_fit"] == "all_atoms_default"
    assert result["changes"][0]["relative_geometry_at_risk"] is False


def test_masked_align_is_replaced_and_marked_as_relative_geometry_risk():
    result = rewrite_alignment("align first :1-628\n")
    assert result["text"] == "rms first :1-628\n"
    assert result["changes"][0]["relative_geometry_at_risk"] is True


def test_move_all_mass_comment_and_template_placeholder_are_preserved_correctly():
    result = rewrite_alignment("autoimage MASK\nalign first MASK mass move * # whole system\n")
    assert result["text"] == "autoimage MASK\nrms first MASK mass # whole system\n"
    assert result["changes"][0]["relative_geometry_at_risk"] is False


def test_no_mask_is_never_guessed_without_anchor():
    with pytest.raises(ValueError, match="mask"):
        rewrite_alignment("align first\n")


def test_patch_is_audited_backed_up_and_idempotent(tmp_path):
    source = tmp_path / "align_10ns.in"
    old = "align first :1-494\n"
    source.write_text(old)
    plan = patch_inputs([source], tmp_path / "plan", apply=False)
    assert source.read_text() == old
    assert not (tmp_path / "plan").exists()
    assert plan["changed_files"] == 1
    report = patch_inputs([source], tmp_path / "backup", apply=True)
    assert source.read_text() == "rms first :1-494\n"
    assert json.loads((tmp_path / "backup/manifest.json").read_text())["status"] == "APPLIED"
    backup = tmp_path / "backup" / report["files"][0]["backup"]
    assert backup.read_text() == old
    assert patch_inputs([source], tmp_path / "second", apply=False)["changed_files"] == 0
