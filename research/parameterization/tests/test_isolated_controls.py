import json

from rdkit import Chem

from go_param import mixed_panel
from go_param.modular import partition_model


def test_isolated_controls_are_pristine_or_single_basal_epoxide_records():
    controls = mixed_panel.build_isolated_control_records()

    assert len(controls) == 4
    assert {record["metadata"]["parent_carbon_count"] for _, record in controls} == {
        24,
        54,
    }
    for mol, record in controls:
        metadata = record["metadata"]
        assert record["formal_charge"] == 0
        assert record["multiplicity"] == 1
        assert record["electron_count"] % 2 == 0
        assert sum(atom["formal_charge"] for atom in record["atoms"]) == 0
        assert all(atom["charge"] is None for atom in record["atoms"])
        assert (
            metadata["geometry_status"]
            == "initial_constrained_force_field_geometry_not_QM_optimized"
        )
        assert metadata["multiplicity_basis"].endswith(
            "QM wavefunction stability unverified"
        )
        partition_model(record)
        assert mol.GetNumAtoms() == len(record["atoms"])

    pristine = [
        record
        for _, record in controls
        if record["metadata"]["control_kind"] == "pristine_parent_reference"
    ]
    isolated_epx = [
        record
        for _, record in controls
        if record["metadata"]["control_kind"] == "single_basal_EPX"
    ]
    assert len(pristine) == 2
    assert len(isolated_epx) == 2
    assert all(record["metadata"]["functional_groups"] == [] for record in pristine)
    for record in isolated_epx:
        groups = record["metadata"]["functional_groups"]
        assert len(groups) == 1
        assert groups[0]["type"] == "EPX"
        assert record["metadata"]["pi_matching_complete"] is True


def test_single_epoxide_controls_have_no_pair_self_distance():
    controls = mixed_panel.build_isolated_control_records()

    for _, record in controls:
        metadata = record["metadata"]
        assert metadata["pairwise_distance_applicable"] is False
        assert metadata["measured_pair_motif_indices"] is None
        assert metadata.get("pair_graph_distance") is None
        assert metadata.get("pair_centroid_distance_angstrom") is None
        assert metadata.get("pair_oxygen_distance_angstrom") is None


def test_supplemental_export_writes_four_cataloged_sdf_xyz_and_models(tmp_path):
    catalog = mixed_panel.export_isolated_controls(tmp_path)

    assert len(catalog["models"]) == 4
    assert catalog["status"] == "initial_models_only_no_QM_results_or_fitted_charges"
    assert all(
        entry["pairwise_distance_applicable"] is False for entry in catalog["models"]
    )
    for entry in catalog["models"]:
        folder = tmp_path / entry["model_id"]
        record = json.loads((folder / "model.json").read_text())
        assert (folder / "initial.sdf").is_file()
        assert (folder / "initial.xyz").is_file()
        molecule = Chem.SDMolSupplier(str(folder / "initial.sdf"), removeHs=False)[0]
        assert molecule is not None
        assert molecule.GetNumAtoms() == len(record["atoms"])


def test_supplemental_export_does_not_modify_existing_panel_catalog(tmp_path):
    panel_dir = tmp_path / "generatedmodels_panel"
    panel_dir.mkdir()
    existing_catalog = panel_dir / "catalog.json"
    existing_catalog.write_text('{"models": ["preexisting"]}\n')
    before = existing_catalog.read_bytes()

    mixed_panel.export_isolated_controls(tmp_path / "generatedmodels_isolated_controls")

    assert existing_catalog.read_bytes() == before
