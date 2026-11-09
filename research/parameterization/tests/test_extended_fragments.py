"""Graph gates catch motif substitutions; no fitted charges are fabricated."""

import importlib
import json

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors


def api():
    assert importlib.util.find_spec("go_param.extended_fragments") is not None, "extended fragment API missing"
    return importlib.import_module("go_param.extended_fragments")


@pytest.mark.parametrize("name,formula", [
    ("dibenzofuran_ether_edge", "C12H8O"),
    ("diphenyl_ether_bridge_control", "C12H10O"),
    ("phthalide_lactone5_edge", "C8H6O2"),
    ("coumarin_lactone6_edge", "C9H6O2"),
    ("oxolan2ol_lactol_control", "C4H8O2"),
    ("phenol_aromatic_OH_control", "C6H6O"),
    ("naphthoquinone14_edge", "C10H6O2"),
])
def test_explicit_reference_graphs_are_neutral_closed_shell_and_resp_serializable(name, formula):
    from go_param.resp import _validate_model

    mol, record = api().build_candidate(name, seed=20261105)
    Chem.SanitizeMol(mol)
    assert rdMolDescriptors.CalcMolFormula(mol) == record["formula"] == formula
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 0
    assert Chem.GetFormalCharge(mol) == record["formal_charge"] == 0
    assert record["electron_count"] == sum(a.GetAtomicNum() for a in mol.GetAtoms())
    assert record["electron_count"] % 2 == 0 and record["multiplicity"] == 1
    assert all(a.GetNumImplicitHs() == 0 for a in mol.GetAtoms())
    assert sum(a.GetSymbol() == "H" for a in mol.GetAtoms()) > 0
    assert [a.GetAtomMapNum() for a in mol.GetAtoms()] == list(range(1, mol.GetNumAtoms()+1))
    assert all(a["role"].startswith("extended_") for a in record["atoms"])
    assert all("charge_e" not in a for a in record["atoms"])
    assert record["metadata"]["parameter_status"] == "UNPARAMETERIZED_RESEARCH_REFERENCE"
    assert record["metadata"]["geometry_status"] == "MMFF94s_starting_geometry_not_QM_accepted"
    assert record["metadata"]["wavefunction_stability_validated"] is False
    _validate_model(record)


@pytest.mark.parametrize("smiles,motif,ring_size", [
    ("C1CO1", "non_epoxide_ether", None),  # epoxide must not be relabelled ETH
    ("COC(C)=O", "lactone", 5),  # acyclic ester is not a lactone
    ("O=C1OCc2ccccc21", "lactone", 6),  # wrong ester ring size
    ("O=C1OCCC1", "lactol", None),  # cyclic ester is not a hemiacetal
    ("C1CCOC1", "lactol", None),  # cyclic ether lacks anomeric OH
    ("CCO", "phenolic_OH", None),  # aliphatic OH is not phenol
    ("O=C1CCCCC1", "quinone", None),  # one ketone is not a quinone
    ("O=C1CCC(=O)CC1", "quinone", None),  # two saturated ketones lack conjugation
])
def test_wrong_real_graphs_fail_motif_gate(smiles, motif, ring_size):
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    with pytest.raises(ValueError, match="motif|ring|epoxide"):
        api().validate_motif(mol, motif, ring_size=ring_size)


def test_ether_lactone_and_lactol_have_distinct_connectivity():
    module = api()
    ether, _ = module.build_candidate("dibenzofuran_ether_edge")
    assert not any(len(r) == 3 for r in ether.GetRingInfo().AtomRings())
    for name, size in [("phthalide_lactone5_edge", 5), ("coumarin_lactone6_edge", 6)]:
        mol, record = module.build_candidate(name)
        carbonyl_c = {b.GetOtherAtom(o).GetIdx() for o in mol.GetAtoms() if o.GetSymbol() == "O" for b in o.GetBonds() if b.GetBondTypeAsDouble() == 2}
        assert len(carbonyl_c) == 1
        assert any(len(r) == size and carbonyl_c <= set(r) and any(mol.GetAtomWithIdx(i).GetSymbol() == "O" for i in r) for r in mol.GetRingInfo().AtomRings())
        assert record["metadata"]["motif_graph_gate"]["ring_size"] == size
    lactol, record = module.build_candidate("oxolan2ol_lactol_control")
    assert not any(b.GetBondTypeAsDouble() == 2 for b in lactol.GetBonds())
    centers = [a for a in lactol.GetAtoms() if a.GetSymbol() == "C" and sum(n.GetSymbol() == "O" for n in a.GetNeighbors()) == 2]
    assert len(centers) == 1 and centers[0].GetHybridization() == Chem.HybridizationType.SP3
    assert any(n.GetSymbol() == "O" and any(h.GetSymbol() == "H" for h in n.GetNeighbors()) for n in centers[0].GetNeighbors())
    assert record["metadata"]["reference_scope"] == "aliphatic_connectivity_control"


def test_seed_reproduces_geometry_and_unknown_candidate_is_rejected():
    module = api()
    _, first = module.build_candidate("diphenyl_ether_bridge_control", seed=71)
    _, second = module.build_candidate("diphenyl_ether_bridge_control", seed=71)
    np.testing.assert_allclose([a["position"] for a in first["atoms"]], [a["position"] for a in second["atoms"]], atol=1e-12)
    with pytest.raises(ValueError, match="Unknown"):
        module.build_candidate("made_up_fragment")


def test_lactone_carbonyl_and_ring_oxygens_have_distinct_unparameterized_roles():
    module = api()
    _, record = module.build_candidate("phthalide_lactone5_edge")
    roles = [a["role"] for a in record["atoms"] if a["element"] == "O"]
    assert set(roles) == {"extended_lactone_carbonyl_O", "extended_lactone_ring_O"}


def test_coumarin_json_and_sanitized_sdf_keep_the_same_aromatic_graph(tmp_path):
    api().export_candidates(tmp_path / "panel")
    folder = tmp_path / "panel" / "coumarin_lactone6_edge"
    record = json.loads((folder / "model.json").read_text())
    mol = Chem.SDMolSupplier(str(folder / "initial.sdf"), removeHs=False)[0]
    assert mol is not None
    expected = {tuple(sorted((b["a"], b["b"]))): b["order"] for b in record["bonds"]}
    observed = {tuple(sorted((b.GetBeginAtomIdx()+1, b.GetEndAtomIdx()+1))):
                b.GetBondTypeAsDouble() for b in mol.GetBonds()}
    assert observed == expected
    for row, atom in zip(record["atoms"], mol.GetAtoms()):
        if row["role"] == "extended_background_C":
            assert not atom.GetIsAromatic()


def test_export_preserves_mapping_prepares_only_qm_and_refuses_overwrite(tmp_path):
    module = api()
    catalog = module.export_candidates(tmp_path / "panel", seed=20261105, cpus=2, mem_gb=4)
    assert len(catalog["models"]) == 7
    assert catalog["qm_submitted"] is False and catalog["fitted_partial_charges_available"] is False
    for entry in catalog["models"]:
        folder = tmp_path / "panel" / entry["model_id"]
        record = json.loads((folder / "model.json").read_text())
        sdf = next(m for m in Chem.SDMolSupplier(str(folder / "initial.sdf"), removeHs=False) if m is not None)
        assert sdf.GetNumAtoms() == len(record["atoms"])
        assert [a.GetAtomMapNum() for a in sdf.GetAtoms()] == [a["id"] for a in record["atoms"]]
        assert int((folder / "initial.xyz").read_text().splitlines()[0]) == len(record["atoms"])
        dft = (folder / "geometry_opt" / "reference.gjf").read_text()
        hf = (folder / "HF_RESP_template" / "reference.gjf").read_text()
        assert "wB97XD/6-311++G(d,p)" in dft and "Freq" in dft and "Stable" in dft
        assert "HF/6-31G*" in hf and "Pop=MK" in hf and "Opt" not in hf
        assert record["metadata"]["parameter_status"] == "UNPARAMETERIZED_RESEARCH_REFERENCE"
    with pytest.raises(FileExistsError):
        module.export_candidates(tmp_path / "panel")
