"""Graph identity, image conventions and lossless exchange regressions."""

import copy
import json

import pytest

from go_builder.model import load_model, save_model, validate_model, write_pdb


def pair():
    return {
        "schema_version": 1,
        "name": "pair",
        "atoms": [
            {"id": i, "element": "C", "position": p, "formal_charge": 0}
            for i, p in [(4, [0, 0, 0]), (19, [1.42, 0, 0])]
        ],
        "bonds": [{"a": 4, "b": 19, "order": 1.5, "image": [0, 0, 0]}],
        "cell": [[10, 0, 0], [0, 10, 0], [0, 0, 40]],
        "pbc": [False] * 3,
        "metadata": {},
    }


def test_json_preserves_sparse_ids_and_explicit_bonds(tmp_path):
    m = pair()
    save_model(m, tmp_path / "m.json")
    assert load_model(tmp_path / "m.json") == m
    with pytest.raises(FileExistsError):
        save_model(m, tmp_path / "m.json")
    write_pdb(m, tmp_path / "m.pdb", allow_undercoordinated=True)
    pdb_text = (tmp_path / "m.pdb").read_text()
    assert "CONECT    4   19" in pdb_text
    assert "CONECT   19    4" in pdb_text


@pytest.mark.parametrize("mutation", ["duplicate", "dangling", "image", "nan", "charge"])
def test_invalid_graph_is_rejected(mutation):
    m = pair()
    if mutation == "duplicate":
        m["atoms"].append(copy.deepcopy(m["atoms"][0]))
    if mutation == "dangling":
        m["bonds"][0]["b"] = 42
    if mutation == "image":
        m["bonds"][0]["image"] = [1, 0, 0]
    if mutation == "nan":
        m["atoms"][0]["position"][0] = float("nan")
    if mutation == "charge":
        m["atoms"][0]["charge"] = float("inf")
    with pytest.raises(ValueError):
        validate_model(m)


def test_validation_does_not_change_or_renumber_graph():
    m = pair()
    before = json.dumps(m)
    report = validate_model(m)
    assert report["valid"] and report["atom_count"] == 2
    assert json.dumps(m) == before


def test_structurally_undercoordinated_import_is_flagged_for_review():
    # The exchange schema admits fragments, but never silently calls two
    # singly bonded aromatic carbons a saturated chemical structure.
    assert validate_model(pair())["undercoordinated_atom_ids"] == [4, 19]


def test_fused_aromatic_carbon_is_not_mistaken_for_hypervalent_carbon():
    from go_builder.build import build_sheet

    m = build_sheet({"nx": 2, "ny": 2})
    for atom in m["atoms"]:
        atom["metadata"] = {}
    report = validate_model(m)
    assert report["valid"]
    assert report["undercoordinated_atom_ids"] == []


def test_geometry_integrity_gate_catches_clashes_bad_bonds_and_undercoordination(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import (
        assert_export_ready,
        check_geometry_integrity,
        write_mol2,
        write_xyz,
    )

    m = build_sheet({"nx": 2, "ny": 2, "seed": 1})
    rep = check_geometry_integrity(m)
    assert rep["passed"] is True
    assert rep["clash_count"] == 0
    assert rep["undercoordinated_atom_ids"] == []
    assert rep["formal_charge"] == 0

    # Synthetic H clash < 0.90 A
    h_clash = copy.deepcopy(m)
    h_atom = next(a for a in h_clash["atoms"] if a["element"] == "H")
    other_c = next(a for a in h_clash["atoms"] if a["element"] == "C" and a["id"] != h_atom["metadata"]["anchor_id"])
    h_atom["position"] = [other_c["position"][0] + 0.85, other_c["position"][1], other_c["position"][2]]
    h_rep = check_geometry_integrity(h_clash)
    assert h_rep["passed"] is False
    assert h_rep["clash_count"] >= 1
    assert h_rep["min_nonbonded_h_distance"] < 0.90
    for fn in (assert_export_ready, lambda x: write_xyz(x, tmp_path / "c.xyz"), lambda x: write_mol2(x, tmp_path / "c.mol2"), lambda x: write_pdb(x, tmp_path / "c.pdb")):
        with pytest.raises(ValueError, match="clash"):
            fn(h_clash)

    # Synthetic anomalous bond length > 1.85 A
    bad_bond = copy.deepcopy(m)
    bad_bond["atoms"][0]["position"][2] += 2.0
    b_rep = check_geometry_integrity(bad_bond)
    assert b_rep["passed"] is False
    assert b_rep["bond_lengths_valid"] is False
    with pytest.raises(ValueError, match="bond length"):
        assert_export_ready(bad_bond)

    # Undercoordinated fragment fails write_pdb by default (passes only with allow_undercoordinated=True) and fails strict XYZ/MOL2
    frag = pair()
    assert check_geometry_integrity(frag)["passed"] is False
    with pytest.raises(ValueError, match="undercoordinated"):
        write_pdb(frag, tmp_path / "frag_strict.pdb")
    write_pdb(frag, tmp_path / "frag_ok.pdb", allow_undercoordinated=True)
    with pytest.raises(ValueError, match="undercoordinated"):
        write_xyz(frag, tmp_path / "frag.xyz")
    with pytest.raises(ValueError, match="undercoordinated"):
        write_mol2(frag, tmp_path / "frag.mol2")


def test_write_xyz_and_write_mol2_formats_and_periodic_rejection(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import write_cif, write_extxyz, write_mol2, write_sdf, write_xyz

    m = build_sheet(
        {
            "nx": 4,
            "ny": 4,
            "seed": 17,
            "name": "mixed_go",
            "groups": {"OH": 1, "epoxide": 1, "carbonyl": 1, "COOH": 1, "COO": 1},
        }
    )
    xyz_path = tmp_path / "go.xyz"
    mol2_path = tmp_path / "go.mol2"
    write_xyz(m, xyz_path)
    write_mol2(m, mol2_path)

    xyz_lines = xyz_path.read_text().splitlines()
    assert int(xyz_lines[0]) == len(m["atoms"])
    assert "name=mixed_go" in xyz_lines[1]
    assert "formal_charge=-1" in xyz_lines[1]
    assert f"electron_count={m['metadata']['electron_count']}" in xyz_lines[1]
    assert "geometry_status=constructed_unrelaxed" in xyz_lines[1]
    assert len(xyz_lines) == len(m["atoms"]) + 2

    mol2_text = mol2_path.read_text()
    assert "@<TRIPOS>MOLECULE" in mol2_text
    assert "@<TRIPOS>ATOM" in mol2_text
    assert "@<TRIPOS>BOND" in mol2_text
    assert "NO_CHARGES" in mol2_text
    assert "0.0000" in mol2_text
    assert "C.ar" in mol2_text
    assert "C.3" in mol2_text
    assert "C.2" in mol2_text
    assert "O.co2" in mol2_text
    assert "O.2" in mol2_text
    assert "O.3" in mol2_text
    assert " ar\n" in mol2_text
    assert " 2\n" in mol2_text

    # Periodic sheet: standard XYZ/MOL2/SDF rejected, while EXTXYZ and CIF succeed
    pbc_m = build_sheet({"nx": 3, "ny": 3, "pbc": [True, True, False], "groups": {"OH": 1}})
    with pytest.raises(ValueError, match="periodic image bonds"):
        write_xyz(pbc_m, tmp_path / "pbc.xyz")
    with pytest.raises(ValueError, match="periodic image bonds"):
        write_mol2(pbc_m, tmp_path / "pbc.mol2")
    with pytest.raises(ValueError, match="periodic image bonds"):
        write_sdf(pbc_m, tmp_path / "pbc.sdf")
    write_extxyz(pbc_m, tmp_path / "pbc.extxyz")
    write_cif(pbc_m, tmp_path / "pbc.cif")
    assert 'pbc="T T F"' in (tmp_path / "pbc.extxyz").read_text()
    assert "_symmetry_space_group_name_H-M 'P 1'" in (tmp_path / "pbc.cif").read_text()


def test_electron_parity_diagnostics_and_eph_ept_mol2_tripos_types(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import check_geometry_integrity, write_mol2

    odd_sheet = build_sheet({"nx": 2, "ny": 2, "seed": 1, "groups": {"OH": 1}})
    rep_odd = check_geometry_integrity(odd_sheet)
    assert rep_odd["passed"] is True
    assert rep_odd["electron_parity"] == "odd"
    assert rep_odd["electron_count"] % 2 == 1
    assert any("open-shell" in w or "Odd electron count" in w for w in rep_odd["warnings"])

    even_sheet = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {"OH": 2, "epoxide": 1, "epoxide_eph": 2, "epoxide_ept": 1},
        }
    )
    rep_even = check_geometry_integrity(even_sheet)
    assert rep_even["passed"] is True
    assert rep_even["electron_parity"] == "even"
    assert rep_even["electron_count"] % 2 == 0
    assert rep_even["warnings"] == []

    mol2_path = tmp_path / "p2.mol2"
    write_mol2(even_sheet, mol2_path)
    atom_block = mol2_path.read_text().split("@<TRIPOS>ATOM\n")[1].split("@<TRIPOS>BOND")[0]
    tripos_by_res_name = {}
    for line in atom_block.strip().splitlines():
        parts = line.split()
        # id name x y z type resid resname charge
        tripos_by_res_name[(parts[7], parts[1])] = parts[5]

    for name, expected_type in [("CZ", "C.3"), ("CK", "C.3"), ("OE", "O.3"), ("HK", "H")]:
        assert tripos_by_res_name[("EPH", name)] == expected_type
    for name, expected_type in [("CP", "C.3"), ("CK", "C.3"), ("OE", "O.3"), ("HP", "H"), ("HK", "H")]:
        assert tripos_by_res_name[("EPT", name)] == expected_type


def test_14plus_steric_gate_and_motif_bond_length_checks():
    from go_builder.build import build_sheet
    from go_builder.model import check_geometry_integrity

    # Two 1-4+ hydrogens placed at 1.10 A (< 1.20 A threshold) must fail with kind="14plus_h"
    m = build_sheet({"nx": 3, "ny": 3, "seed": 5})
    h_atoms = [a for a in m["atoms"] if a["element"] == "H"]
    # Pick two H atoms bonded to different carbons (1-4+ graph distance >= 3)
    h1, h2 = h_atoms[0], h_atoms[-1]
    c_nbrs = {
        b["a"] if b["b"] == h1["id"] else b["b"]
        for b in m["bonds"]
        if h1["id"] in (b["a"], b["b"])
    }
    assert h2["id"] not in c_nbrs
    h_clash = copy.deepcopy(m)
    ha = next(a for a in h_clash["atoms"] if a["id"] == h1["id"])
    hb = next(a for a in h_clash["atoms"] if a["id"] == h2["id"])
    # Move hb to be 1.10 A from ha along +z so it doesn't distort hb's own C-H? Wait: keep hb bonded length valid
    # by placing ha and hb across a synthetic 4-atom chain!
    synth = {
        "schema_version": 1,
        "name": "synth_14",
        "atoms": [
            {"id": 1, "element": "C", "position": [0.0, 0.0, 0.0], "formal_charge": 0, "metadata": {"hybridization": "sp3"}},
            {"id": 2, "element": "C", "position": [1.50, 0.0, 0.0], "formal_charge": 0, "metadata": {"hybridization": "sp3"}},
            {"id": 3, "element": "H", "position": [0.0, 1.09, 0.0], "formal_charge": 0, "metadata": {}},
            {"id": 4, "element": "H", "position": [1.10, 1.09, 0.0], "formal_charge": 0, "metadata": {}},
        ],
        "bonds": [
            {"a": 1, "b": 2, "order": 1.0, "image": [0, 0, 0]},
            {"a": 1, "b": 3, "order": 1.0, "image": [0, 0, 0]},
            {"a": 2, "b": 4, "order": 1.0, "image": [0, 0, 0]},
        ],
        "cell": [[20.0, 0.0, 0.0], [0.0, 20.0, 0.0], [0.0, 0.0, 40.0]],
        "pbc": [False, False, False],
        "metadata": {},
    }
    # Check C-H bond lengths: 1-3 is 1.09 A, 2-4 is sqrt(0.4^2 + 1.09^2) = 1.161 A (in [1.00, 1.18] A)
    rep14 = check_geometry_integrity(synth)
    assert rep14["passed"] is False
    assert rep14["clash_count"] >= 1
    c_14 = [c for c in rep14["clashes"] if c["kind"] == "14plus_h"]
    assert len(c_14) == 1
    assert {c_14[0]["a"], c_14[0]["b"]} == {3, 4}
    assert c_14[0]["distance"] == pytest.approx(1.10, abs=1e-4)
    assert c_14[0]["threshold"] == 1.20

    # Motif-specific C-H bond length validation: 0.894 A and 1.256 A must fail even if inside [0.85, 1.85]
    for bad_ch in (0.894, 1.256):
        bad = copy.deepcopy(m)
        ha = next(a for a in bad["atoms"] if a["id"] == h1["id"])
        anc = next(a for a in bad["atoms"] if a["id"] == ha["metadata"]["anchor_id"])
        dx = ha["position"][0] - anc["position"][0]
        dy = ha["position"][1] - anc["position"][1]
        dz = ha["position"][2] - anc["position"][2]
        cur = (dx * dx + dy * dy + dz * dz) ** 0.5
        scale = bad_ch / cur
        ha["position"] = [
            anc["position"][0] + dx * scale,
            anc["position"][1] + dy * scale,
            anc["position"][2] + dz * scale,
        ]
        rep_bad = check_geometry_integrity(bad)
        assert rep_bad["passed"] is False
        assert rep_bad["bond_lengths_valid"] is False
        assert any(
            tuple(sorted(ab["elements"])) == ("C", "H")
            for ab in rep_bad["bond_length_violations"]
        )


def test_compute_chemistry_summary_reports_domains_environments_and_qm():
    from go_builder.build import build_sheet
    from go_builder.model import compute_chemistry_summary

    m = build_sheet(
        {
            "nx": 6,
            "ny": 6,
            "seed": 17,
            "groups": {"OH": 4, "epoxide": 2, "carbonyl": 2, "COOH": 1, "COO": 1},
        }
    )
    summary = compute_chemistry_summary(m)
    assert summary["formula"].startswith("C")
    assert summary["atom_counts"]["C"] == sum(1 for a in m["atoms"] if a["element"] == "C")
    assert summary["atom_counts"]["O"] == sum(1 for a in m["atoms"] if a["element"] == "O")
    assert summary["atom_counts"]["H"] == sum(1 for a in m["atoms"] if a["element"] == "H")
    assert sum(summary["atomic_fractions"].values()) == pytest.approx(1.0)
    assert summary["formal_charge"] == -1
    assert summary["c_over_o"] == pytest.approx(
        summary["atom_counts"]["C"] / summary["atom_counts"]["O"]
    )
    assert len(summary["dimensions_nm"]["cell_nm"]) == 3
    assert len(summary["dimensions_nm"]["extent_nm"]) == 3
    assert summary["sp2_carbons"] + summary["sp3_carbons"] == summary["scaffold_carbons"]
    assert summary["sp2_domains"]["count"] >= 1
    assert summary["sp2_domains"]["largest_domain_carbons"] <= summary["sp2_carbons"]

    envs = summary["chemical_environments"]
    assert envs["basis"] == "structural_graph_carbon_environments_not_xps_intensities"
    assert envs["counts"]["C=O_carbonyl"] == 2
    assert envs["counts"]["O-C=O_carboxyl"] == 2
    assert envs["carboxyl_breakdown"] == {"carboxyl_C4": 2, "anchor_CY": 2}
    assert sum(envs["counts"].values()) == summary["atom_counts"]["C"]
    assert envs["counts"]["C-O_sp3"] == 4 + 2 * 2  # 4 OH + 2 epoxides (2 C each)

    qm = summary["qm_diagnostics"]
    assert qm["electron_count"] == m["metadata"]["electron_count"]
    assert qm["electron_parity"] in ("even", "odd")
    sub = qm["sp2_sublattice_counts"]
    assert sub is not None and sub["A"] + sub["B"] <= summary["sp2_carbons"]
    assert qm["sp2_sublattice_imbalance"] == abs(sub["A"] - sub["B"])
    assert "min_spin_multiplicity_hint" not in qm


def test_write_pdb_reciprocal_conect_formal_charge_and_bounds(tmp_path):
    from go_builder.build import build_sheet

    m = build_sheet({"nx": 4, "ny": 4, "seed": 0, "groups": {"COO": 1}})
    pdb_path = tmp_path / "coo.pdb"
    write_pdb(m, pdb_path)
    lines = pdb_path.read_text().splitlines()

    hetatm_lines = [ln for ln in lines if ln.startswith("HETATM")]
    assert len(hetatm_lines) == len(m["atoms"])
    # Check wwPDB columns 77-78 (element) and 79-80 (formal charge): CBD anionic oxygen gets '1-'
    charged_lines = [ln for ln in hetatm_lines if "CBD" in ln and ln[78:80] == "1-"]
    assert len(charged_lines) == 1
    assert charged_lines[0][76:78].strip() == "O"
    neutral_lines = [ln for ln in hetatm_lines if ln[78:80] == "  "]
    assert len(neutral_lines) == len(m["atoms"]) - 1

    # Verify reciprocal CONECT records (every bond a-b appears under CONECT a and CONECT b, max 4 nbrs/line)
    conect_adj: dict[int, set[int]] = {}
    for ln in lines:
        if not ln.startswith("CONECT"):
            continue
        fields = [int(ln[6 + 5 * k : 11 + 5 * k]) for k in range((len(ln) - 6) // 5) if ln[6 + 5 * k : 11 + 5 * k].strip()]
        assert 2 <= len(fields) <= 5
        src = fields[0]
        conect_adj.setdefault(src, set()).update(fields[1:])
    for bond in m["bonds"]:
        a_id, b_id = bond["a"], bond["b"]
        assert b_id in conect_adj[a_id]
        assert a_id in conect_adj[b_id]

    # Strict wwPDB coordinate bounds: -1000.000 must fail (takes 9 chars in %8.3f)
    oob = pair()
    oob["atoms"][0]["position"] = [-1000.000, 0.0, 0.0]
    oob["atoms"][1]["position"] = [-998.580, 0.0, 0.0]
    with pytest.raises(ValueError, match="coordinate range exceeded"):
        write_pdb(oob, tmp_path / "oob.pdb", allow_undercoordinated=True)

    # -999.999 is accepted
    ok = pair()
    ok["atoms"][0]["position"] = [-999.999, 0.0, 0.0]
    ok["atoms"][1]["position"] = [-998.579, 0.0, 0.0]
    write_pdb(ok, tmp_path / "ok_bound.pdb", allow_undercoordinated=True)


def test_write_mol2_contiguous_residues_parmed_and_rdkit_roundtrip(tmp_path):
    import parmed as pmd
    from rdkit import Chem

    from go_builder.build import build_sheet
    from go_builder.model import write_mol2

    m = build_sheet({"nx": 4, "ny": 4, "seed": 0, "groups": {"COO": 1}})
    mol2_path = tmp_path / "go_coo.mol2"
    write_mol2(m, mol2_path)

    # ParmEd with structure=True must read exact atom, bond, and contiguous residue counts
    struct = pmd.load_file(str(mol2_path), structure=True)
    assert len(struct.atoms) == 76
    assert len(struct.bonds) == 93
    assert len(struct.residues) == 54

    # RDKit must load all 76 atoms and preserve -1 formal charge
    mol = Chem.MolFromMol2File(str(mol2_path), sanitize=False, removeHs=False)
    assert mol is not None
    assert mol.GetNumAtoms() == 76
    assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == -1


def test_write_xyz_recomputes_electron_count_live(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import electronic_summary, write_xyz

    m = build_sheet({"nx": 3, "ny": 3, "seed": 2, "groups": {"OH": 1}})
    expected_elec = electronic_summary(m)["electron_count"]
    m["metadata"]["electron_count"] = 999999
    xyz_path = tmp_path / "live.xyz"
    write_xyz(m, xyz_path)
    comment = xyz_path.read_text().splitlines()[1]
    assert f"electron_count={expected_elec}" in comment
    assert "999999" not in comment


def test_write_extxyz_periodic_and_finite(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import write_extxyz

    pbc_m = build_sheet({"nx": 3, "ny": 3, "pbc": [True, True, False], "name": "pbc_sheet", "groups": {"OH": 1}})
    pbc_path = tmp_path / "pbc.extxyz"
    write_extxyz(pbc_m, pbc_path)
    pbc_lines = pbc_path.read_text().splitlines()
    assert int(pbc_lines[0]) == len(pbc_m["atoms"])
    assert len(pbc_lines) == len(pbc_m["atoms"]) + 2
    assert 'Lattice="' in pbc_lines[1]
    assert "Properties=species:S:1:pos:R:3" in pbc_lines[1]
    assert 'pbc="T T F"' in pbc_lines[1]
    assert 'name="pbc_sheet"' in pbc_lines[1]
    assert "formal_charge=0" in pbc_lines[1]
    assert "geometry_status=" in pbc_lines[1]

    fin_m = build_sheet({"nx": 3, "ny": 3, "name": "fin_sheet"})
    fin_m["cell"] = None
    fin_path = tmp_path / "fin.extxyz"
    write_extxyz(fin_m, fin_path)
    fin_header = fin_path.read_text().splitlines()[1]
    assert "Lattice=" not in fin_header
    assert 'pbc="F F F"' in fin_header


def test_write_sdf_v3000_rdkit_roundtrip_preserves_cbd_charge(tmp_path):
    from rdkit import Chem

    from go_builder.build import build_sheet
    from go_builder.model import write_sdf

    m = build_sheet({"nx": 4, "ny": 4, "seed": 0, "groups": {"COO": 1}})
    sdf_path = tmp_path / "go_coo.sdf"
    write_sdf(m, sdf_path)
    sdf_text = sdf_path.read_text()
    assert "V3000" in sdf_text
    assert "M  V30 BEGIN CTAB" in sdf_text
    assert "CHG=-1" in sdf_text
    assert ">  <formal_charge>" in sdf_text
    assert ">  <electron_count>" in sdf_text
    assert ">  <geometry_status>" in sdf_text

    suppl = Chem.SDMolSupplier(str(sdf_path), sanitize=False, removeHs=False)
    mol = suppl[0]
    assert mol is not None
    assert mol.GetNumAtoms() == 76
    assert mol.GetNumBonds() == 93
    assert sum(a.GetFormalCharge() for a in mol.GetAtoms()) == -1
    conf = mol.GetConformer()
    for idx, atom in enumerate(m["atoms"]):
        pos = conf.GetAtomPosition(idx)
        assert [pos.x, pos.y, pos.z] == pytest.approx(atom["position"], abs=1e-4)


def test_write_cif_p1_periodic_and_finite(tmp_path):
    from go_builder.build import build_sheet
    from go_builder.model import write_cif

    pbc_m = build_sheet({"nx": 3, "ny": 3, "pbc": [True, True, False], "name": "cif pbc"})
    pbc_cif = tmp_path / "pbc.cif"
    write_cif(pbc_m, pbc_cif)
    text_pbc = pbc_cif.read_text()
    assert "data_cif_pbc" in text_pbc
    assert "_cell_length_a" in text_pbc
    assert "_symmetry_space_group_name_H-M 'P 1'" in text_pbc
    site_lines = [
        ln.split()
        for ln in text_pbc.splitlines()
        if ln and ln[0] in ("C", "O", "H") and not ln.startswith("_")
    ]
    assert len(site_lines) == len(pbc_m["atoms"])
    # Periodic x, y fractional coordinates must be wrapped in [0, 1)
    for cols in site_lines:
        fx, fy = float(cols[2]), float(cols[3])
        assert 0.0 <= fx < 1.0
        assert 0.0 <= fy < 1.0

    fin_m = build_sheet({"nx": 3, "ny": 3, "name": "cif_fin"})
    fin_m["cell"] = None
    fin_cif = tmp_path / "fin.cif"
    write_cif(fin_m, fin_cif)
    text_fin = fin_cif.read_text()
    assert "_cell_length_a" in text_fin
    assert "_atom_site_Cartn_z" in text_fin
    xs = [a["position"][0] for a in fin_m["atoms"]]
    a_line = next(ln for ln in text_fin.splitlines() if ln.startswith("_cell_length_a"))
    assert float(a_line.split()[1]) == pytest.approx((max(xs) - min(xs)) + 20.0, abs=1e-4)
