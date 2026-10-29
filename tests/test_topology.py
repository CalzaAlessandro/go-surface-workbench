from copy import deepcopy
from pathlib import Path

import numpy as np
import pytest

from go_builder.topology import export_gromacs, gbsa_readiness, import_amber, prepare_amber

REFERENCE = Path(__file__).resolve().parents[1] / "resources/reference/current"


def graph():
    return {
        "schema_version": 1,
        "name": "finite",
        "atoms": [
            {
                "id": 10,
                "element": "C",
                "name": "CH",
                "resname": "GRT",
                "resid": 8,
                "position": [0.0, 0.0, 0.0],
                "atom_type": None,
                "charge": None,
                "formal_charge": 0,
                "metadata": {},
            },
            {
                "id": 30,
                "element": "H",
                "name": "HC",
                "resname": "GRT",
                "resid": 8,
                "position": [1.08, 0.0, 0.0],
                "atom_type": None,
                "charge": None,
                "formal_charge": 0,
                "metadata": {},
            },
        ],
        "bonds": [{"a": 10, "b": 30, "order": 1, "image": [0, 0, 0]}],
        "cell": None,
        "pbc": [False] * 3,
        "metadata": {},
    }


def test_leap_uses_every_explicit_bond_no_distance_and_preserves_source(tmp_path):
    model = graph()
    original = deepcopy(model)
    result = prepare_amber(model, tmp_path / "amber", reference_dir=REFERENCE)
    leap = Path(result["leap_input"]).read_text()
    assert "bondbydistance" not in leap.lower()
    assert "bond mol.1.CH mol.1.HC" in leap
    assert "source leaprc.gaff2" in leap
    assert "set default PBRadii mbondi2" in leap
    assert result["status"] == "inputs_prepared"
    assert not result["parameterized"]
    assert model == original
    with pytest.raises(FileExistsError):
        prepare_amber(model, tmp_path / "amber", reference_dir=REFERENCE)


def test_nonzero_images_amber_gate_and_periodic_solvent_distinction(tmp_path):
    model = graph()
    model["pbc"] = [True, True, True]
    model["cell"] = [[30, 0, 0], [0, 30, 0], [0, 0, 30]]
    assert "periodic" not in " ".join(gbsa_readiness(model)["blockers"]).lower()
    prepare_amber(model, tmp_path / "finite_in_box", reference_dir=REFERENCE)
    model["bonds"][0]["image"] = [1, 0, 0]
    with pytest.raises(ValueError, match="periodic|image|winding"):
        prepare_amber(model, tmp_path / "winding", reference_dir=REFERENCE)
    assert any("periodic" in x.lower() for x in gbsa_readiness(model)["blockers"])


def test_gbsa_needs_parameterization_radii_and_partition():
    model = graph()
    report = gbsa_readiness(model)
    assert not report["ready"]
    assert any("charge" in x.lower() for x in report["blockers"])
    for atom in model["atoms"]:
        atom.update(charge=0, atom_type="ca", gb_radius=1.7, gb_screen=0.72)
    model["metadata"]["parameterized"] = True
    report = gbsa_readiness(model)
    assert report["ready"]
    assert report["masks"]["receptor"] == ":1"
    assert report["binding_ready"] is False


def parameterized_pair(tmp_path):
    import parmed as pmd

    s = pmd.Structure()
    c = pmd.Atom(name="C1", type="c3", atomic_number=6, mass=12.01, charge=-0.1)
    h = pmd.Atom(name="H1", type="h1", atomic_number=1, mass=1.008, charge=0.1)
    c.atom_type = pmd.topologyobjects.AtomType("c3", 1, 12.01, 6)
    c.atom_type.set_lj_params(0.1094, 1.908)
    h.atom_type = pmd.topologyobjects.AtomType("h1", 2, 1.008, 1)
    h.atom_type.set_lj_params(0.0157, 1.387)
    c.solvent_radius = 1.7
    c.screen = 0.72
    h.solvent_radius = 1.2
    h.screen = 0.85
    s.add_atom(c, "TST", 1)
    s.add_atom(h, "TST", 1)
    bt = pmd.BondType(340.0, 1.09)
    s.bond_types.append(bt)
    s.bonds.append(pmd.Bond(c, h, type=bt))
    s.coordinates = np.array([[0.0, 0.0, 0.0], [1.09, 0.0, 0.0]])
    # Close nonbonded coordinates are deliberately not bonded on import.
    extra = pmd.Atom(name="C2", type="c3", atomic_number=6, mass=12.01, charge=0)
    extra.atom_type = c.atom_type
    extra.solvent_radius = 1.7
    extra.screen = 0.72
    s.add_atom(extra, "TST", 2)
    s.coordinates = np.array([[0.0, 0.0, 0.0], [1.09, 0.0, 0.0], [0.2, 0.2, 0.0]])
    top = tmp_path / "input.prmtop"
    crd = tmp_path / "input.rst7"
    s.save(str(top), overwrite=False)
    s.save(str(crd), overwrite=False)
    return top, crd


def test_import_and_gromacs_roundtrip_preserve_parameterization(tmp_path):
    import parmed as pmd

    top, crd = parameterized_pair(tmp_path)
    model = import_amber(top, crd)
    assert len(model["atoms"]) == 3
    assert len(model["bonds"]) == 1
    assert model["atoms"][0]["charge"] == pytest.approx(-0.1)
    assert model["atoms"][0]["gb_radius"] == pytest.approx(1.7)
    assert model["metadata"]["parameterized"] is True
    out = export_gromacs(top, crd, tmp_path / "converted")
    assert out["parameterized"]
    assert Path(out["itp"]).exists()
    imported = pmd.load_file(out["top"], xyz=out["gro"])
    original = pmd.load_file(str(top), xyz=str(crd))
    assert len(imported.bonds) == 1
    assert np.array([a.charge for a in imported.atoms]) == pytest.approx(
        [a.charge for a in original.atoms]
    )
    assert imported.bonds[0].type.k == pytest.approx(original.bonds[0].type.k)
    assert imported.bonds[0].type.req == pytest.approx(original.bonds[0].type.req)
    assert imported.atoms[0].epsilon == pytest.approx(original.atoms[0].epsilon)
    assert imported.atoms[0].rmin == pytest.approx(original.atoms[0].rmin)
    assert imported.coordinates == pytest.approx(original.coordinates, abs=0.006)
    with pytest.raises(FileExistsError):
        export_gromacs(top, crd, tmp_path / "converted")


def test_leap_bond_order_uses_current_documented_string(tmp_path):
    model = graph()
    model["bonds"][0]["order"] = 1.5
    result = prepare_amber(model, tmp_path / "aromatic", reference_dir=REFERENCE)
    assert 'bond mol.1.CH mol.1.HC "A"' in Path(result["leap_input"]).read_text()


def test_leap_canonical_order_records_stable_ids_for_interleaved_residues(tmp_path):
    import json

    model = graph()
    c = deepcopy(model["atoms"][0])
    c.update(id=40, resid=20)
    h = deepcopy(model["atoms"][1])
    h.update(id=50, resid=20)
    model["atoms"] = [model["atoms"][0], c, model["atoms"][1], h]
    model["bonds"].append({"a": 40, "b": 50, "image": [0, 0, 0], "order": 1})
    result = prepare_amber(model, tmp_path / "interleaved", reference_dir=REFERENCE)
    assigned = json.loads(Path(result["assigned_model"]).read_text())
    assert [a["id"] for a in assigned["atoms"]] == [10, 30, 40, 50]
    assert result["atom_mapping"] == {"1": 10, "2": 30, "3": 40, "4": 50}


def test_itp_contains_reusable_molecule_definition(tmp_path):
    top, crd = parameterized_pair(tmp_path)
    out = export_gromacs(top, crd, tmp_path / "molecule")
    assert "[ moleculetype ]" in Path(out["itp"]).read_text()
    assert '#include "molecule.itp"' in Path(out["top"]).read_text()


def test_gromacs_preserves_angles_torsions_and_14_scaling(tmp_path):
    import parmed as pmd

    structure = pmd.Structure()
    atype = pmd.topologyobjects.AtomType("c3", 1, 12.01, 6)
    atype.set_lj_params(0.1094, 1.908)
    for index in range(4):
        atom = pmd.Atom(
            name=f"C{index}",
            type="c3",
            atomic_number=6,
            mass=12.01,
            charge=[-0.1, 0.05, 0.1, -0.05][index],
        )
        atom.atom_type = atype
        atom.solvent_radius = 1.7
        atom.screen = 0.72
        structure.add_atom(atom, "NUM", 1)
    bt = pmd.BondType(310.0, 1.526)
    structure.bond_types.append(bt)
    for a, b in ((0, 1), (1, 2), (2, 3)):
        structure.bonds.append(pmd.Bond(structure.atoms[a], structure.atoms[b], type=bt))
    angle = pmd.AngleType(40.0, 109.5)
    structure.angle_types.append(angle)
    for a, b, c in ((0, 1, 2), (1, 2, 3)):
        structure.angles.append(pmd.Angle(*(structure.atoms[x] for x in (a, b, c)), type=angle))
    for k, period, phase in ((0.155, 3, 0.0), (0.25, 1, 180.0)):
        term = pmd.DihedralType(k, period, phase, scee=1.2, scnb=2.0)
        structure.dihedral_types.append(term)
        structure.dihedrals.append(pmd.Dihedral(*structure.atoms, type=term))
    structure.coordinates = np.array([[0, 0, 0], [1.526, 0, 0], [2.0, 1.45, 0], [3.0, 1.8, 1.0]])
    top = tmp_path / "torsion.prmtop"
    crd = tmp_path / "torsion.rst7"
    structure.save(str(top))
    structure.save(str(crd))
    out = export_gromacs(top, crd, tmp_path / "torsion_gmx")
    amber = pmd.load_file(str(top), xyz=str(crd))
    gmx = pmd.load_file(out["top"], xyz=out["gro"])
    assert [(a.type.k, a.type.theteq) for a in gmx.angles] == pytest.approx(
        np.asarray([(a.type.k, a.type.theteq) for a in amber.angles])
    )

    def terms(s):
        result = []
        for dihedral in s.dihedrals:
            values = dihedral.type if isinstance(dihedral.type, list) else [dihedral.type]
            result.extend((t.phi_k, t.per, t.phase, t.scee, t.scnb) for t in values)
        return sorted(result)

    assert terms(gmx) == pytest.approx(np.asarray(terms(amber)))
    assert out["roundtrip_verified"]["torsion_terms"] == 2


def test_actual_mixed_sheet_fragment_assignment_and_stable_amber_order(tmp_path):
    import json

    from go_builder.build import build_sheet

    model = build_sheet(
        {
            "nx": 10,
            "ny": 8,
            "seed": 17,
            "groups": {"OH": 4, "epoxide": 3, "carbonyl": 2, "COOH": 2, "COO": 2},
        }
    )
    result = prepare_amber(model, tmp_path / "mixed", reference_dir=REFERENCE)
    assigned = json.loads(Path(result["assigned_model"]).read_text())
    assert len(assigned["atoms"]) == len(model["atoms"])
    assert set(result["atom_mapping"].values()) == {a["id"] for a in model["atoms"]}
    assert result["net_charge"] == pytest.approx(-2.0)
    by_residue = {}
    for atom in assigned["atoms"]:
        by_residue.setdefault(atom["resid"], []).append(atom)
    assert [a["id"] for atoms in by_residue.values() for a in atoms] == [
        a["id"] for a in assigned["atoms"]
    ]
    assert not result["parameterized"]


def test_gbsa_compresses_masks_and_rejects_covalent_binding_partition():
    model = graph()
    for atom in model["atoms"]:
        atom.update(charge=0.0, atom_type="ca", gb_radius=1.7, gb_screen=0.72)
    extra = deepcopy(model["atoms"][0])
    extra.update(id=50, resid=9)
    ligand = deepcopy(extra)
    ligand.update(id=70, resid=10, resname="LIG")
    model["atoms"].extend([extra, ligand])
    model["metadata"].update(parameterized=True, ligand_resids=[10])
    report = gbsa_readiness(model)
    assert report["masks"]["receptor"] == ":1-2"
    assert report["masks"]["ligand"] == ":3"
    assert report["binding_ready"]
    model["bonds"].append({"a": 50, "b": 70, "order": 1, "image": [0, 0, 0]})
    assert not gbsa_readiness(model)["binding_ready"]


def test_external_missing_parameter_error_cannot_be_reported_as_parameterized(tmp_path):
    import json

    tleap = tmp_path / "tleap_missing_terms"
    tleap.write_text('#!/bin/sh\nprintf "Could not find bond parameter for: cx - n3\\n"\nexit 0\n')
    tleap.chmod(0o755)
    out = tmp_path / "failed"
    with pytest.raises(ValueError, match="parameterization failed"):
        prepare_amber(graph(), out, reference_dir=REFERENCE, run=True, tleap=str(tleap))
    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert not manifest["parameterized"]
    assert manifest["topology"] is None
    assert "cx - n3" in (out / "tleap.log").read_text()


def test_preparing_new_inputs_invalidates_prior_parameterized_claim(tmp_path):
    model = graph()
    model["metadata"]["parameterized"] = True
    for atom in model["atoms"]:
        atom.update(charge=0.0, atom_type="ca", gb_radius=1.7, gb_screen=0.72)
    out = prepare_amber(model, tmp_path / "new_assignment", reference_dir=REFERENCE)
    assert not out["parameterized"]
    assert not out["gbsa"]["ready"]


@pytest.mark.parametrize("name", ["gox", "rgo"])
def test_supplied_old_amber_topologies_roundtrip_without_refitting(name, tmp_path):
    import parmed as pmd

    old = REFERENCE.parent / "old"
    top = old / f"{name}.parm7"
    crd = old / f"{name}.rst7"
    if not top.exists() or not crd.exists():
        pytest.skip("Historical topology snapshot not supplied")
    original = pmd.load_file(str(top), xyz=str(crd))
    model = import_amber(top, crd)
    assert len(model["atoms"]) == len(original.atoms)
    assert len(model["bonds"]) == len(original.bonds)
    assert [a["charge"] for a in model["atoms"]] == pytest.approx(
        [a.charge for a in original.atoms]
    )
    out = export_gromacs(top, crd, tmp_path / name)
    assert out["roundtrip_verified"]["atoms"] == len(original.atoms)
    assert out["roundtrip_verified"]["torsion_terms"] >= len(original.dihedrals)
