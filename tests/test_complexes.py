import json
import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import parmed as pmd
import pytest
from parmed.tools.changeradii import ChRad

from go_builder.complexes import assemble_complex


def methane(tmp_path, name, offset=(0.0, 0.0, 0.0), epsilon=0.1094, radii="mbondi"):
    """Complete numerical methane fixture with independently assigned terms."""
    structure = pmd.Structure()
    carb = pmd.topologyobjects.AtomType("c3", 1, 12.01, 6)
    carb.set_lj_params(epsilon, 1.908)
    hydr = pmd.topologyobjects.AtomType("hc", 2, 1.008, 1)
    hydr.set_lj_params(0.0157, 1.487)
    for i in range(5):
        atom = pmd.Atom(
            name="C1" if i == 0 else f"H{i}",
            type="c3" if i == 0 else "hc",
            atomic_number=6 if i == 0 else 1,
            mass=12.01 if i == 0 else 1.008,
            charge=-0.4 if i == 0 else 0.1,
        )
        atom.atom_type = carb if i == 0 else hydr
        structure.add_atom(atom, "MET", 1)
    bt = pmd.BondType(340.0, 1.09)
    structure.bond_types.append(bt)
    for i in range(1, 5):
        structure.bonds.append(pmd.Bond(structure.atoms[0], structure.atoms[i], type=bt))
    at = pmd.AngleType(33.0, 109.4712)
    structure.angle_types.append(at)
    for i, j in combinations(range(1, 5), 2):
        structure.angles.append(
            pmd.Angle(structure.atoms[i], structure.atoms[0], structure.atoms[j], type=at)
        )
    structure.coordinates = (
        np.asarray(
            [
                [0, 0, 0],
                [0.629, 0.629, 0.629],
                [0.629, -0.629, -0.629],
                [-0.629, 0.629, -0.629],
                [-0.629, -0.629, 0.629],
            ]
        )
        + offset
    )
    parm = pmd.amber.AmberParm.from_structure(structure)
    ChRad(parm, radii)
    top = tmp_path / f"{name}.parm7"
    crd = tmp_path / f"{name}.rst7"
    parm.save(str(top))
    parm.save(str(crd))
    return top, crd


def test_single_ligand_assembly_preserves_parameters_and_adds_no_cross_bonds(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    before = {p: p.read_bytes() for p in (rec, rc, lig, lc)}
    result = assemble_complex(rec, rc, [lig], [lc], tmp_path / "complex")
    assert result["parameterized"]
    assert result["masks"] == {
        "receptor": ":1",
        "ligand": ":2",
        "complex": ":1-2",
        "ligands": [":2"],
    }
    assert result["gbsa"]["binding_ready"]
    complex_ = pmd.load_file(result["complex_topology"], xyz=result["complex_coordinates"])
    assert len(complex_.atoms) == 10
    assert len(complex_.bonds) == 8
    assert len(complex_.angles) == 12
    assert all((b.atom1.idx < 5) == (b.atom2.idx < 5) for b in complex_.bonds)
    assert [a.charge for a in complex_.atoms] == pytest.approx([-0.4, 0.1, 0.1, 0.1, 0.1] * 2)
    assert [a.solvent_radius for a in complex_.atoms] == pytest.approx(
        [1.7] + [1.3] * 4 + [1.7] + [1.3] * 4
    )
    assert [b.type.k for b in complex_.bonds] == pytest.approx([340.0] * 8)
    assert all(p.read_bytes() == data for p, data in before.items())
    assert Path(result["complex_pdb"]).read_text().count("CONECT") == 8
    assert Path(result["receptor_pdb"]).exists()
    assert (tmp_path / "complex/lig.parm7").exists()
    with pytest.raises(FileExistsError):
        assemble_complex(rec, rc, [lig], [lc], tmp_path / "complex")


def test_multiple_simultaneous_ligands_keep_independent_masks(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    l1, c1 = methane(tmp_path, "l1", (0, 0, 5))
    l2, c2 = methane(tmp_path, "l2", (0, 0, -5))
    out = assemble_complex(rec, rc, [l1, l2], [c1, c2], tmp_path / "multi")
    assert out["masks"]["ligands"] == [":2", ":3"]
    assert out["masks"]["ligand"] == ":2-3"
    assert len(pmd.load_file(out["ligand_topology"]).atoms) == 10
    assert len(out["ligands"]) == 2
    for entry in out["ligands"]:
        assert len(pmd.load_file(entry["topology"]).atoms) == 5


def test_pose_must_be_complete_and_formats_explicit(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    with pytest.raises(ValueError, match="length|one|count"):
        assemble_complex(rec, rc, [lig], [], tmp_path / "mismatch")
    unknown = tmp_path / "pose.pdbqt"
    unknown.write_text("not a supported full pose")
    with pytest.raises(ValueError, match="format|PDBQT|pose"):
        assemble_complex(rec, rc, [lig], [unknown], tmp_path / "unknown")
    short = pmd.amber.Rst7(natom=1)
    short.coordinates = np.zeros((1, 3))
    short.write(str(tmp_path / "partial.rst7"))
    with pytest.raises(ValueError, match="atom|count"):
        assemble_complex(rec, rc, [lig], [tmp_path / "partial.rst7"], tmp_path / "short")


def test_strict_pdb_pose_atom_identity_order_element_and_partial_detection(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    s = pmd.load_file(str(lig), xyz=str(lc))
    pose = tmp_path / "pose.pdb"
    s.save(str(pose))
    out = assemble_complex(rec, rc, [lig], [pose], tmp_path / "valid_pdb")
    assert out["ligands"][0]["coordinate_format"] == "pdb_checked"
    lines = pose.read_text().splitlines()
    index = next(i for i, line in enumerate(lines) if line.startswith(("ATOM  ", "HETATM")))
    original = lines[index]
    lines[index] = original[:12] + "XXXX" + original[16:]
    bad = tmp_path / "wrong_name.pdb"
    bad.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="identity|name|order"):
        assemble_complex(rec, rc, [lig], [bad], tmp_path / "bad_name")
    lines[index] = original[:76] + " N" + original[78:]
    bad2 = tmp_path / "wrong_element.pdb"
    bad2.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="element|identity"):
        assemble_complex(rec, rc, [lig], [bad2], tmp_path / "bad_element")


@pytest.mark.parametrize("case", ["radii", "lj", "overlap"])
def test_incompatible_inputs_are_rejected_without_output(tmp_path, case):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(
        tmp_path,
        "lig",
        (0, 0, 0 if case == "overlap" else 5),
        epsilon=0.12 if case == "lj" else 0.1094,
        radii="mbondi2" if case == "radii" else "mbondi",
    )
    out = tmp_path / "rejected"
    with pytest.raises(ValueError, match="radii|Lennard|overlap|contact"):
        assemble_complex(rec, rc, [lig], [lc], out)
    assert not out.exists()


def test_old_receptor_with_independently_parameterized_ligand_preserves_all_terms(tmp_path):
    old = Path(__file__).resolve().parents[1] / "resources/reference/old"
    top = old / "gox.parm7"
    crd = old / "gox.rst7"
    if not top.exists():
        pytest.skip("Historical receptor snapshot not supplied")
    receptor = pmd.load_file(str(top), xyz=str(crd))
    offset = np.max(receptor.coordinates, axis=0) + [5, 5, 5]
    lig, lc = methane(tmp_path, "lig", offset)
    out = assemble_complex(top, crd, [lig], [lc], tmp_path / "old_complex")
    combined = pmd.load_file(out["complex_topology"], xyz=out["complex_coordinates"])
    assert len(combined.atoms) == 945
    assert len(combined.dihedrals) == len(receptor.dihedrals)
    assert len(combined.bonds) == len(receptor.bonds) + 4
    assert np.asarray([a.charge for a in combined.atoms[:940]]) == pytest.approx(
        [a.charge for a in receptor.atoms]
    )
    assert np.asarray(combined.coordinates[:940]) == pytest.approx(receptor.coordinates, abs=1e-6)
    assert out["masks"]["receptor"] == ":1"
    assert out["masks"]["ligand"] == ":2"
    assert out["roundtrip_verified"]["complex"]["torsion_terms"] == len(receptor.dihedrals)


def test_cli_assembles_complete_pose_and_reports_dynamic_masks(tmp_path):
    rec, rc = methane(tmp_path, "cli_rec")
    lig, lc = methane(tmp_path, "cli_lig", (0, 0, 5))
    root = Path(__file__).resolve().parents[1]
    output = tmp_path / "cli_complex"
    run = subprocess.run(
        [
            sys.executable,
            str(root / "06_assemble_docked_complex.py"),
            "--receptor-topology",
            str(rec),
            "--receptor-coordinates",
            str(rc),
            "--ligand-topology",
            str(lig),
            "--ligand-coordinates",
            str(lc),
            "--output-dir",
            str(output),
            "--json",
        ],
        text=True,
        capture_output=True,
    )
    assert run.returncode == 0, run.stdout + run.stderr
    report = json.loads(run.stdout)
    assert report["masks"]["receptor"] == ":1"
    assert report["masks"]["ligand"] == ":2"
    assert report["cross_component_bonds"] == 0
    assert (output / "compl.parm7").is_file()


def pdbqt_pose(structure, path, order=None, second_model=False):
    order = list(range(len(structure.atoms))) if order is None else order
    records = []
    for serial, index in enumerate(order, 1):
        atom = structure.atoms[index]
        x, y, z = structure.coordinates[index]
        ad = {1: "H", 6: "C", 7: "N", 8: "OA", 9: "F", 17: "Cl", 35: "Br", 53: "I"}[
            atom.atomic_number
        ]
        records.append(
            f"ATOM  {serial:5d} {atom.name:<4s} MET A   1    {x:8.3f}{y:8.3f}{z:8.3f}  0.00  0.00    0.123 {ad:>2s}\n"
        )
    path.write_text(
        "MODEL 1\nROOT\n"
        + "".join(records)
        + "ENDROOT\nTORSDOF 0\nENDMDL\n"
        + ("MODEL 2\nATOM invalid second pose ignored\nENDMDL\n" if second_model else "")
    )


def sdf_pose(structure, path, order=None, names=False, invert=False):
    from rdkit import Chem

    rw = Chem.RWMol()
    for atom in structure.atoms:
        rd = Chem.Atom(atom.atomic_number)
        rd.SetNoImplicit(True)
        rw.AddAtom(rd)
    for bond in structure.bonds:
        rw.AddBond(bond.atom1.idx, bond.atom2.idx, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    conf = Chem.Conformer(len(structure.atoms))
    for index, xyz in enumerate(structure.coordinates):
        xyz = list(map(float, xyz))
        xyz[0] = -xyz[0] if invert else xyz[0]
        conf.SetAtomPosition(index, xyz)
    mol.AddConformer(conf)
    if order is not None:
        mol = Chem.RenumberAtoms(mol, order)
    if names:
        mol.SetProp(
            "GO_BUILDER_ATOM_NAMES",
            json.dumps(
                [
                    structure.atoms[index].name
                    for index in (order or list(range(len(structure.atoms))))
                ]
            ),
        )
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()
    return mol


def test_complete_pdbqt_first_model_maps_names_after_tree_permutation(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    template = pmd.load_file(str(lig), xyz=str(lc))
    pose = tmp_path / "vina.pdbqt"
    pdbqt_pose(template, pose, [3, 0, 4, 2, 1], second_model=True)
    out = assemble_complex(rec, rc, [lig], [pose], tmp_path / "pdbqt")
    coordinates = pmd.load_file(out["ligand_topology"], xyz=out["ligand_coordinates"]).coordinates
    assert coordinates == pytest.approx(template.coordinates, abs=0.0006)
    assert out["ligands"][0]["coordinate_format"] == "pdbqt_first_model_named"


@pytest.mark.parametrize("problem", ["missing_h", "duplicate_name", "wrong_element"])
def test_pdbqt_rejects_partial_or_ambiguous_atom_identity(tmp_path, problem):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    template = pmd.load_file(str(lig), xyz=str(lc))
    pose = tmp_path / "bad.pdbqt"
    if problem == "missing_h":
        pdbqt_pose(template, pose, [0, 1, 2, 3])
    else:
        pdbqt_pose(template, pose)
        lines = pose.read_text().splitlines()
        index = next(i for i, line in enumerate(lines) if line.startswith("ATOM"))
        if problem == "duplicate_name":
            lines[index] = lines[index][:12] + "H1  " + lines[index][16:]
        else:
            lines[index] = lines[index][:-2] + " N"
        pose.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="PDBQT|identity|hydrogen|atom"):
        assemble_complex(rec, rc, [lig], [pose], tmp_path / "rejected")
    assert not (tmp_path / "rejected").exists()


def test_sdf_graph_mapping_complete_hydrogens_and_equivalent_symmetry(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    template = pmd.load_file(str(lig), xyz=str(lc))
    pose = tmp_path / "gnina.sdf"
    sdf_pose(template, pose, [3, 0, 4, 2, 1])
    out = assemble_complex(rec, rc, [lig], [pose], tmp_path / "sdf")
    reread = pmd.load_file(out["ligand_topology"], xyz=out["ligand_coordinates"])
    assert [a.charge for a in reread.atoms] == pytest.approx([a.charge for a in template.atoms])
    assert out["ligands"][0]["coordinate_format"] == "sdf_first_pose_explicit_graph"
    assert sorted(map(tuple, np.round(reread.coordinates, 3))) == sorted(
        map(tuple, np.round(template.coordinates, 3))
    )


def test_sdf_missing_h_and_changed_connectivity_are_rejected(tmp_path):
    from rdkit import Chem

    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    template = pmd.load_file(str(lig), xyz=str(lc))
    pose = tmp_path / "incomplete.sdf"
    mol = sdf_pose(template, pose)
    rw = Chem.RWMol(mol)
    rw.RemoveAtom(4)
    writer = Chem.SDWriter(str(pose))
    writer.write(rw.GetMol())
    writer.close()
    with pytest.raises(ValueError, match="atom|hydrogen|graph"):
        assemble_complex(rec, rc, [lig], [pose], tmp_path / "partial")
    rw = Chem.RWMol(mol)
    rw.RemoveBond(0, 1)
    rw.AddBond(1, 2, Chem.BondType.SINGLE)
    writer = Chem.SDWriter(str(pose))
    writer.write(rw.GetMol())
    writer.close()
    with pytest.raises(ValueError, match="graph|hydrogen|connectivity"):
        assemble_complex(rec, rc, [lig], [pose], tmp_path / "changed")


def test_sdf_distinct_resp_hydrogens_require_explicit_name_map(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    template = pmd.load_file(str(lig), xyz=str(lc))
    template.atoms[1].charge += 0.01
    template.atoms[2].charge -= 0.01
    changed = tmp_path / "resp.parm7"
    template.save(str(changed))
    pose = tmp_path / "resp.sdf"
    sdf_pose(template, pose, [4, 3, 2, 1, 0])
    with pytest.raises(ValueError, match="symmetr|equivalent|ambig|name map"):
        assemble_complex(rec, rc, [changed], [pose], tmp_path / "ambiguous")
    sdf_pose(template, pose, [4, 3, 2, 1, 0], names=True)
    out = assemble_complex(rec, rc, [changed], [pose], tmp_path / "named")
    reread = pmd.load_file(out["ligand_topology"], xyz=out["ligand_coordinates"])
    assert reread.coordinates == pytest.approx(template.coordinates, abs=0.0001)
    assert [a.charge for a in reread.atoms] == pytest.approx([a.charge for a in template.atoms])


def chiral_template():
    structure = pmd.Structure()
    for i, element in enumerate([6, 9, 17, 35, 53]):
        atom = pmd.Atom(
            name=["C1", "F1", "CL1", "BR1", "I1"][i],
            atomic_number=element,
            type=f"X{i}",
            mass=float(element + 1),
            charge=0.1 * i,
        )
        structure.add_atom(atom, "CHI", 1)
    for i in range(1, 5):
        structure.bonds.append(pmd.Bond(structure.atoms[0], structure.atoms[i]))
    structure.coordinates = np.asarray(
        [[0, 0, 5], [1, 1, 6], [1, -1, 4], [-1, 1, 4], [-1, -1, 6]], dtype=float
    )
    return structure


def test_sdf_tetrahedral_stereoinversion_requires_and_checks_reference(tmp_path):
    from go_builder.complexes import _read_pose

    template = chiral_template()
    pose = tmp_path / "chiral.sdf"
    sdf_pose(template, pose)
    # A prmtop alone carries no reliable reference geometry.
    bare = template.copy(pmd.Structure)
    bare.coordinates = None
    with pytest.raises(ValueError, match="reference|stereo"):
        _read_pose(pose, bare)
    coords, _ = _read_pose(pose, template)
    assert coords == pytest.approx(template.coordinates)
    sdf_pose(template, pose, invert=True)
    with pytest.raises(ValueError, match="stereo|chiral|inver"):
        _read_pose(pose, template)


@pytest.mark.parametrize("difference", ["charge", "mass", "type"])
def test_sdf_heavy_symmetry_cannot_exchange_unequal_amber_atoms(tmp_path, difference):
    from rdkit import Chem
    from rdkit.Chem import AllChem

    from go_builder.complexes import _read_pose

    molecule = Chem.AddHs(Chem.MolFromSmiles("CC"))
    assert AllChem.EmbedMolecule(molecule, randomSeed=19) == 0
    structure = pmd.Structure()
    for atom in molecule.GetAtoms():
        index = atom.GetIdx()
        number = atom.GetAtomicNum()
        structure.add_atom(
            pmd.Atom(
                name=f"A{index}",
                atomic_number=number,
                type="c3" if number == 6 else "hc",
                mass=12.01 if number == 6 else 1.008,
                charge=0.0,
            ),
            "ETH",
            1,
        )
    for bond in molecule.GetBonds():
        structure.bonds.append(
            pmd.Bond(structure.atoms[bond.GetBeginAtomIdx()], structure.atoms[bond.GetEndAtomIdx()])
        )
    structure.coordinates = molecule.GetConformer().GetPositions()
    setattr(
        structure.atoms[1], difference, {"charge": 0.1, "mass": 13.01, "type": "cx"}[difference]
    )
    pose = tmp_path / "ethane.sdf"
    sdf_pose(structure, pose)
    with pytest.raises(ValueError, match="symmetry|unequal|name map"):
        _read_pose(pose, structure)


def test_sdf_first_pose_and_mol_keep_graph_without_hydrogen_autofill(tmp_path):
    from rdkit import Chem

    from go_builder.complexes import _read_pose

    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    structure = pmd.load_file(str(lig), xyz=str(lc))
    sdf = tmp_path / "poses.sdf"
    mol = sdf_pose(structure, sdf)
    sdf.write_text(sdf.read_text() + "broken later pose ignored\n$$$$\n")
    coords, kind = _read_pose(sdf, structure)
    assert kind == "sdf_first_pose_explicit_graph"
    assert len(coords) == len(structure.atoms)
    mol_path = tmp_path / "single.mol"
    Chem.MolToMolFile(mol, str(mol_path))
    assert _read_pose(mol_path, structure)[0] == pytest.approx(coords)


def test_reference_coordinate_option_is_length_checked(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    with pytest.raises(ValueError, match="reference"):
        assemble_complex(
            rec, rc, [lig], [lc], tmp_path / "bad_reference", ligand_reference_coordinates=[]
        )
    pose = tmp_path / "lig.sdf"
    sdf_pose(pmd.load_file(str(lig), xyz=str(lc)), pose)
    out = assemble_complex(
        rec, rc, [lig], [pose], tmp_path / "reference", ligand_reference_coordinates=[lc]
    )
    assert out["ligands"][0]["reference_coordinates"] == str(lc.resolve())
    assert out["ligands"][0]["reference_coordinates_sha256"]


def test_corrupt_amber_topology_is_explanatory_value_error_before_outputs(tmp_path):
    rec, rc = methane(tmp_path, "rec")
    lig, lc = methane(tmp_path, "lig", (0, 0, 5))
    lines = lig.read_text().splitlines()
    index = lines.index("%FLAG ATOM_NAME") + 2
    lines[index] = lines[index][:-4]
    corrupted = tmp_path / "corrupt.parm7"
    corrupted.write_text("\n".join(lines) + "\n")
    with pytest.raises(ValueError, match="corrupt Amber|Invalid.*Amber|FLAG"):
        assemble_complex(rec, rc, [corrupted], [lc], tmp_path / "bad_topology")
    assert not (tmp_path / "bad_topology").exists()
