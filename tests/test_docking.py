import hashlib
import json

import pytest

from go_builder.docking import _gnina_score, _validate_pose_heavy_atom_count, dock_boxes


def _write_sdf(path, *, affinity=None, molecule="C"):
    from rdkit import Chem

    mol = Chem.MolFromSmiles(molecule)
    conf = Chem.Conformer(mol.GetNumAtoms())
    conf.Set3D(True)
    for index in range(mol.GetNumAtoms()):
        conf.SetAtomPosition(index, (float(index), 0.0, 0.0))
    mol.AddConformer(conf)
    if affinity is not None:
        mol.SetProp("minimizedAffinity", str(affinity))
    writer = Chem.SDWriter(str(path))
    writer.write(mol)
    writer.close()


def _pdb_atom(serial, name, x, y, z, element, record="ATOM"):
    return (
        f"{record:<6}{serial:5d} {name:>4s} ALA A{serial:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {element:>2s}  \n"
    )


def _pdbqt_atom(serial, name, x, y, z, atom_type):
    return (
        f"ATOM  {serial:5d} {name:<4s} LIG     1    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}    0.000 {atom_type}\n"
    )


def _inputs(tmp_path, *, coords=((1.0, 2.0, 3.0, "C"), (4.0, 5.0, 6.0, "O"))):
    tmp_path.mkdir(parents=True, exist_ok=True)
    receptor = tmp_path / "receptor.pdb"
    receptor.write_text(
        "".join(
            _pdb_atom(i, name, x, y, z, element)
            for i, (name, x, y, z, element) in enumerate(
                (("C1", *coords[0]), ("O1", *coords[1])), 1
            )
        )
    )
    pdbqt = tmp_path / "receptor.pdbqt"
    # Deliberately reorder atoms and omit hydrogen: identity is element+coordinate, not order.
    pdbqt.write_text(
        _pdbqt_atom(2, "O1", *coords[1][:3], "OA") + _pdbqt_atom(1, "C1", *coords[0][:3], "C")
    )
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text(_pdbqt_atom(1, "C1", 0.0, 0.0, 0.0, "C"))
    boxes_file = tmp_path / "boxes.json"
    (tmp_path / "boxes").mkdir()
    config = tmp_path / "boxes" / "box_001.txt"
    config.write_text(
        "center_x = 2\ncenter_y = 3\ncenter_z = 4\nsize_x = 10\nsize_y = 11\nsize_z = 12\n"
    )
    payload = {
        "schema_version": 1,
        "receptor_pdb": "receptor.pdb",
        "receptor_pdb_sha256": hashlib.sha256(receptor.read_bytes()).hexdigest(),
        "boxes": [
            {
                "id": "box_001",
                "center": [2, 3, 4],
                "size": [10, 11, 12],
                "config": "boxes/box_001.txt",
            }
        ],
        "source_frame": 1,
    }
    boxes_file.write_text(json.dumps(payload))
    return boxes_file, pdbqt, ligand, payload


def test_prepare_only_writes_snapshot_matched_configs_and_literal_commands(tmp_path):
    boxes, receptor, ligand, payload = _inputs(tmp_path)
    out = tmp_path / "dock"

    manifest = dock_boxes(boxes, receptor, ligand, out, cpus=3, runs=2, run=False)

    assert manifest["status"] == "inputs_prepared"
    assert manifest["docking_completed"] is False
    assert len(json.loads((out / "commands.json").read_text())["commands"]) == 2
    assert (out / "dock_boxes.sh").stat().st_mode & 0o111
    cfg = (out / "boxes" / "box_001.txt").read_text()
    assert "center_z = 4" in cfg and "size_y = 11" in cfg
    assert payload["receptor_pdb_sha256"] in (out / "manifest.json").read_text()


def test_run_uses_allocated_cpus_and_incrementing_seeds_and_parses_best_scores(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)
    fake = tmp_path / "fake-vina"
    log = tmp_path / "argv.jsonl"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import json, os, sys\n"
        "from pathlib import Path\n"
        "args=sys.argv[1:]\n"
        f"Path({str(log)!r}).open('a').write(json.dumps(args)+'\\n')\n"
        "out=Path(args[args.index('--out')+1])\n"
        "out.write_text('REMARK VINA RESULT: -7.5 0.0 0.0\\nATOM      1 C    LIG     1       0.000   0.000   0.000  0.00  0.00    0.000 C\\n')\n"
    )
    fake.chmod(0o755)
    out = tmp_path / "dock"

    manifest = dock_boxes(
        boxes,
        receptor,
        ligand,
        out,
        vina=str(fake),
        cpus=4,
        exhaustiveness=8,
        num_modes=5,
        runs=2,
        seed=41,
        run=True,
    )

    commands = [json.loads(line) for line in log.read_text().splitlines()]
    assert [args[args.index("--cpu") + 1] for args in commands] == ["4", "4"]
    assert [args[args.index("--seed") + 1] for args in commands] == ["41", "42"]
    assert manifest["status"] == "docking_completed"
    assert manifest["docking_completed"] is True
    assert (out / "summary.csv").read_text().splitlines()[1].split(",")[1] == "-7.5"


def test_rejects_truncated_or_different_receptor_before_creating_output(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)
    receptor.write_text(_pdbqt_atom(1, "C1", 1.0, 2.0, 3.0, "C"))
    out = tmp_path / "must-not-exist"

    with pytest.raises(ValueError, match="heavy[- ]atoms?"):
        dock_boxes(boxes, receptor, ligand, out)

    assert not out.exists()


def test_rejects_stale_config_and_bad_box_geometry_without_output(tmp_path):
    boxes, receptor, ligand, payload = _inputs(tmp_path)
    (tmp_path / "boxes" / "box_001.txt").write_text("center_x = 99\n")
    out = tmp_path / "dock"

    with pytest.raises(ValueError, match="config"):
        dock_boxes(boxes, receptor, ligand, out)
    assert not out.exists()

    payload["boxes"][0]["size"][0] = 0
    (tmp_path / "boxes" / "box_001.txt").write_text(
        "center_x = 2\ncenter_y = 3\ncenter_z = 4\nsize_x = 10\nsize_y = 11\nsize_z = 12\n"
    )
    boxes.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="size"):
        dock_boxes(boxes, receptor, ligand, out)
    assert not out.exists()

    payload["boxes"][0]["size"][0] = 10
    (tmp_path / "boxes" / "box_001.txt").write_text(
        "center_x = 2\ncenter_y = 3\ncenter_z = 4\nsize_x = 10\n"
        "size_y = 11\nsize_z = 12\nreceptor = secret.pdbqt\n"
    )
    boxes.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="unexpected or duplicate key"):
        dock_boxes(boxes, receptor, ligand, out)
    assert not out.exists()


def test_rejects_hash_mismatch_and_nonempty_output_directory(tmp_path):
    boxes, receptor, ligand, payload = _inputs(tmp_path)
    payload["receptor_pdb_sha256"] = "0" * 64
    boxes.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="SHA-256"):
        dock_boxes(boxes, receptor, ligand, tmp_path / "out")

    boxes, receptor, ligand, _ = _inputs(tmp_path / "other")
    occupied = tmp_path / "occupied"
    occupied.mkdir()
    (occupied / "keep.txt").write_text("original")
    with pytest.raises(ValueError, match="empty directory"):
        dock_boxes(boxes, receptor, ligand, occupied)
    assert (occupied / "keep.txt").read_text() == "original"


def test_run_failure_marks_manifest_failed(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)
    fake = tmp_path / "fake-failure"
    fake.write_text("#!/bin/sh\necho deliberate failure >&2\nexit 7\n")
    fake.chmod(0o755)
    out = tmp_path / "dock"

    with pytest.raises(RuntimeError, match="Vina failed"):
        dock_boxes(boxes, receptor, ligand, out, vina=str(fake), run=True)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["docking_completed"] is False
    assert "deliberate failure" in (out / "box_001_run_001.stderr.txt").read_text()


def test_run_success_without_pose_is_marked_failed(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)
    fake = tmp_path / "fake-no-pose"
    fake.write_text("#!/bin/sh\necho completed without writing output\n")
    fake.chmod(0o755)
    out = tmp_path / "dock"

    with pytest.raises(RuntimeError, match="no nonempty output"):
        dock_boxes(boxes, receptor, ligand, out, vina=str(fake), run=True)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["docking_completed"] is False


def test_missing_vina_binary_leaves_failed_manifest(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)
    out = tmp_path / "dock"

    with pytest.raises(FileNotFoundError, match="Vina executable not found"):
        dock_boxes(boxes, receptor, ligand, out, vina="missing-vina-executable", run=True)

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["docking_completed"] is False


def test_large_receptor_identity_check_is_iterative_and_bijective(tmp_path):
    count = 1200
    coords = [
        ((index % 20) * 4.0, ((index // 20) % 20) * 4.0, (index // 400) * 4.0)
        for index in range(count)
    ]
    receptor_pdb = tmp_path / "receptor.pdb"
    receptor_pdb.write_text(
        "".join(_pdb_atom(index + 1, "C", *xyz, "C") for index, xyz in enumerate(coords))
    )
    receptor_pdbqt = tmp_path / "receptor.pdbqt"
    receptor_pdbqt.write_text(
        "".join(
            _pdbqt_atom(index + 1, "C", *xyz, "C")
            for index, xyz in reversed(list(enumerate(coords)))
        )
    )
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text(_pdbqt_atom(1, "C", 0.0, 0.0, 0.0, "C"))
    box_dir = tmp_path / "boxes"
    box_dir.mkdir()
    (box_dir / "box_001.txt").write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\nsize_x = 10\nsize_y = 10\nsize_z = 10\n"
    )
    boxes = tmp_path / "boxes.json"
    boxes.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "receptor_pdb": "receptor.pdb",
                "receptor_pdb_sha256": hashlib.sha256(receptor_pdb.read_bytes()).hexdigest(),
                "boxes": [
                    {
                        "id": "box_001",
                        "center": [0, 0, 0],
                        "size": [10, 10, 10],
                        "config": "boxes/box_001.txt",
                    }
                ],
            }
        )
    )

    manifest = dock_boxes(boxes, receptor_pdbqt, ligand, tmp_path / "dock")

    assert manifest["status"] == "inputs_prepared"


def test_receptor_matching_uses_augmenting_path_for_ambiguous_coordinates(tmp_path):
    receptor_pdb = tmp_path / "receptor.pdb"
    receptor_pdb.write_text(
        _pdb_atom(1, "C1", 0.000, 0.0, 0.0, "C") + _pdb_atom(2, "C2", 0.019, 0.0, 0.0, "C")
    )
    receptor_pdbqt = tmp_path / "receptor.pdbqt"
    # Expected atom 1 can match either PDBQT atom; atom 2 can match only atom 1.
    receptor_pdbqt.write_text(
        _pdbqt_atom(1, "C1", 0.005, 0.0, 0.0, "C") + _pdbqt_atom(2, "C2", -0.015, 0.0, 0.0, "C")
    )
    ligand = tmp_path / "ligand.pdbqt"
    ligand.write_text(_pdbqt_atom(1, "C", 0.0, 0.0, 0.0, "C"))
    box_dir = tmp_path / "boxes"
    box_dir.mkdir()
    (box_dir / "box_001.txt").write_text(
        "center_x = 0\ncenter_y = 0\ncenter_z = 0\nsize_x = 10\nsize_y = 10\nsize_z = 10\n"
    )
    boxes = tmp_path / "boxes.json"
    boxes.write_text(
        json.dumps(
            {
                "schema_version": 1,
                "receptor_pdb": "receptor.pdb",
                "receptor_pdb_sha256": hashlib.sha256(receptor_pdb.read_bytes()).hexdigest(),
                "boxes": [
                    {
                        "id": "box_001",
                        "center": [0, 0, 0],
                        "size": [10, 10, 10],
                        "config": "boxes/box_001.txt",
                    }
                ],
            }
        )
    )

    manifest = dock_boxes(boxes, receptor_pdbqt, ligand, tmp_path / "dock")

    assert manifest["status"] == "inputs_prepared"


def test_gnina_prepare_uses_snapshot_pdb_sdf_and_expected_flags(tmp_path):
    boxes, _, _, _ = _inputs(tmp_path)
    ligand_sdf = tmp_path / "ligand.sdf"
    _write_sdf(ligand_sdf)

    manifest = dock_boxes(
        boxes,
        tmp_path / "receptor.pdb",
        ligand_sdf,
        tmp_path / "gnina-plan",
        engine="gnina",
        gnina="gnina",
        cpus=2,
        seed=7,
    )

    command = json.loads((tmp_path / "gnina-plan" / "commands.json").read_text())["commands"][0][
        "argv"
    ]
    assert manifest["status"] == "inputs_prepared"
    assert command[command.index("--receptor") + 1].endswith("receptor.pdb")
    assert command[command.index("--ligand") + 1].endswith("ligand.sdf")
    assert command[command.index("--cpu") + 1] == "2"
    assert command[command.index("--seed") + 1] == "7"
    assert command[command.index("--scoring") + 1] == "vina"
    assert command[command.index("--cnn_scoring") + 1] == "none"
    assert command[command.index("--pose_sort_order") + 1] == "Energy"
    assert "--no_gpu" in command
    assert command[command.index("--out") + 1].endswith(".sdf")


def test_gnina_accepts_snapshot_matched_pdbqt_inputs(tmp_path):
    boxes, receptor, ligand, _ = _inputs(tmp_path)

    manifest = dock_boxes(
        boxes,
        receptor,
        ligand,
        tmp_path / "gnina-pdbqt-plan",
        engine="gnina",
    )

    assert manifest["engine"] == "gnina"
    assert manifest["receptor_input_format"] == "pdbqt"
    assert manifest["ligand_input_format"] == "pdbqt"


def test_gnina_run_parses_minimized_affinity_and_validates_output_atoms(tmp_path):
    boxes, _, _, _ = _inputs(tmp_path)
    ligand_sdf = tmp_path / "ligand.sdf"
    _write_sdf(ligand_sdf)
    fake = tmp_path / "fake-gnina"
    log = tmp_path / "argv.json"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import json, sys\n"
        "from pathlib import Path\n"
        "from rdkit import Chem\n"
        "args=sys.argv[1:]\n"
        f"Path({str(log)!r}).write_text(json.dumps(args))\n"
        "mol=Chem.MolFromSmiles('C')\n"
        "conf=Chem.Conformer(1)\n"
        "conf.Set3D(True)\n"
        "conf.SetAtomPosition(0,(0.0,0.0,0.0))\n"
        "mol.AddConformer(conf)\n"
        "mol.SetProp('minimizedAffinity','-0.63861')\n"
        "writer=Chem.SDWriter(args[args.index('--out')+1])\n"
        "writer.write(mol)\n"
        "writer.close()\n"
    )
    fake.chmod(0o755)

    manifest = dock_boxes(
        boxes,
        tmp_path / "receptor.pdb",
        ligand_sdf,
        tmp_path / "gnina-run",
        engine="gnina",
        gnina=str(fake),
        cpus=3,
        runs=1,
        seed=31,
        run=True,
    )

    args = json.loads(log.read_text())
    assert args[args.index("--cpu") + 1] == "3"
    assert args[args.index("--seed") + 1] == "31"
    assert args[args.index("--cnn_scoring") + 1] == "none"
    assert manifest["status"] == "docking_completed"
    assert manifest["global_best"]["score"] == pytest.approx(-0.63861)
    assert manifest["score_field"] == "minimizedAffinity"


def test_gnina_score_parser_reads_captured_sdf_pose_record(tmp_path):
    # Captured GNINA result record, used only to verify offline SDF score parsing.
    output = tmp_path / "recorded_pose.sdf"
    output.write_text(
        "\n\n\n  4  3  0  0  0  0  0  0  0  0999 V2000\n"
        "    0.5439    0.4039    3.6412 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.6330   -0.5637    3.6619 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "   -0.1850   -1.8803    3.4882 O   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "    0.0218   -1.9848    2.5227 H   0  0  0  0  0  0  0  0  0  0  0  0\n"
        "  1  2  1  0\n  2  3  1  0\n  3  4  1  0\nM  END\n"
        "> <minimizedAffinity>\n-0.63861\n\n$$$$\n"
    )

    assert _gnina_score(output) == pytest.approx(-0.63861)
    _validate_pose_heavy_atom_count(output, "gnina", expected_count=3)


def test_gnina_rejects_changed_receptor_and_invalid_or_multirecord_sdf_before_output(tmp_path):
    boxes, _, _, _ = _inputs(tmp_path)
    ligand_sdf = tmp_path / "ligand.sdf"
    _write_sdf(ligand_sdf)
    changed = tmp_path / "changed.pdb"
    changed.write_text(
        _pdb_atom(1, "C1", 100.0, 2.0, 3.0, "C") + _pdb_atom(2, "O1", 4.0, 5.0, 6.0, "O")
    )
    out = tmp_path / "invalid-receptor"
    with pytest.raises(ValueError, match="heavy[- ]atoms?"):
        dock_boxes(boxes, changed, ligand_sdf, out, engine="gnina")
    assert not out.exists()

    multi = tmp_path / "multi.sdf"
    _write_sdf(multi)
    with multi.open("a") as handle:
        handle.write(ligand_sdf.read_text())
    out = tmp_path / "invalid-ligand"
    with pytest.raises(ValueError, match="exactly one valid molecule"):
        dock_boxes(boxes, tmp_path / "receptor.pdb", multi, out, engine="gnina")
    assert not out.exists()


def test_gnina_rejects_sdf_without_three_dimensional_coordinates(tmp_path):
    boxes, _, _, _ = _inputs(tmp_path)
    from rdkit import Chem

    mol = Chem.MolFromSmiles("C")
    conformer = Chem.Conformer(mol.GetNumAtoms())
    conformer.Set3D(False)
    mol.AddConformer(conformer)
    ligand_sdf = tmp_path / "two_dimensional.sdf"
    writer = Chem.SDWriter(str(ligand_sdf))
    writer.write(mol)
    writer.close()
    out = tmp_path / "invalid-2d"

    with pytest.raises(ValueError, match="3D coordinates"):
        dock_boxes(boxes, tmp_path / "receptor.pdb", ligand_sdf, out, engine="gnina")
    assert not out.exists()


def test_gnina_rejects_pose_with_missing_ligand_heavy_atoms(tmp_path):
    boxes, _, _, _ = _inputs(tmp_path)
    ligand_sdf = tmp_path / "ligand.sdf"
    _write_sdf(ligand_sdf, molecule="CO")
    fake = tmp_path / "fake-incomplete-gnina"
    fake.write_text(
        "#!/usr/bin/env python3\n"
        "import sys\n"
        "from rdkit import Chem\n"
        "mol=Chem.MolFromSmiles('C')\n"
        "conf=Chem.Conformer(1)\n"
        "conf.Set3D(True)\n"
        "conf.SetAtomPosition(0,(0.0,0.0,0.0))\n"
        "mol.AddConformer(conf)\n"
        "mol.SetProp('minimizedAffinity','-0.5')\n"
        "writer=Chem.SDWriter(sys.argv[sys.argv.index('--out')+1])\n"
        "writer.write(mol)\n"
        "writer.close()\n"
    )
    fake.chmod(0o755)
    out = tmp_path / "missing-pose-atoms"

    with pytest.raises(ValueError, match="heavy-atom count does not match"):
        dock_boxes(
            boxes,
            tmp_path / "receptor.pdb",
            ligand_sdf,
            out,
            engine="gnina",
            gnina=str(fake),
            runs=1,
            run=True,
        )

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["docking_completed"] is False
