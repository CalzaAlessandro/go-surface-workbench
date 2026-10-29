"""Equilibrated conformations keep their graph and receive full 3D box coverage."""

import json
from itertools import product

import numpy as np
import parmed as pmd
import pytest
from test_complexes import methane

from go_builder.snapshot import cover_boxes, make_whole, prepare_docking_receptor


def test_wrapped_finite_receptor_is_made_whole_without_flattening_or_retyping(tmp_path):
    top, coordinates = methane(tmp_path, "source", (9.7, 0, 0))
    original = pmd.load_file(str(top), xyz=str(coordinates))
    restart = pmd.amber.Rst7(natom=5)
    restart.box = [10, 10, 10, 90, 90, 90]
    restart.coordinates = np.mod(original.coordinates, 10)
    wrapped = tmp_path / "wrapped.rst7"
    restart.write(str(wrapped))
    output = tmp_path / "docking"
    result = prepare_docking_receptor(
        top, wrapped, output, receptor_mask=":1", box_size=[12, 12, 12], overlap=4, padding=3
    )
    receptor = pmd.load_file(result["topology"], xyz=result["coordinates"])
    assert receptor.coordinates == pytest.approx(original.coordinates, abs=1e-6)
    assert [a.charge for a in receptor.atoms] == pytest.approx([-0.4, 0.1, 0.1, 0.1, 0.1])
    assert len(receptor.bonds) == 4
    assert (output / "receptor.pdb").read_text().count("CONECT") == 4
    assert json.loads((output / "boxes.json").read_text())["receptor_pdb"] == "receptor.pdb"
    assert result["flattened"] is False


def test_selection_cannot_cut_a_covalent_graft_or_surface_bond(tmp_path):
    top, coordinates = methane(tmp_path, "source")
    output = tmp_path / "rejected"
    with pytest.raises(ValueError, match="covalent|bond|cut"):
        prepare_docking_receptor(top, coordinates, output, receptor_mask="@1")
    assert not output.exists()


def test_3d_boxes_cover_folded_extent_including_z_and_padding():
    xyz = np.asarray([[-10, 0, 0], [0, 20, 15], [10, 0, 30]], dtype=float)
    boxes = cover_boxes(xyz, box_size=[12, 12, 12], overlap=3, padding=2)
    assert len({round(b["center"][2], 6) for b in boxes}) > 1
    for point in product([-12, 12], [-2, 22], [-2, 32]):
        assert any(
            np.all(np.abs(np.asarray(point) - b["center"]) <= np.asarray(b["size"]) / 2 + 1e-8)
            for b in boxes
        )


def test_covalently_winding_graph_cannot_be_unwrapped_as_finite():
    xyz = np.asarray([[0, 0, 0], [3.333333, 0, 0], [6.666667, 0, 0]])
    with pytest.raises(ValueError, match="winding|infinite|periodic"):
        make_whole(xyz, [(0, 1), (1, 2), (0, 2)], np.diag([10.0, 10.0, 10.0]))


def test_netcdf_frame_is_one_based_and_preserves_actual_three_dimensional_pose(tmp_path):
    from netCDF4 import Dataset

    top, coordinates = methane(tmp_path, "source")
    original = pmd.load_file(str(top), xyz=str(coordinates))
    trajectory = tmp_path / "production.nc"
    with Dataset(trajectory, "w") as nc:
        nc.createDimension("frame", 2)
        nc.createDimension("atom", 5)
        nc.createDimension("spatial", 3)
        nc.createVariable("coordinates", "f8", ("frame", "atom", "spatial"))[:] = [
            original.coordinates,
            original.coordinates + [0, 0, 7],
        ]
        nc.createVariable("time", "f8", ("frame",))[:] = [10, 20]
    result = prepare_docking_receptor(top, trajectory, tmp_path / "frame2", frame=2)
    receptor = pmd.load_file(result["topology"], xyz=result["coordinates"])
    assert receptor.coordinates == pytest.approx(original.coordinates + [0, 0, 7], abs=1e-6)
    assert result["source_frame"] == 2
    assert result["source_time_ps"] == 20
    with pytest.raises(ValueError, match="frame"):
        prepare_docking_receptor(top, trajectory, tmp_path / "badframe", frame=3)


@pytest.mark.parametrize(
    "options", [{"overlap": 12}, {"padding": -1}, {"box_size": [0, 12, 12]}, {"max_boxes": 1}]
)
def test_invalid_or_excessive_box_requests_are_explicit(options):
    with pytest.raises(ValueError):
        cover_boxes(
            np.asarray([[0, 0, 0], [40, 40, 40]]),
            box_size=options.pop("box_size", [12, 12, 12]),
            **options,
        )
