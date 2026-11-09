"""Unparameterized ether/lactone/lactol references; no production assignment.

Graphs are explicit neutral C/H/O molecules. MMFF94s geometries and Gaussian
templates are starting artifacts, never accepted QM, fitted charges, or MD data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

from .protocols import DFT_METHOD
from .qm import gaussian_input
from .reference_qm import make_sp_input
from .resp import _validate_model

SEED = 20261105
STATUS = "UNPARAMETERIZED_RESEARCH_REFERENCE"
CANDIDATES = (
    {"model_id": "dibenzofuran_ether_edge", "name": "dibenzofuran", "smiles": "c1ccc2c(c1)oc1ccccc12", "motif": "non_epoxide_ether", "ring_size": 5, "formula": "C12H8O", "inchikey": "TXCDCPKCNAJMEE-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:28145", "reference_scope": "fused_aromatic_edge_motif_reference"},
    {"model_id": "diphenyl_ether_bridge_control", "name": "diphenyl ether", "smiles": "c1ccc(Oc2ccccc2)cc1", "motif": "non_epoxide_ether", "ring_size": None, "formula": "C12H10O", "inchikey": "USIUVYZYUHIAEV-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:39258", "reference_scope": "diaryl_bridge_connectivity_control"},
    {"model_id": "phthalide_lactone5_edge", "name": "phthalide", "smiles": "O=C1OCc2ccccc21", "motif": "lactone", "ring_size": 5, "formula": "C8H6O2", "inchikey": "WNZQDUSMALZDQF-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:38085", "reference_scope": "fused_aromatic_edge_motif_reference"},
    {"model_id": "coumarin_lactone6_edge", "name": "coumarin", "smiles": "O=c1ccc2ccccc2o1", "motif": "lactone", "ring_size": 6, "formula": "C9H6O2", "inchikey": "ZYGHJZDHTFUPRJ-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:28794", "reference_scope": "fused_aromatic_edge_motif_reference"},
    {"model_id": "oxolan2ol_lactol_control", "name": "oxolan-2-ol", "smiles": "OC1CCCO1", "motif": "lactol", "ring_size": 5, "formula": "C4H8O2", "inchikey": "JNODDICFTDYODH-UHFFFAOYSA-N", "source": "https://jglobal.jst.go.jp/en/detail?JGLOBAL_ID=200907060640520321", "reference_scope": "aliphatic_connectivity_control"},
    {"model_id": "phenol_aromatic_OH_control", "name": "phenol", "smiles": "Oc1ccccc1", "motif": "phenolic_OH", "ring_size": None, "formula": "C6H6O", "inchikey": "ISWSIDIOOBJBQZ-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:15882", "reference_scope": "monocyclic_aromatic_OH_control"},
    {"model_id": "naphthoquinone14_edge", "name": "1,4-naphthoquinone", "smiles": "O=C1C=CC(=O)c2ccccc21", "motif": "quinone", "ring_size": 6, "formula": "C10H6O2", "inchikey": "FRASJONUBLZVQX-UHFFFAOYSA-N", "source": "https://www.ebi.ac.uk/chebi/CHEBI:27418", "reference_scope": "fused_aromatic_edge_motif_reference"},
)


def validate_motif(mol, motif, *, ring_size=None):
    """Fail closed on wrong connectivity; return explicit one-based motif IDs.

    Aromatic C/O are allowed for dibenzofuran and coumarin; RDKit can encode
    conjugated C–O bonds as aromatic rather than literal single bonds.
    """
    Chem.SanitizeMol(mol)
    rings = mol.GetRingInfo().AtomRings()
    carbonyls = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
                 for b in mol.GetBonds() if b.GetBondTypeAsDouble() == 2
                 and {b.GetBeginAtom().GetSymbol(), b.GetEndAtom().GetSymbol()} == {"C", "O"}]
    carbonyl_c = {next(i for i in pair if mol.GetAtomWithIdx(i).GetSymbol() == "C") for pair in carbonyls}
    matches = []
    if motif == "non_epoxide_ether":
        for oxygen in mol.GetAtoms():
            heavy = [n for n in oxygen.GetNeighbors() if n.GetAtomicNum() > 1]
            if oxygen.GetSymbol() != "O" or oxygen.GetFormalCharge() or len(heavy) != 2 or any(n.GetSymbol() != "C" for n in heavy):
                continue
            if any(n.GetIdx() in carbonyl_c for n in heavy):
                continue
            sizes = [len(r) for r in rings if oxygen.GetIdx() in r]
            if 3 in sizes or (ring_size is not None and ring_size not in sizes):
                continue
            matches.append([oxygen.GetIdx(), *(n.GetIdx() for n in heavy)])
    elif motif == "lactone":
        for c in carbonyl_c:
            oxygen_double = next(i for pair in carbonyls if c in pair for i in pair if i != c)
            for ring in rings:
                if c not in ring or (ring_size is not None and len(ring) != ring_size):
                    continue
                for n in mol.GetAtomWithIdx(c).GetNeighbors():
                    if n.GetSymbol() == "O" and n.GetIdx() in ring and sum(a.GetSymbol() == "C" for a in n.GetNeighbors()) == 2:
                        matches.append([c, oxygen_double, n.GetIdx()])
    elif motif == "lactol":
        if not carbonyl_c:
            for carbon in mol.GetAtoms():
                if carbon.GetSymbol() != "C" or carbon.GetHybridization() != Chem.HybridizationType.SP3 or not carbon.IsInRing():
                    continue
                oxygens = [n for n in carbon.GetNeighbors() if n.GetSymbol() == "O"]
                hydroxyl = [o for o in oxygens if o.GetTotalNumHs(includeNeighbors=True) == 1]
                ether = [o for o in oxygens if o.IsInRing() and sum(n.GetSymbol() == "C" for n in o.GetNeighbors()) == 2]
                if len(hydroxyl) == len(ether) == 1 and hydroxyl[0].GetIdx() != ether[0].GetIdx() and all(mol.GetBondBetweenAtoms(carbon.GetIdx(), o.GetIdx()).GetBondTypeAsDouble() == 1 for o in oxygens):
                    if ring_size is None or any(len(r) == ring_size and carbon.GetIdx() in r and ether[0].GetIdx() in r for r in rings):
                        matches.append([carbon.GetIdx(), hydroxyl[0].GetIdx(), ether[0].GetIdx()])
    elif motif == "phenolic_OH":
        for oxygen in mol.GetAtoms():
            if oxygen.GetSymbol() == "O" and oxygen.GetTotalNumHs(includeNeighbors=True) == 1:
                for n in oxygen.GetNeighbors():
                    if n.GetSymbol() == "C" and n.GetIsAromatic():
                        matches.append([n.GetIdx(), oxygen.GetIdx()])
    elif motif == "quinone":
        for ring in rings:
            cs = carbonyl_c.intersection(ring)
            if len(ring) == 6 and len(cs) == 2 and all(mol.GetAtomWithIdx(i).GetSymbol() == "C" for i in ring):
                # The catalog uses para quinone: opposite sites in its six-cycle.
                a, b = sorted(cs)
                if min((ring.index(a)-ring.index(b)) % 6, (ring.index(b)-ring.index(a)) % 6) == 3:
                    if sum(mol.GetBondBetweenAtoms(ring[i], ring[(i+1) % 6]).GetBondTypeAsDouble() >= 1.5 for i in range(6)) >= 2:
                        matches.append(sorted(cs | {i for pair in carbonyls for i in pair}))
    else:
        raise ValueError(f"Unknown motif: {motif}")
    if not matches:
        raise ValueError(f"Required motif/ring gate failed: {motif}, ring_size={ring_size}; epoxide is not general ether")
    return {"passed": True, "motif": motif, "ring_size": ring_size,
            "matches_atom_ids": [[i+1 for i in match] for match in matches],
            "carbonyl_carbon_atom_ids": sorted(i+1 for i in carbonyl_c)}


def build_candidate(name, *, seed=SEED):
    """Return (explicit-H RDKit molecule, serializable research model)."""
    try:
        spec = next(s for s in CANDIDATES if s["model_id"] == name)
    except StopIteration as exc:
        raise ValueError(f"Unknown extended reference: {name}") from exc
    if type(seed) is not int or not 0 <= seed < 2**31:
        raise ValueError("Seed must be a nonnegative signed 32-bit integer")
    mol = Chem.MolFromSmiles(spec["smiles"])
    if mol is None:
        raise ValueError("Invalid catalog SMILES")
    Chem.SanitizeMol(mol)
    if rdMolDescriptors.CalcMolFormula(mol) != spec["formula"] or Chem.MolToInchiKey(mol) != spec["inchikey"]:
        raise ValueError("Catalog source identity mismatch")
    mol = Chem.AddHs(mol)
    if Chem.GetFormalCharge(mol) != 0 or any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError("Only neutral zero-radical reference graphs are permitted")
    electrons = sum(a.GetAtomicNum() for a in mol.GetAtoms())
    if electrons % 2:
        raise ValueError("Closed-shell singlet electron parity gate failed")
    gate = validate_motif(mol, spec["motif"], ring_size=spec["ring_size"])
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    if AllChem.EmbedMolecule(mol, params) < 0:
        raise ValueError("ETKDG embedding failed")
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError("Missing MMFF94s starting-geometry parameters")
    status = AllChem.MMFFOptimizeMolecule(mol, mmffVariant="MMFF94s", maxIters=3000)
    if status != 0:
        raise ValueError("MMFF94s starting minimization did not converge")
    # MMFF assigns its own aromaticity. Restore the standard RDKit graph before
    # classifying roles/ranks so JSON and sanitized SDF encode the same bonds.
    Chem.SanitizeMol(mol)
    gate = validate_motif(mol, spec["motif"], ring_size=spec["ring_size"])
    Chem.AssignAtomChiralTagsFromStructure(mol, confId=0, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    unlabelled = Chem.Mol(mol)
    ranks = list(Chem.CanonicalRankAtoms(unlabelled, breakTies=False))
    motifids = {i for match in gate["matches_atom_ids"] for i in match}
    carbonylids = set(gate["carbonyl_carbon_atom_ids"])
    positions = mol.GetConformer().GetPositions()
    if not all(math.isfinite(float(x)) for point in positions for x in point):
        raise ValueError("Nonfinite starting coordinates")
    roles = {}
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()+1
        if atom.GetSymbol() == "H":
            continue
        if idx in motifids:
            suffix = "carbonyl_C" if idx in carbonylids else atom.GetSymbol()
            if atom.GetSymbol() == "O" and spec["motif"] == "lactone":
                suffix = "ring_O" if atom.IsInRing() else "carbonyl_O"
            roles[idx] = f"extended_{spec['motif']}_{suffix}"
        else:
            roles[idx] = "extended_aromatic_C" if atom.GetIsAromatic() else f"extended_background_{atom.GetSymbol()}"
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == "H":
            roles[atom.GetIdx()+1] = "extended_H_on_"+roles[atom.GetNeighbors()[0].GetIdx()+1].removeprefix("extended_")
        atom.SetAtomMapNum(atom.GetIdx()+1)
        atom.SetProp("role", roles[atom.GetIdx()+1])
    atoms = [{"id": a.GetIdx()+1, "name": f"{a.GetSymbol()}{a.GetIdx()+1}", "element": a.GetSymbol(), "position": positions[a.GetIdx()].tolist(), "role": roles[a.GetIdx()+1], "roles": [roles[a.GetIdx()+1]], "formal_charge": a.GetFormalCharge(), "symmetry_class": int(ranks[a.GetIdx()]), "chiral_tag": str(a.GetChiralTag()), "cip_code": a.GetProp("_CIPCode") if a.HasProp("_CIPCode") else None} for a in mol.GetAtoms()]
    clean = Chem.RemoveHs(Chem.Mol(mol))
    for atom in clean.GetAtoms():
        atom.SetAtomMapNum(0)
    metadata = {**spec, "source_access_date": "2026-11-05", "source_identity_verification": "exact formula and nonstereospecific InChIKey from documented SMILES before embedding", "parameter_status": STATUS, "modular_partition_supported": False, "builder_assignment_supported": False, "geometry_status": "MMFF94s_starting_geometry_not_QM_accepted", "initial_geometry_method": "ETKDGv3 + MMFF94s", "optimization_status": "initial_force_field_minimization_converged", "coordinate_units": "angstrom", "seed": seed, "rdkit_version": rdBase.rdkitVersion, "motif_graph_gate": gate, "canonical_isomeric_SMILES": Chem.MolToSmiles(clean), "charge_status": "formal integers only; no fitted partial charges", "wavefunction_stability_validated": False, "multiplicity_basis": "neutral zero-radical even-electron graph; closed-shell singlet starting ansatz, not ground-state validation", "stereo_status": "input unspecified where applicable; one deterministic 3D stereoisomer recorded, no population claim", "scope": "unparameterized_small_molecule_reference_not_full_sheet_or_size_transfer_validation", "size_transfer_requirement": "preserve graph classes and molecular/group charge constraints from about 4 nm to about 10 nm; held-out larger-sheet validation required; no whole-sheet QM required for each production size", "experimental_presence_or_abundance_claimed": False}
    record = {"schema_version": 1, "model_id": name, "formula": spec["formula"], "formal_charge": 0, "multiplicity": 1, "electron_count": electrons, "atoms": atoms, "bonds": [{"a": b.GetBeginAtomIdx()+1, "b": b.GetEndAtomIdx()+1, "order": b.GetBondTypeAsDouble()} for b in mol.GetBonds()], "metadata": metadata}
    _validate_model(record)
    mol.SetProp("_Name", name)
    mol.SetProp("parameter_status", STATUS)
    mol.SetProp("input_smiles", spec["smiles"])
    return mol, record


def export_candidates(output_dir, *, seed=SEED, cpus=8, mem_gb=6):
    """Create a fresh complete reference panel; refuse an existing output path."""
    folder = Path(output_dir)
    if folder.exists():
        raise FileExistsError(f"Refusing to overwrite reference panel: {folder}")
    # Build and validate everything before writing any candidate files.
    prepared = []
    for spec in CANDIDATES:
        mol, record = build_candidate(spec["model_id"], seed=seed)
        dft = gaussian_input(record, cpus=cpus, mem_gb=mem_gb, method=DFT_METHOD)
        hf = make_sp_input(record, "hf_631gstar", cpus=cpus, mem_gb=mem_gb)
        hf = hf.replace("fixed-geometry reference hf_631gstar", "TEMPLATE ONLY: regenerate using accepted DFT geometry before HF MK/RESP")
        prepared.append((mol, record, dft, hf))
    folder.mkdir(parents=True)
    entries = []
    for mol, record, dft, hf in prepared:
        case = folder / record["model_id"]
        case.mkdir()
        (case / "model.json").write_text(json.dumps(record, indent=2, ensure_ascii=False)+"\n")
        with Chem.SDWriter(str(case / "initial.sdf")) as writer:
            writer.write(mol)
        lines = [str(len(record["atoms"])), record["model_id"]+": MMFF94s starting geometry; NOT QM accepted"]
        lines.extend(a["element"]+" "+" ".join(f"{x:.12f}" for x in a["position"]) for a in record["atoms"])
        (case / "initial.xyz").write_text("\n".join(lines)+"\n")
        for dirname, text in (("geometry_opt", dft), ("HF_RESP_template", hf)):
            job = case / dirname
            job.mkdir()
            (job / "reference.gjf").write_text(text)
        (case / "PREPARED_ONLY.json").write_text(json.dumps({"status": "inputs_prepared_not_submitted", "DFT_geometry_accepted": False, "HF_ESP_template_requires_regeneration_from_accepted_DFT_geometry": True, "orientations_prepared": False, "RESP_fit_completed": False, "force_field_validated": False}, indent=2)+"\n")
        entries.append({k: record["metadata"][k] for k in ("model_id", "name", "smiles", "motif", "ring_size", "formula", "inchikey", "source", "reference_scope", "parameter_status")})
        entries[-1].update({"atom_count": len(record["atoms"]), "formal_charge": 0, "multiplicity": 1, "electron_count": record["electron_count"], "model_path": record["model_id"]+"/model.json"})
    catalog = {"schema_version": 1, "status": STATUS, "models": entries, "qm_submitted": False, "fitted_partial_charges_available": False, "seed": seed, "rdkit_version": rdBase.rdkitVersion, "smiles_source_access_date": "2026-11-05", "hiergo_scope": "ETH ether exists in HierGO; lactone and lactol candidates are separate literature-inspired extensions, not established HierGO blocks", "modular_partition_supported": False, "size_transfer_validated": False, "size_transfer_target_nm": [4, 10], "production_size_policy": "same validated local graph classes and charge constraints; held-out size validation, not fresh whole-sheet QM at every production size"}
    (folder / "catalog.json").write_text(json.dumps(catalog, indent=2, ensure_ascii=False)+"\n")
    paths = sorted(p for p in folder.rglob("*") if p.is_file())
    (folder / "artifact_manifest.json").write_text(json.dumps({"algorithm": "SHA-256", "files": [{"path": str(p.relative_to(folder)), "sha256": hashlib.sha256(p.read_bytes()).hexdigest(), "bytes": p.stat().st_size} for p in paths]}, indent=2)+"\n")
    return catalog


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--cpus", type=int, default=8)
    parser.add_argument("--mem-gb", type=int, default=6)
    args = parser.parse_args()
    result = export_candidates(args.output_dir, seed=args.seed, cpus=args.cpus, mem_gb=args.mem_gb)
    print(json.dumps({"models": len(result["models"]), "status": result["status"], "output_dir": args.output_dir}))
