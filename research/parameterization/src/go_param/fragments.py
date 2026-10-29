"""Explicit closed-shell finite PAH models for a prospective GO QM campaign.

All geometries are force-field starting structures, never QM results. Atom maps
are stable within each model. The hydrogen-terminated pyrene boundary is an
intrinsic finite-flake boundary, not a representation of bulk graphene.
"""
from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

PYRENE_SMILES = 'c1cc2ccc3cccc4ccc(c1)c2c34'
MODEL_NAMES = (
    'pyrene_pristine', 'pyrene_epoxide', 'pyrene_diol_cis',
    'pyrene_diol_trans', 'pyrene_dione', 'pyrene_COOH', 'pyrene_COO',
    'pyrene_diol_epoxide',
)
SEED = 29102026


def _set_role(atom, role):
    atom.SetProp('role', role)


def _base(name):
    mol = Chem.MolFromSmiles(PYRENE_SMILES)
    if name == 'pyrene_diol_epoxide':
        # A different Kekule matching provides two disjoint basal C=C bonds.
        candidates = Chem.ResonanceMolSupplier(mol, Chem.KEKULE_ALL)
        for candidate in candidates:
            basal = [b for b in candidate.GetBonds()
                     if b.GetBondType() == Chem.BondType.DOUBLE
                     and all(a.GetDegree() == 3 for a in (b.GetBeginAtom(), b.GetEndAtom()))]
            if len(basal) >= 2 and not ({basal[0].GetBeginAtomIdx(), basal[0].GetEndAtomIdx()}
                                     & {basal[1].GetBeginAtomIdx(), basal[1].GetEndAtomIdx()}):
                mol = Chem.Mol(candidate)
                for atom in mol.GetAtoms():
                    atom.SetIsAromatic(False)
                for bond in mol.GetBonds():
                    bond.SetIsAromatic(False)
                break
        else:
            raise ValueError('No disjoint basal double bonds in pyrene')
    else:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    for atom in mol.GetAtoms():
        _set_role(atom, 'edge_C' if atom.GetDegree() == 2 else 'scaffold_C')
    return Chem.RWMol(mol)


def _oxygen(rw, role, charge=0):
    atom = Chem.Atom('O')
    atom.SetFormalCharge(charge)
    _set_role(atom, role)
    return rw.AddAtom(atom)


def _basal_pairs(rw):
    return [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rw.GetBonds()
            if b.GetBondType() == Chem.BondType.DOUBLE
            and all(a.GetDegree() == 3 and a.GetSymbol() == 'C'
                    for a in (b.GetBeginAtom(), b.GetEndAtom()))]


def _add_basal_group(rw, pair, kind):
    a, b = pair
    rw.GetBondBetweenAtoms(a, b).SetBondType(Chem.BondType.SINGLE)
    if kind == 'epoxide':
        oxygen = _oxygen(rw, 'epoxide_O')
        for carbon in pair:
            rw.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
            _set_role(rw.GetAtomWithIdx(carbon), 'epoxide_C')
    else:
        for carbon in pair:
            oxygen = _oxygen(rw, 'hydroxyl_O')
            rw.AddBond(carbon, oxygen, Chem.BondType.SINGLE)
            _set_role(rw.GetAtomWithIdx(carbon), 'hydroxyl_C')


def _graph(name):
    rw = _base(name)
    if name in ('pyrene_epoxide', 'pyrene_diol_cis', 'pyrene_diol_trans', 'pyrene_diol_epoxide'):
        pairs = _basal_pairs(rw)
        if name == 'pyrene_diol_epoxide':
            _add_basal_group(rw, pairs[0], 'diol')
            _add_basal_group(rw, pairs[1], 'epoxide')
        else:
            _add_basal_group(rw, pairs[0], 'epoxide' if name == 'pyrene_epoxide' else 'diol')
    elif name == 'pyrene_dione':
        pair = next((b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in rw.GetBonds()
                    if b.GetBondType() == Chem.BondType.DOUBLE
                    and all(a.GetDegree() == 2 for a in (b.GetBeginAtom(), b.GetEndAtom())))
        rw.GetBondBetweenAtoms(*pair).SetBondType(Chem.BondType.SINGLE)
        for carbon in pair:
            oxygen = _oxygen(rw, 'carbonyl_O')
            rw.AddBond(carbon, oxygen, Chem.BondType.DOUBLE)
            _set_role(rw.GetAtomWithIdx(carbon), 'carbonyl_C')
    elif name in ('pyrene_COOH', 'pyrene_COO'):
        prefix = 'COOH' if name == 'pyrene_COOH' else 'COO'
        anchor = next(a.GetIdx() for a in rw.GetAtoms() if a.GetDegree() == 2)
        _set_role(rw.GetAtomWithIdx(anchor), prefix + '_anchor')
        carbon = Chem.Atom('C')
        _set_role(carbon, prefix + '_C')
        carbon_idx = rw.AddAtom(carbon)
        rw.AddBond(anchor, carbon_idx, Chem.BondType.SINGLE)
        carbonyl = _oxygen(rw, prefix + '_O')
        hydroxyl = _oxygen(rw, 'COOH_OH' if prefix == 'COOH' else 'COO_O', 0 if prefix == 'COOH' else -1)
        rw.AddBond(carbon_idx, carbonyl, Chem.BondType.DOUBLE)
        rw.AddBond(carbon_idx, hydroxyl, Chem.BondType.SINGLE)
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    mol = Chem.AddHs(mol)
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'H':
            parent_role = atom.GetNeighbors()[0].GetProp('role')
            _set_role(atom, {'hydroxyl_O': 'hydroxyl_H', 'COOH_OH': 'COOH_H'}.get(parent_role, 'edge_H'))
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    if any(a.GetNumRadicalElectrons() for a in mol.GetAtoms()):
        raise ValueError(f'Radical fragment rejected: {name}')
    return mol


def _diol_face(mol):
    xyz = mol.GetConformer().GetPositions()
    scaffold = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'C'
                and a.GetProp('role') not in ('hydroxyl_C', 'epoxide_C')]
    oxygens = [a.GetIdx() for a in mol.GetAtoms() if a.GetProp('role') == 'hydroxyl_O']
    center = xyz[scaffold].mean(axis=0)
    normal = np.linalg.svd(xyz[scaffold] - center)[2][-1]
    heights = (xyz[oxygens] - center) @ normal
    return float(np.prod(heights)), heights.tolist()


def _embed(mol, seed):
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(mol, params) != 0:
        return False
    if AllChem.MMFFHasAllMoleculeParams(mol):
        status = AllChem.MMFFOptimizeMolecule(mol, maxIters=3000)
        method = 'MMFF94'
    elif AllChem.UFFHasAllMoleculeParams(mol):
        status = AllChem.UFFOptimizeMolecule(mol, maxIters=3000)
        method = 'UFF'
    else:
        raise ValueError('No supported initial force field')
    if status != 0:
        return False
    mol.SetProp('initial_geometry_method', method)
    mol.SetIntProp('embedding_seed', seed)
    return True


def build_fragment(name):
    """Return a validated molecule with explicit H, stable maps and one conformer."""
    if name not in MODEL_NAMES:
        raise ValueError(f'Unknown model: {name}')
    graph = _graph(name)
    desired = -1 if name == 'pyrene_diol_trans' else 1
    hydroxyl_c = [a.GetIdx() for a in graph.GetAtoms() if a.GetProp('role') == 'hydroxyl_C']
    tag_pairs = itertools.product((Chem.ChiralType.CHI_TETRAHEDRAL_CW,
                                   Chem.ChiralType.CHI_TETRAHEDRAL_CCW), repeat=len(hydroxyl_c))
    for attempt, tags in enumerate(tag_pairs):
        mol = Chem.Mol(graph)
        for idx, tag in zip(hydroxyl_c, tags):
            mol.GetAtomWithIdx(idx).SetChiralTag(tag)
        if not _embed(mol, SEED + attempt):
            continue
        if hydroxyl_c:
            product, heights = _diol_face(mol)
            if np.sign(product) != desired or min(abs(h) for h in heights) < 0.3:
                continue
        # Capture actual tetrahedral stereochemistry of the generated geometry.
        Chem.AssignAtomChiralTagsFromStructure(mol, confId=0, replaceExistingTags=True)
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        xyz = mol.GetConformer().GetPositions()
        if not np.isfinite(xyz).all():
            raise ValueError(f'Nonfinite coordinates: {name}')
        for bond in mol.GetBonds():
            length = np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()])
            if not 0.85 < length < 1.9:
                raise ValueError(f'Invalid initial bond length: {name}, {length}')
        mol.SetProp('_Name', name)
        return mol
    raise ValueError(f'Unable to construct requested stereochemistry: {name}')


def model_record(name, mol):
    """Serialize the explicit QM graph; symmetry classes are graph equivalences."""
    unlabelled = Chem.Mol(mol)
    for atom in unlabelled.GetAtoms():
        atom.SetAtomMapNum(0)
    classes = list(Chem.CanonicalRankAtoms(unlabelled, breakTies=False))
    equivalent_o = [a.GetIdx() for a in mol.GetAtoms() if a.GetProp('role') == 'COO_O']
    if equivalent_o:
        # Localized integer valence graph: force only the chemically justified
        # resonance oxygen pair equivalent for later constrained charge fitting.
        shared = max(classes) + 1
        for idx in equivalent_o:
            classes[idx] = shared
    xyz = mol.GetConformer().GetPositions()
    atoms = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        role = atom.GetProp('role')
        roles = [role]
        if atom.GetSymbol() == 'C' and idx < 16:
            roles.append('scaffold_C')
            carbon_degree = sum(n.GetSymbol() == 'C' and n.GetIdx() < 16 for n in atom.GetNeighbors())
            if carbon_degree == 2:
                roles.append('edge_C')
        atoms.append({'id': idx + 1, 'name': f'{atom.GetSymbol()}{idx + 1}',
                      'element': atom.GetSymbol(), 'position': xyz[idx].tolist(),
                      'role': role, 'roles': sorted(set(roles)),
                      'formal_charge': atom.GetFormalCharge(), 'symmetry_class': int(classes[idx]),
                      'chiral_tag': str(atom.GetChiralTag()),
                      'cip_code': atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None})
    metadata = {
        'parent_name': 'pyrene', 'parent_smiles': PYRENE_SMILES,
        'boundary': 'intrinsic finite hydrogen-terminated PAH; edge H are standard caps',
        'interpretation_limit': 'small finite PAH model; size, edge and loading effects require held-out validation',
        'geometry_status': 'initial_force_field_geometry_not_QM_optimized',
        'initial_geometry_method': mol.GetProp('initial_geometry_method'),
        'seed': mol.GetIntProp('embedding_seed'), 'rdkit_version': rdBase.rdkitVersion,
        'coordinate_units': 'angstrom', 'charge_status': 'formal integer charges only; no fitted partial charges',
        'carboxylate_resonance_equivalence': [idx + 1 for idx in equivalent_o],
        'symmetry_class_basis': 'unmapped graph automorphism; COO oxygen resonance pair additionally tied',
        'stereo_status': 'assigned from explicit 3D geometry',
        'multiplicity_basis': 'closed-shell singlet starting ansatz; QM wavefunction stability unverified',
        'optimization_status': 'initial force-field minimization converged',
    }
    if any(a['role'] == 'hydroxyl_C' for a in atoms):
        product, heights = _diol_face(mol)
        metadata.update({'diol_face': 'cis' if product > 0 else 'trans',
                         'hydroxyl_O_plane_heights_angstrom': heights,
                         'face_reference': 'least squares plane of unfunctionalized scaffold carbon atoms'})
    return {'schema_version': 1, 'model_id': name,
            'formula': rdMolDescriptors.CalcMolFormula(mol),
            'formal_charge': Chem.GetFormalCharge(mol), 'multiplicity': 1,
            'electron_count': sum(a.GetAtomicNum() for a in mol.GetAtoms()) - Chem.GetFormalCharge(mol),
            'atoms': atoms,
            'bonds': [{'a': b.GetBeginAtomIdx() + 1, 'b': b.GetEndAtomIdx() + 1,
                       'order': b.GetBondTypeAsDouble()} for b in mol.GetBonds()],
            'metadata': metadata}


def export_models(output_dir, names=None):
    """Write model.json, SDF and XYZ files plus a machine-readable catalog."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for name in MODEL_NAMES if names is None else names:
        mol = build_fragment(name)
        record = model_record(name, mol)
        folder = output_dir / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'model.json').write_text(json.dumps(record, indent=2) + '\n')
        with Chem.SDWriter(str(folder / 'initial.sdf')) as writer:
            writer.write(mol)
        lines = [str(len(record['atoms'])), f'{name}: initial force-field geometry; charge={record["formal_charge"]}; multiplicity=1']
        lines.extend(f'{a["element"]:2s} ' + ' '.join(f'{x:.10f}' for x in a['position']) for a in record['atoms'])
        (folder / 'initial.xyz').write_text('\n'.join(lines) + '\n')
        records.append({'model_id': name, 'path': name + '/model.json',
                        'formula': record['formula'], 'formal_charge': record['formal_charge'],
                        'multiplicity': 1, 'atom_count': len(record['atoms'])})
    catalog = {'schema_version': 1, 'models': records,
               'status': 'initial_models_only_no_QM_results_or_fitted_parameters'}
    (output_dir / 'catalog.json').write_text(json.dumps(catalog, indent=2) + '\n')
    return catalog


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir', type=Path)
    parser.add_argument('--models', nargs='+', choices=MODEL_NAMES)
    args = parser.parse_args()
    print(json.dumps(export_models(args.output_dir, args.models), indent=2))
