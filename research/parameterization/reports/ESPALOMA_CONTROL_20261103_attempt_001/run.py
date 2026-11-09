"""Bounded EspalomaCharge benchmark; original libraries/QM data are read-only."""
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import csv
import os
from pathlib import Path
import random
import statistics
import sys
import time
import traceback

HERE = Path(__file__).resolve().parent

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save(path, data):
    with Path(path).open('x') as handle:
        json.dump(data, handle, indent=2, allow_nan=False)
        handle.write('\n')

def table(path, rows):
    with Path(path).open('x', newline='') as handle:
        w = csv.DictWriter(handle, fieldnames=list(rows[0]))
        w.writeheader(); w.writerows(rows)

def molecule(model):
    from rdkit import Chem
    from rdkit.Geometry import Point3D
    rw = Chem.RWMol(); index = {}
    for atom in model['atoms']:
        a = Chem.Atom(atom['element'])
        a.SetFormalCharge(atom['formal_charge']); a.SetNoImplicit(True)
        a.SetIntProp('source_atom_id', atom['id'])
        index[atom['id']] = rw.AddAtom(a)
    orders = {1: Chem.BondType.SINGLE, 2: Chem.BondType.DOUBLE,
              3: Chem.BondType.TRIPLE, 1.5: Chem.BondType.AROMATIC}
    for bond in model['bonds']:
        rw.AddBond(index[bond['a']], index[bond['b']], orders[bond['order']])
    mol = rw.GetMol(); Chem.SanitizeMol(mol)
    assert mol.GetNumAtoms() == len(model['atoms'])
    assert Counter(a.GetSymbol() for a in mol.GetAtoms()) == Counter(a['element'] for a in model['atoms'])
    assert Chem.GetFormalCharge(mol) == model['formal_charge']
    assert sum(a.GetNumImplicitHs() for a in mol.GetAtoms()) == 0
    assert sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()) == 0
    assert {tuple(sorted((b.GetBeginAtomIdx(), b.GetEndAtomIdx()))) for b in mol.GetBonds()} == {
        tuple(sorted((index[b['a']], index[b['b']]))) for b in model['bonds']}
    conf = Chem.Conformer(mol.GetNumAtoms()); conf.Set3D(True)
    for atom in model['atoms']:
        conf.SetAtomPosition(index[atom['id']], Point3D(*atom['position']))
    mol.AddConformer(conf)
    Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    stereo_discrepancies = []
    for atom in model['atoms']:
        if atom.get('cip_code') in ('R', 'S'):
            a = mol.GetAtomWithIdx(index[atom['id']])
            actual = a.GetProp('_CIPCode') if a.HasProp('_CIPCode') else None
            if actual != atom['cip_code']:
                stereo_discrepancies.append({'atom_id': atom['id'], 'stored_CIP': atom['cip_code'],
                                             'coordinate_derived_CIP': actual})
    # Upstream features contain no stereo tags or coordinates. Retain these
    # discrepancies as a separate source-QC issue, without claiming a stereo
    # validation or changing the reference graph, coordinates or charge state.
    mol.SetProp('stereo_metadata_discrepancies', json.dumps(stereo_discrepancies))
    assert [a.GetIntProp('source_atom_id') for a in mol.GetAtoms()] == [a['id'] for a in model['atoms']]
    return mol

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--validate-input-only', action='store_true')
    parser.add_argument('--backend', choices=['dgl', 'pyg'], default='dgl')
    args = parser.parse_args()
    source = json.loads((HERE/'inputs.json').read_text())
    from rdkit import Chem
    prepared = [(case, molecule(case['model'])) for case in source['cases']]
    assert len(prepared) == 16
    if args.validate_input_only:
        print(json.dumps({'graphs_passing_atom_order_composition_charge_valence_connectivity': len(prepared),
                          'stereo_metadata_discrepancies': {case['model']['model_id']: json.loads(mol.GetProp('stereo_metadata_discrepancies'))
                              for case, mol in prepared if json.loads(mol.GetProp('stereo_metadata_discrepancies'))},
                          'stereochemical_validation_claimed': False}, indent=2))
        return
    import numpy as np
    import torch
    import espaloma_charge.app as app
    import espaloma_charge.models
    from espaloma_charge.utils import from_rdkit_mol
    from urllib.request import urlretrieve
    torch.set_num_threads(2); torch.set_num_interop_threads(1)
    torch.manual_seed(0); torch.use_deterministic_algorithms(True)
    output = HERE/'results'; output.mkdir(exist_ok=False)
    checkpoint = HERE/('model_release.pt' if args.backend == 'dgl' else 'model.pt')
    url = ('https://github.com/choderalab/espaloma-charge/releases/download/v0.0.8/model.pt'
           if args.backend == 'dgl' else app.MODEL_URL.strip())
    if not checkpoint.exists(): urlretrieve(url, checkpoint)
    if args.backend == 'dgl':
        os.chdir(HERE)
        cache = HERE/'.model.pt'
        if not cache.exists(): cache.write_bytes(checkpoint.read_bytes())
        assert digest(cache) == digest(checkpoint)
    else:
        app.MODEL_PATH = str(checkpoint)
    model = torch.load(checkpoint, map_location='cpu', weights_only=False)
    original_training = model.training
    model.eval()
    def predict(mol):
        graph = from_rdkit_mol(mol)
        if args.backend == 'dgl':
            assert graph.num_nodes() == mol.GetNumAtoms()
            # Use the published readout's q_ref hook to supply molecular Q.
            # The v0.0.8 featurizer otherwise provides only formal_charge.
            graph.ndata['q_ref'] = graph.ndata['formal_charge'].clone()
        else:
            assert graph.num_nodes == mol.GetNumAtoms()
        with torch.inference_mode():
            result = model(graph)
            values = result.ndata['q'] if args.backend == 'dgl' else result.q
            return values.detach().cpu().numpy().reshape(-1).astype(float)
    probes = []
    for smiles in ['CO', 'CC(=O)[O-]', 'N#N']:
        mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
        q = predict(mol); api = app.charge(mol)
        assert len(q) == mol.GetNumAtoms() and np.isfinite(q).all()
        assert abs(q.sum()-Chem.GetFormalCharge(mol)) <= 5e-5
        probes.append({'smiles': smiles, 'formal_charge': Chem.GetFormalCharge(mol),
                       'sum_charge_e': float(q.sum()), 'unadapted_API_sum_charge_e': float(api.sum()),
                       'unadapted_API_vs_benchmark_max_difference_e': float(np.max(abs(q-api)))})
    save(output/'CONTROL_PROBES.json', probes)
    records, atoms, roles, failures = [], [], [], []
    for case, mol in prepared:
        info = case['model']; mid = info['model_id']
        try:
            started = time.perf_counter(); q = predict(mol); elapsed = time.perf_counter()-started
            assert len(q) == len(info['atoms']) and np.isfinite(q).all()
            assert abs(q.sum()-info['formal_charge']) <= 5e-5
            repeated = predict(mol)
            repeat_difference = float(np.max(abs(q-repeated)))
            assert repeat_difference <= 5e-6
            order = list(range(len(q))); random.Random(97).shuffle(order)
            permuted = predict(Chem.RenumberAtoms(mol, order))
            permutation_difference = float(np.max(abs(q[np.array(order)]-permuted)))
            assert permutation_difference <= 5e-5
            refs = {key: np.array(value, dtype=float) if value is not None else None
                    for key, value in case['references'].items()}
            values = dict(espaloma=q, **refs)
            comparisons = {}
            for a, b in [('espaloma', 'native_AM1_BCC'), ('espaloma', 'current_library'),
                         ('espaloma', 'molecular_HF_RESP'), ('native_AM1_BCC', 'current_library'),
                         ('native_AM1_BCC', 'molecular_HF_RESP'), ('current_library', 'molecular_HF_RESP')]:
                if values[a] is not None and values[b] is not None:
                    delta = values[a]-values[b]
                    comparisons[a+'__'+b] = {'charge_RMS_difference_e': float(np.sqrt(np.mean(delta**2))),
                        'charge_MAE_e': float(np.mean(abs(delta))), 'max_absolute_difference_e': float(np.max(abs(delta)))}
            row = {'model_id': mid, 'status': 'technical_gates_passed', 'atoms': len(q),
                   'formal_charge': info['formal_charge'], 'Espaloma_charge_sum_e': float(q.sum()),
                   'repeat_max_difference_e': repeat_difference, 'permutation_max_difference_e': permutation_difference,
                   'inference_seconds': elapsed, 'native_BCC_primary_available': case['native_BCC_primary_available'],
                   'native_BCC_exclusion': case['native_BCC_exclusion'], 'comparisons': comparisons,
                   'canonical_isomeric_SMILES': Chem.MolToSmiles(mol, isomericSmiles=True),
                   'stereo_metadata_discrepancies': json.loads(mol.GetProp('stereo_metadata_discrepancies')),
                   'stereochemical_validation_claimed': False,
                   'final_parameter_validation': False}
            records.append(row)
            case_rows = []
            for i, atom in enumerate(info['atoms']):
                entry = {'model_id': mid, 'atom_id': atom['id'], 'element': atom['element'], 'role': atom['role'],
                         'espaloma_charge_e': float(q[i]),
                         **{key+'_charge_e': None if val is None else float(val[i]) for key, val in refs.items()}}
                atoms.append(entry); case_rows.append(entry)
            for role in dict.fromkeys(a['role'] for a in info['atoms']):
                selected = [r for r in case_rows if r['role'] == role]
                roles.append({'model_id': mid, 'role': role, 'atoms': len(selected),
                              **{method+'_mean_charge_e': None if selected[0][method+'_charge_e'] is None else
                                 statistics.mean(r[method+'_charge_e'] for r in selected)
                                 for method in ['espaloma', *refs]}})
            save(output/(mid+'.json'), {'summary': row, 'atoms': case_rows, 'source_SHA256': case['source_SHA256']})
            print(mid, 'passed', flush=True)
        except Exception as exc:
            failure = {'model_id': mid, 'error': repr(exc), 'traceback': traceback.format_exc()}
            failures.append(failure); save(output/(mid+'.FAILED.json'), failure)
            print(mid, 'FAILED', repr(exc), flush=True)
    if atoms:
        table(output/'per_atom_charge_comparison.csv', atoms)
        table(output/'per_context_role_comparison.csv', roles)
    matched = [r for r in records if r['native_BCC_primary_available']]
    aggregate = {}
    for key in ['espaloma__native_AM1_BCC', 'espaloma__current_library', 'native_AM1_BCC__current_library',
                'espaloma__molecular_HF_RESP', 'native_AM1_BCC__molecular_HF_RESP', 'current_library__molecular_HF_RESP']:
        vals = [r['comparisons'][key]['charge_RMS_difference_e'] for r in matched]
        if vals: aggregate[key] = {'common_contexts': len(vals), 'median_context_charge_RMS_difference_e': statistics.median(vals),
                                   'min_context_RMS_difference_e': min(vals), 'max_context_RMS_difference_e': max(vals)}
    packages = {name: importlib.metadata.version(name) for name in ['torch','dgl' if args.backend == 'dgl' else 'torch-geometric','rdkit','numpy','espaloma_charge']}
    summary = {'status': 'BENCHMARK_COMPLETE' if len(records) == 16 else 'BENCHMARK_PARTIAL_WITH_FAILURES',
        'completed_utc': datetime.now(timezone.utc).isoformat(), 'successful_contexts': len(records),
        'failed_contexts': len(failures), 'matched_BCC_contexts': len(matched), 'aggregate_common_subset': aggregate,
        'model_source_commit': 'b32649c8df2029fbd25dd18716fe838e0c62d6ae' if args.backend == 'dgl' else source['source_commit'],
        'initial_development_commit_not_run': source['source_commit'],
        'backend': args.backend, 'checkpoint_URL': url, 'checkpoint_SHA256': digest(checkpoint),
        'package_versions': packages, 'Python': sys.version,
        'upstream_serialized_training_mode': original_training, 'benchmark_explicit_eval_mode': True,
        'charge_adapter': 'formal_charge to supported q_ref field' if args.backend == 'dgl' else 'native q_ref',
        'API_probe_agreement_recorded': True, 'graph_only_input': True,
        'upstream_features_use_stereochemistry_or_coordinates': False,
        'contexts_with_stereo_metadata_discrepancies': sum(bool(r['stereo_metadata_discrepancies']) for r in records),
        'stereochemical_validation_claimed': False,
        'weights_refitted': False, 'new_Gaussian_or_RESP_fit': False,
        'charge_similarity_is_physical_validation': False, 'final_parameter_validation': False,
        'inputs_SHA256': digest(HERE/'inputs.json'), 'worker_SHA256': digest(__file__),
        'artifacts_SHA256': {str(p.relative_to(output)): digest(p) for p in output.rglob('*') if p.is_file()}}
    save(output/'BENCHMARK_COMPLETE.json', summary)
    print(json.dumps(summary, indent=2), flush=True)

if __name__ == '__main__':
    main()
