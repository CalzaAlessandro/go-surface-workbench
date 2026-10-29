"""Audited charge-only controls; never fabricate missing coverage or renormalize.

This module does not claim complete force-field transfer (bonded/LJ parameters
and original typing rules are outside these charge-only mappings).
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CANDIDATES = ('current_gaff2', 'hiergo_opls', 'gopy_raw', 'stauffer_esp',
              'stauffer_sesp', 'stauffer_esp_hybrid', 'stauffer_sesp_hybrid',
              'old_whole_sheet', 'pinto_gaff2_aim')
BOHR_ANGSTROM = 0.529177210903
CHARGE_TOLERANCE_E = 1e-6


def _read(value):
    return json.loads(Path(value).read_text()) if isinstance(value, (str, Path)) else value


def _tables(candidate):
    """Return role -> (charge, source, locator), classification, and caveats."""
    libraries = _read(ROOT / 'references/library_charges.json')
    table = {}

    def lib(role, filename, name):
        rows = [r for r in libraries if r['file'] == filename]
        atoms = [a for r in rows for a in r['atoms'] if a['name'] == name]
        if len(atoms) != 1:
            raise ValueError(f'Ambiguous/missing library charge: {filename}:{name}')
        table[role] = (atoms[0]['charge_e'], filename, name)

    current = {
        'edge_C': ('GRT', 'CH'), 'edge_H': ('GRT', 'HC'), 'scaffold_C': ('GGG', 'CX'),
        'hydroxyl_C': ('HYD', 'CL'), 'hydroxyl_O': ('HYD', 'OL'), 'hydroxyl_H': ('HYD', 'HL'),
        'epoxide_C': ('EPX', 'CY'), 'epoxide_O': ('EPX', 'OE'),
        'carbonyl_C': ('CAB', 'CN'), 'carbonyl_O': ('CAB', 'ON'),
        'COOH_anchor': ('CBX', 'CY'), 'COOH_C': ('CBX', 'C4'),
        'COOH_O': ('CBX', 'OJ'), 'COOH_OH': ('CBX', 'OK'), 'COOH_H': ('CBX', 'HK'),
        'COO_anchor': ('CBD', 'CY'), 'COO_C': ('CBD', 'C4'), 'COO_O': ('CBD', 'OJ'),
    }
    if candidate == 'current_gaff2':
        for role, (unit, atom) in current.items():
            lib(role, f'local_snapshot/current/{unit}_gaff2.lib', atom)
        return table, 'user_block_charge_control', ['Local block-charge mapping; not a newly fitted RESP model.']
    if candidate == 'hiergo_opls':
        rows = _read(ROOT / 'references/hiergo_17_charge_rows.json')
        role_keys = {
            'hydroxyl_C': ('HYD', 'CD'), 'hydroxyl_O': ('HYD', 'OH1'), 'hydroxyl_H': ('HYD', 'HH'),
            'epoxide_C': ('EPX', 'CE'), 'epoxide_O': ('EPX', 'OE'),
            'ether_C': ('ETH', 'CT'), 'ether_O': ('ETH', 'OT'),
            'carbonyl_C': ('CAB', 'CN'), 'carbonyl_O': ('CAB', 'ON'),
            'COOH_anchor': ('CARB', 'CC'), 'COOH_C': ('CARB', 'CA'),
            'COOH_O': ('CARB', 'OC'), 'COOH_OH': ('CARB', 'OH'), 'COOH_H': ('CARB', 'HC1'),
            'edge_C': ('GRAT', 'CH'), 'edge_H': ('GRAT', 'HC'), 'scaffold_C': ('GRA', 'C')}
        for role, (residue, name) in role_keys.items():
            row = next(r for r in rows if (r['residue'], r['atom_name']) == (residue, name))
            table[role] = (row['charge_e'], 'hiergo_17_charge_rows.json', f'{residue}:{name}')
        return table, 'published_charge_table_control', ['Charge-only role mapping, not execution of the complete HierGO typing/force field.', 'COO− is not supplied by the 17-row table.']
    if candidate == 'gopy_raw' or candidate.endswith('_hybrid'):
        for role in ('edge_C', 'edge_H', 'scaffold_C'):
            unit, atom = current[role]
            lib(role, f'local_snapshot/current/{unit}_gaff2.lib', atom)
    if candidate == 'gopy_raw':
        mappings = {'hydroxyl_C': ('H1A', 'CY'), 'hydroxyl_O': ('H1A', 'OL'),
                    'hydroxyl_H': ('H1A', 'HK'), 'epoxide_C': ('E1A', 'CY'),
                    'epoxide_O': ('E1A', 'OE'), 'COOH_anchor': ('C1A', 'CY'),
                    'COOH_C': ('C1A', 'C4'), 'COOH_O': ('C1A', 'OJ'),
                    'COOH_OH': ('C1A', 'OK'), 'COOH_H': ('C1A', 'HK')}
        for role, (unit, atom) in mappings.items():
            lib(role, f'gopy/GO_tutorial/{unit}.lib', atom)
        lib('scaffold_C', 'gopy/GO_tutorial/GGG.lib', 'CX')
        return table, 'hybrid_control', ['Original GOPY motif charges retained without correction; current GRT supplies finite CH caps.',
                                         'Nonzero nominally neutral motif sums are reported and excluded from primary matched-charge ESP scoring.']
    if candidate.startswith('stauffer_'):
        source = _read(ROOT / 'references/stauffer2014_table_I_GO_charges.json')
        method = 'SESP' if 'sesp' in candidate else 'ESP'
        for prefix, group in (('hydroxyl', 'hydroxyl'), ('epoxide', 'epoxy')):
            for element, charge in source['models'][method][group].items():
                table[f'{prefix}_{element}'] = (charge, 'stauffer2014_table_I_GO_charges.json', f'Table I:{method}:{group}:{element}')
        hybrid = candidate.endswith('_hybrid')
        return table, 'hybrid_control' if hybrid else 'published_partial_charge_table', [
            'Table I indexed primary text; full PDF download unavailable. SESP is not RESP.',
            'Current GRT caps and GGG scaffold complete the basal table; this is a hybrid control.' if hybrid else
            'Only basal OH/epoxy table entries are supplied; finite caps and scaffold are unsupported.']
    if candidate == 'old_whole_sheet':
        return {}, 'descriptive_reference_only', ['OLD charges belong to their original whole sheets. No exact mapping to new fragments; no averaging transfer.']
    if candidate == 'pinto_gaff2_aim':
        return {}, 'descriptive_reference_only', ['Published periodic GAFF2/AIM DDEC6 charges and AIM LJ depend on source sheet/context. No automatic transfer onto finite fragments.']
    raise ValueError(f'Unknown candidate: {candidate}')


def assign_candidate(model_path, candidate):
    """Assign covered roles without padding or changing total charge.

    Accept a model.json path or parsed record. Uncovered atoms retain None;
    complete totals are reported only when every atom is covered. Atom IDs and
    order are preserved. The primary role overrides edge/scaffold annotations.
    """
    model = _read(model_path)
    ids = [a['id'] for a in model['atoms']]
    if len(set(ids)) != len(ids) or not ids:
        raise ValueError('Model atom identity must be nonempty and unique')
    table, kind, notes = _tables(candidate)
    atoms = []
    for atom in model['atoms']:
        role = atom['role']
        entry = table.get(role)
        charge = None if entry is None else float(entry[0])
        if charge is not None and not np.isfinite(charge):
            raise ValueError('Assigned charges must be finite')
        atoms.append({'atom_id': atom['id'], 'element': atom['element'], 'role': role,
                      'charge_e': charge, 'source': None if entry is None else entry[1],
                      'locator': None if entry is None else entry[2]})
    missing = [a for a in atoms if a['charge_e'] is None]
    total = None if missing else float(sum(a['charge_e'] for a in atoms))
    error = None if total is None else total - model['formal_charge']
    return {'model_id': model['model_id'], 'candidate': candidate, 'model_kind': kind,
            'formal_charge_e': model['formal_charge'], 'atoms': atoms,
            'coverage': {'complete': not missing, 'covered_atoms': len(atoms)-len(missing),
                         'total_atoms': len(atoms), 'unsupported_atom_ids': [a['atom_id'] for a in missing],
                         'unsupported_roles': sorted({a['role'] for a in missing})},
            'covered_charge_sum_e': float(sum(a['charge_e'] for a in atoms if a['charge_e'] is not None)),
            'total_charge_e': total, 'charge_error_e': error,
            'charge_matches_formal': error is not None and abs(error) <= CHARGE_TOLERANCE_E,
            'primary_esp_eligible': error is not None and abs(error) <= CHARGE_TOLERANCE_E,
            'notes': notes}


def compare_esp(model, esp_ref, candidates, *, chunk_size=1024,
                allow_charge_mismatch_diagnostic=False,
                charge_tolerance_e=CHARGE_TOLERANCE_E):
    """Score a shared grid/geometry in atomic units with bounded pair storage.

    esp_ref must explicitly contain atom_ids, atom_positions_angstrom,
    formal_charge, points_angstrom, potential_au and dipole_origin_angstrom.
    Candidate dictionaries allow later fitted charges, but are revalidated here.
    Grid and atom chunks bound pair arrays to chunk_size²; zero-distance points
    are rejected. No potential offset removal or charge rescaling is performed.
    """
    if (isinstance(charge_tolerance_e, bool) or not isinstance(charge_tolerance_e, (int, float))
            or not np.isfinite(charge_tolerance_e) or charge_tolerance_e <= 0):
        raise ValueError('charge_tolerance_e must be a finite positive number')
    model, ref = _read(model), _read(esp_ref)
    if not isinstance(chunk_size, int) or chunk_size < 1:
        raise ValueError('chunk_size must be a positive integer')
    ids = [a['id'] for a in model['atoms']]
    xyz = np.asarray([a['position'] for a in model['atoms']], dtype=float)
    ref_xyz = np.asarray(ref['atom_positions_angstrom'], dtype=float)
    points = np.asarray(ref['points_angstrom'], dtype=float)
    values = np.asarray(ref['potential_au'], dtype=float)
    origin = np.asarray(ref['dipole_origin_angstrom'], dtype=float)
    if not all(np.isfinite(x).all() for x in (xyz, ref_xyz, points, values, origin)):
        raise ValueError('ESP inputs and geometry must be finite')
    if not ids or len(set(ids)) != len(ids) or ids != ref['atom_ids']:
        raise ValueError('ESP atom identity/order mismatch')
    if xyz.shape != (len(ids), 3) or ref_xyz.shape != xyz.shape or not np.allclose(xyz, ref_xyz, atol=1e-8, rtol=0):
        raise ValueError('ESP geometry mismatch')
    if model['formal_charge'] != ref['formal_charge']:
        raise ValueError('ESP formal charge mismatch')
    if points.ndim != 2 or points.shape[1] != 3 or not len(points) or values.shape != (len(points),) or origin.shape != (3,):
        raise ValueError('Invalid ESP grid, potential or dipole origin shape')
    results = []
    for candidate in candidates:
        assigned = assign_candidate(model, candidate) if isinstance(candidate, str) else candidate
        result = {'candidate': assigned['candidate'], 'coverage': assigned['coverage']}
        results.append(result)
        if [a['atom_id'] for a in assigned['atoms']] != ids:
            raise ValueError('Candidate atom identity/order mismatch')
        if any(a['charge_e'] is None for a in assigned['atoms']):
            result.update(status='skipped', reason='unsupported_chemistry')
            continue
        q = np.asarray([a['charge_e'] for a in assigned['atoms']], dtype=float)
        if not np.isfinite(q).all():
            raise ValueError('Candidate charges must be finite')
        error = float(q.sum()) - model['formal_charge']
        result['charge_error_e'] = error
        matched = abs(error) <= charge_tolerance_e
        if not matched and not allow_charge_mismatch_diagnostic:
            result.update(status='skipped', reason='charge_mismatch')
            continue
        prediction = np.zeros(len(points))
        for start in range(0, len(points), chunk_size):
            block = points[start:start+chunk_size]
            for atom_start in range(0, len(xyz), chunk_size):
                distance = np.linalg.norm(block[:, None, :] - xyz[None, atom_start:atom_start+chunk_size, :], axis=2)
                if np.any(distance <= 1e-12):
                    raise ValueError('ESP grid intersects an atom')
                prediction[start:start+len(block)] += (BOHR_ANGSTROM / distance) @ q[atom_start:atom_start+chunk_size]
        residual = prediction - values
        mse = float(np.mean(residual**2))
        reference_rms = float(np.sqrt(np.mean(values**2)))
        result.update(status='matched_charge_control' if matched else 'diagnostic_charge_mismatch',
                      rmse_au=float(np.sqrt(mse)), mae_au=float(np.mean(abs(residual))),
                      bias_au=float(np.mean(residual)), max_abs_error_au=float(np.max(abs(residual))),
                      relative_rmse=None if reference_rms == 0 else float(np.sqrt(mse))/reference_rms,
                      dipole_e_angstrom=(q @ (xyz-origin)).tolist())
    return {'model_id': model['model_id'], 'n_grid_points': len(points),
            'potential_unit': 'atomic unit (e/bohr)', 'dipole_unit': 'e angstrom',
            'dipole_origin_angstrom': origin.tolist(),
            'charge_tolerance_e': float(charge_tolerance_e), 'results': results}


def export_reference_summary(root=ROOT):
    """Write auditable charge-only baseline tables; never synthesize QM values."""
    root = Path(root)
    records = [assign_candidate(path, candidate) for path in sorted((root/'models').glob('*/model.json'))
               for candidate in CANDIDATES]
    catalog = _read(root/'references/reference_catalog.json')
    urls = {r['id']: r.get('url', r.get('repository', r.get('path'))) for r in catalog['sources']}
    role_charges = []
    for candidate in CANDIDATES:
        table, _, _ = _tables(candidate)
        for role, (charge, source, locator) in sorted(table.items()):
            source_id = ('hiergo2023' if source.startswith('hiergo') else
                         'gopy2020' if source.startswith('gopy/') else
                         'stauffer2014' if source.startswith('stauffer') else None)
            role_charges.append({'candidate': candidate, 'role': role, 'charge_e': charge,
                                 'source': source, 'locator': locator,
                                 'source_url': urls.get(source_id)})
    payload = {'status': 'baseline_charge_mapping_only_no_new_QM_or_RESP_results',
               'charge_unit': 'elementary charge', 'sources': urls,
               'role_charges': role_charges, 'records': records}
    target = root/'reports/reference_charge_summary.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(payload, indent=2, allow_nan=False)+'\n')
    with target.with_suffix('.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['model_id', 'candidate', 'model_kind', 'complete',
                               'covered_atoms', 'total_atoms', 'total_charge_e', 'charge_error_e',
                               'primary_esp_eligible', 'unsupported_roles'])
        writer.writeheader()
        for r in records:
            writer.writerow({key: r[key] for key in ('model_id', 'candidate', 'model_kind', 'total_charge_e',
                            'charge_error_e', 'primary_esp_eligible')} | {
                            key: r['coverage'][key] for key in ('complete', 'covered_atoms', 'total_atoms')} | {
                            'unsupported_roles': ';'.join(r['coverage']['unsupported_roles'])})
    with target.with_name('reference_charge_summary.roles.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['candidate', 'role', 'charge_e', 'source', 'locator', 'source_url'])
        writer.writeheader()
        writer.writerows(role_charges)
    return payload


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Export honest baseline charge coverage (no QM fit).')
    parser.add_argument('--root', type=Path, default=ROOT)
    args = parser.parse_args()
    export_reference_summary(args.root)
