"""Small, explicit distance/face panel on compact finite PAH scaffolds.

Graph construction uses a perfect matching of the residual bipartite pi network.
Geometries are constrained force-field starting models, never QM evidence.
"""
from __future__ import annotations

import argparse
import copy
import itertools
import json
import math
from collections import deque
from functools import cache
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdMolDescriptors

CATEGORIES = ('OH_OH', 'EPX_EPX', 'OH_EPX')
BINS = {'near': (1, 1), 'intermediate': (2, 3), 'far': (4, 1000)}
SEED = 30102026


@cache
def _lattice(shell):
    if type(shell) is not int or shell < 2:
        raise ValueError('shell must be an integer >= 2')
    vertices, edges = set(), set()
    length = 1.42
    radius = shell - 1
    for q in range(-radius, radius + 1):
        for r in range(-radius, radius + 1):
            if max(abs(q), abs(r), abs(q + r)) > radius:
                continue
            center = (math.sqrt(3) * length * (q + r / 2), 1.5 * length * r)
            ring = [(round(center[0] + length * math.cos(math.radians(30 + 60 * i)), 6),
                     round(center[1] + length * math.sin(math.radians(30 + 60 * i)), 6))
                    for i in range(6)]
            vertices.update(ring)
            edges.update(tuple(sorted((ring[i], ring[(i + 1) % 6]))) for i in range(6))
    points = sorted(vertices, key=lambda p: (p[1], p[0]))
    lookup = {point: idx for idx, point in enumerate(points)}
    pairs = sorted(tuple(sorted((lookup[a], lookup[b]))) for a, b in edges)
    neighbors = {i: [] for i in range(len(points))}
    for a, b in pairs:
        neighbors[a].append(b)
        neighbors[b].append(a)
    colors = {0: 0}
    queue = deque([0])
    while queue:
        a = queue.popleft()
        for b in neighbors[a]:
            if b in colors and colors[b] == colors[a]:
                raise ValueError('Scaffold is not bipartite')
            if b not in colors:
                colors[b] = 1 - colors[a]
                queue.append(b)
    distances = []
    for source in range(len(points)):
        distance = {source: 0}
        queue = deque([source])
        while queue:
            a = queue.popleft()
            for b in neighbors[a]:
                if b not in distance:
                    distance[b] = distance[a] + 1
                    queue.append(b)
        distances.append(distance)
    return points, pairs, neighbors, colors, distances


def _matching(shell, removed):
    points, _, neighbors, colors, _ = _lattice(shell)
    remaining = set(range(len(points))) - set(removed)
    left = sorted(i for i in remaining if colors[i] == 0)
    if 2 * len(left) != len(remaining):
        return None
    matched = {}

    def augment(a, visited):
        for b in sorted(neighbors[a]):
            if b not in remaining or b in visited:
                continue
            visited.add(b)
            if b not in matched or augment(matched[b], visited):
                matched[b] = a
                return True
        return False

    if not all(augment(a, set()) for a in left):
        return None
    return {tuple(sorted((a, b))) for b, a in matched.items()}


def _graph(shell, motifs, allow_edge=False):
    points, pairs, neighbors, _, _ = _lattice(shell)
    removed = [identifier - 1 for motif in motifs for identifier in motif['carbon_ids']]
    if len(removed) != len(set(removed)):
        raise ValueError('Functional motifs must use distinct, nonoverlapping carbon sites')
    if any(i not in neighbors or (len(neighbors[i]) != 3 and not allow_edge) for i in removed):
        raise ValueError('Functional sites must be existing basal carbons')
    for motif in motifs:
        ids = [i - 1 for i in motif['carbon_ids']]
        if motif['type'] not in ('OH', 'EPX') or len(ids) != (1 if motif['type'] == 'OH' else 2):
            raise ValueError('OH needs one site; EPX needs two bonded sites')
        if motif['type'] == 'EPX' and tuple(sorted(ids)) not in pairs:
            raise ValueError('Epoxide sites must share an existing C-C bond')
    matching = _matching(shell, removed)
    if matching is None:
        raise ValueError('No perfect matching of remaining pi network; radical ansatz rejected')
    rw = Chem.RWMol()
    for idx in range(len(points)):
        atom = Chem.Atom('C')
        atom.SetProp('role', 'edge_C' if len(neighbors[idx]) == 2 else 'scaffold_C')
        rw.AddAtom(atom)
    for pair in pairs:
        rw.AddBond(*pair, Chem.BondType.DOUBLE if pair in matching else Chem.BondType.SINGLE)
    for motif in motifs:
        oxygen = Chem.Atom('O')
        oxygen.SetProp('role', 'hydroxyl_O' if motif['type'] == 'OH' else 'epoxide_O')
        oxygen_idx = rw.AddAtom(oxygen)
        for identifier in motif['carbon_ids']:
            rw.AddBond(identifier - 1, oxygen_idx, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(identifier - 1).SetProp(
                'role', 'hydroxyl_C' if motif['type'] == 'OH' else 'epoxide_C')
    mol = rw.GetMol()
    Chem.SanitizeMol(mol)
    # Keep the selected explicit Kekule graph; aromaticity perception can give
    # an alternate localized form, which is equivalent but obscures the gate.
    Chem.Kekulize(mol, clearAromaticFlags=True)
    return mol


def _map_hydrogens(mol):
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            parent = atom.GetNeighbors()[0]
            atom.SetProp('role', 'hydroxyl_H' if parent.GetProp('role') == 'hydroxyl_O' else 'edge_H')
        atom.SetAtomMapNum(atom.GetIdx() + 1)
    if any(atom.GetNumRadicalElectrons() for atom in mol.GetAtoms()):
        raise ValueError('Radical graph rejected')
    return mol


def build_parent(shell):
    """Hydrogen-terminated coronene (2), circumcoronene (3), or a larger shell."""
    mol = _graph(shell, [])
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for idx, point in enumerate(_lattice(shell)[0]):
        conformer.SetAtomPosition(idx, (*point, 0.0))
    mol.AddConformer(conformer)
    return _map_hydrogens(Chem.AddHs(mol, addCoords=True))


@cache
def _candidates(shell, category):
    if category not in CATEGORIES:
        raise ValueError(f'Unknown category: {category}')
    _, pairs, neighbors, colors, distances = _lattice(shell)
    basal = sorted(i for i, ns in neighbors.items() if len(ns) == 3)
    basal_edges = [pair for pair in pairs if all(i in basal for i in pair)]
    if category == 'OH_OH':
        raw = [([('OH', (a,)), ('OH', (b,))], (0, 1))
               for a, b in itertools.combinations(basal, 2) if colors[a] != colors[b]]
    elif category == 'EPX_EPX':
        raw = [([('EPX', a), ('EPX', b)], (0, 1))
               for a, b in itertools.combinations(basal_edges, 2) if not set(a) & set(b)]
    else:
        # Keep the second OH as an adjacent cis balancing partner. Its full
        # placement is reported: it is not an invisible electron-count fix.
        raw = [([('OH', (a[0],)), ('OH', (a[1],)), ('EPX', b)], (0, 2))
               for a in basal_edges for b in basal_edges if not set(a) & set(b)]
    result = []
    for groups, measured in raw:
        removed = [i for _, sites in groups for i in sites]
        if _matching(shell, removed) is None:
            continue
        first, second = [groups[i][1] for i in measured]
        distance = min(distances[a][b] for a in first for b in second)
        result.append({
            'motifs': [{'type': kind, 'carbon_ids': [i + 1 for i in sites]} for kind, sites in groups],
            'measured_pair_motif_indices': list(measured), 'graph_distance': distance,
            'sublattice_imbalance': sum(1 if colors[i] == 0 else -1 for i in removed),
            'pi_matching_complete': True,
        })
    return result


def enumerate_candidates(shell, category):
    """Return valid, nonoverlapping motif placements with original-graph distances."""
    return copy.deepcopy(_candidates(shell, category))


def _initial_geometry(mol, shell, motifs, faces):
    points = _lattice(shell)[0]
    conformer = Chem.Conformer(mol.GetNumAtoms())
    for idx, point in enumerate(points):
        conformer.SetAtomPosition(idx, (*point, 0.0))
    for offset, (motif, face) in enumerate(zip(motifs, faces, strict=True)):
        indices = [i - 1 for i in motif['carbon_ids']]
        for idx in indices:
            conformer.SetAtomPosition(idx, (*points[idx], 0.2 * face))
        center = np.mean([points[idx] for idx in indices], axis=0)
        conformer.SetAtomPosition(len(points) + offset,
                                  (*center, (1.63 if motif['type'] == 'OH' else 1.43) * face))
    mol.AddConformer(conformer)
    mol = _map_hydrogens(Chem.AddHs(mol, addCoords=True))
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        raise ValueError('MMFF parameters unavailable for initial panel geometry')
    properties = AllChem.MMFFGetMoleculeProperties(mol)
    forcefield = AllChem.MMFFGetMoleculeForceField(mol, properties)
    for idx in range(len(points)):
        forcefield.MMFFAddPositionConstraint(idx, 0.35, 100.0)
    for idx in range(len(points), len(points) + len(motifs)):
        forcefield.MMFFAddPositionConstraint(idx, 0.4, 100.0)
    forcefield.Initialize()
    status = forcefield.Minimize(maxIts=3000)
    Chem.AssignAtomChiralTagsFromStructure(mol, replaceExistingTags=True)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return mol, status


def build_configuration(shell, category, candidate, relative_face, seed=SEED):
    """Build one explicitly gated graph and constrained MMFF starting geometry."""
    if category not in (*CATEGORIES, 'EPH', 'EPT') or relative_face not in ('cis', 'trans'):
        raise ValueError('Provide a supported category and cis/trans relative face')
    if type(seed) is not int or seed < 0:
        raise ValueError('seed must be a nonnegative integer')
    motifs = copy.deepcopy(candidate['motifs'])
    expected = {'OH_OH': ['OH', 'OH'], 'EPX_EPX': ['EPX', 'EPX'],
                'OH_EPX': ['OH', 'OH', 'EPX'], 'EPH': ['EPX'], 'EPT': ['EPX']}[category]
    if [group['type'] for group in motifs] != expected:
        raise ValueError('Category does not match explicit motif counts')
    mol = _graph(shell, motifs, allow_edge=category in ('EPH', 'EPT'))
    measured = [0, 0] if category in ('EPH', 'EPT') else [0, 2] if category == 'OH_EPX' else [0, 1]
    faces = [1] * len(motifs)
    faces[measured[1]] = 1 if relative_face == 'cis' else -1
    mol, status = _initial_geometry(mol, shell, motifs, faces)
    xyz = mol.GetConformer().GetPositions()
    ncarbon = len(_lattice(shell)[0])
    for offset, (motif, face) in enumerate(zip(motifs, faces, strict=True)):
        motif.update(oxygen_id=ncarbon + offset + 1, requested_face=face)
    unmodified = [i for i in range(ncarbon) if mol.GetAtomWithIdx(i).GetProp('role')
                  not in ('hydroxyl_C', 'epoxide_C')]
    center = xyz[unmodified].mean(axis=0)
    normal = np.linalg.svd(xyz[unmodified] - center)[2][-1]
    if normal[2] < 0:
        normal *= -1
    heights = [(xyz[group['oxygen_id'] - 1] - center) @ normal for group in motifs]
    for motif, height in zip(motifs, heights, strict=True):
        motif['oxygen_plane_height_angstrom'] = float(height)
        if height * motif['requested_face'] < 0.35:
            raise ValueError('Constrained geometry did not retain requested face')
    a, b = [motifs[idx] for idx in measured]
    distance = min(_lattice(shell)[4][i - 1][j - 1] for i in a['carbon_ids'] for j in b['carbon_ids'])
    center_a = xyz[[i - 1 for i in a['carbon_ids']]].mean(axis=0)
    center_b = xyz[[i - 1 for i in b['carbon_ids']]].mean(axis=0)
    if not np.isfinite(xyz).all():
        raise ValueError('Nonfinite geometry')
    for bond in mol.GetBonds():
        length = np.linalg.norm(xyz[bond.GetBeginAtomIdx()] - xyz[bond.GetEndAtomIdx()])
        if not 0.85 < length < 1.95:
            raise ValueError(f'Implausible initial bond length {length:.5f} A')
    unlabelled = Chem.Mol(mol)
    for atom in unlabelled.GetAtoms():
        atom.SetAtomMapNum(0)
    classes = list(Chem.CanonicalRankAtoms(unlabelled, breakTies=False))
    atoms = []
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        role = atom.GetProp('role')
        roles = {role}
        if idx < ncarbon:
            roles.add('scaffold_C')
            if len(_lattice(shell)[2][idx]) == 2:
                roles.add('edge_C')
        atoms.append({'id': idx + 1, 'name': f'{atom.GetSymbol()}{idx + 1}',
                      'element': atom.GetSymbol(), 'position': xyz[idx].tolist(),
                      'role': role, 'roles': sorted(roles), 'formal_charge': atom.GetFormalCharge(),
                      'charge': None, 'symmetry_class': int(classes[idx]),
                      'chiral_tag': str(atom.GetChiralTag()),
                      'cip_code': atom.GetProp('_CIPCode') if atom.HasProp('_CIPCode') else None})
    metadata = {
        'parent_name': 'coronene' if shell == 2 else 'circumcoronene' if shell == 3 else f'hexflake_shell{shell}',
        'parent_carbon_count': ncarbon, 'parent_formula': f'C{ncarbon}H{6 * shell}',
        'scaffold_definition': 'union of hexagonal honeycomb cells of axial radius shell-1',
        'scaffold_carbon_ids': list(range(1, ncarbon + 1)), 'shell': shell,
        'category': category, 'functional_groups': motifs,
        'measured_pair_motif_indices': measured, 'relative_face': relative_face,
        'observed_relative_face': 'cis' if heights[measured[0]] * heights[measured[1]] > 0 else 'trans',
        'face_reference': 'least-squares unfunctionalized scaffold carbon plane, normal toward +z',
        'pair_graph_distance': distance,
        'graph_distance_definition': 'minimum number of original scaffold C-C edges between motif carbon sets',
        'pair_centroid_distance_angstrom': float(np.linalg.norm(center_a - center_b)),
        'pair_oxygen_distance_angstrom': float(np.linalg.norm(xyz[a['oxygen_id'] - 1] - xyz[b['oxygen_id'] - 1])),
        'coordinate_units': 'angstrom', 'geometry_status': 'initial_constrained_force_field_geometry_not_QM_optimized',
        'initial_geometry_method': 'analytic honeycomb + MMFF94 with carbon/oxygen position constraints',
        'position_constraints': {'carbon_max_displacement_angstrom': 0.35,
                                 'oxygen_max_displacement_angstrom': 0.4,
                                 'force_constant_kcal_mol_A2': 100.0},
        'optimization_status': 'converged' if status == 0 else 'iteration_limit',
        'seed': seed, 'seed_usage': 'provenance only; geometry and matching are deterministic without random sampling',
        'rdkit_version': rdBase.rdkitVersion,
        'pi_matching_complete': True, 'removed_site_sublattice_imbalance': 0,
        'multiplicity_basis': 'closed-shell singlet starting ansatz; QM wavefunction stability unverified',
        'charge_status': 'formal integer charges only; no fitted partial charges',
        'boundary': 'finite hydrogen-terminated PAH with intrinsic edge hydrogens',
        'stereo_status': 'assigned from constrained initial coordinates',
        'interpretation_limit': 'context-dependence test only; finite-size/edge/loading effects and constrained geometry need QM validation',
    }
    record = {'schema_version': 1, 'model_id': f'C{ncarbon}_{category}_{distance}_{relative_face}',
              'formula': rdMolDescriptors.CalcMolFormula(mol), 'formal_charge': Chem.GetFormalCharge(mol),
              'multiplicity': 1, 'electron_count': sum(a.GetAtomicNum() for a in mol.GetAtoms()),
              'atoms': atoms, 'bonds': [{'a': b.GetBeginAtomIdx() + 1, 'b': b.GetEndAtomIdx() + 1,
                                       'order': b.GetBondTypeAsDouble()} for b in mol.GetBonds()],
              'metadata': metadata}
    return mol, record


def build_edge_control(shell, subtype):
    """Terminal epoxide controls classified by actual attached-H counts."""
    if subtype not in ('EPH', 'EPT'):
        raise ValueError('Terminal epoxide control must be EPH or EPT')
    _, pairs, neighbors, _, _ = _lattice(shell)
    wanted = [2, 3] if subtype == 'EPH' else [2, 2]
    for pair in pairs:
        if sorted(len(neighbors[i]) for i in pair) != wanted or _matching(shell, pair) is None:
            continue
        candidate = {'motifs': [{'type': 'EPX', 'carbon_ids': [i + 1 for i in pair]}]}
        mol, record = build_configuration(shell, subtype, candidate, 'cis')
        group = record['metadata']['functional_groups'][0]
        counts = sorted(sum(n.GetAtomicNum() == 1 for n in mol.GetAtomWithIdx(i - 1).GetNeighbors())
                        for i in group['carbon_ids'])
        actual = 'EPH' if counts == [0, 1] else 'EPT' if counts == [1, 1] else 'EPX' if counts == [0, 0] else 'unsupported'
        if actual != subtype:
            raise ValueError('Terminal epoxide H subtype disagrees with explicit graph')
        for identifier in group['carbon_ids']:
            atom = mol.GetAtomWithIdx(identifier - 1)
            hs = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 1]
            role = ('EPH_CH' if hs else 'EPH_C') if subtype == 'EPH' else 'EPT_C'
            atom.SetProp('role', role)
            record['atoms'][identifier - 1]['role'] = role
            record['atoms'][identifier - 1]['roles'] = sorted({role, 'scaffold_C'} | ({'edge_C'} if hs else set()))
            for hydrogen in hs:
                hydrogen.SetProp('role', subtype + '_H')
                row = record['atoms'][hydrogen.GetIdx()]
                row.update(role=subtype + '_H', roles=[subtype + '_H'])
        oxygen = record['atoms'][group['oxygen_id'] - 1]
        oxygen.update(role=subtype + '_O', roles=[subtype + '_O'])
        mol.GetAtomWithIdx(group['oxygen_id'] - 1).SetProp('role', subtype + '_O')
        metadata = record['metadata']
        for key in list(metadata):
            if key.startswith('pair_') or key in ('measured_pair_motif_indices', 'relative_face', 'observed_relative_face'):
                del metadata[key]
        metadata.update(control_class=subtype, epoxide_hydrogen_subtype=actual,
                        attached_hydrogens_per_epoxide_carbon=counts,
                        interpretation_limit='edge subtype control; not a basal EPX replacement')
        record['model_id'] = f'C{metadata["parent_carbon_count"]}_{subtype}_edge_control'
        return mol, record
    raise ValueError(f'No matching-valid {subtype} edge pair on shell {shell}')


def _select(shell, category, distance_bin):
    candidates = _candidates(shell, category)
    lower, upper = BINS[distance_bin]
    pool = [candidate for candidate in candidates if lower <= candidate['graph_distance'] <= upper]
    if not pool:
        return None
    # A canonical placement, centered as far from the finite boundary as possible.
    points = np.array(_lattice(shell)[0])
    return copy.deepcopy(min(pool, key=lambda c: (
        c['graph_distance'],
        float(np.linalg.norm(points[[i - 1 for motif in c['motifs'] for i in motif['carbon_ids']]].mean(axis=0))),
        tuple(i for motif in c['motifs'] for i in motif['carbon_ids']))))


def export_panel(output_dir):
    """Export all feasible factorial cells; QM-stage selection is a later action."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    all_cells = list(itertools.product((2, 3), CATEGORIES, BINS, ('cis', 'trans')))
    entries, omissions, warnings = [], [], []
    for cell in all_cells:
        shell, category, distance_bin, face = cell
        candidate = _select(shell, category, distance_bin)
        if candidate is None:
            reason = 'unattainable_with_nonoverlapping_basal_sites_and_perfect_pi_matching'
        else:
            try:
                mol, record = build_configuration(shell, category, candidate, face)
            except ValueError as exc:
                reason = f'initial_geometry_rejected: {exc}'
                warnings.append({'cell': list(cell), 'reason': reason})
            else:
                metadata = record['metadata']
                metadata.update(distance_bin=distance_bin,
                                dataset_split='held_out_size' if shell == 2 else 'training_context',
                                paired_face_comparison_key=f'{shell}_{category}_{distance_bin}')
                name = f'C{metadata["parent_carbon_count"]}_{category}_{distance_bin}_{face}'
                record['model_id'] = name
                mol.SetProp('_Name', name)
                folder = output_dir / name
                folder.mkdir(parents=True, exist_ok=True)
                (folder / 'model.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
                with Chem.SDWriter(str(folder / 'initial.sdf')) as writer:
                    writer.write(mol)
                lines = [str(len(record['atoms'])), name + ': constrained initial MMFF geometry, not QM']
                lines += [a['element'] + ' ' + ' '.join(f'{v:.10f}' for v in a['position']) for a in record['atoms']]
                (folder / 'initial.xyz').write_text('\n'.join(lines) + '\n')
                entries.append({'model_id': name, 'path': name + '/model.json', 'formula': record['formula'],
                                'atom_count': len(record['atoms']), 'formal_charge': 0, 'multiplicity': 1,
                                'parent_carbon_count': metadata['parent_carbon_count'], 'category': category,
                                'distance_bin': distance_bin, 'relative_face': face,
                                'graph_distance': metadata['pair_graph_distance'],
                                'centroid_distance_angstrom': metadata['pair_centroid_distance_angstrom'],
                                'dataset_split': metadata['dataset_split']})
                continue
        omissions.append({'shell': shell, 'category': category, 'distance_bin': distance_bin,
                          'relative_face': face, 'reason': reason})
    controls = []
    from .fragments import build_fragment, model_record
    for control_class, parent_model in (('EPH', None), ('EPT', None), ('CAB', 'pyrene_dione'),
                                       ('COOH', 'pyrene_COOH'), ('COO', 'pyrene_COO')):
        if parent_model is None:
            mol, record = build_edge_control(2, control_class)
        else:
            mol = build_fragment(parent_model)
            record = model_record(parent_model, mol)
            record['model_id'] += '_edge_control'
            record['metadata'].update(control_class=control_class,
                                      control_provenance='existing validated go_param.fragments builder, unchanged')
        name = record['model_id']
        folder = output_dir / name
        folder.mkdir(parents=True, exist_ok=True)
        (folder / 'model.json').write_text(json.dumps(record, indent=2, allow_nan=False) + '\n')
        mol.SetProp('_Name', name)
        with Chem.SDWriter(str(folder / 'initial.sdf')) as writer:
            writer.write(mol)
        lines = [str(len(record['atoms'])), name + ': initial force-field edge control, not QM']
        lines += [a['element'] + ' ' + ' '.join(f'{v:.10f}' for v in a['position']) for a in record['atoms']]
        (folder / 'initial.xyz').write_text('\n'.join(lines) + '\n')
        controls.append({'model_id': name, 'path': name + '/model.json', 'control_class': control_class,
                         'formula': record['formula'], 'formal_charge': record['formal_charge'],
                         'dataset_split': 'edge_control_not_pair_distance_training'})
    warnings += [dict(cell) for cell in omissions if 'unattainable' in cell['reason']]
    catalog = {'schema_version': 1, 'models': entries, 'edge_controls': controls,
               'omitted_cells': omissions, 'warnings': warnings,
               'distance_bins_graph_edges': {key: [lo, None if hi == 1000 else hi] for key, (lo, hi) in BINS.items()},
               'design': 'all feasible factorial cells; compare within parent/category/loading/protonation; omitted cells are explicit',
               'parent_formulas': {'2': 'C24H12', '3': 'C54H18'},
               'status': 'initial_models_only_no_QM_results_or_fitted_charges'}
    (output_dir / 'catalog.json').write_text(json.dumps(catalog, indent=2, allow_nan=False) + '\n')
    return catalog


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('output_dir', type=Path)
    args = parser.parse_args()
    print(json.dumps(export_panel(args.output_dir), indent=2))
