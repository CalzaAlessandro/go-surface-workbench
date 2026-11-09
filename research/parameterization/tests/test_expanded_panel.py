"""Expanded closed-shell charge references retain requested motif counts."""
from collections import Counter

import numpy as np
import pytest
from rdkit import Chem


def test_multiple_carboxylates_are_distinct_resonance_pairs_with_correct_total():
    from go_param.expanded_panel import build_expanded_record
    from go_param.modular import partition_model

    mol,record=build_expanded_record(2,{'COO':2},distance_bin='far',seed=41)
    assert record['formal_charge']==-2
    assert record['electron_count']%2==0
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())
    assert record['metadata']['pair_graph_distances'][0]['distance']>=4
    assert record['metadata']['actual_group_counts']=={'COO':2}
    groups=[g for g in record['metadata']['functional_groups'] if g['type']=='COO']
    pairs=[{record['atoms'][i-1]['symmetry_class'] for i in g['oxygen_ids']} for g in groups]
    assert all(len(p)==1 for p in pairs) and pairs[0]!=pairs[1]
    partition_model(record)


def test_mixed_loading_and_seed_are_reproducible_without_radicals():
    from go_param.expanded_panel import build_expanded_record
    from go_param.modular import partition_model

    counts={'OH':2,'EPX':1,'COOH':1}
    _,a=build_expanded_record(3,counts,distance_bin='near',seed=9)
    _,b=build_expanded_record(3,counts,distance_bin='near',seed=9)
    assert a['metadata']['actual_group_counts']==counts
    assert a['formal_charge']==0
    assert a['metadata']['canonical_isomeric_SMILES']==b['metadata']['canonical_isomeric_SMILES']
    np.testing.assert_allclose([x['position'] for x in a['atoms']],[x['position'] for x in b['atoms']],atol=1e-10)
    assert Counter(x['role'] for x in a['atoms'])['hydroxyl_O']==2
    partition_model(a)


def test_impossible_counts_fail_without_silently_reducing_loading():
    from go_param.expanded_panel import build_expanded_record

    with pytest.raises(ValueError,match='place|sites|count'):
        build_expanded_record(2,{'OH':99},seed=1)
