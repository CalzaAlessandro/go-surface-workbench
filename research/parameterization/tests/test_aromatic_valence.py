"""Aromatic representation must admit the original atom-wise bond valences."""
import json
from pathlib import Path

import pytest

from go_param import alternatives as a


def test_native_c24_partial_aromatic_output_is_a_valid_kekule_representation():
    data=json.loads((Path(__file__).parent/'data/c24_single_epx_native_bond_representation.json').read_text())
    bonds={tuple(pair):value for pair,value in data['actual_bonds']}
    report=a.validate_bond_representation(data['model'],bonds,data['output_types'])
    assert report['fixed_sigma_and_functional_group_bonds_preserved']


def benzene():
    model={'atoms':[{'id':i,'element':'C','formal_charge':0} for i in range(1,7)],
           'bonds':[{'a':i,'b':i+1 if i<6 else 1,'order':2 if i%2 else 1} for i in range(1,7)]}
    model['metadata']={'am1_bcc_resonance':{'graph_sha256':a.graph_digest(model),
                        'aromatic_bond_pairs':[sorted([b['a'],b['b']]) for b in model['bonds']]}}
    return model


def test_aromatic_metadata_does_not_allow_changing_every_ring_bond_to_double():
    model=benzene()
    actual={tuple(sorted((b['a'],b['b']))):'2' for b in model['bonds']}
    with pytest.raises(ValueError,match='valence'):
        a.validate_bond_representation(model,actual,[])


def test_full_aromatic_encoding_has_an_integer_valence_preserving_assignment():
    model=benzene()
    actual={tuple(sorted((b['a'],b['b']))):'ar' for b in model['bonds']}
    assert a.validate_bond_representation(model,actual,[])['validated_aromatic_metadata']
