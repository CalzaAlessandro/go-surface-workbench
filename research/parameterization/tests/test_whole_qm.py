"""Whole-sheet QM input and output gates; a normal exit is not a spin proof."""
from pathlib import Path

import pytest


def test_electron_parity_and_finite_whole_sheet_input():
    from go_param.whole_qm import make_orca_input

    model={'name':'odd','formal_charge':-1,'atoms':[
        {'element':'C','atomic_number':6,'position':[0,0,0]}]}
    with pytest.raises(ValueError,match='parity'):
        make_orca_input(model,multiplicity=1,processes=64,maxcore_mb=3000)
    text=make_orca_input(model,multiplicity=2,processes=64,maxcore_mb=3000)
    assert 'ma-def2-SVP' in text and 'CPCM(Water)' in text
    assert 'nprocs 64' in text and '* xyz -1 2' in text
    assert 'Opt' not in text and 'HF' not in text


def test_real_orca_probe_requires_convergence_energy_and_atom_mapped_populations():
    from go_param.whole_qm import inspect_orca

    text=(Path(__file__).parent/'data/orca_pbe_water_runtime.out').read_text()
    result=inspect_orca(text,['O','H','H'],formal_charge=0,multiplicity=1)
    assert result['accepted_single_point_diagnostic']
    assert not result['ground_spin_state_established']
    assert len(result['hirshfeld_atoms'])==3
    bad=text.replace('SCF CONVERGED','SCF UNCONVERGED')
    assert not inspect_orca(bad,['O','H','H'],formal_charge=0,multiplicity=1)['accepted_single_point_diagnostic']
    assert not inspect_orca(text,['C','H','H'],formal_charge=0,multiplicity=1)['accepted_single_point_diagnostic']
