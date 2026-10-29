"""AM1-BCC transport/safety fixtures; no test is an actual charge calculation."""
import importlib
import json
from pathlib import Path
import subprocess

import pytest


def api():
    return importlib.import_module('go_param.alternatives')


@pytest.fixture
def model(tmp_path):
    m = {'model_id': 'explicit_water_fixture', 'formal_charge': 0, 'multiplicity': 1,
         'metadata': {'geometry_status': 'test_fixture_not_QM'},
         'atoms': [{'id': 1, 'element': 'O', 'position': [0., 0., 0.], 'formal_charge': 0, 'role': 'test'},
                   {'id': 2, 'element': 'H', 'position': [.95, 0., 0.], 'formal_charge': 0, 'role': 'test'},
                   {'id': 3, 'element': 'H', 'position': [-.3, .9, 0.], 'formal_charge': 0, 'role': 'test'}],
         'bonds': [{'a': 1, 'b': 2, 'order': 1}, {'a': 1, 'b': 3, 'order': 1}]}
    p = tmp_path/'model.json'; p.write_text(json.dumps(m))
    return p, m


def mol2(model, q=None):
    q = q or [-.8, .4, .4]
    return ('@<TRIPOS>MOLECULE\nGOF\n3 2 1 0 0\nSMALL\nUSER_CHARGES\n'
            '@<TRIPOS>ATOM\n' + ''.join(f"{a['id']} {a['element']}{a['id']:03d} {a['position'][0]} {a['position'][1]} {a['position'][2]} {a['element'].lower()} 1 GOF {v}\n" for a,v in zip(model['atoms'],q))
            + '@<TRIPOS>BOND\n1 1 2 1\n2 1 3 1\n@<TRIPOS>SUBSTRUCTURE\n1 GOF 1\n')


def fixtures(monkeypatch, model, *, q=None, failure=None):
    p, m = model
    calls = []
    params=p.parent/'TEST_BCCPARM.DAT'; params.write_text('1 31 91 1 0.0000\n')
    monkeypatch.setattr(api(), '_bcc_parameter_file', lambda programs: params)
    monkeypatch.setattr(api(), '_program', lambda name: {'path': str(name), 'sha256': 'TEST_STUB_EXECUTABLE'})
    def run(args, **kw):
        calls.append(args)
        out = Path(kw['cwd'])
        if failure == 'exit':
            return subprocess.CompletedProcess(args, 1)
        text = mol2(m, q)
        if failure == 'mapping':
            text = text.replace('1 O001', '1 H002', 1)
        target = out/args[args.index('-o')+1]
        target.write_text(text)
        (out/'sqm.in').write_text(" &qmmm\n qm_theory='AM1', maxcyc=9999, qmcharge=0\n /\n")
        values=q or [-.8,.4,.4]
        for fname in ('ANTECHAMBER_AM1BCC_PRE.AC','ANTECHAMBER_AM1BCC.AC'):
            (out/fname).write_text(''.join(f"ATOM{a['id']:7d}  {a['element']}{a['id']:03d}GOF 1 0 0 0 {v:.6f} {31 if a['element']=='O' else 91}\n" for a,v in zip(m['atoms'],values))+'BOND 1 1 2 1\nBOND 2 1 3 1\n')
        (out/'sqm.out').write_text(('SCF failed to converge\n' if failure == 'sqm' else ' ... geometry converged !\n')
                                  + ' Atom Element Mulliken Charge\n'+''.join(f" {a['id']} {a['element']} {v:.3f}\n" for a,v in zip(m['atoms'],values))+' Total Mulliken Charge = 0.000\n'
                                  + ' Final Structure\n QMMM: QM_NO. MM_NO. ATOM X Y Z\n'
                                  + ''.join(f" QMMM: {a['id']} {a['id']} {a['element']} {' '.join(map(str,a['position']))}\n" for a in m['atoms'])
                                  + ' --------- Calculation Completed ----------\n')
        return subprocess.CompletedProcess(args, 0)
    monkeypatch.setattr(api().subprocess, 'run', run)
    return calls


def test_sdf_preserves_index_graph_formal_charge_and_coordinates(model):
    _, m = model
    m['formal_charge'] = -1; m['atoms'][0]['formal_charge'] = -1
    text = api().sdf_text(m)
    assert text.splitlines()[3][:6] == '  3  2'
    assert 'M  CHG  1   1  -1' in text
    assert text.endswith('M  END\n$$$$\n')
    atom = text.splitlines()[4]
    assert atom[31:34].strip() == 'O'
    assert int(atom[60:63]) == 1


def test_real_external_command_contract_and_candidate_outputs(tmp_path, monkeypatch, model):
    calls = fixtures(monkeypatch, model)
    result = api().run_am1_bcc(model[0], tmp_path/'fits')
    assert calls and '-dr' not in calls[0]
    assert calls[0][calls[0].index('-c')+1] == 'bcc'
    assert calls[0][calls[0].index('-at')+1] == 'gaff2'
    candidate = json.loads(Path(result['candidate_json']).read_text())
    assert [a['atom_id'] for a in candidate['atoms']] == [1,2,3]
    assert [a['charge_e'] for a in candidate['atoms']] == [-.8,.4,.4]
    assert candidate['coverage']['complete']
    assert result['geometry']['sqm_optimization_converged']
    assert result['final_parameter_validation'] is False
    assert (Path(result['attempt_dir'])/'sqm.out').exists()
    again = api().run_am1_bcc(model[0], tmp_path/'fits')
    assert again['attempt_dir'] != result['attempt_dir']


@pytest.mark.parametrize('failure,q', [('sqm',None),('mapping',None),('exit',None),(None,[-.7,.4,.4]),(None,[float('nan'),.4,.4])])
def test_bad_external_results_never_publish_candidate(tmp_path, monkeypatch, model, failure, q):
    fixtures(monkeypatch, model, failure=failure, q=q)
    with pytest.raises((ValueError, RuntimeError)):
        api().run_am1_bcc(model[0], tmp_path/'fits')
    assert list((tmp_path/'fits').glob('*/AM1_BCC_FAILED.json'))
    assert not list((tmp_path/'fits').glob('*/candidate.json'))
    assert not list((tmp_path/'fits').glob('*/AM1_BCC_COMPLETE.json'))


def test_model_invalid_charge_or_ids_cannot_reach_program(tmp_path, monkeypatch, model):
    p,m=model
    m['atoms'][1]['id']=1; p.write_text(json.dumps(m))
    monkeypatch.setattr(api(), '_program', lambda *a: pytest.fail('invalid model executed'))
    with pytest.raises(ValueError,match='IDs'):
        api().run_am1_bcc(p,tmp_path/'fits')


def test_diagnostic_reference_scope_is_retained(tmp_path,monkeypatch,model):
    p,m=model
    m['metadata'].update(scope='rhf_instability_diagnostic_control',wavefunction_stability_validated=False)
    p.write_text(json.dumps(m)); fixtures(monkeypatch,model)
    result=api().run_am1_bcc(p,tmp_path/'fits')
    assert result['source_reference_scope']=='rhf_instability_diagnostic_control'
    candidate=json.loads(Path(result['candidate_json']).read_text())
    assert candidate['source_reference_scope']=='rhf_instability_diagnostic_control'


def test_am1_geometry_change_cannot_hide_a_broken_bond(tmp_path,monkeypatch,model):
    fixtures(monkeypatch,model)
    original=api().subprocess.run
    def broken(*a,**kw):
        result=original(*a,**kw)
        path=Path(kw['cwd'])/'sqm.out'
        path.write_text(path.read_text().replace('QMMM: 2 2 H 0.95 0.0 0.0','QMMM: 2 2 H 9.5 0.0 0.0'))
        return result
    monkeypatch.setattr(api().subprocess,'run',broken)
    with pytest.raises(ValueError,match='geometry|bond'):
        api().run_am1_bcc(model[0],tmp_path/'fits')
    assert not list((tmp_path/'fits').glob('*/candidate.json'))


def test_slurm_submits_identical_worker_without_computing_science(tmp_path,monkeypatch,model):
    calls=[]
    def submit(args,**kw):
        calls.append(args)
        return subprocess.CompletedProcess(args,0,stdout='12345\n',stderr='')
    monkeypatch.setattr(api().subprocess,'run',submit)
    api().main(['--model',str(model[0]),'--output-dir',str(tmp_path/'bcc'),'--slurm'])
    script=(tmp_path/'bcc_am1bcc_worker.slurm').read_text()
    assert 'module load amber/25' in script
    assert '-m go_param.alternatives' in script
    assert '--slurm' not in script
    assert calls[0][0]=='sbatch'
    assert not (tmp_path/'bcc').exists()


def test_six_decimal_amber_charge_rounding_is_not_given_five_decimal_slack(tmp_path,monkeypatch,model):
    fixtures(monkeypatch,model,q=[-.79999,.4,.4])
    with pytest.raises(ValueError,match='charge sum'):
        api().run_am1_bcc(model[0],tmp_path/'fits')


def test_named_working_input_uses_unique_element_compatible_four_character_names(tmp_path,model):
    path=tmp_path/'input.mol2'; path.write_text(api().input_mol2(model[1]))
    atoms,_=api().read_mol2(path)
    assert [a['name'] for a in atoms]==['O001','H002','H003']
    assert all(len(a['name'])==4 for a in atoms)


def test_five_character_old_names_and_duplicate_output_names_are_rejected(tmp_path,model):
    path=tmp_path/'out.mol2'
    path.write_text(mol2(model[1]).replace('O001','A0001'))
    with pytest.raises(ValueError,match='name|mapping'):
        api().validate_mol2(path,model[1])
    path.write_text(mol2(model[1]).replace('H003','H002'))
    with pytest.raises(ValueError,match='unique|duplicate'):
        api().validate_mol2(path,model[1])


def test_name_range_cannot_silently_exceed_four_characters():
    assert api().atom_name({'id':999,'element':'C'})=='C999'
    with pytest.raises(ValueError,match='999|name'):
        api().atom_name({'id':1000,'element':'C'})


def test_explicit_rounding_policy_retains_raw_and_is_minimum_l2_change():
    raw=[-.127,-.114,.133,.11]
    fixed,proof=api().constrain_rounding(raw,0,source_decimals=3)
    residual=sum(raw)
    assert proof['raw_charge_sum_e']==pytest.approx(residual)
    assert sum(fixed)==pytest.approx(0,abs=1e-12)
    assert all(a-b==pytest.approx(-residual/len(raw)) for a,b in zip(fixed,raw))
    assert raw==[-.127,-.114,.133,.11]
    with pytest.raises(ValueError):api().constrain_rounding([1.,1.],0,source_decimals=3)
    with pytest.raises(ValueError):api().constrain_rounding(raw,0,source_decimals=2)


def test_alternating_carbon_cycle_resonance_preserves_valence_and_fixed_fg():
    model={'atoms':[{'id':i,'element':'C'} for i in range(1,7)],
           'bonds':[{'a':i,'b':i+1 if i<6 else 1,'order':2 if i%2 else 1} for i in range(1,7)]}
    expected={tuple(sorted((b['a'],b['b']))):b['order'] for b in model['bonds']}
    swapped={p:str(3-o) for p,o in expected.items()}
    assert api().validate_bond_representation(model,swapped,[])['alternating_cycle_atom_ids']==[[1,2,3,4,5,6]]
    bad=dict(swapped);bad[(1,2)]='2'
    with pytest.raises(ValueError,match='bond|valence|cycle'):api().validate_bond_representation(model,bad,[])
    model['atoms'][0]['element']='O'
    with pytest.raises(ValueError):api().validate_bond_representation(model,swapped,[])


def test_native_bcc_missing_stage_cannot_be_accepted_even_with_small_charge_error(tmp_path,model):
    with pytest.raises((ValueError,FileNotFoundError)):
        api().validate_bcc_evidence(tmp_path,model[1],[-.8,.4,.4])


def test_active_bcc_table_must_explain_native_per_atom_corrections(tmp_path,monkeypatch,model):
    fixtures(monkeypatch,model)
    result=api().run_am1_bcc(model[0],tmp_path/'fits')
    out=Path(result['attempt_dir'])
    (out/'BCCPARM.DAT').write_text('1 31 91 1 0.1000\n')
    post=(out/'ANTECHAMBER_AM1BCC.AC').read_text().replace('-0.800000','-0.600000').replace('0.400000','0.300000')
    (out/'ANTECHAMBER_AM1BCC.AC').write_text(post)
    proof=api().validate_bcc_evidence(out,model[1],[-.6,.3,.3])
    assert proof['native_BCC_verified'] and proof['nonzero_BCC_bonds']==2
    assert proof['max_abs_BCC_delta_e']==pytest.approx(.2)
    (out/'BCCPARM.DAT').write_text('1 31 91 1 0.0900\n')
    with pytest.raises(ValueError,match='BCC correction'):
        api().validate_bcc_evidence(out,model[1],[-.6,.3,.3])


def test_explicit_policy_exports_distinct_candidate_and_preserves_raw(tmp_path,monkeypatch,model):
    fixtures(monkeypatch,model,q=[-.799,.4,.4])
    result=api().run_am1_bcc(model[0],tmp_path/'fits',charge_policy='constrained-rounding')
    out=Path(result['attempt_dir'])
    candidate=json.loads((out/'candidate.json').read_text())
    raw=json.loads((out/'raw_am1_bcc.json').read_text())
    assert candidate['candidate']=='new_am1_bcc_charge_constrained_vsraw'
    assert raw['candidate']=='raw_native_am1_bcc'
    assert raw['charge_sum_e']==pytest.approx(.001)
    assert candidate['total_charge_e']==pytest.approx(0,abs=1e-12)
    assert candidate['charge_policy']['max_abs_correction_e']==pytest.approx(.001/3)
    assert api().read_mol2(out/'am1_bcc.mol2')[0][0]['charge_e']==-.799


def test_carboxylate_exception_is_limited_to_declared_terminal_resonance_pair():
    path=Path(__file__).resolve().parents[1]/'models/pyrene_COO/model.json'
    m=json.loads(path.read_text())
    actual={tuple(sorted((b['a'],b['b']))):('ar' if b['order']==1.5 else str(int(b['order']))) for b in m['bonds']}
    actual[(17,18)]='1'
    types=[{'id':a['id'],'type':'o' if a['element']=='O' else 'ca' if a['element']=='C' else 'ha'} for a in m['atoms']]
    result=api().validate_bond_representation(m,actual,types)
    assert result['representation_changes'][0]['kind']=='terminal_carboxylate_delocalized_MOL2_representation'
    m['atoms'][18]['symmetry_class']=-999
    with pytest.raises(ValueError):api().validate_bond_representation(m,actual,types)


def test_local_aromatic_metadata_is_graph_bound_and_preserves_source_coordinates(tmp_path,model):
    p,m=model
    prepared=api().prepare_resonance_model(p,tmp_path/'prepared.json')
    assert [a['position'] for a in prepared['atoms']]==[a['position'] for a in m['atoms']]
    assert prepared['metadata']['am1_bcc_resonance']['graph_sha256']==api().graph_digest(m)
    actual={(1,2):'1',(1,3):'1'}
    prepared['atoms'][0]['formal_charge']=-1
    with pytest.raises(ValueError,match='digest'):
        api().validate_bond_representation(prepared,actual,[])
