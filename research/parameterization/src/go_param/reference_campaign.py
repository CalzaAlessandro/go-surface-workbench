"""Isolated QM method comparison worker. No fitting or MD parameter changes."""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import subprocess
import sys
import time

from .reference_qm import make_sp_input, inspect_sp_output, property_job_text, finite_field_dipole
from .qm import gaussian_input, inspect_reference
from .protocols import DFT_METHOD


def write_json(path,value):
    Path(path).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')


def run_gaussian(directory):
    scratch=directory/'scratch';scratch.mkdir(exist_ok=False)
    start=time.time()
    with (directory/'reference.com').open('rb') as inp,(directory/'reference.log').open('wb') as out:
        result=subprocess.run(['g16'],stdin=inp,stdout=out,stderr=subprocess.STDOUT,
                              cwd=directory,env=dict(os.environ,GAUSS_SCRDIR=str(scratch)))
    write_json(directory/'PROCESS_STATUS.json',{'returncode':result.returncode,'seconds':time.time()-start})
    return (directory/'reference.log').read_text(errors='replace')


def compare_properties(directory,model,results):
    """Pairwise direct QM properties; reject unmatched grids rather than zip them."""
    import numpy as np
    from .resp import read_esp
    rows=[]
    accepted=[r for r in results if r.get('accepted') and r.get('esp_extracted')]
    for i,a in enumerate(accepted):
        for b in accepted[i+1:]:
            xa,ga,pa=read_esp(directory/a['method_key']/'reference.esp')
            xb,gb,pb=read_esp(directory/b['method_key']/'reference.esp')
            row={'method_a':a['method_key'],'method_b':b['method_key'],
                 'dipole_vector_difference_debye':math.dist(a['dipole_debye'],b['dipole_debye']),
                 'origin_angstrom':[0,0,0],'formal_charge':model['formal_charge']}
            matched=xa.shape==xb.shape and ga.shape==gb.shape and np.allclose(xa,xb,atol=1e-8,rtol=0) and np.allclose(ga,gb,atol=1e-8,rtol=0)
            row['matched_grid']=bool(matched)
            if matched:
                rms=float(np.sqrt(np.mean((pa-pb)**2)))
                row.update(points=len(pa),potential_rmse_au=rms,
                           rrms_to_a_with_floor=rms/max(float(np.sqrt(np.mean(pa**2))),1e-4),
                           rrms_to_b_with_floor=rms/max(float(np.sqrt(np.mean(pb**2))),1e-4),
                           potential_rms_floor_au=1e-4)
            rows.append(row)
    write_json(directory/'PAIRWISE_QM_PROPERTIES.json',rows)


def run_model(root,index):
    from .resp import read_esp, BOHR_ANGSTROM
    import numpy as np
    manifest=json.loads((root/'INPUT_MANIFEST.json').read_text())
    for rel,sha in manifest['files_sha256'].items():
        if hashlib.sha256((root/rel).read_bytes()).hexdigest()!=sha:
            raise ValueError('frozen input hash mismatch: '+rel)
    cases=json.loads((root/'CASES.json').read_text());case=cases[index]
    directory=root/'outputs'/case['model_id'];directory.mkdir(parents=True,exist_ok=False)
    model=json.loads((root/case['model_file']).read_text())
    if case['geometry_mode']=='new_optimization':
        geom=directory/'geometry';geom.mkdir()
        (geom/'reference.com').write_text((root/case['geometry_input']).read_text())
        text=run_gaussian(geom)
        report=inspect_reference(model,text,method=DFT_METHOD)
        write_json(geom/'GEOMETRY_GATE.json',report)
        if not report['accepted']:
            write_json(directory/'COMPARISON_BLOCKED.json',{'stage':'geometry','issues':report['issues']})
            return
        model=report['optimized_model']
    write_json(directory/'fixed_model.json',model)
    results=[]
    for key in case['methods']:
        d=directory/key;d.mkdir()
        try:
            frozen=root/'fixed_inputs'/case['model_id']/(key+'.com')
            text=frozen.read_text() if frozen.is_file() else make_sp_input(model,key,cpus=16,mem_gb=54)
            (d/'reference.com').write_text(text)
            log=run_gaussian(d)
            result=inspect_sp_output(log,model,key)
            result['esp_extracted']=False
            if result['accepted'] and result['dipole_debye'] is not None:
                (d/'property.log').write_text(property_job_text(log))
                proc=subprocess.run(['espgen','-i','property.log','-o','reference.esp'],cwd=d,capture_output=True,text=True)
                (d/'espgen.log').write_text(proc.stdout+proc.stderr)
                if proc.returncode:raise ValueError('native espgen failed')
                xyz,grid,phi=read_esp(d/'reference.esp')
                if xyz.shape!=(len(model['atoms']),3) or not np.allclose(xyz*BOHR_ANGSTROM,[a['position'] for a in model['atoms']],atol=3e-6,rtol=0):
                    raise ValueError('espgen coordinates changed')
                if not np.isfinite(phi).all() or not np.isfinite(grid).all():raise ValueError('nonfinite potential/grid')
                result.update(esp_extracted=True,grid_points=len(grid))
        except Exception as exc:
            result={'accepted':False,'issues':[str(exc)],'method_key':key,'esp_extracted':False}
        write_json(d/'PROPERTY_GATE.json',result);results.append(result)
        print(json.dumps({'model':model['model_id'],'method':key,'accepted':result['accepted'],'issues':result['issues']}),flush=True)
    if case.get('finite_field') and any(r['method_key']=='ccsd_avdz' and r['accepted'] for r in results):
        field_root=directory/'ccsdt_finite_field';field_root.mkdir()
        zero=field_root/'zero';zero.mkdir()
        (zero/'reference.com').write_text(make_sp_input(model,'ccsdt_avdz',cpus=16,mem_gb=54))
        gate=inspect_sp_output(run_gaussian(zero),model,'ccsdt_avdz')
        write_json(zero/'PROPERTY_GATE.json',gate)
        if gate['accepted']:
            energies={a:{} for a in 'XYZ'};failures=[]
            for axis in 'XYZ':
                for n in [-2,-1,1,2]:
                    d=field_root/f'{axis}{n:+d}';d.mkdir()
                    (d/'reference.com').write_text(make_sp_input(model,'ccsdt_avdz',cpus=16,mem_gb=54,field=(axis,n)))
                    r=inspect_sp_output(run_gaussian(d),model,'ccsdt_avdz',field=(axis,n),require_stability=False)
                    r['separate_zero_field_stability_gate']='../zero/PROPERTY_GATE.json'
                    write_json(d/'PROPERTY_GATE.json',r)
                    if r['accepted']:energies[axis][n]=r['energy_hartree']
                    else:failures.append({'field':[axis,n],'issues':r['issues']})
            if failures:write_json(field_root/'FIELD_FAILURES.json',failures)
            else:
                r=finite_field_dipole(energies,operator_sign=1)
                r.update(origin_angstrom=[0,0,0],ionic_dipole_origin_dependent=model['formal_charge']!=0,
                         runtime_sign_probe_job='1469071',step_test_passed=r['max_step_difference_debye']<.005)
                write_json(field_root/'FINITE_FIELD_DIPOLE.json',r)
    compare_properties(directory,model,results)
    write_json(directory/'COMPARISON_STATUS.json',{'model_id':model['model_id'],'methods':results,
               'reference_selected':False,'production_charges_changed':False})


if __name__=='__main__':
    run_model(Path(sys.argv[1]).resolve(),int(os.environ['SLURM_ARRAY_TASK_ID']))
