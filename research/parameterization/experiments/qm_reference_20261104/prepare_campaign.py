"""Freeze isolated comparison inputs. Run once; refuses existing directories."""
import hashlib,json,pathlib,shutil
from rdkit import Chem
from rdkit.Chem import AllChem,rdMolDescriptors
from go_param.reference_qm import make_sp_input
from go_param.qm import gaussian_input
from go_param.protocols import DFT_METHOD

ROOT=pathlib.Path(__file__).resolve().parents[2]
CALC=ROOT/'calculations'
DFT=['wb97xd_pople','wb97xd_avtz','pbe0_avtz','pw6b95_avtz','pw6b95_avdz']

def dump(path,value):path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def campaign(name,cases,models,*,cpus,mem,concurrency,geometry_inputs=False):
    root=CALC/name;root.mkdir(exist_ok=False)
    code=root/'code/go_param';code.mkdir(parents=True)
    for file in ['reference_qm.py','reference_campaign.py','qm.py','protocols.py','resp.py']:
        shutil.copyfile(ROOT/'src/go_param'/file,code/file)
    for m,case in zip(models,cases):
        d=root/'models'/m['model_id'];d.mkdir(parents=True)
        dump(d/'model.json',m);case['model_file']=str((d/'model.json').relative_to(root))
        if case['geometry_mode']=='new_optimization':
            (d/'geometry.com').write_text(gaussian_input(m,cpus=cpus,mem_gb=mem-10,method=DFT_METHOD))
            case['geometry_input']=str((d/'geometry.com').relative_to(root))
        else:
            d=root/'fixed_inputs'/m['model_id'];d.mkdir(parents=True)
            for key in case['methods']:(d/(key+'.com')).write_text(make_sp_input(m,key,cpus=cpus,mem_gb=mem-10 if mem>10 else mem-2))
    dump(root/'CASES.json',cases)
    (root/'run.slurm').write_text(f'''#!/bin/bash -l
#SBATCH --job-name=GO-QM-method-validation
#SBATCH --partition=ciamician_cpu
#SBATCH --qos=ciamician_cpu_caps
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task={cpus}
#SBATCH --mem={mem}G
#SBATCH --time=5-00:00:00
#SBATCH --array=0-{len(cases)-1}%{concurrency}
#SBATCH --output=reference-%A_%a.out
#SBATCH --error=reference-%A_%a.err
set -eo pipefail
module load gaussian/16
module load amber/25
set -u
cd "${{SLURM_SUBMIT_DIR}}"
export PYTHONPATH="$PWD/code"
export OPENBLAS_NUM_THREADS=1
for candidate in "$AMBERHOME/bin/amber.python" "$AMBERHOME/miniconda/bin/python" python python3; do
 if "$candidate" -c 'import numpy' >/dev/null 2>&1; then exec "$candidate" -m go_param.reference_campaign "$PWD"; fi
done
exit 5
''')
    dump(root/'INPUT_MANIFEST.json',{'purpose':'Validate QM methods; not a charge/force-field preset','probe_job_id':'1469071','methods':'Explicit gas-phase fixed-geometry method/basis comparison; no DFT-derived MD charge replacement','resources':{'cpus':cpus,'mem_gb':mem,'walltime':'5-00:00:00','concurrency':concurrency},'files_sha256':{str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in sorted(root.rglob('*')) if f.is_file()}})
    print(root)

water=json.loads((CALC/'qm_reference_runtime_20261104_attempt_001/model.json').read_text())
campaign('qm_reference_integration_20261104_attempt_001',[{'model_id':water['model_id'],'geometry_mode':'already_fixed','methods':['pw6b95_avtz','mp2_avdz','ccsd_avdz']}],[water],cpus=4,mem=8,concurrency=1)

ids=['pyrene_pristine','pyrene_diol_trans','pyrene_diol_epoxide','pyrene_dione','pyrene_COOH','pyrene_COO','C24_OH_EPX_near_trans','C24_OH_EPX_far_trans']
models=[];cases=[]
for model_id in ids:
    source=CALC/'complementary_results_20261103_snapshot_001/atomic_inputs'/model_id/'DFT/model.json'
    model=json.loads(source.read_text())
    assert model['metadata']['wavefunction_stability_validated']
    model['metadata']['comparison_source']={'path':str(source.relative_to(ROOT)),'sha256':hashlib.sha256(source.read_bytes()).hexdigest()}
    models.append(model);cases.append({'model_id':model_id,'geometry_mode':'already_fixed','methods':['hf_631gstar']+DFT})
campaign('qm_reference_go_controls_20261104_attempt_001',cases,models,cpus=16,mem=64,concurrency=2)

models=[];cases=[]
anchors=[('benzene','c1ccccc1'),('phenol','Oc1ccccc1'),('benzene_oxide','C1=CC2OC2C=C1'),('p_benzoquinone','O=C1C=CC(=O)C=C1'),('acetic_acid','CC(=O)O'),('acetate','CC(=O)[O-]')]
for label,smiles in anchors:
    mol=Chem.AddHs(Chem.MolFromSmiles(smiles));params=AllChem.ETKDGv3();params.randomSeed=20261104
    assert AllChem.EmbedMolecule(mol,params)==0
    assert AllChem.MMFFHasAllMoleculeParams(mol)
    AllChem.MMFFOptimizeMolecule(mol,maxIters=2000)
    conf=mol.GetConformer()
    q=sum(a.GetFormalCharge() for a in mol.GetAtoms());model_id='anchor_'+label
    model={'schema_version':1,'model_id':model_id,'formula':rdMolDescriptors.CalcMolFormula(mol),'formal_charge':q,'multiplicity':1,'atoms':[{'id':a.GetIdx()+1,'element':a.GetSymbol(),'position':list(conf.GetAtomPosition(a.GetIdx()))} for a in mol.GetAtoms()],'bonds':[{'a':b.GetBeginAtomIdx()+1,'b':b.GetEndAtomIdx()+1,'order':b.GetBondTypeAsDouble()} for b in mol.GetBonds()],'metadata':{'smiles':smiles,'purpose':'Small chemical anchor for independent correlated QM controls, not GO library training','initial_geometry':'RDKit ETKDGv3 seed 20261104 + MMFF; NOT an accepted QM geometry'}}
    models.append(model);cases.append({'model_id':model_id,'geometry_mode':'new_optimization','methods':['hf_631gstar']+DFT+['mp2_avdz','mp2_avtz','ccsd_avdz'],'finite_field':label in ['benzene_oxide','acetate']})
campaign('qm_reference_anchors_20261104_attempt_001',cases,models,cpus=16,mem=64,concurrency=2)
