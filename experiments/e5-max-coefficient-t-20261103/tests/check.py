#!/usr/bin/env python3
import sys
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
from inputs import (OLD, PROJECT, old_records, reconstruct, name, cv_seed, array_sha,
                    DESIGNS, generator, oracle_copies, g_copies)
from statistic import coefficient_t, max_coefficient_t
import ast
import hashlib
import json
import os
import subprocess
import numpy as np

work = ROOT/'validation'
work.mkdir(exist_ok=True)
checks = {}

# Existing manifests and a before-state inventory, all read-only.
inventory = {}
for manifest in ('CODE_SHA256SUMS', 'SHA256SUMS'):
    lines = (OLD/manifest).read_text().splitlines()
    for line in lines:
        expected, relative = line.split(None, 1)
        relative = relative.strip()
        path = OLD/relative
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        assert actual == expected, (manifest, relative)
        inventory[str(path.relative_to(PROJECT))] = actual
    checks[manifest] = len(lines)
for path in [OLD/'SHA256SUMS', OLD/'results/complete/raw.jsonl',
             PROJECT/'experiments/readable-experiments/R/gaussian_tests.R',
             PROJECT/'experiments/readable-experiments/R/linear_tests.R']:
    inventory[str(path.relative_to(PROJECT))] = hashlib.sha256(path.read_bytes()).hexdigest()
before = work/'original-files-before.json'
if before.exists():
    assert json.loads(before.read_text()) == inventory
else:
    before.write_text(json.dumps(inventory, indent=2, sort_keys=True)+'\n')

tree = ast.parse((OLD/'run_pilot.py').read_text())
designs = [node for node in tree.body if isinstance(node, ast.Assign)
           and any(isinstance(t, ast.Name) and t.id == 'DESIGNS' for t in node.targets)]
assert ast.literal_eval(designs[0].value) == DESIGNS
rows = old_records()
seeds = [cv_seed(row) for row in rows]
assert len(set(seeds)) == len(rows)
checks['unique_cv_seeds'] = len(seeds)
selected = [r for r in rows if r['rep'] == 10000]
for row in selected:
    reconstruct(row)
checks['original_data_hashes'] = len(selected)

row = selected[0]
x,y,target,s,om,graph,key = reconstruct(row)
candidates = np.concatenate((x[None,:,target],
    oracle_copies(x,target,om,6,generator(row['design_seed'],'oracle',*key)),
    g_copies(x,target,graph,5,3,generator(row['design_seed'],'g',3,*key))[:,:,target]),axis=0)
with (work/'reference-input.bin').open('wb') as output:
    np.asarray(y,dtype='<f8').tofile(output)
    np.asarray(x[:,s],dtype='<f8').ravel(order='F').tofile(output)
    for candidate in candidates:
        np.asarray(candidate,dtype='<f8').ravel(order='F').tofile(output)
env = dict(os.environ, TMPDIR=str(work))
result = subprocess.run(['Rscript','--vanilla',str(ROOT/'tests/reference.R'),str(ROOT),str(work)],
                        check=True,capture_output=True,text=True,cwd=ROOT,env=env)
(work/'R-reference.log').write_text(result.stdout+result.stderr)
residual = np.fromfile(work/'R-residual.bin',dtype='<f8')
expected = np.loadtxt(work/'R-coefficient-t.csv',delimiter=',')
before_inputs = array_sha(candidates,residual)
actual = coefficient_t(candidates,residual)
error = float(np.max(np.abs(actual-expected)))
assert np.allclose(actual,expected,atol=2e-11,rtol=2e-11),error
checks['R_glm_max_absolute_t_error'] = error
assert array_sha(candidates,residual) == before_inputs
assert np.array_equal(max_coefficient_t(candidates,residual),np.max(np.abs(actual),axis=1))
order = np.array([4,0,9,1,2,3,5,6,7,8,10,11])
assert np.array_equal(coefficient_t(candidates[order],residual),actual[order])
assert np.allclose(coefficient_t(candidates[0],residual),actual[0],atol=1e-12)
scales = np.array([.1,-2,3,.3,-.5,4,2,10])
assert np.allclose(coefficient_t(candidates*scales,residual),actual*np.sign(scales),atol=2e-11)
assert np.allclose(coefficient_t(candidates,3*residual+7),actual,atol=2e-11)
assert np.allclose(max_coefficient_t(candidates[:,:,::-1],residual),np.max(abs(actual),axis=1),atol=2e-11)
marginal = np.array([coefficient_t(candidates[0,:,j:j+1],residual)[0] for j in range(8)])
assert np.max(np.abs(marginal-actual[0])) > .1
checks['invariances_and_joint_adjustment'] = True
for bad_x,bad_y in [(np.ones((80,8)),residual),
                    (candidates[0,:9],residual[:9]),
                    (candidates[0],np.zeros(80)),
                    (candidates[0]*np.nan,residual)]:
    try:
        max_coefficient_t(bad_x,bad_y)
    except ValueError:
        pass
    else:
        raise AssertionError('Invalid input was silently scored')
checks['rank_df_variance_finiteness_failures'] = True
checks['old_methods_evaluated'] = 0
(work/'implementation-checks.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
