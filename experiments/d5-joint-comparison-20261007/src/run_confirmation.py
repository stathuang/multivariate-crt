#!/usr/bin/env python3
"""Persist paired CRT experiments. Every run has a new, non-overwritten directory.

Raw JSONL is flushed after each dataset. Completed replicate identities, p-values,
diagnostics, paired indicators, seeds, sources and failures remain inspectable.
"""
import os
for _name in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[_name] = '1'
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import hashlib
import itertools
import json
import platform
import shlex
import time
import traceback
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import scipy
from scipy.stats import beta
import geometry
from geometry import reference, maxpartial, bonferroni_two
from mscrt.core import RidgeStatistic, conditional_parameters, coordinate_pvalues, g_copies, oracle_copies, rank_pvalue
from mscrt.experiments import generator, graph_design

ROOT = Path(__file__).resolve().parents[1]
_DESIGN_CACHE = {}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def array_sha(*arrays):
    h = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        h.update(str((a.shape, a.dtype.str)).encode())
        h.update(a.tobytes())
    return h.hexdigest()


def one_dataset(job):
    cfg, design, arrangement, pattern, theta, rep = job
    start = time.perf_counter()
    spec = reference.DESIGNS[design]
    n, p, t, m = cfg['n'], cfg['p'], cfg['t'], cfg['M']
    target = np.arange(t)
    cache_key = (design, arrangement, n, p, t)
    if cache_key not in _DESIGN_CACHE:
        om, sigma, graph = graph_design(p, t, arrangement, spec['bandwidth'], spec['offdiag'])
        s, b, v = conditional_parameters(om, target)
        _DESIGN_CACHE[cache_key] = om, np.linalg.cholesky(sigma), graph, s, b, v
    om, chol, graph, s, b, v = _DESIGN_CACHE[cache_key]
    direction = np.ones(t)/np.sqrt(t) if pattern != 'sparse' else np.eye(t)[0]
    direction = direction/np.sqrt(direction @ v @ direction)
    gamma = np.zeros(len(s)); gamma[:4] = [1., -.7, .5, .3]
    key = (arrangement, pattern, float(theta), int(rep))
    rng = generator(spec['seed'], 'data', *key)
    x = rng.normal(size=(n, p)) @ chol.T
    eps = rng.normal(size=n)
    y = x[:, s] @ gamma + theta*((x[:, target]-x[:, s]@b.T)@direction) + eps
    stat = RidgeStatistic.fit(y, x[:, s], n*cfg['nuisance_penalty'], cfg['target_penalty'])
    pcoord, pcoord_oracle = coordinate_pvalues(x, target, graph, stat, om)
    partial0, detail = maxpartial(x[:, target], x, target, graph, stat, details=True)
    obs = dict(Ridge=stat(x[:, target]), MaxPartial=float(partial0[0]),
               MaxCorr=float(reference.maxcorr_scores(x[:, target], stat)[0]))
    pvalues = {'CoordG':pcoord, 'CoordOracle':pcoord_oracle}
    diagnostics = {}
    scores_saved = {}
    for sampler in ('Oracle', 'G'):
        if sampler == 'Oracle':
            xt = oracle_copies(x, target, om, m, generator(spec['seed'], 'oracle', *key))
        else:
            full = g_copies(x, target, graph, m, cfg['L'], generator(spec['seed'], 'g', cfg['L'], *key))
            xt = full[:, :, target]
        scores = dict(Ridge=stat(xt), MaxPartial=maxpartial(xt, x, target, graph, stat),
                      MaxCorr=reference.maxcorr_scores(xt, stat))
        for name, values in scores.items():
            method = sampler+name
            pvalues[method] = rank_pvalue(obs[name], values)
            rounded = np.round(values, 12)
            diagnostics[method] = dict(observed=obs[name], minimum=float(values.min()), maximum=float(values.max()),
                ties=int(np.sum(rounded == np.round(obs[name], 12))),
                greater=int(np.sum(rounded > np.round(obs[name], 12))),
                score_sha256=array_sha(np.r_[obs[name],values]))
            if rep == cfg['rep_start']:
                scores_saved[method] = np.r_[obs[name], values]
        pvalues[sampler+'AdaptiveBonf'] = bonferroni_two(pvalues[sampler+'Ridge'],pvalues[sampler+'MaxPartial'])
    row = dict(design=design, arrangement=arrangement, pattern=pattern, theta=theta, rep=rep,
        design_seed=spec['seed'], data_stream=['data',*key], M=m, L=cfg['L'],
        data_sha256=array_sha(x,y), pvalues=pvalues,
        reject={name:int(value<=cfg['alpha']) for name,value in pvalues.items()},
        observed_local_p=detail['local_p'][0].tolist(), observed_w=detail['squared_partial_correlation'][0].tolist(),
        observed_df=detail['df'][0].tolist(), score_diagnostics=diagnostics,
        conditional_effects=(np.sqrt(n)*theta*np.abs(direction)/np.sqrt(np.diag(om)[:t])).tolist(),
        seconds=time.perf_counter()-start)
    # Dataset-level paired differences retained explicitly, as well as in summaries.
    row['paired_reject_differences'] = {
        a+' minus '+b: row['reject'][a]-row['reject'][b]
        for a,b in comparisons()
    }
    return row, scores_saved


def comparisons():
    for sampler in ('G', 'Oracle'):
        for a,b in [('MaxPartial','Ridge'),('AdaptiveBonf','Ridge'),('AdaptiveBonf','MaxPartial'),('MaxPartial','MaxCorr')]:
            yield sampler+a, sampler+b
        for name in ('Ridge','MaxPartial','AdaptiveBonf'):
            yield sampler+name, 'CoordG'
    for name in ('Ridge','MaxPartial','AdaptiveBonf'):
        yield 'G'+name, 'Oracle'+name


def write_csv(path, rows):
    with path.open('w', newline='') as f:
        writer=csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def summarize(out, records):
    groups = {}
    for r in records:
        groups.setdefault(tuple(r[k] for k in ('design','arrangement','pattern','theta')), []).append(r)
    summaries, pairs = [], []
    for key, rows in sorted(groups.items()):
        common=dict(zip(('design','arrangement','pattern','theta'), key))
        count=len(rows)
        for method in rows[0]['pvalues']:
            values=np.array([r['reject'][method] for r in rows])
            k=int(values.sum()); rate=k/count
            lo=0. if k==0 else float(beta.ppf(.025,k,count-k+1))
            hi=1. if k==count else float(beta.ppf(.975,k+1,count-k))
            summaries.append(dict(**common,method=method,reps=count,rejections=k,rate=rate,
                mcse=float(values.std(ddof=1)/np.sqrt(count)) if count>1 else None,exact_low=lo,exact_high=hi))
        for a,b in comparisons():
            d=np.array([r['reject'][a]-r['reject'][b] for r in rows])
            pairs.append(dict(**common,method=a,reference=b,reps=count,difference=float(d.mean()),
                paired_mcse=float(d.std(ddof=1)/np.sqrt(count)) if count>1 else None,
                wins=int(np.sum(d==1)),losses=int(np.sum(d==-1))))
    write_csv(out/'summary.csv',summaries)
    write_csv(out/'paired_differences.csv',pairs)


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--config',required=True)
    ap.add_argument('--workers',type=int,default=2)
    args=ap.parse_args()
    cfg=json.loads(Path(args.config).read_text())
    if cfg['rep_stop']<=cfg['rep_start'] or cfg['M']<1 or cfg['L']<0 or args.workers<1:
        raise ValueError('Invalid replicate range, copies, sweeps or workers.')
    stamp=datetime.now(ZoneInfo('Asia/Singapore')).strftime('%Y%m%d-%H%M%S')
    out=ROOT/'runs'/(stamp+'-'+cfg['name'])
    out.mkdir(parents=True,exist_ok=False)
    (out/'scores').mkdir()
    (out/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
    sources=[Path(__file__),Path(geometry.__file__),ROOT/'run_pilot.py',Path(reference.V17_CODE)/'mscrt/core.py',Path(reference.V17_CODE)/'mscrt/experiments.py']
    sources += sorted((ROOT/'tests').glob('test_*.py'))
    sources += [p for p in (ROOT/'src').rglob('*.py') if p not in sources]
    meta_sources_relative = {str(p.relative_to(ROOT)):sha(p) for p in sources}
    meta=dict(started=stamp,status='running',python=sys.version,numpy=np.__version__,scipy=scipy.__version__,
        platform=platform.platform(),cwd=str(Path.cwd()),command=shlex.join([sys.executable,*sys.argv]),workers=args.workers,
        source_sha256={str(p):sha(p) for p in sources},source_relative_sha256=meta_sources_relative,config_sha256=sha(args.config),
        thread_limits={k:os.environ[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS')})
    (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
    print(out,flush=True)
    jobs=[(cfg,d,g,cell['pattern'],cell['theta'],r)
        for d,g,cell,r in itertools.product(cfg['designs'],cfg['arrangements'],cfg['cells'],range(cfg['rep_start'],cfg['rep_stop']))]
    records=[]; start=time.perf_counter()
    try:
        with (out/'raw.jsonl').open('x') as raw, (out/'progress.log').open('x') as log, ProcessPoolExecutor(max_workers=args.workers) as pool:
            for row,scores in pool.map(one_dataset,jobs,chunksize=1):
                raw.write(json.dumps(row,sort_keys=True)+'\n');raw.flush()
                records.append(row)
                if scores:
                    ident=f"{row['design']}-{row['arrangement']}-{row['pattern']}-{row['theta']:g}-r{row['rep']}"
                    np.savez_compressed(out/'scores'/(ident+'.npz'),**scores)
                if len(records)%25==0 or len(records)==len(jobs):
                    message=f'{len(records)}/{len(jobs)} datasets; {time.perf_counter()-start:.1f} seconds'
                    print(message,flush=True);log.write(message+'\n');log.flush()
        summarize(out,records)
        meta.update(status='complete',datasets=len(records),method_records=sum(len(r['pvalues']) for r in records),
                    elapsed_seconds=time.perf_counter()-start,finished=datetime.now(ZoneInfo('Asia/Singapore')).isoformat())
    except BaseException:
        meta.update(status='failed',datasets=len(records),elapsed_seconds=time.perf_counter()-start,exception=traceback.format_exc())
        (out/'failure.txt').write_text(meta['exception'])
        if records:
            summarize(out,records)
        raise
    finally:
        (out/'metadata.json').write_text(json.dumps(meta,indent=2)+'\n')
        paths=sorted(p for p in out.rglob('*') if p.is_file() and p.name!='SHA256SUMS')
        (out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in paths))


if __name__=='__main__':
    main()
