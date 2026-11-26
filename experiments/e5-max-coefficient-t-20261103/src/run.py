#!/usr/bin/env python3
"""Resume only new MaxCoefT scoring; reuse every old p-value verbatim."""
from inputs import (ROOT, OLD, PROJECT, MATCHED_RAW, CONFIG, DESIGNS, KEYS, old_records,
                    reconstruct, identity, name, cv_seed, array_sha, generator,
                    g_copies, oracle_copies, rank_pvalue)
import argparse
import atexit
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback
from zoneinfo import ZoneInfo
import numpy as np
import scipy
from statistic import max_coefficient_t

_R = None
_WORK = None
_OUT = None


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    temporary = path.with_suffix(path.suffix + '.tmp')
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n')
    temporary.replace(path)


def stamp():
    return datetime.now(ZoneInfo('Asia/Singapore')).isoformat()


def init_worker(out):
    global _R, _WORK, _OUT
    _OUT = Path(out)
    _WORK = _OUT / 'work' / str(os.getpid())
    _WORK.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, TMPDIR=str(_WORK), R_USER_CACHE_DIR=str(_WORK))
    stderr = (_WORK / 'R-stderr.log').open('a')
    _R = subprocess.Popen(['Rscript', '--vanilla', str(ROOT/'src/background_server.R'), str(ROOT/'src')],
                          stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr,
                          text=True, bufsize=1, cwd=ROOT, env=env)
    atexit.register(_R.terminate)


def background(x, y, seed):
    input_path, output_path = _WORK/'background-input.bin', _WORK/'background-output.bin'
    np.asarray(np.column_stack((x, y)), dtype='<f8').ravel(order='F').tofile(input_path)
    task = dict(input=str(input_path), output=str(output_path), n=len(y), p=x.shape[1], seed=seed)
    _R.stdin.write(json.dumps(task)+'\n'); _R.stdin.flush()
    response = _R.stdout.readline()
    if not response:
        raise RuntimeError(f'R background process exited; see {_WORK}/R-stderr.log')
    metadata = json.loads(response)
    if not metadata['ok']:
        raise RuntimeError(metadata)
    residual = np.fromfile(output_path, dtype='<f8')
    if residual.shape != y.shape or not np.all(np.isfinite(residual)):
        raise ValueError('Invalid R residual output')
    return residual, metadata


def one_replication(old):
    start = time.perf_counter()
    ident = name(old)
    try:
        x, y, target, s, omega, graph, key = reconstruct(old)
        seed = cv_seed(old)
        residual, lasso = background(x[:, s], y, seed)
        background_seconds = time.perf_counter() - start
        observed = max_coefficient_t(x[:, target], residual)
        row = {k: old[k] for k in KEYS}
        row.update(data_sha256=array_sha(x, y), data_hash_matches_old=True,
                   design_seed=old['design_seed'], data_stream=old['data_stream'],
                   cv_seed=seed, cv_rng=['Mersenne-Twister', 'Inversion', 'Rejection'],
                   lasso=lasso, residual_sha256=array_sha(residual), M=999, L=3,
                   alpha=.05, pvalues={}, reject={}, diagnostics={},
                   reused_pvalues=old['pvalues'], reused_reject=old['reject'])
        saved = {'residual_response': residual}
        for sampler in ('Oracle', 'G'):
            clock = time.perf_counter()
            stream = ['oracle', *key] if sampler == 'Oracle' else ['g', 3, *key]
            rng = generator(old['design_seed'], *stream)
            if sampler == 'Oracle':
                xt = oracle_copies(x, target, omega, 999, rng)
            else:
                full = g_copies(x, target, graph, 999, 3, rng)
                xt = full[:, :, target]
                del full
            candidates_hash = array_sha(xt)
            values = max_coefficient_t(xt, residual)
            del xt
            method = sampler + 'MaxCoefT'
            scores = np.r_[observed, values]
            pvalue = rank_pvalue(observed, values)
            row['pvalues'][method] = pvalue
            row['pvalues'][sampler+'RidgeMaxCoefTBonf'] = min(1., 2*min(old['pvalues'][sampler+'Ridge'], pvalue))
            rounded = np.round(scores, 12)
            row['diagnostics'][method] = dict(observed=observed, minimum=float(values.min()),
                maximum=float(values.max()), ties=int(np.sum(rounded[1:] == rounded[0])),
                greater=int(np.sum(rounded[1:] > rounded[0])), score_sha256=array_sha(scores),
                candidates_sha256=candidates_hash, stream=stream,
                seconds=time.perf_counter()-clock)
            saved[method] = scores
        row['reject'] = {k: int(v <= .05) for k, v in row['pvalues'].items()}
        row['seconds'] = time.perf_counter() - start
        row['background_seconds'] = background_seconds
        score_path = _OUT/'scores'/f'{ident}.npz'
        temporary = score_path.with_suffix('.tmp.npz')
        np.savez_compressed(temporary, **saved)
        temporary.replace(score_path)
        row['scores_file_sha256'] = sha(score_path)
        write_json(_OUT/'records'/f'{ident}.json', row)
        return dict(ok=True, name=ident, seconds=row['seconds'])
    except Exception:
        failure = dict(ok=False, name=ident, identity=identity(old), time=stamp(),
                       traceback=traceback.format_exc())
        # Every failed attempt survives retries.
        write_json(_OUT/'failures'/f'{ident}-{time.time_ns()}.json', failure)
        return failure


def source_hashes():
    paths = sorted((ROOT/'src').glob('*.py')) + sorted((ROOT/'src').glob('*.R'))
    paths += [OLD/'src/mscrt/core.py', OLD/'src/mscrt/experiments.py',
              OLD/'src/run_confirmation.py', OLD/'run_pilot.py', OLD/'configs/full.json',
              MATCHED_RAW,
              PROJECT/'experiments/readable-experiments/R/gaussian_tests.R',
              PROJECT/'experiments/readable-experiments/R/linear_tests.R']
    return {str(p.relative_to(PROJECT)): sha(p) for p in paths}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', choices=['pilot', 'full'], required=True)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--out', default='results')
    args = parser.parse_args()
    if not 1 <= args.workers <= 6:
        raise ValueError('Use 1 to 6 workers')
    out = (ROOT/args.out).resolve()
    if not out.is_relative_to(ROOT) or out == ROOT:
        raise ValueError('All output must be inside the new experiment directory')
    for folder in ('records', 'scores', 'failures', 'work', 'logs'):
        (out/folder).mkdir(parents=True, exist_ok=True)
    sources = source_hashes()
    binding = out/'source_binding.json'
    if binding.exists():
        if json.loads(binding.read_text()) != sources:
            raise ValueError('Source changed since this run; use a new output directory')
    else:
        write_json(binding, sources)
    rows = old_records()
    if args.stage == 'pilot':
        rows = [r for r in rows if r['rep'] < 10002]
    pending = []
    for row in rows:
        path = out/'records'/f'{name(row)}.json'
        if path.exists():
            done = json.loads(path.read_text())
            scores = out/'scores'/f'{name(row)}.npz'
            if (identity(done) != identity(row) or done['data_sha256'] != row['data_sha256']
                    or done['scores_file_sha256'] != sha(scores)):
                raise ValueError(f'Invalid checkpoint {path}')
        else:
            pending.append(row)
    # Interleave designs/cells so the timing pilot and progress cover the grid.
    pending.sort(key=lambda r: (r['rep'], r['design'], r['pattern'], r['theta']))
    start = time.perf_counter()
    run_id = datetime.now(ZoneInfo('Asia/Singapore')).strftime('%Y%m%d-%H%M%S')
    metadata = dict(stage=args.stage, started=stamp(), status='running', expected=len(rows),
        resumed=len(rows)-len(pending), workers=args.workers, python=sys.version,
        numpy=np.__version__, scipy=scipy.__version__, platform=platform.platform(),
        command=[sys.executable, *sys.argv], config=CONFIG)
    meta_path = out/'logs'/f'{run_id}-{args.stage}.json'
    write_json(meta_path, metadata)
    completed = failures = 0
    print(json.dumps(dict(output=str(out), pending=len(pending), workers=args.workers)), flush=True)
    try:
        with ProcessPoolExecutor(max_workers=args.workers, initializer=init_worker, initargs=(str(out),)) as pool:
            futures = [pool.submit(one_replication, row) for row in pending]
            for future in as_completed(futures):
                result = future.result()
                completed += int(result['ok']); failures += int(not result['ok'])
                if not result['ok']:
                    print(json.dumps(result), flush=True)
                if (completed+failures) % 25 == 0 or completed+failures == len(pending):
                    progress = dict(time=stamp(), completed=completed, failed=failures,
                        total=len(pending), elapsed_seconds=time.perf_counter()-start)
                    with (out/'logs'/f'{run_id}-progress.jsonl').open('a') as log:
                        log.write(json.dumps(progress)+'\n')
                    print(json.dumps(progress), flush=True)
        metadata['status'] = 'complete' if failures == 0 else 'incomplete'
    except BaseException:
        metadata.update(status='interrupted', traceback=traceback.format_exc())
        raise
    finally:
        metadata.update(finished=stamp(), completed=completed, failures=failures,
                        elapsed_seconds=time.perf_counter()-start)
        write_json(meta_path, metadata)
    if failures:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
