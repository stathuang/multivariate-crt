"""Read-only reuse of the original generators; no old method is evaluated."""
import os
import sys
from pathlib import Path

for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ[name] = "1"
os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
sys.dont_write_bytecode = True

import hashlib
import json
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parents[1]
OLD = ROOT.parent / "d5-joint-comparison-20261007"
MATCHED_RAW = Path(os.environ.get('REPLICATION_MATCHED_RAW', str(OLD/'results/complete/raw.jsonl'))).resolve()
if not MATCHED_RAW.is_relative_to(PROJECT):
    raise ValueError('Matched input must be inside this replication package')
sys.path.insert(0, str(OLD / "src"))
from mscrt.core import conditional_parameters, g_copies, oracle_copies, rank_pvalue
from mscrt.experiments import generator, graph_design

# Source: old run_pilot.py DESIGNS. The values are checked against its AST.
DESIGNS = {"A": {"seed": 27092611, "bandwidth": 6, "offdiag": .20},
           "B": {"seed": 27092612, "bandwidth": 1, "offdiag": -.45}}
CONFIG = json.loads((OLD / "configs/full.json").read_text())
KEYS = ("design", "arrangement", "pattern", "theta", "rep")
_CACHE = {}


def identity(row):
    return tuple(row[k] for k in KEYS)


def name(row):
    return f"{row['design']}-{row['arrangement']}-{row['pattern']}-{row['theta']:g}-r{row['rep']}"


def array_sha(*arrays):
    digest = hashlib.sha256()
    for a in arrays:
        a = np.ascontiguousarray(a)
        digest.update(str((a.shape, a.dtype.str)).encode())
        digest.update(a.tobytes())
    return digest.hexdigest()


def cv_seed(row):
    # Independent R stream, stable across worker counts and execution order.
    keys = [row['design_seed'], 'e5-max-coefficient-t-background-cv',
            row['arrangement'], row['pattern'], float(row['theta']), int(row['rep'])]
    digest = hashlib.sha256(json.dumps(keys, sort_keys=True).encode()).digest()
    return 1 + int.from_bytes(digest[:8], 'little') % 2147483646


def reconstruct(row):
    """The data-only lines of old run_confirmation.one_dataset, unchanged."""
    cfg = CONFIG
    design, arrangement, pattern, theta, rep = identity(row)
    spec = DESIGNS[design]
    n, p, t = cfg['n'], cfg['p'], cfg['t']
    target = np.arange(t)
    cache_key = (design, arrangement, n, p, t)
    if cache_key not in _CACHE:
        om, sigma, graph = graph_design(p, t, arrangement, spec['bandwidth'], spec['offdiag'])
        s, b, v = conditional_parameters(om, target)
        _CACHE[cache_key] = om, np.linalg.cholesky(sigma), graph, s, b, v
    om, chol, graph, s, b, v = _CACHE[cache_key]
    direction = np.ones(t)/np.sqrt(t) if pattern != 'sparse' else np.eye(t)[0]
    direction = direction/np.sqrt(direction @ v @ direction)
    gamma = np.zeros(len(s)); gamma[:4] = [1., -.7, .5, .3]
    key = (arrangement, pattern, float(theta), int(rep))
    rng = generator(spec['seed'], 'data', *key)
    x = rng.normal(size=(n, p)) @ chol.T
    eps = rng.normal(size=n)
    y = x[:, s] @ gamma + theta*((x[:, target]-x[:, s]@b.T)@direction) + eps
    if array_sha(x, y) != row['data_sha256']:
        raise ValueError(f"Original data hash mismatch: {name(row)}")
    return x, y, target, s, om, graph, key


def old_records():
    rows = [json.loads(line) for line in MATCHED_RAW.read_text().splitlines()]
    expected = {(d, 'connected', c['pattern'], float(c['theta']), r)
                for d in CONFIG['designs'] for c in CONFIG['cells']
                for r in range(10000, 10500)}
    if len(rows) != 7000 or {identity(r) for r in rows} != expected:
        raise ValueError("Old records do not match the authorised 14 x 500 grid")
    for row in rows:
        if row['M'] != 999 or row['L'] != 3 or row['design_seed'] != DESIGNS[row['design']]['seed']:
            raise ValueError("Old metadata differs from fixed protocol")
    return sorted(rows, key=identity)
