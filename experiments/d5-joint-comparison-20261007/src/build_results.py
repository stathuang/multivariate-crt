#!/usr/bin/env python3
"""Combine verified, disjoint runs into the complete 14-cell D.5 grid."""
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path
from run_confirmation import summarize

ROOT = Path(__file__).resolve().parents[1]
KEYS = ('design', 'arrangement', 'pattern', 'theta', 'rep')
SHARED = ('n', 'p', 't', 'M', 'L', 'alpha', 'nuisance_penalty', 'target_penalty',
          'designs', 'arrangements', 'rep_start', 'rep_stop')
METHODS = ('CoordG', 'GRidge', 'GMaxPartial', 'GAdaptiveBonf',
           'OracleRidge', 'OracleMaxPartial', 'OracleAdaptiveBonf')


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def identity(row):
    return (row['design'], row['arrangement'], row['pattern'],
            float(row['theta']), int(row['rep']))


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--runs', nargs='+', required=True)
    ap.add_argument('--out', required=True, help='A new output directory; existing outputs are never overwritten.')
    args = ap.parse_args()
    full = json.loads((ROOT/'configs/full.json').read_text())
    expected = {(d, a, c['pattern'], float(c['theta']), r)
                for d, a, c, r in itertools.product(full['designs'], full['arrangements'],
                                                  full['cells'], range(full['rep_start'], full['rep_stop']))}
    records, seen, cells, sources = [], set(), {}, []
    for name in args.runs:
        run = Path(name).resolve()
        cfg = json.loads((run/'config.json').read_text())
        meta = json.loads((run/'metadata.json').read_text())
        if meta['status'] != 'complete':
            raise ValueError(f'Incomplete run: {run}')
        for key in SHARED:
            if cfg[key] != full[key]:
                raise ValueError(f'Configuration mismatch: {run}: {key}')
        for line in (run/'SHA256SUMS').read_text().splitlines():
            digest, path = line.split('  ', 1)
            if sha(run/path) != digest:
                raise ValueError(f'Output hash mismatch: {run/path}')
        source_rows = [json.loads(line) for line in (run/'raw.jsonl').read_text().splitlines()]
        expected_run = {(d, a, c['pattern'], float(c['theta']), r)
                        for d, a, c, r in itertools.product(cfg['designs'], cfg['arrangements'],
                            cfg['cells'], range(cfg['rep_start'], cfg['rep_stop']))}
        actual_run = [identity(row) for row in source_rows]
        if len(actual_run) != len(expected_run) or set(actual_run) != expected_run:
            raise ValueError(f'Missing or repeated datasets: {run}')
        if seen.intersection(actual_run):
            raise ValueError(f'Overlapping dataset identities: {run}')
        seen.update(actual_run)
        records.extend(source_rows)
        for cell in sorted({ident[:-1] for ident in actual_run}):
            cells[cell] = dict(source_index=len(sources), datasets=sum(ident[:-1] == cell for ident in actual_run))
        try:
            location = str(run.relative_to(ROOT))
        except ValueError:
            location = str(run)
        sources.append(dict(run=location, datasets=len(source_rows), raw_sha256=sha(run/'raw.jsonl'),
                            metadata_sha256=sha(run/'metadata.json'), config_sha256=sha(run/'config.json')))
    if seen != expected or len(records) != len(expected):
        raise ValueError(f'Incomplete or unexpected grid: missing={len(expected-seen)}, extra={len(seen-expected)}')
    if len({r['data_sha256'] for r in records}) != len(records):
        raise ValueError('Repeated generated datasets in the complete grid.')
    methods = set(records[0]['pvalues'])
    if any(set(row['pvalues']) != methods for row in records):
        raise ValueError('Inconsistent method coverage.')
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=False)
    records.sort(key=identity)
    with (out/'raw.jsonl').open('x') as stream:
        for row in records:
            stream.write(json.dumps(row, sort_keys=True)+'\n')
    summarize(out, records)
    summaries = list(csv.DictReader((out/'summary.csv').open()))
    rates = {(s['design'], s['pattern'], float(s['theta']), s['method']): float(s['rate']) for s in summaries}
    ordered_cells = [(d, p, theta) for d in full['designs'] for p, theta in
                     [('null', 0.), ('dense', .5), ('dense', .75), ('dense', 1.),
                      ('sparse', .5), ('sparse', .75), ('sparse', 1.)]]
    with (out/'comparison.csv').open('x', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['design', 'pattern', 'theta', 'datasets', *METHODS])
        for d, p, theta in ordered_cells:
            writer.writerow([d, p, theta, 500, *[rates[d, p, theta, m] for m in METHODS]])
    lines = ['# Complete D.5 comparison', '',
             'Each row contains 500 paired datasets. The nominal level is 0.05. G-CRT uses 999 randomized copies and 3 sweeps.', '',
             '| Design | Pattern | θ | CoordG + Bonferroni | G-CRT + ridge | G-CRT + MaxPartial | G-CRT p-value combination |',
             '|---|---|---:|---:|---:|---:|---:|']
    for d, p, theta in ordered_cells:
        lines.append('| '+ ' | '.join([d, p, f'{theta:g}', *[f'{rates[d,p,theta,m]:.3f}' for m in METHODS[:4]]])+' |')
    lines += ['', 'The CSV files also retain oracle comparisons, MaxCorr diagnostics, exact binomial intervals, Monte Carlo standard errors and paired rejection differences.',
              r'The two G-CRT p-values are combined as $\min\{1,2\min(p_{\mathrm{ridge}},p_{\mathrm{MaxPartial}})\}$.',
              'Source runs and per-cell provenance are recorded in `provenance.json`.', '']
    (out/'RESULTS.md').write_text('\n'.join(lines))
    latex = ['% Generated by src/build_results.py; complete 14-cell grid.',
             '% Columns: design, pattern, theta, CoordG, GRidge, GMaxPartial, GAdaptiveBonf, OracleRidge, OracleMaxPartial, OracleAdaptiveBonf.']
    for d, p, theta in ordered_cells:
        latex.append(' & '.join([d, p, f'{theta:g}', *[f'{rates[d,p,theta,m]:.3f}' for m in METHODS]])+r' \\')
    (out/'table_rows.tex').write_text('\n'.join(latex)+'\n')
    manifest = dict(status='complete', datasets=len(records), method_records=len(records)*len(methods),
                    cells=len(cells), datasets_per_cell=500, sources=sources,
                    cell_sources=[dict(zip(KEYS[:-1], cell), **value) for cell, value in sorted(cells.items())],
                    full_config_sha256=sha(ROOT/'configs/full.json'), builder_sha256=sha(Path(__file__)))
    (out/'provenance.json').write_text(json.dumps(manifest, indent=2)+'\n')
    paths = sorted(p for p in out.rglob('*') if p.is_file())
    (out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.relative_to(out)}\n' for p in paths))
    print(json.dumps({k:manifest[k] for k in ('status', 'datasets', 'method_records', 'cells', 'datasets_per_cell')}, indent=2))


if __name__ == '__main__':
    main()
