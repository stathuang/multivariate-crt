"""Read-only recomputation of saved runs, with a separate verification output."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import numpy as np
from scipy.stats import binomtest

PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def source_hashes(meta):
    """Resolve source evidence inside this package, including historical snapshots."""
    if 'source_relative_sha256' in meta:
        return {PACKAGE_ROOT / name: value for name, value in meta['source_relative_sha256'].items()}
    provenance = json.loads((PACKAGE_ROOT/'provenance/SOURCES.json').read_text())
    snapshots = {item['original_path']: PACKAGE_ROOT/item['snapshot'] for item in provenance['sources']}
    return {snapshots[name]: value for name, value in meta['source_sha256'].items()}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('run')
    ap.add_argument('--out',required=True)
    ap.add_argument('--allow-test-source-changes',action='store_true',help='Only for smoke runs preceding expanded tests; implementation hashes remain strict.')
    args=ap.parse_args()
    run=Path(args.run)
    cfg=json.loads((run/'config.json').read_text())
    meta=json.loads((run/'metadata.json').read_text())
    assert meta['status']=='complete'
    checked_hashes=0
    for line in (run/'SHA256SUMS').read_text().splitlines():
        expected,name=line.split('  ',1)
        assert hashlib.sha256((run/name).read_bytes()).hexdigest()==expected,name
        checked_hashes+=1
    changed_sources=[]
    for path,expected in source_hashes(meta).items():
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            changed_sources.append(str(path))
            assert args.allow_test_source_changes and '/tests/test_' in str(path),path
    rows=[json.loads(line) for line in (run/'raw.jsonl').read_text().splitlines()]
    ident=lambda r:(r['design'],r['arrangement'],r['pattern'],float(r['theta']),int(r['rep']))
    expected={(d,g,c['pattern'],float(c['theta']),r) for d,g,c,r in itertools.product(cfg['designs'],cfg['arrangements'],cfg['cells'],range(cfg['rep_start'],cfg['rep_stop']))}
    assert len(rows)==len(expected)==len({ident(r) for r in rows})
    assert {ident(r) for r in rows}==expected
    groups={}
    coord_errors=[]
    for row in rows:
        groups.setdefault(ident(row)[:-1],[]).append(row)
        p=row['pvalues']
        assert all(0<=v<=1 and math.isfinite(v) for v in p.values())
        assert row['reject']=={k:int(v<=cfg['alpha']) for k,v in p.items()}
        coord_errors.append(abs(p['CoordG']-min(1.,cfg['t']*min(row['observed_local_p']))))
        assert coord_errors[-1]<2e-11
        for sampler in ('Oracle','G'):
            assert p[sampler+'AdaptiveBonf']==min(1.,2*min(p[sampler+'Ridge'],p[sampler+'MaxPartial']))
            for score in ('Ridge','MaxPartial','MaxCorr'):
                diag=row['score_diagnostics'][sampler+score]
                assert p[sampler+score]==(1+diag['greater']+diag['ties'])/(cfg['M']+1)
        for score in ('Ridge','MaxPartial','MaxCorr'):
            assert row['score_diagnostics']['G'+score]['observed']==row['score_diagnostics']['Oracle'+score]['observed']
        for comparison,value in row['paired_reject_differences'].items():
            a,b=comparison.split(' minus ')
            assert value==row['reject'][a]-row['reject'][b]
    maximum_summary_error=0.
    summaries=list(csv.DictReader((run/'summary.csv').open()))
    expected_summaries = {(*key, method) for key, records in groups.items() for method in records[0]['pvalues']}
    actual_summaries = [(s['design'],s['arrangement'],s['pattern'],float(s['theta']),s['method']) for s in summaries]
    assert len(actual_summaries) == len(set(actual_summaries)) == len(expected_summaries)
    assert set(actual_summaries) == expected_summaries
    for s in summaries:
        g=groups[(s['design'],s['arrangement'],s['pattern'],float(s['theta']))]
        n=len(g);k=sum(r['reject'][s['method']] for r in g)
        ci=binomtest(k,n).proportion_ci(confidence_level=.95,method='exact')
        checks={'reps':n,'rejections':k,'rate':k/n,'mcse':math.sqrt((k-k*k/n)/(n*(n-1))),
                'exact_low':ci.low,'exact_high':ci.high}
        for key,value in checks.items():
            error=abs(float(s[key])-value)
            maximum_summary_error=max(maximum_summary_error,error)
            assert error<2e-10,(s,key,error)
    pairs=list(csv.DictReader((run/'paired_differences.csv').open()))
    expected_pairs = {(*key, *comparison.split(' minus ')) for key, records in groups.items() for comparison in records[0]['paired_reject_differences']}
    actual_pairs = [(p['design'],p['arrangement'],p['pattern'],float(p['theta']),p['method'],p['reference']) for p in pairs]
    assert len(actual_pairs) == len(set(actual_pairs)) == len(expected_pairs)
    assert set(actual_pairs) == expected_pairs
    for pair in pairs:
        g=groups[(pair['design'],pair['arrangement'],pair['pattern'],float(pair['theta']))]
        d=[r['reject'][pair['method']]-r['reject'][pair['reference']] for r in g]
        n=len(d);total=sum(d);squares=sum(v*v for v in d)
        assert int(pair['reps'])==n
        assert abs(float(pair['difference'])-total/n)<1e-14
        assert abs(float(pair['paired_mcse'])-math.sqrt((squares-total*total/n)/(n*(n-1))))<1e-13
        assert int(pair['wins'])==d.count(1) and int(pair['losses'])==d.count(-1)
    score_files=0
    for row in rows:
        if row['rep']!=cfg['rep_start']:
            continue
        name=f"{row['design']}-{row['arrangement']}-{row['pattern']}-{row['theta']:g}-r{row['rep']}.npz"
        with np.load(run/'scores'/name) as scores:
            for method in scores.files:
                values=np.round(scores[method],12)
                assert len(values)==cfg['M']+1
                assert row['pvalues'][method]==np.mean(values>=values[0])
        score_files+=1
    result=dict(status='passed',run=str(run.resolve()),datasets=len(rows),method_records=sum(len(r['pvalues']) for r in rows),
        distinct_cells=len(groups),marginal_summary_rows=len(summaries),paired_summary_rows=len(pairs),
        saved_candidate_score_files=score_files,checked_output_hashes=checked_hashes,
        source_hashes_unchanged=not changed_sources,changed_test_sources=changed_sources,
        implementation_hashes_unchanged=True,maximum_coordinate_tail_error=max(coord_errors),maximum_summary_error=maximum_summary_error)
    assert meta['datasets'] == len(rows)
    assert meta['method_records'] == result['method_records']
    result['source_verification_basis'] = 'package_relative_sources' if 'source_relative_sha256' in meta else 'preserved_original_source_snapshots'
    Path(args.out).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    main()
