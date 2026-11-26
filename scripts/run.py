#!/usr/bin/env python3
"""Explicit study dispatcher; full protocols remain in the sectioned source files."""
from pathlib import Path
import argparse
import datetime
import json
import os
import shutil
import subprocess
import sys

ROOT=Path(__file__).resolve().parents[1]


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('study',choices=['gaussian-size','gaussian-power','direction','centering','connected','nongaussian','matched','maxcoeft','bcr-observed','bcr-null','stock']+[f'E{i:02}' for i in range(3,16)]+['stock-graph','stock-gof'])
    parser.add_argument('--trial',action='store_true')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--scenario',default='stock_three_stat_replay_T1')
    parser.add_argument('--matched-raw',help='Fresh complete matched-study raw.jsonl, relative to the package.')
    args=parser.parse_args()
    assert 1<=args.workers<=6
    env=os.environ.copy()
    env.update(REPLICATION_ROOT=str(ROOT),OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',VECLIB_MAXIMUM_THREADS='1',MKL_NUM_THREADS='1',PYTHONDONTWRITEBYTECODE='1')
    env['REPLICATION_TRIAL']='1' if args.trial else '0'
    if args.matched_raw:
        source=(ROOT/args.matched_raw).resolve()
        if not source.is_relative_to(ROOT) or not source.is_file():
            parser.error('--matched-raw must name an existing file inside this package')
        env['REPLICATION_MATCHED_RAW']=str(source)
    stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
    out=ROOT/'outputs'/f'{args.study}-{stamp}'
    out.mkdir(parents=True)
    commands=[]
    def run(command,cwd=ROOT):
        commands.append(dict(argv=list(map(str,command)),cwd=str(Path(cwd).relative_to(ROOT))))
        (out/'commands.json').write_text(json.dumps(commands,indent=2)+'\n')
        print('RUN',subprocess.list2cmdline(list(map(str,command))),flush=True)
        with (out/'run.log').open('a') as log:
            subprocess.run(list(map(str,command)),cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
    def clone(name):
        target=out/name
        shutil.copytree(ROOT/'experiments'/name,target,ignore=shutil.ignore_patterns('outputs','runs','__pycache__','*.pyc'))
        return target
    if args.study.startswith('E') or args.study in ['stock-graph','stock-gof']:
        if args.study=='stock-gof' and args.trial:
            parser.error('Stock goodness of fit has a fixed 1,000-copy protocol; omit --trial.')
        entry={'stock-graph':'E19','stock-gof':'E20'}.get(args.study,args.study)
        run(['Rscript','--vanilla',ROOT/'experiments/readable-experiments/run.R',entry]+(['--trial'] if args.trial else []))
    elif args.study.startswith('gaussian-'):
        entry='E01' if args.study.endswith('size') else 'E02'
        d=clone('appendix-b-lw-20261001')
        run(['Rscript','--vanilla','-e',f'options(readable.trial={"TRUE" if args.trial else "FALSE"},lw.workers={args.workers}L); source("{entry}.R")'],d)
    elif args.study in ['direction','centering']:
        d=clone({'direction':'gcrt-direction-centering-m999-20261008','centering':'gcrt-symmetric-centering-20261004'}[args.study])
        # These inexpensive fixed studies have no reduced protocol.
        if args.trial: raise SystemExit('Use the full fixed protocol for direction/centering; see README runtime.')
        run(['Rscript','--vanilla','run.R'],d)
    elif args.study=='connected':
        d=clone('gcrt-centering-connected-v20')
        run(['Rscript','--vanilla','validate.R'],d)
        if not args.trial:
            run(['Rscript','--vanilla','run.R'],d)
            control=out/'controls';(control/'results').mkdir(parents=True)
            (control/'INPUTS.json').write_text(json.dumps(dict(source_directory=str(d))))
            run(['Rscript','--vanilla',ROOT/'experiments/gcrt-centering-connected-controls-20261009/analyse.R',control])
    elif args.study=='nongaussian':
        d=clone('mvn-crt-nongaussian-m3-20261018')
        for folder in ['logs','checks','results']: (d/folder).mkdir(exist_ok=True)
        for stage in ['prepare','check','pilot']+([] if args.trial else ['formal','summarize','verify']):
            run(['Rscript','--vanilla','run.R',stage],d)
    elif args.study=='matched':
        d=ROOT/'experiments/d5-joint-comparison-20261007'
        run([sys.executable,'-B',d/'src/run_confirmation.py','--config','configs/smoke.json' if args.trial else 'configs/full.json','--workers',args.workers],d)
    elif args.study=='maxcoeft':
        d=ROOT/'experiments/e5-max-coefficient-t-20261103'
        run([sys.executable,'-B',d/'src/run.py','--stage','pilot' if args.trial else 'full','--workers',args.workers,'--out','runs/'+stamp],d)
        records=[json.loads(p.read_text()) for p in sorted((d/'runs'/stamp/'records').glob('*.json'))]
        assert len(records)==(28 if args.trial else 7000)
        (d/'runs'/stamp/'raw.jsonl').write_text(''.join(json.dumps(r,sort_keys=True)+'\n' for r in records))
    elif args.study.startswith('bcr-'):
        run(['Rscript','--vanilla',ROOT/'scripts/bcr-adopted.R',args.study,out]+(['--trial'] if args.trial else []))
    elif args.study=='stock':
        d=ROOT/'experiments/stock-rf1000/engine'
        command=['Rscript','--vanilla',d/'scripts/run.R','--scenario',args.scenario,'--campaign','replication-'+stamp,'--workers','1','--backend','serial']
        if args.trial: command+=['--purpose','timing','--limit','1']
        run(command,d)
    (out/'COMPLETE.json').write_text(json.dumps(dict(study=args.study,trial=args.trial,commands=commands),indent=2)+'\n')
    print('Completed. Logs:',out)


if __name__=='__main__': main()
