"""Copy an explicit completed run into a separate rebuild input set.
The bundled reference is immutable. INPUTS.json identifies every replacement;
files not listed there still come from the bundled reference (a mixed rebuild).
"""
from pathlib import Path
import argparse,json,shutil,hashlib,datetime
ROOT=Path(__file__).resolve().parents[1]
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('kind',choices=['gaussian-size','gaussian-power','graphical','matched','maxcoeft','direction','centering','connected','nongaussian','stock','stock-gof','stock-graph','bcr-observed','bcr-null'])
p.add_argument('run_dir',type=Path,help='Directory containing the completed result files, as printed by the study runner.')
p.add_argument('--output',default='outputs/fresh-reference')
a=p.parse_args();src=a.run_dir.resolve();dest=(ROOT/a.output).resolve()
if not dest.is_relative_to(ROOT/'outputs'):p.error('Output must be below outputs/')
if not src.is_dir():p.error('Run directory is missing')
if not dest.exists():
 shutil.copytree(ROOT/'reference',dest)
 (dest/'matched').mkdir();shutil.copyfile(ROOT/'experiments/d5-joint-comparison-20261007/results/complete/raw.jsonl',dest/'matched/raw.jsonl')
 (dest/'INPUTS.json').write_text('{}\n')
log=json.loads((dest/'INPUTS.json').read_text())
def put(source,target):
 if not source.is_file():raise FileNotFoundError(source)
 d=dest/target;d.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(source,d)
 log[target]=dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),kind=a.kind,
                 imported_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
if a.kind.startswith('gaussian-'):put(src/'results.csv','gaussian/'+('E01' if a.kind.endswith('size') else 'E02')+'/results.csv')
elif a.kind in ['graphical','stock']:
 import csv
 with (src/'results.csv').open() as f:ids={r['scenario'] for r in csv.DictReader(f)}
 if len(ids)!=1:raise ValueError('Choose one scenario directory')
 put(src/'results.csv',('graphical' if a.kind=='graphical' else 'stock')+'/'+ids.pop()+'.csv')
elif a.kind in ['matched','maxcoeft']:put(src/'raw.jsonl',a.kind+'/raw.jsonl')
elif a.kind=='direction':put(src/'direction-tests.csv','direction/direction-tests.csv')
elif a.kind=='centering':put(src/'tests.csv','centering/tests.csv')
elif a.kind=='connected':put(src/'tests.csv','connected/tests.csv')
elif a.kind=='nongaussian':
 files=sorted(src.glob('cell-*/replications.csv'))
 if len(files)!=12:raise ValueError('Choose results/formal with all 12 cells')
 for f in files:put(f,'nongaussian/'+f.parent.name+'.csv')
elif a.kind=='stock-gof':
 put(src/'result.rds','stock/gof.rds');put(src/'results.csv','stock/Stock_GoF.csv')
elif a.kind=='stock-graph':put(src/'graph.rds','stock/graph.rds')
elif a.kind.startswith('bcr-'):put(src/'results.csv','bcr/'+('observed.csv' if a.kind.endswith('observed') else 'UNRESTRICTED_RESULTS.csv'))
(dest/'INPUTS.json').write_text(json.dumps(log,indent=2)+'\n')
print('Rebuild inputs:',dest,'; replacements:',len(log),'; untouched files retain bundled provenance.')
