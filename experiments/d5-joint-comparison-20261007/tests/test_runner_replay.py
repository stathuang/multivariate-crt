"""The persisted runner must reproduce original-pilot method decisions."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key]='1'
import sys
sys.dont_write_bytecode=True
import unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
import run_confirmation as runner


class ReplayTests(unittest.TestCase):
    def test_unchanged_pilot_decisions_in_both_designs(self):
        cfg=dict(n=80,p=120,t=8,M=199,L=3,nuisance_penalty=.2,target_penalty=1.,alpha=.05,rep_start=0)
        mapping={'CoordG':'CoordG-Bonf','OracleRidge':'OracleRidge','GRidge':'GRidge-L',
            'OracleMaxCorr':'OracleMaxCorr','GMaxCorr':'GMaxCorr-L',
            'OracleMaxPartial':'OracleMaxPartial','GMaxPartial':'GMaxPartial-L'}
        for design in ('A','B'):
            original=runner.reference.run_design(design,2,199,3,['dense','sparse'],[.75])
            for pattern,theta in [('null',0.),('dense',.75),('sparse',.75)]:
                records=[runner.one_dataset((cfg,design,'connected',pattern,theta,rep))[0] for rep in range(2)]
                for method,old in mapping.items():
                    self.assertEqual(sum(r['reject'][method] for r in records),
                        original['cells'][f'{pattern}:{theta:g}']['methods'][old]['rejections'])
                # Reordering sampler-specific computations cannot affect their
                # data stream: deterministic rerun gives identical data and p-values.
                replay=runner.one_dataset((cfg,design,'connected',pattern,theta,0))[0]
                self.assertEqual(replay['data_sha256'],records[0]['data_sha256'])
                self.assertEqual(replay['pvalues'],records[0]['pvalues'])


if __name__=='__main__':
    unittest.main(verbosity=2)
