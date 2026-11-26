"""Check local imports and exact replay against already saved formal data."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import sys
sys.dont_write_bytecode = True
import json
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))
import run_confirmation as runner
import mscrt.core
import mscrt.experiments


class PortabilityTests(unittest.TestCase):
    def test_all_implementation_imports_are_local(self):
        for module in [runner, runner.geometry, runner.reference, mscrt.core, mscrt.experiments]:
            self.assertTrue(Path(module.__file__).resolve().is_relative_to(ROOT))

    def test_exact_replay_of_saved_500_replication_confirmation(self):
        baseline = ROOT/'reference_results/initial-confirmation'
        cfg = json.loads((baseline/'config.json').read_text())
        wanted = {('A','dense',0.75,10000), ('A','sparse',0.75,10000),
                  ('B','dense',0.75,10000), ('B','sparse',0.75,10000)}
        checked = set()
        with (baseline/'raw.jsonl').open() as stream:
            for line in stream:
                old = json.loads(line)
                key = (old['design'],old['pattern'],old['theta'],old['rep'])
                if key not in wanted:
                    continue
                new, _ = runner.one_dataset((cfg,old['design'],old['arrangement'],old['pattern'],old['theta'],old['rep']))
                for field in ['data_sha256','pvalues','reject','observed_local_p','observed_df',
                              'score_diagnostics','conditional_effects','paired_reject_differences']:
                    self.assertEqual(new[field],old[field],(key,field))
                checked.add(key)
        self.assertEqual(checked,wanted)

    def test_supplement_identity_is_disjoint_from_existing_grid(self):
        base = json.loads((ROOT/'reference_results/initial-confirmation/config.json').read_text())
        added = json.loads((ROOT/'configs/sparse-completion.json').read_text())
        identity = lambda c: {(d,a,x['pattern'],float(x['theta']),r)
            for d in c['designs'] for a in c['arrangements'] for x in c['cells']
            for r in range(c['rep_start'],c['rep_stop'])}
        old, new = identity(base), identity(added)
        self.assertEqual(len(new),2000)
        self.assertFalse(old & new)
        self.assertEqual(len(old | new),7000)


if __name__ == '__main__':
    unittest.main(verbosity=2)
