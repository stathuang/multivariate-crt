"""Independent scalar checks of the unchanged handoff pilot and new runner."""
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))
import run_pilot as pilot
import geometry
from mscrt.core import RidgeStatistic, coordinate_pvalues, rank_pvalue, spherical_pvalue


def scalar_local(xt, x, graph, statistic):
    """Least-squares residuals, independent of the pilot's batched QR code."""
    n, t = xt.shape
    full = x.copy()
    full[:, :t] = xt
    local = []
    for j in range(t):
        a = np.column_stack((np.ones(n), statistic.h @ xt[:, np.arange(t) != j]))
        response = statistic.yres - a @ np.linalg.lstsq(a, statistic.yres, rcond=None)[0]
        d = np.column_stack((np.ones(n), full[:, np.flatnonzero(graph[j])]))
        r = xt[:, j] - d @ np.linalg.lstsq(d, xt[:, j], rcond=None)[0]
        u = response - d @ np.linalg.lstsq(d, response, rcond=None)[0]
        df = n - np.linalg.matrix_rank(d)
        den = (r @ r) * (u @ u)
        w = np.clip((r @ u)**2 / den, 0, 1) if den > 1e-26 else 0.
        local.append(beta.sf(w, .5, (df-1)/2) if df > 1 else 1.)
    return np.asarray(local)


class GeometryTests(unittest.TestCase):
    def fixture(self, seed=4107, n=24, p=9, t=3):
        rng = np.random.default_rng(seed)
        x = rng.normal(size=(n, p))
        y = rng.normal(size=n)
        graph = np.abs(np.arange(p)[:, None] - np.arange(p)) <= 2
        np.fill_diagonal(graph, False)
        stat = RidgeStatistic.fit(y, x[:, t:], .2*n, 1.)
        xt = np.concatenate((x[None, :, :t], rng.normal(size=(6, n, t))))
        return x, xt, graph, stat

    def test_scalar_agreement_and_retained_local_quantities(self):
        for seed in (7, 41, 4107):
            x, xt, graph, stat = self.fixture(seed)
            target = np.arange(xt.shape[-1])
            expected = np.array([scalar_local(z, x, graph, stat) for z in xt])
            actual_local = []
            sf = beta.sf
            def capture(w, a, b):
                result = sf(w, a, b)
                actual_local.append(result.copy())
                return result
            with patch.object(pilot.beta, 'sf', side_effect=capture):
                scores = pilot.maxpartial_scores(xt, x, target, graph, stat)
            np.testing.assert_allclose(np.array(actual_local).T, expected, atol=2e-12, rtol=2e-10)
            np.testing.assert_allclose(scores, -np.log(expected.min(axis=1)), atol=2e-12)
            # Compare every observed-coordinate tail with retained SVD implementation.
            for j in target:
                others = target[target != j]
                a = np.column_stack((np.ones(len(x)), stat.h @ x[:, others]))
                ry = stat.yres - a @ np.linalg.lstsq(a, stat.yres, rcond=None)[0]
                d = np.column_stack((np.ones(len(x)), x[:, np.flatnonzero(graph[j])]))
                self.assertAlmostEqual(expected[0, j], spherical_pvalue(x[:, j], d, ry), places=12)
            pg, _ = coordinate_pvalues(x, target, graph, stat)
            self.assertAlmostEqual(pg, min(1., len(target)*expected[0].min()), places=12)

    def test_permutation_batching_fixed_data_and_observed_target_independence(self):
        x, xt, graph, stat = self.fixture()
        target = np.arange(xt.shape[-1])
        backups = [a.copy() for a in (x, xt, graph, stat.h, stat.yres)]
        scores = pilot.maxpartial_scores(xt, x, target, graph, stat)
        order = np.array([3, 0, 6, 1, 5, 2, 4])
        np.testing.assert_allclose(pilot.maxpartial_scores(xt[order], x, target, graph, stat), scores[order], atol=1e-13)
        single = [pilot.maxpartial_scores(z, x, target, graph, stat)[0] for z in xt]
        np.testing.assert_allclose(single, scores, atol=1e-13)
        other_x = x.copy()
        other_x[:, target] = 1e8  # Candidate targets must replace these in neighbors.
        np.testing.assert_array_equal(pilot.maxpartial_scores(xt, other_x, target, graph, stat), scores)
        for actual, before in zip((x, xt, graph, stat.h, stat.yres), backups):
            np.testing.assert_array_equal(actual, before)

    def test_rank_ties(self):
        self.assertEqual(rank_pvalue(1., np.ones(19)), 1.)
        self.assertEqual(rank_pvalue(2., np.ones(19)), .05)
        self.assertEqual(rank_pvalue(1.+1e-14, np.ones(19)), 1.)

    def test_new_implementation_matches_pilot(self):
        for seed in (7, 41, 4107):
            x, xt, graph, stat = self.fixture(seed)
            target = np.arange(xt.shape[-1])
            expected = pilot.maxpartial_scores(xt, x, target, graph, stat)
            score, detail = geometry.maxpartial(xt, x, target, graph, stat, details=True)
            np.testing.assert_array_equal(score, expected)
            for b in range(len(xt)):
                np.testing.assert_allclose(detail['local_p'][b], scalar_local(xt[b], x, graph, stat), atol=2e-12)
            order = np.array([3, 0, 6, 1, 5, 2, 4])
            before=[a.copy() for a in (x,xt,graph,stat.h,stat.yres)]
            np.testing.assert_array_equal(geometry.maxpartial(xt[order],x,target,graph,stat),score[order])
            for a,b in zip((x,xt,graph,stat.h,stat.yres),before):
                np.testing.assert_array_equal(a,b)

    def test_rank_deficiency_and_zero_response(self):
        x, xt, graph, stat = self.fixture()
        xt[:, :, 1] = xt[:, :, 2]
        target=np.arange(xt.shape[-1])
        score, detail=geometry.maxpartial(xt,x,target,graph,stat,details=True)
        expected=np.array([scalar_local(z,x,graph,stat) for z in xt])
        np.testing.assert_allclose(detail['local_p'],expected,atol=2e-12)
        np.testing.assert_allclose(score,-np.log(expected.min(axis=1)),atol=2e-12)
        zero_stat=RidgeStatistic(stat.h,np.zeros_like(stat.yres),stat.penalty)
        np.testing.assert_array_equal(geometry.maxpartial(xt,x,target,graph,zero_stat),np.zeros(len(xt)))

    def test_two_test_correction(self):
        self.assertEqual(geometry.bonferroni_two(.04,.08),.08)
        self.assertEqual(geometry.bonferroni_two(.01,.08),.02)
        self.assertEqual(geometry.bonferroni_two(.6,.9),1.)
        # Enumerate exchangeable candidate positions; union bound holds for every
        # pair of score arrays including tied and perfectly dependent scores.
        for b in (np.arange(40), -np.arange(40), np.zeros(40)):
            a=np.arange(40)
            ps=[]
            for j in range(40):
                keep=np.arange(40)!=j
                ps.append(geometry.bonferroni_two(rank_pvalue(a[j],a[keep]),rank_pvalue(b[j],b[keep])))
            self.assertLessEqual(np.mean(np.array(ps)<=.05),.05)

    def test_smoke_both_designs_null_dense_sparse(self):
        for design in ('A', 'B'):
            result = pilot.run_design(design, 2, 9, 3, ['dense', 'sparse'], [.75])
            self.assertEqual(set(result['cells']), {'null:0', 'dense:0.75', 'sparse:0.75'})
            for cell in result['cells'].values():
                self.assertEqual(len(cell['methods']), 7)
                for method in cell['methods'].values():
                    self.assertEqual(method['reps'], 2)
                    self.assertTrue(0 <= method['rate'] <= 1)


if __name__ == '__main__':
    unittest.main(verbosity=2)
