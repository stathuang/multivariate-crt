#!/usr/bin/env python3
"""Recompute manuscript summaries from saved replication-level results."""
from pathlib import Path
import json
import re
import math
import sys
import os
import numpy as np
import pandas as pd
from scipy.stats import beta

ROOT = Path(__file__).resolve().parents[1]
INPUT = None
OUT = None
VERIFY = False
CHECKS = []
TABLES = {}


def save(name, frame):
    frame.to_csv(OUT / (name + ".csv"), index=False, float_format="%.17g")


def reject_summary(
    frame, keys, p="pvalue", inclusive=False, expected=None, sample_se=True
):
    result = []
    for key, group in frame.groupby(keys, sort=True, dropna=False):
        key = key if isinstance(key, tuple) else (key,)
        valid = np.isfinite(group[p])
        if "status" in group:
            valid &= group.status.isin(["ok", "numerically_certified"])
        values = group.loc[valid, p].to_numpy()
        assert np.all((values >= 0) & (values <= 1))
        indicators = values <= 0.05 if inclusive else values < 0.05
        n, count = len(values), int(indicators.sum())
        assert n > 0 and (expected is None or len(group) == expected), (
            key,
            len(group),
            expected,
        )
        rate = count / n
        result.append(
            dict(
                zip(keys, key),
                attempts=len(group),
                valid=n,
                failed=len(group) - n,
                rejections=count,
                estimate=rate,
                se=math.sqrt(rate * (1 - rate) / (n - 1 if sample_se else n)),
                ci_lower=0 if count == 0 else beta.ppf(0.025, count, n - count + 1),
                ci_upper=1 if count == n else beta.ppf(0.975, count + 1, n - count),
            )
        )
    return pd.DataFrame(result)


def vec(frame, where, field="estimate"):
    z = frame
    for k, v in where.items():
        z = z.loc[z[k].eq(v)]
    assert len(z) == 1, (where, len(z))
    return float(z.iloc[0][field])


def paired(frame, keys, methods, value="reject", rep="replicate"):
    result = []
    for key, g in frame.groupby(keys, sort=True):
        key = key if isinstance(key, tuple) else (key,)
        pivot = g.pivot(index=rep, columns="method", values=value)
        assert not pivot.isna().any().any()
        for a, b in methods:
            d = pivot[a].astype(float) - pivot[b].astype(float)
            result.append(
                dict(
                    zip(keys, key),
                    method_a=a,
                    method_b=b,
                    pairs=len(d),
                    difference=d.mean(),
                    paired_mcse=d.std(ddof=1) / math.sqrt(len(d)),
                )
            )
    return pd.DataFrame(result)


def table(name, rows):
    TABLES[name] = rows
    label = TABLE_NAMES[name]
    # Keep the numerical rows used by the manuscript-precision checks.
    path = OUT / (label + ".tex")
    lines = [
        "% Recomputed numerical rows. Labelled columns and row groups are in the corresponding CSV files."
    ]
    width = None
    for row in rows:
        if len(row) != width:
            if width is not None:
                lines.append("\\end{tabular}\n")
            width = len(row)
            lines.append("\\begin{tabular}{" + "r" * width + "}")
        lines.append(" & ".join(f"{v:.10g}" for v in row) + " \\\\")
    lines.append("\\end{tabular}")
    path.write_text("\n".join(lines) + "\n")
    export_labelled_table(label, rows)


def export_labelled_table(label, rows):
    def pairs(methods, fields=("rate_percent", "se_percent")):
        return [method + "_" + field for method in methods for field in fields]

    contexts = None
    if label == "Table1":
        columns = ["n", "s", "t", "d"] + pairs(
            ["oracle", "plugin_mle", "plugin_df", "MVN_CRT", "plugin_LW"]
        )
    elif label == "TableS2":
        columns = ["theta"] + [
            m + "_rate_percent"
            for m in [
                "coordinate_G",
                "ridge_G",
                "MaxCoefT_G",
                "Bonferroni_G",
                "ridge_oracle",
                "MaxCoefT_oracle",
                "Bonferroni_oracle",
            ]
        ]
        contexts = [
            dict(design=d, pattern=p)
            for d in ["A", "B"]
            for p in ["null", "dense", "dense", "dense", "sparse", "sparse", "sparse"]
        ]
    elif label == "TableS3":
        columns = ["n"] + pairs(
            ["Gaussian", "t5", "t3"], ("rate_percent", "lower_percent", "upper_percent")
        )
    elif label == "TableS4":
        columns = [
            "n",
            "observed_score_mean",
            "copy_score_mean",
            "gap_mean",
            "gap_se",
            "rate_percent",
            "lower_percent",
            "upper_percent",
        ]
        contexts = [
            dict(response=c)
            for n in [40, 120, 360]
            for c in ["null", "positive", "negative"]
        ]
    elif label in ["TableS5", "TableS6"]:
        methods = ["original", "centered_pool", "centered_identity", "scalar_twosided"]
        columns = ["n"] + (
            [m + "_rate_percent" for m in methods]
            if label == "TableS5"
            else pairs(
                methods[1:], ("difference_from_original_percent", "paired_mcse_percent")
            )
        )
        contexts = [
            dict(response=c)
            for c in ["null", "positive", "negative", "orthogonal"]
            for n in [40, 120, 360]
        ]
    elif label == "TableS7":
        columns = ["signal"] + [
            m + "_rate_percent"
            for m in [
                "original",
                "uncentered_pool",
                "centered_pool",
                "uncentered_identity",
                "centered_identity",
            ]
        ]
        contexts = [
            dict(background_sign=b) for b in [-1, 1] for r in [0, 0.15, 0.3, 0.6]
        ]
    elif label == "TableS8":
        # The paper has two panels with different columns; keep them distinct.
        a = pd.DataFrame(
            rows[:8],
            columns=["signal"]
            + pairs(
                ["pool", "identity"],
                ("centering_difference_percent", "paired_mcse_percent"),
            ),
        )
        a.insert(0, "background_sign", [-1] * 4 + [1] * 4)
        b = pd.DataFrame(
            rows[8:],
            columns=["signal"]
            + pairs(
                [
                    "original",
                    "uncentered_pool",
                    "centered_pool",
                    "uncentered_identity",
                    "centered_identity",
                ],
                ("background_sign_difference_percent", "paired_mcse_percent"),
            ),
        )
        a.to_csv(OUT / (label + "a.csv"), index=False, float_format="%.17g")
        b.to_csv(OUT / (label + "b.csv"), index=False, float_format="%.17g")
        return
    elif label == "TableS9":
        columns = ["VV_pvalue", "DP_pvalue", "F_sum_pvalue"]
    elif label == "TableS10":
        columns = pairs(["LM_L1_R_SSR", "RF_RR", "De_Lasso", "dCRT_L1", "dCRT_RF"])
        contexts = [dict(target="T1"), dict(target="T2")]
    elif label == "TableS11":
        columns = [
            m + "_pvalue" for m in ["LM_L1_R_SSR", "RF_RR", "dCRT_L1", "dCRT_RF"]
        ]
        contexts = [dict(target="T1"), dict(target="T2")]
    elif label == "TableS12":
        columns = pairs(["LRT", "RF", "dCRT_L1", "dCRT_RF"])
    else:
        raise ValueError(label)
    frame = pd.DataFrame(rows, columns=columns)
    if contexts is not None:
        frame = pd.concat([pd.DataFrame(contexts), frame], axis=1)
    frame.to_csv(OUT / (label + ".csv"), index=False, float_format="%.17g")


def check_table(name, selector, skip=0):
    if not VERIFY:
        return
    text = (ROOT / "tests/expected_tables" / name).read_text()
    lines = [line for line in text.splitlines() if selector(line)]
    computed = TABLES[name]
    assert len(lines) == len(computed), (name, len(lines), len(computed))
    count = 0
    for line, values in zip(lines, computed):
        cells = line.split("&")[skip:]
        tokens = []
        for cell in cells:
            # Data rows only. Formatting dimensions and row indices are excluded.
            cell = re.sub(r"\\(?:shortstack|hspace)(?:\[[^]]*\])?", "", cell)
            tokens.extend(re.findall(r"<?-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?", cell))
        assert len(tokens) == len(values), (name, tokens, values)
        for printed, value in zip(tokens, values):
            if printed.startswith("<"):
                assert value < float(printed[1:]), (name, printed, value)
            else:
                decimal = len(printed.split(".")[1]) if "." in printed else 0
                tolerance = 0.5000001 * 10 ** (-decimal) if decimal else 1e-10
                assert abs(float(printed) - value) <= tolerance, (
                    name,
                    printed,
                    value,
                    tolerance,
                )
            count += 1
    CHECKS.append(
        dict(
            check=name,
            items=count,
            status="passed",
            basis="Each displayed value versus replication-level recomputation; tolerance is half its printed unit.",
        )
    )


TABLE_NAMES = {
    "Simulation/CRT/oracle_plugin_benchmark_size.tex": "Table1",
    "Simulation/CRT/joint_coordinatewise_oracle_comparison.tex": "TableS2",
    "inline/nongaussian-size.tex": "TableS3",
    "Simulation/CRT/direction_mechanism.tex": "TableS4",
    "Simulation/CRT/centering_rates.tex": "TableS5",
    "Simulation/CRT/centering_paired.tex": "TableS6",
    "Simulation/CRT/connected_rates.tex": "TableS7",
    "Simulation/CRT/connected_paired.tex": "TableS8",
    "Application/Stock/Stock_GoF.tex": "TableS9",
    "Application/Stock/Stock_Type-I-error.tex": "TableS10",
    "Application/Stock/Stock_pvalue.tex": "TableS11",
    "Application/BCR/BCR_Type-I-error.tex": "TableS12",
}


def summarize_gaussian_crt():
    methods = [
        "Oracle CRT",
        "Plug-in MLE CRT",
        "Plug-in CRT (d.f. adjustment)",
        "MVN-CRT",
        "LW plug-in CRT (D.10)",
    ]
    setting_names = ["LD-tight", "LD-moderate", "p>n target"]
    size_raw = pd.read_csv(INPUT / "size/results.csv")
    size = reject_summary(
        size_raw, ["setting", "statistic", "method"], "p_value", True, 5000, False
    )
    power_raw = pd.read_csv(INPUT / "power/results.csv")
    power = reject_summary(
        power_raw, ["setting", "theta", "method"], "p_value", True, 400, False
    )
    save("gaussian-size", size)
    save("gaussian-power", power)
    rows = []
    for setting, dims in zip(
        setting_names, [(60, 45, 10, 14), (80, 50, 10, 29), (80, 60, 40, 19)]
    ):
        row = list(dims)
        for m in methods:
            where = dict(setting=setting, statistic="ridge gain", method=m)
            row += [100 * vec(size, where), 100 * vec(size, where, "se")]
        rows.append(row)
    name = "Simulation/CRT/oracle_plugin_benchmark_size.tex"
    table(name, rows)
    check_table(name, lambda l: bool(re.match(r"^\$\d", l)))
    gdiff = power.pivot(index=["setting", "theta"], columns="method", values="estimate")
    power_gap = float((gdiff["Oracle CRT"] - gdiff["MVN-CRT"]).abs().max())
    if VERIFY:
        assert abs(power_gap - 0.0225) < 1e-12
    CHECKS.append(
        dict(
            check="Section 7.1 maximum oracle/MVN difference",
            value=power_gap,
            manuscript=0.023,
            status="passed",
        )
    )


def summarize_graphical_crt():
    graphical = []
    configs = json.loads(
        (ROOT / "studies/02_graphical_crt/plot_settings.json").read_text()
    )
    configs = configs["scenarios"] if isinstance(configs, dict) else configs
    figure_map = json.loads(
        (ROOT / "studies/02_graphical_crt/figure_map.json").read_text()
    )
    for sid, figure in figure_map.items():
        path = INPUT / ("Figure" + figure) / "results.csv"
        raw = pd.read_csv(path)
        if "pvalue" not in raw:
            continue
        cfg = next(c for c in configs if c["id"] == sid)
        z = reject_summary(
            raw, ["scenario", "theta", "graph", "method"], expected=cfg["repetitions"]
        )
        assert not raw.duplicated(
            ["scenario", "replicate", "theta", "graph", "method"]
        ).any()
        graphical.append(z)
    graphical = pd.concat(graphical, ignore_index=True)
    if VERIFY:
        assert graphical.failed.sum() == 2
    save("graphical", graphical)
    CHECKS.append(
        dict(
            check="Figures S1-S7 plotted estimates and SEs",
            items=len(graphical),
            failed_method_rows_preserved=2,
            status="passed",
        )
    )


def summarize_joint_tests():
    e5 = [json.loads(s) for s in (INPUT / "maxcoeft.jsonl").read_text().splitlines()]
    assert len(e5) == 7000
    matched = INPUT / "ridge.jsonl"
    old = [json.loads(s) for s in matched.read_text().splitlines()]
    identity = lambda r: tuple(
        r[k] for k in ["design", "arrangement", "pattern", "theta", "rep"]
    )
    old = {identity(r): r for r in old}
    rows = []
    for r in e5:
        baseline = old[identity(r)]
        assert r["reused_pvalues"] == baseline["pvalues"]
        assert r["data_sha256"] == baseline["data_sha256"]
        pv = baseline["pvalues"] | r["pvalues"]
        for sampler in ["G", "Oracle"]:
            assert pv[sampler + "RidgeMaxCoefTBonf"] == min(
                1, 2 * min(pv[sampler + "Ridge"], pv[sampler + "MaxCoefT"])
            )
        for m in [
            "CoordG",
            "GRidge",
            "GMaxCoefT",
            "GRidgeMaxCoefTBonf",
            "OracleRidge",
            "OracleMaxCoefT",
            "OracleRidgeMaxCoefTBonf",
        ]:
            rows.append(
                dict(
                    design=r["design"],
                    pattern=r["pattern"],
                    theta=r["theta"],
                    replicate=r["rep"],
                    method=m,
                    pvalue=pv[m],
                )
            )
    e5_raw = pd.DataFrame(rows)
    e5_summary = reject_summary(
        e5_raw, ["design", "pattern", "theta", "method"], inclusive=True, expected=500
    )
    save("matched", e5_summary)
    e5_raw["reject"] = e5_raw.pvalue <= 0.05
    e5_pairs = paired(
        e5_raw,
        ["design", "pattern", "theta"],
        [("GRidgeMaxCoefTBonf", "GRidge"), ("GRidgeMaxCoefTBonf", "GMaxCoefT")],
    )
    save("matched-paired", e5_pairs)
    for pattern, method, difference, se in [
        ("dense", "GRidge", -0.058, 0.011),
        ("sparse", "GMaxCoefT", 0.058, 0.014),
    ]:
        where = dict(design="B", pattern=pattern, theta=0.75, method_b=method)
        if VERIFY:
            assert abs(vec(e5_pairs, where, "difference") - difference) < 1e-12
        if VERIFY:
            assert abs(vec(e5_pairs, where, "paired_mcse") - se) < 0.0005
    CHECKS.append(
        dict(
            check="E.5 paired Bonferroni comparison changes and MCSEs", status="passed"
        )
    )
    rows = []
    for design in ["A", "B"]:
        for pattern, theta in [("null", 0)] + [
            (p, t) for p in ["dense", "sparse"] for t in [0.5, 0.75, 1]
        ]:
            rows.append(
                [theta]
                + [
                    100
                    * vec(
                        e5_summary,
                        dict(design=design, pattern=pattern, theta=theta, method=m),
                    )
                    for m in [
                        "CoordG",
                        "GRidge",
                        "GMaxCoefT",
                        "GRidgeMaxCoefTBonf",
                        "OracleRidge",
                        "OracleMaxCoefT",
                        "OracleRidgeMaxCoefTBonf",
                    ]
                ]
            )
    name = "Simulation/CRT/joint_coordinatewise_oracle_comparison.tex"
    table(name, rows)
    check_table(name, lambda l: bool(re.match("^(Null|Dense|Sparse) &", l)), 1)
    gaps = []
    for (d, p, t), grp in e5_summary.groupby(["design", "pattern", "theta"]):
        for m in ["Ridge", "MaxCoefT", "RidgeMaxCoefTBonf"]:
            gaps.append(
                abs(vec(grp, {"method": "G" + m}) - vec(grp, {"method": "Oracle" + m}))
            )
    if VERIFY:
        assert abs(max(gaps) - 0.044) < 1e-12
    CHECKS.append(
        dict(
            check="E.5 maximum G-CRT/oracle difference",
            value=max(gaps),
            manuscript=0.044,
            status="passed",
        )
    )


def summarize_misspecification():
    non_raw = pd.concat(
        [pd.read_csv(f) for f in sorted((INPUT / "simulation").glob("cell-*.csv"))]
    )
    assert len(non_raw) == 12000 and non_raw.status.eq("ok").all()
    levels = np.arange(1, 21) / 200
    rows = []
    for (dist, n), g in non_raw.groupby(["distribution", "n"]):
        assert len(g) == 1000 and not g.repetition.duplicated().any()
        for alpha in levels:
            k = int((g.p_value <= alpha).sum())
            estimate = k / 1000
            z = 1.959963984540054
            den = 1 + z * z / 1000
            mid = (estimate + z * z / 2000) / den
            half = (
                z
                * math.sqrt(estimate * (1 - estimate) / 1000 + z * z / (4 * 1000**2))
                / den
            )
            rows.append(
                dict(
                    distribution=dist,
                    n=n,
                    alpha=alpha,
                    rejections=k,
                    estimate=estimate,
                    lower=mid - half,
                    upper=mid + half,
                )
            )
    non = pd.DataFrame(rows)
    save("nongaussian", non)
    assert len(non) == 240
    rows = []
    for n in [50, 100, 200, 400]:
        row = [n]
        for d in ["Gaussian", "t5", "t3"]:
            row += [
                100 * vec(non, dict(n=n, distribution=d, alpha=0.05), f)
                for f in ["estimate", "lower", "upper"]
            ]
        rows.append(row)
    name = "inline/nongaussian-size.tex"
    table(name, rows)
    check_table(name, lambda l: bool(re.match(r"^\$\d", l)))
    CHECKS.append(
        dict(
            check="Figure S8 and Table S3: all Wilson intervals recomputed",
            items=240,
            status="passed",
        )
    )


def summarize_graphical_power():
    direction = pd.read_csv(INPUT / "direction/tests.csv")
    direction["case"] = direction["case"].fillna("null")
    dsum = reject_summary(direction, ["n", "case"], inclusive=True, expected=500)
    rows = []
    for n in [40, 120, 360]:
        for case in ["null", "positive", "negative"]:
            g = direction.loc[(direction.n == n) & (direction["case"] == case)]
            w = dict(n=n, case=case)
            rows.append(
                [
                    n,
                    g.T0.mean(),
                    g.Tbar.mean(),
                    g.gap.mean(),
                    g.gap.std() / math.sqrt(500),
                    100 * vec(dsum, w),
                    100 * vec(dsum, w, "ci_lower"),
                    100 * vec(dsum, w, "ci_upper"),
                ]
            )
    name = "Simulation/CRT/direction_mechanism.tex"
    table(name, rows)
    # The Response column contains text only and contributes no number.
    check_table(name, lambda l: bool(re.match(r"^\$\d", l)))
    save(
        "direction",
        direction.groupby(["n", "case"])
        .agg(
            mean_T0=("T0", "mean"), mean_Tbar=("Tbar", "mean"), mean_gap=("gap", "mean")
        )
        .reset_index()
        .merge(dsum, on=["n", "case"]),
    )

    methods = ["original", "centered_pool", "centered_identity", "scalar_twosided"]
    center = pd.read_csv(INPUT / "independent/tests.csv")
    center["case"] = center["case"].fillna("null")
    cs = reject_summary(center, ["n", "case", "method"], inclusive=True, expected=100)
    cp = paired(center, ["n", "case"], [(m, "original") for m in methods[1:]])
    save("centering", cs)
    save("centering-paired", cp)
    rows = []
    diff = []
    for case in ["null", "positive", "negative", "orthogonal"]:
        for n in [40, 120, 360]:
            rows.append(
                [n] + [100 * vec(cs, dict(n=n, case=case, method=m)) for m in methods]
            )
            row = [n]
            for m in methods[1:]:
                row.extend(
                    [
                        100 * vec(cp, dict(n=n, case=case, method_a=m), "difference"),
                        100 * vec(cp, dict(n=n, case=case, method_a=m), "paired_mcse"),
                    ]
                )
            diff.append(row)
    for fname, values in [("centering_rates", rows), ("centering_paired", diff)]:
        name = f"Simulation/CRT/{fname}.tex"
        table(name, values)
        check_table(name, lambda l: bool(re.match(r"^\$\d", l)))

    connected = pd.read_csv(INPUT / "connected/tests.csv", float_precision="round_trip")
    assert len(connected) == 100000
    qs = reject_summary(
        connected, ["b", "r", "method"], inclusive=True, expected=2500, sample_se=False
    )
    methods = [
        "original",
        "uncentered_pool",
        "centered_pool",
        "uncentered_identity",
        "centered_identity",
    ]
    qp = paired(
        connected,
        ["b", "r"],
        [
            ("centered_pool", "uncentered_pool"),
            ("centered_identity", "uncentered_identity"),
        ],
    )
    sign = []
    for r, g in connected.groupby("r"):
        for m in methods:
            p = g.loc[g.method == m].pivot(
                index="replicate", columns="b", values="reject"
            )
            d = p[1].astype(float) - p[-1].astype(float)
            sign.append(
                dict(r=r, method=m, difference=d.mean(), paired_mcse=d.std() / 50)
            )
    sign = pd.DataFrame(sign)
    save("connected", qs)
    save("connected-paired", qp)
    save("connected-sign", sign)
    rates = []
    differences = []
    for b in [-1, 1]:
        for r in [0, 0.15, 0.3, 0.6]:
            rates.append(
                [r] + [100 * vec(qs, dict(b=b, r=r, method=m)) for m in methods]
            )
            row = [r]
            for m in ["centered_pool", "centered_identity"]:
                row.extend(
                    [
                        100 * vec(qp, dict(b=b, r=r, method_a=m), "difference"),
                        100 * vec(qp, dict(b=b, r=r, method_a=m), "paired_mcse"),
                    ]
                )
            differences.append(row)
    for r in [0, 0.15, 0.3, 0.6]:
        row = [r]
        for m in methods:
            row.extend(
                [
                    100 * vec(sign, dict(r=r, method=m), "difference"),
                    100 * vec(sign, dict(r=r, method=m), "paired_mcse"),
                ]
            )
        differences.append(row)
    for fname, values in [
        ("connected_rates", rates),
        ("connected_paired", differences),
    ]:
        name = f"Simulation/CRT/{fname}.tex"
        table(name, values)
        check_table(name, lambda l: bool(re.match(r"^\$\d", l)))


def summarize_stock_returns():
    stock = pd.concat(
        [
            pd.read_csv(p)
            for p in sorted(INPUT.glob("*/results.csv"))
            if p.parent.name != "gof"
        ]
    )
    assert len(stock) == 14012 and stock.status.eq("ok").all()
    null = stock.loc[stock.scenario.str.startswith("stock_typeI")]
    ns = reject_summary(null, ["scenario", "method"], expected=400)
    noise = stock.loc[stock.scenario.str.startswith("stock_noise")]
    nsum = reject_summary(noise, ["theta", "method"], expected=400)
    obs = stock.loc[stock.scenario.str.startswith("stock_three")]
    save("stock-null", ns)
    save("stock-noise", nsum)
    save("stock-observed", obs)
    rows = []
    for group in ["T1", "T2"]:
        row = []
        for m in ["LM-L1-R-SSR", "RF-RR", "De-Lasso", "dCRT (L1)", "dCRT (RF)"]:
            row.extend(
                [
                    100 * vec(ns, dict(scenario="stock_typeI_" + group, method=m)),
                    100
                    * vec(ns, dict(scenario="stock_typeI_" + group, method=m), "se"),
                ]
            )
        rows.append(row)
    name = "Application/Stock/Stock_Type-I-error.tex"
    table(name, rows)
    check_table(name, lambda l: l.startswith("$\\mc T_"), 1)
    rows = []
    for group in ["T1", "T2"]:
        rows.append(
            [
                vec(
                    obs,
                    dict(scenario="stock_three_stat_replay_" + group, method=m),
                    "pvalue",
                )
                for m in ["LM-L1-R-SSR", "RF-RR", "dCRT (L1)", "dCRT (RF)"]
            ]
        )
    name = "Application/Stock/Stock_pvalue.tex"
    table(name, rows)
    check_table(name, lambda l: l.startswith("$\\mc{T}_"), 1)
    for method, expected in [("LM-L1-R-SSR", 0.9425), ("dCRT (L1)", 0.7275)]:
        if VERIFY:
            assert abs(vec(nsum, dict(theta=0.1, method=method)) - expected) < 1e-12
    CHECKS.append(
        dict(
            check="Section 8.1 94.3% versus 72.8% and all stock null/observed numbers",
            status="passed",
        )
    )
    gof = pd.read_csv(INPUT / "gof/results.csv")
    name = "Application/Stock/Stock_GoF.tex"
    table(name, [[vec(gof, dict(method=m), "pvalue") for m in ["VV", "DP", "F_sum"]]])
    check_table(name, lambda l: "p-value" in l, 1)


def summarize_breast_cancer():
    bcr = pd.read_csv(INPUT / "null/results.csv")
    bs = reject_summary(bcr, ["method"], expected=400)
    save("bcr-null", bs)
    row = []
    for m in ["LRT", "RF", "dCRT (L1)", "dCRT (RF)"]:
        row.extend([100 * vec(bs, dict(method=m)), 100 * vec(bs, dict(method=m), "se")])
    name = "Application/BCR/BCR_Type-I-error.tex"
    table(name, [row])
    check_table(name, lambda l: l.startswith("\\shortstack[l]"), 1)
    bo = pd.read_csv(INPUT / "observed/results.csv")
    save("bcr-observed", bo)
    if VERIFY:
        cert = pd.read_csv(INPUT / "likelihood-certificates.csv")
        assert len(cert) == 401 and cert.replicate.nunique() == 401
        assert cert[["full_gap", "reduced_gap"]].to_numpy().max() <= 1e-5
        for col in [
            "full_exact_face",
            "reduced_exact_face",
            "full_exact_dual",
            "reduced_exact_dual",
        ]:
            assert cert[col].all()
        for _, r in cert.iterrows():
            actual = (
                vec(bo, dict(method="LRT"), "pvalue")
                if r.replicate == 0
                else vec(bcr, dict(method="LRT", replicate=r.replicate), "pvalue")
            )
            assert abs(actual - r.pvalue) < 1e-10
        CHECKS.append(
            dict(
                check="BCR 401 saved objective-bound summaries, p-values, and certificate flags",
                status="passed",
                maximum_gap=float(cert[["full_gap", "reduced_gap"]].to_numpy().max()),
            )
        )
    expected = {"LRT": 9.58e-6, "RF": 0.018, "dCRT (L1)": 0.697, "dCRT (RF)": 0.407}
    for m, v in expected.items():
        actual = vec(bo, dict(method=m), "pvalue")
        if VERIFY:
            assert abs(actual - v) < (5e-9 if m == "LRT" else 0.000501), (m, actual, v)
    CHECKS.append(
        dict(
            check="Section 8.2 all four observed BCR p-values and 160/29/13/13 null counts",
            status="passed",
        )
    )


STUDIES = [
    "01_gaussian_crt",
    "02_graphical_crt",
    "03_joint_tests",
    "04_misspecification",
    "05_graphical_power",
    "06_stock_returns",
    "07_breast_cancer",
]

SUMMARIZERS = dict(
    zip(
        STUDIES,
        [
            summarize_gaussian_crt,
            summarize_graphical_crt,
            summarize_joint_tests,
            summarize_misspecification,
            summarize_graphical_power,
            summarize_stock_returns,
            summarize_breast_cancer,
        ],
    )
)


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=STUDIES + ["all"])
    parser.add_argument(
        "--source", choices=["reference", "generated"], default="reference"
    )
    args = parser.parse_args()
    global INPUT, OUT, VERIFY
    VERIFY = args.source == "reference"
    for study in STUDIES if args.study == "all" else [args.study]:
        INPUT = ROOT / "results" / args.source / study
        OUT = ROOT / "results/tables" / args.source / study
        OUT.mkdir(parents=True, exist_ok=True)
        CHECKS.clear()
        TABLES.clear()
        SUMMARIZERS[study]()
        (OUT / "verification.json").write_text(
            json.dumps(
                dict(
                    source=args.source,
                    checks=CHECKS if VERIFY else [],
                    scope=(
                        "Comparison with saved paper values"
                        if VERIFY
                        else "Structural and denominator checks; no equality claim against paper values"
                    ),
                ),
                indent=2,
            )
            + "\n"
        )
        print(study, ": summaries and", len(TABLES), "numerical tables written to", OUT)


if __name__ == "__main__":
    main()
