#!/usr/bin/env python3
"""Run every calculation in one paper study, or the complete paper."""
import argparse, json, os, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
commands = json.loads((ROOT / "studies/commands.json").read_text())
p = argparse.ArgumentParser(description=__doc__)
p.add_argument("study", choices=list(commands) + ["all"])
p.add_argument(
    "--quick",
    action="store_true",
    help="Use fewer datasets, retaining formal method/copy settings",
)
a = p.parse_args()
env = dict(os.environ)
for key in [
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "MKL_NUM_THREADS",
]:
    env[key] = "1"
for study in commands if a.study == "all" else [a.study]:
    for command in commands[study]:
        command = [sys.executable if item == "PYTHON" else item for item in command]
        if a.quick:
            command += ["--quick"]
        print(" ".join(command), flush=True)
        subprocess.run(command, cwd=ROOT, env=env, check=True)
    if not a.quick:
        subprocess.run(
            [sys.executable, "scripts/summarize.py", study, "--source", "generated"],
            cwd=ROOT,
            env=env,
            check=True,
        )
        subprocess.run(
            ["Rscript", "--vanilla", "scripts/figures.R", study, "generated"],
            cwd=ROOT,
            env=env,
            check=True,
        )
