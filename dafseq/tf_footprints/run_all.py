#!/usr/bin/env python3
"""Run the four TF-footprint panels."""
import subprocess
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent

RUNS = [
    ["make_FigS8A.py"],
    ["make_FigS8B.py", "--extra-plot-tf", "IRF1", "--extra-plot-threshold", "0.70"],
    ["make_FigS8C.py"],
    ["make_FigS8D.py"],
]

for run in RUNS:
    print(f"=== {run[0]}")
    subprocess.run([sys.executable, str(HERE / run[0]), *run[1:]], check=True)
