#!/usr/bin/env python3
"""Run the four numbered stages, then every panel script, in order."""
import subprocess
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent

STAGES = ["00_build_universe.py",
          "01_resolve_coordinates.py",
          "02_build_feature_matrix.py",
          "03_cluster.py"]

PANELS = ["make_pc_scatters.py",        # Fig 5B, Fig 5D, Fig S9C, Fig S10B
          "make_Fig5C_bars.py",
          "make_Fig5C_pies.py",         # reads the bar geometry, so it runs after the bars
          "make_Fig5A_FigS9A.py",       # matches the bar axes width, so it runs after the bars
          "make_Fig5E.py",
          "make_Fig5F.py",
          "make_FigS9B.py",
          "make_FigS10A.py",
          "make_FigS10C.py",
          "make_FigS10D.py"]

for script in STAGES + PANELS:
    print(f"\n=== {script} ===", flush=True)
    subprocess.run([sys.executable, str(HERE / script)], check=True)
