#!/usr/bin/env python3
"""Run the CRISPRa RNA-seq stages and draw the panels."""
import subprocess
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
STAGES = ["00_pydeseq2.py", "01_compute_tpm.py",
          "make_Fig2I.py", "make_Fig2K.py", "make_FigS5A.py"]

for stage in STAGES:
    print(f"== {stage}", flush=True)
    subprocess.run([sys.executable, str(HERE / stage)], check=True)
