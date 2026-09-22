#!/usr/bin/env python3
"""Run the accessibility stages and panels in order."""
import subprocess
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
HERE = Path(__file__).resolve().parent

SCRIPTS = [
    "00_msp_summary.py",
    "01_msp_area_stats.py",
    "02_tss_accessibility.py",
    "make_Fig4A.py",
    "make_Fig4B_4C.py",
]

for name in SCRIPTS:
    print(f"=== {name}")
    subprocess.run([sys.executable, str(HERE / name)], check=True)
