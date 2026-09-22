"""Run the numbered rnaseq stages in order, then the panel scripts."""
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

STAGES = sorted(p.name for p in HERE.glob("[0-9][0-9]_*") if p.suffix in (".py", ".R"))
PANELS = sorted(p.name for p in HERE.glob("make_*") if p.suffix in (".py", ".R"))


def run(name):
    cmd = [sys.executable, name] if name.endswith(".py") else ["Rscript", name]
    print(f"\n=== {name} ===", flush=True)
    subprocess.run(cmd, cwd=HERE, check=True)


if __name__ == "__main__":
    for name in STAGES + PANELS:
        run(name)
