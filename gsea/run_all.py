import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
TABLES = cfg.out_dir("gsea_tables")

CONTRASTS = [
    ("deseq2_naive_vs_stim", "GSEA_GO_BP_Naive_vs_Stim_results.csv"),
    ("deseq2_naive_vs_rested", "GSEA_GO_BP_Naive_vs_Rested_results.csv"),
]


def run(cmd):
    print("+", " ".join(str(c) for c in cmd))
    subprocess.run([str(c) for c in cmd], check=True)


for key, table in CONTRASTS:
    run(["Rscript", HERE / "00_run_gsea.R", cfg.path("rnaseq", key), TABLES / table])

run(["Rscript", HERE / "make_FigS2A.R",
     HERE / "data" / "GSEA_top5_IFNg.csv",
     cfg.out_dir("FigS2A") / "FigS2A.pdf"])

run([sys.executable, HERE / "make_FigS2E.py"])
