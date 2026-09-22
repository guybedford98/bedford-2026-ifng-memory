"""Fig S6D: CXCL9 and CXCL10 TPM across the four IFNg conditions.

Drawn together with CXCL11 and IDO1 so the four panels form one matched set: the
H3K27me3 derepression domain spans CXCL11 together with its adjacent, co-regulated
neighbours CXCL9 and CXCL10, and all three are potently IFNg-induced memory genes.

That claim is verified rather than assumed. All three are checked against the
protein-coding positive-memory list at run time and the script aborts if any of them is
missing, so the sentence cannot drift away from the gene list. Induction and potentiation
are annotated from the DESeq2 tables, not from the TPM ratios: Naive->Stim for induction,
Stim->Rech for the memory bracket.

Per-replicate TPM, S1-S3 Naive, S4-S6 Stimulated, S7-S9 Rested, S10-S12 Rechallenged.
"""
import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs
import plotstats as ps

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

OUT = cfg.out_dir("FigS6D")
MEMLIST = cfg.out_dir("rnaseq_tables") / "Positive_memory_genes_maxTPM5_proteincoding.csv"

CONDS = ["Naive", "Stim", "Rested", "Rech"]
COLORS = [fs.COLOR["naive"], fs.COLOR["stimulated"],
          fs.COLOR["rested"], fs.COLOR["rechallenge"]]

DATA = {
    "CXCL9":  {"Naive": [0.0, 0.05, 0.02], "Stim": [16.95, 15.64, 14.77],
               "Rested": [0.61, 0.58, 0.52], "Rech": [59.52, 45.81, 80.77]},
    "CXCL10": {"Naive": [0.0, 0.56, 1.01], "Stim": [295.77, 342.85, 344.73],
               "Rested": [3.12, 4.42, 6.22], "Rech": [595.03, 773.34, 574.24]},
    "CXCL11": {"Naive": [0.0, 1.59, 2.00], "Stim": [593.10, 641.29, 637.17],
               "Rested": [5.99, 10.85, 17.46], "Rech": [1413.72, 1676.76, 1411.91]},
    "IDO1":   {"Naive": [0.0, 1.58, 1.77], "Stim": [512.55, 491.75, 512.80],
               "Rested": [4.49, 9.78, 11.94], "Rech": [1479.19, 1619.31, 1513.05]},
}

DESEQ = {                       # gene: (padj Naive->Stim, padj Stim->Rech)
    "CXCL9":  (6.470866e-28, 1.45e-14),
    "CXCL10": (1.750347e-20, 3.31e-19),
    "CXCL11": (4.132650e-10, 2.14e-56),
    "IDO1":   (4.519928e-09, 2.22e-126),
}

YAXIS = {                       # gene: (ymax, ticks)
    "CXCL9":  (100, [0, 25, 50, 75, 100]),
    "CXCL10": (1000, [0, 250, 500, 750, 1000]),   # headroom: top replicate is 773
    "CXCL11": (2000, [0, 500, 1000, 1500, 2000]),
    "IDO1":   (2000, [0, 500, 1000, 1500, 2000]),
}

ORDER = ["CXCL9", "CXCL10", "CXCL11", "IDO1"]

# CXCL9 and CXCL10 are the Fig S6D panel; CXCL11 and IDO1 are the matched references.
STEM = {"CXCL9": "FigS6D_CXCL9", "CXCL10": "FigS6D_CXCL10",
        "CXCL11": "TPM_CXCL11", "IDO1": "TPM_IDO1"}

mem = set(pd.read_csv(MEMLIST)["gene_id"])
missing = [g for g in ("CXCL9", "CXCL10", "CXCL11") if g not in mem]
if missing:
    raise SystemExit(
        "these genes are NOT in the protein-coding positive-memory list, so the "
        "manuscript sentence is wrong: %s" % ", ".join(missing))
print("memory-list check: CXCL9, CXCL10 and CXCL11 are all among the %d "
      "protein-coding memory genes" % len(mem))

fs.setup()


def panel(gene, ax):
    d = DATA[gene]
    series = [d[c] for c in CONDS]
    info = fs.bars(ax, series, colors=COLORS)
    ymax, ticks = YAXIS[gene]
    ax.set_ylim(0, ymax)
    ax.set_yticks(ticks)
    ax.set_ylabel("TPM", fontweight="bold")
    ax.set_xticks(range(len(CONDS)))
    ax.set_xticklabels(CONDS, rotation=45, ha="right", fontweight="bold")
    ax.set_title(gene, fontsize=8, fontstyle="italic")

    # Memory bracket: Stimulated vs Rechallenged, asterisks from DESeq2.
    p_stim, p_mem = DESEQ[gene]
    fs.sig_bracket(ax, 1, 3, info["ytop"][3] + ymax * 0.06, ps.pstars(p_mem))
    fs.prism(ax)
    print("  %-7s Stim %7.1f -> Rech %7.1f TPM   induction padj %.2g   memory padj %.2g %s"
          % (gene, info["means"][1], info["means"][3], p_stim, p_mem, ps.pstars(p_mem)))
    return info


print()
for gene in ORDER:
    fig, ax = plt.subplots(figsize=fs.panel_mm(32, 36))
    panel(gene, ax)
    fs.save(fig, OUT / STEM[gene])

fig, axes = plt.subplots(1, len(ORDER), figsize=fs.panel_mm(32 * len(ORDER), 38))
for ax, gene in zip(axes, ORDER):
    panel(gene, ax)
fig.tight_layout(w_pad=1.6)
fs.save(fig, OUT / "FigS6D_CXCL_cluster_row")

print("\nfigures in:", OUT)
for p in sorted(OUT.glob("*.png")):
    print("  ", p.name)
