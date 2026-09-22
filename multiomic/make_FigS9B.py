"""Fig S9B: contribution of each of the 11 input features to PC1, PC2 and PC3. Three facets,
one horizontal lollipop per feature, markers coloured by sign, x limits derived from the data.

Output: results/FigS9B/FigS9B.{pdf,png,svg}
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

fs.setup()
np.random.seed(42)
DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("FigS9B")

load = pd.read_csv(DATA / "pc_loadings.tsv", sep="\t", index_col=0)
ve = pd.read_csv(DATA / "variance_explained.tsv", sep="\t", index_col=0)
pct = {pc: round(100 * float(ve.loc[pc, "variance_explained"])) for pc in ("PC1", "PC2", "PC3")}

# n the PCA was fit on -- reported, not drawn
n_genes = int(pd.read_csv(DATA / "feature_matrix.tsv", sep="\t", index_col=0).shape[0])

ROWS = [
    ("RNA_Stim_vs_Naive",                 "RNA induction (Stim)"),
    ("RNA_Rested_vs_Naive",               "RNA persistence (Rested)"),
    ("RNA_Rechallenged_vs_Stim",          "RNA priming (Rechal)"),
    ("K27me3_TSS3kb_Stim_vs_Naive",       "H3K27me3 (Stim)"),
    ("K27me3_TSS3kb_Rested_vs_Naive",     "H3K27me3 (Rested)"),
    ("K27me3_TSS3kb_Rechallenged_vs_Naive","H3K27me3 (Rechal)"),
    ("K4me3_TSS3kb_Stim_vs_Naive",        "H3K4me3 (Stim)"),
    ("K4me3_TSS3kb_Rested_vs_Naive",      "H3K4me3 (Rested)"),
    ("K4me3_TSS3kb_Rechallenged_vs_Naive","H3K4me3 (Rechal)"),
    ("K36me2_TSS3kb_Stim_vs_Naive",       "H3K36me2 (Stim)"),
    ("K36me3_geneBody_Stim_vs_Naive",     "H3K36me3 gene body (Stim)"),
]
POS, NEG = "#b2182b", "#2a81ba"     # red = positive loading, blue = negative (RdBu endpoints)

feats = [r[0] for r in ROWS]; labels = [r[1] for r in ROWS]
missing = [f for f in feats if f not in load.index]
assert not missing, f"features absent from pc_loadings.tsv: {missing}"
y = np.arange(len(ROWS))[::-1]

# x limits from the data, rounded out to the next 0.2, so the panel is not clipped
vmax = float(np.abs(load.loc[feats, ["PC1", "PC2", "PC3"]].values).max())
lim = np.ceil((vmax * 1.12) / 0.2) * 0.2
ticks = [t for t in np.arange(-1.2, 1.21, 0.4) if abs(t) <= lim + 1e-9]

fig, axes = plt.subplots(1, 3, figsize=fs.panel_mm(104, 58), sharey=True,
                         gridspec_kw=dict(wspace=0.16, left=0.305, right=0.99, top=0.86, bottom=0.15))
for ax, pc in zip(axes, ("PC1", "PC2", "PC3")):
    vals = load.loc[feats, pc].values
    ax.axvline(0, color="#BBBBBB", lw=0.7, zorder=1)
    for yi, v in zip(y, vals):
        c = POS if v >= 0 else NEG
        ax.plot([0, v], [yi, yi], color=c, lw=1.3, zorder=2, solid_capstyle="round")
        ax.plot(v, yi, "o", ms=3.8, color=c, mec="white", mew=0.4, zorder=3)
    fs.prism(ax, lw=1.0)
    ax.set_title(f"{pc} ({pct[pc]}% var)", fontsize=8, fontweight="bold", pad=3)
    ax.set_xlim(-lim, lim); ax.set_xticks(ticks)
    ax.tick_params(labelsize=6)
    ax.set_xlabel("loading", fontsize=8)
    ax.set_ylim(-0.6, len(ROWS) - 0.4)

axes[0].set_yticks(y)
axes[0].set_yticklabels(labels, fontsize=6)

fs.save(fig, OUT / "FigS9B")

print(f"PC loadings written | n = {n_genes} genes x {len(ROWS)} features")
print("var%:", pct, "| x limit:", round(lim, 2), "| max |loading|:", round(vmax, 3))
print(load.loc[feats, ["PC1", "PC2", "PC3"]].round(3).to_string())
pv = pd.read_csv(DATA / "per_view_pc_contribution.tsv", sep="\t", index_col=0)
print("\nper-view share of squared loading:")
print(pv.round(3).to_string())
