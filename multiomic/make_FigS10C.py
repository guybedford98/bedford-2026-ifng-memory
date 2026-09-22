"""Fig S10C: per-cluster mean log2FC heatmap over the 11 features (chromatin + RNA).

Output: results/FigS10C/FigS10C.{pdf,png,svg}, data/cluster_feature_means.tsv
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

fs.setup()
np.random.seed(42)
DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("FigS10C")

CLUSTER_COLORS = {1: "#20566b", 2: "#c9a227", 3: "#2e8b57", 4: "#b0532b"}
COND_COLORS = {"Stim": fs.COLOR["stimulated"], "Rested": fs.COLOR["rested"],
               "Rechallenged": fs.COLOR["rechallenge"]}
CELL_FONTSIZE = 8            # in-square numbers

feat = pd.read_csv(DATA / "feature_matrix.tsv", sep="\t").set_index("gene")
cl = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t").set_index("gene")
n_universe = sum(1 for _ in open(DATA / "universe.tsv", encoding="utf-8")) - 1

FEATURES = [
    ("RNA_Stim_vs_Naive",                  "RNA",                "Stim"),
    ("RNA_Rested_vs_Naive",                "RNA",                "Rested"),
    ("RNA_Rechallenged_vs_Stim",           "RNA",                "Rechallenged"),
    ("K27me3_TSS3kb_Stim_vs_Naive",        "H3K27me3",           "Stim"),
    ("K27me3_TSS3kb_Rested_vs_Naive",      "H3K27me3",           "Rested"),
    ("K27me3_TSS3kb_Rechallenged_vs_Naive","H3K27me3",           "Rechallenged"),
    ("K4me3_TSS3kb_Stim_vs_Naive",         "H3K4me3",            "Stim"),
    ("K4me3_TSS3kb_Rested_vs_Naive",       "H3K4me3",            "Rested"),
    ("K4me3_TSS3kb_Rechallenged_vs_Naive", "H3K4me3",            "Rechallenged"),
    ("K36me2_TSS3kb_Stim_vs_Naive",        "H3K36me2",           "Stim"),
    ("K36me3_geneBody_Stim_vs_Naive",      "H3K36me3\ngene body","Stim"),
]
cols = [c for c, _, _ in FEATURES]; marks = [m for _, m, _ in FEATURES]; conds = [d for _, _, d in FEATURES]
is_rna = np.array([m == "RNA" for m in marks])

df = feat.loc[feat.index.intersection(cl.index), cols].join(cl["cluster"])
CL = sorted(df["cluster"].unique(), key=lambda x: int(x.split("_")[1]))
n_by = df.groupby("cluster").size()
H = df.groupby("cluster")[cols].mean().loc[CL]
V = H.values

# colour scale from chromatin cells only (preserve PTM contrast); RNA clips
chrom_vmax = float(np.nanmax(np.abs(V[:, ~is_rna])))
norm = TwoSlopeNorm(vmin=-chrom_vmax, vcenter=0.0, vmax=chrom_vmax)

# every column the same width, one white separator between each mark group
POS = list(range(len(FEATURES))); NDISP = len(FEATURES)
SEPS = [2.5, 5.5, 8.5, 9.5]

D = np.full((len(CL), NDISP), np.nan)
for jc, jp in enumerate(POS):
    D[:, jp] = V[:, jc]
Dm = np.ma.masked_invalid(D)
cmap = plt.get_cmap("RdBu_r").copy(); cmap.set_bad("white", 0.0)

# 170 mm wide so the 8 pt numbers clear their cell walls
fig, ax = plt.subplots(figsize=fs.panel_mm(170, 48))
ax.set_facecolor("white")
im = ax.imshow(Dm, cmap=cmap, norm=norm, aspect="auto")
for i in range(len(CL)):
    for jc, jp in enumerate(POS):
        v = V[i, jc]
        ax.text(jp, i, f"{v:+.2f}", ha="center", va="center", fontsize=CELL_FONTSIZE,
                color="white" if abs(v) > 0.62 * chrom_vmax else "black")

ax.set_yticks(range(len(CL)))
ax.set_yticklabels([f"{fs.clabel(cc)}  (n={int(n_by[cc])})" for cc in CL], fontsize=6)
for tick, cc in zip(ax.get_yticklabels(), CL):
    tick.set_color(CLUSTER_COLORS[int(fs.clabel(cc))]); tick.set_fontweight("bold")
ax.set_ylabel("cluster", fontsize=8, fontweight="bold")

short = {"Stim": "Stim", "Rested": "Rest", "Rechallenged": "Rech"}
ax.set_xticks(POS)
ax.set_xticklabels([short[c] for c in conds], rotation=45, ha="right", fontsize=6)
for tick, c in zip(ax.get_xticklabels(), conds):
    tick.set_color(COND_COLORS[c])
for xb in SEPS:
    ax.axvline(xb, color="white", lw=1.8)

groups = {}
for jc, m in enumerate(marks):
    groups.setdefault(m, []).append(POS[jc])
for m, ps in groups.items():
    # H3K36me2 and H3K36me3 are single-column groups sitting side by side, so their headers
    # get slightly less room than the three-column groups
    ax.text(float(np.mean(ps)), -0.72, m, ha="center", va="bottom",
            fontsize=8, fontweight="bold", color="black")

ax.tick_params(length=0)
for s in ("top", "right", "left", "bottom"):
    ax.spines[s].set_visible(False)
ax.set_xlim(-0.5, NDISP - 0.5)
cb = fig.colorbar(im, ax=ax, fraction=0.018, pad=0.012, extend="both")
cb.ax.tick_params(labelsize=6, length=2.0, width=0.8)
cb.outline.set_linewidth(0.6)
cb.set_label("mean log2FC vs naive\n(RNA Rechal = vs Stim)", fontsize=8)

n_clustered = int(n_by.sum())
ax.text(-0.5, len(CL) + 0.22,
        f"n = genes in the cluster, {n_clustered:,} of {n_universe:,} universe genes clustered",
        ha="left", va="top", fontsize=8, color="#555", clip_on=False)

fig.subplots_adjust(top=0.76)
fs.save(fig, OUT / "FigS10C")
H.round(3).to_csv(DATA / "cluster_feature_means.tsv", sep="\t")
print("heatmap written | chrom vmax=%.2f | per cluster n=%s" % (chrom_vmax, dict(n_by.loc[CL])))
print("clustered %d of universe %d | features %d x clusters %d" % (n_clustered, n_universe, len(cols), len(CL)))
print(H.round(2).to_string())
