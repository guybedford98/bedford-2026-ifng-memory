"""Fig 5E: per-cluster RNA response, memory vs non-memory violins with one dot per gene.

log2FC vs naive for Stimulated, Rested and Rechallenged; Naive is the reference and is not
plotted. Two violins per condition, memory and non-memory, each with its genes overlaid.

Refractory genes are excluded entirely (not folded into non-memory), so the n and the
"% memory" in each title are over memory + non-memory only. That percentage therefore differs
from the Fig 5C composition pies', which keep refractory in the denominator.

Contrasts come from the same feature matrix the clustering used, with Rechallenged
reconstructed as `RNA_Stim_vs_Naive + RNA_Rechallenged_vs_Stim` (the matrix stores rechallenge
relative to stim).

Output: results/Fig5E/Fig5E.{pdf,png,svg}, data/rna_violin_medians.tsv
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

fs.setup()
DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("Fig5E")

CLUSTER_COLORS = {1: "#20566b", 2: "#c9a227", 3: "#2e8b57", 4: "#b0532b"}
FILL = {"memory": "#b2182b", "non_memory": "#C7CBD1"}
DOT = {"memory": "#8f1322", "non_memory": "#8A9099"}   # darker, so points read against the fill
CLAB = {"memory": "memory", "non_memory": "non-memory"}
ORDER = ["memory", "non_memory"]
OFF = {"memory": -0.19, "non_memory": 0.19}
CONDS = ["Stim", "Rested", "Rechallenged"]
SHORT = ["Stim", "Rested", "Rech"]
RNG = np.random.RandomState(42)

cl = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t").set_index("gene")
fm = pd.read_csv(DATA / "feature_matrix.tsv", sep="\t").set_index("gene")
df = fm.join(cl[["cluster", "memory_class"]], how="inner")
df["Stim"] = df["RNA_Stim_vs_Naive"]
df["Rested"] = df["RNA_Rested_vs_Naive"]
df["Rechallenged"] = df["RNA_Stim_vs_Naive"] + df["RNA_Rechallenged_vs_Stim"]
print(f"clustered genes joined to features: {len(df)}")
n_refr = int((df.memory_class == "refractory").sum())
df = df[df.memory_class.isin(ORDER)].copy()
CL = sorted(df["cluster"].unique(), key=lambda x: int(x.split("_")[1]))
print(f"genes plotted: {len(df)} (refractory excluded: {n_refr})")

ylo = float(df[CONDS].min().min()); yhi = float(df[CONDS].max().max())
pad = 0.06 * (yhi - ylo); YL = (ylo - pad, yhi + pad)

# Marker sizes are in points^2 so the dots keep their physical size; only the panel is shorter.
PANEL_W_MM, PANEL_H_MM = 36, 46 * 0.75
fig, axs = plt.subplots(1, len(CL),
                        figsize=fs.panel_mm(len(CL) * PANEL_W_MM, PANEL_H_MM), squeeze=False)
for j, c in enumerate(CL):
    ax = axs[0][j]; sub = df[df.cluster == c]
    ax.axhline(0, color="#CFCFCF", lw=0.6, zorder=1)
    for k, cond in enumerate(CONDS):
        for klass in ORDER:
            v = sub.loc[sub.memory_class == klass, cond].values
            if not len(v): continue
            pos = k + OFF[klass]
            if len(v) >= 5:
                p = ax.violinplot([v], positions=[pos], widths=0.32, showextrema=False,
                                  showmedians=True)
                for b in p["bodies"]:
                    b.set_facecolor(FILL[klass]); b.set_edgecolor(FILL[klass])
                    b.set_linewidth(0.7); b.set_alpha(0.32); b.set_zorder(2)
                p["cmedians"].set_color("#2b2b2b")
                p["cmedians"].set_linewidth(0.8); p["cmedians"].set_zorder(5)
            big = len(v) > 200
            ax.scatter(np.full(len(v), pos) + RNG.uniform(-0.085, 0.085, len(v)), v,
                       s=2.2 if big else 4.0, color=DOT[klass], linewidths=0.0,
                       alpha=0.45 if big else 0.85, zorder=3 if big else 4)
    ax.set_xticks(np.arange(len(CONDS))); ax.set_xticklabels(SHORT, fontsize=6)
    ax.set_xlim(-0.55, len(CONDS) - 0.45); ax.set_ylim(*YL)
    fs.prism(ax, lw=1.0); ax.tick_params(labelsize=6); ax.tick_params(axis="x", length=0)
    pct_mem = 100.0 * (sub.memory_class == "memory").sum() / len(sub)
    ax.set_title(f"cluster {fs.clabel(c)}\n(n={len(sub)}, {pct_mem:.1f}% memory)", fontsize=8,
                 fontweight="bold", color=CLUSTER_COLORS[int(fs.clabel(c))])
    if j == 0: ax.set_ylabel("RNA log2FC\nvs naive", fontsize=8, fontweight="bold")
    else: ax.tick_params(labelleft=False)

handles = [Line2D([0], [0], marker="o", linestyle="none", markersize=3.4,
                  markerfacecolor=FILL[k], markeredgecolor="none", label=CLAB[k]) for k in ORDER]
# Margins are held in MILLIMETRES, not as fractions of the figure. The two-line cluster titles
# and the legend need a fixed amount of room; expressing them fractionally means any change to
# the panel height reflows them, which is what collided the legend with the titles.
TITLE_MM, BOTTOM_MM, LEG_GAP_MM = 10.1, 5.5, 4.6
fig.legend(handles=handles, loc="upper center",
           bbox_to_anchor=(0.5, 1 + LEG_GAP_MM / PANEL_H_MM), ncol=2, frameon=False,
           fontsize=8, handlelength=1.0, columnspacing=1.4)
fig.subplots_adjust(left=0.12, right=0.995, wspace=0.16,
                    top=1 - TITLE_MM / PANEL_H_MM, bottom=BOTTOM_MM / PANEL_H_MM)
fs.save(fig, OUT / "Fig5E")

rows = []
for c in CL:
    sub = df[df.cluster == c]
    for klass in ORDER:
        g = sub[sub.memory_class == klass]
        if not len(g): continue
        rows.append({"cluster": fs.clabel(c), "memory_class": klass, "n": len(g),
                     **{f"median_{k}": round(float(g[k].median()), 2) for k in CONDS}})
med = pd.DataFrame(rows)
med.to_csv(DATA / "rna_violin_medians.tsv", sep="\t", index=False)
print(med.to_string(index=False))
print("Fig5E written")
