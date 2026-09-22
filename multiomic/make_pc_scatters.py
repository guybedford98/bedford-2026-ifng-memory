"""PCA scatters on the ISG universe, four panels meant to sit together:

  Fig5B     PC1 x PC2, coloured by memory class
  Fig5D     PC1 x PC2, coloured by GMM cluster
  FigS9C    PC1 x PC3, coloured by memory class
  FigS10B   PC1 x PC3, coloured by GMM cluster

One marker size everywhere (s=10, white hairline edge), with the small classes kept legible
through draw order rather than size; data-derived axis limits; only IDO1 and CXCL11 ringed and
labelled, anchor rings s=70.

Every dot is one clustered gene; no gene is excluded from any panel.
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
np.random.seed(42)

MEM = "#b2182b"      # memory (red)
REF = "#2a81ba"      # refractory (blue)
NON = "#C4C8CE"      # non-memory (grey)
CLUSTER_COLORS = {1: "#20566b", 2: "#c9a227", 3: "#2e8b57", 4: "#b0532b"}
ANCHORS = ["IDO1", "CXCL11"]      # only the two memory flagships are called out
DOT_S = 10
DOT_EW = 0.25
LEG_MS = 4.4
ANCHOR_S = 70

df = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t")
ve = pd.read_csv(DATA / "variance_explained.tsv", sep="\t", index_col=0)
pct = {pc: round(100 * float(ve.loc[pc, "variance_explained"])) for pc in ("PC1", "PC2", "PC3")}
df["disp"] = df["cluster"].map(fs.clabel).astype(int)
anchors = df[df["is_anchor"].isin(ANCHORS)].copy()
print("clustered genes plotted:", len(df))
print({k: int((df.memory_class == k).sum()) for k in ("memory", "refractory", "non_memory")})
print("cluster sizes:", df.disp.value_counts().sort_index().to_dict())
print("variance explained (%):", {k: round(100 * float(ve.loc[k, "variance_explained"]), 2)
                                  for k in ("PC1", "PC2", "PC3")})
print("anchors:", anchors[["gene", "disp", "PC1", "PC2", "PC3"]].to_string(index=False))

# Area scales as the square of the linear size, so the linear factor is sqrt(0.8): 70 x 60 mm
# -> 62.6 x 53.7 mm. Marker sizes are in points^2 and so are absolute -- they do NOT scale with
# the figure: the dots keep their physical size and occupy a larger fraction of the panel.
AREA_SCALE = 0.8
LIN = AREA_SCALE ** 0.5
PANEL_W_MM, PANEL_H_MM = 70 * LIN, 60 * LIN


def draw(ycol, colour_by, panel, legend_loc="upper right"):
    fig, ax = plt.subplots(figsize=fs.panel_mm(PANEL_W_MM, PANEL_H_MM))

    if colour_by == "memoryclass":
        non = df[df.memory_class == "non_memory"]
        ref = df[df.memory_class == "refractory"]
        mem = df[df.memory_class == "memory"]
        # same size for all three; the small classes stay legible through draw order, not size
        ax.scatter(non["PC1"], non[ycol], s=DOT_S, c=NON, edgecolors="white", linewidths=DOT_EW,
                   alpha=0.55, zorder=1, rasterized=True)
        ax.scatter(ref["PC1"], ref[ycol], s=DOT_S, c=REF, edgecolors="white", linewidths=DOT_EW, zorder=2)
        ax.scatter(mem["PC1"], mem[ycol], s=DOT_S, c=MEM, edgecolors="white", linewidths=DOT_EW, zorder=3)
        handles = [
            Line2D([0], [0], marker="o", ls="", mfc=MEM, mec="white", mew=DOT_EW, ms=LEG_MS, label=f"memory (n={len(mem)})"),
            Line2D([0], [0], marker="o", ls="", mfc=REF, mec="white", mew=DOT_EW, ms=LEG_MS, label=f"refractory (n={len(ref)})"),
            Line2D([0], [0], marker="o", ls="", mfc=NON, mec="white", mew=DOT_EW, ms=LEG_MS, label=f"non-memory (n={len(non)})"),
        ]
    else:
        handles = []
        for d in sorted(df.disp.unique()):
            g = df[df.disp == d]
            ax.scatter(g["PC1"], g[ycol], s=DOT_S, c=CLUSTER_COLORS[d], edgecolors="white",
                       linewidths=DOT_EW, alpha=0.9, zorder=2 + d, rasterized=(len(g) > 400))
            handles.append(Line2D([0], [0], marker="o", ls="", mfc=CLUSTER_COLORS[d], mec="white",
                                  mew=DOT_EW, ms=LEG_MS, label=f"cluster {d} (n={len(g)})"))

    ax.scatter(anchors["PC1"], anchors[ycol], s=ANCHOR_S, facecolors="none", edgecolors="black",
               linewidths=1.0, zorder=8)

    fs.prism(ax, lw=1.0)
    xpad = 0.06 * (df.PC1.max() - df.PC1.min()); ypad = 0.06 * (df[ycol].max() - df[ycol].min())
    ax.set_xlim(df.PC1.min() - xpad, df.PC1.max() + xpad)
    ax.set_ylim(df[ycol].min() - ypad, df[ycol].max() + ypad)
    ax.set_xlabel(f"PC1 ({pct['PC1']}% var)", fontsize=8, fontweight="bold")
    ax.set_ylabel(f"{ycol} ({pct[ycol]}% var)", fontsize=8, fontweight="bold")
    ax.tick_params(labelsize=6)
    ax.axhline(0, color="#DDDDDD", lw=0.6, zorder=0); ax.axvline(0, color="#DDDDDD", lw=0.6, zorder=0)

    # label offsets as a fraction of the axis range, so they behave on both PC2 and PC3
    xr = ax.get_xlim()[1] - ax.get_xlim()[0]; yr = ax.get_ylim()[1] - ax.get_ylim()[0]
    off = {"IDO1": (-0.016 * xr, 0.075 * yr, "right"), "CXCL11": (0.022 * xr, -0.01 * yr, "left")}
    for _, r in anchors.iterrows():
        g = r["is_anchor"]; dx, dy, ha = off[g]
        ax.annotate(g, xy=(r["PC1"], r[ycol]), xytext=(r["PC1"] + dx, r[ycol] + dy),
                    ha=ha, va="center", fontsize=8, fontweight="bold", zorder=9,
                    bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))

    ax.legend(handles=handles, loc=legend_loc, frameon=False, fontsize=8,
              handletextpad=0.4, labelspacing=0.35, borderaxespad=0.3)

    fs.save(fig, cfg.out_dir(panel) / panel)
    print("wrote", panel)


# The legend sits upper right on three of the four panels. On PC1 x PC3 coloured by cluster it
# cannot: IDO1 sits high on PC3 and the cluster legend is the only one with four rows, so its
# bottom row runs into IDO1's callout and the callout's white bbox half-erases the word
# "cluster". The label cannot move instead -- below its ring it would cover CXCL10 and GBP5.
# Lower right is empty on that panel (0 genes at PC1 >= 8 and PC3 <= -4).
draw("PC3", "memoryclass", "FigS9C")
draw("PC3", "clusters", "FigS10B", legend_loc="lower right")
draw("PC2", "memoryclass", "Fig5B")
draw("PC2", "clusters", "Fig5D")
