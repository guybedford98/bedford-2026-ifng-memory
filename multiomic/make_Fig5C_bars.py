"""Fig 5C bars: how the clustered ISGs split across the four GMM clusters, so the small size
of the memory-enriched cluster is explicit before any memory overlay. All three classes
(memory / refractory / non-memory) are included; nothing is excluded here.

The per-cluster composition pies go on top of these bars (make_Fig5C_pies.py), which reads the
bar geometry this script writes to Fig5C_bars.json.

Output: results/Fig5C/Fig5C_bars.{pdf,png,svg,json}, data/cluster_sizes.tsv
"""
import sys, json
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

fs.setup()
DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("Fig5C")

CLUSTER_COLORS = {1: "#20566b", 2: "#c9a227", 3: "#2e8b57", 4: "#b0532b"}

g = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t")
u = pd.read_csv(DATA / "universe.tsv", sep="\t")
n_universe = int(len(u))
order = sorted(g["cluster"].unique(), key=lambda c: int(c.split("_")[1]))
disp = [int(fs.clabel(c)) for c in order]
n = np.array([int((g.cluster == c).sum()) for c in order])
total = int(n.sum())
pct = 100.0 * n / total
# The annotation prints BOTH n and %, and the JSON keeps both series, so nothing is lost
# either way.
Y_MODE = "count"                                  # "count" | "percent"
height = n.astype(float) if Y_MODE == "count" else pct
YLABEL = "number of ISGs" if Y_MODE == "count" else "% of clustered ISGs"
print("n per cluster:", n.tolist(), "| total:", total, "| pct:", np.round(pct, 1).tolist())
print("y axis:", Y_MODE)
print("universe:", n_universe, "| unclustered:", sorted(set(u.gene) - set(g.gene)))

mm = fs.MM
# 72 mm wide so the per-cluster composition pies can sit over the bars at the bars' own width.
FIG_W_MM = 72

# The AXES height is matched to the PC scatter panels so the two sit level in a figure.
# Measured from a throwaway figure at the scatters' own figsize rather than hard-coded, so it
# cannot drift if those panels are resized again. Source of the 70 x 60 and the sqrt(0.8):
# make_pc_scatters.py.
SCATTER_W_MM, SCATTER_H_MM = 70 * 0.8 ** 0.5, 60 * 0.8 ** 0.5
_f, _a = plt.subplots(figsize=fs.panel_mm(SCATTER_W_MM, SCATTER_H_MM))
_f.canvas.draw()
AXES_H_MM = _a.get_window_extent().height / _f.dpi * 25.4
plt.close(_f)

# Margins in millimetres, so the axes height is exactly AXES_H_MM whatever the annotation needs.
# Only top/bottom are set: left/right keep matplotlib's defaults so the axes WIDTH -- and hence
# the bar width the pies are cut to -- is unchanged.
TOP_MM, BOTTOM_MM = 8.5, 9.5
FIG_H_MM = AXES_H_MM + TOP_MM + BOTTOM_MM
print(f"matching scatter axes height: {AXES_H_MM:.2f} mm -> figure {FIG_W_MM} x {FIG_H_MM:.2f} mm")
fig, ax = plt.subplots(figsize=fs.panel_mm(FIG_W_MM, FIG_H_MM))
fig.subplots_adjust(top=1 - TOP_MM / FIG_H_MM, bottom=BOTTOM_MM / FIG_H_MM)
x = np.arange(len(order))
ax.bar(x, height, width=0.74, color=[CLUSTER_COLORS[d] for d in disp],
       edgecolor=fs.BAR_EDGE, linewidth=fs.BAR_EDGE_LW, zorder=3)
_gap = max(height) * 0.015
for xi, h, p, ni in zip(x, height, pct, n):
    ax.text(xi, h + _gap, f"n={ni}\n{p:.1f}%", ha="center", va="bottom", fontsize=8,
            color="#333", linespacing=1.15)

fs.prism(ax, lw=1.0)
ax.set_xticks(x); ax.set_xticklabels([str(d) for d in disp])
ax.set_xlabel("cluster", fontsize=8, fontweight="bold")
ax.set_ylabel(YLABEL, fontsize=8, fontweight="bold")
ax.set_ylim(0, max(height) * 1.22)
ax.tick_params(labelsize=6); ax.tick_params(axis="x", length=0)
ax.set_xlim(-0.62, len(order) - 0.38)
ax.set_title(f"{total:,} of {n_universe:,} ISGs\n(max TPM ≥ 5) with multi-omic features",
             fontsize=8, color="#555", pad=4, linespacing=1.2)

# Record the bar geometry in millimetres so the composition pies can be emitted at exactly the
# bar width and dropped straight onto this panel with no rescaling.
fig.canvas.draw()
_bb = ax.get_window_extent()
_ax_w_mm = _bb.width / fig.dpi * 25.4
_xspan = ax.get_xlim()[1] - ax.get_xlim()[0]
_bar_w_mm = _ax_w_mm * (0.74 / _xspan)              # 0.74 = the bar width in data units
_bar_pitch_mm = _ax_w_mm * (1.0 / _xspan)           # centre-to-centre spacing
_bar_x_mm = [(_bb.x0 + (xi - ax.get_xlim()[0]) / _xspan * _bb.width) / fig.dpi * 25.4
             for xi in x]
_bar_top_mm = [(_bb.y0 + (h_ / (ax.get_ylim()[1] - ax.get_ylim()[0])) * _bb.height)
               / fig.dpi * 25.4 for h_ in height]

fs.save(fig, OUT / "Fig5C_bars")

with open(OUT / "Fig5C_bars.json", "w", encoding="utf-8") as _fh:
    json.dump({"panel": "Fig5C_bars",
               "figure_mm": [FIG_W_MM, FIG_H_MM],
               "axes_width_mm": round(_ax_w_mm, 3),
               "bar_width_mm": round(_bar_w_mm, 3),
               "bar_pitch_mm": round(_bar_pitch_mm, 3),
               "bar_centre_x_mm_from_page_left": [round(v, 2) for v in _bar_x_mm],
               "bar_top_y_mm_from_axes_bottom": [round(v, 2) for v in _bar_top_mm],
               "clusters_display": disp,
               "n_genes": n.tolist(),
               "pct_of_clustered": [round(float(v), 2) for v in pct],
               "note": ("bar_width_mm is what make_Fig5C_pies.py uses as the pie diameter, so a "
                        "pie drops onto a bar at 100% scale. x positions are measured on the "
                        "UNCROPPED canvas; bbox_inches='tight' crops the saved file, so use "
                        "them as relative spacing, not absolute page coords.")},
              _fh, indent=2)
print(f"bar width {_bar_w_mm:.2f} mm, pitch {_bar_pitch_mm:.2f} mm -> pies sized to match")

pd.DataFrame({"cluster": disp, "n_genes": n, "pct_of_clustered": np.round(pct, 2)}).to_csv(
    DATA / "cluster_sizes.tsv", sep="\t", index=False)
print("cluster size distribution written")
