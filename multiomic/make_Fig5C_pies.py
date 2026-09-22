"""Fig 5C pies: one small composition pie per cluster, placed over the bars of Fig5C_bars.

Emitted as ONE FILE PER CLUSTER (not a composite) so each can be placed independently, with a
transparent background and no frame, axis or title. Diameter is read from Fig5C_bars.json
(`bar_width_mm`), so a pie drops onto its bar at 100% scale with no rescaling. Run
make_Fig5C_bars.py first; this script reads its geometry.

Within a cluster, memory + refractory + non-memory = 100%, and refractory is its own class,
never folded into non-memory.

A pie cannot show a sliver. The memory wedges are cluster 1 3/782 = 0.4% -> 1.4 degrees and
cluster 2 2/138 = 1.5% -> 5.3 degrees, at or below the width of the wedge outline itself at
9.7 mm. The wedge is drawn with a thin outline so it is at least detectable, and the exact
angles are written to the JSON, but a reader cannot judge those two proportions from the pie.

Output: results/Fig5C/Fig5C_pie_cluster{N}.{pdf,png,svg} (transparent, for placement)
        results/Fig5C/Fig5C_pie_cluster{N}_labelled.{pdf,png,svg}
        results/Fig5C/Fig5C_pies_contact_sheet.{pdf,png,svg}, Fig5C_pies.json
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
mm = fs.MM

MEM, REF, NON = "#b2182b", "#2a81ba", "#C7CBD1"
NON_TXT = "#6b7078"   # the non-memory fill is too pale to read as text
WEDGE_EDGE = "white"
WEDGE_LW = 0.5
FALLBACK_D_MM = 9.7

# diameter = the bar width of the Fig 5C bars, so the pie drops on at 100%
sizes_json = OUT / "Fig5C_bars.json"
if sizes_json.exists():
    D_MM = float(json.loads(sizes_json.read_text(encoding="utf-8"))["bar_width_mm"])
    print(f"pie diameter {D_MM:.2f} mm, taken from Fig5C_bars.json bar_width_mm")
else:
    D_MM = FALLBACK_D_MM
    print(f"Fig5C_bars.json not found; using fallback diameter {D_MM} mm "
          "(run make_Fig5C_bars.py first for an exact match)")

g = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t")
order = sorted(g["cluster"].unique(), key=lambda c: int(c.split("_")[1]))
disp = [int(fs.clabel(c)) for c in order]

rows = []
for c, d in zip(order, disp):
    sub = g[g.cluster == c]
    n_tot = int(len(sub))
    n_mem = int((sub.memory_class == "memory").sum())
    n_ref = int((sub.memory_class == "refractory").sum())
    n_non = n_tot - n_mem - n_ref
    rows.append(dict(cluster=d, n_total=n_tot, n_memory=n_mem, n_refractory=n_ref,
                     n_nonmemory=n_non,
                     pct_memory=100 * n_mem / n_tot, pct_refractory=100 * n_ref / n_tot,
                     pct_nonmemory=100 * n_non / n_tot))


def draw_pie(ax, r):
    """Wedges in a fixed order, starting at 12 o'clock, clockwise, so the four are comparable."""
    vals = [r["n_memory"], r["n_refractory"], r["n_nonmemory"]]
    cols = [MEM, REF, NON]
    keep = [(v, c) for v, c in zip(vals, cols) if v > 0]   # a 0-count wedge draws nothing
    ax.pie([v for v, _ in keep], colors=[c for _, c in keep],
           startangle=90, counterclock=False,
           wedgeprops=dict(edgecolor=WEDGE_EDGE, linewidth=WEDGE_LW, antialiased=True))
    ax.set_aspect("equal")
    # ax.pie draws at radius 1 but leaves a default margin, which shrinks the circle to 83% of
    # the canvas. Pin the limits to the circle so the drawn diameter IS the file size and it
    # drops on at 100%. The half-linewidth allowance keeps the wedge outline from being clipped.
    pad = (WEDGE_LW / 72) / 2 / (D_MM * mm) * 2
    ax.set_xlim(-1 - pad, 1 + pad); ax.set_ylim(-1 - pad, 1 + pad)


def save_transparent(fig, stem, **kw):
    """These pies are placed at 100% scale, so the canvas must not be cropped or padded by
    figstyle.save()'s tight bounding box."""
    for ext in ("pdf", "png", "svg"):
        fig.savefig(f"{stem}.{ext}", dpi=600, transparent=True, **kw)
    plt.close(fig)


for r in rows:
    fig, ax = plt.subplots(figsize=fs.panel_mm(D_MM, D_MM))
    draw_pie(ax, r)
    ax.set_axis_off()
    fig.subplots_adjust(left=0, right=1, top=1, bottom=0)   # pie fills the canvas exactly
    save_transparent(fig, OUT / f"Fig5C_pie_cluster{r['cluster']}",
                     pad_inches=0, bbox_inches=None)
    deg_mem = 360 * r["pct_memory"] / 100
    flag = "  <- sliver, not readable as a proportion" if deg_mem < 10 else ""
    print(f"  Fig5C_pie_cluster{r['cluster']}: memory {r['n_memory']}/{r['n_total']} "
          f"({r['pct_memory']:.1f}%, {deg_mem:.1f} deg){flag}")

# The placement pies above are 9.7 mm and cannot legibly carry text. These are the same pies at
# a size that can, with ALL THREE percentages spelled out, colour-coded, under the pie. Wedge
# area is unreliable for the small classes (cluster 1's memory wedge is 1.4 deg), so the numbers
# are the thing to read; the pie is the shape.
LBL_D_MM = 26.0
for r in rows:
    fig = plt.figure(figsize=fs.panel_mm(LBL_D_MM, LBL_D_MM + 13))
    ax = fig.add_axes([0.0, 13 / (LBL_D_MM + 13), 1.0, LBL_D_MM / (LBL_D_MM + 13)])
    draw_pie(ax, r)
    ax.set_axis_off()
    for i, (lab, key, col) in enumerate([("memory", "pct_memory", MEM),
                                         ("refractory", "pct_refractory", REF),
                                         ("non-memory", "pct_nonmemory", NON_TXT)]):
        fig.text(0.5, 0.24 - i * 0.075, f"{lab} {r[key]:.1f}%", ha="center", va="top",
                 fontsize=8, color=col, fontweight="bold")
    save_transparent(fig, OUT / f"Fig5C_pie_cluster{r['cluster']}_labelled",
                     bbox_inches="tight")
print(f"  labelled variants written at {LBL_D_MM:.0f} mm with all three percentages")

fig, axs = plt.subplots(1, len(rows), figsize=fs.panel_mm(len(rows) * 24, 40))
for ax, r in zip(np.atleast_1d(axs), rows):
    draw_pie(ax, r)
    ax.set_axis_off()
    # all three percentages, not just memory: with refractory kept as its own class the pie has
    # three wedges, and quoting only "% memory" is what makes the numbers look inconsistent
    # against panels that exclude refractory. Stacked, not one line: on one line the three
    # values are wider than the tile and run into the neighbouring panel's title.
    ax.set_title(f"cluster {r['cluster']}  (n={r['n_total']})\n"
                 f"mem {r['pct_memory']:.1f}%\n"
                 f"ref {r['pct_refractory']:.1f}%\n"
                 f"non {r['pct_nonmemory']:.1f}%",
                 fontsize=8, linespacing=1.3)
fig.subplots_adjust(left=0.01, right=0.99, top=0.60, bottom=0.02, wspace=0.30)
fs.save(fig, OUT / "Fig5C_pies_contact_sheet")

(OUT / "Fig5C_pies.json").write_text(json.dumps({
    "panel": "Fig5C_pies",
    "purpose": ("per-cluster composition pies for placement over the bars of Fig5C_bars; "
                "one file per cluster, transparent, no frame"),
    "diameter_mm": round(D_MM, 3),
    "diameter_source": "Fig5C_bars.json bar_width_mm" if sizes_json.exists() else "fallback",
    "wedge_order": ["memory", "refractory", "non-memory"],
    "start_angle_deg": 90, "direction": "clockwise",
    "colours": {"memory": MEM, "refractory": REF, "non_memory": NON},
    "n_denominates": ("each pie is that CLUSTER's composition, so the three wedges sum to 100% "
                      "and the denominator is the cluster, not the universe. DIFFERENT from the "
                      "Fig5C bars, whose heights are shares of the whole clustered set."),
    "why_this_differs_from_the_rna_violin_panel": (
        "Fig 5E EXCLUDES refractory genes by design, so its '% memory' uses "
        "memory/(memory+non-memory) as the denominator. Cluster 3 is 39.0% here (32/82, "
        "refractory included) but 50.8% there (32/63, refractory excluded). Both are correct; "
        "they answer different questions. Quote the one whose denominator matches the panel."),
    "clusters": [{**r,
                  "memory_wedge_deg": round(360 * r["pct_memory"] / 100, 2),
                  "refractory_wedge_deg": round(360 * r["pct_refractory"] / 100, 2),
                  "readable_memory_wedge": (360 * r["pct_memory"] / 100) >= 10}
                 for r in rows],
    "caveat": ("A pie cannot render a sliver. Memory wedges below ~10 degrees (clusters with "
               "under ~3% memory) are detectable but not measurable by eye."),
}, indent=2), encoding="utf-8")

print(f"\nwrote {len(rows)} pies + contact sheet + Fig5C_pies.json to {OUT}")
