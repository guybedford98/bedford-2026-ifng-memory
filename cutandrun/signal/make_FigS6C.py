#!/usr/bin/env python3
"""
make_FigS6C.py  --  naive H3K27me3 distribution across mappable ISG promoters, with the
four manuscript anchors marked: memory genes (IDO1, CXCL11) in red, non-memory controls
(STAT1, IRF1) in grey.

Data: data/isg_promoter_focal_cpm.csv, column focal_Naive, mean CPM over TSS +/- 3 kb for
1,270 mappable ISG promoters (UCSC refGene coordinates).

The x axis is log10(focal CPM + 1e-4). The distribution is heavily right-skewed (median
0.0101, max 0.5705 CPM), so on a linear axis 95% of the genes squeeze into the left 39%
of the range and all four anchors land inside the leftmost 15%, with STAT1 and IRF1
indistinguishable. The log axis spreads the same four anchors over 1.34 log units.
What quantifies the separation is the percentile, so each anchor carries one, printed as
a compact key in the corner rather than left to the eye.

Run:  py make_FigS6C.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
import config as cfg  # noqa: E402
import figstyle as fs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

PANEL = "FigS6C"
CSV = ROOT / "data" / "isg_promoter_focal_cpm.csv"

PSEUDO = 1e-4
MEMORY = "#b2182b"     # IDO1, CXCL11
CONTROL = fs.NEUTRAL   # STAT1, IRF1
BAR = "#cccccc"
EDGE = fs.BAR_EDGE

# order left to right on the axis
ANCHORS = [("STAT1", CONTROL), ("IRF1", CONTROL), ("IDO1", MEMORY), ("CXCL11", MEMORY)]

df = pd.read_csv(CSV)
cpm = df["focal_Naive"].to_numpy(float)
n = len(cpm)

vals, pcts = {}, {}
for g, _ in ANCHORS:
    row = df.loc[df["gene"] == g, "focal_Naive"]
    if not len(row):
        raise SystemExit("anchor %s is not in the focal table" % g)
    vals[g] = float(row.iloc[0])
    pcts[g] = float((cpm < vals[g]).mean() * 100)

print("n = %d mappable ISG promoters" % n)
for g, _ in ANCHORS:
    print("  %-7s %.5f CPM   log10 %+.3f   percentile %.0f"
          % (g, vals[g], np.log10(vals[g] + PSEUDO), pcts[g]))
print("  IDO1/IRF1 %.1fx   CXCL11/STAT1 %.1fx"
      % (vals["IDO1"] / vals["IRF1"], vals["CXCL11"] / vals["STAT1"]))


def panel():
    x = np.log10(cpm + PSEUDO)
    pos = {g: np.log10(v + PSEUDO) for g, v in vals.items()}
    bins = np.linspace(x.min(), x.max(), 46)

    fig, ax = plt.subplots(figsize=(3.0, 2.15))
    ax.hist(x, bins=bins, color=BAR, edgecolor=EDGE, linewidth=0.4)

    ymax = ax.get_ylim()[1]
    ax.set_ylim(0, ymax * 1.34)          # headroom for the rotated gene labels
    # STAT1 and IRF1 sit 0.12 log units apart, closer than two 8 pt rotated labels are
    # wide, so their labels are nudged apart while the lines stay on the true values.
    span = max(pos.values()) - min(pos.values())
    nudge = {"STAT1": -0.055 * span, "IRF1": 0.055 * span}
    for g, col in ANCHORS:
        ax.axvline(pos[g], color=col, lw=1.0, zorder=3)
        ax.text(pos[g] + nudge.get(g, 0.0), ymax * 1.02, g,
                rotation=90, ha="center", va="bottom",
                fontsize=8, color=col, weight="bold")

    # Percentile is what quantifies the separation, so it is printed rather than left
    # to the eye. Stacked in the upper-left, which the distribution leaves empty.
    for i, (g, col) in enumerate(ANCHORS):
        ax.text(0.02, 0.97 - i * 0.075, "%-7s %.0fth" % (g, pcts[g]),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=6, color=col, family="monospace")

    ax.set_xlabel("log10 focal Naive CPM (+1e-4)", fontsize=8)
    ax.set_ylabel("number of ISGs", fontsize=8)
    ax.set_title("Naive H3K27me3 at mappable ISG promoters (n = %d)" % n, fontsize=8)
    fs.prism(ax, ticklen=2.4)
    fs.save(fig, cfg.out_dir(PANEL) / PANEL)
    print("  wrote %s" % PANEL)


fs.setup()
panel()
