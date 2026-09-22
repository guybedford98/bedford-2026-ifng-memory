#!/usr/bin/env python3
"""Memory trajectory panels from the MSP summary table written by 00_msp_summary.py.

  Fig 4B : mean per-fiber MSP-covered area (%), all three loci.
  Fig 4C : fraction of MSP-rich fibers (%) for the two bistable loci IDO1 and CXCL11
           (canonical balanced k=2 model). STAT1 is the broadly-open foil (not bistable)
           and is shown in Fig 4B only.

Markers are colored by condition to match the Fig 4A rails.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import figstyle as fs

WIDTH_IN = 2.6
PANEL_H = 1.7

COND_ORDER = ["Naive", "Stimulated", "Rested", "Rechallenged"]
COND_SHORT = ["Naive", "Stim", "Rested", "Rech"]
COND_COLORS = {"Naive": "#909499", "Stimulated": "#4D9AD4", "Rested": "#4A4A4B", "Rechallenged": "#0B519D"}
LOCUS_MARKER = {"IDO1": "o", "CXCL11": "s", "STAT1": "^"}
LOCUS_LS = {"IDO1": "-", "CXCL11": "--", "STAT1": ":"}


def make_panel(df, loci, col, ylab, yticks, ylim, panel):
    x = np.arange(len(COND_ORDER))
    fig, ax = plt.subplots(figsize=(WIDTH_IN, PANEL_H))
    end_y = {}
    for locus in loci:
        sub = df[df["locus"] == locus].set_index("condition").reindex(COND_ORDER)
        y = sub[col].to_numpy()
        end_y[locus] = float(y[-1])
        ax.plot(x, y, color="#888888", lw=0.9, ls=LOCUS_LS[locus], zorder=1)
        for xi, cond in zip(x, COND_ORDER):
            ax.scatter(xi, sub.loc[cond, col], s=18, marker=LOCUS_MARKER[locus],
                       color=COND_COLORS[cond], edgecolor="#222222", linewidth=0.4, zorder=3)
    # de-collide the line-end locus labels (greedy minimum spacing)
    min_gap = (ylim[1] - ylim[0]) * 0.075
    label_y, prev = {}, -1e9
    for locus in sorted(loci, key=lambda L: end_y[L]):
        ly = max(end_y[locus], prev + min_gap)
        label_y[locus] = ly
        prev = ly
    for locus in loci:
        ax.text(x[-1] + 0.12, label_y[locus], f" {locus}", va="center", ha="left",
                fontsize=6, fontstyle="italic", color="#222222")

    ax.set_ylabel(ylab, fontsize=8)
    ax.set_ylim(*ylim)
    ax.set_yticks(yticks)
    ax.set_xticks(x)
    ax.set_xticklabels(COND_SHORT, rotation=30, ha="right", fontsize=6)
    ax.set_xlim(-0.3, len(COND_ORDER) - 1 + 1.05)
    ax.tick_params(axis="y", labelsize=6, length=2, pad=1)
    ax.tick_params(axis="x", length=2, pad=1)
    ax.grid(axis="y", color="#ebebeb", lw=0.4)
    ax.set_axisbelow(True)
    fs.prism(ax, lw=0.6, ticklen=2)
    fs.save(fig, cfg.out_dir(panel) / panel)
    print(f"{panel}: {WIDTH_IN}x{PANEL_H} in")


def main() -> None:
    fs.setup()
    tsv = cfg.out_dir("tables") / "msp_heatmap_summary.tsv"
    df = pd.read_csv(tsv, sep="\t")
    df["area_pct"] = pd.to_numeric(df["selected_mean_msp_covered_fiber_area_fraction"], errors="coerce") * 100
    df["rich_pct"] = pd.to_numeric(df["msp_rich_fraction_of_eligible"], errors="coerce") * 100   # NaN for STAT1 (foil, unclustered)

    make_panel(df, ["IDO1", "CXCL11", "STAT1"], "area_pct", "Mean MSP area (%)", [0, 20, 40], (0, 42), "Fig4B")
    make_panel(df, ["IDO1", "CXCL11"], "rich_pct", "MSP-rich fibers (%)", [0, 50, 100], (0, 105), "Fig4C")


if __name__ == "__main__":
    main()
