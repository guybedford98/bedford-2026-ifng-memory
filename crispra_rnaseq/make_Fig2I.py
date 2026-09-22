#!/usr/bin/env python3
"""Fig 2I: CXCL11 TPM at day 10 and day 21, four CRISPRa conditions.

Bar = mean, dots = replicates, brackets carry the DESeq2 BH-FDR of the matching
guide contrast. Reads cxcl11_tpm.json and rnaseq_results.json.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))

import json
import numpy as np
import matplotlib.pyplot as plt

import config as cfg
import figstyle as fs
import plotstats as ps

RES = cfg.out_dir("crispra_rnaseq") / "rnaseq_results.json"
TPMJSON = cfg.out_dir("crispra_rnaseq") / "cxcl11_tpm.json"

GENE = "CXCL11"
COND = ["NT_rest", "NT_IFNg", "sg1_rest", "sg1_IFNg"]
COLOR = {"NT_rest": fs.COLOR["naive"], "NT_IFNg": fs.COLOR["crispra_nt"],
         "sg1_rest": fs.COLOR["rested"], "sg1_IFNg": fs.COLOR["crispra_tgt"]}
LAB1 = {"NT_rest": "NT rest", "NT_IFNg": "NT +IFNγ",
        "sg1_rest": "sg1 rest", "sg1_IFNg": "sg1 +IFNγ"}
DOT_JITTER = 0.15
DAYS = ["day10", "day21"]
DAYLAB = {"day10": "Day 10", "day21": "Day 21"}
# guide comparisons -> which contrast padj annotates the bracket
PAIRS = [("NT_rest", "sg1_rest", "priming_at_rest"),
         ("NT_IFNg", "sg1_IFNg", "memory_rechallenge")]

data = json.loads(RES.read_text())
tpmdata = json.loads(TPMJSON.read_text())


def counts(tp):
    """Per-condition replicate values for GENE, in TPM."""
    return {c: np.array(tpmdata[tp]["targets"][GENE][c]) for c in COND}


def padj(tp, key):
    return data["datasets"][tp]["contrasts"][key]["genes"][GENE]["padj"]


def draw_panel(ax, tp, show_ylab=True):
    """One timepoint: 4 bars (mean) + replicate dots + 2 sig brackets."""
    v = counts(tp)
    x = np.arange(4)
    fs.bars(ax, [v[c] for c in COND], x=x, width=0.66,
            colors=[COLOR[c] for c in COND], sem=False, jitter=DOT_JITTER)
    tops = {c: float(v[c].max()) for c in COND}
    vmax = max(tops.values())
    h = 0.02 * vmax
    placed = []
    xmap = {c: xi for xi, c in zip(x, COND)}
    for a, b, key in PAIRS:
        xa, xb = xmap[a], xmap[b]; lo, hi = min(xa, xb), max(xa, xb)
        base = max(tops[c] for c in COND if lo <= xmap[c] <= hi)
        y = base + 0.06 * vmax
        for plo, phi, py in placed:
            if not (hi < plo or lo > phi): y = max(y, py + 0.13 * vmax)
        stars = ps.pstars(ps.deseq_padj(padj(tp, key))["p"])
        fs.sig_bracket(ax, xa, xb, y, stars, h=h)
        placed.append((lo, hi, y))
    top = max(vmax * 1.10, max((p[2] for p in placed), default=0) + 0.16 * vmax)
    ax.set_ylim(0, top); ax.set_xlim(-0.62, 3.62)
    ax.set_xticks(x)
    ax.set_xticklabels([LAB1[c] for c in COND], rotation=45, ha="right")
    ax.set_title(DAYLAB[tp])
    if show_ylab:
        ax.set_ylabel(f"{GENE} TPM")
    fs.prism(ax)


def main():
    fs.setup()
    fig, axes = plt.subplots(1, 2, figsize=fs.PANEL_HALF)
    draw_panel(axes[0], "day10", show_ylab=True)
    draw_panel(axes[1], "day21", show_ylab=False)
    fig.subplots_adjust(left=0.185, right=0.985, top=0.86, bottom=0.30, wspace=0.42)
    fs.save(fig, cfg.out_dir("Fig2I") / "Fig2I")


if __name__ == "__main__":
    main()
