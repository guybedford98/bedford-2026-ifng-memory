"""Fig S4G (IDO1, CXCL11) and Fig S5C (CXCR3): TPM across the four IFNg conditions.

Bar = mean TPM, overlaid individual replicate points (n=3).
Condition -> sample mapping:
    Naive        = NTC          (S1-S3)
    Stimulated   = Naive+IFNg   (S4-S6)
    Rested       = Trained      (S7-S9)
    Rechallenged = Trained+IFNg (S10-S12)
"""
import os
import csv
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

for fp in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\arialbd.ttf"):
    if os.path.exists(fp):
        font_manager.fontManager.addfont(fp)
plt.rcParams["font.family"] = "Arial"
plt.rcParams["svg.fonttype"] = "none"   # keep text editable in Illustrator
plt.rcParams["pdf.fonttype"] = 42       # embed TrueType (editable) in PDF
plt.rcParams["ps.fonttype"] = 42

CONDITIONS = ["Naive", "Stimulated", "Rested", "Rechallenged"]
COLORS = ["#DDDDDD", "#4D9BD0", "#7F7F7F", "#1C5C9F"]  # lt gray, med blue, dk gray, dk blue

# Which panel each gene belongs to, in draft0012 numbering.
PANEL = {
    "IDO1":   ("FigS4G", "FigS4G_IDO1"),
    "CXCL11": ("FigS4G", "FigS4G_CXCL11"),
    "CXCR3":  ("FigS5C", "FigS5C"),
}

YAXIS = {
    "IDO1":   (2000, [0, 500, 1000, 1500, 2000]),
    "CXCL11": (2000, [0, 500, 1000, 1500, 2000]),
    "CXCR3":  (2000, [0, 500, 1000, 1500, 2000]),  # flat zero vs ligand axis
}

DATA = {
    "IDO1": {
        "Naive":        [0.0,    1.58,   1.77],
        "Stimulated":   [512.55, 491.75, 512.80],
        "Rested":       [4.49,   9.78,   11.94],
        "Rechallenged": [1479.19, 1619.31, 1513.05],
    },
    "CXCL11": {
        "Naive":        [0.0,    1.59,   2.00],
        "Stimulated":   [593.10, 641.29, 637.17],
        "Rested":       [5.99,   10.85,  17.46],
        "Rechallenged": [1413.72, 1676.76, 1411.91],
    },
    # CXCR3 receptor: not expressed in MSCs (0 TPM in every sample); shown on the
    # same 0-2000 axis as its ligand CXCL11 to convey receptor absence vs ligand induction.
    "CXCR3": {
        "Naive":        [0.0, 0.0, 0.0],
        "Stimulated":   [0.0, 0.0, 0.0],
        "Rested":       [0.0, 0.0, 0.0],
        "Rechallenged": [0.0, 0.0, 0.0],
    },
}

# horizontal offsets for the 3 replicate dots within a bar
DOT_OFFSETS = [-0.15, 0.0, 0.15]


def make_panel(gene):
    vals = DATA[gene]
    x = np.arange(len(CONDITIONS))
    means = [float(np.mean(vals[c])) for c in CONDITIONS]

    fig, ax = plt.subplots(figsize=(1.75, 2.0))

    ax.bar(x, means, width=0.65, color=COLORS, edgecolor="black",
           linewidth=1.0, zorder=2)

    for i, c in enumerate(CONDITIONS):
        ys = vals[c]
        xs = [x[i] + o for o in DOT_OFFSETS]
        ax.plot(xs, ys, linestyle="None", marker="o", markersize=4.5,
                markerfacecolor=COLORS[i], markeredgecolor="black",
                markeredgewidth=0.8, zorder=3, clip_on=False)

    ymax, yticks = YAXIS[gene]
    ax.set_ylim(0, ymax)
    ax.set_yticks(yticks)
    ax.set_ylabel(f"{gene} TPM", fontsize=8, fontweight="bold")
    ax.tick_params(axis="y", labelsize=6, length=3.5, width=1.2, direction="out")

    ax.set_xticks(x)
    ax.set_xticklabels(CONDITIONS, fontsize=8, rotation=45, ha="right")
    ax.tick_params(axis="x", length=3.5, width=1.2, direction="out")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.set_xlim(-0.6, len(CONDITIONS) - 0.4)

    fig.tight_layout(pad=0.4)
    panel_dir, stem = PANEL[gene]
    fs.save(fig, cfg.out_dir(panel_dir) / stem)
    print(f"  {gene}: means = " + ", ".join(f"{c} {m:.1f}" for c, m in zip(CONDITIONS, means)))


def write_csv():
    out = cfg.out_dir("rnaseq_tables") / "IDO1_CXCL11_CXCR3_TPM.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["gene", "condition", "replicate", "sample", "TPM"])
        sample_map = {
            "Naive": ["S1 NTC", "S2 NTC", "S3 NTC"],
            "Stimulated": ["S4 IFNg", "S5 IFNg", "S6 IFNg"],
            "Rested": ["S7 Trained", "S8 Trained", "S9 Trained"],
            "Rechallenged": ["S10 Rechallenge", "S11 Rechallenge", "S12 Rechallenge"],
        }
        for gene, conds in DATA.items():
            for c in CONDITIONS:
                for r, (s, v) in enumerate(zip(sample_map[c], conds[c]), start=1):
                    w.writerow([gene, c, r, s, v])
    print(f"  wrote {out}")


if __name__ == "__main__":
    print("Rendering supplementary TPM panels:")
    for gene in ("IDO1", "CXCL11", "CXCR3"):
        make_panel(gene)
    write_csv()
