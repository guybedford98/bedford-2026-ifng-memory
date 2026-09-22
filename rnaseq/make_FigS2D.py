"""Fig S2D: naive-vs-rested volcano for the IFNg training set, with its TNFa companion.

Both panels are drawn here from one set of limits. The manuscript claims the TNFa drift
is smaller than the IFNg drift, and that comparison is only visual if the two panels
share axes and style. The TNFa panel and the side-by-side pair are not manuscript
panels; they carry the TNFa drift counts quoted in the text.

DEG cutoff: FDR < 0.05 and |log2FC| > 0.5. Both tables come from the four-condition
DESeq2 model (stages 05 and 06). The two-condition run in stage 04 filters differently
and gives a different IFNg count at the same cutoff, which is the count quoted in the
Results text.

Shared axes: the x range covers every significant gene in BOTH panels plus padding, so
no DEG is ever hidden. A handful of non-significant outliers can fall outside; the count
is printed for each panel and must be reported in the legend if it is not zero.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass

TABLES = cfg.out_dir("rnaseq_tables")
OUT = cfg.out_dir("FigS2D")

G = "\u03b3"  # gamma
A = "\u03b1"  # alpha
PADJ, CUT = 0.05, 0.5
C_NS, C_SIG = "#cccccc", "#c0392b"

PANELS = [
    ("IFNg", TABLES / "DEG_Naive_vs_Rested_full_model.csv",
     "IFN%s training: Rested vs Naive" % G, "FigS2D"),
    ("TNFa", TABLES / "DEG_TNFa_Naive_vs_Rested.csv",
     "TNF%s training: Rested vs Naive" % A, "TNFa_NaiveVsRested_volcano"),
]

for fp in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\arialbd.ttf"):
    if Path(fp).exists():
        fm.fontManager.addfont(fp)
FONT = "Arial" if "Arial" in {f.name for f in fm.fontManager.ttflist} else "DejaVu Sans"
plt.rcParams.update({
    "font.family": FONT, "font.size": 8, "axes.titlesize": 9,
    "axes.labelsize": 8, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "figure.dpi": 150, "savefig.dpi": 600,
    "svg.fonttype": "none", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.8})


def load(path):
    d = pd.read_csv(path)
    d.columns = [c.strip() for c in d.columns]
    x = d["log2FoldChange"].to_numpy(float)
    p = d["padj"].to_numpy(float)
    m = np.isfinite(p) & np.isfinite(x)      # a volcano needs both coordinates
    x, p = x[m], p[m]
    y = -np.log10(np.clip(p, 1e-300, 1.0))
    sig = (p < PADJ) & (np.abs(x) > CUT)
    return x, y, sig


sets = {tag: load(path) for tag, path, _, _ in PANELS}

# Shared limits: every significant gene in either panel stays in view.
sig_lo = min(x[s].min() for x, _, s in sets.values())
sig_hi = max(x[s].max() for x, _, s in sets.values())
XLIM = (float(np.floor(sig_lo - 0.5)), float(np.ceil(sig_hi + 0.5)))
YLIM = (0.0, max(y.max() for _, y, _ in sets.values()) * 1.08 + 1)
print("shared xlim %s   shared ylim (0, %.0f)" % (XLIM, YLIM[1]))


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(length=3, width=0.8)


def draw(ax, tag, title, shared=True):
    x, y, sig = sets[tag]
    ax.scatter(x[~sig], y[~sig], s=3.5, c=C_NS, linewidths=0, rasterized=True)
    ax.scatter(x[sig], y[sig], s=5, c=C_SIG, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(PADJ), color="#888", lw=0.6, ls="--")
    ax.axvline(CUT, color="#888", lw=0.6, ls="--")
    ax.axvline(-CUT, color="#888", lw=0.6, ls="--")

    if shared:
        ax.set_xlim(*XLIM)
        ax.set_ylim(*YLIM)
    else:
        ax.set_xlim(float(np.floor(x[sig].min() - 0.5)), float(np.ceil(x[sig].max() + 0.5)))
        ax.set_ylim(0, float(y.max()) * 1.08 + 1)

    lo, hi = ax.get_xlim()
    off = int(((x < lo) | (x > hi)).sum())
    off_sig = int((sig & ((x < lo) | (x > hi))).sum())

    n, up = int(sig.sum()), int((sig & (x > 0)).sum())
    ax.text(0.97, 0.97, "%d DEG\n(%d up, %d down)" % (n, up, n - up),
            transform=ax.transAxes, ha="right", va="top", fontsize=5.5, color="#555")
    ax.set_xlabel("log2 fold-change", fontsize=7)
    ax.set_title(title, fontsize=8, weight="bold", pad=4)
    style(ax)
    print("  %-5s DEG %4d (up %3d, down %3d)   points outside the axes: %d (significant: %d)"
          % (tag, n, up, n - up, off, off_sig))
    return n


def save(fig, stem):
    for ext in ("pdf", "svg", "png"):
        fig.savefig(OUT / ("%s.%s" % (stem, ext)), bbox_inches="tight", pad_inches=0.03)
    plt.close(fig)


print("\nstandalone panels (shared axes):")
for tag, _, title, stem in PANELS:
    fig, ax = plt.subplots(figsize=(2.6, 2.7))
    draw(ax, tag, title)
    ax.set_ylabel("\u2212log10 FDR", fontsize=7.5)
    fig.tight_layout()
    save(fig, stem)

fig, axes = plt.subplots(1, 2, figsize=(5.2, 2.7))
for ax, (tag, _, title, _) in zip(axes, PANELS):
    draw(ax, tag, title)
axes[0].set_ylabel("\u2212log10 FDR", fontsize=7.5)
fig.tight_layout()
save(fig, "volcano_pair_IFNg_TNFa")

print("\nwrote to:", OUT)
for p in sorted(OUT.glob("*")):
    print("  ", p.name)
