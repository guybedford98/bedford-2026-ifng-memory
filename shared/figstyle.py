"""Manuscript figure style: Arial, two font sizes, the condition and cluster palettes,
Prism axes, the bar convention and vector export.

    import sys; sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
    import figstyle as fs
    fs.setup()
    fig, ax = plt.subplots(figsize=fs.PANEL_QUARTER)
    info = fs.bars(ax, [ctrl, treat], colors=[fs.COLOR["naive"], fs.COLOR["stimulated"]])
    fs.prism(ax)
    fs.save(fig, out_dir / "Fig1E")

Two sizes only: 8 pt for titles and axis labels, 6 pt for everything smaller. Do not
introduce in-between sizes. pdf.fonttype=42 and svg.fonttype="none" keep all text as
editable Arial glyphs; use save() and never rasterize text.
"""
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# Condition palette for the four IFNg training conditions plus the two CRISPRa arms.
COLOR = {
    "naive":       "#90959a",
    "rested":      "#4b4b4c",
    "stimulated":  "#4f9bcf",
    "rechallenge": "#1e5d9f",
    "crispra_nt":  "#9e2064",   # CRISPRa non-targeting control
    "crispra_tgt": "#462f87",   # CRISPRa targeting (dCas9-VP64 + guide)
}
# Synonyms so a panel can index by whatever condition label its table carries.
COLOR.update({
    "stim":         COLOR["stimulated"],
    "rest":         COLOR["rested"],
    "rechallenged": COLOR["rechallenge"],
    "rech":         COLOR["rechallenge"],
    "nt":           COLOR["crispra_nt"],
    "ntc":          COLOR["crispra_nt"],
    "targeting":    COLOR["crispra_tgt"],
    "target":       COLOR["crispra_tgt"],
})
LAB = {"naive": "Naive", "stimulated": "Stim", "rested": "Rested",
       "rechallenge": "Rech", "crispra_nt": "NT", "crispra_tgt": "dCas9-VP64"}
NEUTRAL = "#6f6f6f"   # generic descriptive-panel gray (identity markers etc.)

# Cluster palette keyed by DISPLAY number (1-4); clusters are 0-indexed internally.
CLUSTER_COLORS = {1: "#4C72B0", 2: "#DD8452", 3: "#55A868", 4: "#C44E52"}

CONDITIONS = ["naive", "stimulated", "rested", "rechallenge"]  # IFNg training order

MM = 1 / 25.4
PANEL_HALF    = (55 * MM, 44 * MM)
PANEL_QUARTER = (40 * MM, 40 * MM)
PANEL_WIDE    = (85 * MM, 45 * MM)


def panel_mm(w_mm, h_mm):
    """Figure size in inches from millimetres."""
    return (w_mm * MM, h_mm * MM)


# Bar anatomy (notebook house rule).
BAR_EDGE = "#2b2b2b"
BAR_EDGE_LW = 0.6
MARKER_S = 17
MARKER_LW = 0.5
JITTER = 0.075
ERR_LW = 0.6
CAPSIZE = 1.6


def setup(base=8):
    """rcParams: Arial, 8 pt titles/axis labels, 6 pt everything smaller, vector text."""
    for fp in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\arialbd.ttf"):
        if Path(fp).exists():
            fm.fontManager.addfont(fp)
    font = "Arial" if "Arial" in {f.name for f in fm.fontManager.ttflist} else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": [font, "Liberation Sans", "Helvetica", "DejaVu Sans"],
        "font.size": base, "axes.titlesize": base, "axes.labelsize": base,
        "legend.fontsize": base - 2, "legend.title_fontsize": base - 2,
        "figure.titlesize": base,
        "xtick.labelsize": base - 2, "ytick.labelsize": base - 2,
        "svg.fonttype": "none",   # keep SVG text as editable text
        "pdf.fonttype": 42,       # keep PDF text as editable Arial glyphs
        "ps.fonttype": 42,
        "axes.linewidth": 0.8,
        "xtick.major.width": 0.8, "ytick.major.width": 0.8,
        "lines.linewidth": 0.8,
        "savefig.dpi": 600,
    })
    return font


def prism(ax, lw=0.8, ticklen=2.6):
    """Drop the top and right spines, ticks out."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(lw)
    ax.spines["bottom"].set_linewidth(lw)
    ax.tick_params(length=ticklen, width=lw, direction="out")


def clabel(cluster_label):
    """Internal 0-indexed label ('cluster_2') to display number string ('3')."""
    return str(int(str(cluster_label).split("_")[1]) + 1)


def save(fig, stem, svg=True, png=True):
    """Write PNG at 600 dpi plus PDF and SVG. stem = path without extension."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    out = []
    p = stem.with_suffix(".pdf")
    fig.savefig(p, bbox_inches="tight", pad_inches=0.02)
    out.append(p)
    if svg:
        p = stem.with_suffix(".svg")
        fig.savefig(p, bbox_inches="tight", pad_inches=0.02)
        out.append(p)
    if png:
        p = stem.with_suffix(".png")
        fig.savefig(p, dpi=600, bbox_inches="tight", pad_inches=0.02)
        out.append(p)
    plt.close(fig)
    return out


def summarize(x):
    """mean, SD, SEM, n for one replicate array (ddof=1)."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    n = x.size
    sd = x.std(ddof=1) if n > 1 else 0.0
    return {"mean": x.mean() if n else np.nan, "sd": sd,
            "sem": sd / np.sqrt(n) if n > 1 else 0.0, "n": n}


def bars(ax, series, x=None, width=0.72, colors=None, dots=True, sem=True,
         dot_size=MARKER_S, dot_edge=BAR_EDGE, dot_face=None, dot_lw=MARKER_LW,
         bar_edge=BAR_EDGE, bar_lw=BAR_EDGE_LW, jitter=JITTER, cap=CAPSIZE,
         err_lw=ERR_LW, zero_floor=True):
    """Draw one bar per replicate array.

    Bar height = group mean; whisker = mean +/- SEM (set sem=False for mean only,
    e.g. descriptive n=2 CPM panels); individual replicate dots overlaid with
    deterministic jitter (reproducible, no RNG). Each dot is filled with its
    bar's color (dot_face=None, the default) so dots match the bars.

    series : list of 1D replicate arrays (one per bar).
    x      : bar centres (default 0..k-1). colors : list of hex (default gray).
    Returns dict{x, means, sems, ns, ytop} where ytop[i] is the top of bar i's
    error whisker (use to anchor significance brackets)."""
    k = len(series)
    if x is None:
        x = np.arange(k, dtype=float)
    x = np.asarray(x, float)
    if colors is None:
        colors = [NEUTRAL] * k
    means, sems, ns, ytop = [], [], [], []
    for i, arr in enumerate(series):
        s = summarize(arr)
        m, se, n = s["mean"], s["sem"], s["n"]
        ax.bar(x[i], m, width=width, color=colors[i], edgecolor=bar_edge,
               linewidth=bar_lw, zorder=2)
        top = m
        if sem and n > 1 and se > 0:
            ax.errorbar(x[i], m, yerr=se, fmt="none", ecolor=bar_edge,
                        elinewidth=err_lw, capsize=cap, capthick=err_lw, zorder=3)
            top = m + se
        if dots:
            arr = np.asarray(arr, float)
            arr = arr[np.isfinite(arr)]
            nn = arr.size
            offs = np.linspace(-jitter, jitter, nn) if nn > 1 else np.array([0.0])
            face = colors[i] if dot_face is None else dot_face
            ax.scatter(x[i] + offs, arr, s=dot_size, facecolor=face,
                       edgecolor=dot_edge, linewidth=dot_lw, zorder=4, clip_on=False)
            top = max(top, float(arr.max()) if nn else top)
        means.append(m)
        sems.append(se)
        ns.append(n)
        ytop.append(top)
    if zero_floor:
        lo, hi = ax.get_ylim()
        if lo > 0:
            ax.set_ylim(0, hi)
    return {"x": x, "means": np.array(means), "sems": np.array(sems),
            "ns": ns, "ytop": np.array(ytop)}


def sig_bracket(ax, x1, x2, y, text, h=None, lw=0.9, color="black",
                fs_text=6, pad_frac=0.012):
    """Significance bracket from x1 to x2 at height y with `text` (e.g. '**', 'ns').
    h = tick height (default 2% of the y-range). Returns the y of the label."""
    lo, hi = ax.get_ylim()
    rng = hi - lo
    if h is None:
        h = rng * 0.02
    ax.plot([x1, x1, x2, x2], [y, y + h, y + h, y], lw=lw, color=color,
            clip_on=False, zorder=5)
    ax.text((x1 + x2) / 2, y + h + rng * pad_frac, text, ha="center",
            va="bottom", fontsize=fs_text, color=color, clip_on=False, zorder=6)
    return y + h + rng * pad_frac


def headroom(ax, frac=0.16):
    """Add vertical headroom above the tallest element for brackets and asterisks."""
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo, hi + (hi - lo) * frac)
