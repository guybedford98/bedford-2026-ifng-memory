#!/usr/bin/env python3
"""
tracklib.py  --  plotting engine for the multi-mark CUT&RUN track panels.

One module renders every layout so IDO1 and CXCL11 come out visually identical.
Signal is read with pybigtools; gene models are baked hg38 coordinates from
config/gene_models.json.

Design notes:
  * Manuscript style: Arial, 6 to 8 pt text only, thin Prism axes, top/right
    spines removed. Exports vector PDF + SVG (text stays editable) plus PNG.
  * Per-mark y axes always: marks have very different dynamic ranges (a sharp
    H3K4me3 promoter peak vs broad low H3K27me3), so no cross-mark shared axis.
"""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

# --------------------------------------------------------------------- style
def setup(base=7):
    for fp in (r"C:\Windows\Fonts\arial.ttf", r"C:\Windows\Fonts\arialbd.ttf"):
        if Path(fp).exists():
            fm.fontManager.addfont(fp)
    font = "Arial" if "Arial" in {f.name for f in fm.fontManager.ttflist} else "DejaVu Sans"
    plt.rcParams.update({
        "font.family": font, "font.size": base,
        "axes.titlesize": base + 1, "axes.labelsize": base,
        "xtick.labelsize": base - 1, "ytick.labelsize": base - 1,
        "legend.fontsize": base - 1,
        "svg.fonttype": "none",     # keep SVG text as text
        "pdf.fonttype": 42,         # keep PDF text as editable Arial glyphs
        "ps.fonttype": 42,
        "axes.linewidth": 0.7, "savefig.dpi": 600,
        "lines.antialiased": True,
    })
    return font

TITLE_FS = 8
LABEL_FS = 7
SMALL_FS = 6

# --------------------------------------------------------------------- signal
def read_signal(bw_path, chrom, start, end, nbins):
    """Binned mean coverage over the window as a length-nbins array. Missing -> 0."""
    import pybigtools
    bw = pybigtools.open(str(bw_path))
    try:
        vals = bw.values(chrom, int(start), int(end), bins=int(nbins),
                         summary="mean", exact=True, missing=0.0)
    finally:
        bw.close()
    return np.asarray(vals, dtype=float)


def positions(start, end, nbins):
    # bin centres in genomic coordinates
    edges = np.linspace(start, end, nbins + 1)
    return 0.5 * (edges[:-1] + edges[1:])


# --------------------------------------------------------------------- helpers
def condition_color(cond, palette):
    return cond["cb_color"] if palette == "cb_safe" else cond["color"]


def _lum(hex_color):
    r = int(hex_color[1:3], 16) / 255
    g = int(hex_color[3:5], 16) / 255
    b = int(hex_color[5:7], 16) / 255
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def mark_ymax(mark, signals):
    peak = 0.0
    for cond in mark["conditions"]:
        v = signals.get((mark["key"], cond))
        if v is not None and v.size:
            peak = max(peak, float(np.nanmax(v)))
    return peak if peak > 0 else 1.0


def _fmt_tick(v):
    if v >= 10:
        return f"{v:.0f}"
    if v >= 1:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _clean_signal_ax(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_visible(False)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", length=2.2, width=0.6, pad=1.5)
    ax.margins(x=0)


def _draw_signal(ax, pos, vals, color, fill, ymax, lw=0.8):
    if fill:
        ax.fill_between(pos, vals, color=color, alpha=0.32, linewidth=0)
        ax.plot(pos, vals, color=color, linewidth=lw, solid_capstyle="round")
    else:
        ax.plot(pos, vals, color=color, linewidth=lw + 0.1, solid_capstyle="round")
    ax.set_ylim(0, ymax * 1.15)


def _yticks(ax, ymax):
    ax.set_yticks([0, ymax])
    ax.set_yticklabels(["0", _fmt_tick(ymax)], fontsize=SMALL_FS)


# --------------------------------------------------------------------- gene model
def draw_gene_model(ax, gm, window, label=True):
    """Exon/intron gene model with CDS vs UTR heights, TSS marker, strand chevrons,
    a scale bar, and light kb ticks. Genomic coords increase left to right."""
    start, end = window["start"], window["end"]
    span = end - start
    ax.set_xlim(start, end)
    ax.set_ylim(-1.0, 1.35)
    for sp in ("top", "right", "left"):
        ax.spines[sp].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.6)
    ax.tick_params(left=False, labelleft=False)

    tx0, tx1 = gm["tx_start"], gm["tx_end"]
    cds0, cds1 = gm.get("cds_start", tx0), gm.get("cds_end", tx1)
    strand = gm["strand"]

    # intron line across the transcript
    ax.plot([tx0, tx1], [0, 0], color="#333333", linewidth=0.9,
            solid_capstyle="round", zorder=2)

    # strand chevrons along the body
    marker = ">" if strand == "+" else "<"
    for xp in np.linspace(tx0, tx1, 9)[1:-1]:
        ax.plot(xp, 0, marker=marker, markersize=3, color="#666666",
                markeredgewidth=0, zorder=3, alpha=0.8)

    # exon boxes: full height inside CDS, half height for UTR
    es, ee = gm.get("exon_starts", [tx0]), gm.get("exon_ends", [tx1])
    for s, e in zip(es, ee):
        # UTR portions (outside CDS)
        for a, b in ((s, min(e, cds0)), (max(s, cds1), e)):
            if b > a:
                ax.add_patch(mpatches.Rectangle((a, -0.16), b - a, 0.32,
                             facecolor="#333333", edgecolor="none", zorder=4))
        # CDS portion
        a, b = max(s, cds0), min(e, cds1)
        if b > a:
            ax.add_patch(mpatches.Rectangle((a, -0.34), b - a, 0.68,
                         facecolor="#333333", edgecolor="none", zorder=4))

    # TSS marker (small arrow bent in the direction of transcription)
    tss = gm["tss"]
    ax.plot([tss, tss], [0.0, 0.95], color="#c1272d", linewidth=0.8, zorder=5)
    dx = span * 0.028 * (1 if strand == "+" else -1)
    ax.annotate("", xy=(tss + dx, 0.95), xytext=(tss, 0.95),
                arrowprops=dict(arrowstyle="-|>", color="#c1272d", lw=0.8,
                                mutation_scale=6), zorder=5)
    ax.text(tss, 1.12, "TSS", ha="center", va="bottom", fontsize=SMALL_FS,
            color="#c1272d", fontweight="bold")

    # gene name
    if label:
        cen = (tx0 + tx1) / 2
        cen = min(max(cen, start + span * 0.05), end - span * 0.05)
        ax.text(cen, -0.72, gm.get("name_label", ""), ha="center", va="top",
                fontsize=SMALL_FS + 0.5, fontstyle="italic", color="#222222")

    # scale bar (round number near 1/5 of the span), bottom right
    cand = [200, 500, 1000, 2000, 5000, 10000]
    bar = min(cand, key=lambda x: abs(x - span / 5))
    x1 = end - span * 0.03
    x0 = x1 - bar
    ax.plot([x0, x1], [-0.9, -0.9], color="#222222", linewidth=1.4,
            solid_capstyle="butt", zorder=5)
    lbl = f"{bar // 1000} kb" if bar >= 1000 else f"{bar} bp"
    ax.text((x0 + x1) / 2, -1.02, lbl, ha="center", va="top",
            fontsize=SMALL_FS, color="#222222")

    # a couple of kb ticks for orientation
    ticks = np.linspace(start, end, 4)
    ax.set_xticks(ticks)
    ax.set_xticklabels([f"{t/1e6:.3f}" for t in ticks], fontsize=SMALL_FS)
    ax.tick_params(axis="x", length=2, width=0.5, pad=1.5)
    ax.set_xlabel(f"{window['chrom']} (Mb)", fontsize=SMALL_FS, labelpad=1)


def _empty_xaxis_scalebar(ax, window):
    """When gene model is off, draw only a scale bar + minimal axis on the bottom track."""
    start, end = window["start"], window["end"]
    span = end - start
    cand = [200, 500, 1000, 2000, 5000, 10000]
    bar = min(cand, key=lambda x: abs(x - span / 5))
    x1 = end - span * 0.03
    x0 = x1 - bar
    y = ax.get_ylim()[1] * 0.9
    ax.plot([x0, x1], [y, y], color="#222222", linewidth=1.4, solid_capstyle="butt", zorder=6)
    lbl = f"{bar // 1000} kb" if bar >= 1000 else f"{bar} bp"
    ax.text((x0 + x1) / 2, y * 1.02, lbl, ha="center", va="bottom",
            fontsize=SMALL_FS, color="#222222")


# --------------------------------------------------------------------- legend
def _condition_legend(fig, conditions, palette, fill, anchor=(0.995, 0.5), loc="center left"):
    handles = []
    for c in conditions:
        col = condition_color(c, palette)
        if fill:
            handles.append(mpatches.Patch(facecolor=col, alpha=0.55, edgecolor=col,
                                          linewidth=0.8, label=c["label"]))
        else:
            handles.append(plt.Line2D([0], [0], color=col, linewidth=1.4, label=c["label"]))
    leg = fig.legend(handles=handles, loc=loc, bbox_to_anchor=anchor,
                     frameon=False, fontsize=SMALL_FS, handlelength=1.2,
                     handletextpad=0.5, labelspacing=0.4, borderaxespad=0.2)
    return leg


def _suptitle(fig, title, subtitle):
    fig.text(0.01, 0.995, title, ha="left", va="top", fontsize=TITLE_FS, fontweight="bold")
    if subtitle:
        fig.text(0.01, 0.965, subtitle, ha="left", va="top", fontsize=SMALL_FS, color="#555555")


# --------------------------------------------------------------------- layouts
def _overlay(gm, marks, conditions, signals, window, *, fill, palette,
             show_gene_model, title, subtitle, grouped=False):
    """Rows = marks, all conditions overlaid on one shared per-mark axis."""
    nrows = len(marks) + (1 if show_gene_model else 0)
    ratios = [3.0] * len(marks) + ([1.7] if show_gene_model else [])
    fig_h = 0.62 * len(marks) + (0.75 if show_gene_model else 0.25) + 0.55
    left = 0.205 if grouped else 0.145
    mlx = -0.04 if grouped else -0.03
    top = 0.90
    fig = plt.figure(figsize=(3.6, fig_h))
    gs = GridSpec(nrows, 1, height_ratios=ratios, hspace=0.14 if grouped else 0.08,
                  left=left, right=0.80, top=top, bottom=0.11)
    pos = positions(window["start"], window["end"], signals["_nbins"])

    axes = []
    for i, mark in enumerate(marks):
        ax = fig.add_subplot(gs[i, 0])
        axes.append(ax)
        ax.set_xlim(window["start"], window["end"])
        ymax = mark_ymax(mark, signals)
        for cond in conditions:
            v = signals.get((mark["key"], cond["key"]))
            if v is None:
                continue
            _draw_signal(ax, pos, v, condition_color(cond, palette), fill, ymax)
        _clean_signal_ax(ax)
        _yticks(ax, ymax)
        ax.text(mlx, 0.5, mark["label"], transform=ax.transAxes, ha="right",
                va="center", fontsize=LABEL_FS, fontweight="bold", rotation=0)
        # faint availability note for the 2-condition marks
        navail = sum(1 for c in conditions if (mark["key"], c["key"]) in signals)
        if navail < len(conditions):
            ax.text(0.985, 0.86, "Naive, Stim only", transform=ax.transAxes,
                    ha="right", va="top", fontsize=SMALL_FS - 0.5, color="#999999")

    if grouped:
        _annotate_state_groups(fig, axes, marks, left, 0.80)

    if show_gene_model:
        axg = fig.add_subplot(gs[len(marks), 0])
        draw_gene_model(axg, gm, window)
    else:
        axes[-1].tick_params(axis="x", bottom=True, labelbottom=True, length=2, width=0.5)
        axes[-1].spines["bottom"].set_visible(True)
        axes[-1].spines["bottom"].set_linewidth(0.6)
        ticks = np.linspace(window["start"], window["end"], 4)
        axes[-1].set_xticks(ticks)
        axes[-1].set_xticklabels([f"{t/1e6:.3f}" for t in ticks], fontsize=SMALL_FS)
        axes[-1].set_xlabel(f"{window['chrom']} (Mb)", fontsize=SMALL_FS, labelpad=1)
        _empty_xaxis_scalebar(axes[0], window)

    _condition_legend(fig, conditions, palette, fill)
    _suptitle(fig, title, subtitle)
    return fig


def _annotate_state_groups(fig, axes, marks, left, right):
    """Active vs repressive: a dashed separator at the group boundary plus a
    rotated group label on a far-left rail, clear of the mark names."""
    import itertools
    n_active = sum(1 for m in marks if m["state"] == "active")
    if 0 < n_active < len(marks):
        yb = (axes[n_active - 1].get_position().y0 + axes[n_active].get_position().y1) / 2
        fig.add_artist(plt.Line2D([left, right], [yb, yb], transform=fig.transFigure,
                                  color="#bbbbbb", linewidth=0.8, linestyle=(0, (4, 2))))
    xl = max(left - 0.185, 0.02)
    for state, group in itertools.groupby(enumerate(marks), key=lambda t: t[1]["state"]):
        idx = [i for i, _ in group]
        gtop = axes[idx[0]].get_position().y1
        gbot = axes[idx[-1]].get_position().y0
        ymid = (gtop + gbot) / 2
        fig.add_artist(plt.Line2D([xl, xl], [gbot, gtop], transform=fig.transFigure,
                                  color="#888888", linewidth=1.0))
        fig.text(xl - 0.008, ymid, "Active" if state == "active" else "Repressive",
                 ha="right", va="center", rotation=90, fontsize=SMALL_FS,
                 color="#555555", fontweight="bold")


def _stacked(gm, marks, conditions, signals, window, *, fill, palette, y_scaling,
             show_gene_model, title, subtitle):
    """Rows = marks, each condition as its own stacked subtrack (no overlay)."""
    blocks = [(m, [c for c in conditions if (m["key"], c["key"]) in signals]) for m in marks]
    nsub = sum(len(cs) for _, cs in blocks)
    nrows = nsub + (1 if show_gene_model else 0)
    ratios = [1.0] * nsub + ([1.7] if show_gene_model else [])
    fig_h = 0.30 * nsub + (0.75 if show_gene_model else 0.25) + 0.55
    fig = plt.figure(figsize=(3.5, fig_h))
    gs = GridSpec(nrows, 1, height_ratios=ratios, hspace=0.06,
                  left=0.145, right=0.80, top=0.92, bottom=0.09)
    pos = positions(window["start"], window["end"], signals["_nbins"])

    row = 0
    first_axes = []
    for mark, conds in blocks:
        ymark = mark_ymax(mark, signals)
        block_axes = []
        for j, cond in enumerate(conds):
            ax = fig.add_subplot(gs[row, 0]); row += 1
            block_axes.append(ax)
            ax.set_xlim(window["start"], window["end"])
            v = signals.get((mark["key"], cond["key"]))
            ymax = ymark if y_scaling == "per_mark" else (float(np.nanmax(v)) if v.size and np.nanmax(v) > 0 else 1.0)
            _draw_signal(ax, pos, v, condition_color(cond, palette), fill, ymax)
            _clean_signal_ax(ax)
            _yticks(ax, ymax)
            col = condition_color(cond, palette)
            lblc = col if _lum(col) < 0.62 else "#555555"
            ax.text(0.012, 0.82, cond["label"], transform=ax.transAxes, ha="left",
                    va="top", fontsize=SMALL_FS, color=lblc, fontweight="bold",
                    bbox=dict(facecolor="white", alpha=0.55, edgecolor="none", pad=0.4))
        # mark bracket label spanning the block
        top = block_axes[0].get_position().y1
        bot = block_axes[-1].get_position().y0
        fig.text(0.02, (top + bot) / 2, mark["label"], ha="center", va="center",
                 rotation=90, fontsize=LABEL_FS, fontweight="bold")
        first_axes.append(block_axes[0])

    if show_gene_model:
        axg = fig.add_subplot(gs[nsub, 0])
        draw_gene_model(axg, gm, window)
    else:
        last = fig.axes[nsub - 1]
        last.tick_params(axis="x", bottom=True, labelbottom=True, length=2, width=0.5)
        last.spines["bottom"].set_visible(True); last.spines["bottom"].set_linewidth(0.6)
        ticks = np.linspace(window["start"], window["end"], 4)
        last.set_xticks(ticks)
        last.set_xticklabels([f"{t/1e6:.3f}" for t in ticks], fontsize=SMALL_FS)
        last.set_xlabel(f"{window['chrom']} (Mb)", fontsize=SMALL_FS, labelpad=1)

    ysc = "shared per-mark y" if y_scaling == "per_mark" else "auto-scaled per track"
    _suptitle(fig, title, (subtitle + "   |   " + ysc) if subtitle else ysc)
    return fig


def _grid(gm, marks, conditions, signals, window, *, fill, palette, y_scaling,
          show_gene_model, title, subtitle):
    """Small multiples: rows = marks, columns = conditions."""
    ncol = len(conditions)
    nmark = len(marks)
    nrows = nmark + (1 if show_gene_model else 0)
    ratios = [1.0] * nmark + ([0.9] if show_gene_model else [])
    fig_w = 6.6
    fig_h = 0.52 * nmark + (0.62 if show_gene_model else 0.2) + 0.5
    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = GridSpec(nrows, ncol, height_ratios=ratios, hspace=0.12, wspace=0.12,
                  left=0.085, right=0.995, top=0.88, bottom=0.10)
    pos = positions(window["start"], window["end"], signals["_nbins"])

    for r, mark in enumerate(marks):
        ymark = mark_ymax(mark, signals)
        for c, cond in enumerate(conditions):
            ax = fig.add_subplot(gs[r, c])
            ax.set_xlim(window["start"], window["end"])
            v = signals.get((mark["key"], cond["key"]))
            if v is None:
                ax.set_ylim(0, 1)
                ax.text(0.5, 0.5, "n.a.", transform=ax.transAxes, ha="center",
                        va="center", fontsize=SMALL_FS, color="#bbbbbb")
                for sp in ax.spines.values():
                    sp.set_visible(False)
                ax.set_xticks([]); ax.set_yticks([])
                continue
            ymax = ymark if y_scaling == "per_mark" else (float(np.nanmax(v)) if v.size and np.nanmax(v) > 0 else 1.0)
            _draw_signal(ax, pos, v, condition_color(cond, palette), fill, ymax)
            _clean_signal_ax(ax)
            if c == 0:
                _yticks(ax, ymax)
                ax.text(-0.30, 0.5, mark["label"], transform=ax.transAxes, ha="right",
                        va="center", fontsize=LABEL_FS, fontweight="bold")
            else:
                if y_scaling == "per_mark":
                    ax.set_yticks([])
                    ax.spines["left"].set_visible(False)
                else:
                    _yticks(ax, ymax)
            if r == 0:
                ax.set_title(cond["label"], fontsize=LABEL_FS, fontweight="bold", pad=3)

    if show_gene_model:
        axg = fig.add_subplot(gs[nmark, :])
        draw_gene_model(axg, gm, window)

    ysc = "shared per-mark y" if y_scaling == "per_mark" else "auto-scaled per cell"
    _suptitle(fig, title, (subtitle + "   |   " + ysc) if subtitle else ysc)
    return fig


def render_region_stacked(gm, blocks, conditions, variant, toggle, subtitle=""):
    """Stacked subtracks like Variant 2, but each mark is shown over its own
    genomic window (TSS +/- 3 kb or the gene body). Marks that share a window are
    grouped, and each group gets its own gene-model track and coordinate axis.

    blocks: ordered list of {mark, region, window, signals:{cond_key: vals}}.
    """
    fill = variant["fill"]
    palette = variant["palette"]
    y_scaling = toggle.get("y_scaling", "per_mark")
    show_gm = toggle.get("gene_model", True)
    title = variant["title"]
    region_lab = {"tss3kb": "TSS +/- 3 kb", "genebody": "gene body"}
    label_by_key = {b["mark"]["key"]: b["mark"]["label"] for b in blocks}

    # group blocks that share the same window (same region focus)
    groups = []
    for b in blocks:
        if groups and groups[-1]["win"] == b["window"]:
            groups[-1]["blocks"].append(b)
        else:
            groups.append({"win": b["window"], "region": b["region"], "blocks": [b]})

    # flat row plan: one row per condition subtrack, one gene-model row per group,
    # a spacer row between groups.
    rows = []
    for gi, g in enumerate(groups):
        for b in g["blocks"]:
            mark = b["mark"]
            conds = [c for c in conditions if c["key"] in mark["conditions"]]
            ymark = max((float(np.nanmax(b["signals"][c["key"]]))
                         for c in conds if b["signals"][c["key"]].size), default=0.0)
            ymark = ymark if ymark > 0 else 1.0
            for c in conds:
                rows.append({"t": "trk", "mark": mark, "cond": c, "win": g["win"],
                             "vals": b["signals"][c["key"]], "ymark": ymark, "gi": gi})
        if show_gm:
            rows.append({"t": "gm", "win": g["win"], "region": g["region"], "gi": gi})
        if gi < len(groups) - 1:
            rows.append({"t": "gap", "gi": gi})

    ratios = [1.0 if r["t"] == "trk" else (1.7 if r["t"] == "gm" else 0.5) for r in rows]
    ntrk = sum(1 for r in rows if r["t"] == "trk")
    fig_h = 0.30 * ntrk + 0.78 * len(groups) + 0.5 * (len(groups) - 1) + 0.6
    left = 0.225
    fig = plt.figure(figsize=(3.7, fig_h))
    gs = GridSpec(len(rows), 1, height_ratios=ratios, hspace=0.06,
                  left=left, right=0.80, top=0.93, bottom=0.06)

    axes = []
    axes_by_mark = {}
    for i, r in enumerate(rows):
        ax = fig.add_subplot(gs[i, 0])
        axes.append(ax)
        if r["t"] == "gap":
            ax.axis("off")
            continue
        if r["t"] == "gm":
            draw_gene_model(ax, gm, r["win"], label=(r["region"] == "genebody"))
            continue
        win = r["win"]
        pos = positions(win["start"], win["end"], r["vals"].size)
        ax.set_xlim(win["start"], win["end"])
        v = r["vals"]
        ymax = (r["ymark"] if y_scaling == "per_mark"
                else (float(np.nanmax(v)) if v.size and np.nanmax(v) > 0 else 1.0))
        _draw_signal(ax, pos, v, condition_color(r["cond"], palette), fill, ymax)
        _clean_signal_ax(ax)
        _yticks(ax, ymax)
        col = condition_color(r["cond"], palette)
        lblc = col if _lum(col) < 0.62 else "#555555"
        ax.text(0.012, 0.82, r["cond"]["label"], transform=ax.transAxes, ha="left",
                va="top", fontsize=SMALL_FS, color=lblc, fontweight="bold",
                bbox=dict(facecolor="white", alpha=0.55, edgecolor="none", pad=0.4))
        axes_by_mark.setdefault(r["mark"]["key"], []).append(ax)

    # mark bracket labels spanning each mark's subtracks
    for mkey, maxes in axes_by_mark.items():
        top = maxes[0].get_position().y1
        bot = maxes[-1].get_position().y0
        fig.text(0.108, (top + bot) / 2, label_by_key[mkey], ha="center", va="center",
                 rotation=90, fontsize=LABEL_FS, fontweight="bold")

    # region rail on the far left, spanning each group's tracks + gene model
    for gi in range(len(groups)):
        ga = [axes[i] for i, r in enumerate(rows) if r["gi"] == gi and r["t"] in ("trk", "gm")]
        top = max(a.get_position().y1 for a in ga)
        bot = min(a.get_position().y0 for a in ga)
        xr = 0.02
        fig.add_artist(plt.Line2D([xr, xr], [bot, top], transform=fig.transFigure,
                                  color="#888888", linewidth=1.0))
        fig.text(xr - 0.008, (top + bot) / 2, region_lab[groups[gi]["region"]],
                 ha="center", va="center", rotation=90, fontsize=SMALL_FS,
                 color="#555555", fontweight="bold")

    _condition_legend(fig, conditions, palette, fill)
    ysc = "shared per-mark y" if y_scaling == "per_mark" else "auto-scaled per track"
    _suptitle(fig, title, (subtitle + "   |   " + ysc) if subtitle else ysc)
    return fig


# --------------------------------------------------------------------- dispatch
def render(gene_model, marks, conditions, signals, window, variant, toggle,
           subtitle=""):
    layout = variant["layout"]
    fill = variant["fill"]
    palette = variant["palette"]
    y_scaling = toggle.get("y_scaling", "per_mark")
    show_gm = toggle.get("gene_model", True)
    title = variant["title"]

    if layout == "overlay":
        return _overlay(gene_model, marks, conditions, signals, window, fill=fill,
                        palette=palette, show_gene_model=show_gm, title=title,
                        subtitle=subtitle)
    if layout == "grouped_state":
        return _overlay(gene_model, marks, conditions, signals, window, fill=fill,
                        palette=palette, show_gene_model=show_gm, title=title,
                        subtitle=subtitle, grouped=True)
    if layout == "stacked":
        return _stacked(gene_model, marks, conditions, signals, window, fill=fill,
                        palette=palette, y_scaling=y_scaling, show_gene_model=show_gm,
                        title=title, subtitle=subtitle)
    if layout == "grid":
        return _grid(gene_model, marks, conditions, signals, window, fill=fill,
                     palette=palette, y_scaling=y_scaling, show_gene_model=show_gm,
                     title=title, subtitle=subtitle)
    raise ValueError(f"unknown layout {layout!r}")


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 300})):
        try:
            fig.savefig(stem.with_suffix(ext), bbox_inches="tight", pad_inches=0.02, **kw)
        except PermissionError:
            print(f"    [locked, kept existing] {stem.name}{ext}")
    plt.close(fig)
    return stem.with_suffix(".pdf")
