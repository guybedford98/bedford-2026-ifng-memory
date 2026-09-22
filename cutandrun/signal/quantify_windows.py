#!/usr/bin/env python3
"""
quantify_windows.py  --  mean CUT&RUN CPM per mark per condition: Fig 3B (CXCL11) and
Fig S6B (IDO1).

For each mark, compute the mean CPM within that mark's DISPLAYED track window
(H3K4me3, H3K36me2 -> TSS +/- 3 kb; H3K27me3, H3K36me3 -> gene body, with the
CXCL11 gene-centered 15 kb override), per condition, from the PER-REPLICATE
bigwigs (n = 2, R1/R2). Writes tidy CSVs, then renders the bar-graph panels
(mean + both replicate dots).

The Stim_H3K36me3_SPARQ_R1 methods-comparison library is deliberately excluded
to keep n = 2 balanced.

Windows come from ../tracks/config, so the bars match the rendered tracks exactly.

Run:  py quantify_windows.py            (quantify if needed, then plot both genes)
      py quantify_windows.py --force    (re-read bigwigs even if the CSV exists)
"""
import csv
import json
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
sys.path.insert(0, str(ROOT.parent / "tracks"))
import config as cfg  # noqa: E402
import common as cm  # noqa: E402
import tracklib as tl  # noqa: E402

TRACKS_CONFIG = ROOT.parent / "tracks" / "config"
CONFIG = json.loads((TRACKS_CONFIG / "config.json").read_text(encoding="utf-8"))
MODELS = json.loads((TRACKS_CONFIG / "gene_models.json").read_text(encoding="utf-8"))

DATA_DIR = ROOT / "data"

# condition -> manuscript color
COND_COLOR = {"Naive": "#90959a", "Stim": "#4f9bcf", "Rested": "#4b4b4c", "Rechallenged": "#1e5d9f"}
COND_LABEL = {"Naive": "Naive", "Stim": "Stim", "Rested": "Rested", "Rechallenged": "Rechal"}
REPS = ["R1", "R2"]

# gene -> manuscript panel
PANELS = {"CXCL11": "Fig3B", "IDO1": "FigS6B"}

# per-mark quantification spec: which window (region), per-replicate and merged bigwig
# filename templates (the merged value is a cross-check), and conditions.
MARKS = [
    {"key": "H3K4me3",  "region": "tss3kb",
     "perrep_tmpl": "{cond}_H3K4me3_{rep}_cpm.bw",
     "merged_tmpl": "{cond}_H3K4me3_merged_cpm.bw",
     "conditions": ["Naive", "Stim", "Rested", "Rechallenged"]},
    {"key": "H3K36me2", "region": "tss3kb",
     "perrep_tmpl": "{cond}_H3K36me2_{rep}_cpm.bw",
     "merged_tmpl": "{cond}_K36me2_merged_cpm.bw",
     "conditions": ["Naive", "Stim"]},
    {"key": "H3K27me3", "region": "genebody",
     "perrep_tmpl": "{cond}_H3K27me3_{rep}_cpm.bw",
     "merged_tmpl": "{cond}_H3K27me3_merged_cpm.bw",
     "conditions": ["Naive", "Stim", "Rested", "Rechallenged"]},
    {"key": "H3K36me3", "region": "genebody",
     "perrep_tmpl": "{cond}_H3K36me3_{rep}_cpm.bw",
     "merged_tmpl": "{cond}_K36me3_merged_cpm.bw",
     "conditions": ["Naive", "Stim"]},
]
MARK_ORDER = [m["key"] for m in MARKS]


# --------------------------------------------------------------------- windows
def region_window(gm, region):
    """Identical to run_tracks.region_window so bars match the displayed tracks."""
    rw = CONFIG.get("region_windows", {})
    tss = gm["tss"]
    if region == "tss3kb":
        hw = rw.get("tss3kb_halfwidth", 3000)
        return {"chrom": gm["chrom"], "start": int(tss - hw), "end": int(tss + hw)}
    ov = rw.get("genebody_override", {}).get(gm.get("name_label", ""))
    if ov and ov.get("mode") == "gene_center":
        center = (gm["tx_start"] + gm["tx_end"]) // 2
        half = int(ov["span"]) // 2
        return {"chrom": gm["chrom"], "start": int(center - half), "end": int(center + half)}
    return dict(gm["window"])


def window_mean(bw_path, win):
    """Mean CPM over the whole window (missing -> 0, matching the track render)."""
    v = tl.read_signal(bw_path, win["chrom"], win["start"], win["end"], 1)
    return float(v[0]) if v.size else float("nan")


# --------------------------------------------------------------------- quantify
def quantify():
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    long_rows = []      # per replicate
    summ_rows = []      # per condition (mean/sem + merged cross-check)

    for gene in ("IDO1", "CXCL11"):
        gm = dict(MODELS["genes"][gene])
        gm["name_label"] = gene
        for mark in MARKS:
            win = region_window(gm, mark["region"])
            span = win["end"] - win["start"]
            print(f"\n[{gene} {mark['key']}] {mark['region']}  "
                  f"{win['chrom']}:{win['start']:,}-{win['end']:,}  ({span/1000:.1f} kb)")
            for cond in mark["conditions"]:
                vals = []
                for rep in REPS:
                    p = cm.bigwig(mark["perrep_tmpl"].format(cond=cond, rep=rep))
                    if not p.exists():
                        print(f"    MISSING {p.name}")
                        continue
                    mv = window_mean(p, win)
                    vals.append(mv)
                    long_rows.append({"gene": gene, "mark": mark["key"], "region": mark["region"],
                                      "chrom": win["chrom"], "start": win["start"], "end": win["end"],
                                      "condition": cond, "replicate": rep,
                                      "mean_cpm": f"{mv:.6f}", "bigwig": p.name})
                    print(f"    {cond:13s} {rep}: {mv:.4f}")
                arr = np.array(vals, dtype=float)
                mean = float(np.nanmean(arr)) if arr.size else float("nan")
                sd = float(np.nanstd(arr, ddof=1)) if arr.size > 1 else float("nan")
                sem = sd / np.sqrt(arr.size) if arr.size > 1 else float("nan")
                # merged cross-check (mean of merged bigwig over same window)
                mp = cm.bigwig(mark["merged_tmpl"].format(cond=cond))
                merged = window_mean(mp, win) if mp.exists() else float("nan")
                summ_rows.append({"gene": gene, "mark": mark["key"], "region": mark["region"],
                                  "condition": cond, "n": arr.size,
                                  "mean_cpm": f"{mean:.6f}", "sem": f"{sem:.6f}",
                                  "rep_values": ";".join(f"{v:.4f}" for v in vals),
                                  "merged_cpm_check": f"{merged:.6f}"})

    with open(DATA_DIR / "window_quant_perrep.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(long_rows[0].keys()))
        w.writeheader(); w.writerows(long_rows)
    with open(DATA_DIR / "window_quant_summary.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summ_rows[0].keys()))
        w.writeheader(); w.writerows(summ_rows)
    print(f"\nwrote {DATA_DIR/'window_quant_perrep.csv'} ({len(long_rows)} rows)")
    print(f"wrote {DATA_DIR/'window_quant_summary.csv'} ({len(summ_rows)} rows)")
    return summ_rows, long_rows


def load_cached():
    def rd(p):
        with open(p, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    return rd(DATA_DIR / "window_quant_summary.csv"), rd(DATA_DIR / "window_quant_perrep.csv")


# --------------------------------------------------------------------- plotting
import matplotlib.pyplot as plt  # noqa: E402

# Fixed plotting-box aspect (height / width) for every stacked bar sub-panel, so
# all marks share one panel shape regardless of how many conditions they carry
# (H3K4me3/H3K27me3 have 4 bars, H3K36me2/H3K36me3 only Naive+Stim). Pinning the
# box means a single-mark, paired, or all-four export all drop into Illustrator at
# the identical tile size (previously the 2-condition K36 panels landed wider).
BOX_ASPECT_H = 0.4056   # horizontal bars
BOX_ASPECT_V = 0.4489   # vertical bars

# Every panel is drawn over the SAME number of category slots (the most any mark
# has: Naive/Stim/Rested/Rechallenged = 4). A mark with fewer conditions (the K36
# marks: Naive/Stim only) has its bars centred in that field, so the bar THICKNESS
# and spacing are identical across marks instead of the 2-bar panels rendering fat.
# The unused slots are left as blank whitespace, not labelled empty conditions.
BAR_THICK = 0.65        # bar thickness as a fraction of one slot (both orientations)
MAX_SLOTS = max(len(m["conditions"]) for m in MARKS)


def _mark_data(gene, mark_key, summ, long):
    conds = next(m for m in MARKS if m["key"] == mark_key)["conditions"]
    means, sems, pts = [], [], []
    for cond in conds:
        s = next((r for r in summ if r["gene"] == gene and r["mark"] == mark_key
                  and r["condition"] == cond), None)
        means.append(float(s["mean_cpm"]) if s else np.nan)
        sems.append(float(s["sem"]) if s and s["sem"] not in ("", "nan") else 0.0)
        pv = [float(r["mean_cpm"]) for r in long if r["gene"] == gene
              and r["mark"] == mark_key and r["condition"] == cond]
        pts.append(pv)
    return conds, means, sems, pts


def _dot_offsets(n):
    """Symmetric within-bar offsets for replicate dots."""
    if n <= 1:
        return [0.0]
    return list(np.linspace(-0.15, 0.15, n))


def _prism(ax):
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.spines["left"].set_linewidth(1.2)
    ax.spines["bottom"].set_linewidth(1.2)
    ax.tick_params(length=3.0, width=1.0)


def _bar_subpanel(ax, gene, mark_key, summ, long, *, horizontal=False):
    """Solid bars with a black edge; replicate dots filled in the bar colour with a
    black edge; no SEM whiskers (the dots carry the spread)."""
    conds, means, sems, pts = _mark_data(gene, mark_key, summ, long)
    colors = [COND_COLOR[c] for c in conds]
    labels = [COND_LABEL[c] for c in conds]
    n = len(conds)
    # Centre this mark's n bars in a fixed field of MAX_SLOTS category slots, so the
    # bar thickness is the same on every panel regardless of how many conditions it has.
    x = np.arange(n) + (MAX_SLOTS - n) / 2.0
    cat_lim = (-0.5, MAX_SLOTS - 0.5)
    flat = [v for pv in pts for v in pv]
    top = _round_top(max(means + flat + [1e-9]))

    if horizontal:
        ax.barh(x, means, color=colors, height=BAR_THICK, edgecolor="black",
                linewidth=1.0, zorder=2)
        for xi, c, pv in zip(x, conds, pts):
            if pv:
                ys = [xi + o for o in _dot_offsets(len(pv))]
                ax.plot(pv, ys, linestyle="None", marker="o", markersize=4.5,
                        markerfacecolor=COND_COLOR[c], markeredgecolor="black",
                        markeredgewidth=0.8, zorder=3, clip_on=False)
        ax.set_yticks(x); ax.set_yticklabels(labels, fontsize=tl.SMALL_FS)
        ax.set_ylim(*cat_lim)
        ax.invert_yaxis()
        ax.set_xlim(0, top * 1.02)
        ax.set_xticks([0, top])
        ax.set_xticklabels(["0", _fmt_top(top)], fontsize=tl.SMALL_FS)
        ax.set_xlabel("Mean CPM", fontsize=tl.SMALL_FS, labelpad=1)
    else:
        ax.bar(x, means, color=colors, width=BAR_THICK, edgecolor="black",
               linewidth=1.0, zorder=2)
        for xi, c, pv in zip(x, conds, pts):
            if pv:
                xs = [xi + o for o in _dot_offsets(len(pv))]
                ax.plot(xs, pv, linestyle="None", marker="o", markersize=4.5,
                        markerfacecolor=COND_COLOR[c], markeredgecolor="black",
                        markeredgewidth=0.8, zorder=3, clip_on=False)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=tl.SMALL_FS, rotation=40, ha="right")
        ax.set_xlim(*cat_lim)
        ax.set_ylim(0, top * 1.02)
        ax.set_yticks([0, top])
        ax.set_yticklabels(["0", _fmt_top(top)], fontsize=tl.SMALL_FS)
        ax.set_ylabel("Mean CPM", fontsize=tl.SMALL_FS, labelpad=1)
    _prism(ax)
    ax.set_title(mark_key, fontsize=tl.LABEL_FS, fontweight="bold", pad=2)
    # Pin the plotting box to one fixed aspect so every mark panel is the same
    # tile shape (the K36 marks no longer render wider than the 4-condition marks).
    ax.set_box_aspect(BOX_ASPECT_H if horizontal else BOX_ASPECT_V)


def _round_top(v):
    """Top-of-axis tick that sits just above the largest plotted point.

    A coarse ladder (1, 1.5, 2, 2.5, ...) overshoots badly on the marks whose
    maximum lands just past a step: CXCL11 H3K27me3 peaks at 0.101 and would be drawn
    on a 0.15 axis, and H3K4me3 peaks at 1.14 on a 1.5 axis, so those bars would stop
    two thirds of the way across while the K36 panels run the full width. Rounding
    up to two significant figures instead ends every axis at the highest point
    (<= 10 percent of headroom) and still gives a readable tick label."""
    if v <= 0:
        return 1.0
    import math
    step = 10 ** (math.floor(math.log10(v)) - 1)
    return round(math.ceil(v / step - 1e-9) * step, 10)


def _fmt_top(v):
    """Tick label for _round_top: enough decimals to print both significant digits
    (tracklib._fmt_tick is fixed at 2 dp and would render 0.032 as '0.03')."""
    import math
    if v <= 0:
        return "0"
    decimals = max(0, 1 - math.floor(math.log10(v)))
    return f"{v:.{decimals}f}"


def _stack_figure(gene, mark_keys, summ, long, *, horizontal=False, width=1.7,
                  per_h=0.95, tag=""):
    n = len(mark_keys)
    fig_h = per_h * n + 0.35
    fig, axes = plt.subplots(n, 1, figsize=(width, fig_h))
    if n == 1:
        axes = [axes]
    for ax, mk in zip(axes, mark_keys):
        _bar_subpanel(ax, gene, mk, summ, long, horizontal=horizontal)
    fig.subplots_adjust(left=0.30, right=0.97, top=0.93, bottom=0.16, hspace=0.75)
    return fig


def save(fig, stem):
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext, kw in ((".pdf", {}), (".svg", {}), (".png", {"dpi": 300})):
        fig.savefig(stem.with_suffix(ext), bbox_inches="tight", pad_inches=0.02, **kw)
    plt.close(fig)


def plot_all(summ, long):
    tl.setup()
    for gene, panel in PANELS.items():
        fig = _stack_figure(gene, MARK_ORDER, summ, long, horizontal=True, per_h=0.85)
        save(fig, cfg.out_dir(panel) / panel)
        print(f"    rendered {panel} ({gene})")


def main():
    force = "--force" in sys.argv[1:]
    if force or not (DATA_DIR / "window_quant_summary.csv").exists():
        summ, long = quantify()
    else:
        print("[cache] reading existing window_quant CSVs (use --force to re-read bigwigs)")
        summ, long = load_cached()
    plot_all(summ, long)
    print("\ndone.")


if __name__ == "__main__":
    main()
