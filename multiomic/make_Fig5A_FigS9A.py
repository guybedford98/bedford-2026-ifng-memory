#!/usr/bin/env python3
"""Fig 5A and Fig S9A: Naive H3K27me3 by memory class, all three classes in one run.

Fig 5A   TSS +/- 10 kb metagene of Naive H3K27me3, mean CPM per class, no spread shown.
Fig S9A  Naive promoter level, mean CPM over TSS +/- 3 kb, box and strip on a log10 axis.

Classes come from data/universe.tsv and coordinates from data/coordinates.tsv, so both panels
sit on the same three-class assignment and the same TSS frame as the rest of this folder.
Signal is the merged Naive H3K27me3 CPM bigWig. Per-gene profiles and focal values are cached
to data/naive_k27me3_3class.npz, keyed on the class gene lists, so restyling never re-reads it.

Class order in both panels is memory, refractory, non-memory, the order of decreasing signal.

Output: results/Fig5A/Fig5A.{pdf,png,svg} + .json
        results/FigS9A/FigS9A.{pdf,png,svg} + .json
        data/naive_k27me3_3class.npz, data/k27me3_focal_naive_3class.tsv
"""
import sys, json, random
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd, pybigtools
import matplotlib.pyplot as plt
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import common as cm
import figstyle as fs

fs.setup()

DATA = Path(__file__).resolve().parent / "data"
CACHE = DATA / "naive_k27me3_3class.npz"
FOCAL_TSV = DATA / "k27me3_focal_naive_3class.tsv"

MEM_COLOR, REFR_COLOR, NM_COLOR = "#b2182b", "#2a81ba", "#8a8a8a"
META_HW, META_NBINS = 10000, 200          # TSS +/- 10 kb in 100 bp bins
FOCAL_HW, FOCAL_NBINS = 3000, 1000        # promoter window, TSS +/- 3 kb
EPS = 1e-4
BW = "Naive_H3K27me3_merged_cpm.bw"
UNIVERSE = "protein-coding ISG universe"
# Random-genome null, 90th percentile focal CPM: the dashed reference line on Fig S9A.
# Carried over unchanged, not recomputed on this universe.
NULL_P90_CPM = 0.14494111715232333

CLASSES = [("memory", "memory", MEM_COLOR),
           ("refractory", "refractory", REFR_COLOR),
           ("non_memory", "non-mem", NM_COLOR)]

u = pd.read_csv(DATA / "universe.tsv", sep="\t")
cls = dict(zip(u.gene.astype(str), u.memory_class))
genes = {k: [g for g in u.gene.astype(str) if cls[g] == k] for k, _, _ in CLASSES}
print("%s: %d genes | " % (UNIVERSE, len(cls))
      + ", ".join("%s %d" % (lab, len(genes[k])) for k, lab, _ in CLASSES))

co = pd.read_csv(DATA / "coordinates.tsv", sep="\t")
coords = {str(r.gene): (r.chrom, int(r.tss), r.strand) for r in co.itertuples(index=False)}
print("coordinates.tsv: %d genes placed" % len(coords))

BW_PATH = cm.bigwig(BW)
_bw = pybigtools.open(str(BW_PATH))
try:
    CHROMS = set(_bw.chroms())
finally:
    _bw.close()


def metagene(chrom, tss, strand):
    """Strand-oriented Naive profile over TSS +/- META_HW. None if any bin is non-finite."""
    v = np.asarray(cm.read_signal(BW_PATH, chrom, tss - META_HW, tss + META_HW, META_NBINS), float)
    if not np.all(np.isfinite(v)):
        return None
    return v[::-1] if strand == "-" else v


def focal(chrom, tss):
    return float(np.nanmean(cm.read_signal(BW_PATH, chrom, tss - FOCAL_HW, tss + FOCAL_HW,
                                           FOCAL_NBINS)))


prof, kept, fvals, fgenes = {}, {}, {}, {}
cache_ok = False
if CACHE.exists():
    z = np.load(CACHE, allow_pickle=False)
    cache_ok = all("universe_%s" % k in z.files and list(z["universe_%s" % k]) == genes[k]
                   for k, _, _ in CLASSES)
    if cache_ok:
        for k, _, _ in CLASSES:
            prof[k], kept[k] = z["prof_%s" % k], list(z["kept_%s" % k])
            fvals[k], fgenes[k] = z["focal_%s" % k], list(z["focal_genes_%s" % k])
        print("loaded %s" % CACHE.name)
    else:
        print("%s was built on different gene sets; recomputing" % CACHE.name)

if not cache_ok:
    for k, lab, _ in CLASSES:
        mat, keep, fv, fg, dropped = [], [], [], [], []
        for g in genes[k]:
            c = coords.get(g)
            if c is None or c[0] not in CHROMS:
                dropped.append(g)
                continue
            v = metagene(*c)
            if v is None:
                dropped.append(g)
            else:
                mat.append(v)
                keep.append(g)
            f = focal(c[0], c[1])
            if np.isfinite(f):
                fv.append(f)
                fg.append(g)
        prof[k], kept[k] = np.array(mat), keep
        fvals[k], fgenes[k] = np.array(fv), fg
        print("  %-11s metagene %d, focal %d%s"
              % (lab, len(keep), len(fg), (", dropped %s" % dropped) if dropped else ""))
    np.savez_compressed(
        CACHE,
        **{"prof_%s" % k: prof[k] for k, _, _ in CLASSES},
        **{"kept_%s" % k: np.array(kept[k]) for k, _, _ in CLASSES},
        **{"focal_%s" % k: fvals[k] for k, _, _ in CLASSES},
        **{"focal_genes_%s" % k: np.array(fgenes[k]) for k, _, _ in CLASSES},
        **{"universe_%s" % k: np.array(genes[k]) for k, _, _ in CLASSES})
    print("cached per-gene profiles and focal CPM -> %s" % CACHE.name)

pd.concat([pd.DataFrame({"gene": fgenes[k], "memory_class": k, "naive": fvals[k]})
           for k, _, _ in CLASSES]).to_csv(FOCAL_TSV, sep="\t", index=False)

# ================================================================= Fig S9A: promoter level
kw = sps.kruskal(*[fvals[k] for k, _, _ in CLASSES])
pairs = {}
for a, b in (("memory", "refractory"), ("memory", "non_memory"), ("refractory", "non_memory")):
    pairs["%s vs %s" % (a, b)] = float(
        sps.mannwhitneyu(fvals[a], fvals[b], alternative="two-sided").pvalue)

fig, ax = plt.subplots(figsize=(2.5, 2.4))
bp = ax.boxplot([np.log10(fvals[k] + EPS) for k, _, _ in CLASSES],
                positions=range(len(CLASSES)), widths=0.5, showfliers=False, patch_artist=True)
for patch, (_, _, c) in zip(bp["boxes"], CLASSES):
    patch.set_facecolor(c)
    patch.set_alpha(0.45)
    patch.set_edgecolor("#333333")
for med in bp["medians"]:
    med.set_color("k")
for w in bp["whiskers"] + bp["caps"]:
    w.set_color("#333333")
rng = random.Random(7)
for i, (k, _, c) in enumerate(CLASSES):
    xs = [i + (rng.random() - 0.5) * 0.28 for _ in fvals[k]]
    ax.scatter(xs, np.log10(fvals[k] + EPS), s=4, color=c, alpha=0.45, edgecolors="none")
ax.axhline(np.log10(NULL_P90_CPM + EPS), color="#333333", lw=0.8, ls="--")
ax.set_xticks(range(len(CLASSES)))
ax.set_xticklabels(["%s\n(n=%s)" % (lab, format(len(fvals[k]), ",")) for k, lab, _ in CLASSES],
                   fontsize=6)
ax.set_ylabel("Avg. counts in promoter\nlog10 (CPM)")
ax.set_title("Naive H3K27me3 level", fontsize=8)
ax.spines[["top", "right"]].set_visible(False)
fig.tight_layout()
OUT_S9A = cfg.out_dir("FigS9A")
fs.save(fig, OUT_S9A / "FigS9A")

print("\n--- Fig S9A: Naive H3K27me3 level (focal, TSS +/- %d kb)" % (FOCAL_HW // 1000))
for k, lab, _ in CLASSES:
    print("  %-11s n=%4d  median %.4f  mean %.4f"
          % (lab, len(fvals[k]), np.median(fvals[k]), fvals[k].mean()))
print("  Kruskal-Wallis H=%.1f p=%.3g" % (kw.statistic, kw.pvalue))
for k, p in pairs.items():
    print("  Mann-Whitney %-24s p=%.3g" % (k, p))

(OUT_S9A / "FigS9A.json").write_text(json.dumps({
    "universe": UNIVERSE,
    "metric": "Naive focal H3K27me3 CPM, TSS +/- %d bp" % FOCAL_HW,
    "classes": "memory, refractory, non-memory",
    "n": {lab: int(len(fvals[k])) for k, lab, _ in CLASSES},
    "median": {lab: float(np.median(fvals[k])) for k, lab, _ in CLASSES},
    "mean": {lab: float(fvals[k].mean()) for k, lab, _ in CLASSES},
    "kruskal_wallis_p": float(kw.pvalue),
    "mann_whitney_pairwise_p": pairs,
    "dashed_line_null_p90_cpm": NULL_P90_CPM,
}, indent=2), encoding="utf-8")

# ================================================================= Fig 5A: metagene
xkb = np.linspace(-META_HW / 1000, META_HW / 1000, META_NBINS)

# Both panels of Fig 5 sit side by side, so the PLOT BOXES must be the same width. Matching
# figure widths would not do it, because the two carry different left margins; the target is
# read from the Fig 5C bar geometry rather than hardcoded.
_sizes = cfg.out_dir("Fig5C") / "Fig5C_bars.json"
if _sizes.exists():
    AXES_W_MM = float(json.loads(_sizes.read_text(encoding="utf-8"))["axes_width_mm"])
    print("\naxes width matched to Fig5C_bars.json: %.2f mm" % AXES_W_MM)
else:
    AXES_W_MM = 55.8
    print("\nFig5C_bars.json not found; using fallback axes width %.1f mm "
          "(run make_Fig5C_bars.py first for an exact match)" % AXES_W_MM)
LEFT_MM, RIGHT_MM = 14.0, 2.5
FIG_W_MM = AXES_W_MM + LEFT_MM + RIGHT_MM
FIG_H_MM = 2.4 * 25.4

fig2, ax2 = plt.subplots(figsize=fs.panel_mm(FIG_W_MM, FIG_H_MM))
fig2.subplots_adjust(left=LEFT_MM / FIG_W_MM, right=1 - RIGHT_MM / FIG_W_MM)
peak = 0.0
for k, lab, c in CLASSES:
    m = np.nanmean(prof[k], axis=0)
    peak = max(peak, float(m.max()))
    ax2.plot(xkb, m, color=c, lw=1.3, label="%s (n=%d)" % (lab, prof[k].shape[0]), zorder=3)
ax2.axvline(0, color="#666666", lw=0.6, ls=":", zorder=1)
ax2.set_xlim(xkb[0], xkb[-1])
Y_TOP = 0.175
if peak > Y_TOP:
    raise SystemExit("Y_TOP %.3f clips the data (peak %.4f); raise it" % (Y_TOP, peak))
ax2.set_ylim(0, Y_TOP)
ax2.set_xticks([-10, -5, 0, 5, 10])
ax2.set_xlabel("distance from TSS (kb)")
ax2.set_ylabel("Naive H3K27me3 (mean CPM)")
ax2.legend(fontsize=8, frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22),
           ncol=3, handlelength=1.4, columnspacing=1.2, borderpad=0.1)
ax2.set_title(UNIVERSE, fontsize=8, color="#555", pad=3)
ax2.text(0.985, 0.80, "line, mean\nno spread shown", transform=ax2.transAxes,
         ha="right", va="top", fontsize=8, color="#555", linespacing=1.25)
ax2.spines[["top", "right"]].set_visible(False)
OUT_5A = cfg.out_dir("Fig5A")
fs.save(fig2, OUT_5A / "Fig5A")

focal_mask = np.abs(xkb) <= FOCAL_HW / 1000
print("\n--- Fig 5A: Naive H3K27me3 metagene (mean line)")
meta_stats = {}
for k, lab, _ in CLASSES:
    m = np.nanmean(prof[k], axis=0)
    meta_stats[lab] = {"n": int(prof[k].shape[0]),
                       "focal_mean_CPM": float(np.mean(m[focal_mask])),
                       "peak_mean_CPM": float(m.max())}
    print("  %-11s n=%4d  focal mean %.4f  peak %.4f"
          % (lab, prof[k].shape[0], meta_stats[lab]["focal_mean_CPM"],
             meta_stats[lab]["peak_mean_CPM"]))
(OUT_5A / "Fig5A.json").write_text(json.dumps({
    "universe": UNIVERSE, "band": "none; mean line only",
    "condition": "Naive", "mark": "H3K27me3",
    "halfwidth_bp": META_HW, "nbins": META_NBINS, "bin_bp": 2 * META_HW // META_NBINS,
    "focal_halfwidth_bp": FOCAL_HW,
    "classes": "memory, refractory, non-memory",
    "stats": meta_stats}, indent=2), encoding="utf-8")

print("\nFig5A and FigS9A written")
