#!/usr/bin/env python3
"""
make_Fig3C.py  --  zoomed-out, two-condition H3K27me3 view over CXCL11 showing that the
loss is FOCAL at the gene and absent at the DISTAL Polycomb domain.

Two conditions only (Naive vs Rested) keep the wide view readable. The figure has a
signal overlay (Naive + trained, shared y), a difference track (Naive minus trained)
that peaks at the gene and is flat at the domain, and a packed gene track. It also
prints a quantification (mean H3K27me3 at the gene vs at the domain, per condition) so
the focal-vs-distal claim is backed by numbers, not just the picture.

Run:  py make_Fig3C.py
"""
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.gridspec import GridSpec

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
import config as cfg  # noqa: E402
import common as cm  # noqa: E402
import tracklib as tl  # noqa: E402

MODELS = json.loads((ROOT / "config" / "gene_models.json").read_text(encoding="utf-8"))
PROBE = json.loads((ROOT / "config" / "probe_summary.json").read_text(encoding="utf-8"))
CCOLOR = {"Naive": "#90959a", "Stim": "#4f9bcf", "Rested": "#4b4b4c", "Rechallenged": "#1e5d9f"}
CLABEL = {"Naive": "Naive", "Stim": "Stimulated", "Rested": "Rested", "Rechallenged": "Rechallenged"}
LOSS = "#c1272d"
LOSS_FILL = "#dd888c"   # opaque equivalent of LOSS at 55% over white (no transparency to flatten)
GAIN_FILL = "#b8e2ca"   # opaque equivalent of the gain green
SMALL, TIT = 6, 8
FOCAL_HALF = 4000   # focal window half width around TSS
NEIGHBORS = {}      # {gene: [neighbor gene models]} loaded in main()

PANEL = "Fig3C"
GENE = "CXCL11"
TRAINED = "Rested"
ZOOM_HALFWIDTH = 1000000


def bwp(c):
    return cm.bigwig(f"{c}_H3K27me3_merged_cpm.bw")


def load(gene, s, e, smooth_bp):
    chrom = MODELS["genes"][gene]["chrom"]
    span = e - s
    nb = max(200, int(span / max(100, span / 3000)))
    binbp = span / nb
    edges = np.linspace(s, e, nb + 1)
    x = 0.5 * (edges[:-1] + edges[1:])
    prof = {}
    for c in CCOLOR:
        v = tl.read_signal(bwp(c), chrom, s, e, nb)
        k = max(1, int(round(smooth_bp / binbp)))
        if k > 1:
            v = np.convolve(v, np.ones(k) / k, mode="same")
        prof[c] = v
    return x, prof


def mean_in(prof, x, s, e):
    msk = (x >= s) & (x <= e)
    return {c: float(np.nanmean(prof[c][msk])) for c in prof}


def quantify(gene):
    gm = MODELS["genes"][gene]
    tss = gm["tss"]
    dom = PROBE["genes"][gene]["domain"]
    s = min(tss - FOCAL_HALF, dom["start"]) - 5000
    e = max(tss + FOCAL_HALF, dom["end"]) + 5000
    x, prof = load(gene, s, e, 1000)
    focal = mean_in(prof, x, tss - FOCAL_HALF, tss + FOCAL_HALF)
    distal = mean_in(prof, x, dom["start"], dom["end"])
    out = {"gene": gene, "focal_window": [tss - FOCAL_HALF, tss + FOCAL_HALF],
           "domain_window": [dom["start"], dom["end"]],
           "focal_means": focal, "distal_means": distal, "contrasts": {}}
    print(f"\n=== {gene} focal (gene TSS +/- {FOCAL_HALF//1000}kb) vs distal (Polycomb domain) ===")
    print(f"  focal  means: " + "  ".join(f"{c}={focal[c]:.4f}" for c in CCOLOR))
    print(f"  distal means: " + "  ".join(f"{c}={distal[c]:.4f}" for c in CCOLOR))
    for trained in ("Rested", "Rechallenged"):
        f_drop = focal["Naive"] - focal[trained]
        d_drop = distal["Naive"] - distal[trained]
        f_rel = f_drop / (focal["Naive"] + 1e-9)
        d_rel = d_drop / (distal["Naive"] + 1e-9)
        out["contrasts"][trained] = {"focal_drop": f_drop, "focal_rel": f_rel,
                                     "distal_drop": d_drop, "distal_rel": d_rel,
                                     "specificity_ratio": f_rel / (abs(d_rel) + 1e-9)}
        print(f"  Naive vs {trained:12s}: focal drop {f_rel*100:5.0f}%  |  distal drop {d_rel*100:5.0f}%"
              f"  |  focal/distal specificity {f_rel/(abs(d_rel)+1e-9):5.1f}x")
    return out


def _scalebar_bp(span):
    cand = [5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000]
    bar = min(cand, key=lambda c: abs(c - span / 6))
    lbl = f"{bar//1000} kb" if bar < 1_000_000 else f"{bar//1_000_000} Mb"
    return bar, lbl


def _parse_field(field):
    if isinstance(field, list):
        return [int(x) for x in field if str(x).strip() != ""]
    return [int(x) for x in str(field).split(",") if x.strip() != ""]


def fetch_refgene_window(chrom, lo, hi):
    url = (f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=refGene;"
           f"chrom={chrom};start={int(lo)};end={int(hi)}")
    with urllib.request.urlopen(url, timeout=25) as r:
        recs = json.loads(r.read()).get("refGene", [])
    bysym = {}
    for g in recs:
        sym = g.get("name2") or g.get("name")
        if not sym:
            continue
        span = int(g["txEnd"]) - int(g["txStart"])
        if sym not in bysym or span > bysym[sym]["_span"]:
            bysym[sym] = {"name": sym, "chrom": chrom, "tx_start": int(g["txStart"]),
                          "tx_end": int(g["txEnd"]), "strand": g.get("strand", "+"),
                          "cds_start": int(g.get("cdsStart", g["txStart"])),
                          "cds_end": int(g.get("cdsEnd", g["txEnd"])),
                          "exon_starts": _parse_field(g.get("exonStarts", [])),
                          "exon_ends": _parse_field(g.get("exonEnds", [])), "_span": span}
    genes = list(bysym.values())
    for g in genes:
        g.pop("_span", None)
    return genes


def load_neighbor_genes():
    """Genes within +/- 1.2 Mb of each focal gene, from UCSC refGene, cached to config/."""
    cache = ROOT / "config" / "neighbor_genes.json"
    if cache.exists():
        try:
            return json.loads(cache.read_text(encoding="utf-8"))
        except Exception:
            pass
    out = {}
    for gene, gm in MODELS["genes"].items():
        center = (gm["tx_start"] + gm["tx_end"]) // 2
        try:
            out[gene] = fetch_refgene_window(gm["chrom"], center - 1_200_000, center + 1_200_000)
            print(f"  fetched {len(out[gene])} genes near {gene}")
        except Exception as exc:
            print(f"  [neighbor fetch failed for {gene}: {exc!r}] focal gene only")
            out[gene] = [{"name": gene, "chrom": gm["chrom"], "tx_start": gm["tx_start"],
                          "tx_end": gm["tx_end"], "strand": gm["strand"],
                          "cds_start": gm.get("cds_start", gm["tx_start"]),
                          "cds_end": gm.get("cds_end", gm["tx_end"]),
                          "exon_starts": gm.get("exon_starts", [gm["tx_start"]]),
                          "exon_ends": gm.get("exon_ends", [gm["tx_end"]])}]
    cache.write_text(json.dumps(out, indent=2), encoding="utf-8")
    return out


def pack_genes(genes, s, e):
    """Pack genes into non-overlapping rows (genome-browser style). Each gene's footprint
    is the wider of its transcript extent and its label, so neither glyphs nor names
    overlap. Every gene is labeled. Returns (list of (gene, row, center), n_rows)."""
    span = e - s
    vis = sorted((g for g in genes if g["tx_end"] >= s and g["tx_start"] <= e),
                 key=lambda g: g["tx_start"])
    pad = span * 0.008
    rows_end = []
    assign = []
    for g in vis:
        cen = (g["tx_start"] + g["tx_end"]) / 2
        lh = len(g["name"]) * span * 0.0046   # approx half label width in bp
        d0 = min(g["tx_start"], cen - lh)
        d1 = max(g["tx_end"], cen + lh)
        row = None
        for ri, rend in enumerate(rows_end):
            if d0 > rend + pad:
                row = ri
                rows_end[ri] = d1
                break
        if row is None:
            rows_end.append(d1)
            row = len(rows_end) - 1
        assign.append((g, row, cen))
    return assign, max(1, len(rows_end))


def _gene_track(ax, assign, nrows, s, e, focal_name, simple=False):
    """Draw packed genes with the italic name below each. Default: exon boxes (tall CDS,
    short UTR), intron line, and strand chevrons. simple=True: just a single line spanning
    each gene (no exon/intron detail), for narrow panels where the structure is too small
    to read. Focal gene bold. No TSS arrow."""
    span = e - s
    ax.set_xlim(s, e)
    ax.set_ylim(-(nrows - 1) - 0.62, 0.34)
    ax.axis("off")
    for g, ri, cen in assign:
        y = -float(ri)
        if simple:
            ax.plot([g["tx_start"], g["tx_end"]], [y, y], color="#333333", lw=1.6,
                    solid_capstyle="round", zorder=2)
        else:
            cds0 = g.get("cds_start", g["tx_start"])
            cds1 = g.get("cds_end", g["tx_end"])
            ax.plot([g["tx_start"], g["tx_end"]], [y, y], color="#333333", lw=0.9,
                    solid_capstyle="round", zorder=2)
            for a, b in zip(g.get("exon_starts", [g["tx_start"]]), g.get("exon_ends", [g["tx_end"]])):
                for lo, hi, h in ((a, min(b, cds0), 0.13), (max(a, cds1), b, 0.13),
                                  (max(a, cds0), min(b, cds1), 0.26)):
                    if hi > lo:
                        ax.add_patch(mpatches.Rectangle((lo, y - h / 2), hi - lo, h,
                                     facecolor="#333333", edgecolor="none", zorder=3))
            marker = ">" if g["strand"] == "+" else "<"
            gs_, ge_ = max(g["tx_start"], s), min(g["tx_end"], e)
            if ge_ > gs_:
                for xp in np.linspace(gs_, ge_, 4)[1:-1]:
                    ax.plot(xp, y, marker=marker, markersize=2.2, color="#999999",
                            markeredgewidth=0, zorder=4)
        cen_c = min(max(cen, s + span * 0.01), e - span * 0.01)
        ax.text(cen_c, y - 0.32, g["name"], ha="center", va="top", fontsize=SMALL - 0.5,
                fontstyle="italic", color="#222222",
                fontweight=("bold" if g["name"] == focal_name else "normal"))


def zoom_win(gene, hw):
    gm = MODELS["genes"][gene]
    return (gm["chrom"], gm["tss"] - hw, gm["tss"] + hw)


def diff_dmax_at(gene, trained, hw):
    """Max abs(Naive minus trained) over a TSS +/- hw window, for fixing a shared axis."""
    gm = MODELS["genes"][gene]
    s, e = gm["tss"] - hw, gm["tss"] + hw
    smooth_bp = min(8000, max(2500, int((e - s) / 150)))
    _, prof = load(gene, s, e, smooth_bp)
    return float(np.nanmax(np.abs(prof["Naive"] - prof[trained])))


def combo(gene, trained, stem, window=None, diff_ymax=None, fig_width=7.2, simple_genes=False):
    chrom, s, e = window
    span = e - s
    q = quant_lookup[gene]["contrasts"][trained]
    smooth_bp = min(8000, max(2500, int(span / 150)))
    x, prof = load(gene, s, e, smooth_bp)
    naive, tr = prof["Naive"], prof[trained]
    diff = naive - tr

    assign, nrows = pack_genes(NEIGHBORS.get(gene, []), s, e)
    fig = plt.figure(figsize=(fig_width, 3.15 + 0.30 * nrows))
    gs = GridSpec(4, 1, height_ratios=[1.6, 2.4, 0.55 * nrows + 0.15, 0.5],
                  left=0.10, right=0.985, top=0.90, bottom=0.08, hspace=0.28)

    # row 1: two-condition overlay, solid lines only (no fill)
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(s, e)
    ymax = max(np.nanmax(naive), np.nanmax(tr))
    for c, v in ((trained, tr), ("Naive", naive)):
        ax1.plot(x, v, color=CCOLOR[c], lw=0.7, zorder=3, solid_capstyle="round")
    ax1.set_ylim(0, ymax * 1.1)
    ax1.set_yticks([0, ymax])
    ax1.set_yticklabels(["0", f"{ymax:.2f}"], fontsize=SMALL)
    for sp in ("top", "right"):
        ax1.spines[sp].set_visible(False)
    ax1.spines["bottom"].set_visible(True)
    ax1.spines["bottom"].set_linewidth(0.6)
    ax1.spines["bottom"].set_color("#cccccc")
    ax1.tick_params(axis="x", bottom=False, labelbottom=False)
    ax1.tick_params(axis="y", length=2.2, width=0.6, labelsize=SMALL)
    ax1.set_ylabel("H3K27me3\nCPM", fontsize=SMALL)
    lblbox = dict(facecolor="white", edgecolor="none", pad=0.5)
    ax1.text(0.995, 0.94, "Naive", transform=ax1.transAxes, ha="right", va="top",
             color=CCOLOR["Naive"], fontsize=SMALL + 0.5, fontweight="bold", bbox=lblbox)
    ax1.text(0.995, 0.76, CLABEL[trained], transform=ax1.transAxes, ha="right", va="top",
             color=CCOLOR[trained], fontsize=SMALL + 0.5, fontweight="bold", bbox=lblbox)

    # row 2: difference track (Naive minus trained)
    ax2 = fig.add_subplot(gs[1, 0])
    ax2.set_xlim(s, e)
    ax2.axhline(0, color="#999", lw=0.6, zorder=1)
    ax2.fill_between(x, 0, diff, where=(diff >= 0), color=LOSS_FILL, linewidth=0, zorder=2)
    ax2.fill_between(x, 0, diff, where=(diff < 0), color=GAIN_FILL, linewidth=0, zorder=2)
    ax2.plot(x, diff, color="#333", lw=0.5, zorder=3)
    dmax = float(np.nanmax(np.abs(diff)))
    tickmax = diff_ymax if diff_ymax is not None else dmax
    ax2.set_ylim(-tickmax * 1.1, tickmax * 1.1)
    ax2.set_yticks([round(-tickmax, 2), 0, round(tickmax, 2)])
    ax2.tick_params(axis="y", length=2.2, width=0.6, labelsize=SMALL)
    for sp in ("top", "right", "bottom"):
        ax2.spines[sp].set_visible(False)
    ax2.tick_params(axis="x", bottom=False, labelbottom=False)
    ax2.set_ylabel(f"Naive minus\n{CLABEL[trained]} (CPM)", fontsize=SMALL)
    ax2.text(0.995, 0.92, "H3K27me3 lost (Naive higher)", transform=ax2.transAxes,
             ha="right", va="top", color=LOSS, fontsize=SMALL, fontweight="bold")

    # row 3: multi-gene track (every gene labeled, packed into rows so names do not overlap)
    axg = fig.add_subplot(gs[2, 0])
    _gene_track(axg, assign, nrows, s, e, gene, simple=simple_genes)

    # row 4: coordinate axis + scale bar
    ax3 = fig.add_subplot(gs[3, 0])
    ax3.set_xlim(s, e)
    ax3.set_ylim(-1.0, 0.55)
    ax3.axis("off")
    ax3.plot([s, e], [0.32, 0.32], color="#999", lw=0.5, clip_on=False)
    for t in np.linspace(s, e, 7):
        ax3.plot([t, t], [0.24, 0.32], color="#666", lw=0.5, clip_on=False)
        ax3.text(t, 0.16, f"{t/1e6:.2f}", ha="center", va="top", fontsize=SMALL, color="#444")
    ax3.text(s, -0.42, f"{chrom} (Mb)", ha="left", va="top", fontsize=SMALL, color="#444")
    bar, blbl = _scalebar_bp(span)
    x1 = e - span * 0.02
    x0 = x1 - bar
    ax3.plot([x0, x1], [-0.42, -0.42], color="#222", lw=1.6, solid_capstyle="butt")
    ax3.text((x0 + x1) / 2, -0.58, blbl, ha="center", va="top", fontsize=SMALL, color="#222")

    matched = "   |   difference axis matched to CXCL11" if diff_ymax is not None else ""
    coord = f"{chrom}:{s:,}-{e:,}  ({span/1000:.0f} kb)  hg38  CPM"
    title = f"H3K27me3 is lost focally at {gene}, not across the distal Polycomb domain"
    sub = f"Naive vs {CLABEL[trained]}  |  {coord}   |   focal loss {q['focal_rel']*100:.0f}%{matched}"
    fig.text(0.01, 0.985, title, ha="left", va="top", fontsize=TIT, fontweight="bold")
    fig.text(0.01, 0.935, sub, ha="left", va="top", fontsize=SMALL, color="#555")

    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".svg"), bbox_inches="tight", pad_inches=0.02)
    fig.savefig(stem.with_suffix(".png"), dpi=200, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print(f"  wrote {stem.name}")


quant_lookup = {}


def main():
    tl.setup()
    quant_lookup[GENE] = quantify(GENE)
    (cfg.out_dir(PANEL) / "focal_distal_quant.json").write_text(
        json.dumps(quant_lookup, indent=2), encoding="utf-8")
    global NEIGHBORS
    NEIGHBORS = load_neighbor_genes()

    # the difference-track axis is fixed to CXCL11's own maximum over the 1 Mb window
    ref = diff_dmax_at(GENE, TRAINED, ZOOM_HALFWIDTH)
    print(f"  1 Mb difference axis (Naive vs {TRAINED}) fixed to {GENE}: +/-{ref*1.1:.3f} CPM")
    combo(GENE, TRAINED, cfg.out_dir(PANEL) / PANEL,
          window=zoom_win(GENE, ZOOM_HALFWIDTH), diff_ymax=ref)


if __name__ == "__main__":
    main()
