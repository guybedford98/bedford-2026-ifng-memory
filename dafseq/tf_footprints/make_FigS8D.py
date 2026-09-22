#!/usr/bin/env python3
"""Fig S8D: motif-annotated binding map for the CXCL11 promoter-proximal element,
TSS +225 to +425 bp (element centre +325 +/- 100 bp).

The scan reuses make_FigS8B.py (imported, not copied) so the motif calls, occupancy
statistics, expression annotation and candidate ranking are computed exactly as they were at
the IDO1 TSS. Two local substitutions: pysam is stubbed out, the reference is read through the
.fai index below, and the spanning-fiber sets come from the pure-Python BAM reader in
make_FigS8C.py. Scanning the TSS +/-100 bp window this way reproduces the published
cxcl11_tss_pm100_motif_hits.tsv exactly (853/853 hits, identical motif IDs, coordinates and
strands), and the fiber counts match naive_eligible_fibers in the published ranking table.

Coordinate note: hg38 is read 0-based half-open, matching pysam's fetch.

Strand note: CXCL11 is on the MINUS strand. x is genomic bp relative to the TSS, so POSITIVE
x = genomic downstream = biologically upstream (the promoter).
"""

from __future__ import annotations

import importlib.util
import sys
import types
from collections import OrderedDict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "shared"))

import config as cfg
import figstyle as fs

HG38 = cfg.path("reference", "hg38_fasta")
JASPAR = cfg.path("reference", "jaspar_motifs")
TPM_XLSX = cfg.path("rnaseq", "tpm_xlsx")
COUNTS_XLSX = cfg.path("rnaseq", "counts")

LOCUS = "CXCL11"
CHROM = "chr4"
TSS = 76_036_070
INDUCIBLE_LO, INDUCIBLE_HI = 290, 360
CENTRE = 325
CONSTITUTIVE_PEAK = 460
PWM_THRESHOLD = 0.85

# Curated representatives, resolved against the ranking table each run writes.
# Taking the raw top-N by rank gives near-identical STAT calls tiled across the
# window (the family bonus plus overlapping PWMs), which hides everything else -
# the same reason the TSS panel is curated.  One representative per distinct
# footprint instead, chosen from the ranking and the family table.
# Each entry is (tf_name, relative_start, note).
WINDOWS = [
    {
        "lo": CENTRE - 100,
        "hi": CENTRE + 100,
        "prefix": "cxcl11_centre_pm100",
        "panel": "FigS8D",
        "title": f"{LOCUS} memory element (centre +{CENTRE} ± 100 bp): "
                 "FiberHMM protection density with top candidate motifs",
        "figsize": (5.8, 4.8),
        "height_ratios": [2.5, 2.6],
        "ymax": 0.092,
        "xticks": [225, 275, 325, 375, 425],
        "show_constitutive": False,
        "centre_axis": True,
        "curated": [
            ("STAT1", 231, "high-scoring GAS-like site, NOT inducible"),
            ("STAT3", 291, "inducible cluster"),
            ("ETV7", 294, "inducible cluster"),
            ("STAT1::STAT2", 305, "ISGF3 heterodimer, inducible cluster"),
            ("ETS2", 306, "inducible cluster"),
            ("Irf1", 307, "co-located with ISGF3 - ISRE-like composite"),
            ("KLF6", 327, "inducible cluster"),
            ("NR2C2", 343, "inducible cluster, 3' edge"),
            ("STAT3", 366, "+397 shoulder"),
            ("STAT1::STAT2", 393, "second ISGF3 site at the +397 peak"),
            ("IRF7", 393, "co-located with the second ISGF3 site"),
        ],
    },
]

FONT_LABEL = 8
FONT_NUMBER = 6


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module  # dataclass decorator needs this
    spec.loader.exec_module(module)
    return module


# pysam has no Windows build; nothing we call from `ident` touches it.
sys.modules.setdefault("pysam", types.ModuleType("pysam"))
ident = load_module("ident_cxcl11", HERE / "make_FigS8B.py")
panels = load_module("cxcl11_panels", HERE / "make_FigS8C.py")

# Point the reused pipeline functions at CXCL11.
ident.TSS = TSS
ident.LOCUS_LABEL = LOCUS
ident.SAMPLES = panels.SAMPLES

COLORS = panels.COLORS
FAMILY_COLORS = ident.FAMILY_COLORS


def fetch_reference(chrom: str, start: int, end: int) -> str:
    """0-based half-open sequence fetch via the .fai index (pysam fetch semantics)."""
    fai = HG38.parent / (HG38.name + ".fai")
    for line in fai.read_text().splitlines():
        fields = line.split("\t")
        if fields[0] == chrom:
            _length, offset, line_bases, line_width = (int(x) for x in fields[1:5])
            break
    else:
        raise KeyError(f"{chrom} not in {fai}")
    span = end - start
    with HG38.open("rb") as handle:
        handle.seek(offset + start // line_bases * line_width + start % line_bases)
        raw = handle.read(span + span // line_bases + 2)
    return raw.decode().replace("\n", "").replace("\r", "")[:span].upper()


def draw_motif_track(ax, motif_rows: list[dict], lo: int, hi: int) -> None:
    """One motif per row, box at its true coordinates, horizontal label beside it."""
    n = len(motif_rows)
    ax.set_xlim(lo, hi)
    ax.set_ylim(0, n)
    ax.set_yticks([])
    ax.spines[["left", "right", "top", "bottom"]].set_visible(False)
    ax.tick_params(axis="x", length=0)
    bp_per_char = (hi - lo) / 62.0
    for i, row in enumerate(motif_rows):
        y = n - i - 0.5
        x0, x1 = row["relative_start"], row["relative_end"]
        ax.add_patch(
            plt.Rectangle(
                (x0, y - 0.26),
                max(x1 - x0, 1),
                0.52,
                facecolor=FAMILY_COLORS.get(row["family"], FAMILY_COLORS["Other"]),
                edgecolor="none",
                alpha=0.85,
            )
        )
        label = row["tf_name"] if row["family"] == "Other" else f"{row['family']} {row['tf_name']}"
        pad = 0.012 * (hi - lo)
        if hi - x1 - pad >= len(label) * bp_per_char:
            ax.text(x1 + pad, y, label, fontsize=FONT_NUMBER, ha="left",
                    va="center", color="#222222")
        else:
            ax.text(x0 - pad, y, label, fontsize=FONT_NUMBER, ha="right",
                    va="center", color="#222222")


def build_window(window: dict, motifs, tpm, counts) -> None:
    lo, hi = window["lo"], window["hi"]
    window_start, window_end = TSS + lo, TSS + hi
    out_dir = cfg.out_dir(window["panel"])
    prefix = window["prefix"]
    ident.OUTPUT_PREFIX = prefix
    print(f"\n=== {CHROM}:{window_start:,}-{window_end:,} ({lo:+d} to {hi:+d} bp) ===")
    out_dir.mkdir(parents=True, exist_ok=True)

    # ---- motif scan (reused verbatim from make_FigS8B.py) ----
    seq = fetch_reference(CHROM, window_start, window_end)
    hits = ident.scan_motifs(motifs, seq, window_start, PWM_THRESHOLD)
    for hit in hits:
        hit["relative_start"] = hit["abs_start"] - TSS
        hit["relative_end"] = hit["abs_end"] - TSS
        hit["relative_center"] = ((hit["abs_start"] + hit["abs_end"]) // 2) - TSS
    print(f"  motif hits at PWM score fraction >= {PWM_THRESHOLD}: {len(hits)}")

    # ---- expression support + occupancy + ranking ----
    expr_by_tf = ident.annotate_hits_with_expression(hits, tpm, counts)
    condition_data = panels.load_condition_data(CHROM, window_start, window_end)
    occupancy_rows, summary_rows = ident.compute_occupancy(hits, condition_data)
    candidate_rows = ident.add_candidate_scores(summary_rows, expr_by_tf)
    family_rows = ident.collapse_families(candidate_rows)

    ident.save_tsv(out_dir / f"{prefix}_motif_hits.tsv", hits, ident.motif_hit_fields())
    ident.save_tsv(out_dir / f"{prefix}_motif_occupancy_by_condition.tsv", occupancy_rows)
    ident.save_tsv(out_dir / f"{prefix}_candidate_tf_ranking.tsv", candidate_rows,
                   ident.candidate_fields())
    ident.save_tsv(out_dir / f"{prefix}_candidate_families.tsv", family_rows)
    print(f"  wrote 4 tables to {out_dir}")

    # ---- density (identical computation to Fig S8B and Fig S8C) ----
    density = panels.per_base_density(condition_data, window_start, window_end)
    x = np.arange(window_start, window_end) - TSS
    smoothed = {c: gaussian_filter1d(density[c], sigma=2) for c in panels.SAMPLES}

    # ---- per-position TSV so the panel never has to be reverse-engineered ----
    rows = []
    for condition, sample in panels.SAMPLES.items():
        n_eligible = condition_data[condition]["n_eligible"]
        raw = density[condition]
        for i, pos in enumerate(range(window_start, window_end)):
            rows.append({
                "locus": LOCUS, "chrom": CHROM, "position": pos,
                "position_relative_to_tss": pos - TSS,
                "position_relative_to_element_centre": pos - TSS - CENTRE,
                "condition": condition, "sample": sample,
                "eligible_fibers": n_eligible,
                "fibers_with_tf_call": int(round(raw[i] * n_eligible)),
                "tf_call_fraction": f"{raw[i]:.6f}",
                "tf_call_fraction_smoothed_sigma2": f"{smoothed[condition][i]:.6f}",
            })
    panels.write_tsv(out_dir / f"{window['panel']}_per_position.tsv", rows)
    print(f"  wrote {window['panel']}_per_position.tsv ({len(rows)} rows)")

    by_key = {}
    for row in candidate_rows:
        by_key.setdefault((row["tf_name"], row["relative_start"]), row)
    plot_rows = []
    for tf_name, rel_start, note in window["curated"]:
        row = by_key.get((tf_name, rel_start))
        if row is None:
            raise ValueError(f"no ranking row for {tf_name} at {rel_start:+d}")
        plot_rows.append({**row, "curation_note": note})
    plot_rows.sort(key=lambda r: r["relative_start"])

    fig, (ax_track, ax) = plt.subplots(
        2, 1, figsize=window["figsize"], sharex=True,
        gridspec_kw={"height_ratios": window["height_ratios"], "hspace": 0.08},
    )
    draw_motif_track(ax_track, plot_rows, lo, hi)
    ax_track.axvspan(INDUCIBLE_LO, INDUCIBLE_HI, color="#0B519D", alpha=0.07, lw=0)
    ax_track.set_ylabel("Top\nmotifs", fontsize=FONT_LABEL, rotation=0, labelpad=20, va="center")
    ax_track.set_title(window["title"], fontsize=FONT_LABEL, pad=6)

    ymax = window["ymax"]
    ax.axvspan(INDUCIBLE_LO, INDUCIBLE_HI, color="#0B519D", alpha=0.07, lw=0)
    for condition in panels.SAMPLES:
        ax.plot(x, smoothed[condition], color=COLORS[condition], lw=1.1, label=condition)
    ax.axvline(CENTRE, color="#222222", lw=0.7, ls="--")
    ax.annotate(f"centre {CENTRE:+d}", xy=(CENTRE - 5, ymax * 0.82), fontsize=FONT_NUMBER,
                ha="right", va="top", color="#222222")
    ax.annotate(f"inducible cluster {INDUCIBLE_LO:+d} to {INDUCIBLE_HI:+d}",
                xy=(CENTRE + 6, ymax * 0.96), fontsize=FONT_NUMBER, ha="left", va="top",
                color="#0B519D")
    if window["show_constitutive"]:
        ax.annotate(f"constitutive peak {CONSTITUTIVE_PEAK:+d}",
                    xy=(CONSTITUTIVE_PEAK + 6, ymax * 0.89), fontsize=FONT_NUMBER,
                    ha="left", va="center", color="#555555")
    ax.set_ylim(0, ymax)
    ax.set_xlim(lo, hi)
    if window["xticks"]:
        ax.set_xticks(window["xticks"])
    ax.set_xlabel(f"Position relative to {LOCUS} TSS (bp)")
    ax.set_ylabel("Fraction of spanning fibers\nwith FiberHMM TF call")
    ax.legend(frameon=False, ncols=1, loc="upper left", handlelength=1.4,
              columnspacing=1.1, borderpad=0, labelspacing=0.3)
    ax.spines[["top", "right"]].set_visible(False)

    note_y = -0.28
    if window["centre_axis"]:
        # Second axis in element-centred coordinates.  The primary axis stays in bp
        # relative to the TSS - the convention every other CXCL11 output uses - so the
        # "+/-100 of the centre" framing is added, never substituted.
        sec = ax.secondary_xaxis(-0.26, functions=(lambda v: v - CENTRE, lambda v: v + CENTRE))
        sec.set_xlabel(f"bp from element centre (+{CENTRE})")
        sec.set_xticks([-100, -50, 0, 50, 100])
        sec.spines["bottom"].set_linewidth(0.6)
        sec.tick_params(width=0.6, length=2.4, labelsize=FONT_NUMBER)
        note_y = -0.62

    ax.annotate(panels.STRAND_NOTE, xy=(0.5, note_y), xycoords="axes fraction",
                ha="center", va="top", fontsize=FONT_NUMBER, color="#555555")

    fs.save(fig, out_dir / window["panel"])

    # ---- what is actually in the inducible cluster ----
    print(f"\n  curated candidates plotted (window {lo:+d}..{hi:+d})")
    for row in plot_rows:
        print(f"    rank {row['rank']:4d}  {row['tf_name']:14s} {row['family']:11s}"
              f" {row['relative_start']:+4d}..{row['relative_end']:+4d}"
              f"  score {row['score_frac']:.3f}"
              f"  stim-naive {row['stim_vs_naive_delta']:+.3f}"
              f"  rechal-rested {row['rechallenge_vs_rested_delta']:+.3f}"
              f"  {row['expression_tier']:22s} {row['curation_note']}")

    print(f"\n  top 12 candidates INSIDE the inducible cluster"
          f" ({INDUCIBLE_LO:+d}..{INDUCIBLE_HI:+d})")
    inside = [r for r in candidate_rows
              if r["relative_center"] >= INDUCIBLE_LO and r["relative_center"] <= INDUCIBLE_HI]
    for row in inside[:12]:
        print(f"    rank {row['rank']:4d}  {row['tf_name']:14s} {row['family']:11s}"
              f" {row['relative_start']:+4d}..{row['relative_end']:+4d}"
              f"  score {row['score_frac']:.3f}"
              f"  stim-naive {row['stim_vs_naive_delta']:+.3f}"
              f"  rechal-rested {row['rechallenge_vs_rested_delta']:+.3f}"
              f"  q(rechal) {row['rechallenge_vs_rested_q']:.2g}  {row['expression_tier']}")

    print("\n  families represented in the inducible cluster")
    fam = OrderedDict()
    for row in inside:
        fam.setdefault(row["family"], []).append(row)
    for family, rows in sorted(fam.items(), key=lambda kv: kv[1][0]["rank"]):
        best = rows[0]
        print(f"    {family:11s} n={len(rows):3d}  best rank {best['rank']:4d}"
              f" {best['tf_name']:14s} at {best['relative_start']:+d}..{best['relative_end']:+d}"
              f"  {best['expression_tier']}")

    if window["show_constitutive"]:
        print(f"\n  candidates at the constitutive peak ({CONSTITUTIVE_PEAK:+d} +/- 15)")
        at_const = [r for r in candidate_rows
                    if abs(r["relative_center"] - CONSTITUTIVE_PEAK) <= 15][:6]
        for row in at_const:
            print(f"    rank {row['rank']:4d}  {row['tf_name']:14s} {row['family']:11s}"
                  f" {row['relative_start']:+4d}..{row['relative_end']:+4d}"
                  f"  stim-naive {row['stim_vs_naive_delta']:+.3f}"
                  f"  rechal-rested {row['rechallenge_vs_rested_delta']:+.3f}")


def main() -> None:
    fs.setup()
    motifs = ident.parse_meme_motifs(JASPAR)
    tpm = ident.load_tpm_expression(TPM_XLSX)
    counts = ident.load_count_expression(COUNTS_XLSX)
    for window in WINDOWS:
        build_window(window, motifs, tpm, counts)


if __name__ == "__main__":
    main()
