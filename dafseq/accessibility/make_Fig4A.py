#!/usr/bin/env python3
"""Fig 4A: continuous stacked single-fiber MSP rasters for IDO1, CXCL11 and STAT1.

One raster per locus, four conditions stacked with spacer rows, gene track below.
"""
import argparse
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib
import matplotlib.colors as mcolors
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import figstyle as fs
import msp_matrices as mm

ILLUSTRATOR_SAFE_RCPARAMS = {
    "pdf.compression": 0,
    "image.composite_image": False,
    "svg.image_inline": True,
    "savefig.transparent": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
}
CONDENSED_EXPORT_DPI = 600
MSP_GRAY = "#4A4A4B"
SPACER_ROWS = 28


def matrix_to_rgb(matrix: np.ndarray) -> np.ndarray:
    """Render binary MSP matrix as white background with neutral-gray calls."""
    color = np.array(mcolors.to_rgb(MSP_GRAY), dtype=float)
    rgb = np.ones((matrix.shape[0], matrix.shape[1], 3), dtype=float)
    if matrix.size:
        mask = matrix > 0
        rgb[mask] = color
    return rgb


def blank_rgb(rows: int, cols: int, color: str = "#f1f1f1") -> np.ndarray:
    arr = np.ones((rows, cols, 3), dtype=float)
    arr[:] = np.array(mcolors.to_rgb(color), dtype=float)
    return arr


def build_matrices(args: argparse.Namespace) -> Dict[str, List[object]]:
    configs = mm.load_locus_configs()
    sample_sets = mm.load_sample_configs()
    out: Dict[str, List[object]] = {}

    for locus in args.loci:
        base_cfg = configs[locus]
        locus_cfg = base_cfg.trimmed(args.trim_left_bp, args.trim_right_bp)
        matrices = []
        for sample_cfg in sample_sets[locus]:
            matrices.append(
                mm.build_sample_matrix(
                    locus_cfg,
                    sample_cfg,
                    "both",
                    args.bin_size,
                    args.min_coverage_frac,
                )
            )
        mm.apply_global_msp_state_order(
            matrices,
            max_fibers_per_condition=0,
            sort_window_bp=args.sort_window_bp,
            display_depth_mode="equal_subsample",
            random_state=42,
            cluster_mode="unclustered" if locus in args.no_cluster_loci else "clustered",
        )
        out[locus] = matrices
    return out


def combined_locus_rgb(mats: Sequence[object], spacer_rows: int = SPACER_ROWS) -> Tuple[np.ndarray, List[float], List[Tuple[str, int, int, str, str]]]:
    blocks = []
    boundaries = []
    labels = []
    row_cursor = 0
    n_cols = mats[0].selected_matrix.shape[1]
    for idx, sample in enumerate(mats):
        rgb = matrix_to_rgb(sample.selected_matrix)
        blocks.append(rgb)
        if sample.msp_state_mode == "clustered":
            value = (sample.high_state_count / sample.selected_reads * 100.0) if sample.selected_reads else 0.0
            labels.append((sample.condition, row_cursor, row_cursor + sample.selected_reads, f"{value:.1f}%", "rich"))
        else:
            value = sample.selected_mean * 100.0
            labels.append((sample.condition, row_cursor, row_cursor + sample.selected_reads, f"{value:.1f}%", "MSP-covered area"))
        row_cursor += sample.selected_reads
        if idx != len(mats) - 1:
            blocks.append(blank_rgb(spacer_rows, n_cols, "#f5f5f5"))
            row_cursor += spacer_rows
            boundaries.append(row_cursor - spacer_rows / 2)
    return np.vstack(blocks), boundaries, labels


def add_condition_rails(ax, labels: Sequence[Tuple[str, int, int, str, str]], total_rows: int) -> None:
    """Condition-color rails outside the neutral MSP raster."""
    for condition, y0, y1, _value, _metric_label in labels:
        y_bottom = 1.0 - y1 / total_rows
        height = (y1 - y0) / total_rows
        ax.add_patch(
            mpatches.Rectangle(
                (-0.045, y_bottom),
                0.022,
                height,
                transform=ax.transAxes,
                clip_on=False,
                facecolor=mm.CONDITION_COLORS.get(condition, "#2166ac"),
                edgecolor="none",
                zorder=6,
            )
        )


def plot_stacked_rasters(matrices_by_locus: Dict[str, List[object]], out_stem: Path, title: str) -> None:
    n_loci = len(matrices_by_locus)
    rendered = {locus: combined_locus_rgb(mats) for locus, mats in matrices_by_locus.items()}

    # Size the raster so every fiber gets its own pixel row at the export dpi.
    max_rows = max(rgb.shape[0] for rgb, _boundaries, _labels in rendered.values())
    raster_h = max_rows / CONDENSED_EXPORT_DPI
    gene_h = 0.12 * raster_h
    fig_w = max(7.6, 3.6 * n_loci)
    fig_h = raster_h + gene_h + 1.0

    fig = plt.figure(figsize=(fig_w, fig_h))
    gs = fig.add_gridspec(
        2,
        n_loci,
        height_ratios=[1, 0.12],
        width_ratios=[1] * n_loci,
        top=1 - 0.5 / fig_h,
        bottom=0.5 / fig_h,
        hspace=0.02,
        wspace=0.36,
    )

    for col, locus in enumerate(matrices_by_locus):
        mats = matrices_by_locus[locus]
        rgb, boundaries, labels = rendered[locus]
        locus_cfg = mats[0].config
        total_rows = rgb.shape[0]
        ax = fig.add_subplot(gs[0, col])
        ax.imshow(
            rgb,
            aspect="auto",
            interpolation="nearest",
            origin="upper",
            extent=[locus_cfg.start, locus_cfg.end, total_rows, 0],
            rasterized=False,
        )
        ax.axvline(locus_cfg.tss, color="#111111", lw=0.6, ls=":")
        for b in boundaries:
            ax.axhline(b, color="#333333", lw=0.55)
        add_condition_rails(ax, labels, total_rows)
        for condition, y0, y1, value, metric_label in labels:
            y_frac = 1.0 - ((y0 + y1) / 2 / total_rows)
            ax.text(1.01, y_frac, f"{condition}\n{value} {metric_label}", transform=ax.transAxes, ha="left", va="center", fontsize=6)
        ax.set_title(locus, fontsize=8, pad=5)
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ax.spines.values():
            spine.set_linewidth(0.7)
        ax_gene = fig.add_subplot(gs[1, col])
        mm.draw_gene_track(ax_gene, locus, locus_cfg.chrom, locus_cfg.start, locus_cfg.end, fig_width=3.5)

    fig.suptitle(title, fontsize=8, y=0.995)
    fig.text(0.5, 0.5 / fig_h * 0.4, f"All MSP calls use the same neutral gray ({MSP_GRAY}); condition colors are categorical rails only.", ha="center", fontsize=6, color="#444444")
    fs.save(fig, out_stem)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fig 4A stacked single-fiber MSP rasters.")
    parser.add_argument("--loci", nargs="+", default=["IDO1", "CXCL11", "STAT1"])
    parser.add_argument("--no-cluster-loci", nargs="*", default=["STAT1"])
    parser.add_argument("--trim-left-bp", type=int, default=0)
    parser.add_argument("--trim-right-bp", type=int, default=100)
    parser.add_argument("--bin-size", type=int, default=10)
    parser.add_argument("--min-coverage-frac", type=float, default=0.60)
    parser.add_argument("--sort-window-bp", type=int, default=250)
    parser.add_argument("--title", default="Continuous stacked MSP rasters: IDO1, CXCL11, STAT1")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.loci = [x.upper() for x in args.loci]
    args.no_cluster_loci = [x.upper() for x in args.no_cluster_loci]

    fs.setup()
    matplotlib.rcParams.update(ILLUSTRATOR_SAFE_RCPARAMS)

    matrices = build_matrices(args)
    out_stem = cfg.out_dir("Fig4A") / "Fig4A"
    plot_stacked_rasters(matrices, out_stem, args.title)
    print(out_stem)


if __name__ == "__main__":
    main()
