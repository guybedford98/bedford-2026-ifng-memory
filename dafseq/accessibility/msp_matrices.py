"""Per-fiber binary MSP matrices from FiberHMM MSP calls, and the balanced k=2
MSP-rich/MSP-low classification used for Fig 4A-C.

Imported by the numbered stages and panel scripts in this folder.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
import zlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.patches as mpatches
import numpy as np
import pysam

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))
import config as cfg

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
LOCUS_CONFIGS = DATA_DIR / "locus_configs.tsv"
LOCUS_SAMPLE_SETS = DATA_DIR / "locus_sample_sets.tsv"
SAMPLE_MANIFEST = DATA_DIR / "sample_manifest.tsv"

CONDITION_ORDER = ["Naive", "Stimulated", "Rested", "Rechallenged"]
CONDITION_COLORS = {
    "Naive": "#909499",
    "Stimulated": "#4D9AD4",
    "Rested": "#4A4A4B",
    "Rechallenged": "#0B519D",
}

GENE_PRESETS = {
    "CXCL11": {
        "name": "CXCL11",
        "strand": "-",
        "tss": 76036070,
        "exons": [
            (76033681, 76034816),
            (76035046, 76035119),
            (76035215, 76035342),
            (76035926, 76036070),
        ],
        "tx_start": 76033681,
        "tx_end": 76036070,
    },
    "STAT1": {
        "name": "STAT1",
        "strand": "-",
        "tss": 191014171,
        "exons": [
            (191013524, 191013678),
            (191014017, 191014171),
        ],
        "tx_start": 190969148,
        "tx_end": 191014171,
    },
    "IDO1": {
        "name": "IDO1",
        "strand": "+",
        "tss": 39913890,
        "exons": [
            (39913890, 39914009),
            (39917874, 39917970),
            (39918087, 39918207),
            (39918814, 39918933),
            (39920099, 39920114),
            (39922551, 39922651),
            (39923468, 39923586),
            (39924720, 39924772),
            (39925222, 39925371),
            (39927829, 39928790),
        ],
        "tx_start": 39913890,
        "tx_end": 39928790,
    },
}


@dataclass(frozen=True)


class LocusConfig:
    locus: str
    chrom: str
    start: int
    end: int
    tss: int
    strand: str
    edge_trim_bp: int

    def trimmed(self, left_bp: int = 0, right_bp: int = 0) -> "LocusConfig":
        new_start = self.start + max(0, left_bp)
        new_end = self.end - max(0, right_bp)
        if new_start >= new_end:
            raise ValueError(
                f"Trimmed region for {self.locus} is empty: "
                f"{self.chrom}:{new_start}-{new_end}"
            )
        return LocusConfig(
            locus=self.locus,
            chrom=self.chrom,
            start=new_start,
            end=new_end,
            tss=self.tss,
            strand=self.strand,
            edge_trim_bp=self.edge_trim_bp,
        )


@dataclass(frozen=True)


class SampleConfig:
    locus: str
    condition: str
    sample: str
    order: int
    decorated_bam: Path
    msp_bed: Path


@dataclass(frozen=True)


class ReadInfo:
    name: str
    start: int
    end: int
    strand_tag: str


@dataclass


class SampleMatrix:
    locus: str
    condition: str
    sample: str
    strand_view: str
    config: LocusConfig
    sample_config: SampleConfig
    bin_edges: np.ndarray
    bin_centers: np.ndarray
    matrix: np.ndarray
    selected_matrix: np.ndarray
    selected_indices: np.ndarray
    selected_names: List[str]
    reads: List[ReadInfo]
    eligible_reads: int
    normalized_depth: int
    selected_reads: int
    total_msp_blocks: int
    matrix_mean: float
    selected_mean: float
    high_state_count: int
    low_state_count: int
    msp_state_mode: str = "clustered"
    global_cluster_train_reads_per_condition: int = 0


def parse_region(region: str) -> Tuple[str, int, int]:
    chrom, rest = region.split(":", 1)
    start_s, end_s = rest.replace(",", "").split("-", 1)
    return chrom, int(start_s), int(end_s)


def read_tsv(path: Path) -> List[Dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_locus_configs() -> Dict[str, LocusConfig]:
    rows = read_tsv(LOCUS_CONFIGS)
    configs: Dict[str, LocusConfig] = {}
    for row in rows:
        chrom, start, end = parse_region(row["region"])
        locus = row["locus"]
        configs[locus] = LocusConfig(
            locus=locus,
            chrom=chrom,
            start=start,
            end=end,
            tss=int(row["tss"]),
            strand=row["strand"],
            edge_trim_bp=int(row.get("edge_trim_bp") or 0),
        )
    return configs


def load_sample_configs() -> Dict[str, List[SampleConfig]]:
    bam_dir = cfg.path("dafseq", "bam_dir")
    fiberhmm_dir = cfg.path("dafseq", "fiberhmm_dir")
    manifest = {row["sample"]: row for row in read_tsv(SAMPLE_MANIFEST)}
    sample_rows = read_tsv(LOCUS_SAMPLE_SETS)

    by_locus: Dict[str, List[SampleConfig]] = {}
    for row in sample_rows:
        sample = row["sample"]
        if sample not in manifest:
            raise ValueError(f"Sample {sample} is missing from sample_manifest.tsv")
        cfg_row = SampleConfig(
            locus=row["locus"],
            condition=row["condition"],
            sample=sample,
            order=int(row["order"]),
            decorated_bam=bam_dir / manifest[sample]["decorated_bam"],
            msp_bed=fiberhmm_dir / f"{sample}.fiberhmm.recalled_msp.bed",
        )
        by_locus.setdefault(row["locus"], []).append(cfg_row)

    for locus in by_locus:
        by_locus[locus].sort(
            key=lambda c: (CONDITION_ORDER.index(c.condition), c.order)
            if c.condition in CONDITION_ORDER
            else (99, c.order)
        )
    return by_locus


def parse_bed12_blocks(
    bed_path: Path,
    chrom: str,
    start: int,
    end: int,
) -> Tuple[Dict[str, List[Tuple[int, int]]], int]:
    by_read: Dict[str, List[Tuple[int, int]]] = {}
    total_blocks = 0
    with bed_path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 12:
                continue
            bed_chrom = fields[0]
            if bed_chrom != chrom:
                continue
            bed_start = int(fields[1])
            name = fields[3]
            block_count = int(fields[9])
            sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
            rel_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
            for size, rel in zip(sizes[:block_count], rel_starts[:block_count]):
                block_start = bed_start + rel
                block_end = block_start + size
                if block_end <= start or block_start >= end:
                    continue
                clipped = (max(block_start, start), min(block_end, end))
                by_read.setdefault(name, []).append(clipped)
                total_blocks += 1

    for name in by_read:
        by_read[name].sort()
    return by_read, total_blocks


def load_eligible_reads(
    bam_path: Path,
    chrom: str,
    start: int,
    end: int,
    strand_view: str,
    min_coverage_frac: float,
) -> List[ReadInfo]:
    reads: List[ReadInfo] = []
    region_len = end - start
    strand_view = strand_view.upper()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for read in bam.fetch(chrom, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            read_start = read.reference_start
            read_end = read.reference_end
            if read_start is None or read_end is None:
                continue
            overlap = max(0, min(read_end, end) - max(read_start, start))
            if overlap / region_len < min_coverage_frac:
                continue
            try:
                strand_tag = read.get_tag("st")
            except KeyError:
                strand_tag = "NA"
            if strand_view != "BOTH" and strand_tag.upper() != strand_view:
                continue
            reads.append(
                ReadInfo(
                    name=read.query_name,
                    start=max(read_start, start),
                    end=min(read_end, end),
                    strand_tag=strand_tag,
                )
            )

    reads.sort(key=lambda r: r.name)
    return reads


def make_bin_edges(start: int, end: int, bin_size: int) -> Tuple[np.ndarray, np.ndarray]:
    n_bins = int(math.ceil((end - start) / bin_size))
    edges = start + np.arange(n_bins + 1) * bin_size
    edges[-1] = end
    centers = (edges[:-1] + edges[1:]) / 2.0
    return edges.astype(int), centers


def fill_interval_matrix(
    reads: Sequence[ReadInfo],
    blocks_by_read: Dict[str, List[Tuple[int, int]]],
    bin_edges: np.ndarray,
) -> np.ndarray:
    n_reads = len(reads)
    n_bins = len(bin_edges) - 1
    mat = np.zeros((n_reads, n_bins), dtype=np.float32)
    if n_reads == 0 or n_bins == 0:
        return mat

    start = int(bin_edges[0])
    bin_size = int(bin_edges[1] - bin_edges[0]) if n_bins > 1 else int(bin_edges[-1] - bin_edges[0])

    for row_i, read in enumerate(reads):
        for block_start, block_end in blocks_by_read.get(read.name, []):
            if block_end <= bin_edges[0] or block_start >= bin_edges[-1]:
                continue
            first = max(0, int((block_start - start) // bin_size))
            last = min(n_bins - 1, int((block_end - 1 - start) // bin_size))
            for bin_i in range(first, last + 1):
                ov = max(
                    0,
                    min(block_end, int(bin_edges[bin_i + 1]))
                    - max(block_start, int(bin_edges[bin_i])),
                )
                denom = max(1, int(bin_edges[bin_i + 1] - bin_edges[bin_i]))
                mat[row_i, bin_i] += ov / denom

    np.clip(mat, 0.0, 1.0, out=mat)
    return mat


def row_focus_scores(
    matrix: np.ndarray,
    bin_centers: np.ndarray,
    tss: int,
    sort_window_bp: int,
) -> np.ndarray:
    if matrix.shape[0] == 0:
        return np.zeros(0, dtype=np.float32)
    mask = (bin_centers >= tss - sort_window_bp) & (bin_centers <= tss + sort_window_bp)
    if not mask.any():
        return matrix.sum(axis=1)
    return matrix[:, mask].sum(axis=1)


def select_evenly(indices: np.ndarray, target_n: int) -> np.ndarray:
    if target_n <= 0 or len(indices) == 0:
        return np.zeros(0, dtype=int)
    if len(indices) <= target_n:
        return indices
    keep_pos = np.linspace(0, len(indices) - 1, target_n).round().astype(int)
    keep_pos = np.unique(keep_pos)
    if len(keep_pos) < target_n:
        missing = target_n - len(keep_pos)
        extras = [i for i in range(len(indices)) if i not in set(keep_pos)]
        keep_pos = np.sort(np.concatenate([keep_pos, np.array(extras[:missing], dtype=int)]))
    return indices[keep_pos[:target_n]]


def build_sample_matrix(
    locus_cfg: LocusConfig,
    sample_cfg: SampleConfig,
    strand_view: str,
    bin_size: int,
    min_coverage_frac: float,
) -> SampleMatrix:
    blocks_by_read, total_blocks = parse_bed12_blocks(
        sample_cfg.msp_bed, locus_cfg.chrom, locus_cfg.start, locus_cfg.end
    )
    reads = load_eligible_reads(
        sample_cfg.decorated_bam,
        locus_cfg.chrom,
        locus_cfg.start,
        locus_cfg.end,
        strand_view,
        min_coverage_frac,
    )
    bin_edges, bin_centers = make_bin_edges(locus_cfg.start, locus_cfg.end, bin_size)
    matrix = fill_interval_matrix(reads, blocks_by_read, bin_edges)
    matrix = (matrix > 0).astype(np.float32)
    return SampleMatrix(
        locus=locus_cfg.locus,
        condition=sample_cfg.condition,
        sample=sample_cfg.sample,
        strand_view=strand_view,
        config=locus_cfg,
        sample_config=sample_cfg,
        bin_edges=bin_edges,
        bin_centers=bin_centers,
        matrix=matrix,
        selected_matrix=matrix,
        selected_indices=np.arange(matrix.shape[0], dtype=int),
        selected_names=[r.name for r in reads],
        reads=list(reads),
        eligible_reads=matrix.shape[0],
        normalized_depth=matrix.shape[0],
        selected_reads=matrix.shape[0],
        total_msp_blocks=total_blocks,
        matrix_mean=float(matrix.mean()) if matrix.size else 0.0,
        selected_mean=float(matrix.mean()) if matrix.size else 0.0,
        high_state_count=0,
        low_state_count=0,
    )


def fit_balanced_global_kmeans(
    matrices: List[SampleMatrix],
    random_state: int,
) -> Tuple[Dict[str, np.ndarray], int]:
    """
    Fit one k=2 MSP-state model across all conditions for this locus/strand view.

    Each condition contributes the same number of rows to the fitting set, so a
    deeper sample cannot define the MSP-rich/low boundary by itself.
    """
    nonempty = [m for m in matrices if m.matrix.shape[0] > 0]
    if not nonempty:
        return {m.condition: np.zeros(0, dtype=bool) for m in matrices}, 0

    train_n = min(m.matrix.shape[0] for m in nonempty)
    rng = np.random.default_rng(random_state)
    train_rows = []
    for m in nonempty:
        if m.matrix.shape[0] == train_n:
            idx = np.arange(train_n)
        else:
            idx = rng.choice(m.matrix.shape[0], train_n, replace=False)
        train_rows.append(m.matrix[idx, :])

    training = np.vstack(train_rows)
    if training.shape[0] < 2 or np.allclose(training, training[0, :]):
        labels = {
            m.condition: np.ones(m.matrix.shape[0], dtype=bool)
            for m in matrices
        }
        return labels, train_n

    try:
        from sklearn.cluster import KMeans

        km = KMeans(n_clusters=2, random_state=random_state, n_init=20)
        km.fit(training)
        high_cluster = int(np.argmax(km.cluster_centers_.mean(axis=1)))
        labels = {
            m.condition: (km.predict(m.matrix) == high_cluster)
            for m in matrices
        }
    except Exception:
        # Lightweight fallback: split by the midpoint between two 1D k-means
        # centers on total MSP burden. This keeps the workflow usable if sklearn
        # is not available in a future environment.
        scores = training.mean(axis=1)
        c0, c1 = np.quantile(scores, [0.25, 0.75])
        for _ in range(50):
            d0 = np.abs(scores - c0)
            d1 = np.abs(scores - c1)
            assign = d1 < d0
            if assign.any():
                c1 = scores[assign].mean()
            if (~assign).any():
                c0 = scores[~assign].mean()
        low_c, high_c = sorted([c0, c1])
        threshold = (low_c + high_c) / 2
        labels = {
            m.condition: (m.matrix.mean(axis=1) >= threshold)
            for m in matrices
        }

    return labels, train_n


def apply_global_msp_state_order(
    matrices: List[SampleMatrix],
    max_fibers_per_condition: int,
    sort_window_bp: int,
    display_depth_mode: str,
    random_state: int,
    cluster_mode: str = "clustered",
) -> None:
    counts = [m.eligible_reads for m in matrices]
    positive_counts = [c for c in counts if c > 0]
    if display_depth_mode == "equal_subsample":
        target_n = min(positive_counts) if positive_counts else 0
        if max_fibers_per_condition > 0:
            target_n = min(target_n, max_fibers_per_condition)
    else:
        target_n = max(positive_counts) if positive_counts else 0
        if max_fibers_per_condition > 0:
            target_n = min(target_n, max_fibers_per_condition)

    if cluster_mode == "clustered":
        high_by_condition, train_n = fit_balanced_global_kmeans(matrices, random_state)
    else:
        high_by_condition = {
            m.condition: np.zeros(m.matrix.shape[0], dtype=bool)
            for m in matrices
        }
        train_n = 0

    for m in matrices:
        names = [read.name for read in m.reads]
        high = high_by_condition.get(m.condition, np.zeros(m.matrix.shape[0], dtype=bool))
        focus = row_focus_scores(m.matrix, m.bin_centers, m.config.tss, sort_window_bp)
        total = m.matrix.sum(axis=1)
        name_order = np.array([zlib.crc32(name.encode("utf-8")) for name in names], dtype=np.uint64)

        if cluster_mode == "clustered":
            high_idx = np.where(high)[0]
            low_idx = np.where(~high)[0]
            high_idx = high_idx[np.lexsort((name_order[high_idx], -total[high_idx], -focus[high_idx]))]
            low_idx = low_idx[np.lexsort((name_order[low_idx], -total[low_idx], -focus[low_idx]))]
            order = np.concatenate([high_idx, low_idx])
        else:
            order = np.lexsort((name_order, -total, -focus))

        if display_depth_mode == "equal_subsample":
            selected = select_evenly(order, target_n)
        elif max_fibers_per_condition > 0 and len(order) > max_fibers_per_condition:
            selected = select_evenly(order, max_fibers_per_condition)
        else:
            selected = order

        selected_mat = m.matrix[selected, :] if len(selected) else np.empty((0, m.matrix.shape[1]))
        selected_high = high[selected] if len(selected) else np.zeros(0, dtype=bool)

        m.selected_indices = selected
        m.selected_matrix = selected_mat
        m.selected_names = [m.reads[i].name for i in selected] if len(selected) else []
        m.selected_reads = int(selected_mat.shape[0])
        m.selected_mean = float(selected_mat.mean()) if selected_mat.size else 0.0
        m.normalized_depth = int(max(target_n, m.selected_reads))
        if cluster_mode == "clustered":
            m.high_state_count = int(selected_high.sum())
            m.low_state_count = int((~selected_high).sum())
        else:
            m.high_state_count = 0
            m.low_state_count = 0
        m.msp_state_mode = cluster_mode
        m.global_cluster_train_reads_per_condition = int(train_n)


def draw_gene_track(ax, locus: str, chrom: str, region_start: int, region_end: int, fig_width: float) -> None:
    gene_info = GENE_PRESETS.get(locus)
    if gene_info is None:
        ax.set_visible(False)
        return

    name = gene_info["name"]
    strand = gene_info["strand"]
    exons = sorted(gene_info["exons"])
    tx_s = int(gene_info["tx_start"])
    tx_e = int(gene_info["tx_end"])
    region_w = region_end - region_start
    scale = min(1.0, fig_width / 13.0)
    fs_label = max(5.5, 8 * scale)

    ax.set_xlim(region_start, region_end)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)

    body_s = max(tx_s, region_start)
    body_e = min(tx_e, region_end)
    if body_s < body_e:
        ax.plot([body_s, body_e], [0.48, 0.48], color="#555555", lw=1.0, solid_capstyle="butt", zorder=1)

    arrow_step = max(region_w / 14, 1)
    arrow_dx = arrow_step * 0.28
    arrowprops = dict(arrowstyle="-|>", color="#888888", mutation_scale=5, lw=0.6)
    x = body_s + arrow_step * 0.5
    while x < body_e - arrow_step * 0.1:
        in_exon = any(max(es, region_start) <= x <= min(ee, region_end) for es, ee in exons)
        if not in_exon:
            if strand == "-":
                src, dst = x + arrow_dx, x - arrow_dx
            else:
                src, dst = x - arrow_dx, x + arrow_dx
            ax.annotate("", xy=(dst, 0.48), xytext=(src, 0.48), arrowprops=arrowprops, zorder=2)
        x += arrow_step

    exon_h = 0.20
    for ex_s, ex_e in exons:
        vis_s = max(ex_s, region_start)
        vis_e = min(ex_e, region_end)
        if vis_s >= vis_e:
            continue
        ax.add_patch(
            mpatches.FancyBboxPatch(
                (vis_s, 0.48 - exon_h / 2),
                vis_e - vis_s,
                exon_h,
                boxstyle="square,pad=0",
                facecolor="#2166ac",
                edgecolor="none",
                zorder=3,
            )
        )

    direction = "<" if strand == "-" else ">"
    visible_lefts = [
        max(es, region_start)
        for es, ee in exons
        if min(ee, region_end) > max(es, region_start)
    ]
    label_x = min(visible_lefts) if visible_lefts else body_s
    ha = "right" if strand == "-" else "left"
    ax.text(label_x, 0.70, f"{name} {direction}", ha=ha, va="bottom", fontsize=fs_label, color="#333333", fontstyle="italic")

    ax.set_xticks([region_start, region_end])
    ax.set_xticklabels([f"{chrom}:{region_start:,}", f"{chrom}:{region_end:,}"], fontsize=6)
    ax.tick_params(axis="x", length=4, width=0.8)


def write_summary_tsv(
    path: Path,
    matrices: Sequence[SampleMatrix],
    bin_size: int,
    min_coverage_frac: float,
    max_fibers: int,
    trim_left_bp: int,
    trim_right_bp: int,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "locus",
        "condition",
        "sample",
        "strand_view",
        "region",
        "trim_left_bp",
        "trim_right_bp",
        "tss",
        "eligible_reads",
        "normalized_display_depth",
        "selected_reads",
        "msp_state_mode",
        "high_msp_state_reads",
        "low_msp_state_reads",
        "msp_rich_fraction_of_eligible",
        "global_cluster_train_reads_per_condition",
        "total_msp_blocks_in_region",
        "all_eligible_mean_msp_covered_fiber_area_fraction",
        "selected_mean_msp_covered_fiber_area_fraction",
        "bin_size_bp",
        "min_coverage_frac",
        "max_fibers_per_condition",
        "decorated_bam",
        "msp_bed",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for m in matrices:
            is_clustered = m.msp_state_mode == "clustered"
            writer.writerow(
                {
                    "locus": m.locus,
                    "condition": m.condition,
                    "sample": m.sample,
                    "strand_view": m.strand_view,
                    "region": f"{m.config.chrom}:{m.config.start}-{m.config.end}",
                    "trim_left_bp": trim_left_bp,
                    "trim_right_bp": trim_right_bp,
                    "tss": m.config.tss,
                    "eligible_reads": m.eligible_reads,
                    "normalized_display_depth": m.normalized_depth,
                    "selected_reads": m.selected_reads,
                    "msp_state_mode": m.msp_state_mode,
                    "high_msp_state_reads": m.high_state_count if is_clustered else "",
                    "low_msp_state_reads": m.low_state_count if is_clustered else "",
                    "msp_rich_fraction_of_eligible": f"{(m.high_state_count / m.selected_reads) if m.selected_reads else 0.0:.6f}" if is_clustered else "",
                    "global_cluster_train_reads_per_condition": m.global_cluster_train_reads_per_condition if is_clustered else "",
                    "total_msp_blocks_in_region": m.total_msp_blocks,
                    "all_eligible_mean_msp_covered_fiber_area_fraction": f"{m.matrix_mean:.6f}",
                    "selected_mean_msp_covered_fiber_area_fraction": f"{m.selected_mean:.6f}",
                    "bin_size_bp": bin_size,
                    "min_coverage_frac": min_coverage_frac,
                    "max_fibers_per_condition": max_fibers,
                    "decorated_bam": m.sample_config.decorated_bam.name,
                    "msp_bed": m.sample_config.msp_bed.name,
                }
            )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="MSP matrix build options shared by the accessibility stages.")
    parser.add_argument("--loci", nargs="+", default=["IDO1", "CXCL11", "STAT1"])
    parser.add_argument("--strand-views", nargs="+", default=["both"])
    parser.add_argument("--trim-left-bp", type=int, default=0)
    parser.add_argument("--trim-right-bp", type=int, default=100)
    parser.add_argument("--bin-size", type=int, default=10)
    parser.add_argument("--min-coverage-frac", type=float, default=0.60)
    parser.add_argument("--display-depth-mode", choices=["equal_subsample", "all"], default="equal_subsample")
    parser.add_argument("--max-fibers-per-condition", type=int, default=0)
    parser.add_argument("--sort-window-bp", type=int, default=250)
    parser.add_argument("--no-cluster-loci", nargs="*", default=["STAT1"])
    parser.add_argument("--n-permutations", type=int, default=5000)
    parser.add_argument("--n-bootstraps", type=int, default=3000)
    parser.add_argument("--random-seed", type=int, default=42)
    return parser.parse_args()


def build_matrices(args: argparse.Namespace) -> List[SampleMatrix]:
    configs = load_locus_configs()
    sample_sets = load_sample_configs()
    matrices_all: List[SampleMatrix] = []
    for locus in args.loci:
        base_cfg = configs[locus]
        locus_cfg = base_cfg.trimmed(args.trim_left_bp, args.trim_right_bp)
        for strand_view in args.strand_views:
            matrices: List[SampleMatrix] = []
            for sample_cfg in sample_sets[locus]:
                matrices.append(
                    build_sample_matrix(
                        locus_cfg,
                        sample_cfg,
                        strand_view,
                        args.bin_size,
                        args.min_coverage_frac,
                    )
                )
            apply_global_msp_state_order(
                matrices,
                args.max_fibers_per_condition,
                args.sort_window_bp,
                args.display_depth_mode,
                random_state=args.random_seed,
                cluster_mode="unclustered" if locus in args.no_cluster_loci else "clustered",
            )
            matrices_all.extend(matrices)
    return matrices_all
