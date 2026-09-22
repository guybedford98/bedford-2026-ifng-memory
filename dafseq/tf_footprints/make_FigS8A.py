#!/usr/bin/env python3
"""Fig S8A: aggregate FiberHMM TF-sized footprint density across the IDO1 amplicon."""
import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pysam

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import figstyle as fs

FIBERHMM_DIR = cfg.path("dafseq", "fiberhmm_dir")

CHROM = "chr8"
REGION_START = 39_910_187
REGION_END = 39_914_251
TSS = 39_913_890
BIN_SIZE = 50
EDGE_TRIM_BP = 150
EDGE_TRIM_START = REGION_START + EDGE_TRIM_BP
EDGE_TRIM_END = REGION_END - EDGE_TRIM_BP

SAMPLES = {
    "Naive": "GB-270_IDO1_Naive_merged",
    "Stimulated": "GB-270_IDO1_Stim_merged",
    "Rested": "GB-270_IDO1_Rested_merged",
    "Rechallenged": "GB-270_IDO1_Rechallenge_merged",
}

COLORS = {
    "Naive": "#909499",
    "Stimulated": "#4D9AD4",
    "Rested": "#4A4A4B",
    "Rechallenged": "#0B519D",
}


def sample_paths(sample: str) -> dict[str, Path]:
    return {
        "recalled_bam": FIBERHMM_DIR / f"{sample}.fiberhmm.recalled.bam",
        "tf_bed": FIBERHMM_DIR / f"{sample}.fiberhmm.recalled_tf.bed",
    }


def parse_bed12_blocks(path: Path):
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0]
            start = int(f[1])
            read = f[3]
            score = int(f[4])
            strand = f[5]
            if len(f) >= 12:
                block_count = int(f[9])
                sizes = [int(x) for x in f[10].rstrip(",").split(",") if x]
                starts = [int(x) for x in f[11].rstrip(",").split(",") if x]
                for size, rel in zip(sizes[:block_count], starts[:block_count]):
                    s = start + rel
                    yield {
                        "chrom": chrom,
                        "start": s,
                        "end": s + size,
                        "read": read,
                        "score": score,
                        "strand": strand,
                        "size": size,
                    }
            else:
                end = int(f[2])
                yield {
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "read": read,
                    "score": score,
                    "strand": strand,
                    "size": end - start,
                }


def load_blocks(sample: str) -> list[dict]:
    paths = sample_paths(sample)
    return [
        b
        for b in parse_bed12_blocks(paths["tf_bed"])
        if b["chrom"] == CHROM and b["end"] > REGION_START and b["start"] < REGION_END
    ]


def density_by_bin(sample: str, blocks: list[dict], bin_size: int = BIN_SIZE, start: int = REGION_START, end: int = REGION_END):
    paths = sample_paths(sample)
    bins = np.arange(start, end + bin_size, bin_size)
    if bins[-1] < end:
        bins = np.append(bins, end)
    mids = (bins[:-1] + bins[1:]) / 2
    spanning = np.zeros(len(mids), dtype=int)
    covered_bp = np.zeros(len(mids), dtype=int)
    tf_reads = [set() for _ in mids]
    protected_bp = np.zeros(len(mids), dtype=int)

    with pysam.AlignmentFile(str(paths["recalled_bam"]), "rb") as bam:
        for read in bam.fetch(CHROM, start, end):
            if read.is_unmapped or read.is_secondary or read.is_supplementary:
                continue
            rs = max(read.reference_start, start)
            re = min(read.reference_end, end)
            if rs >= re:
                continue
            i0 = max(0, np.searchsorted(bins, rs, side="right") - 1)
            i1 = min(len(mids) - 1, np.searchsorted(bins, re - 1, side="right") - 1)
            for i in range(i0, i1 + 1):
                ov = max(0, min(re, bins[i + 1]) - max(rs, bins[i]))
                if ov:
                    spanning[i] += 1
                    covered_bp[i] += ov

    for block in blocks:
        bs = max(block["start"], start)
        be = min(block["end"], end)
        if bs >= be:
            continue
        i0 = max(0, np.searchsorted(bins, bs, side="right") - 1)
        i1 = min(len(mids) - 1, np.searchsorted(bins, be - 1, side="right") - 1)
        for i in range(i0, i1 + 1):
            ov = max(0, min(be, bins[i + 1]) - max(bs, bins[i]))
            if ov:
                tf_reads[i].add(block["read"])
                protected_bp[i] += ov

    rows = []
    for i in range(len(mids)):
        rows.append(
            {
                "chrom": CHROM,
                "start": int(bins[i]),
                "end": int(bins[i + 1]),
                "mid": round(float(mids[i]), 1),
                "mid_relative_to_tss": round(float(mids[i] - TSS), 1),
                "spanning_reads": int(spanning[i]),
                "reads_with_tf": len(tf_reads[i]),
                "covered_bp": int(covered_bp[i]),
                "protected_bp": int(protected_bp[i]),
                "tf_read_fraction": round(len(tf_reads[i]) / max(int(spanning[i]), 1), 6),
                "tf_bp_fraction": round(int(protected_bp[i]) / max(int(covered_bp[i]), 1), 6),
            }
        )
    return rows


def write_tsv(path: Path, rows: list[dict]):
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def plot_amplicon_density(density_rows: list[dict], out_stem: Path):
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    for label in SAMPLES:
        rows = [
            r
            for r in density_rows
            if r["condition"] == label and EDGE_TRIM_START <= r["mid"] <= EDGE_TRIM_END
        ]
        x = np.array([r["mid_relative_to_tss"] for r in rows])
        y = np.array([r["tf_read_fraction"] for r in rows])
        ax.plot(x, y, lw=1.9, color=COLORS[label], label=label)
    ax.axvline(0, color="#222222", lw=1.0, ls="--", label="IDO1 TSS")
    ax.axvspan(-75, 75, color="#777777", alpha=0.08, lw=0)
    ax.set_xlim(EDGE_TRIM_START - TSS, EDGE_TRIM_END - TSS)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Position relative to IDO1 TSS (bp)")
    ax.set_ylabel("Fraction of spanning reads\nwith TF-sized call")
    ax.legend(frameon=False, ncols=5, fontsize=6)
    fs.prism(ax)
    fs.save(fig, out_stem)


def main() -> None:
    fs.setup()
    out_dir = cfg.out_dir("FigS8A")
    density_rows = []
    for label, sample in SAMPLES.items():
        blocks = load_blocks(sample)
        for row in density_by_bin(sample, blocks):
            row["condition"] = label
            row["sample"] = sample
            density_rows.append(row)

    edge_trimmed_density_rows = [
        row for row in density_rows if EDGE_TRIM_START <= row["mid"] <= EDGE_TRIM_END
    ]
    write_tsv(out_dir / "ido1_edge_trimmed_tf_density_50bp.tsv", edge_trimmed_density_rows)
    plot_amplicon_density(density_rows, out_dir / "FigS8A")
    print(out_dir)


if __name__ == "__main__":
    main()
