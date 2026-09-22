#!/usr/bin/env python3
"""Fig S8C: aggregate FiberHMM TF-call density across the CXCL11 amplicon, 50 bp bins.

Also holds the CXCL11 locus config, the pure-Python BAM reader and the per-base density
helpers that make_FigS8D.py imports.

pysam has no Windows build, so the spanning-fiber denominator is read from the recalled BAM
with a small pure-Python BGZF/BAM reader. It reproduces the pipeline's counts exactly: Naive
full-span over chr4:76,035,970-76,036,170 = 4485, which is the naive_eligible_fibers value in
the published ranking table.

Strand note: CXCL11 is on the MINUS strand. The panel plots genomic bp relative to the TSS, so
POSITIVE x = genomic downstream = biologically upstream (the promoter). The axis is deliberately
NOT mirrored to match IDO1.
"""
import csv
import gzip
import struct
import sys
from collections import OrderedDict, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import figstyle as fs

FIBERHMM_DIR = cfg.path("dafseq", "fiberhmm_dir")
SCOUT_DENSITY = Path(__file__).resolve().parents[1] / "data" / "cxcl11_scout_tf_density_50bp.tsv"

# data/locus_configs.tsv
LOCUS = "CXCL11"
CHROM = "chr4"
REGION_START = 76_034_865
REGION_END = 76_038_554
TSS = 76_036_070
STRAND = "-"
EDGE_TRIM_BP = 150
EDGE_TRIM_START = REGION_START + EDGE_TRIM_BP
EDGE_TRIM_END = REGION_END - EDGE_TRIM_BP

# data/locus_sample_sets.tsv
SAMPLES = OrderedDict(
    [
        ("Naive", "GB-270_CXCL11_Naive"),
        ("Stimulated", "GB-270_CXCL11_STIM"),
        ("Rested", "GB-270_CXCL11_REST"),
        ("Rechallenged", "GB-270_CXCL11_RECHALLENGE"),
    ]
)

COLORS = {
    "Naive": "#909499",
    "Stimulated": "#4D9AD4",
    "Rested": "#4A4A4B",
    "Rechallenged": "#0B519D",
}

FONT_LABEL = 8
FONT_NUMBER = 6

STRAND_NOTE = (
    f"{LOCUS} is on the minus strand; x is genomic bp relative to the TSS "
    "(+ = genomic downstream = biologically upstream / promoter)"
)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def tidy_axes(ax) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


CIGAR_CONSUMES_REF = {0, 2, 3, 7, 8}  # M, D, N, =, X


def iter_bam_alignments(bam_path: Path):
    """Yield (read_name, ref_name, ref_start, ref_end, flag) for every BAM record."""
    with gzip.open(str(bam_path), "rb") as handle:
        if handle.read(4) != b"BAM\x01":
            raise ValueError(f"{bam_path} is not a BAM file")
        (l_text,) = struct.unpack("<i", handle.read(4))
        handle.read(l_text)
        (n_ref,) = struct.unpack("<i", handle.read(4))
        refs = []
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", handle.read(4))
            refs.append(handle.read(l_name)[:-1].decode())
            handle.read(4)  # l_ref
        while True:
            head = handle.read(4)
            if len(head) < 4:
                break
            (block_size,) = struct.unpack("<i", head)
            block = handle.read(block_size)
            ref_id, pos, l_read_name, _mapq, _bin, n_cigar, flag = struct.unpack(
                "<iiBBHHH", block[:16]
            )
            offset = 32
            read_name = block[offset : offset + l_read_name - 1].decode()
            offset += l_read_name
            ref_len = 0
            for i in range(n_cigar):
                (op_len,) = struct.unpack("<I", block[offset + 4 * i : offset + 4 * i + 4])
                if (op_len & 0xF) in CIGAR_CONSUMES_REF:
                    ref_len += op_len >> 4
            yield read_name, refs[ref_id] if ref_id >= 0 else None, pos, pos + ref_len, flag


def eligible_reads(bam_path: Path, chrom: str, start: int, end: int) -> set[str]:
    """Primary reads that fully span [start, end) - mirrors the pysam version."""
    reads = set()
    for name, ref, ref_start, ref_end, flag in iter_bam_alignments(bam_path):
        if flag & 0x4 or flag & 0x100 or flag & 0x800:  # unmapped/secondary/supplementary
            continue
        if ref != chrom:
            continue
        if ref_start <= start and ref_end >= end:
            reads.add(name)
    return reads


def parse_bed12_blocks(path: Path):
    """Same BED12 block walk as make_FigS8B.py."""
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            chrom = fields[0]
            start = int(fields[1])
            read = fields[3]
            if len(fields) >= 12:
                block_count = int(fields[9])
                sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
                starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
                for size, rel in zip(sizes[:block_count], starts[:block_count]):
                    block_start = start + rel
                    yield {
                        "chrom": chrom,
                        "start": block_start,
                        "end": block_start + size,
                        "read": read,
                    }
            else:
                yield {
                    "chrom": chrom,
                    "start": start,
                    "end": int(fields[2]),
                    "read": read,
                }


def load_condition_data(chrom: str, start: int, end: int) -> dict[str, dict]:
    """Mirrors make_FigS8B.py::load_condition_data."""
    data = {}
    for condition, sample in SAMPLES.items():
        bam = FIBERHMM_DIR / f"{sample}.fiberhmm.recalled.bam"
        tf_bed = FIBERHMM_DIR / f"{sample}.fiberhmm.recalled_tf.bed"
        eligible = eligible_reads(bam, chrom, start, end)
        blocks_by_read = defaultdict(list)
        for block in parse_bed12_blocks(tf_bed):
            if block["chrom"] != chrom or block["end"] <= start or block["start"] >= end:
                continue
            if block["read"] in eligible:
                blocks_by_read[block["read"]].append(block)
        data[condition] = {
            "sample": sample,
            "blocks_by_read": blocks_by_read,
            "n_eligible": len(eligible),
        }
        print(f"  {condition:13s} {sample:26s} spanning fibers = {len(eligible)}")
    return data


def per_base_density(condition_data: dict, start: int, end: int) -> dict[str, np.ndarray]:
    """Verbatim from make_FigS8B.py::per_base_density."""
    width = end - start
    densities = {}
    for condition, data in condition_data.items():
        counts = np.zeros(width, dtype=float)
        for _read_name, blocks in data["blocks_by_read"].items():
            for block in blocks:
                s = max(block["start"], start) - start
                e = min(block["end"], end) - start
                if s < e:
                    counts[s:e] += 1
        densities[condition] = counts / max(data["n_eligible"], 1)
    return densities

def build_panel_a(out_stem: Path) -> None:
    rows = read_tsv(SCOUT_DENSITY)

    fig, ax = plt.subplots(figsize=(6.5, 2.1))
    for condition in SAMPLES:
        kept = [
            r
            for r in rows
            if r["condition"] == condition
            and EDGE_TRIM_START <= float(r["mid"]) <= EDGE_TRIM_END
        ]
        kept.sort(key=lambda r: float(r["mid"]))
        x = np.array([float(r["mid_relative_to_tss"]) for r in kept])
        y = np.array([float(r["tf_read_fraction"]) for r in kept])
        ax.plot(x, y, lw=1.0, color=COLORS[condition], label=condition)
        if condition == "Naive":
            print(f"  bins kept after {EDGE_TRIM_BP} bp edge trim: {len(kept)}"
                  f" ({x.min():+.0f} to {x.max():+.0f} bp)")

    ax.axvline(0, color="#222222", lw=0.7, ls="--", label=f"{LOCUS} TSS")
    ax.axvspan(-75, 75, color="#777777", alpha=0.08, lw=0)
    ax.set_xlim(EDGE_TRIM_START - TSS, EDGE_TRIM_END - TSS)
    ax.set_ylim(bottom=0)
    ax.set_xlabel(f"Position relative to {LOCUS} TSS (bp)")
    ax.set_ylabel("Fraction of spanning reads\nwith TF-sized call")
    ax.legend(frameon=False, ncols=5, loc="upper center", bbox_to_anchor=(0.5, 1.22),
              handlelength=1.4, columnspacing=1.1, borderpad=0)
    ax.annotate(STRAND_NOTE, xy=(0.5, -0.42), xycoords="axes fraction",
                ha="center", va="top", fontsize=FONT_NUMBER, color="#555555")
    tidy_axes(ax)
    fs.save(fig, out_stem)


def main() -> None:
    fs.setup()
    build_panel_a(cfg.out_dir("FigS8C") / "FigS8C")


if __name__ == "__main__":
    main()
