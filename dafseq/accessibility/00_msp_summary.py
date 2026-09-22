#!/usr/bin/env python3
"""Build the per-fiber MSP matrices for IDO1, CXCL11 and STAT1 and write the
summary table that Fig 4B and Fig 4C are drawn from.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import msp_matrices as mm


def main() -> None:
    args = mm.parse_args()
    args.loci = [x.upper() for x in args.loci]
    args.no_cluster_loci = [x.upper() for x in args.no_cluster_loci]

    matrices = mm.build_matrices(args)
    summary_tsv = cfg.out_dir("tables") / "msp_heatmap_summary.tsv"
    mm.write_summary_tsv(
        summary_tsv,
        matrices,
        args.bin_size,
        args.min_coverage_frac,
        args.max_fibers_per_condition,
        args.trim_left_bp,
        args.trim_right_bp,
    )
    for m in matrices:
        print(f"{m.locus}\t{m.condition}\teligible={m.eligible_reads}\tselected={m.selected_reads}")
    print(f"Summary TSV: {summary_tsv}")


if __name__ == "__main__":
    main()
