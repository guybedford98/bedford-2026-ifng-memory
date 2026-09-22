#!/usr/bin/env python3
"""TSS-localized MSP frequency, the percentages quoted in the Results.

For each locus x condition, take the equal-depth selected binary MSP matrix and a TSS-centered
window (TSS +/- W bp). Report:
  mean_tss_freq_pct : mean fraction of the TSS window covered by MSP, averaged across fibers
                      ("MSP frequency at the TSS")
  frac_open_pct     : fraction of fibers whose TSS window is majority (>=50%) MSP-covered
                      ("fraction of fibers with an accessible TSS region")
"""
import csv
import sys
from pathlib import Path

import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import msp_matrices as mm


def main() -> None:
    args = mm.parse_args()
    args.loci = [x.upper() for x in args.loci]
    args.no_cluster_loci = [x.upper() for x in args.no_cluster_loci]

    mats = mm.build_matrices(args)
    rows = []
    for W in (150, 250):
        print(f"\n=== TSS window +/- {W} bp ===")
        print("locus\tcondition\tn_fibers\tmean_tss_freq_pct\tfrac_open_pct")
        for m in mats:
            centers = np.asarray(m.bin_centers, dtype=float)
            mat = np.asarray(m.selected_matrix, dtype=float)
            tss = float(m.config.tss)
            mask = np.abs(centers - tss) <= W
            if mask.sum() == 0 or mat.shape[0] == 0:
                continue
            per_fiber = mat[:, mask].mean(axis=1)        # fraction of TSS window accessible, per fiber
            mean_freq = 100.0 * per_fiber.mean()
            frac_open = 100.0 * (per_fiber >= 0.5).mean()
            print(f"{m.locus}\t{m.condition}\t{mat.shape[0]}\t{mean_freq:.1f}\t{frac_open:.1f}")
            rows.append(
                {
                    "tss_window_bp": W,
                    "locus": m.locus,
                    "condition": m.condition,
                    "n_fibers": mat.shape[0],
                    "mean_tss_freq_pct": f"{mean_freq:.1f}",
                    "frac_open_pct": f"{frac_open:.1f}",
                }
            )

    path = cfg.out_dir("tables") / "tss_msp_accessibility.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    print(f"\n{path}")


if __name__ == "__main__":
    main()
