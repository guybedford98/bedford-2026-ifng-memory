#!/usr/bin/env python3
"""Molecule-level statistics on per-fiber MSP-covered area.

The observation being tested is one selected fiber's MSP-covered area:
    selected_matrix.mean(axis=1) * 100

These are not biological-replicate inference.
"""
import csv
import math
import sys
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np
from scipy import stats

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "shared"))

import config as cfg
import msp_matrices as mm
from msp_matrices import CONDITION_ORDER, SampleMatrix

TARGETED_CONTRASTS = [
    ("Naive", "Stimulated"),
    ("Naive", "Rested"),
    ("Stimulated", "Rechallenged"),
    ("Rested", "Rechallenged"),
]


def bh_qvalues(p_values: Sequence[float]) -> List[float]:
    indexed = sorted(enumerate(p_values), key=lambda item: item[1])
    n = len(indexed)
    q = [math.nan] * n
    running = 1.0
    for rank_from_end, (idx, p) in enumerate(reversed(indexed), start=1):
        rank = n - rank_from_end + 1
        val = min(running, p * n / rank)
        running = val
        q[idx] = min(1.0, val)
    return q


def cliffs_delta_from_u(u_stat: float, n_x: int, n_y: int) -> float:
    if n_x == 0 or n_y == 0:
        return math.nan
    return (2.0 * u_stat / (n_x * n_y)) - 1.0


def permutation_mean_diff(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> Tuple[float, float]:
    observed = float(np.mean(y) - np.mean(x))
    combined = np.concatenate([x, y])
    n_x = len(x)
    if n_x == 0 or len(y) == 0:
        return observed, math.nan
    more_extreme = 0
    for _ in range(n_perm):
        shuffled = rng.permutation(combined)
        diff = float(np.mean(shuffled[n_x:]) - np.mean(shuffled[:n_x]))
        if abs(diff) >= abs(observed):
            more_extreme += 1
    p_value = (more_extreme + 1) / (n_perm + 1)
    return observed, p_value


def bootstrap_mean_diff_ci(
    x: np.ndarray,
    y: np.ndarray,
    n_boot: int,
    rng: np.random.Generator,
) -> Tuple[float, float, float]:
    observed = float(np.mean(y) - np.mean(x))
    if len(x) == 0 or len(y) == 0:
        return observed, math.nan, math.nan
    diffs = np.empty(n_boot, dtype=float)
    for i in range(n_boot):
        bx = rng.choice(x, size=len(x), replace=True)
        by = rng.choice(y, size=len(y), replace=True)
        diffs[i] = np.mean(by) - np.mean(bx)
    lo, hi = np.percentile(diffs, [2.5, 97.5])
    return observed, float(lo), float(hi)


def per_fiber_values(matrix: SampleMatrix) -> np.ndarray:
    if matrix.selected_matrix.size == 0:
        return np.array([], dtype=float)
    return matrix.selected_matrix.mean(axis=1) * 100.0


def summarize_values(matrices: Sequence[SampleMatrix], out_dir: Path) -> Path:
    path = out_dir / "per_condition_msp_area_summary.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            delimiter="\t",
            fieldnames=[
                "locus",
                "condition",
                "strand_view",
                "selected_fibers",
                "mean_msp_covered_fiber_area_pct",
                "median_msp_covered_fiber_area_pct",
                "sd_pct",
                "iqr_pct",
            ],
        )
        writer.writeheader()
        for m in matrices:
            v = per_fiber_values(m)
            q25, q75 = np.percentile(v, [25, 75]) if len(v) else (math.nan, math.nan)
            writer.writerow(
                {
                    "locus": m.locus,
                    "condition": m.condition,
                    "strand_view": m.strand_view,
                    "selected_fibers": len(v),
                    "mean_msp_covered_fiber_area_pct": f"{np.mean(v):.4f}" if len(v) else "NA",
                    "median_msp_covered_fiber_area_pct": f"{np.median(v):.4f}" if len(v) else "NA",
                    "sd_pct": f"{np.std(v, ddof=1):.4f}" if len(v) > 1 else "NA",
                    "iqr_pct": f"{q75 - q25:.4f}" if len(v) else "NA",
                }
            )
    return path


def run_kruskal(locus_values: Dict[str, Dict[str, np.ndarray]], out_dir: Path) -> Path:
    rows = []
    pvals = []
    for locus, values_by_condition in locus_values.items():
        arrays = [values_by_condition[c] for c in CONDITION_ORDER if c in values_by_condition]
        if len(arrays) >= 2 and all(len(a) for a in arrays):
            h_stat, p_value = stats.kruskal(*arrays)
        else:
            h_stat, p_value = math.nan, math.nan
        rows.append({"locus": locus, "test": "Kruskal-Wallis", "h_statistic": h_stat, "p_value": p_value})
        pvals.append(p_value)
    qvals = bh_qvalues([p if not math.isnan(p) else 1.0 for p in pvals])
    path = out_dir / "overall_kruskal_wallis.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=["locus", "test", "h_statistic", "p_value", "q_value"])
        writer.writeheader()
        for row, q in zip(rows, qvals):
            writer.writerow(
                {
                    "locus": row["locus"],
                    "test": row["test"],
                    "h_statistic": f"{row['h_statistic']:.6g}" if not math.isnan(row["h_statistic"]) else "NA",
                    "p_value": f"{row['p_value']:.8g}" if not math.isnan(row["p_value"]) else "NA",
                    "q_value": f"{q:.8g}",
                }
            )
    return path


def run_pairwise_mannwhitney(locus_values: Dict[str, Dict[str, np.ndarray]], out_dir: Path) -> Tuple[Path, Dict[Tuple[str, str, str], float]]:
    rows = []
    pvals = []
    for locus, values_by_condition in locus_values.items():
        for i, cond_a in enumerate(CONDITION_ORDER):
            for cond_b in CONDITION_ORDER[i + 1 :]:
                x = values_by_condition.get(cond_a, np.array([], dtype=float))
                y = values_by_condition.get(cond_b, np.array([], dtype=float))
                if len(x) and len(y):
                    result = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
                    u_stat = float(result.statistic)
                    p_value = float(result.pvalue)
                    delta_ab = cliffs_delta_from_u(u_stat, len(x), len(y))
                    mean_diff = float(np.mean(y) - np.mean(x))
                    median_diff = float(np.median(y) - np.median(x))
                else:
                    u_stat = p_value = delta_ab = mean_diff = median_diff = math.nan
                rows.append(
                    {
                        "locus": locus,
                        "condition_a": cond_a,
                        "condition_b": cond_b,
                        "n_a": len(x),
                        "n_b": len(y),
                        "mean_b_minus_a_pct": mean_diff,
                        "median_b_minus_a_pct": median_diff,
                        "mannwhitney_u_for_a": u_stat,
                        "p_value": p_value,
                        "cliffs_delta_a_vs_b": delta_ab,
                    }
                )
                pvals.append(p_value)
    qvals = bh_qvalues([p if not math.isnan(p) else 1.0 for p in pvals])
    q_lookup: Dict[Tuple[str, str, str], float] = {}
    path = out_dir / "pairwise_mannwhitney.tsv"
    with path.open("w", newline="") as handle:
        fieldnames = [
            "locus",
            "condition_a",
            "condition_b",
            "n_a",
            "n_b",
            "mean_b_minus_a_pct",
            "median_b_minus_a_pct",
            "mannwhitney_u_for_a",
            "p_value",
            "q_value",
            "cliffs_delta_a_vs_b",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row, q in zip(rows, qvals):
            q_lookup[(row["locus"], row["condition_a"], row["condition_b"])] = q
            writer.writerow(
                {
                    **{k: row[k] for k in ["locus", "condition_a", "condition_b", "n_a", "n_b"]},
                    "mean_b_minus_a_pct": f"{row['mean_b_minus_a_pct']:.6g}" if not math.isnan(row["mean_b_minus_a_pct"]) else "NA",
                    "median_b_minus_a_pct": f"{row['median_b_minus_a_pct']:.6g}" if not math.isnan(row["median_b_minus_a_pct"]) else "NA",
                    "mannwhitney_u_for_a": f"{row['mannwhitney_u_for_a']:.6g}" if not math.isnan(row["mannwhitney_u_for_a"]) else "NA",
                    "p_value": f"{row['p_value']:.8g}" if not math.isnan(row["p_value"]) else "NA",
                    "q_value": f"{q:.8g}",
                    "cliffs_delta_a_vs_b": f"{row['cliffs_delta_a_vs_b']:.6g}" if not math.isnan(row["cliffs_delta_a_vs_b"]) else "NA",
                }
            )
    return path, q_lookup


def run_targeted_permutation_bootstrap(
    locus_values: Dict[str, Dict[str, np.ndarray]],
    out_dir: Path,
    n_perm: int,
    n_boot: int,
    random_seed: int,
) -> Path:
    rng = np.random.default_rng(random_seed)
    rows = []
    pvals = []
    for locus, values_by_condition in locus_values.items():
        for cond_a, cond_b in TARGETED_CONTRASTS:
            x = values_by_condition.get(cond_a, np.array([], dtype=float))
            y = values_by_condition.get(cond_b, np.array([], dtype=float))
            mean_diff_perm, p_perm = permutation_mean_diff(x, y, n_perm, rng)
            mean_diff_boot, ci_lo, ci_hi = bootstrap_mean_diff_ci(x, y, n_boot, rng)
            rows.append(
                {
                    "locus": locus,
                    "condition_a": cond_a,
                    "condition_b": cond_b,
                    "n_a": len(x),
                    "n_b": len(y),
                    "mean_b_minus_a_pct": mean_diff_perm,
                    "permutation_p_value": p_perm,
                    "bootstrap_mean_b_minus_a_pct": mean_diff_boot,
                    "bootstrap_ci95_low_pct": ci_lo,
                    "bootstrap_ci95_high_pct": ci_hi,
                    "n_permutations": n_perm,
                    "n_bootstraps": n_boot,
                }
            )
            pvals.append(p_perm)
    qvals = bh_qvalues([p if not math.isnan(p) else 1.0 for p in pvals])
    path = out_dir / "targeted_permutation_bootstrap.tsv"
    with path.open("w", newline="") as handle:
        fieldnames = [
            "locus",
            "condition_a",
            "condition_b",
            "n_a",
            "n_b",
            "mean_b_minus_a_pct",
            "permutation_p_value",
            "permutation_q_value",
            "bootstrap_ci95_low_pct",
            "bootstrap_ci95_high_pct",
            "n_permutations",
            "n_bootstraps",
        ]
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row, q in zip(rows, qvals):
            writer.writerow(
                {
                    "locus": row["locus"],
                    "condition_a": row["condition_a"],
                    "condition_b": row["condition_b"],
                    "n_a": row["n_a"],
                    "n_b": row["n_b"],
                    "mean_b_minus_a_pct": f"{row['mean_b_minus_a_pct']:.6g}",
                    "permutation_p_value": f"{row['permutation_p_value']:.8g}" if not math.isnan(row["permutation_p_value"]) else "NA",
                    "permutation_q_value": f"{q:.8g}",
                    "bootstrap_ci95_low_pct": f"{row['bootstrap_ci95_low_pct']:.6g}" if not math.isnan(row["bootstrap_ci95_low_pct"]) else "NA",
                    "bootstrap_ci95_high_pct": f"{row['bootstrap_ci95_high_pct']:.6g}" if not math.isnan(row["bootstrap_ci95_high_pct"]) else "NA",
                    "n_permutations": row["n_permutations"],
                    "n_bootstraps": row["n_bootstraps"],
                }
            )
    return path


def main() -> None:
    args = mm.parse_args()
    args.loci = [x.upper() for x in args.loci]
    args.no_cluster_loci = [x.upper() for x in args.no_cluster_loci]
    out_dir = cfg.out_dir("tables")

    matrices = mm.build_matrices(args)
    summary_tsv = summarize_values(matrices, out_dir)

    locus_values: Dict[str, Dict[str, np.ndarray]] = {}
    for m in matrices:
        if m.strand_view != "both":
            continue
        locus_values.setdefault(m.locus, {})[m.condition] = per_fiber_values(m)

    kruskal_tsv = run_kruskal(locus_values, out_dir)
    pairwise_tsv, _q_lookup = run_pairwise_mannwhitney(locus_values, out_dir)
    permutation_tsv = run_targeted_permutation_bootstrap(
        locus_values,
        out_dir,
        args.n_permutations,
        args.n_bootstraps,
        args.random_seed,
    )
    for path in (summary_tsv, kruskal_tsv, pairwise_tsv, permutation_tsv):
        print(path)


if __name__ == "__main__":
    main()
