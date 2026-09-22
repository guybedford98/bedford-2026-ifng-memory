"""Rebuild the IFNg-stimulated gene (ISG) universe under a max-TPM floor.

The published universe (1,336 ISGs) was built from RSEM TPM with two rules:
  1. drop genes with SUMMED TPM < 5 across the four conditions
  2. keep genes with log2(IFNg / NTC) >= 0.5      [raw ratio, no pseudocount]
Memory classes were then called within that universe on
     log2(Rechallenge / IFNg)                     [raw ratio, no pseudocount]
     >= +0.5  positive memory        <= -0.5  refractory        else non-memory

The summed-TPM floor of 5 is a mean of 1.25 TPM per condition, so it admits a large
tail of genes that are never meaningfully expressed. This script replaces that floor
with a per-condition one (max TPM across the four conditions >= MAX_TPM_FLOOR) and
reclassifies.

Because max(TPM) >= 5 implies sum(TPM) >= 5, the new universe is a strict SUBSET of
the original 1,336. No genes outside the published list can enter.

The classification rule itself is unchanged, so a gene's class never flips; only
membership in the universe changes.
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

MAX_TPM_FLOOR = 5.0
ISG_L2FC = 0.5          # log2(IFNg / NTC) threshold defining an ISG
MEMORY_L2FC = 0.5       # |log2(Rechallenge / IFNg)| threshold for memory classes

DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("rnaseq_tables")

CONDS = ["NTC", "IFNg", "Trained", "Rechallenge"]


def log2_ratio(num, den):
    """Raw log2 ratio, reproducing the spreadsheet's divide-by-zero behaviour.

    A zero denominator with a positive numerator went to +inf (called positive
    memory); a zero numerator with a positive denominator went to -inf (called
    refractory). Both zero is left as NaN.
    """
    num = np.asarray(num, dtype=float)
    den = np.asarray(den, dtype=float)
    out = np.full(num.shape, np.nan)
    ok = den > 0
    out[ok] = np.log2(np.where(num[ok] > 0, num[ok], np.nan) / den[ok])
    out[(den == 0) & (num > 0)] = np.inf
    out[(den > 0) & (num == 0)] = -np.inf
    out[(den == 0) & (num == 0)] = np.nan
    return out


def classify(memory_l2fc):
    cls = np.full(len(memory_l2fc), "non-memory", dtype=object)
    cls[memory_l2fc >= MEMORY_L2FC] = "positive memory"
    cls[memory_l2fc <= -MEMORY_L2FC] = "refractory"
    return cls


u = pd.read_csv(DATA / "IFNg stimulated genes.csv")
published = {
    "positive memory": set(pd.read_csv(DATA / "Positive memory genes.csv").gene_id),
    "refractory": set(pd.read_csv(DATA / "negitive memory genes.csv").gene_id),
    "non-memory": set(pd.read_csv(DATA / "Non-Memory genes.csv").gene_id),
}

u["maxTPM"] = u[CONDS].max(axis=1)
u["sumTPM"] = u[CONDS].sum(axis=1)
u["isg_l2fc"] = log2_ratio(u.IFNg, u.NTC)
u["memory_l2fc"] = log2_ratio(u.Rechallenge, u.IFNg)
u["class"] = classify(u.memory_l2fc.values)

# a gene called an ISG despite being flat-zero in both naive and stimulated
u["zero_naive_and_stim"] = (u.NTC == 0) & (u.IFNg == 0)
u["div_by_zero"] = np.isinf(u.memory_l2fc) | (u.NTC == 0)

print("=" * 70)
print("REPRODUCTION CHECK against the published lists")
print("=" * 70)
ok = True
for cls, members in published.items():
    got = set(u.gene_id[u["class"] == cls])
    match = got == members
    ok &= match
    print(f"  {cls:<16} published {len(members):>5}   reproduced {len(got):>5}   "
          f"identical: {match}")
if not ok:
    raise SystemExit("Classification rule does not reproduce the published lists; stop.")
print("  -> rule reproduced exactly, safe to re-filter\n")

keep = u.maxTPM >= MAX_TPM_FLOOR
new = u[keep].copy()
dropped = u[~keep].copy()

print("=" * 70)
print(f"NEW UNIVERSE: max TPM across the four conditions >= {MAX_TPM_FLOOR:g}")
print("=" * 70)
print(f"  original universe : {len(u):>5}")
print(f"  dropped           : {len(dropped):>5}  ({100 * len(dropped) / len(u):.1f}%)")
print(f"  new universe      : {len(new):>5}\n")

rows = []
for cls in ["positive memory", "refractory", "non-memory"]:
    before = int((u["class"] == cls).sum())
    after = int((new["class"] == cls).sum())
    rows.append({
        "class": cls,
        "n_original": before,
        "n_maxTPM5": after,
        "n_dropped": before - after,
        "pct_lost": round(100 * (before - after) / before, 1) if before else 0.0,
    })
summary = pd.DataFrame(rows)
totals = {
    "class": "TOTAL",
    "n_original": len(u),
    "n_maxTPM5": len(new),
    "n_dropped": len(dropped),
    "pct_lost": round(100 * len(dropped) / len(u), 1),
}
summary = pd.concat([summary, pd.DataFrame([totals])], ignore_index=True)
print(summary.to_string(index=False))

print("\n--- divide-by-zero genes that SURVIVE the max-TPM floor ---")
surv = new[new.div_by_zero]
print(f"  {len(surv)} genes; of these {int(surv.zero_naive_and_stim.sum())} are zero in "
      f"BOTH naive and stimulated yet still called ISGs")
if len(surv):
    cols = ["gene_id"] + CONDS + ["class"]
    print(surv[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))

outcols = ["gene_id"] + CONDS + ["maxTPM", "sumTPM", "isg_l2fc", "memory_l2fc",
                                 "class", "div_by_zero", "zero_naive_and_stim"]

new[outcols].to_csv(OUT / "ISG_universe_maxTPM5.csv", index=False)
dropped[outcols].to_csv(OUT / "dropped_genes_maxTPM_under5.csv", index=False)
summary.to_csv(OUT / "class_counts_summary.csv", index=False)

fname = {
    "positive memory": "Positive_memory_genes_maxTPM5.csv",
    "refractory": "Refractory_genes_maxTPM5.csv",
    "non-memory": "Non_memory_genes_maxTPM5.csv",
}
for cls, f in fname.items():
    sub = new[new["class"] == cls]
    sub[outcols].to_csv(OUT / f, index=False)
    print(f"\nwrote {f:<40} {len(sub):>5} genes")

# The max-TPM floor does not remove the divide-by-zero genes: a gene can be 0 TPM
# in naive and/or stimulated and still peak above 5 in a later condition. Report
# what the counts would be if a computable fold change were also required
# (NTC > 0 AND IFNg > 0). This is a side calculation, not the deliverable above.
strict = new[(new.NTC > 0) & (new.IFNg > 0)]
print("\n" + "=" * 70)
print("SENSITIVITY (not applied to the files above): additionally require")
print("NTC > 0 AND IFNg > 0 so the ISG fold change is computable")
print("=" * 70)
for cls in ["positive memory", "refractory", "non-memory"]:
    print(f"  {cls:<16} {int((new['class'] == cls).sum()):>5} -> "
          f"{int((strict['class'] == cls).sum()):>5}")
print(f"  {'TOTAL':<16} {len(new):>5} -> {len(strict):>5}")
strict[outcols].to_csv(OUT / "sensitivity_ISG_universe_maxTPM5_nonzero_NTC_IFNg.csv",
                       index=False)

with open(OUT / "summary_isg_universe.json", "w", encoding="utf-8") as fh:
    json.dump({
        "max_tpm_floor": MAX_TPM_FLOOR,
        "isg_l2fc_threshold": ISG_L2FC,
        "memory_l2fc_threshold": MEMORY_L2FC,
        "original_universe": len(u),
        "new_universe": len(new),
        "dropped": len(dropped),
        "counts": summary.to_dict(orient="records"),
        "div_by_zero_surviving": int(len(surv)),
        "zero_naive_and_stim_surviving": int(surv.zero_naive_and_stim.sum()),
    }, fh, indent=2)

print(f"\nwrote {'ISG_universe_maxTPM5.csv':<40} {len(new):>5} genes")
print(f"wrote {'dropped_genes_maxTPM_under5.csv':<40} {len(dropped):>5} genes")
print("wrote class_counts_summary.csv, summary_isg_universe.json")
