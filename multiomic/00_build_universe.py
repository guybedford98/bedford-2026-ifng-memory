#!/usr/bin/env python3
"""Stage 00 -- build the gene universe.

Start from the protein-coding ISG universe written by rnaseq/ and apply the
zero-in-both-naive-and-stim exclusion. Nothing else. This is the only place the universe is
defined; every later stage reads universe.tsv and must never re-filter.

Output: data/universe.tsv, data/universe_excluded.tsv, data/universe_summary.json
"""
import sys, json
from pathlib import Path
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

UNIVERSE_CSV = cfg.path("multiomic", "isg_universe")
EXTRA_COLS = ["resolved_symbol", "ensembl_gene_id", "biotype", "biotype_source"]
OUT = Path(__file__).resolve().parent / "data"
OUT.mkdir(parents=True, exist_ok=True)

# The four genes this removes are MIR1307, PAPPA-AS1, RPL37P18, SNHG9: 0.00 TPM in BOTH
# Naive and IFNg, so the ISG call is log2(0/0) and the memory call is log2(x/0) = +inf.
# Do NOT widen to "any zero in naive or stim" -- that also removes LGALS17A
# (0 -> 9.21 -> 0.06 -> 15.86), a textbook memory gene whose naive zero is the whole point.
EXCLUDE_ZERO_IN_BOTH = True

u = pd.read_csv(UNIVERSE_CSV).set_index("gene_id")
u.index = u.index.astype(str)
print(f"ISG universe as published: {len(u)}")
for c in ["positive memory", "refractory", "non-memory"]:
    print(f"  {c:<16} {int((u['class'] == c).sum())}")

if EXCLUDE_ZERO_IN_BOTH:
    drop = u[u.zero_naive_and_stim]
    keep = u[~u.zero_naive_and_stim]
else:
    drop = u.iloc[0:0]
    keep = u

print(f"\nexcluded (zero in both Naive and IFNg): {len(drop)}")
if len(drop):
    print(drop[["NTC", "IFNg", "Trained", "Rechallenge", "class"]].round(2).to_string())

keep = keep.rename(columns={"NTC": "Naive", "IFNg": "Stim",
                            "Trained": "Rested", "Rechallenge": "Rechallenged"})
keep = keep.rename(columns={"class": "memory_class"})
keep["memory_class"] = keep["memory_class"].map(
    {"positive memory": "memory", "refractory": "refractory", "non-memory": "non_memory"})

cols = ["Naive", "Stim", "Rested", "Rechallenged", "maxTPM", "sumTPM",
        "isg_l2fc", "memory_l2fc", "memory_class", "div_by_zero"]
cols += [c for c in EXTRA_COLS if c in keep.columns]      # biotype provenance, when present
keep[cols].to_csv(OUT / "universe.tsv", sep="\t", index_label="gene")
drop.to_csv(OUT / "universe_excluded.tsv", index_label="gene")

summary = {
    "source": str(UNIVERSE_CSV),
    "exclusion_rule": "zero TPM in both Naive and Stim" if EXCLUDE_ZERO_IN_BOTH else "none",
    "n_published": int(len(u)),
    "n_excluded": int(len(drop)),
    "excluded_genes": sorted(drop.index.tolist()),
    "n_universe": int(len(keep)),
    "by_class": {c: int((keep.memory_class == c).sum())
                 for c in ["memory", "refractory", "non_memory"]},
}
(OUT / "universe_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

print(f"\nuniverse: {summary['n_universe']}  " +
      "  ".join(f"{k} {v}" for k, v in summary["by_class"].items()))
print(f"wrote {OUT/'universe.tsv'}")

if summary["n_universe"] != 1063:
    print(f"\nNOTE: got {summary['n_universe']}, not the expected 1063 protein-coding ISGs. "
          "Confirm before proceeding.")
