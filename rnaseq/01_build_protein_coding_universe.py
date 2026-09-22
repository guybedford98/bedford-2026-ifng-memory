"""Restrict the max-TPM>=5 ISG universe to protein-coding genes.

Stage 00 rebuilt the IFNg-stimulated gene universe on a per-condition floor (max TPM
across the four training conditions >= 5), giving 1,122 genes: 57 positive memory, 44
refractory, 1,021 non-memory. The classification rule itself was untouched, so classes
never flip; only membership changes.

That universe still carries pseudogenes, lncRNAs, miRNAs and misc RNAs. Several of the
memory-class members are pseudogenes and snoRNAs, which are awkward to interpret and
hard to follow up. This script applies a second membership filter, protein-coding
biotype, on top of the max-TPM floor. Again the classification rule is unchanged and no
gene changes class.

Biotype source: GENCODE v44 gene records reduced to a symbol -> gene_type table.
Fourteen universe symbols have no v44 gene record; they were resolved once against HGNC
REST and NCBI Gene and the verdicts are frozen in data/symbol_resolution.tsv so this
script needs no network. Two of them (ARNTL2 -> BMAL2, DDX58 -> RIGI) are protein-coding
genes that were simply renamed and are recovered here; the rest are ncRNAs or
pseudogenes, except LOC100419583, which has no record anywhere and is excluded as
unresolved.

Only one universe symbol (GBP1P1) maps to more than one v44 gene record; both records
are non-coding, so the ambiguity does not affect the filter.
"""

import sys
import json
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("rnaseq_tables")

CONDS = ["NTC", "IFNg", "Trained", "Rechallenge"]
CLASSES = ["positive memory", "refractory", "non-memory"]

u = pd.read_csv(OUT / "ISG_universe_maxTPM5.csv")

gtf = pd.read_csv(cfg.path("reference", "gencode_v44_biotypes"), sep="\t", header=None,
                  names=["symbol", "biotype", "ensembl_gene_id"])
res = pd.read_csv(DATA / "symbol_resolution.tsv", sep="\t")

# GENCODE carries a handful of repeated symbols (PAR_Y copies and multi-locus names).
# Collapse to one record per symbol, preferring a protein_coding record if one exists
# so the filter is never made stricter by an arbitrary row order.
gtf["is_pc"] = gtf.biotype == "protein_coding"
gtf = (gtf.sort_values("is_pc", ascending=False)
          .drop_duplicates("symbol")
          .set_index("symbol"))

u["biotype"] = u.gene_id.map(gtf.biotype)
u["ensembl_gene_id"] = u.gene_id.map(gtf.ensembl_gene_id)
u["biotype_source"] = u.biotype.notna().map({True: "gencode_v44_gtf", False: None})
u["resolved_symbol"] = u.gene_id

fallback = res.set_index("query_symbol")
miss = u.biotype.isna()
u.loc[miss, "biotype"] = u.loc[miss, "gene_id"].map(fallback.biotype)
u.loc[miss, "biotype_source"] = u.loc[miss, "gene_id"].map(fallback.source)
u.loc[miss, "resolved_symbol"] = u.loc[miss, "gene_id"].map(fallback.resolved_symbol)

unannotated = u[u.biotype.isna()]
if len(unannotated):
    raise SystemExit(f"{len(unannotated)} symbols have no biotype call: "
                     f"{list(unannotated.gene_id)}")

u["protein_coding"] = u.biotype == "protein_coding"

print("=" * 72)
print("BIOTYPE ANNOTATION of the max-TPM>=5 universe")
print("=" * 72)
print(f"  universe in            : {len(u)}")
print(f"  matched in GENCODE v44 : {(u.biotype_source == 'gencode_v44_gtf').sum()}")
print(f"  resolved by fallback   : {(u.biotype_source != 'gencode_v44_gtf').sum()}")
print()
print(u.biotype.value_counts().to_string())

new = u[u.protein_coding].copy()
dropped = u[~u.protein_coding].copy()

print("\n" + "=" * 72)
print("PROTEIN-CODING UNIVERSE")
print("=" * 72)

rows = []
for cls in CLASSES:
    before = int((u["class"] == cls).sum())
    after = int((new["class"] == cls).sum())
    rows.append({
        "class": cls,
        "n_maxTPM5": before,
        "n_maxTPM5_protein_coding": after,
        "n_dropped": before - after,
        "pct_lost": round(100 * (before - after) / before, 1) if before else 0.0,
    })
summary = pd.DataFrame(rows)
summary = pd.concat([summary, pd.DataFrame([{
    "class": "TOTAL",
    "n_maxTPM5": len(u),
    "n_maxTPM5_protein_coding": len(new),
    "n_dropped": len(dropped),
    "pct_lost": round(100 * len(dropped) / len(u), 1),
}])], ignore_index=True)
print(summary.to_string(index=False))

for cls in ["positive memory", "refractory"]:
    lost = dropped[dropped["class"] == cls]
    print(f"\n--- {cls} genes dropped as non-coding ({len(lost)}) ---")
    if len(lost):
        print(lost[["gene_id", "biotype", "maxTPM", "memory_l2fc"]]
              .to_string(index=False, float_format=lambda x: f"{x:.2f}"))

print("\n--- divide-by-zero genes surviving both filters ---")
surv = new[new.div_by_zero]
print(f"  {len(surv)} genes; {int(surv.zero_naive_and_stim.sum())} zero in BOTH naive "
      f"and stimulated")
if len(surv):
    print(surv[["gene_id"] + CONDS + ["class"]]
          .to_string(index=False, float_format=lambda x: f"{x:.2f}"))

outcols = (["gene_id", "resolved_symbol", "ensembl_gene_id", "biotype",
            "biotype_source"] + CONDS +
           ["maxTPM", "sumTPM", "isg_l2fc", "memory_l2fc", "class",
            "div_by_zero", "zero_naive_and_stim"])

new[outcols].to_csv(OUT / "ISG_universe_maxTPM5_proteincoding.csv", index=False)
dropped[outcols].to_csv(OUT / "dropped_genes_non_protein_coding.csv", index=False)
summary.to_csv(OUT / "class_counts_summary_proteincoding.csv", index=False)
u[outcols].to_csv(OUT / "ISG_universe_maxTPM5_biotype_annotated.csv", index=False)

fname = {
    "positive memory": "Positive_memory_genes_maxTPM5_proteincoding.csv",
    "refractory": "Refractory_genes_maxTPM5_proteincoding.csv",
    "non-memory": "Non_memory_genes_maxTPM5_proteincoding.csv",
}
for cls, f in fname.items():
    sub = new[new["class"] == cls]
    sub[outcols].to_csv(OUT / f, index=False)
    print(f"\nwrote {f:<52} {len(sub):>5} genes")

# Same side calculation stage 00 reported: additionally require a computable ISG
# fold change (NTC > 0 AND IFNg > 0).
strict = new[(new.NTC > 0) & (new.IFNg > 0)]
print("\n" + "=" * 72)
print("SENSITIVITY (not applied above): also require NTC > 0 AND IFNg > 0")
print("=" * 72)
for cls in CLASSES:
    print(f"  {cls:<16} {int((new['class'] == cls).sum()):>5} -> "
          f"{int((strict['class'] == cls).sum()):>5}")
print(f"  {'TOTAL':<16} {len(new):>5} -> {len(strict):>5}")
strict[outcols].to_csv(
    OUT / "sensitivity_ISG_universe_maxTPM5_proteincoding_nonzero_NTC_IFNg.csv",
    index=False)

sens = {cls: int((strict["class"] == cls).sum()) for cls in CLASSES}
sens["TOTAL"] = len(strict)

with open(OUT / "summary_isg_universe_proteincoding.json", "w", encoding="utf-8") as fh:
    json.dump({
        "parent_universe": "ISG_universe_maxTPM5.csv",
        "max_tpm_floor": 5.0,
        "isg_l2fc_threshold": 0.5,
        "memory_l2fc_threshold": 0.5,
        "biotype_annotation": "GENCODE v44 gene records; 14 symbols resolved via "
                              "HGNC REST / NCBI Gene (data/symbol_resolution.tsv)",
        "universe_in": len(u),
        "universe_protein_coding": len(new),
        "dropped_non_protein_coding": len(dropped),
        "counts": summary.to_dict(orient="records"),
        "biotypes_dropped": dropped.biotype.value_counts().to_dict(),
        "renamed_symbols_recovered": ["ARNTL2 -> BMAL2", "DDX58 -> RIGI"],
        "unresolved_excluded": ["LOC100419583"],
        "div_by_zero_surviving": int(len(surv)),
        "zero_naive_and_stim_surviving": int(surv.zero_naive_and_stim.sum()),
        "sensitivity_nonzero_NTC_IFNg": sens,
    }, fh, indent=2)

print(f"\nwrote {'ISG_universe_maxTPM5_proteincoding.csv':<52} {len(new):>5} genes")
print(f"wrote {'dropped_genes_non_protein_coding.csv':<52} {len(dropped):>5} genes")
print(f"wrote {'ISG_universe_maxTPM5_biotype_annotated.csv':<52} {len(u):>5} genes")
print("wrote class_counts_summary_proteincoding.csv, summary_isg_universe_proteincoding.json")
