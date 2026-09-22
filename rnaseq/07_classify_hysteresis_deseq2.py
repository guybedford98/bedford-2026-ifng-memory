"""A DESeq2 test for hysteresis, run alongside the TPM-ratio memory classes.

Every ISG definition in stages 00-01 is a raw TPM ratio with no statistical test: the
universe is log2(IFNg / NTC) >= 0.5, and the memory classes are log2(Rechallenge /
IFNg) >= +0.5 or <= -0.5 on mean TPM. This is the DESeq2 cross-reference.

What this adds, and what it does NOT change
-------------------------------------------
Universe   : the 1,063 protein-coding ISGs, unchanged.
Hysteresis : DESeq2 Stimulated vs Rechallenged, padj < 0.05, either direction. NEW.
Memory     : TPM-ratio log2(Rech / IFNg) >= +0.5. UNCHANGED -- it arrives in the
             `class` column of the input universe and is never rewritten here.
Refractory : TPM-ratio log2(Rech / IFNg) <= -0.5. UNCHANGED, likewise.

The two definitions are deliberately kept SEPARATE. Hysteresis is a test ("did anything
change between the first and second stimulation"); memory and refractory stay a
magnitude call on mean TPM. This script adds a `hysteretic` flag to every gene.

It also computes a fully DESeq2-based `deseq2_class` (padj < 0.05 plus a +/-0.5 cut on
the DESeq2 log2FoldChange) purely as a CROSS-REFERENCE, so the two rules can be
compared. That column is not the manuscript classification.

Which log2FC the 0.5 cut applies to -- a trap
---------------------------------------------
The memory rule thresholds the TPM-ratio log2FC (`memory_l2fc`), NOT the DESeq2
`log2FoldChange`. DESeq2 shrinks its estimate, so the two are close but not equal: five
memory genes sit at 0.51-0.55 on the TPM ratio and 0.47-0.49 in DESeq2 (RTKN, CD40,
CRISPLD2, PARP8, LYPD1). Applying the cut to the wrong column silently drops them. That
is exactly the difference between the 50 and the 41 in the comparison table below.

Sign convention
---------------
The DESeq2 contrast is Naive+IFNg vs Trained+IFNg, i.e. Stimulated vs Rechallenged, and
POSITIVE log2FoldChange means higher in Rechallenged. The script asserts this against
IDO1 and CXCL11 before classifying anything, so a re-run against a re-exported DESeq2
table cannot silently flip the up and down lists.
"""

import sys
import json
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

OUT = cfg.out_dir("rnaseq_tables")
UNIVERSE = OUT / "ISG_universe_maxTPM5_proteincoding.csv"
DESEQ = cfg.path("rnaseq", "deseq2_stim_vs_rechallenge")

PADJ = 0.05
LFC = 0.5

CLASS_MEM = "memory"
CLASS_REF = "refractory"
CLASS_SUB = "sub-threshold hysteresis"
CLASS_NONE = "no hysteresis"

uni = pd.read_csv(UNIVERSE)
de = pd.read_csv(DESEQ)
de = de.rename(columns={de.columns[0]: "gene"})
de = de[["gene", "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj"]]

print("universe            : %d" % len(uni))
print("DESeq2 rows         : %d" % len(de))

if de["gene"].duplicated().any():
    raise SystemExit("DESeq2 table has duplicate gene symbols; resolve before classifying")

m = uni.merge(de, left_on="gene_id", right_on="gene", how="left").drop(columns="gene")
n_hit = int(m["log2FoldChange"].notna().sum())
n_miss = len(m) - n_hit
n_padj_na = int((m["log2FoldChange"].notna() & m["padj"].isna()).sum())
print("matched to DESeq2   : %d   unmatched: %d" % (n_hit, n_miss))
print("padj NA (independent filtering): %d" % n_padj_na)
if n_miss:
    m.loc[m["log2FoldChange"].isna(), ["gene_id"]].to_csv(
        OUT / "unmatched_in_deseq2.csv", index=False)

# Both anchors are strongly higher on rechallenge in TPM; if the contrast ran the other
# way round their log2FoldChange would be negative.
anchors = m.set_index("gene_id").loc[["IDO1", "CXCL11"]]
for g, row in anchors.iterrows():
    if not (row["Rechallenge"] > row["IFNg"] and row["log2FoldChange"] > 0):
        raise SystemExit(
            "sign check failed on %s: TPM %.1f -> %.1f but log2FoldChange %+.3f. "
            "The contrast direction is not Stimulated -> Rechallenged; do not trust "
            "the class lists." % (g, row["IFNg"], row["Rechallenge"],
                                  row["log2FoldChange"]))
print("sign check          : OK (positive log2FoldChange = higher in Rechallenged)")

sig = m["padj"].notna() & (m["padj"] < PADJ)
m["hysteretic"] = sig
n_hyst = int(sig.sum())
n_up = int((sig & (m["log2FoldChange"] > 0)).sum())

# Cross-reference only: what a fully DESeq2-based rule would call. The incoming `class`
# column is left exactly as it arrived.
m["deseq2_class"] = CLASS_NONE
m.loc[sig & (m["log2FoldChange"] >= LFC), "deseq2_class"] = CLASS_MEM
m.loc[sig & (m["log2FoldChange"] <= -LFC), "deseq2_class"] = CLASS_REF
m.loc[sig & (m["log2FoldChange"].abs() < LFC), "deseq2_class"] = CLASS_SUB

order = [CLASS_MEM, CLASS_REF, CLASS_SUB, CLASS_NONE]
counts = {c: int((m["deseq2_class"] == c).sum()) for c in order}

print()
print("--- DESeq2 Stimulated vs Rechallenged on %d protein-coding ISGs" % len(m))
print("HYSTERESIS, padj < %.2f            : %d  (%.1f%%)   up %d / down %d"
      % (PADJ, n_hyst, 100.0 * n_hyst / len(m), n_up, n_hyst - n_up))
print()
print("manuscript classes, carried through unchanged from the input universe:")
for c, n in m["class"].value_counts().items():
    print("  %-26s: %d   of which hysteretic: %d"
          % (c, n, int((sig & (m["class"] == c)).sum())))
print()
print("cross-reference only -- padj plus a +/-%.1f cut on the DESeq2 log2FoldChange:" % LFC)
for c in order:
    print("  %-26s: %d" % (c, counts[c]))

xt = pd.crosstab(m["class"], m["deseq2_class"]).reindex(columns=order, fill_value=0)
print()
print("--- TPM-ratio class (rows) vs DESeq2 cross-reference class (columns)")
print(xt.to_string())
xt.to_csv(OUT / "concordance_tpmratio_vs_deseq2.csv")

# Same padj, but the +/-0.5 cut applied to each of the two log2FC columns in turn.
rules = pd.DataFrame([
    {"rule": "TPM-ratio log2FC only (the manuscript rule)",
     "memory": int((m["memory_l2fc"] >= LFC).sum()),
     "refractory": int((m["memory_l2fc"] <= -LFC).sum())},
    {"rule": "padj + TPM-ratio log2FC",
     "memory": int((sig & (m["memory_l2fc"] >= LFC)).sum()),
     "refractory": int((sig & (m["memory_l2fc"] <= -LFC)).sum())},
    {"rule": "padj + DESeq2 log2FC",
     "memory": counts[CLASS_MEM], "refractory": counts[CLASS_REF]},
    {"rule": "padj + direction only, no magnitude cut",
     "memory": n_up, "refractory": n_hyst - n_up},
])
rules.to_csv(OUT / "rule_comparison.csv", index=False)
print()
print("--- what each rule would give")
print(rules.to_string(index=False))

rows = []
for padj_thr in (0.05, 0.01, 0.001):
    s = m["padj"].notna() & (m["padj"] < padj_thr)
    for lfc_thr in (0.25, 0.5, 1.0):
        rows.append({
            "padj_threshold": padj_thr,
            "lfc_threshold": lfc_thr,
            "hysteretic": int(s.sum()),
            "memory_deseq2_lfc": int((s & (m["log2FoldChange"] >= lfc_thr)).sum()),
            "refractory_deseq2_lfc": int((s & (m["log2FoldChange"] <= -lfc_thr)).sum()),
        })
sens = pd.DataFrame(rows)
sens.to_csv(OUT / "sensitivity_thresholds.csv", index=False)
print()
print("--- threshold sensitivity (cross-reference rule)")
print(sens.to_string(index=False))

outcols = ["gene_id", "resolved_symbol", "ensembl_gene_id", "biotype",
           "NTC", "IFNg", "Trained", "Rechallenge", "maxTPM",
           "isg_l2fc", "memory_l2fc", "class",
           "baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj",
           "hysteretic", "deseq2_class"]
outcols = [c for c in outcols if c in m.columns]
m = m.sort_values("padj", na_position="last")
m[outcols].to_csv(OUT / "ISG_universe_proteincoding_deseq2_classified.csv", index=False)

# The manuscript classes, now carrying the hysteresis flag.
for cls, fname in [("positive memory", "Memory_genes_with_hysteresis_flag.csv"),
                   ("refractory", "Refractory_genes_with_hysteresis_flag.csv")]:
    sub = m[m["class"] == cls].sort_values("memory_l2fc",
                                           ascending=(cls == "refractory"))
    sub[outcols].to_csv(OUT / fname, index=False)

# The cross-reference lists.
for cls, fname in [(CLASS_MEM, "xref_deseq2_memory.csv"),
                   (CLASS_REF, "xref_deseq2_refractory.csv"),
                   (CLASS_SUB, "xref_deseq2_subthreshold_hysteresis.csv"),
                   (CLASS_NONE, "xref_deseq2_no_hysteresis.csv")]:
    sub = m[m["deseq2_class"] == cls].sort_values(
        "log2FoldChange", ascending=(cls == CLASS_REF))
    sub[outcols].to_csv(OUT / fname, index=False)

summary_tbl = pd.DataFrame(
    [{"quantity": "universe", "n": len(m)},
     {"quantity": "hysteretic (DESeq2 padj < %.2f)" % PADJ, "n": n_hyst},
     {"quantity": "  of which up on rechallenge", "n": n_up},
     {"quantity": "  of which down on rechallenge", "n": n_hyst - n_up}]
    + [{"quantity": "TPM-ratio %s" % k, "n": int(v)}
       for k, v in m["class"].value_counts().items()]
    + [{"quantity": "xref DESeq2 %s" % c, "n": counts[c]} for c in order])
summary_tbl.to_csv(OUT / "class_counts_summary_deseq2.csv", index=False)

mem = m[m["class"] == "positive memory"]
top = m[m["hysteretic"] & (m["log2FoldChange"] > 0)].nlargest(15, "log2FoldChange")
with open(OUT / "summary_hysteresis.json", "w", encoding="utf-8") as fh:
    json.dump({
        "universe": "ISG_universe_maxTPM5_proteincoding.csv",
        "deseq2_table": str(DESEQ),
        "contrast": "Stimulated vs Rechallenged; positive log2FoldChange = higher in Rechallenged",
        "padj_threshold": PADJ,
        "lfc_threshold_crossref_only": LFC,
        "n_universe": len(m),
        "n_matched_in_deseq2": n_hit,
        "n_unmatched": n_miss,
        "n_padj_na": n_padj_na,
        "hysteresis": {"n": n_hyst, "up": n_up, "down": n_hyst - n_up,
                       "pct_of_universe": round(100.0 * n_hyst / len(m), 1)},
        "manuscript_classes_unchanged": {
            k: int(v) for k, v in m["class"].value_counts().items()},
        "manuscript_memory_hysteretic": int((mem["padj"] < PADJ).sum()),
        "manuscript_memory_not_hysteretic": mem.loc[
            mem["padj"] >= PADJ, "gene_id"].tolist(),
        "manuscript_refractory_not_hysteretic": m.loc[
            (m["class"] == "refractory") & (m["padj"] >= PADJ), "gene_id"].tolist(),
        "crossref_deseq2_counts": counts,
        "rule_comparison": rules.to_dict("records"),
        "anchors": {g: {"IFNg_TPM": float(anchors.loc[g, "IFNg"]),
                        "Rechallenge_TPM": float(anchors.loc[g, "Rechallenge"]),
                        "log2FoldChange": float(anchors.loc[g, "log2FoldChange"]),
                        "padj": float(anchors.loc[g, "padj"])} for g in anchors.index},
        "top15_up_by_log2fc": top[["gene_id", "IFNg", "Rechallenge",
                                   "log2FoldChange", "padj"]].to_dict("records"),
    }, fh, indent=2)

print()
print("wrote:")
for p in sorted(OUT.glob("*")):
    print("  ", p.name)
