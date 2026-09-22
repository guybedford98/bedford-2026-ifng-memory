"""Build the TNFa counterpart of the IFNg ISG universe.

Applies the same two membership filters stages 00 and 01 applied to the IFNg universe,
in the same order, to the published TNFa gene set:

    1. max TPM across the four training conditions >= 5
    2. protein-coding only, GENCODE v44 biotypes

and re-derives the memory classification with the same raw-ratio rules. Before filtering
anything it checks that those rules reproduce the published TNFa lists on disk, exactly
as stage 00 checked itself against the published IFNg lists.

Symbols with no v44 gene record are resolved from the frozen table in data/, so the
script re-runs with no network.
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

DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("rnaseq_tables")

UNIVERSE_CSV = DATA / "TNF stimulated genes.csv"
BIOTYPES_TSV = cfg.path("reference", "gencode_v44_biotypes")
RESOLUTION_TSV = DATA / "symbol_resolution_tnfa.tsv"

CONDS = ["Naive", "Stimulated", "Rested", "Rechallenged"]
MAXTPM_FLOOR = 5.0
ISG_CUT = 0.5       # log2(Stimulated / Naive) >= 0.5 defines the stimulated universe
MEM_CUT = 0.5       # log2(Rechallenged / Stimulated) >= +0.5 / <= -0.5

# The published TNFa class lists on disk, for the reproduction check
ONDISK = {
    "positive memory": "positive memory.csv",
    "refractory": "Negative memory.csv",
    "non-memory": "non memory TNF stimulated genes.csv",
}


def log2_ratio(num, den):
    """Raw ratio, no pseudocount. 0 denominators go to +/-inf, 0/0 to NaN.

    This is how the original spreadsheet behaved and is required to reproduce the
    published counts; the IFNg stages document the same choice.
    """
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.log2(num / den)


def classify(mem_l2fc):
    return np.where(mem_l2fc >= MEM_CUT, "positive memory",
                    np.where(mem_l2fc <= -MEM_CUT, "refractory", "non-memory"))


def main():
    print("=" * 78)
    print("TNFa universe, max TPM >= 5 and protein-coding")
    print("=" * 78)

    u = pd.read_csv(UNIVERSE_CSV).rename(columns={"Gene": "gene_id"})
    if list(u.columns) != ["gene_id"] + CONDS:
        raise SystemExit(f"unexpected columns in {UNIVERSE_CSV.name}: {list(u.columns)}")
    if u["gene_id"].duplicated().any():
        raise SystemExit("duplicate gene_id in the published TNFa universe")
    print(f"published TNFa universe: {len(u)} genes")

    u["maxTPM"] = u[CONDS].max(axis=1)
    u["sumTPM"] = u[CONDS].sum(axis=1)
    u["isg_l2fc"] = log2_ratio(u["Stimulated"], u["Naive"])
    u["memory_l2fc"] = log2_ratio(u["Rechallenged"], u["Stimulated"])
    u["class"] = classify(u["memory_l2fc"])
    u["div_by_zero"] = (u["Naive"] == 0) | (u["Stimulated"] == 0)
    u["zero_naive_and_stim"] = (u["Naive"] == 0) & (u["Stimulated"] == 0)

    # the universe should already satisfy the stimulated-gene call
    below = int((u["isg_l2fc"] < ISG_CUT).sum())
    print(f"  log2(Stim/Naive) >= {ISG_CUT}: holds for {len(u) - below} of {len(u)} "
          f"(min {np.nanmin(u['isg_l2fc']):.4f})")

    print("\n" + "-" * 78)
    print("reproduction check against the published TNFa lists")
    print("-" * 78)
    ondisk = {k: set(pd.read_csv(DATA / f)["Gene"]) for k, f in ONDISK.items()}
    overlaps = {}
    keys = list(ondisk)
    for i, a in enumerate(keys):
        for b in keys[i + 1:]:
            both = ondisk[a] & ondisk[b]
            if both:
                overlaps[(a, b)] = sorted(both)
    for (a, b), genes in overlaps.items():
        print(f"  NOTE: {len(genes)} genes appear in BOTH the on-disk '{a}' and '{b}' "
              f"lists: {', '.join(genes)}")

    mismatches = {}
    for cls, on in ondisk.items():
        calc = set(u.loc[u["class"] == cls, "gene_id"])
        missing, extra = sorted(on - calc), sorted(calc - on)
        mismatches[cls] = (missing, extra)
        flag = "exact" if not missing and not extra else f"-{len(missing)} +{len(extra)}"
        print(f"  {cls:16s} rule {len(calc):5d}   on disk {len(on):5d}   {flag}")

    unexplained = []
    for cls, (missing, extra) in mismatches.items():
        for g in missing + extra:
            row = u.loc[u["gene_id"] == g].iloc[0]
            # a 0/0 ratio is unclassifiable and is what put a gene in two lists at once
            if row["zero_naive_and_stim"]:
                continue
            # a value sitting on the threshold to 4 dp is a rounding edge, not a rule change
            if abs(abs(row["memory_l2fc"]) - MEM_CUT) < 1e-3:
                continue
            unexplained.append((cls, g, float(row["memory_l2fc"])))
    if unexplained:
        print("\n  UNEXPLAINED disagreements (not 0/0, not on the threshold):")
        for cls, g, v in unexplained:
            print(f"    {g}  on-disk/rule split for {cls}, memory_l2fc = {v:+.4f}")
        raise SystemExit("the rule does not reproduce the published TNFa lists; stopping")
    print("\n  every disagreement is either a 0/0 ratio or a value sitting exactly on the "
          "+/-0.5 threshold.\n  The rule reproduces the published TNFa classification.")

    print("\n" + "-" * 78)
    print(f"filter 1: max TPM across the four conditions >= {MAXTPM_FLOOR:g}")
    print("-" * 78)
    u["pass_maxtpm"] = u["maxTPM"] >= MAXTPM_FLOOR
    step1 = u[u["pass_maxtpm"]].copy()
    print(f"  {len(u)} -> {len(step1)} genes ({len(u) - len(step1)} dropped)")

    print("\n" + "-" * 78)
    print("filter 2: protein-coding only, GENCODE v44")
    print("-" * 78)
    bt = pd.read_csv(BIOTYPES_TSV, sep="\t", header=None,
                     names=["symbol", "biotype", "ensembl_gene_id"])
    print(f"  biotype table: {len(bt)} v44 gene records")
    # A symbol can map to several v44 gene records. If any of them is protein_coding, the
    # symbol is protein-coding: the RSEM quantification is per symbol, and the coding locus
    # is the one carrying the expression. Only if no record is coding does the symbol fall
    # out, and then the biotype reported is the first one. The IFNg stage collapsed every
    # multi-record symbol to "AMBIGUOUS" and dropped it, which was safe there only because
    # its single ambiguous symbol (GBP1P1) is non-coding in both records. It is not safe
    # here: SPATA13 has a protein_coding record and an lncRNA one.
    def collapse(s):
        if (s == "protein_coding").any():
            return "protein_coding"
        return s.iloc[0] if s.nunique() == 1 else f"AMBIGUOUS:{'/'.join(sorted(set(s)))}"

    grp = bt.groupby("symbol")["biotype"].agg(collapse)
    multi = bt["symbol"].value_counts()
    multi = set(multi[multi > 1].index)
    ens = bt.drop_duplicates("symbol").set_index("symbol")["ensembl_gene_id"]

    res = pd.read_csv(RESOLUTION_TSV, sep="\t")
    res_bt = res.set_index("query_symbol")

    rows = []
    for g in step1["gene_id"]:
        if g in grp.index:
            rows.append((g, g, grp[g], ens.get(g, ""), "gencode_v44_gtf"))
        elif g in res_bt.index:
            r = res_bt.loc[g]
            rows.append((g, r["resolved_symbol"] if pd.notna(r["resolved_symbol"]) else g,
                         r["biotype"], "", r["source"]))
        else:
            raise SystemExit(
                f"'{g}' has no v44 gene record and no entry in {RESOLUTION_TSV.name}. "
                "Resolve it and add a frozen row before re-running.")
    ann = pd.DataFrame(rows, columns=["gene_id", "resolved_symbol", "biotype",
                                      "ensembl_gene_id", "biotype_source"])
    step1 = step1.merge(ann, on="gene_id", how="left")

    n_res = int((step1["biotype_source"] != "gencode_v44_gtf").sum())
    print(f"  {len(step1) - n_res} symbols annotated directly, {n_res} from the frozen table")
    dup = step1[step1["gene_id"].isin(multi)]
    if len(dup):
        print(f"  {len(dup)} symbols map to more than one v44 gene record:")
        for _, r in dup.iterrows():
            print(f"    {r['gene_id']:14s} -> {r['biotype']}")

    step1["protein_coding"] = step1["biotype"] == "protein_coding"
    pc = step1[step1["protein_coding"]].copy()
    print(f"  {len(step1)} -> {len(pc)} genes ({len(step1) - len(pc)} dropped)")

    recovered = step1[(step1["biotype_source"] != "gencode_v44_gtf") & step1["protein_coding"]]
    if len(recovered):
        print("\n  protein-coding genes recovered by the frozen resolution table "
              "(symbol matching alone would have dropped these):")
        for _, r in recovered.iterrows():
            print(f"    {r['gene_id']} -> {r['resolved_symbol']}")

    dropped_bt = step1[~step1["protein_coding"]]["biotype"].value_counts()
    print("\n  biotypes removed:")
    for b, n in dropped_bt.items():
        print(f"    {n:4d}  {b}")

    print("\n" + "-" * 78)
    print("class counts through both filters")
    print("-" * 78)
    order = ["positive memory", "refractory", "non-memory"]
    tab = pd.DataFrame({
        "published": u["class"].value_counts().reindex(order),
        "maxTPM>=5": step1["class"].value_counts().reindex(order),
        "+ protein-coding": pc["class"].value_counts().reindex(order),
    }).fillna(0).astype(int)
    tab.loc["TOTAL"] = tab.sum()
    tab["% lost"] = ((1 - tab["+ protein-coding"] / tab["published"]) * 100).round(1)
    print(tab.to_string())

    for cls in ("positive memory", "refractory"):
        lost = sorted(set(u.loc[u["class"] == cls, "gene_id"]) - set(pc.loc[pc["class"] == cls, "gene_id"]))
        print(f"\n  {cls} dropped ({len(lost)}): {', '.join(lost)}")

    cols = (["gene_id", "resolved_symbol", "ensembl_gene_id", "biotype", "biotype_source"]
            + CONDS + ["maxTPM", "sumTPM", "isg_l2fc", "memory_l2fc", "class",
                       "div_by_zero", "zero_naive_and_stim"])
    pc[cols].to_csv(OUT / "TNFa_universe_maxTPM5_proteincoding.csv", index=False)
    step1[cols + ["protein_coding"]].to_csv(
        OUT / "TNFa_universe_maxTPM5_biotype_annotated.csv", index=False)
    for cls, name in (("positive memory", "Positive_memory"), ("refractory", "Refractory"),
                      ("non-memory", "Non_memory")):
        pc[pc["class"] == cls][cols].to_csv(
            OUT / f"{name}_genes_TNFa_maxTPM5_proteincoding.csv", index=False)
    dropped = u[~u["gene_id"].isin(pc["gene_id"])].copy()
    dropped["reason"] = np.where(~dropped["pass_maxtpm"], "maxTPM<5", "not protein-coding")
    dropped[["gene_id"] + CONDS + ["maxTPM", "class", "reason"]].to_csv(
        OUT / "dropped_genes_TNFa.csv", index=False)
    tab.to_csv(OUT / "class_counts_summary_TNFa.csv")

    summary = {
        "source": str(UNIVERSE_CSV),
        "biotype_reference": str(BIOTYPES_TSV),
        "n_published": int(len(u)),
        "n_after_maxtpm": int(len(step1)),
        "n_protein_coding": int(len(pc)),
        "maxtpm_floor": MAXTPM_FLOOR,
        "isg_cut": ISG_CUT,
        "memory_cut": MEM_CUT,
        "counts": {k: {c: int(tab.loc[c, k]) for c in tab.index} for k in
                   ("published", "maxTPM>=5", "+ protein-coding")},
        "on_disk_list_overlaps": {f"{a} & {b}": g for (a, b), g in overlaps.items()},
        "protein_coding_recovered_by_resolution": recovered["gene_id"].tolist(),
    }
    (OUT / "summary_tnfa_universe.json").write_text(json.dumps(summary, indent=2),
                                                    encoding="utf-8")
    print(f"\nwritten to {OUT}")
    print("done.")


if __name__ == "__main__":
    main()
