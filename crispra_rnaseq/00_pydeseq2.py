#!/usr/bin/env python3
"""CRISPRa CXCL11 synthetic-memory RNA-seq, differential expression (day 10 & day 21).

Plasmidsaurus featureCounts matrices (raw _count) -> DESeq2 (pydeseq2), one fit per
timepoint on a single 4-level group factor. Writes the per-contrast tables and the
results JSON that the panel scripts read.

2x2 design per timepoint (guide x IFNg), n=3:
  samples 1-3  NT_rest      4-6  NT_IFNg
  samples 7-9  sg1_rest    10-12 sg1_IFNg
Normalization + Wald tests via DESeq2 (median-of-ratios; NB GLM). Low-count
filter: total counts >= 10 (matches paper). DEG = padj<0.05 & |log2FC|>1.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))

import csv, json, warnings
import numpy as np
import pandas as pd
warnings.simplefilter("ignore")
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

import config as cfg

MANIFEST = cfg.path("crispra_rnaseq", "sample_sheet")
DATASETS = {"day10": cfg.path("crispra_rnaseq", "counts_day10"),
            "day21": cfg.path("crispra_rnaseq", "counts_day21")}
PROJ = {"day10": "CF63KH", "day21": "GVSYGD"}
DAY_LABEL = {"day10": "Day 10", "day21": "Day 21"}

OUT = cfg.out_dir("crispra_rnaseq")
TAB = cfg.out_dir("crispra_rnaseq", "tables")

COND_ORDER = ["NT_rest", "NT_IFNg", "sg1_rest", "sg1_IFNg"]
PADJ, LFC_T = 0.05, 1.0

TARGET = ["CXCL11", "IDO1"]
CONTRASTS = [
    ("IFNg_response_NT", "NT_IFNg", "NT_rest", "IFNγ response (NT+IFNγ vs NT rest)"),
    ("priming_at_rest", "sg1_rest", "NT_rest", "CRISPRa priming at rest (sg1 rest vs NT rest)"),
    ("memory_rechallenge", "sg1_IFNg", "NT_IFNg", "Synthetic memory (sg1+IFNγ vs NT+IFNγ)"),
]


def load_manifest():
    g_map = {"Synthego NC #1": "NT", "CXCL11 sg1": "sg1"}
    a_map = {"rested": "rest", "+IFNg": "IFNg"}
    out = {"day10": {}, "day21": {}}
    with open(MANIFEST, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            tp = "day10" if "Day 10" in r["order"] else "day21"
            out[tp][int(r["sample_label"].split("_")[1])] = f"{g_map[r['guide']]}_{a_map[r['arm']]}"
    return out


def load_matrix(path, proj):
    ids, names, biotype, count_cols = [], [], [], {}
    with open(path, newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        header = next(rdr)
        for j, h in enumerate(header):
            if h.endswith("_count"):
                count_cols[int(h[len(proj) + 1:-len("_count")])] = j
        order = sorted(count_cols)
        col_idx = [count_cols[n] for n in order]
        rows = []
        for row in rdr:
            ids.append(row[0]); names.append(row[1]); biotype.append(row[2])
            rows.append([float(row[j]) for j in col_idx])
    return np.array(ids), np.array(names), np.array(biotype), np.asarray(rows, float), order


def size_factors(counts):
    with np.errstate(divide="ignore"):
        logc = np.log(counts)
    mask = np.all(counts > 0, axis=1)
    ref = logc[mask].mean(axis=1, keepdims=True)
    return np.exp(np.median(logc[mask] - ref, axis=0))


def run_deseq(tp, ids, names, counts, cond, order):
    keep = counts.sum(1) >= 10
    gid = ids[keep]
    id2name = dict(zip(ids, names))
    cnt = np.rint(counts[keep]).astype(int).T          # samples x genes
    samples = [f"s{n}" for n in order]
    cdf = pd.DataFrame(cnt, index=samples, columns=gid)
    meta = pd.DataFrame({"group": list(cond)}, index=samples)
    try:
        dds = DeseqDataSet(counts=cdf, metadata=meta, design="~group")
    except TypeError:
        dds = DeseqDataSet(counts=cdf, metadata=meta, design_factors="group")
    dds.deseq2()
    pos = {g: i for i, g in enumerate(ids)}
    de = {}
    for key, a, b, title in CONTRASTS:
        ds = DeseqStats(dds, contrast=["group", a, b])
        ds.summary()
        r = ds.results_df
        rr = r.copy(); rr.insert(0, "gene_name", [id2name.get(g, "") for g in rr.index])
        rr.sort_values("padj").to_csv(TAB / f"{tp}_{key}.csv")
        lfc = np.full(len(ids), np.nan); padj = np.full(len(ids), np.nan)
        pval = np.full(len(ids), np.nan); base = np.full(len(ids), np.nan)
        for g, l, pa, pv, bm in zip(r.index.to_numpy(), r["log2FoldChange"].to_numpy(),
                                    r["padj"].to_numpy(), r["pvalue"].to_numpy(),
                                    r["baseMean"].to_numpy()):
            i = pos.get(g)
            if i is not None:
                lfc[i] = l; padj[i] = pa; pval[i] = pv; base[i] = bm
        de[key] = dict(lfc=lfc, padj=padj, pval=pval, base=base)
    return de, keep


def analyze():
    man = load_manifest()
    res = {"datasets": {}, "thresholds": {"padj": PADJ, "lfc": LFC_T},
           "contrasts": {c[0]: c[3] for c in CONTRASTS}}
    for tp, path in DATASETS.items():
        ids, names, biotype, counts, order = load_matrix(path, PROJ[tp])
        cond = np.array([man[tp][n] for n in order])
        cols = {c: np.where(cond == c)[0] for c in COND_ORDER}
        sf = size_factors(counts)
        norm = counts / sf
        name2i = {}
        for i, nm in enumerate(names):
            name2i.setdefault(nm, i)
        print(f"[{tp}] DESeq2 ...", flush=True)
        de, keep = run_deseq(tp, ids, names, counts, cond, order)

        ds = {"n_genes": int(len(ids)), "n_expressed": int(keep.sum()), "contrasts": {}}
        for key, a, b, title in CONTRASTS:
            lfc, padj = de[key]["lfc"], de[key]["padj"]
            ok = ~np.isnan(padj)
            sig = ok & (padj < PADJ) & (np.abs(lfc) > LFC_T)
            sig05 = ok & (padj < PADJ) & (np.abs(lfc) > 0.5)
            order_idx = np.where(ok)[0]
            order_idx = order_idx[np.argsort(padj[order_idx])]
            top = [{"name": str(names[i]), "lfc": float(lfc[i]), "padj": float(padj[i])}
                   for i in order_idx[:15]]
            ds["contrasts"][key] = {
                "title": title, "a": a, "b": b,
                "n_sig": int(sig.sum()), "n_sig_lfc05": int(sig05.sum()),
                "n_up": int((sig & (lfc > 0)).sum()), "n_down": int((sig & (lfc < 0)).sum()),
                "top": top,
                "genes": {nm: {"lfc": float(lfc[name2i[nm]]),
                               "padj": (None if np.isnan(padj[name2i[nm]]) else float(padj[name2i[nm]]))}
                          for nm in TARGET if nm in name2i}}

        def gene_vals(sym):
            if sym not in name2i: return None
            i = name2i[sym]
            return {c: [float(x) for x in norm[i, cols[c]]] for c in COND_ORDER}
        ds["targets"] = {s: gene_vals(s) for s in TARGET}
        res["datasets"][tp] = ds

    (OUT / "rnaseq_results.json").write_text(json.dumps(res, indent=2))
    print_summary(res)
    return res


def print_summary(res):
    print("=" * 64)
    for tp in ["day10", "day21"]:
        ds = res["datasets"][tp]
        print(f"\n{DAY_LABEL[tp]}: {ds['n_expressed']}/{ds['n_genes']} expressed")
        for k, c in ds["contrasts"].items():
            cx = c["genes"].get("CXCL11", {})
            pa = cx.get("padj"); pa = float("nan") if pa is None else pa
            print(f"  {k:20s} DEG={c['n_sig']:5d} (up {c['n_up']}, dn {c['n_down']})"
                  f"  |log2FC|>0.5: {c['n_sig_lfc05']:5d}"
                  f"  CXCL11 lfc={cx.get('lfc', float('nan')):.2f} padj={pa:.1e}")
    print("=" * 64)


if __name__ == "__main__":
    analyze()
