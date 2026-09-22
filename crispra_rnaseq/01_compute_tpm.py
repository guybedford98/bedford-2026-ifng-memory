#!/usr/bin/env python3
"""Compute CXCL11 (and IDO1) TPM for the CRISPRa RNA-seq.

Plasmidsaurus exports counts + CPM only (no length, no TPM). TPM needs gene
lengths; TPM is invariant to a global length rescale, so we use the paper's own
per-gene lengths (RSEM, tabulated in the IFNg-training combined-TPM workbook)
applied to the CRISPRa featureCounts counts:  TPM_i = (c_i/L_i) / sum_j(c_j/L_j) * 1e6.
This keeps the CRISPRa data on the same length basis as the IFNg training TPM panels.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))

import csv, json, re
import numpy as np
import openpyxl

import config as cfg

TPMXLSX = cfg.path("crispra_rnaseq", "gene_lengths")
MATS = {"day10": (cfg.path("crispra_rnaseq", "counts_day10"), "CF63KH"),
        "day21": (cfg.path("crispra_rnaseq", "counts_day21"), "GVSYGD")}
OUT = cfg.out_dir("crispra_rnaseq")

COND = {**{n: "NT_rest" for n in (1,2,3)}, **{n: "NT_IFNg" for n in (4,5,6)},
        **{n: "sg1_rest" for n in (7,8,9)}, **{n: "sg1_IFNg" for n in (10,11,12)}}
COND_ORDER = ["NT_rest","NT_IFNg","sg1_rest","sg1_IFNg"]
GENES = ["CXCL11","IDO1"]


def load_lengths():
    wb = openpyxl.load_workbook(TPMXLSX, read_only=True, data_only=True)
    ws = wb.active
    it = ws.iter_rows(values_only=True)
    next(it)                      # row1 merged category labels
    hdr = list(next(it))          # row2 = real header: gene_id, length, ...
    gi = hdr.index("gene_id"); li = hdr.index("length")
    L = {}
    for row in it:
        sym = row[gi]
        if sym is None: continue
        try: length = float(row[li])
        except (TypeError, ValueError): continue
        if length and length > 0 and sym not in L:
            L[str(sym)] = length
    wb.close()
    return L


def load_counts(path, proj):
    with open(path, newline="", encoding="utf-8") as fh:
        rdr = csv.reader(fh, delimiter="\t")
        hdr = next(rdr)
        name_i = hdr.index("gene_name")
        cnt_cols = {}
        for j, h in enumerate(hdr):
            m = re.match(re.escape(proj) + r"_(\d+)_count$", h)
            if m: cnt_cols[int(m.group(1))] = j
        order = sorted(cnt_cols)
        names = []; counts = []
        for row in rdr:
            names.append(row[name_i])
            counts.append([float(row[cnt_cols[n]]) for n in order])
    return names, np.array(counts), order


def main():
    L = load_lengths()
    print(f"length map: {len(L)} genes; CXCL11={L.get('CXCL11')}, IDO1={L.get('IDO1')}")
    out = {}
    for tp, (path, proj) in MATS.items():
        names, counts, order = load_counts(path, proj)          # counts: genes x samples
        lengths = np.array([L.get(n, np.nan) for n in names])
        have = ~np.isnan(lengths)
        rate = np.zeros_like(counts)
        rate[have] = counts[have] / lengths[have, None]
        denom = rate.sum(axis=0)                                 # per-sample sum of count/length
        tpm = rate / denom * 1e6                                 # genes x samples
        # matched fraction of total counts (quality check)
        frac = counts[have].sum() / counts.sum()
        name2i = {}
        for i, n in enumerate(names): name2i.setdefault(n, i)
        out[tp] = {"matched_count_frac": float(frac), "targets": {}}
        for g in GENES:
            i = name2i[g]
            per_cond = {c: [] for c in COND_ORDER}
            for k, n in enumerate(order):
                per_cond[COND[n]].append(float(tpm[i, k]))
            out[tp]["targets"][g] = per_cond
        print(f"[{tp}] matched {frac*100:.1f}% of counts; genes w/ length {have.sum()}/{len(names)}")
        for g in GENES:
            pc = out[tp]["targets"][g]
            means = {c: np.mean(pc[c]) for c in COND_ORDER}
            print(f"    {g} TPM means: " + "  ".join(f"{c}={means[c]:.1f}" for c in COND_ORDER))
    outp = OUT / "cxcl11_tpm.json"
    outp.write_text(json.dumps(out, indent=2))
    print("wrote", outp)


if __name__ == "__main__":
    main()
