#!/usr/bin/env python3
"""Fig S5A: rested-state volcanoes, CRISPRa CXCL11-targeted vs non-targeted.

Real IFNg training (Rested vs Naive) as the reference panel, then CRISPRa priming
at day 10 and day 21. Shared log2FC x; -log10 FDR auto-scaled per panel, since the
pipelines differ in depth and FDR is not comparable across them. Only CXCL11 is
labelled. Reads the pre-computed DESeq2 tables.
"""
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import config as cfg
import figstyle as fs

TAB = cfg.out_dir("crispra_rnaseq", "tables")
G = "γ"  # gamma

# label, path, name_col
COND = [
    (f"IFN{G}\nat rest",          cfg.path("rnaseq", "deseq2_naive_vs_rested"), None),
    ("CRISPRa d10\nat rest",      TAB / "day10_priming_at_rest.csv",  "gene_name"),
    ("CRISPRa d21\nat rest",      TAB / "day21_priming_at_rest.csv",  "gene_name"),
]

PADJ, CUT = 0.05, 0.5                       # DEG = FDR<0.05 & |log2FC|>0.5
CX_BLUE = "#1a56a4"
C_NS, C_SIG = "#cccccc", "#c0392b"           # volcano non-significant / significant


def load(path, name_col):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip() for c in df.columns]
    name = df[name_col] if name_col else pd.Series(df.index, index=df.index)
    out = pd.DataFrame({"gene_name": np.asarray(name, dtype=object),
                        "lfc":  df["log2FoldChange"].to_numpy(float),
                        "padj": df["padj"].to_numpy(float)})
    return out.dropna(subset=["lfc"]).reset_index(drop=True)


DATA = [(lab, load(p, nc)) for (lab, p, nc) in COND]


def _volcano(ax, df, xlim):
    """One volcano panel: x = log2FC (shared), y = -log10 FDR (auto), only CXCL11 labelled."""
    x = df["lfc"].to_numpy(); padj = df["padj"].to_numpy(); names = df["gene_name"].to_numpy()
    m = np.isfinite(padj)                         # a volcano needs padj on y
    x = x[m]; y = -np.log10(np.clip(padj[m], 1e-300, 1.0)); nm = names[m]
    sig = (padj[m] < PADJ) & (np.abs(x) > CUT)
    ax.scatter(x[~sig], y[~sig], s=3.5, c=C_NS, linewidths=0)
    ax.scatter(x[sig], y[sig], s=5, c=C_SIG, linewidths=0)
    ax.axhline(-np.log10(PADJ), color="#888", lw=0.6, ls="--")
    ax.axvline(CUT, color="#888", lw=0.6, ls="--"); ax.axvline(-CUT, color="#888", lw=0.6, ls="--")
    w = np.where(nm == "CXCL11")[0]
    if len(w):
        j = w[0]; left = x[j] > xlim[1] - 2.0     # near the right edge -> label to the left
        ax.scatter([x[j]], [y[j]], s=18, c=CX_BLUE, zorder=5)
        ax.annotate("CXCL11", (x[j], y[j]), xytext=(-4 if left else 4, 2),
                    textcoords="offset points", ha="right" if left else "left",
                    fontsize=6, weight="bold", color=CX_BLUE, zorder=6)
    ax.set_xlim(*xlim)
    inw = (x >= xlim[0]) & (x <= xlim[1])
    ytop = float(y[inw].max()) if inw.any() else float(y.max())
    ax.set_ylim(0, ytop * 1.08 + 1)
    n = int(sig.sum()); up = int((sig & (x > 0)).sum()); dn = n - up
    ax.text(0.03, 0.97, f"{n} DEG\n({up} up, {dn} down)", transform=ax.transAxes,
            ha="left", va="top", fontsize=6, color="#555")
    fs.prism(ax)


def fig_rested_volcano():
    dfs = [d[1] for d in DATA]
    def cxval(df):
        w = df.index[df["gene_name"] == "CXCL11"]
        return float(df.loc[w[0], "lfc"]) if len(w) else np.nan
    los, his = [], []
    for df in dfs:
        x = df["lfc"].to_numpy(); p = df["padj"].to_numpy()
        sig = (p < PADJ) & (np.abs(x) > CUT)
        v = x[sig] if sig.sum() >= 10 else x[np.isfinite(x)]
        los.append(np.percentile(v, 1.5)); his.append(np.percentile(v, 98.5))
    lo, hi = min(los), max(his)
    for df in dfs:                                          # force every CXCL11 call-out in view
        c = cxval(df)
        if np.isfinite(c):
            lo = min(lo, c); hi = max(hi, c)
    xS = (float(np.floor(lo - 0.6)), float(np.ceil(hi + 0.6)))

    fig, ax = plt.subplots(1, 3, figsize=(7.6, 2.75))
    fig.subplots_adjust(top=0.72, bottom=0.24, left=0.085, right=0.985, wspace=0.3)
    for a, df in zip(ax, dfs):
        _volcano(a, df, xS)
        a.set_xlabel("log2 fold-change")
    ax[0].set_ylabel("−log10 FDR")
    for a, day in zip(ax, ["Day 10", "Day 10", "Day 21"]):
        a.set_title(day, weight="bold", pad=3)
    # spanning system + contrast headers above the per-panel "Day X" titles
    p0, p1, p2 = (ax[i].get_position() for i in range(3))
    fig.text((p0.x0 + p0.x1) / 2, 0.885, f"IFN{G} training: Rested vs Naive",
             ha="center", va="bottom", weight="bold")
    fig.text((p1.x0 + p2.x1) / 2, 0.885, "CRISPRa priming: CXCL11-Targeted vs Non-Targeted",
             ha="center", va="bottom", weight="bold")
    fs.save(fig, cfg.out_dir("FigS5A") / "FigS5A")
    return xS


if __name__ == "__main__":
    fs.setup()
    print(f"DEG cutoff = FDR<{PADJ} & |log2FC|>{CUT}")
    xS = fig_rested_volcano()
    print(f"  rested-volcano shared x = {xS}")
    for lab, df in DATA:
        w = df.index[df["gene_name"] == "CXCL11"]
        if len(w):
            r = df.loc[w[0]]
            print(f"    {lab.replace(chr(10), ' '):22s} CXCL11 lfc={r['lfc']:+.2f} padj={r['padj']:.1e}")
