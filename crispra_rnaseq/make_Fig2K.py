#!/usr/bin/env python3
"""Fig 2K: MA plots, three panels left -> right.

Real IFNg training (the broad-program reference), then CRISPRa synthetic memory
at day 10 and at day 21, showing the decay. IFNg training is a single ~10-day-rested
re-stim, so it has no day-21 counterpart and stays one reference panel.

DEG cutoff is FDR < 0.05 & |log2FC| > 0.5. log2FC axis shared across panels, log10
mean-count x shared out to 5.5. Reads the pre-computed DESeq2 tables, so pydeseq2
is not needed here.
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

# re-challenge tables: CRISPRa (gene_name column) and IFNg training (index = gene_name)
CRISPRA_TAB = (TAB / "day10_memory_rechallenge.csv", "gene_name")
TRAINING_TAB = (cfg.path("rnaseq", "deseq2_stim_vs_rechallenge"), None)

PADJ, LFC_T = 0.05, 0.5
C_NS, C_SIG = "#cccccc", "#c0392b"
CX_BLUE = "#1a56a4"


def load(path, name_col):
    df = pd.read_csv(path, index_col=0)
    df.columns = [c.strip() for c in df.columns]
    name = df[name_col] if name_col else pd.Series(df.index, index=df.index)
    out = pd.DataFrame({"gene_name": np.asarray(name, dtype=object),
                        "base": df["baseMean"].to_numpy(float),
                        "lfc":  df["log2FoldChange"].to_numpy(float),
                        "padj": df["padj"].to_numpy(float)})
    # keep NaN-padj genes (they have a valid baseMean + log2FC -> grey in the MA cloud)
    return out.dropna(subset=["base", "lfc"]).reset_index(drop=True)


TRAINING  = load(*TRAINING_TAB)
CRISPRA10 = load(*CRISPRA_TAB)                                       # day 10 re-challenge
CRISPRA21 = load(TAB / "day21_memory_rechallenge.csv", "gene_name")  # day 21 re-challenge


def callout_lfcs(dfs):
    vals = []
    for d in dfs:
        for g in ("CXCL11", "IDO1"):
            w = d.index[d["gene_name"] == g]
            if len(w):
                vals.append(float(d.loc[w[0], "lfc"]))
    return vals


def shared_lfc_range(dfs, pct=1.5, pad=0.7):
    """Robust shared log2FC range over the panels, forced to include the call-outs."""
    los = [np.percentile(d["lfc"].to_numpy(), pct) for d in dfs]
    his = [np.percentile(d["lfc"].to_numpy(), 100 - pct) for d in dfs]
    lo, hi = min(los), max(his)
    for v in callout_lfcs(dfs):
        lo = min(lo, v); hi = max(hi, v)
    return (float(np.floor(lo - pad)), float(np.ceil(hi + pad)))


def ma_panel(ax, df, title, sublabel, ylim, xlim):
    base = df["base"].to_numpy(); lfc = df["lfc"].to_numpy(); padj = df["padj"].to_numpy()
    names = df["gene_name"].to_numpy()
    ok = base > 0
    x = np.log10(base[ok]); y = lfc[ok]; p = padj[ok]; nm = names[ok]
    sig = (p < PADJ) & (np.abs(y) > LFC_T)
    ax.scatter(x[~sig], y[~sig], s=3.5, c=C_NS, linewidths=0)
    ax.scatter(x[sig], y[sig], s=5, c=C_SIG, linewidths=0)
    ax.axhline(0, color="#888", lw=0.7)
    ax.axhline(LFC_T, color="#888", lw=0.6, ls="--")
    ax.axhline(-LFC_T, color="#888", lw=0.6, ls="--")
    for g, col, ha, off in (("CXCL11", CX_BLUE, "left", (3, 3)),):
        w = np.where(nm == g)[0]
        if len(w):
            j = w[0]
            ax.scatter([x[j]], [y[j]], s=16, c=col, zorder=5)
            ax.annotate(g, (x[j], y[j]), fontsize=6, weight="bold", ha=ha,
                        xytext=off, textcoords="offset points", color=col, zorder=6)
    ax.set_xlim(*xlim); ax.set_ylim(*ylim)
    n_off = int((((y < ylim[0]) | (y > ylim[1])) & sig).sum())   # significant only
    if n_off:
        ax.text(0.97, 0.03, f"+{n_off} DEG off-axis", transform=ax.transAxes,
                ha="right", va="bottom", fontsize=6, color="#999", style="italic")
    n_sig = int(sig.sum()); n_up = int((sig & (y > 0)).sum()); n_dn = n_sig - n_up
    ax.text(0.03, 0.95, f"{n_sig} DEG\n({n_up} up, {n_dn} down)", transform=ax.transAxes,
            ha="left", va="top", fontsize=6, color="#444")
    ax.set_xlabel("log10 mean norm. counts" + (f"\n{sublabel}" if sublabel else ""))
    ax.set_ylabel("log2 fold-change")
    ax.set_title(title)
    fs.prism(ax)
    return n_sig, n_up, n_dn


# panels left -> right: real IFNg-training reference, then CRISPRa memory at day 10 & 21.
C_SUB = f"CXCL11-Targeted+IFN{G}\nvs Non-Targeted+IFN{G}"
PANELS = [
    dict(df=TRAINING,  title=f"IFN{G} training", sub="Rechallenged\nvs Stimulated"),
    dict(df=CRISPRA10, title="CRISPRa day 10",    sub=C_SUB),
    dict(df=CRISPRA21, title="CRISPRa day 21",    sub=C_SUB),
]
DFS = [p["df"] for p in PANELS]


def build_MA():
    yR = shared_lfc_range(DFS)
    def logb(d):
        b = d["base"].to_numpy(); return np.log10(b[b > 0])
    xmin = min(logb(d).min() for d in DFS)
    xS = (float(np.floor(xmin * 10) / 10 - 0.15), 5.5)   # shared x, all extend to 5.5
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.2))
    stats = [ma_panel(a, p["df"], p["title"], p["sub"], ylim=yR, xlim=xS)
             for a, p in zip(ax, PANELS)]
    fig.tight_layout(w_pad=1.4)
    fs.save(fig, cfg.out_dir("Fig2K") / "Fig2K")
    return stats


if __name__ == "__main__":
    fs.setup()
    ma_stats = build_MA()
    for p, s in zip(PANELS, ma_stats):
        print(f"  MA      {p['title'].replace(G, 'g'):16s} DEG={s}")
