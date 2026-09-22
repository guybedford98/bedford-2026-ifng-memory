"""Fig S2E: top enriched GO biological processes for the naive-versus-rested contrast,
as a diverging NES bar chart. Terms elevated in rested cells to the right, terms reduced
in rested cells to the left.

Sign convention: positive NES = higher in RESTED than naive. Verified against IDO1,
which carries log2FoldChange +2.9 in the Naive_vs_Rested DESeq2 table it was ranked on.

Redundant GO terms are collapsed. clusterProfiler reports alias terms that wrap the
identical gene set (e.g. "extracellular matrix organization" / "extracellular structure
organization" / "external encapsulating structure organization", all 261 genes at
NES 2.132). Terms sharing a setSize and an NES to 6 dp are the same set; the shortest
description is kept and the number collapsed is recorded in the console output.

Reads the GSEA table written by 00_run_gsea.R.
"""

import sys
import textwrap
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")

GSEA_TABLE = cfg.out_dir("gsea_tables") / "GSEA_GO_BP_Naive_vs_Rested_results.csv"

N_PER_SIDE = 5
FDR_CUT = 0.05
UP_COLOR = "#b2182b"      # higher in rested
DOWN_COLOR = "#2166ac"    # lower in rested
EDGE = fs.BAR_EDGE

# Terms called out by name in the manuscript text; drawn in bold.
HIGHLIGHT = {"extracellular matrix organization", "cytoplasmic translation"}

# Display-only shortenings for GO names that would otherwise wrap to three or four
# lines and force the panel taller. Standard abbreviations only, no meaning changed.
# The full GO description is preserved in the *_plotted_terms.csv written alongside.
SHORT_NAME = {
    "peroxisome proliferator activated receptor signaling pathway":
        "PPAR signaling pathway",
    "regulation of blood vessel endothelial cell proliferation involved in sprouting angiogenesis":
        "regulation of sprouting angiogenesis",
    "antigen processing and presentation of exogenous peptide antigen via MHC class II":
        "antigen presentation, exogenous peptide (MHC II)",
    "antigen processing and presentation of peptide antigen via MHC class II":
        "antigen presentation, peptide (MHC II)",
    "antigen processing and presentation of exogenous peptide antigen":
        "antigen presentation, exogenous peptide",
    "antigen processing and presentation of exogenous antigen":
        "antigen presentation, exogenous antigen",
    "antigen processing and presentation of peptide antigen":
        "antigen presentation, peptide antigen",
}

PANELS = [
    ("IFNg", GSEA_TABLE, "IFN$\\gamma$: naive vs rested", "FigS2E"),
]


def collapse_aliases(df):
    """Drop GO terms that wrap an identical gene set (same setSize, same NES)."""
    df = df.copy()
    df["_key"] = list(zip(df["setSize"], df["NES"].round(6)))
    df["_len"] = df["Description"].str.len()
    df = df.sort_values(["_key", "_len"])
    kept = df.drop_duplicates("_key", keep="first")
    n_dropped = len(df) - len(kept)
    return kept.drop(columns=["_key", "_len"]), n_dropped


def make_panel(tag, path, title, stem):
    d = pd.read_csv(path)
    sig = d[d["p.adjust"] < FDR_CUT]
    sig, n_alias = collapse_aliases(sig)
    # Top N each way, sized to drop into a supplementary panel slot. A term the
    # manuscript names is never dropped: if it falls outside the top N it displaces the
    # weakest entry on its own side, and the swap is reported.
    def pick(pool, n, largest):
        top = pool.nlargest(n, "NES") if largest else pool.nsmallest(n, "NES")
        forced = pool[pool["Description"].isin(HIGHLIGHT)]
        forced = forced[forced["NES"] > 0] if largest else forced[forced["NES"] < 0]
        missing = forced[~forced["Description"].isin(top["Description"])]
        if len(missing):
            print("  forced in (named in the text, outside the top %d): %s"
                  % (n, ", ".join(missing["Description"])))
            top = pd.concat([top.head(n - len(missing)), missing])
        return top

    up = pick(sig, N_PER_SIDE, True)
    dn = pick(sig, N_PER_SIDE, False)
    sel = pd.concat([dn, up]).sort_values("NES").reset_index(drop=True)

    print("\n=== %s" % tag)
    print("  terms at FDR < %.2f: %d (%d alias duplicates collapsed)"
          % (FDR_CUT, len(sig), n_alias))
    for _, r in sel.iloc[::-1].iterrows():
        print("   %+.2f  FDR %-9.2g  %s" % (r["NES"], r["p.adjust"], r["Description"]))
    missing = HIGHLIGHT - set(sel["Description"])
    if missing:
        print("  NOTE: highlighted term(s) not in the plotted top %d: %s"
              % (N_PER_SIDE, ", ".join(sorted(missing))))

    # Row height set so a 3-line wrapped GO name at 8 pt clears its neighbours.
    fig, ax = plt.subplots(figsize=(3.2, 0.36 * len(sel) + 0.90))
    colors = [UP_COLOR if v > 0 else DOWN_COLOR for v in sel["NES"]]
    y = range(len(sel))
    ax.barh(y, sel["NES"], color=colors, edgecolor=EDGE,
            linewidth=fs.BAR_EDGE_LW, height=0.72)

    for i, r in sel.iterrows():
        off = 0.06 if r["NES"] > 0 else -0.06
        ax.text(r["NES"] + off, i, "%.0e" % r["p.adjust"],
                va="center", ha="left" if r["NES"] > 0 else "right",
                fontsize=6, color="#555555")

    # Long GO names blow out the panel width; wrap rather than truncate.
    display = [SHORT_NAME.get(t, t) for t in sel["Description"]]
    wrapped = ["\n".join(textwrap.wrap(t, 30)) for t in display]
    ax.set_yticks(list(y))
    ax.set_yticklabels(wrapped, fontsize=8)
    # Match on the UNWRAPPED description: a wrapped label no longer equals its GO name.
    for lab, term in zip(ax.get_yticklabels(), sel["Description"]):
        if term in HIGHLIGHT:
            lab.set_fontweight("bold")
            lab.set_color("#000000")

    ax.axvline(0, color=EDGE, linewidth=0.8)
    ax.set_xlabel("Normalized enrichment score", fontsize=8)
    ax.set_title(title, fontsize=8)
    lim = max(abs(sel["NES"])) * 1.30
    ax.set_xlim(-lim, lim)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.tick_params(axis="y", length=0)
    ax.tick_params(axis="x", length=2.4, width=0.8, direction="out")

    # The diverging layout leaves the top-left and bottom-right corners empty.
    ax.text(0.02, 0.985, "← lower in rested", transform=ax.transAxes,
            ha="left", va="top", fontsize=8, color=DOWN_COLOR)
    ax.text(0.98, 0.015, "higher in rested →", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=8, color=UP_COLOR)

    out = cfg.out_dir(stem)
    fs.save(fig, out / stem)
    print("  wrote %s.{pdf,svg,png}" % stem)
    sel.to_csv(out / (stem + "_plotted_terms.csv"), index=False)


fs.setup()
for args in PANELS:
    make_panel(*args)
