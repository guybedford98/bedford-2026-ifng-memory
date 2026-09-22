"""Fig 1F: k-means heatmap of the 1,063-gene protein-coding ISG universe.

Per-gene z-score across the four training conditions -> KMeans -> seaborn clustermap,
row_cluster=False, col_cluster=False, cmap="vlag", categorical row-color strip.

The panel is heatmap + cluster color strip only. Cluster names, brackets, the
"#Genes per Cluster" column and the z-score legend are added on the artboard.

Usage:  py make_Fig1F.py [--k N]
"""

import sys
import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

TABLES = cfg.out_dir("rnaseq_tables")
INPUT_CSV = TABLES / "ISG_universe_maxTPM5_proteincoding.csv"

# k for the manuscript panel: k=5 recovers all 50 positive memory genes and has the
# better silhouette of the candidates considered.
PANEL_K = 5
N_CLUSTERS = int(sys.argv[sys.argv.index("--k") + 1]) if "--k" in sys.argv else PANEL_K
RANDOM_STATE = 42
# The panel k keeps the plain filenames. Any other k is a variant and gets a _kN
# suffix, so trying alternatives never overwrites the panel built for the figure.
SUFFIX = "" if N_CLUSTERS == PANEL_K else f"_k{N_CLUSTERS}"

# Panel geometry, measured off the heatmap as it sits on the artboard: the four
# condition columns occupy 74.04 x 213.96 pt and the categorical colour strip is
# 3.84 pt wide and abuts them with no gap. Building at exactly that size means the
# export drops onto the artboard at 100% instead of being scaled by hand.
HEAT_W_IN = 74.04 / 72
HEAT_H_IN = 213.96 / 72
STRIP_W_IN = 3.84 / 72
# margins only need to be big enough for the labels; bbox_inches="tight" trims the rest
MARGIN = {"left": 0.06, "right": 0.55, "bottom": 0.80, "top": 0.06}
K_SWEEP = range(2, 11)

# TPM column in the universe table -> condition label on the panel, in panel order
CONDITIONS = {
    "NTC": "Naive",
    "IFNg": "Stimulated",
    "Trained": "Rested",
    "Rechallenge": "Rechallenged",
}
COLS = list(CONDITIONS.values())

# Genes labelled on the right of the panel. Every other row is unlabelled.
SELECTED_GENES = [
    "IRF1", "TLR3", "SOCS3", "HK2", "FAS", "CIITA", "STAT1", "IL32", "CXCL11",
    "GBP5", "IDO1", "IL6", "PLOD2", "STING1", "BATF3", "CARD6", "MET", "KRAS",
    "BCL10", "GBP1", "CCL2",
]

# The four grouped response classes, top to bottom as in the published panel
CLASS_ORDER = [
    "Non-Memory Like Response",
    "Memory Like Response",
    "Refractory Like Response",
    "Other Diverse Response",
]

# Published k=6 cluster centroids from the earlier 1,179-gene TPM-sum>=10 heatmap.
# Used only as an independent archetype check on the new clusters, never to assign them.
PUBLISHED_ARCHETYPES = {
    #  name              Naive   Stim    Rested  Rech     n
    "nonmemory":       (-0.992,  1.041, -0.991,  0.943,  514),
    "memory":          (-1.101,  0.581, -0.786,  1.306,  194),
    "additive_memory": (-1.468,  0.059,  0.224,  1.186,   58),
    "refractory":      (-0.756,  1.400, -1.047,  0.404,  208),
    "nonresetting":    (-1.374,  1.089, -0.456,  0.741,  167),
    "undefined":       (-1.285,  1.075,  0.643, -0.434,   38),
}

# Class-assignment rule thresholds. Calibrated so the rule reproduces the published
# cluster naming on BOTH published runs (the 1,179-gene panel and the 1,336-gene
# variant), 12 of 12 clusters, then applied unchanged here.
DELTA_CUT = 0.5     # |Rechallenged - Stimulated| in z units
RESET_CUT = 0.30    # (Rested - Naive) / (Stimulated - Naive)

# Diverging scale. With four conditions and a population sigma, the per-gene z-score is
# bounded at +/- sqrt(n-1) = sqrt(3) = 1.732 by construction (the extreme case is one
# condition away from three equal ones). That is the true full range of the matrix and it
# is symmetric about zero, so fixing the scale here clips nothing and makes the z-score
# legend well defined.
VMAX = float(np.sqrt(len(CONDITIONS) - 1))
VMIN = -VMAX


def classify(delta, reset, reset_cut=RESET_CUT):
    """The calibrated cluster -> response class rule."""
    if delta >= DELTA_CUT:
        return CLASS_ORDER[1], f"delta {delta:+.3f} >= {DELTA_CUT:+.2f}"
    if np.isfinite(reset) and reset >= reset_cut:
        return CLASS_ORDER[3], f"reset {reset:+.3f} >= {reset_cut:+.2f}"
    if delta <= -DELTA_CUT:
        return CLASS_ORDER[2], f"delta {delta:+.3f} <= {-DELTA_CUT:+.2f}"
    return CLASS_ORDER[0], f"delta {delta:+.3f} within +/-{DELTA_CUT:.2f}, reset {reset:+.3f}"


def rule_line(n, s, r, q, reset_cut=RESET_CUT):
    """Return (class, delta, reset, reason) for one cluster centroid."""
    delta = q - s
    span = s - n
    reset = (r - n) / span if abs(span) > 1e-9 else np.nan
    cls, why = classify(delta, reset, reset_cut)
    return cls, delta, reset, why


def spread_labels(pos, min_sep, lo, hi):
    """
    Push a sorted list of label positions apart to at least min_sep, staying inside
    [lo, hi]. Standard two-pass relaxation; returns positions in the same order.
    """
    pos = list(pos)
    for i in range(1, len(pos)):
        pos[i] = max(pos[i], pos[i - 1] + min_sep)
    overflow = pos[-1] - hi
    if overflow > 0:
        pos = [p - overflow for p in pos]
    for i in range(len(pos) - 2, -1, -1):
        pos[i] = min(pos[i], pos[i + 1] - min_sep)
    underflow = lo - pos[0]
    if underflow > 0:
        pos = [p + underflow for p in pos]
    return pos


def main():
    print("=" * 78)
    print("Fig 1F: k-means heatmap, 1,063 protein-coding ISGs")
    print("=" * 78)
    print(f"input: {INPUT_CSV}")

    raw = pd.read_csv(INPUT_CSV)
    n_in = len(raw)
    print(f"genes in:  {n_in}")
    if raw["gene_id"].duplicated().any():
        raise SystemExit("duplicate gene_id in the input universe")

    tpm = raw.set_index("gene_id")[list(CONDITIONS)].rename(columns=CONDITIONS)
    if tpm.isna().to_numpy().any():
        raise SystemExit("missing TPM values in the input universe")

    universe_class = raw.set_index("gene_id")["class"]

    # Per-gene z-score across the four conditions, sigma population (ddof=0).
    mu = tpm.mean(axis=1)
    sigma = tpm.std(axis=1, ddof=0)
    if (sigma == 0).any():
        raise SystemExit(f"{(sigma == 0).sum()} genes have zero variance and cannot be z-scored")
    z = tpm.sub(mu, axis=0).div(sigma, axis=0)
    zmin, zmax = float(z.to_numpy().min()), float(z.to_numpy().max())
    print(f"z-scored:  {len(z)} genes x {z.shape[1]} conditions "
          f"(range {zmin:+.3f} to {zmax:+.3f}; colour scale fixed at "
          f"{VMIN:+.3f} to {VMAX:+.3f}, the theoretical bound, so nothing clips)")
    if zmin < VMIN or zmax > VMAX:
        raise SystemExit("z-scores fall outside the fixed colour scale; values would clip")
    z.to_csv(TABLES / "ISG_1063_zscore.csv")

    X = z.to_numpy()

    print("\n" + "-" * 78)
    print("silhouette sweep (KMeans, random_state=42)")
    print("-" * 78)
    sweep = []
    for k in K_SWEEP:
        km = KMeans(n_clusters=k, random_state=RANDOM_STATE).fit(X)
        score = silhouette_score(X, km.labels_)
        sweep.append({"k": k, "silhouette": round(float(score), 4),
                      "inertia": round(float(km.inertia_), 3)})
        print(f"  k={k:2d}   silhouette = {score:.4f}   inertia = {km.inertia_:9.2f}")
    sweep_df = pd.DataFrame(sweep)
    sweep_df.to_csv(TABLES / "silhouette_sweep_k2_k10.csv", index=False)
    best_k = int(sweep_df.loc[sweep_df["silhouette"].idxmax(), "k"])
    print(f"\n  best silhouette at k={best_k} "
          f"({sweep_df['silhouette'].max():.4f}); k={N_CLUSTERS} {'the panel k' if N_CLUSTERS == PANEL_K else 'a variant, not the panel k'} "
          f"(silhouette "
          f"{float(sweep_df.loc[sweep_df['k'] == N_CLUSTERS, 'silhouette'].iloc[0]):.4f}).")

    km = KMeans(n_clusters=N_CLUSTERS, random_state=RANDOM_STATE).fit(X)
    labels = pd.Series(km.labels_, index=z.index, name="Cluster")

    cent = pd.DataFrame(km.cluster_centers_, columns=COLS)
    cent.index.name = "Cluster"
    cent["n"] = labels.value_counts().sort_index()

    print("\n" + "-" * 78)
    print(f"KMeans k={N_CLUSTERS}, random_state={RANDOM_STATE}: centroids (z units)")
    print("-" * 78)
    print(cent.round(3).to_string())

    # nearest published archetype, Euclidean on the four centroid values
    arche = []
    for c in cent.index:
        v = cent.loc[c, COLS].to_numpy(dtype=float)
        d = {nm: float(np.linalg.norm(v - np.array(a[:4]))) for nm, a in PUBLISHED_ARCHETYPES.items()}
        nearest = min(d, key=d.get)
        arche.append({"Cluster": c, "nearest_published_archetype": nearest,
                      "distance": round(d[nearest], 3)})
    arche_df = pd.DataFrame(arche).set_index("Cluster")
    print("\nnearest published cluster archetype:")
    print(arche_df.to_string())

    xt = pd.crosstab(labels, universe_class.reindex(labels.index))
    xt = xt.reindex(columns=[c for c in ["positive memory", "refractory", "non-memory"]
                             if c in xt.columns])
    print("\n" + "-" * 78)
    print("cluster x universe `class` (TPM-ratio call, independent of the clustering)")
    print("-" * 78)
    print(xt.to_string())
    pct = (xt.T / xt.sum(axis=1)).T * 100
    print("\nrow %:")
    print(pct.round(1).to_string())
    xt.to_csv(TABLES / f"crosstab_cluster_vs_universe_class{SUFFIX}.csv")

    print("\n" + "-" * 78)
    print("cluster -> response class, from the centroids")
    print("-" * 78)
    assign = {}
    rows = []
    for c in cent.index:
        n, s, r, q = (float(cent.loc[c, col]) for col in COLS)
        cls, delta, reset, why = rule_line(n, s, r, q)
        assign[c] = cls
        rows.append({"Cluster": c, "n": int(cent.loc[c, "n"]),
                     "Naive": round(n, 3), "Stimulated": round(s, 3),
                     "Rested": round(r, 3), "Rechallenged": round(q, 3),
                     "delta_Rech_minus_Stim": round(delta, 3),
                     "reset_fraction": round(reset, 3),
                     "nearest_published_archetype": arche_df.loc[c, "nearest_published_archetype"],
                     "response_class": cls, "reason": why})
        print(f"  cluster {c} (n={int(cent.loc[c, 'n']):4d})  ->  {cls:26s}  [{why}]")
    assign_df = pd.DataFrame(rows).set_index("Cluster")
    assign_df.to_csv(TABLES / f"cluster_centroids_and_class_assignment{SUFFIX}.csv")

    # sensitivity of the grouping to the reset threshold, so any borderline cluster is
    # on record rather than being decided silently by the cut-off
    print("\nsensitivity of the grouping to the reset threshold:")
    sens = []
    for cut in (0.20, 0.25, 0.30, 0.35, 0.40):
        alt = {}
        for c in cent.index:
            n, s, r, q = (float(cent.loc[c, col]) for col in COLS)
            alt[c] = rule_line(n, s, r, q, reset_cut=cut)[0]
        tot = {cls: int(sum(int(cent.loc[c, "n"]) for c in cent.index if alt[c] == cls))
               for cls in CLASS_ORDER}
        moved = [int(c) for c in cent.index if alt[c] != assign[c]]
        sens.append({"reset_cut": cut, **{k: v for k, v in tot.items()},
                     "clusters_moved_vs_default": ",".join(str(c) for c in moved)})
        print(f"  reset_cut {cut:.2f}  " +
              "  ".join(f"{cls.split()[0]:11s}{tot[cls]:5d}" for cls in CLASS_ORDER) +
              (f"   moves cluster {moved}" if moved else ""))
    pd.DataFrame(sens).to_csv(TABLES / f"sensitivity_reset_threshold{SUFFIX}.csv", index=False)

    # Row order: classes top to bottom in the published order, clusters within a class
    # by descending size, genes within a cluster left in input order (stable sort).
    block_order = []
    for cls in CLASS_ORDER:
        members = [c for c in cent.index if assign[c] == cls]
        members.sort(key=lambda c: -int(cent.loc[c, "n"]))
        block_order.extend(members)
    if sorted(block_order) != sorted(cent.index.tolist()):
        raise SystemExit("block ordering lost or duplicated a cluster")

    rank = {c: i for i, c in enumerate(block_order)}
    order_key = labels.map(rank)
    z_sorted = z.loc[order_key.sort_values(kind="stable").index]
    lab_sorted = labels.loc[z_sorted.index]

    n_plot = len(z_sorted)
    print("\n" + "-" * 78)
    print(f"genes in {n_in}, genes plotted {n_plot}, "
          f"{'no loss' if n_in == n_plot else 'LOSS ' + str(n_in - n_plot)}")
    print("-" * 78)
    if n_in != n_plot:
        raise SystemExit(f"gene loss between input ({n_in}) and panel ({n_plot})")

    print("\ncluster sizes, top to bottom as plotted:")
    for c in block_order:
        print(f"  block {block_order.index(c)}  cluster {c}  n={int(cent.loc[c, 'n']):4d}  "
              f"{assign[c]}")
    print("\ngrouped totals, published order:")
    class_totals = {}
    for cls in CLASS_ORDER:
        tot = int(sum(int(cent.loc[c, "n"]) for c in cent.index if assign[c] == cls))
        class_totals[cls] = tot
        print(f"  {cls:26s} {tot:5d}")
    print(f"  {'TOTAL':26s} {sum(class_totals.values()):5d}")

    # Row-color strip: tab20, keyed to the plotted block position (top band = tab20[0])
    # so the bands read in order down the strip.
    palette = sns.color_palette("tab20", N_CLUSTERS)
    block_color = {c: palette[i] for i, c in enumerate(block_order)}
    row_colors = pd.Series(lab_sorted.map(block_color).to_numpy(),
                           index=z_sorted.index, name=None)

    colour_map = pd.DataFrame([
        {"block_position": i, "Cluster": c, "n": int(cent.loc[c, "n"]),
         "response_class": assign[c],
         "hex": matplotlib.colors.to_hex(block_color[c])}
        for i, c in enumerate(block_order)
    ])
    print("\nrow-colour strip, top to bottom:")
    print(colour_map.to_string(index=False))
    colour_map.to_csv(TABLES / f"row_color_strip_map{SUFFIX}.csv", index=False)

    matplotlib.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "pdf.fonttype": 42,          # keep text editable
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.linewidth": 0.5,
    })

    fig_w = MARGIN["left"] + STRIP_W_IN + HEAT_W_IN + MARGIN["right"]
    fig_h = MARGIN["bottom"] + HEAT_H_IN + MARGIN["top"]

    g = sns.clustermap(
        z_sorted,
        row_cluster=False,
        col_cluster=False,
        cmap="vlag",
        vmin=VMIN, vmax=VMAX,
        figsize=(fig_w, fig_h),
        yticklabels=False,
        xticklabels=True,
        row_colors=row_colors,
        dendrogram_ratio=(0.001, 0.001),
        colors_ratio=0.035,
        cbar_pos=None,               # the z-score legend is added on the artboard
    )

    ax = g.ax_heatmap
    g.ax_row_dendrogram.set_visible(False)
    g.ax_col_dendrogram.set_visible(False)

    # Override clustermap's gridspec with the measured panel geometry. The strip sits
    # flush against the heatmap, as it does on the artboard.
    left_f = MARGIN["left"] / fig_w
    bot_f = MARGIN["bottom"] / fig_h
    h_f = HEAT_H_IN / fig_h
    g.ax_row_colors.set_position([left_f, bot_f, STRIP_W_IN / fig_w, h_f])
    ax.set_position([left_f + STRIP_W_IN / fig_w, bot_f, HEAT_W_IN / fig_w, h_f])

    # condition labels, 8 pt, rotated as in the published panel
    ax.set_xticklabels(COLS, rotation=90, fontsize=8)
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.tick_params(axis="x", length=2, width=0.5, pad=2)

    # gene labels, 6 pt, only the selected list. Rows that fall too close together are
    # pushed apart to a legible spacing and connected back to their own row by a leader
    # line, so a label always points at the gene it names.
    pos = {gene: i for i, gene in enumerate(z_sorted.index)}
    missing = [gene for gene in SELECTED_GENES if gene not in pos]
    present = sorted((gene for gene in SELECTED_GENES if gene in pos), key=lambda gene: pos[gene])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)

    g.fig.canvas.draw()
    ax_h_px = ax.get_window_extent().height
    min_sep = (7.0 * g.fig.dpi / 72.0) / ax_h_px * n_plot   # 6 pt text + 1 pt gap
    row_y = [pos[gene] + 0.5 for gene in present]
    lab_y = spread_labels(row_y, min_sep, 0.5, n_plot - 0.5)
    n_moved = sum(1 for a, b in zip(row_y, lab_y) if abs(a - b) > 0.5)

    # elbow leaders: a short stub off the row, a diagonal, then a short stub into the
    # label. Straight lines fan into a tangle where several labelled genes sit within a
    # few rows of each other; the elbow keeps them readable at this width.
    x_edge = len(COLS)
    for gene, y_row, y_lab in zip(present, row_y, lab_y):
        ax.add_line(matplotlib.lines.Line2D(
            [x_edge, x_edge + 0.10, x_edge + 0.22, x_edge + 0.30],
            [y_row, y_row, y_lab, y_lab],
            linewidth=0.4, color="#2b2b2b", clip_on=False, solid_joinstyle="miter",
        ))
        ax.text(x_edge + 0.36, y_lab, gene, fontsize=6, va="center", ha="left",
                style="italic", clip_on=False)   # italic gene symbols

    g.ax_row_colors.set_xticks([])
    g.ax_row_colors.set_yticks([])
    for spine in g.ax_row_colors.spines.values():
        spine.set_visible(False)

    print("\nlabelled genes: "
          f"{len(present)} of {len(SELECTED_GENES)} requested "
          f"({n_moved} nudged off their row onto a leader line)")
    if missing:
        print(f"  NOT IN THE 1,063-GENE UNIVERSE, not labelled: {', '.join(missing)}")
    else:
        print("  all requested label genes are present in the universe")

    fs.save(g.fig, cfg.out_dir("Fig1F") / f"Fig1F{SUFFIX}")
    plt.close("all")

    out = z_sorted.copy()
    out["Cluster"] = lab_sorted
    out["block_position"] = lab_sorted.map({c: i for i, c in enumerate(block_order)})
    out["response_class"] = lab_sorted.map(assign)
    out["universe_class"] = universe_class.reindex(out.index)
    out["labelled_on_panel"] = out.index.isin(present)
    out.to_csv(TABLES / f"ISG_1063_zscore_with_clusters_ordered{SUFFIX}.csv")

    counts = pd.DataFrame(
        [{"response_class": cls, "n_genes": class_totals[cls],
          "clusters": ",".join(str(c) for c in block_order if assign[c] == cls)}
         for cls in CLASS_ORDER]
        + [{"response_class": "TOTAL", "n_genes": sum(class_totals.values()), "clusters": ""}]
    )
    counts.to_csv(TABLES / f"class_counts_per_response_class{SUFFIX}.csv", index=False)

    summary = {
        "input_csv": str(INPUT_CSV),
        "n_genes_in": n_in,
        "n_genes_plotted": n_plot,
        "gene_loss": n_in - n_plot,
        "n_clusters": N_CLUSTERS,
        "random_state": RANDOM_STATE,
        "zscore": "per gene across the four conditions, sigma ddof=0",
        "heatmap_scale": {"vmin": round(VMIN, 4), "vmax": round(VMAX, 4), "cmap": "vlag",
                          "note": "sqrt(n_conditions - 1), the theoretical z bound",
                          "observed_z_min": round(zmin, 4), "observed_z_max": round(zmax, 4)},
        "silhouette_sweep": sweep,
        "best_silhouette_k": best_k,
        "cluster_sizes": {int(c): int(cent.loc[c, "n"]) for c in cent.index},
        "block_order_top_to_bottom": [int(c) for c in block_order],
        "cluster_to_class": {int(c): assign[c] for c in cent.index},
        "class_totals_published_order": class_totals,
        "assignment_rule": {
            "delta_cut": DELTA_CUT, "reset_cut": RESET_CUT,
            "delta": "Rechallenged - Stimulated, centroid z units",
            "reset_fraction": "(Rested - Naive) / (Stimulated - Naive), centroid z units",
        },
        "label_genes_requested": SELECTED_GENES,
        "label_genes_missing": missing,
    }
    (TABLES / f"summary_Fig1F{SUFFIX}.json").write_text(json.dumps(summary, indent=2),
                                                        encoding="utf-8")
    print("\ndone.")


if __name__ == "__main__":
    main()
