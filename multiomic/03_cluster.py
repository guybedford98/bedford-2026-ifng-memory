#!/usr/bin/env python3
"""Stage 03 -- PC-count selection, PCA + GMM clustering, memory projection.

    per-feature z-score
    -> PCA(n_components=N_PCS, random_state=42)
    -> GaussianMixture(covariance_type='full', n_init=10, init_params='kmeans',
                       random_state=42, max_iter=500)
    -> k sweep 2..10, BIC-best, relabel clusters by size (largest = cluster_0)

Output: data/gene_clusters_memory.tsv, gene_clusters_bicbest.tsv, gene_clusters_k2..10.tsv,
        pc_values.tsv, pc_loadings.tsv, variance_explained.tsv, per_view_pc_contribution.tsv,
        gmm_bic_sweep.tsv, anchor_placements_by_k.tsv, cluster_summary_memory.tsv,
        n_pcs_selection.tsv, run_meta.json
"""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
from sklearn.decomposition import PCA
from sklearn.mixture import GaussianMixture
from scipy.stats import hypergeom

sys.stdout.reconfigure(encoding="utf-8")

DATA = Path(__file__).resolve().parent / "data"

SEED = 42
ANCHORS = ["IDO1", "CXCL11", "IRF1", "STAT1"]
N_PERM = 500          # Horn's parallel analysis permutations
PA_PCT = 95           # keep components above this percentile of the null

# Report the RNA-priming component as PC2 and the H3K27me3-leaning one as PC3, i.e. swap the
# two. Legitimate only because they are near-degenerate on this universe (12.068% vs 12.378%,
# and they swap in 37% of gene bootstraps), so their variance ordering carries no information.
# The GMM is affine-equivariant, so this cannot change the clustering.
SWAP_PC2_PC3 = True


def parallel_analysis(Xz, n_perm=N_PERM, pct=PA_PCT, seed=SEED):
    """Horn's parallel analysis by column-wise permutation (correlation PCA)."""
    obs = PCA(random_state=seed).fit(Xz).explained_variance_
    rng = np.random.default_rng(seed)
    A = Xz.copy()
    null = np.empty((n_perm, A.shape[1]))
    for i in range(n_perm):
        P = np.column_stack([rng.permutation(A[:, j]) for j in range(A.shape[1])])
        P = (P - P.mean(0)) / np.where(P.std(0) == 0, 1.0, P.std(0))
        null[i] = PCA(random_state=seed).fit(P).explained_variance_
    thr = np.percentile(null, pct, axis=0)
    keep = 0
    for i in range(len(obs)):                 # leading run only, stop at the first failure
        if obs[i] > thr[i]:
            keep += 1
        else:
            break
    return obs, thr, max(1, keep)


def relabel_by_size(labels):
    vals, counts = np.unique(labels, return_counts=True)
    order = list(vals[np.argsort(-counts, kind="stable")])
    remap = {old: i for i, old in enumerate(order)}
    return np.array([remap[l] for l in labels])


def run_pipeline(X, n_pcs):
    std = X.std(axis=0).replace(0, 1.0)
    Xz = (X - X.mean(axis=0)) / std
    pcnames = [f"PC{i+1}" for i in range(n_pcs)]
    pca = PCA(n_components=n_pcs, random_state=SEED)
    Z = pd.DataFrame(pca.fit_transform(Xz.values), index=X.index, columns=pcnames)
    loadings = pd.DataFrame(pca.components_.T, index=X.columns, columns=pcnames)
    ve = pca.explained_variance_ratio_
    if SWAP_PC2_PC3 and n_pcs >= 3:
        # relabel components 2 and 3: reorder the columns, then re-apply the PC1..PCn names
        order = [0, 2, 1] + list(range(3, n_pcs))
        Z = Z.iloc[:, order]; Z.columns = pcnames
        loadings = loadings.iloc[:, order]; loadings.columns = pcnames
        ve = ve[order]
        print("  PC2/PC3 SWAPPED for reporting: what PCA ranks 3rd is reported as PC2. "
              "Legitimate only because they are near-degenerate here.")
    print(f"  variance explained: " +
          "  ".join(f"{p} {v*100:.2f}%" for p, v in zip(pcnames, ve)))

    VIEWS = {"RNA": [c for c in X.columns if c.startswith("RNA")],
             "H3K27me3": [c for c in X.columns if c.startswith("K27me3")],
             "H3K4me3": [c for c in X.columns if c.startswith("K4me3")],
             "H3K36me2": [c for c in X.columns if c.startswith("K36me2")],
             "H3K36me3": [c for c in X.columns if c.startswith("K36me3")]}
    contrib = {}
    for pc in pcnames:
        sq = loadings[pc] ** 2
        contrib[pc] = {v: float(sq[cols].sum() / sq.sum()) for v, cols in VIEWS.items()}
    contrib_df = pd.DataFrame(contrib)

    def fit_gmm(k):
        g = GaussianMixture(n_components=k, covariance_type="full", n_init=10,
                            init_params="kmeans", random_state=SEED, max_iter=500)
        g.fit(Z.values)
        lab = relabel_by_size(g.predict(Z.values))
        post = g.predict_proba(Z.values).max(axis=1)
        return lab, post, float(g.bic(Z.values)), float(g.aic(Z.values))

    sweep, kclust = [], {}
    for k in range(2, 11):
        lab, post, bic, aic = fit_gmm(k)
        df = Z.copy()
        df["cluster"] = [f"cluster_{c}" for c in lab]
        df["posterior"] = post
        kclust[k] = df
        df.to_csv(DATA / f"gene_clusters_k{k}.tsv", sep="\t")
        sweep.append({"k": k, "BIC": bic, "AIC": aic, "n_clusters": len(set(lab))})
    sweep = pd.DataFrame(sweep)
    bic_best_k = int(sweep.loc[sweep.BIC.idxmin(), "k"])
    clusters = kclust[bic_best_k].copy()

    clusters.to_csv(DATA / "gene_clusters_bicbest.tsv", sep="\t")
    Z.to_csv(DATA / "pc_values.tsv", sep="\t")
    loadings.to_csv(DATA / "pc_loadings.tsv", sep="\t")
    sweep.to_csv(DATA / "gmm_bic_sweep.tsv", sep="\t", index=False)
    contrib_df.to_csv(DATA / "per_view_pc_contribution.tsv", sep="\t")
    pd.Series(ve, index=pcnames, name="variance_explained").to_csv(
        DATA / "variance_explained.tsv", sep="\t")
    print(f"  BIC-best k = {bic_best_k}")
    return dict(Z=Z, loadings=loadings, ve=ve, contrib=contrib_df, sweep=sweep,
                kclust=kclust, clusters=clusters, bic_best_k=bic_best_k, pcnames=pcnames)


def project_memory(clusters, u, kclust):
    clusters = clusters.copy()
    clusters["memory_class"] = u.memory_class.reindex(clusters.index)
    CL = sorted(clusters.cluster.unique(), key=lambda c: int(c.split("_")[1]))
    M = len(clusters)
    n_mem = int((clusters.memory_class == "memory").sum())
    rows = []
    for cc in CL:
        g = clusters[clusters.cluster == cc]
        N = len(g)
        k = int((g.memory_class == "memory").sum())
        p = hypergeom.sf(k - 1, M, n_mem, N) if k > 0 else 1.0
        row = {"cluster": cc, "n_genes": N,
               "mean_posterior": round(g.posterior.mean(), 3)}
        for pc in [c for c in g.columns if c.startswith("PC")]:
            row[f"{pc}_mean"] = round(g[pc].mean(), 2)
        row.update({"n_memory": k, "pct_memory": round(100 * k / N, 1),
                    "n_refractory": int((g.memory_class == "refractory").sum()),
                    "fold_enrich": round((k / N) / (n_mem / M), 2) if k else 0.0,
                    "hypergeom_p": p,
                    "anchors": ",".join(a for a in ANCHORS
                                        if a in g.index)})
        rows.append(row)
    summary = pd.DataFrame(rows)
    summary.to_csv(DATA / "cluster_summary_memory.tsv", sep="\t", index=False)

    outg = clusters.join(u[["Naive", "Stim", "Rested", "Rechallenged", "maxTPM",
                            "memory_l2fc", "isg_l2fc"]])
    outg["Naive_TPM"] = outg["Naive"]           # names the panel scripts expect
    outg["Stim_TPM"] = outg["Stim"]
    outg["is_anchor"] = [g if g in ANCHORS else "" for g in outg.index]
    outg.reset_index().rename(columns={"index": "gene"}).to_csv(
        DATA / "gene_clusters_memory.tsv", sep="\t", index=False)

    apr = []
    for k in range(2, 11):
        row = {"k": k}
        for a in ANCHORS:
            row[a] = kclust[k].loc[a, "cluster"] if a in kclust[k].index else "NA"
        row["IDO1_CXCL11_together"] = (row["IDO1"] == row["CXCL11"]
                                       if row["IDO1"] != "NA" and row["CXCL11"] != "NA"
                                       else None)
        apr.append(row)
    apr = pd.DataFrame(apr)
    apr.to_csv(DATA / "anchor_placements_by_k.tsv", sep="\t", index=False)
    return clusters, summary, outg, apr


def main():
    u = pd.read_csv(DATA / "universe.tsv", sep="\t").set_index("gene")
    X = pd.read_csv(DATA / "feature_matrix.tsv", sep="\t").set_index("gene")
    comp = pd.read_csv(DATA / "feature_completeness.tsv", sep="\t").set_index("gene")

    print("--- INPUT RECONCILIATION ---")
    print(f"  universe.tsv            {len(u)}")
    print(f"  feature_matrix.tsv      {len(X)}")
    lost = sorted(set(u.index) - set(X.index))
    print(f"  in universe, not in matrix: {len(lost)}  {lost}")
    for g in lost:
        print(f"     {g}: {comp.at[g, 'exclusion_reason']}")
    extra = sorted(set(X.index) - set(u.index))
    assert not extra, f"matrix has genes outside the universe: {extra}"
    print(f"  -> clustering {len(X)} genes. No further filtering happens below this line.")
    print("  by class: " + "  ".join(
        f"{c} {int((u.memory_class.reindex(X.index) == c).sum())}"
        for c in ["memory", "refractory", "non_memory"]))

    print("\n--- PC-COUNT SELECTION (Horn's parallel analysis, "
          f"{N_PERM} permutations, {PA_PCT}th pct, seed {SEED}) ---")
    std = X.std(axis=0).replace(0, 1.0)
    Xz = ((X - X.mean(axis=0)) / std).values
    obs, thr, n_pcs = parallel_analysis(Xz)
    pa = pd.DataFrame({"component": np.arange(1, len(obs) + 1),
                       "observed_eigenvalue": obs, "null_p95": thr,
                       "retained": obs > thr})
    pa.to_csv(DATA / "n_pcs_selection.tsv", sep="\t", index=False)
    print(pa.round(4).to_string(index=False))
    print(f"\n  parallel analysis retains {n_pcs} component(s)")
    margin = float(pa.observed_eigenvalue.iloc[2] - pa.null_p95.iloc[2])
    rel = margin / float(pa.null_p95.iloc[2])
    print(f"  component-3 margin: {margin:+.4f} ({100*rel:+.1f}% of the null threshold)")
    if abs(rel) < 0.10:
        print("  *** That margin is inside 10% of the threshold, so the PC count is NOT "
              "robustly determined. ***")

    print(f"\n--- RUN, n_pcs = {n_pcs} ---")
    P = run_pipeline(X, n_pcs)
    clusters, summary, outg, apr = project_memory(P["clusters"], u, P["kclust"])

    print("\n--- CLUSTER SUMMARY (display numbers are cluster+1) ---")
    print(summary.to_string(index=False))

    print("\n--- ANCHOR GATE ---")
    print(apr.to_string(index=False))
    together = apr.IDO1_CXCL11_together
    all_k = bool(together.all())
    at_best = bool(apr.loc[apr.k == P["bic_best_k"], "IDO1_CXCL11_together"].iloc[0])
    if all_k:
        print("\n  IDO1 and CXCL11 co-cluster at EVERY k (2..10). Gate PASSED.")
    else:
        bad = apr.loc[~together.astype(bool), "k"].tolist()
        print(f"\n  *** ANCHOR GATE: IDO1 and CXCL11 do NOT co-cluster at k = {bad}. "
              f"At the BIC-best k = {P['bic_best_k']} they are "
              f"{'together' if at_best else 'APART'}. ***")
    print("\n  anchor clusters at BIC-best k (display numbers):")
    for a in ANCHORS:
        if a in clusters.index:
            cl = clusters.at[a, "cluster"]
            print(f"    {a:<7} cluster {int(cl.split('_')[1])+1}  "
                  f"({u.at[a, 'memory_class']})")

    meta = {
        "method": ("z-score -> PCA(n_pcs, seed42) -> GMM(full, n_init10, kmeans, seed42), "
                   "k 2-10, BIC-best, relabel by size"),
        "seed": SEED,
        "rna_source": "tpm", "k27_norm": "cpm", "coordinate_source_mode": "union_local",
        "n_universe": int(len(u)), "n_clustered": int(len(X)),
        "n_lost_universe_to_matrix": len(lost),
        "lost_genes": {g: comp.at[g, "exclusion_reason"] for g in lost},
        "n_by_class_clustered": {c: int((u.memory_class.reindex(X.index) == c).sum())
                                 for c in ["memory", "refractory", "non_memory"]},
        "n_by_class_universe": {c: int((u.memory_class == c).sum())
                                for c in ["memory", "refractory", "non_memory"]},
        "n_pcs": int(n_pcs),
        "n_pcs_selection": {
            "method": f"Horn parallel analysis, {N_PERM} column permutations, "
                      f"{PA_PCT}th percentile",
            "retained_by_parallel_analysis": int(n_pcs),
            "component3_observed": float(pa.observed_eigenvalue.iloc[2]),
            "component3_null_p95": float(pa.null_p95.iloc[2]),
            "component3_margin": margin,
            "component3_margin_pct_of_threshold": round(100 * rel, 2)},
        "swap_pc2_pc3": SWAP_PC2_PC3,
        "variance_explained": {p: float(v) for p, v in zip(P["pcnames"], P["ve"])},
        "per_view_pc_contribution": P["contrib"].to_dict(),
        "bic_best_k": int(P["bic_best_k"]),
        "cluster_sizes": summary.n_genes.tolist(),
        "pct_memory_by_cluster": summary.pct_memory.tolist(),
        "hypergeom_p_by_cluster": summary.hypergeom_p.tolist(),
        "anchor_clusters_display": {a: int(clusters.at[a, "cluster"].split("_")[1]) + 1
                                    for a in ANCHORS if a in clusters.index},
        "anchor_gate_ido1_cxcl11_together_at_every_k": all_k,
        "anchor_gate_at_bic_best_k": at_best,
        "universe_csv": json.loads(
            (DATA / "universe_summary.json").read_text(encoding="utf-8")).get("source", "unknown"),
    }
    (DATA / "run_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"\nwrote {DATA/'run_meta.json'}")


if __name__ == "__main__":
    main()
