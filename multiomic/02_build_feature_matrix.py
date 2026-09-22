#!/usr/bin/env python3
"""Stage 02 -- build the 11-feature matrix for the universe.

Features (11):
    RNA_Stim_vs_Naive, RNA_Rested_vs_Naive, RNA_Rechallenged_vs_Stim
    K27me3_TSS3kb_{Stim,Rested,Rechallenged}_vs_Naive
    K4me3_TSS3kb_{Stim,Rested,Rechallenged}_vs_Naive
    K36me2_TSS3kb_Stim_vs_Naive
    K36me3_geneBody_Stim_vs_Naive

Chromatin recipe:
    mean bigWig signal over the window (nanmean, 0.0 if empty/all-NaN)
    -> log2(sig + 0.01)
    -> subtract the Naive value of the same mark
Windows: TSS +/- 3 kb for K27me3/K4me3/K36me2; gene body for K36me3.

Output: data/feature_matrix.tsv, data/feature_completeness.tsv, data/raw_signal.tsv,
        data/feature_build_meta.json
"""
import sys, json, time, warnings
from pathlib import Path
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
warnings.filterwarnings("ignore", category=RuntimeWarning)
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import common as cm

DATA = Path(__file__).resolve().parent / "data"

RNA_SOURCE = "tpm"
K27_NORM = "cpm"
TPM_PSEUDO = 0.1                    # RNA ratios ONLY; the memory classes are untouched
PSEUDO = 0.01                       # log2(signal + 0.01)
TSS_HALFWIDTH = 3000

# mark -> (conditions available, window kind)
MARKS = {
    "H3K27me3": (["Naive", "Stim", "Rested", "Rechallenged"], "tss"),
    "H3K4me3":  (["Naive", "Stim", "Rested", "Rechallenged"], "tss"),
    "K36me2":   (["Naive", "Stim"], "tss"),
    "K36me3":   (["Naive", "Stim"], "body"),
}
FEATURES = ["RNA_Stim_vs_Naive", "RNA_Rested_vs_Naive", "RNA_Rechallenged_vs_Stim",
            "K27me3_TSS3kb_Stim_vs_Naive", "K27me3_TSS3kb_Rested_vs_Naive",
            "K27me3_TSS3kb_Rechallenged_vs_Naive",
            "K4me3_TSS3kb_Stim_vs_Naive", "K4me3_TSS3kb_Rested_vs_Naive",
            "K4me3_TSS3kb_Rechallenged_vs_Naive",
            "K36me2_TSS3kb_Stim_vs_Naive", "K36me3_geneBody_Stim_vs_Naive"]
CHROM_FEATURES = [f for f in FEATURES if not f.startswith("RNA_")]


def mean_sig(bw, chrom, s, e):
    try:
        v = bw.values(chrom, int(max(0, s)), int(e))
    except Exception:
        return np.nan          # NaN here means "could not read", distinct from a real 0
    if v is None or len(v) == 0:
        return np.nan
    m = np.nanmean(v)
    return float(m) if not np.isnan(m) else 0.0


def extract_signal(co):
    """Raw mean signal per (mark, condition) for every gene. Cached to data/raw_signal.tsv."""
    import pybigtools, hashlib
    cache = DATA / "raw_signal.tsv"
    stamp = DATA / ".raw_signal.key"
    # Key the cache on the COORDINATES, not the gene names, so a coordinate change invalidates
    # it rather than silently reusing pre-fix signal.
    key = hashlib.sha256(
        (co[["chrom", "tss", "body_start", "body_end", "strand"]].to_csv()
         + f"|halfwidth={TSS_HALFWIDTH}|pseudo={PSEUDO}").encode()).hexdigest()
    if cache.exists() and stamp.exists() and stamp.read_text().strip() == key:
        print(f"  raw signal read from cache ({cache.name}; coordinate key matches)")
        return pd.read_csv(cache, sep="\t", index_col=0)
    if cache.exists():
        print("  cache is STALE (the coordinates changed since it was written); recomputing")

    # sort by locus so each bigWig is read in roughly increasing file order
    order = co.sort_values(["chrom", "body_start"]).index
    sig = pd.DataFrame(index=co.index, dtype=float)
    for mark, (conds, kind) in MARKS.items():
        for cond in conds:
            t0 = time.time()
            bw = pybigtools.open(str(cm.bigwig(f"{cond}_{mark}_merged_cpm.bw")))
            col, vals = f"{mark}|{cond}", {}
            for g in order:
                r = co.loc[g]
                if kind == "body":
                    s, e = r.body_start, r.body_end
                else:
                    s, e = r.tss - TSS_HALFWIDTH, r.tss + TSS_HALFWIDTH
                vals[g] = mean_sig(bw, r.chrom, s, e)
            bw.close()
            sig[col] = pd.Series(vals)
            print(f"  {col:<24} {time.time()-t0:6.1f}s  "
                  f"zero {int((sig[col] == 0).sum()):>4}  unreadable {int(sig[col].isna().sum()):>3}")
    sig.to_csv(cache, sep="\t", index_label="gene")
    stamp.write_text(key)
    return sig


def build_chromatin(sig):
    lg = np.log2(sig + PSEUDO)
    F = pd.DataFrame(index=sig.index)
    for cond in ["Stim", "Rested", "Rechallenged"]:
        F[f"K27me3_TSS3kb_{cond}_vs_Naive"] = lg[f"H3K27me3|{cond}"] - lg["H3K27me3|Naive"]
        F[f"K4me3_TSS3kb_{cond}_vs_Naive"] = lg[f"H3K4me3|{cond}"] - lg["H3K4me3|Naive"]
    F["K36me2_TSS3kb_Stim_vs_Naive"] = lg["K36me2|Stim"] - lg["K36me2|Naive"]
    F["K36me3_geneBody_Stim_vs_Naive"] = lg["K36me3|Stim"] - lg["K36me3|Naive"]
    return F


def build_rna(u):
    p = TPM_PSEUDO
    F = pd.DataFrame(index=u.index)
    F["RNA_Stim_vs_Naive"] = np.log2((u.Stim + p) / (u.Naive + p))
    F["RNA_Rested_vs_Naive"] = np.log2((u.Rested + p) / (u.Naive + p))
    F["RNA_Rechallenged_vs_Stim"] = np.log2((u.Rechallenged + p) / (u.Stim + p))
    return F


def main():
    u = pd.read_csv(DATA / "universe.tsv", sep="\t").set_index("gene")
    co = pd.read_csv(DATA / "coordinates.tsv", sep="\t").set_index("gene")
    print(f"universe {len(u)}   coordinates {len(co)}")
    missing_coords = sorted(set(u.index) - set(co.index))
    print(f"universe genes without coordinates: {len(missing_coords)} {missing_coords}")
    co = co.loc[[g for g in u.index if g in co.index]]

    print(f"\nRNA features from {RNA_SOURCE} (pseudocount {TPM_PSEUDO} on the ratios only)")
    # A zero NUMERATOR is exactly as load-bearing as a zero denominator: log2(0/x) = -inf,
    # so the flag tests all four conditions, not just Naive and Stim.
    zero_any = u[(u[["Naive", "Stim", "Rested", "Rechallenged"]] == 0).any(axis=1)]
    print(f"  genes with a zero TPM in any condition (pseudocount LOAD-BEARING, the feature "
          f"would be +/-inf without it): {len(zero_any)}")
    if len(zero_any):
        print(zero_any[["Naive", "Stim", "Rested", "Rechallenged", "memory_class"]]
              .round(2).to_string())
    Frna = build_rna(u)
    # separately: how far does the pseudocount move a feature that was finite anyway?
    raw = pd.DataFrame({
        "RNA_Stim_vs_Naive": np.log2(u.Stim / u.Naive),
        "RNA_Rested_vs_Naive": np.log2(u.Rested / u.Naive),
        "RNA_Rechallenged_vs_Stim": np.log2(u.Rechallenged / u.Stim)})
    shift = (Frna - raw).abs()
    shift_finite = shift[np.isfinite(raw).all(axis=1)]
    n_shift_01 = int((shift_finite > 0.1).any(axis=1).sum())
    n_shift_05 = int((shift_finite > 0.5).any(axis=1).sum())
    print(f"  genes whose feature value it MATERIALLY SHIFTS (finite either way): "
          f"{n_shift_01} by >0.1, {n_shift_05} by >0.5")
    top = shift_finite.max(axis=1).sort_values(ascending=False).head(8)
    print("    largest shifts: " + ", ".join(f"{g} {v:.2f}" for g, v in top.items()))
    print("    (the compression sits at the top of the induction axis, on genes with a very "
          "low Naive TPM, i.e. on PC1's dominant feature)")

    print(f"\nchromatin signal ({K27_NORM} bigWigs), {len(co)} genes x 12 tracks")
    sig = extract_signal(co)
    Fchr = build_chromatin(sig)

    F = Frna.join(Fchr, how="left")[FEATURES]
    F.index.name = "gene"

    # completeness, per gene, every reason written down
    comp = pd.DataFrame(index=u.index)
    comp["memory_class"] = u.memory_class
    comp["has_coordinate"] = comp.index.isin(co.index)
    comp["coordinate_source"] = co.coordinate_source.reindex(comp.index).fillna("")
    comp["n_features_finite"] = F.reindex(comp.index).notna().sum(axis=1)
    comp["complete"] = comp.n_features_finite == len(FEATURES)
    comp["rna_pseudocount_used"] = comp.index.isin(zero_any.index)
    comp["short_model_lt1kb"] = co.short_model.reindex(comp.index).fillna(False)
    # a feature computed from two all-zero tracks is 0.0 by arithmetic, not by measurement
    flat = pd.Series(0, index=comp.index)
    one_sided = pd.Series(0, index=comp.index)
    below = pd.Series(0, index=comp.index)
    for mark, (conds, _) in MARKS.items():
        for cond in conds[1:]:
            a, b = sig[f"{mark}|Naive"], sig[f"{mark}|{cond}"]
            both0 = (a == 0) & (b == 0)
            # exactly one track zero: the feature is a pure pseudocount step, and these are
            # the most extreme values in the K36me3 column
            one0 = ((a == 0) ^ (b == 0))
            # both tracks under the signal pseudocount: the feature is floor-dominated
            bl = (a < PSEUDO) & (b < PSEUDO)
            for tgt, m in ((flat, both0), (one_sided, one0), (below, bl)):
                tgt += m.reindex(comp.index).fillna(False).astype(int)
    comp["n_features_from_zero_signal"] = flat.astype(int)
    comp["n_features_from_onesided_zero"] = one_sided.astype(int)
    comp["n_features_below_signal_floor"] = below.astype(int)
    # a gene with no promoter mark at all in any condition: its 6 TSS-window K27me3/K4me3
    # features are floor noise rather than measurement
    k4max = sig[[f"H3K4me3|{c}" for c in MARKS["H3K4me3"][0]]].max(axis=1)
    comp["no_promoter_k4me3_any_condition"] = (k4max < 0.05).reindex(comp.index).fillna(False)
    comp["all_12_tracks_zero"] = (sig == 0).all(axis=1).reindex(comp.index).fillna(False)
    reason = []
    for g in comp.index:
        if not comp.at[g, "has_coordinate"]:
            reason.append("no coordinate in any source")
        elif not comp.at[g, "complete"]:
            miss = F.loc[g][F.loc[g].isna()].index.tolist()
            reason.append("non-finite: " + ",".join(miss))
        else:
            reason.append("")
    comp["exclusion_reason"] = reason
    comp.to_csv(DATA / "feature_completeness.tsv", sep="\t", index_label="gene")

    Fc = F[comp.complete.reindex(F.index).fillna(False)]
    Fc.to_csv(DATA / "feature_matrix.tsv", sep="\t")

    print("\n--- RECONCILIATION, universe -> feature matrix ---")
    print(f"  universe                       {len(u)}")
    print(f"  with a coordinate              {int(comp.has_coordinate.sum())}"
          f"   (lost {int((~comp.has_coordinate).sum())})")
    print(f"  with all 11 features finite    {int(comp.complete.sum())}"
          f"   (lost {int((comp.has_coordinate & ~comp.complete).sum())})")
    lost = comp[~comp.complete]
    if len(lost):
        print("\n  EVERY GENE LOST, WITH ITS REASON:")
        print(lost[["memory_class", "exclusion_reason"]].to_string())
    else:
        print(f"\n  NO GENE LOST. universe == feature matrix == {len(Fc)}.")
    print("\n  by class:")
    print(pd.crosstab(comp.memory_class, comp.complete).to_string())

    print(f"\n  genes with >=1 feature computed from two all-zero tracks (FLAG, not a drop): "
          f"{int((comp.n_features_from_zero_signal > 0).sum())}")
    print(f"  genes with a <1 kb model (K36me3 gene-body caveat, FLAG not a drop): "
          f"{int(comp.short_model_lt1kb.sum())}")

    meta = {
        "rna_source": RNA_SOURCE, "k27_norm": K27_NORM, "tpm_pseudo": TPM_PSEUDO,
        "signal_pseudo": PSEUDO, "tss_halfwidth": TSS_HALFWIDTH,
        "n_universe": int(len(u)), "n_with_coordinates": int(comp.has_coordinate.sum()),
        "n_complete": int(comp.complete.sum()),
        "n_lost": int((~comp.complete).sum()),
        "lost_genes": {g: comp.at[g, "exclusion_reason"] for g in comp[~comp.complete].index},
        "by_class_complete": {k: int(v) for k, v in
                              comp[comp.complete].memory_class.value_counts().items()},
        "n_rna_pseudocount_loadbearing": int(comp.rna_pseudocount_used.sum()),
        "rna_pseudocount_loadbearing_genes": sorted(comp[comp.rna_pseudocount_used].index.tolist()),
        "n_rna_materially_shifted_gt_0.1": n_shift_01,
        "n_rna_materially_shifted_gt_0.5": n_shift_05,
        "n_short_models": int(comp.short_model_lt1kb.sum()),
        "short_models_by_class": {c: int(comp[comp.memory_class == c].short_model_lt1kb.sum())
                                  for c in ["memory", "refractory", "non_memory"]},
        "n_with_zero_signal_features": int((comp.n_features_from_zero_signal > 0).sum()),
        "n_with_onesided_zero_features": int((comp.n_features_from_onesided_zero > 0).sum()),
        "n_with_below_floor_features": int((comp.n_features_below_signal_floor > 0).sum()),
        "n_no_promoter_k4me3": int(comp.no_promoter_k4me3_any_condition.sum()),
        "genes_all_12_tracks_zero": sorted(comp[comp.all_12_tracks_zero].index.tolist()),
    }
    (DATA / "feature_build_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(f"\nwrote {DATA/'feature_matrix.tsv'}  ({Fc.shape[0]} x {Fc.shape[1]})")


if __name__ == "__main__":
    main()
