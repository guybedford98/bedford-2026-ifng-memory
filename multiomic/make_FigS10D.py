"""Fig S10D: memory vs non-memory within a cluster, the four chromatin marks.
Built for EVERY cluster in the run; the manuscript panel uses clusters 3 and 4.

Refractory genes are excluded (this is a memory vs non-memory contrast), so the column n values
are memory and non-memory only, and the implied "% memory" has a different denominator than the
Fig 5C composition pies'.

Coordinates come from coordinates.tsv (0-based half-open, body = [body_start, body_end),
tss = body_start if '+' else body_end-1, TSS window = [tss-3000, tss+3000)). bigWig extraction
is cached to data/memVnon_metaplot_profiles.npz; delete the file to force a recompute.

Output: results/FigS10D/FigS10D_cluster{N}.{pdf,png,svg} + _n.json
"""
import sys, json
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")
import numpy as np, pandas as pd, pybigtools
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import common as cm
import figstyle as fs

fs.setup()
np.random.seed(42)

DATA = Path(__file__).resolve().parent / "data"
OUT = cfg.out_dir("FigS10D")
CACHE = DATA / "memVnon_metaplot_profiles.npz"
COND_COLORS = {"Naive": fs.COLOR["naive"], "Stim": fs.COLOR["stimulated"],
               "Rested": fs.COLOR["rested"], "Rechallenged": fs.COLOR["rechallenge"]}

CONDS = ["Naive", "Stim", "Rested", "Rechallenged"]
HALF, BIN = 3000, 100
NBIN = 2 * HALF // BIN
FLANK, UPB, BODYB, DNB = 1000, 10, 40, 10
TSS_MARKS = [("H3K27me3", "H3K27me3", CONDS), ("H3K4me3", "H3K4me3", CONDS),
             ("H3K36me2", "K36me2", ["Naive", "Stim"])]
GB_MARK = ("H3K36me3", "K36me3", ["Naive", "Stim"])
GRP_COLOR = {"memory": "#b2182b", "non-memory": "#7a7a7a"}
GRPS = ["memory", "non-memory"]
SMALL_N = 10                       # below this a group mean is not a group profile; say so
cl_all = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t")
coords_tab = pd.read_csv(DATA / "coordinates.tsv", sep="\t")

# every cluster in this run, in display order (1..k), rather than a hard-coded pair
DISPLAY_CLUSTERS = sorted(int(c.split("_")[1]) + 1 for c in cl_all.cluster.unique())
UNIVERSE_NOTE = (f"ISG universe ({len(cl_all):,} clustered genes); "
                 "refractory genes excluded from this panel")


def groups_for(disp_cluster):
    key = f"cluster_{disp_cluster - 1}"            # data files are 0-indexed
    sub = cl_all[cl_all.cluster == key]
    return sub, {"memory": set(sub[sub.memory_class == "memory"].gene),
                 "non-memory": set(sub[sub.memory_class == "non_memory"].gene)}


def binned(arr, n):
    idx = np.linspace(0, len(arr), n + 1).astype(int)
    return np.array([np.nanmean(arr[idx[i]:idx[i+1]]) if idx[i+1] > idx[i] else np.nan for i in range(n)])
def tss_prof(bw, chrom, t, strand):
    s, e = t - HALF, t + HALF
    if s < 0: return None
    try: v = np.asarray(bw.values(chrom, s, e), dtype=float)
    except Exception: return None
    if v.size != 2 * HALF: return None
    p = np.nanmean(v.reshape(NBIN, BIN), axis=1); return p[::-1] if strand == "-" else p
def gene_scaled(bw, chrom, start, end, strand):
    if start - FLANK < 0 or end <= start: return None
    try:
        up = np.asarray(bw.values(chrom, start - FLANK, start), dtype=float)
        body = np.asarray(bw.values(chrom, start, end), dtype=float)
        dn = np.asarray(bw.values(chrom, end, end + FLANK), dtype=float)
    except Exception: return None
    p = np.concatenate([binned(up, UPB), binned(body, BODYB), binned(dn, DNB)]); return p[::-1] if strand == "-" else p


ALL_GENES = sorted(set().union(*[set().union(*groups_for(c)[1].values()) for c in DISPLAY_CLUSTERS]))
_c = coords_tab[coords_tab.gene.isin(ALL_GENES)]
TSS_COORD = {r.gene: (r.chrom, int(r.tss), r.strand) for r in _c.itertuples()}
GB_COORD = {r.gene: (r.chrom, int(r.body_start), int(r.body_end), r.strand) for r in _c.itertuples()}
print(f"genes needed (clusters {DISPLAY_CLUSTERS}, memory + non-memory): {len(ALL_GENES)}; "
      f"with coordinates: {len(TSS_COORD)}")
_missing = [g for g in ALL_GENES if g not in TSS_COORD]
if _missing:
    print("  NO COORDINATE:", _missing)


def extract_all():
    """{f'{mark}|{cond}': (nGene x nBin) array over ALL_GENES}; an all-NaN row = no profile."""
    store = {}
    jobs = [(disp, tok, conds, TSS_COORD, tss_prof, NBIN) for disp, tok, conds in TSS_MARKS]
    jobs.append((GB_MARK[0], GB_MARK[1], GB_MARK[2], GB_COORD, gene_scaled, UPB + BODYB + DNB))
    for disp, tok, conds, coords, fn, nb in jobs:
        for cond in conds:
            bw = pybigtools.open(str(cm.bigwig(f"{cond}_{tok}_merged_cpm.bw")))
            chroms0 = set(bw.chroms().keys())
            mat = np.full((len(ALL_GENES), nb), np.nan)
            ok = 0
            for i, g in enumerate(ALL_GENES):
                cd = coords.get(g)
                if cd is None or cd[0] not in chroms0:
                    continue
                p = fn(bw, *cd)
                if p is not None:
                    mat[i] = p; ok += 1
            bw.close()
            store[f"{disp}|{cond}"] = mat
            print(f"  {disp} {cond}: {ok}/{len(ALL_GENES)} genes with a profile")
    return store


_key = "\n".join(ALL_GENES)
PROF = None
if CACHE.exists():
    _z = np.load(CACHE, allow_pickle=False)
    if str(_z["key"]) == _key:
        PROF = {k: _z[k] for k in _z.files if k not in ("key", "genes")}
        print(f"loaded cached profiles <- {CACHE.name}")
    else:
        print("cache is stale (gene set changed); re-extracting")
if PROF is None:
    PROF = extract_all()
    np.savez_compressed(CACHE, key=np.array(_key), genes=np.array(ALL_GENES), **PROF)
    print(f"cached profiles -> {CACHE.name}")
GIDX = {g: i for i, g in enumerate(ALL_GENES)}


def build(disp_cluster):
    sub, groups = groups_for(disp_cluster)
    n_refr = int((sub.memory_class == "refractory").sum())
    print(f"cluster {disp_cluster}:", {k: len(v) for k, v in groups.items()},
          f"(refractory excluded: {n_refr})")

    def rows_with_profile(disp, cond, gs):
        m = PROF[f"{disp}|{cond}"][[GIDX[g] for g in gs if g in GIDX]]
        return m[~np.all(np.isnan(m), axis=1)]

    res, n_used = {}, {}
    for disp, tok, conds in TSS_MARKS + [GB_MARK]:
        res[disp] = {c: {grp: np.nanmean(rows_with_profile(disp, c, gs), axis=0)
                         for grp, gs in groups.items()} for c in conds}
        n_used[disp] = {grp: int(rows_with_profile(disp, conds[0], gs).shape[0])
                        for grp, gs in groups.items()}
    print(f"  done cluster {disp_cluster}")

    MARK_ROWS = [m[0] for m in TSS_MARKS] + [GB_MARK[0]]
    LAST_TSS_ROW = len(TSS_MARKS) - 1      # where the x unit label belongs: the gene-body row
                                           # is last and carries its own -1kb/TSS/TES/+1kb axis
    xax_tss = np.linspace(-3, 3, NBIN); xgb = np.arange(UPB + BODYB + DNB)
    nrow, ncol = len(MARK_ROWS), len(GRPS)
    fig, axs = plt.subplots(nrow, ncol, figsize=fs.panel_mm(ncol * 40, nrow * 22), squeeze=False)

    for i, disp in enumerate(MARK_ROWS):
        is_gb = (disp == GB_MARK[0])
        ymax = max(np.nanmax(res[disp][cond][g]) for cond in res[disp] for g in GRPS) * 1.06
        for j, grp in enumerate(GRPS):
            ax = axs[i][j]
            for cond in CONDS:
                if cond in res[disp]:
                    ax.plot(xgb if is_gb else xax_tss, res[disp][cond][grp],
                            color=COND_COLORS[cond], lw=1.1, zorder=3)
            ax.set_ylim(0, ymax)
            ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
            ax.tick_params(labelsize=6, length=2.0, width=0.8)
            if is_gb:
                for xv in (UPB, UPB + BODYB): ax.axvline(xv, color="#cccccc", lw=0.6, zorder=0)
                ax.set_xticks([0, UPB, UPB + BODYB, UPB + BODYB + DNB])
                if i == nrow - 1: ax.set_xticklabels(["-1kb", "TSS", "TES", "+1kb"], fontsize=6)
                else: ax.set_xticklabels([])
            else:
                ax.axvline(0, color="#cccccc", lw=0.6, zorder=0); ax.set_xticks([-3, 0, 3])
            if i == 0:
                ax.set_title(f"cluster {disp_cluster} {grp}\n(n={len(groups[grp])})",
                             fontsize=8, fontweight="bold", color=GRP_COLOR[grp])
            if j == 0: ax.set_ylabel(f"{disp}\nCPM", fontsize=8, fontweight="bold")
            else: ax.tick_params(labelleft=False)
            if i == LAST_TSS_ROW and not is_gb:
                ax.set_xlabel("kb from TSS", fontsize=8)

    handles = [Line2D([0], [0], color=COND_COLORS[c], lw=1.6, label=c) for c in CONDS]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.035), ncol=4,
               frameon=False, fontsize=8, handlelength=1.4, columnspacing=1.4)
    # Producing this for EVERY cluster means some panels average a handful of genes -- on the
    # memory side that can be 1-3. Those curves are not interpretable as group means, so say so
    # on the figure rather than letting the reader assume every panel carries equal weight.
    small = {g: len(v) for g, v in groups.items() if len(v) < SMALL_N}
    note = UNIVERSE_NOTE
    if small:
        note += ("\nCAUTION: " + ", ".join(f"{g} n={n}" for g, n in small.items()) +
                 f" -- fewer than {SMALL_N} genes, so this curve is a mean over very few "
                 "loci and should not be read as a group profile")
    fig.text(0.5, -0.012, note, ha="center", va="top", fontsize=8, color="#4a4a4a")
    # hspace 0.60: with no RNA row the x unit label falls on a MIDDLE row (the last TSS row),
    # not the bottom one, and at 0.42 the row below overdrew it.
    fig.subplots_adjust(left=0.135, right=0.99, top=0.88, bottom=0.075, hspace=0.60, wspace=0.18)
    stem = OUT / f"FigS10D_cluster{disp_cluster}"
    fs.save(fig, stem)

    meta = {
        "panel": f"FigS10D_cluster{disp_cluster}",
        "universe_n": int(len(cl_all)),
        "clustered_n": int(len(cl_all)),
        "display_cluster": disp_cluster,
        "internal_cluster": f"cluster_{disp_cluster - 1}",
        "cluster_n_total": int(len(sub)),
        "n_title_memory": len(groups["memory"]),
        "n_title_non_memory": len(groups["non-memory"]),
        "n_refractory_excluded": n_refr,
        "n_denominates": ("column-title n = genes of that memory class assigned to this cluster; "
                          "refractory excluded, so memory + non-memory < cluster total. Each "
                          "chromatin row averages the subset of those with an extractable "
                          "profile (n_with_chromatin_profile), which is what the curve shows."),
        "n_with_chromatin_profile": n_used,
        "small_n_caution": {g: len(v) for g, v in groups.items() if len(v) < SMALL_N},
    }
    (OUT / f"FigS10D_cluster{disp_cluster}_n.json").write_text(
        json.dumps(meta, indent=2), encoding="utf-8")
    print(f"  written -> FigS10D_cluster{disp_cluster} (+ _n.json)")


for c in DISPLAY_CLUSTERS:
    build(c)
