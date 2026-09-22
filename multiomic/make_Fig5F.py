"""Fig 5F: chromatin metaplot grid, mark (row) x cluster (column), conditions as line colour.

H3K27me3, H3K4me3 and H3K36me2 at TSS +/- 3 kb; H3K36me3 as a scaled gene body.
H3K36me2 / H3K36me3 bigWigs exist for Naive and Stim only.
Columns are clusters 1-4 (data files are 0-indexed cluster_0..3; display via figstyle.clabel()).
Per-gene profiles are cached to data/metaplot_profiles.npz so restyling never re-reads bigWigs.

Coordinates come from coordinates.tsv for every clustered gene (0-based half-open;
body = [body_start, body_end); tss = body_start if '+' else body_end-1;
TSS window = [tss-3000, tss+3000)).

Output: results/Fig5F/Fig5F.{pdf,png,svg}, data/metaplot_gene_counts.tsv
"""
import sys, hashlib
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
OUT = cfg.out_dir("Fig5F")
CLUSTER_COLORS = {1: "#20566b", 2: "#c9a227", 3: "#2e8b57", 4: "#b0532b"}
COND_COLORS = {"Naive": fs.COLOR["naive"], "Stim": fs.COLOR["stimulated"],
               "Rested": fs.COLOR["rested"], "Rechallenged": fs.COLOR["rechallenge"]}

CONDS = ["Naive", "Stim", "Rested", "Rechallenged"]
HALF, BIN = 3000, 100
NBIN = 2 * HALF // BIN
FLANK, UPB, BODYB, DNB = 1000, 10, 40, 10
TSS_MARKS = [("H3K27me3", "H3K27me3", CONDS), ("H3K4me3", "H3K4me3", CONDS),
             ("H3K36me2", "K36me2", ["Naive", "Stim"])]
GB_MARK = ("H3K36me3", "K36me3", ["Naive", "Stim"])
CACHE = DATA / "metaplot_profiles.npz"

cl = pd.read_csv(DATA / "gene_clusters_memory.tsv", sep="\t").set_index("gene")
CL = sorted(cl["cluster"].unique(), key=lambda x: int(x.split("_")[1]))
cl_map = cl["cluster"].to_dict()
print("clustered genes: %d  clusters displayed: %s" % (len(cl), [fs.clabel(c) for c in CL]))

co = pd.read_csv(DATA / "coordinates.tsv", sep="\t")
co = co[co["gene"].isin(cl_map)]
tss, gb = {}, {}
for r in co.itertuples(index=False):
    tss[r.gene] = (r.chrom, int(r.tss), r.strand)
    gb[r.gene] = (r.chrom, int(r.body_start), int(r.body_end), r.strand)
missing = sorted(set(cl_map) - set(tss))
print("coords found: TSS %d / gene body %d of %d clustered genes" % (len(tss), len(gb), len(cl)))
if missing:
    print("  clustered genes with no coordinate (%d): %s" % (len(missing), ", ".join(missing)))

genes_by_cl = {c: [g for g in tss if cl_map.get(g) == c] for c in CL}
for c in CL:
    print("  cluster %s: %d genes with TSS coords" % (fs.clabel(c), len(genes_by_cl[c])))


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


# The cache is keyed on a hash of the COORDINATE FRAME plus the window geometry, not on the
# gene names, so a coordinate change invalidates it instead of serving stale profiles.
_key_src = (co[["chrom", "tss", "body_start", "body_end", "strand"]].to_csv()
            + "|HALF=%d|BIN=%d|FLANK=%d|UPB=%d|BODYB=%d|DNB=%d"
            % (HALF, BIN, FLANK, UPB, BODYB, DNB))
CACHE_KEY = hashlib.sha256(_key_src.encode()).hexdigest()

cache = {}
if CACHE.exists():
    z = np.load(CACHE, allow_pickle=True)
    _stored = str(z["_coord_key"]) if "_coord_key" in z.files else None
    if _stored == CACHE_KEY:
        cache = {k: z[k] for k in z.files if k != "_coord_key"}
        print("loaded cached profiles: %s (%d mark x condition sets, coordinate key matches)"
              % (CACHE.name, len(cache) // 2))
    else:
        print("cache %s is STALE (coordinates or window geometry changed); recomputing" % CACHE.name)


def per_gene(disp, tok, cond, coords, fn):
    kp, kg = "%s|%s|prof" % (disp, cond), "%s|%s|genes" % (disp, cond)
    if kp in cache and kg in cache:
        return dict(zip(list(cache[kg]), cache[kp]))
    bw = pybigtools.open(str(cm.bigwig("%s_%s_merged_cpm.bw" % (cond, tok))))
    chroms0 = set(bw.chroms().keys()); per = {}
    for g, cd in coords.items():
        if cd[0] not in chroms0: continue
        p = fn(bw, *cd)
        if p is not None: per[g] = p
    bw.close()
    gl = sorted(per)
    cache[kg] = np.array(gl, dtype=object); cache[kp] = np.vstack([per[g] for g in gl])
    return per


def compute(disp, tok, conds, coords, fn):
    out, used = {}, {}
    for cond in conds:
        per = per_gene(disp, tok, cond, coords, fn)
        out[cond] = {c: np.nanmean([per[g] for g in genes_by_cl[c] if g in per], axis=0) for c in CL}
        used[cond] = {c: sum(1 for g in genes_by_cl[c] if g in per) for c in CL}
    return out, used


res, used_n = {}, {}
for disp, tok, conds in TSS_MARKS:
    res[disp], used_n[disp] = compute(disp, tok, conds, tss, tss_prof); print("done", disp)
res[GB_MARK[0]], used_n[GB_MARK[0]] = compute(GB_MARK[0], GB_MARK[1], GB_MARK[2], gb, gene_scaled)
print("done", GB_MARK[0])
np.savez_compressed(CACHE, _coord_key=np.array(CACHE_KEY), **cache)
print("cached per-gene profiles -> %s" % CACHE)

MARK_ROWS = [m[0] for m in TSS_MARKS] + [GB_MARK[0]]
xax_tss = np.linspace(-3, 3, NBIN); xgb = np.arange(UPB + BODYB + DNB)
nrow, ncol = len(MARK_ROWS), len(CL)
fig, axs = plt.subplots(nrow, ncol, figsize=fs.panel_mm(ncol * 33, nrow * 22), squeeze=False)

for i, disp in enumerate(MARK_ROWS):
    is_gb = (disp == GB_MARK[0]); conds_here = list(res[disp].keys())
    ymax = max(np.nanmax(res[disp][cond][c]) for cond in conds_here for c in CL) * 1.06
    for j, c in enumerate(CL):
        ax = axs[i][j]
        for cond in CONDS:
            if cond in res[disp]:
                ax.plot(xgb if is_gb else xax_tss, res[disp][cond][c], color=COND_COLORS[cond], lw=1.0, zorder=3)
        ax.set_ylim(0, ymax); ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
        ax.tick_params(labelsize=6, length=2.0, width=0.8)
        if i == 0:
            ax.set_title("cluster %s\n(n=%d)" % (fs.clabel(c), len(genes_by_cl[c])), fontsize=8,
                         fontweight="bold", color=CLUSTER_COLORS[int(fs.clabel(c))])
        if is_gb:
            for xv in (UPB, UPB + BODYB): ax.axvline(xv, color="#cccccc", lw=0.6, zorder=0)
            ax.set_xticks([0, UPB, UPB + BODYB, UPB + BODYB + DNB])
            if i == nrow - 1: ax.set_xticklabels(["-1kb", "TSS", "TES", "+1kb"], fontsize=6)
            else: ax.set_xticklabels([])
        else:
            ax.axvline(0, color="#cccccc", lw=0.6, zorder=0); ax.set_xticks([-3, 0, 3])
        if j == 0: ax.set_ylabel("%s\nCPM" % disp, fontsize=8, fontweight="bold")
        else: ax.tick_params(labelleft=False)
        # The x unit belongs on the bottom-most TSS row: MARK_ROWS ends with the gene-body mark,
        # which keeps its own scaled-body axis (-1kb / TSS / TES / +1kb). labelpad=0 because the
        # row gap is 14.9 pt and the tick labels already use 9.8 pt of it; the default 4 pt pad
        # drives the label into the gene-body axes below.
        if not is_gb and i == len(TSS_MARKS) - 1:
            ax.set_xlabel("kb from TSS", fontsize=8, labelpad=0)

handles = [Line2D([0], [0], color=COND_COLORS[c], lw=1.6, label=c) for c in CONDS]
fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.05), ncol=4, frameon=False,
           fontsize=8, handlelength=1.4, columnspacing=1.4)
fig.subplots_adjust(left=0.115, right=0.995, top=0.85, bottom=0.10, hspace=0.42, wspace=0.16)
CAPTION = ("n = clustered genes per cluster ({:,} clustered, all with resolvable coordinates). "
           "All memory classes included.\nLines are per-cluster mean profiles of merged CPM. "
           "H3K36me2 / H3K36me3 have Naive and Stim bigWigs only.").format(len(cl))
fig.text(0.5, -0.035, CAPTION, ha="center", va="top", fontsize=8, color="#444444", linespacing=1.5)
fs.save(fig, OUT / "Fig5F")

counts = pd.DataFrame({
    "cluster": [fs.clabel(c) for c in CL],
    "n_genes_with_TSS_coords": [len(genes_by_cl[c]) for c in CL],
    "n_genes_with_genebody_coords": [sum(1 for g in genes_by_cl[c] if g in gb) for c in CL],
    "n_profiled_TSS_min_over_conds": [min(used_n[d][cd][c] for d, _, cds in TSS_MARKS for cd in cds) for c in CL],
    "n_profiled_genebody_min_over_conds": [min(used_n[GB_MARK[0]][cd][c] for cd in GB_MARK[2]) for c in CL],
    "n_memory": [int((cl.loc[genes_by_cl[c], "memory_class"] == "memory").sum()) for c in CL],
    "n_refractory": [int((cl.loc[genes_by_cl[c], "memory_class"] == "refractory").sum()) for c in CL],
    "n_non_memory": [int((cl.loc[genes_by_cl[c], "memory_class"] == "non_memory").sum()) for c in CL]})
counts.to_csv(DATA / "metaplot_gene_counts.tsv", sep="\t", index=False)
print(counts.to_string(index=False))
print("Fig5F written (4 mark rows)")
