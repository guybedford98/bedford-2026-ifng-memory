#!/usr/bin/env python3
"""Stage 01 -- resolve a TSS and gene body for every universe gene.

Resolution order and the rule used inside each source
-----------------------------------------------------
1. gencode_v49_pc   GENCODE v49 protein-coding consensus (one transcript per gene).
                    Read from the BODY bed6 only -- see the TRAP below.
2. gencode_v44_gtf  GENCODE v44 full GTF, all biotypes. Transcript pick:
                    MANE_Select > Ensembl_canonical > longest transcript.
3. ucsc_refgene     UCSC refGene, resolved through shared/common.py and cached.
4. ncbi_eutils      data/coordinates_ncbi_supplement.tsv -- the 5 LOC/pseudogene cases that
                    no GENCODE/UCSC file carries.

TRAP -- the v49 "TSS" bed6 that ships beside the body bed6 is NOT strand-aware
------------------------------------------------------------------------------
It has TSS == body start for all 19,746 rows on BOTH strands. For a minus-strand gene that is
the 3' end, not the TSS: CXCL11 is off by 2,388 bp, IRF1 by 9,164 bp, STAT1 by 45,022 bp.
That file is therefore not read here. TSS is derived from the v49 BODY bed6, strand-aware.

Coordinate convention used everywhere downstream:
    0-based, half-open.  body = [body_start, body_end).  tss = body_start if '+'
    else body_end - 1.   TSS window = [tss - 3000, tss + 3000).

Output: data/coordinates.tsv, data/coordinates_failed.tsv, data/coordinates_summary.json,
        data/coordinates_outliers.tsv
"""
import sys, json, re
from pathlib import Path
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import common as cm

DATA = Path(__file__).resolve().parent / "data"

SOURCE_MODE = "union_local"

V49_BODY = cfg.path("reference", "protein_coding_bed6")
V44_GTF = cfg.path("reference", "gencode_gtf")
NCBI_SUPP = DATA / "coordinates_ncbi_supplement.tsv"

SHORT_MODEL_BP = 1000                # flag, never drop
BW_FOR_CONTIGS = cm.bigwig("Naive_H3K27me3_merged_cpm.bw")
V44_TSS_DISAGREE_BP = 10000          # flag a v49 pick this far from the v44 canonical model

# Deprecated HGNC symbols the universe still uses, mapped to the current symbol the v49 file
# carries. Without these six the alias-blind symbol match demotes them off the version-matched
# source that comes first. Verified by locus+strand overlap: every pair is the same model
# (reciprocal overlap >= 0.994, TSS delta 0).
V49_ALIASES = {"DDX58": "RIGI", "ARNTL2": "BMAL2", "BVES": "POPDC1",
               "C12orf4": "FERRY3", "KIAA1522": "NHSL3", "ODF3B": "CIMAP1B"}

# Symbols the v49 consensus file attaches to the WRONG model, demoted to v44.
# SCHIP1: the v49 bed6 gives SCHIP1 -> ENST00000638749.2, which in v44 belongs to gene
# IQCJ-SCHIP1 (ENSG00000283154.3, tag readthrough_gene). The same v49 file also carries a
# separate IQCJ-SCHIP1 row (ENST00000485419.7), so the assignment is internally inconsistent
# within that one file. Taking it puts the promoter window 566,617 bp from the SCHIP1 TSS
# that v44 uses, and inflates the K36me3 "gene body" from 57.5 kb to 624 kb. This is the same
# class of defect as the non-strand-aware TSS bed6 above: a bad input file, not a version
# difference. Caught automatically by the outlier guard below; listed here so the override is
# explicit rather than a silent rule.
V49_REJECT = {"SCHIP1": "v49 attaches this symbol to the IQCJ-SCHIP1 readthrough transcript"}

BED_COLS = ["chrom", "start", "end", "name", "score", "strand"]


def usable_contigs():
    """Contigs the bigWigs actually carry. A coordinate on anything else cannot be measured,
    so stage 01 must not accept it -- otherwise the source chain stops early on an unusable
    hit and the failure only surfaces at stage 02, past the layer that reports losses."""
    try:
        import pybigtools
        bw = pybigtools.open(str(BW_FOR_CONTIGS))
        c = set(bw.chroms())
        bw.close()
        return c
    except Exception as e:
        print(f"  (bigWig contig list unavailable: {e}; falling back to primary-assembly rule)")
        return None


def rec(chrom, body_start, body_end, strand, source, model_id=""):
    """One normalised coordinate record. 0-based half-open, strand-aware TSS."""
    body_start, body_end = int(body_start), int(body_end)
    if body_start > body_end:
        body_start, body_end = body_end, body_start
    tss = body_start if strand == "+" else body_end - 1
    return dict(chrom=chrom, body_start=body_start, body_end=body_end, strand=strand,
                tss=int(tss), model_length=int(body_end - body_start),
                coordinate_source=source, model_id=model_id)


def parse_v49_body(wanted):
    """GENCODE v49 protein-coding consensus, one transcript per gene.

    27 symbols are duplicated (mostly chrX/chrY pseudoautosomal pairs). Deterministic rule:
    drop the chrY copy when a non-chrY copy exists, then keep the longest span, then the
    lexicographically smallest transcript id. Every universe gene this rule touches is
    reported.
    """
    d = pd.read_csv(V49_BODY, sep="\t", header=None, names=BED_COLS)
    d["gene"] = d.name.str.split("|").str[0]
    d["tx"] = d.name.str.split("|").str[1]
    # The universe still uses six deprecated HGNC symbols. Matching on the symbol alone
    # demotes them off the version-matched source, so map the current v49 symbol back to the
    # alias the universe carries.
    rev = {v: k for k, v in V49_ALIASES.items() if k in wanted}
    d["gene"] = d.gene.replace(rev)
    d = d[d.gene.isin(wanted) & ~d.gene.isin(V49_REJECT)].copy()
    dup_hits = sorted(d.loc[d.gene.duplicated(keep=False), "gene"].unique().tolist())
    d["_y"] = (d.chrom == "chrY").astype(int)
    d["_span"] = d.end - d.start
    d = d.sort_values(["gene", "_y", "_span", "tx"], ascending=[True, True, False, True])
    d = d.drop_duplicates("gene", keep="first")
    out = {r.gene: rec(r.chrom, r.start, r.end, r.strand, "gencode_v49_pc", r.tx)
           for r in d.itertuples()}
    return out, dup_hits


ATTR = re.compile(r'(\S+) "([^"]*)"')


def _attrs(field):
    return dict(ATTR.findall(field))


def parse_v44_gtf(wanted):
    """v44 transcripts for `wanted`, picking MANE_Select > Ensembl_canonical > longest."""
    best = {}
    with open(V44_GTF, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "transcript":
                continue
            a = _attrs(f[8])
            g = a.get("gene_name")
            if g not in wanted:
                continue
            start0, end = int(f[3]) - 1, int(f[4])          # GTF is 1-based inclusive
            tags = f[8]
            rank = 0 if 'tag "MANE_Select"' in tags else (1 if 'tag "Ensembl_canonical"' in tags else 2)
            key = (rank, -(end - start0), a.get("transcript_id", ""))
            if g not in best or key < best[g][0]:
                best[g] = (key, (f[0], start0, end, f[6], a.get("transcript_id", "")))
    return {g: rec(v[1][0], v[1][1], v[1][2], v[1][3], "gencode_v44_gtf", v[1][4])
            for g, v in best.items()}


def parse_ucsc(wanted):
    """UCSC refGene records through shared/common.py. txStart is 0-based, txEnd exclusive.

    common.gene_coords() reports `tss` as txEnd for minus-strand genes, i.e. one past the last
    base; rec() recomputes it as body_end - 1 so the convention matches the BEDs.
    """
    out = {}
    if not wanted:
        return out
    for g, c in cm.gene_coords(sorted(wanted)).items():
        if not isinstance(c, dict):
            continue
        if not all(c.get(k) is not None for k in ("chrom", "strand", "tx_start", "tx_end")):
            continue
        out[g] = rec(c["chrom"], c["tx_start"], c["tx_end"], c["strand"],
                     "ucsc_refgene", c.get("name", ""))
    return out


def parse_ncbi_supplement(wanted):
    """Checked-in NCBI Gene E-utilities table for the 5 LOC/pseudogene cases.

    Its `tss` column is body_end for minus-strand genes (the exclusive end); rec()
    renormalises to body_end - 1. That 1 bp shift is immaterial for a +/-3 kb window but is
    fixed anyway so one convention holds across all four sources.
    """
    if not NCBI_SUPP.exists():
        return {}, pd.DataFrame()
    d = pd.read_csv(NCBI_SUPP, sep="\t")
    d = d[d.gene.isin(wanted)]
    out = {r.gene: rec(r.chrom, r.body_start, r.body_end, r.strand,
                       "ncbi_eutils", str(r.ncbi_gene_id)) for r in d.itertuples()}
    return out, d


def main():
    u = pd.read_csv(DATA / "universe.tsv", sep="\t")
    genes = u.gene.astype(str).tolist()
    wanted = set(genes)
    print(f"universe in: {len(genes)}")

    contigs = usable_contigs()
    print(f"  usable contigs from the bigWigs: "
          f"{len(contigs) if contigs is not None else 'n/a (primary-assembly rule)'}")
    resolved, prov, rejected = {}, [], []

    def ok_contig(chrom):
        if contigs is not None:
            return chrom in contigs
        return "_" not in chrom          # primary-assembly fallback rule

    def take(name, table, extra=None):
        n0 = len(resolved)
        for g in genes:
            if g not in resolved and g in table:
                if not ok_contig(table[g]["chrom"]):
                    rejected.append({"gene": g, "source": name,
                                     "chrom": table[g]["chrom"],
                                     "reason": "contig absent from the bigWigs; "
                                               "unmeasurable, so the gene falls through "
                                               "to the next source"})
                    continue
                resolved[g] = table[g]
        got = len(resolved) - n0
        print(f"  {name:<16} resolved {got:>5}   cumulative {len(resolved)}/{len(genes)}")
        prov.append({"source": name, "newly_resolved": got, "cumulative": len(resolved),
                     **(extra or {})})

    print("\nresolution, in priority order:")
    v49, v49_dups = parse_v49_body(wanted)
    take("gencode_v49_pc", v49, {"duplicate_symbols_deduped": v49_dups})

    still = wanted - set(resolved)
    print(f"  ({len(still)} unresolved, going to v44)")
    take("gencode_v44_gtf", parse_v44_gtf(still))

    take("ucsc_refgene", parse_ucsc(wanted - set(resolved)))

    ncbi, ncbi_df = parse_ncbi_supplement(wanted - set(resolved))
    take("ncbi_eutils", ncbi)

    failed = [g for g in genes if g not in resolved]

    co = pd.DataFrame.from_dict(resolved, orient="index")
    co.index.name = "gene"
    co = co.loc[[g for g in genes if g in resolved]]
    co["short_model"] = co.model_length < SHORT_MODEL_BP

    # parent_gene rides along from the NCBI supplement (multi-mapping caveat)
    co["parent_gene"] = ""
    if len(ncbi_df):
        supp = pd.read_csv(NCBI_SUPP, sep="\t").set_index("gene")
        for g in supp.index:
            if g in co.index:
                co.loc[g, "parent_gene"] = "" if pd.isna(supp.loc[g, "parent_gene"]) else supp.loc[g, "parent_gene"]

    co = co.join(u.set_index("gene")[["memory_class"]])

    # Cross-source guard: where does the accepted model disagree with v44? A symbol can be
    # attached to a readthrough transcript in one annotation and to the gene proper in another
    # (SCHIP1 / IQCJ-SCHIP1 is the worked case). Silent acceptance puts the chromatin window on
    # the wrong promoter, so every large disagreement is written out.
    v44_all = parse_v44_gtf(wanted)
    rows = []
    for g, r in co.iterrows():
        v = v44_all.get(g)
        if v is None or v["chrom"] != r.chrom:
            continue
        d_tss = int(r.tss) - int(v["tss"])
        ratio = r.model_length / max(1, v["model_length"])
        if abs(d_tss) > V44_TSS_DISAGREE_BP or ratio > 5 or ratio < 0.2:
            rows.append({"gene": g, "memory_class": r.memory_class,
                         "accepted_source": r.coordinate_source, "accepted_model": r.model_id,
                         "accepted_tss": int(r.tss), "accepted_length": int(r.model_length),
                         "v44_tss": int(v["tss"]), "v44_length": int(v["model_length"]),
                         "v44_model": v["model_id"], "tss_delta_vs_v44": d_tss,
                         "length_ratio_vs_v44": round(ratio, 3)})
    outliers = pd.DataFrame(rows).sort_values(
        "tss_delta_vs_v44", key=lambda s: s.abs(), ascending=False) if rows else pd.DataFrame()
    outliers.to_csv(DATA / "coordinates_outliers.tsv", sep="\t", index=False)

    cols = ["chrom", "body_start", "body_end", "strand", "tss", "model_length",
            "coordinate_source", "model_id", "short_model", "parent_gene", "memory_class"]
    co[cols].to_csv(DATA / "coordinates.tsv", sep="\t")

    fdf = u.set_index("gene").loc[failed] if failed else u.set_index("gene").iloc[0:0]
    fdf = fdf.assign(reason="no model in v49 / v44 / UCSC / NCBI supplement")
    fdf.to_csv(DATA / "coordinates_failed.tsv", sep="\t")

    print(f"\nresolved {len(resolved)} / {len(genes)}")
    print("\nby source:")
    print(co.coordinate_source.value_counts().to_string())
    print("\nby source x memory class:")
    print(pd.crosstab(co.coordinate_source, co.memory_class).to_string())

    if rejected:
        print(f"\nsource hits REJECTED on an unmeasurable contig ({len(rejected)}); the gene "
              "then fell through to the next source:")
        print(pd.DataFrame(rejected).to_string(index=False))

    if failed:
        print(f"\nFAILED ({len(failed)}): {', '.join(failed)}")
    else:
        print("\nFAILED: none -- every universe gene has a coordinate.")

    if len(outliers):
        print(f"\ncoordinate OUTLIERS vs the v44 canonical model (|TSS delta| > "
              f"{V44_TSS_DISAGREE_BP} bp, or length ratio outside 0.2-5): {len(outliers)}")
        print(outliers.to_string(index=False))

    short = co[co.short_model]
    print(f"\nmodels < {SHORT_MODEL_BP} bp (FLAGGED, not dropped): {len(short)}")
    sc = co.groupby("memory_class").short_model.agg(["sum", "count"])
    sc["pct"] = (100 * sc["sum"] / sc["count"]).round(1)
    print("  by class -- the K36me3 gene-body metric is degenerate for these, and the flag is "
          "strongly class-asymmetric:")
    print(sc.to_string())
    if len(short):
        print(short[["chrom", "model_length", "coordinate_source", "memory_class"]]
              .sort_values("model_length").to_string())

    print("\nanchors:")
    for a in ["IDO1", "CXCL11", "IRF1", "STAT1"]:
        if a in co.index:
            r = co.loc[a]
            print(f"  {a:<7} {r.chrom}:{r.body_start}-{r.body_end} ({r.strand}) "
                  f"tss={r.tss} len={r.model_length} src={r.coordinate_source} {r.model_id}")

    summary = {
        "source_mode": SOURCE_MODE,
        "n_universe": len(genes),
        "n_resolved": len(resolved),
        "n_failed": len(failed),
        "failed_genes": failed,
        "by_source": co.coordinate_source.value_counts().to_dict(),
        "v49_duplicate_symbols_deduped": v49_dups,
        "n_short_models": int(co.short_model.sum()),
        "short_models": short.model_length.to_dict(),
        "short_models_by_class": {c: {"n_short": int(co[co.memory_class == c].short_model.sum()),
                                      "n_total": int((co.memory_class == c).sum())}
                                  for c in ["memory", "refractory", "non_memory"]},
        "v49_aliases_applied": {k: v for k, v in V49_ALIASES.items() if k in set(co.index)},
        "rejected_unmeasurable_contig": rejected,
        "coordinate_outliers_vs_v44": (outliers.to_dict("records") if len(outliers) else []),
        "provenance": prov,
        "tss_convention": "0-based; body_start if '+' else body_end-1",
        "v49_tss_bed6_rejected": ("the v49 protein-coding TSS bed6 is not strand-aware "
                                  "(TSS==body start on both strands); TSS is derived from the "
                                  "body bed6 instead"),
    }
    (DATA / "coordinates_summary.json").write_text(json.dumps(summary, indent=2),
                                                   encoding="utf-8")
    print(f"\nwrote {DATA/'coordinates.tsv'}")


if __name__ == "__main__":
    main()
