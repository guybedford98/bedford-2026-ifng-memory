"""Shared loaders and helpers used by more than one analysis folder: bigwig reading,
gene-set loading, strand-aware TSS from the protein-coding bed6, and coordinate
resolution against UCSC refGene.

    import common as cm
    prof = cm.read_signal(cm.bigwig("H3K27me3_Naive.bw"), "chr4", 76033070, 76039070, 1000)
    coords = cm.gene_coords(["CXCL11", "IDO1"])

All paths come from config.yaml through config.py.
"""
import csv
import json
import sys
import time
import urllib.request

import numpy as np

import config as cfg

_COORD_CACHE = cfg.REPO_ROOT / "cache" / "gene_coords.json"


def bigwig(filename):
    """Path to a CUT&RUN bigwig inside the configured bigwig directory."""
    return cfg.path("cutandrun", "bigwig_dir") / filename


def bam(filename):
    """Path to a CUT&RUN bam inside the configured bam directory."""
    return cfg.path("cutandrun", "bam_dir") / filename


def read_signal(bw_path, chrom, start, end, nbins):
    """Binned mean CPM over a window. Missing -> 0. Uses pybigtools."""
    import pybigtools
    bw = pybigtools.open(str(bw_path))
    try:
        vals = bw.values(chrom, int(start), int(end), bins=int(nbins),
                         summary="mean", exact=True, missing=0.0)
    finally:
        bw.close()
    return np.asarray(vals, dtype=float)


def load_gene_set(path, column="gene"):
    """Return the list of gene symbols in a one-gene-per-row CSV."""
    out = []
    with open(path, newline="", encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            sym = (row.get(column) or "").strip()
            if sym:
                out.append(sym)
    return out


def load_bed6_genes():
    """Protein-coding consensus bed6 -> {symbol: {chrom, start, end, strand, tss}}.

    The TSS is derived here by strand from the gene body; the companion TSS bed6 that
    ships with the same annotation is not strand-aware and must not be used."""
    genes = {}
    with open(cfg.path("reference", "protein_coding_bed6"), encoding="utf-8") as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            chrom, start, end, name, _score, strand = line.rstrip("\n").split("\t")[:6]
            start, end = int(start), int(end)
            genes[name] = {"chrom": chrom, "start": start, "end": end, "strand": strand,
                           "tss": start if strand == "+" else end - 1}
    return genes


def _ucsc(url):
    with urllib.request.urlopen(url, timeout=25) as r:
        return json.loads(r.read())


def _parse_field(field):
    if isinstance(field, list):
        return [int(x) for x in field if str(x).strip() != ""]
    return [int(x) for x in str(field).split(",") if x.strip() != ""]


def _resolve_one(symbol):
    """Resolve one gene symbol to hg38 coordinates via UCSC search + refGene. None if not found."""
    try:
        d = _ucsc(f"https://api.genome.ucsc.edu/search?search={symbol}&genome=hg38")
        pos = None
        for grp in d.get("positionMatches", []):
            for h in grp.get("matches", []):
                if h.get("position"):
                    pos = h["position"]
                    break
            if pos:
                break
        if not pos:
            return None
        chrom, coords = pos.split(":")
        lo, hi = (int(x) for x in coords.replace(",", "").split("-"))
        recs = _ucsc(f"https://api.genome.ucsc.edu/getData/track?genome=hg38;track=refGene;"
                     f"chrom={chrom};start={max(0, lo - 2000)};end={hi + 2000}").get("refGene", [])
        hits = [g for g in recs if str(g.get("name2", "")).upper() == symbol.upper()]
        if not hits:
            return None
        g = max(hits, key=lambda r: int(r["txEnd"]) - int(r["txStart"]))
        strand = g.get("strand", "+")
        tx0, tx1 = int(g["txStart"]), int(g["txEnd"])
        return {"chrom": chrom, "strand": strand, "tx_start": tx0, "tx_end": tx1,
                "tss": tx0 if strand == "+" else tx1,
                "cds_start": int(g.get("cdsStart", tx0)), "cds_end": int(g.get("cdsEnd", tx1)),
                "exon_starts": _parse_field(g.get("exonStarts", [])),
                "exon_ends": _parse_field(g.get("exonEnds", [])), "name": g.get("name", "")}
    except Exception as exc:
        print(f"    coord resolve failed for {symbol}: {exc!r}", file=sys.stderr)
        return None


def gene_coords(symbols, sleep=0.05):
    """Resolve many symbols to coordinates, cached to cache/gene_coords.json.
    Returns {symbol: coord_dict}; unresolved symbols map to None."""
    cache = {}
    if _COORD_CACHE.exists():
        cache = json.loads(_COORD_CACHE.read_text(encoding="utf-8"))
    changed = False
    for sym in symbols:
        if sym in cache:
            continue
        cache[sym] = _resolve_one(sym)
        changed = True
        if sleep:
            time.sleep(sleep)
    if changed:
        _COORD_CACHE.parent.mkdir(parents=True, exist_ok=True)
        _COORD_CACHE.write_text(json.dumps(cache, indent=2), encoding="utf-8")
    return {s: cache.get(s) for s in symbols}


def tss_window(coord, halfwidth):
    """(chrom, start, end) for TSS +/- halfwidth."""
    return coord["chrom"], coord["tss"] - halfwidth, coord["tss"] + halfwidth


def gene_body_window(coord, pad=0):
    """(chrom, start, end) for the gene body, optionally padded."""
    return coord["chrom"], coord["tx_start"] - pad, coord["tx_end"] + pad


def window_mean(bw_path, chrom, start, end, nbins):
    """Mean signal over one window."""
    return float(np.nanmean(read_signal(bw_path, chrom, start, end, nbins)))


def loss_fraction(naive, trained):
    """Fractional naive-vs-trained loss = (naive - trained) / naive. Positive = lost."""
    if naive is None or naive <= 0:
        return float("nan")
    return (naive - trained) / naive
