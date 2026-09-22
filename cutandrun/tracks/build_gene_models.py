#!/usr/bin/env python3
"""
build_gene_models.py

Fetch hg38 gene models (exon structure, strand, TSS) for the track loci from the
UCSC refGene track, pick one representative transcript per gene, compute a tight
promoter-plus-body plot window, and bake everything into config/gene_models.json.

The models are baked to a static JSON so track plotting is deterministic and needs
no network at render time. If UCSC is unreachable, a hardcoded fallback (verified
against UCSC on build day) is used instead.

Run:  py build_gene_models.py
"""
import json
import os
import sys
import urllib.request

sys.stdout.reconfigure(encoding="utf-8")

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "config", "gene_models.json")
GENOME = "hg38"

# Genes to resolve. The approx window narrows the refGene query so we never pull a
# neighbouring gene: a stale IDO1 coordinate in an earlier locus file pointed at ADAM2,
# about 170 kb away.
GENES = {
    "IDO1":   {"chrom": "chr8", "approx": (39_910_000, 39_932_000)},
    "CXCL11": {"chrom": "chr4", "approx": (76_030_000, 76_040_000)},
    "IRF1":   {"chrom": "chr5", "approx": (132_480_000, 132_495_000)},
    "STAT1":  {"chrom": "chr2", "approx": (190_960_000, 191_020_000)},
}

# Promoter-side and gene-body-side padding (bp) for the plot window.
PROMOTER_PAD = 3000
BODY_PAD = 2500

# Fallback models (verified against UCSC hg38) in case the API is down.
FALLBACK = {
    "IDO1": {
        "chrom": "chr8", "strand": "+", "tx_start": 39_913_891, "tx_end": 39_928_790,
        "name": "NM_002164", "source": "fallback",
    },
    "CXCL11": {
        "chrom": "chr4", "strand": "-", "tx_start": 76_033_682, "tx_end": 76_036_070,
        "name": "NM_005409", "source": "fallback",
    },
    "IRF1": {
        "chrom": "chr5", "strand": "-", "tx_start": 132_481_447, "tx_end": 132_490_976,
        "name": "NM_002198", "source": "fallback",
    },
    "STAT1": {
        "chrom": "chr2", "strand": "-", "tx_start": 190_965_024, "tx_end": 191_014_313,
        "name": "NM_007315", "source": "fallback",
    },
}


def ucsc_json(url):
    with urllib.request.urlopen(url, timeout=20) as r:
        return json.loads(r.read())


def parse_field(field):
    if isinstance(field, list):
        return [int(x) for x in field if str(x).strip() != ""]
    return [int(x) for x in str(field).split(",") if x.strip() != ""]


def fetch_refgene(gene, chrom, lo, hi):
    url = (f"https://api.genome.ucsc.edu/getData/track?genome={GENOME}"
           f";track=refGene;chrom={chrom};start={lo};end={hi}")
    recs = ucsc_json(url).get("refGene", [])
    hits = [g for g in recs if str(g.get("name2", "")).upper() == gene.upper()]
    if not hits:
        return None
    # Representative transcript: the one with the widest span (most complete body).
    hits.sort(key=lambda g: int(g["txEnd"]) - int(g["txStart"]), reverse=True)
    return hits[0]


def build_one(gene, cfg):
    chrom = cfg["chrom"]
    lo, hi = cfg["approx"]
    model = None
    try:
        g = fetch_refgene(gene, chrom, lo, hi)
        if g is not None:
            es = parse_field(g.get("exonStarts", []))
            ee = parse_field(g.get("exonEnds", []))
            model = {
                "chrom": chrom,
                "strand": g.get("strand", "+"),
                "tx_start": int(g["txStart"]),
                "tx_end": int(g["txEnd"]),
                "cds_start": int(g.get("cdsStart", g["txStart"])),
                "cds_end": int(g.get("cdsEnd", g["txEnd"])),
                "exon_starts": es,
                "exon_ends": ee,
                "name": g.get("name", ""),
                "source": "ucsc_refGene",
            }
            print(f"[UCSC] {gene}: {chrom}:{model['tx_start']:,}-{model['tx_end']:,} "
                  f"({model['strand']}) {len(es)} exons  {model['name']}")
    except Exception as exc:
        print(f"[UCSC failed for {gene}: {exc!r}] using fallback", file=sys.stderr)

    if model is None:
        fb = dict(FALLBACK[gene])
        fb.setdefault("cds_start", fb["tx_start"])
        fb.setdefault("cds_end", fb["tx_end"])
        fb.setdefault("exon_starts", [fb["tx_start"]])
        fb.setdefault("exon_ends", [fb["tx_end"]])
        model = fb
        print(f"[fallback] {gene}: {model['chrom']}:{model['tx_start']:,}-{model['tx_end']:,} "
              f"({model['strand']})")

    # TSS is txStart on + strand, txEnd on - strand.
    tss = model["tx_start"] if model["strand"] == "+" else model["tx_end"]
    model["tss"] = tss

    # Tight promoter-plus-body window: pad the promoter side beyond the TSS and the
    # far (gene-end) side beyond the TES, in genomic coordinates.
    if model["strand"] == "+":
        win_start = tss - PROMOTER_PAD
        win_end = model["tx_end"] + BODY_PAD
    else:
        win_start = model["tx_start"] - BODY_PAD
        win_end = tss + PROMOTER_PAD
    model["window"] = {"chrom": model["chrom"], "start": int(win_start), "end": int(win_end)}
    span = win_end - win_start
    print(f"        window {model['chrom']}:{win_start:,}-{win_end:,}  ({span/1000:.1f} kb)")
    return model


def main():
    out = {"genome": GENOME, "promoter_pad": PROMOTER_PAD, "body_pad": BODY_PAD, "genes": {}}
    for gene, cfg in GENES.items():
        out["genes"][gene] = build_one(gene, cfg)
    with open(OUT, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
