#!/usr/bin/env python3
"""
run_tracks.py  --  multi-mark CUT&RUN track panels: Fig 3A (CXCL11) and Fig S6A (IDO1).

Reads config/config.json + config/gene_models.json and renders each locus with the
region-stacked layout: H3K4me3 and H3K36me2 over TSS +/- 3 kb, H3K27me3 and H3K36me3
over the gene body.

Run:  py run_tracks.py            (both loci)
      py run_tracks.py IDO1       (single locus)
"""
import json
import sys
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8")

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1] / "shared"))
import config as cfg  # noqa: E402
import common as cm  # noqa: E402
import tracklib as tl  # noqa: E402

CONFIG = json.loads((ROOT / "config" / "config.json").read_text(encoding="utf-8"))
MODELS = json.loads((ROOT / "config" / "gene_models.json").read_text(encoding="utf-8"))

MARK_BY_KEY = {m["key"]: m for m in CONFIG["marks"]}


def region_window(gm, region):
    """Per-mark window: 'tss3kb' = TSS +/- 3 kb; 'genebody' = the full locus window,
    unless a per-gene override widens it. CXCL11 is a small gene, so its gene-body view
    is widened to a gene-centered window for more IDO1-like context."""
    rw = CONFIG.get("region_windows", {})
    tss = gm["tss"]
    if region == "tss3kb":
        hw = rw.get("tss3kb_halfwidth", 3000)
        return {"chrom": gm["chrom"], "start": int(tss - hw), "end": int(tss + hw)}
    ov = rw.get("genebody_override", {}).get(gm.get("name_label", ""))
    if ov and ov.get("mode") == "gene_center":
        center = (gm["tx_start"] + gm["tx_end"]) // 2
        half = int(ov["span"]) // 2
        return {"chrom": gm["chrom"], "start": int(center - half), "end": int(center + half)}
    return dict(gm["window"])


def build_region_blocks(gm, variant, nbins):
    """Read each mark in the variant over its assigned window."""
    blocks = []
    for mkey, region in variant["regions"].items():
        mark = MARK_BY_KEY[mkey]
        win = region_window(gm, region)
        sig = {}
        for cond in mark["conditions"]:
            fname = mark["file_tmpl"].format(cond=cond)
            sig[cond] = tl.read_signal(cm.bigwig(fname), win["chrom"], win["start"],
                                       win["end"], nbins)
        blocks.append({"mark": mark, "region": region, "window": win, "signals": sig})
    return blocks


def main():
    only = sys.argv[1] if len(sys.argv) > 1 else None
    tl.setup()
    conditions = CONFIG["conditions"]

    loci = [l for l in CONFIG["loci"] if (only is None or l["gene"] == only)]
    for locus in loci:
        gene = locus["gene"]
        panel = locus["panel"]
        print(f"\n=== {gene} ({panel}) ===")
        gm = dict(MODELS["genes"][gene])
        gm["name_label"] = gene

        for variant in CONFIG["variants"]:
            for toggle in variant["toggles"]:
                subtitle = f"{gene}   hg38   per-mark region focus   y-axis: CPM"
                blocks = build_region_blocks(gm, variant, CONFIG["nbins"])
                fig = tl.render_region_stacked(gm, blocks, conditions, variant, toggle,
                                               subtitle=subtitle)
                tl.save(fig, cfg.out_dir(panel) / panel)
                print(f"    rendered {panel}")


if __name__ == "__main__":
    main()
