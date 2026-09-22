# CUT&RUN promoter signal

`quantify_windows.py` writes the per-replicate and per-condition CPM tables in `data/` and renders
Fig 3B (CXCL11) and Fig S6B (IDO1), taking its windows from `../tracks/config` so the bars match
the rendered tracks; pass `--force` to re-read the bigwigs instead of the cached tables.
`make_FigS6C.py` draws the naive promoter H3K27me3 distribution across the mappable ISG promoters
in `data/isg_promoter_focal_cpm.csv`.
