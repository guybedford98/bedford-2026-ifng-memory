# CUT&RUN genome-browser tracks

`run_tracks.py` renders Fig 3A (CXCL11) and Fig S6A (IDO1); `make_Fig3C.py` renders the
1 Mb H3K27me3 focal-versus-distal view at CXCL11. Both read the merged CPM bigwigs named in
`config/config.json` from `cutandrun.bigwig_dir` and the hg38 locus definitions from
`config/gene_models.json`, which `build_gene_models.py` regenerates from the UCSC refGene API
if the loci change.
