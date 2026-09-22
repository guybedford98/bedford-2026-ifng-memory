# rnaseq

Bulk RNA-seq of the IFNg and TNFa training sets: the ISG universe and its memory classification,
the DESeq2 contrasts, and the Fig 1, Fig S2, Fig S3, Fig S4, Fig S5 and Fig S6 RNA-seq
panels. `py run_all.py` runs the numbered stages in order (mixed python and `Rscript`) and then
the panel scripts; stage 00 builds the universe every other number in the paper hangs off, and
stages 03-06 write the DESeq2 tables. Panel scripts read the deposited DESeq2 tables through
`config.yaml`, so they can be run without re-running the R stages.
