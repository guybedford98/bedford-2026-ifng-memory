# gsea

clusterProfiler `gseGO` on GO Biological Process for the IFNg RNA-seq, ranked on the DESeq2 Wald
statistic. `python run_all.py` runs the enrichment for both contrasts (`Rscript` must be on PATH),
writes the tables to `results/gsea_tables/`, then draws Fig S2A and Fig S2E. The five terms in
`data/GSEA_top5_IFNg.csv` were picked by hand as a non-redundant set, so that file is checked in
rather than generated, and its two `pvalue` columns hold clusterProfiler's `p.adjust`.
