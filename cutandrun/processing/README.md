# CUT&RUN processing

Shell pipeline that turns the CUT&RUN FASTQs into blacklist-filtered hg38 BAMs and CPM bigWigs.
Reference and FASTQ paths come from `config.yaml` through `config.sh`, and the tools (bowtie2,
samtools, fastp, fastqc, multiqc, deeptools, bedtools, macs2, trim_galore, cutadapt, SEACR) must
be on PATH.

```
bash run_all.sh [output_root]
```

Run order:

1. `00_crop_fastqs.sh` Trim Galore, then cutadapt crop of both reads to 50 bp.
2. `01_merge_resequenced_fastqs.sh` concatenates the resequenced H3K27me3 FASTQs with the original
   reads of the same library.
3. `02_cutandrun_pipeline.sh <samples.tsv> <out_dir>` once per sample sheet in `data/`.
4. `03_merge_replicates.sh` merges replicate BAMs into CPM bigWigs for H3K27ac, H3K4me3, H3K36me2,
   H3K36me3 and IgG.
5. `04_merge_h3k27me3_deep.sh` does the same for H3K27me3, CPM and spike-in.

`data/samples_GB249.tsv` (H3K27me3, H3K4me3, H3K27ac across the four conditions),
`data/samples_GB263.tsv` (H3K36me2, H3K36me3, naive and stimulated) and `data/samples_GB301.tsv`
(deep-sequenced H3K27me3) map each library to mark, condition, replicate and its FASTQ pair,
relative to `cutandrun.fastq_dir`. `data/spike_in_stats_*.tsv` are the K-MetStat read counts and
scale factors from those runs.
