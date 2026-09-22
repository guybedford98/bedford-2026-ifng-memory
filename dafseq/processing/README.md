# DAF-seq processing

Alignment and footprint calling for the targeted DAF-seq amplicons; both steps need external
tools, so run them before anything in `../accessibility` or `../tf_footprints`.

## 1. Alignment, strand filtering and deamination decoration

Reads were processed with the Stergachis lab DAF-seq QC workflow (`StephanieBohaczuk/DAF-QC-SMK`,
run as `pixi run snakemake --configfile config/config.yaml --cores 8`), which is not
redistributed here. Pinned tool versions were minimap2 2.30, samtools 1.22.1, pbmarkdup 1.2.0 and
rustybam 0.1.34, and the config values used were `platform: ont`, `is_fastq: True`,
`chimera_cutoff: 0.9`, `min_deamination_count: 50`, `end_tolerance: 30`,
`decorated_samplesize: 5000`.

```
minimap2 -t {threads} --MD --secondary=no -Y -y -a -x map-ont {ref.fa} {reads.fq} | samtools sort > {sample}.mapped.reads.bam
```

Read filtering keeps primary alignments (`samtools view -F 2306`) whose QC strand call is CT or
GA; decoration rewrites deaminated positions to IUPAC Y (CT) or R (GA) and sets the `st:Z` strand
tag.

## 2. FiberHMM

FiberHMM 2.10.5, one locus series at a time. `01_run_fiberhmm.sh` calls and extracts each sample,
reading `dafseq.bam_dir`, `dafseq.fiberhmm_dir` and `reference.hg38_fasta` from `config.yaml` and
writing the recalled BAM plus the MSP, TF and footprint BED12 files into `dafseq.fiberhmm_dir`.

```
bash processing/02_run_fiberhmm_locus_series.sh IDO1
bash processing/02_run_fiberhmm_locus_series.sh CXCL11
bash processing/02_run_fiberhmm_locus_series.sh STAT1
```
