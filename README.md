Analysis code for Bedford et al. Mesenchymal stromal cells were trained with IFNg, rested and
rechallenged; this repository holds the RNA-seq, CUT&RUN, DAF-seq and assay analyses that produce
the manuscript figures.

## Contents

- `rnaseq/` bulk RNA-seq of the IFNg and TNFa training sets: ISG universe, memory classification, DESeq2 contrasts, k-means heatmaps, PCA, scatter/volcano/MA and TPM panels
- `gsea/` gene-set enrichment on the IFNg RNA-seq
- `crispra_rnaseq/` CRISPRa CXCL11 RNA-seq, re-analyzed with pydeseq2
- `cutandrun/` CUT&RUN processing, browser tracks and promoter-window quantification
- `multiomic/` multi-omic feature table, PCA and Gaussian mixture clustering of the ISG universe
- `dafseq/` DAF-seq accessibility (MSP) and transcription-factor footprints
- `shared/` config loader, figure style and bar-graph statistics used by the folders above

## Setup

```
conda env create -f environment.yml
conda activate ifng-msc-memory
```

The conda environment covers the command-line tools as well. For the Python side alone, use
`pip install -r requirements.txt` on Python 3.12 or 3.13. The R side needs tximport, DESeq2,
clusterProfiler, org.Hs.eg.db, DOSE, tidyverse, readxl and yaml; versions are in
`r_session_info.txt`.

Every path lives in `config.yaml`. Edit it so each key points at your local copy of the data;
scripts read paths through `shared/config.py`, and outputs land under the `output_dir` key.
Sequencing data are deposited at GEO under [GEO accession]; processed tables and figure source
data are archived at [Zenodo DOI].

## Run order

`rnaseq/run_all.py` first: it writes the protein-coding ISG universe that `multiomic/` starts
from, and the DESeq2 tables that `gsea/run_all.py` and the CRISPRa panels read.
`cutandrun/processing/run_all.sh` has to finish before anything that reads bigwigs, which is
`cutandrun/tracks/`, `cutandrun/signal/` and `multiomic/run_all.py`. `dafseq/processing/` runs
before `dafseq/accessibility/run_all.py` and `dafseq/tf_footprints/run_all.py`.

## Figure panels

The table below maps each manuscript panel to the script that produces it.

| Panel    | What it shows                                                                | Script                                    |
| -------- | ---------------------------------------------------------------------------- | ----------------------------------------- |
| Fig 1E   | Rechallenged vs Stimulated mean TPM per gene, linear axes                    | `rnaseq/make_Fig1E.py`                  |
| Fig 1F   | k-means heatmap of the 1,063-gene protein-coding ISG universe                | `rnaseq/make_Fig1F.py`                  |
| Fig 1G   | PCA of the four IFNg training conditions, VST counts                         | `rnaseq/make_Fig1G.R`                   |
| Fig 2I   | CXCL11 TPM at day 10 and day 21, four CRISPRa conditions                     | `crispra_rnaseq/make_Fig2I.py`          |
| Fig 2K   | Rechallenge MA plots, IFNg training and CRISPRa priming at day 10 and day 21 | `crispra_rnaseq/make_Fig2K.py`          |
| Fig 3A   | Multi-mark CUT&RUN tracks at CXCL11, four conditions                         | `cutandrun/tracks/run_tracks.py`        |
| Fig 3B   | Mean CPM per mark per condition in the displayed CXCL11 windows              | `cutandrun/signal/quantify_windows.py`  |
| Fig 3C   | Wide H3K27me3 view at CXCL11, focal loss vs distal Polycomb domain           | `cutandrun/tracks/make_Fig3C.py`        |
| Fig 4A   | Single-fiber MSP rasters at IDO1, CXCL11 and STAT1                           | `dafseq/accessibility/make_Fig4A.py`    |
| Fig 4B   | Mean per-fiber MSP-covered area by condition, three loci                     | `dafseq/accessibility/make_Fig4B_4C.py` |
| Fig 4C   | Fraction of MSP-rich fibers by condition, IDO1 and CXCL11                    | `dafseq/accessibility/make_Fig4B_4C.py` |
| Fig 5A   | Aggregate Naive H3K27me3 around the TSS, memory / refractory / non-memory    | `multiomic/make_Fig5A_FigS9A.py`        |
| Fig 5B   | PC1 x PC2 scatter coloured by memory class                                   | `multiomic/make_pc_scatters.py`         |
| Fig 5C   | Cluster sizes across the four GMM clusters                                   | `multiomic/make_Fig5C_bars.py`          |
| Fig 5C   | Per-cluster memory / refractory / non-memory composition pies                | `multiomic/make_Fig5C_pies.py`          |
| Fig 5D   | PC1 x PC2 scatter coloured by GMM cluster                                    | `multiomic/make_pc_scatters.py`         |
| Fig 5E   | Per-cluster RNA response, memory vs non-memory violins                       | `multiomic/make_Fig5E.py`               |
| Fig 5F   | Chromatin metaplot grid, mark by cluster                                     | `multiomic/make_Fig5F.py`               |
| Fig S2A  | Five GO biological processes, NES in stimulated and in rested                | `gsea/make_FigS2A.R`                    |
| Fig S2D  | Rested vs Naive volcano, IFNg training set                                   | `rnaseq/make_FigS2D.py`                 |
| Fig S2E  | Diverging NES bars for GO terms, naive vs rested                             | `gsea/make_FigS2E.py`                   |
| Fig S3B  | k-means heatmap of the protein-coding TNFa universe                          | `rnaseq/make_FigS3B.py`                 |
| Fig S3C  | PCA of the four TNFa training conditions, VST counts                         | `rnaseq/make_FigS3C.R`                  |
| Fig S4G  | IDO1 and CXCL11 TPM across the four IFNg conditions                          | `rnaseq/make_FigS4G_FigS5C.py`          |
| Fig S5A  | Rested-state volcanoes, CRISPRa targeted vs non-targeted                     | `crispra_rnaseq/make_FigS5A.py`         |
| Fig S5C  | CXCR3 TPM across the four IFNg conditions                                    | `rnaseq/make_FigS4G_FigS5C.py`          |
| Fig S6A  | Multi-mark CUT&RUN tracks at IDO1, four conditions                           | `cutandrun/tracks/run_tracks.py`        |
| Fig S6B  | Mean CPM per mark per condition in the displayed IDO1 windows                | `cutandrun/signal/quantify_windows.py`  |
| Fig S6C  | Naive promoter H3K27me3 across mappable ISG promoters                        | `cutandrun/signal/make_FigS6C.py`       |
| Fig S6D  | CXCL9 and CXCL10 TPM across the four IFNg conditions                         | `rnaseq/make_FigS6D.py`                 |
| Fig S8A  | Aggregate TF-footprint density across the IDO1 amplicon                      | `dafseq/tf_footprints/make_FigS8A.py`   |
| Fig S8B  | IDO1 TSS +/-100 bp footprint density with motif annotation                   | `dafseq/tf_footprints/make_FigS8B.py`   |
| Fig S8C  | Aggregate TF-call density across the CXCL11 amplicon                         | `dafseq/tf_footprints/make_FigS8C.py`   |
| Fig S8D  | CXCL11 promoter-proximal element, TSS +225 to +425 bp                        | `dafseq/tf_footprints/make_FigS8D.py`   |
| Fig S9A  | Naive promoter H3K27me3 by memory class, log10 mean CPM                      | `multiomic/make_Fig5A_FigS9A.py`        |
| Fig S9B  | Contribution of the 11 input features to PC1, PC2 and PC3                    | `multiomic/make_FigS9B.py`              |
| Fig S9C  | PC1 x PC3 scatter coloured by memory class                                   | `multiomic/make_pc_scatters.py`         |
| Fig S10A | GMM BIC across k = 2 to 10                                                   | `multiomic/make_FigS10A.py`             |
| Fig S10B | PC1 x PC3 scatter coloured by GMM cluster                                    | `multiomic/make_pc_scatters.py`         |
| Fig S10C | Per-cluster mean log2FC heatmap over the 11 features                         | `multiomic/make_FigS10C.py`             |
| Fig S10D | Memory vs non-memory within a cluster, four chromatin marks                  | `multiomic/make_FigS10D.py`             |
