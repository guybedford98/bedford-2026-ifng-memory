# CRISPRa CXCL11 RNA-seq

pydeseq2 re-analysis of the day-10 and day-21 gene-level count matrices, each timepoint fitted
independently on a single four-level group factor (NT_rest, NT_IFNg, sg1_rest, sg1_IFNg). Run
`python run_all.py`, which runs `00_pydeseq2.py` and `01_compute_tpm.py` and then draws the
panels from their outputs. Alignment and counting were done by the sequencing provider, so only
the gene-level matrices are read here.
