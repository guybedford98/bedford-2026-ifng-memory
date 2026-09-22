# multiomic

Builds the 11-feature RNA + chromatin table for the protein-coding ISG universe, runs PCA and
Gaussian mixture clustering on it, and draws Fig 5A-F, Fig S9A-C and Fig S10A-D. `py run_all.py`
runs the four numbered stages in order and then the panel scripts; stage 00 reads the ISG universe
CSV written by `rnaseq/` (`multiomic.isg_universe`), and stages 01 and 02 read the GENCODE
annotations and the CUT&RUN bigwigs, caching their extracted signal under `data/`.

Data files carry 0-indexed `cluster_0..3`; every figure displays 1-4 via `figstyle.clabel()`.
PC2 and PC3 are near-degenerate here (12.07% vs 12.38%) and are reported as a 2D subspace: the
RNA-priming component is labelled PC2 and the H3K27me3-leaning one PC3, the reverse of raw
variance order (`SWAP_PC2_PC3` in `03_cluster.py`). The GMM is affine-equivariant, so the
relabelling does not change the clustering.
