# DESeq2 on RSEM gene-level quantification, Stimulated (Naive+IFNg) vs Rechallenged
# (Trained+IFNg). This is the contrast Fig 2K, the hysteresis count and the Fig S4G
# adjusted P values are taken from.

library(tximport)
library(DESeq2)
library(yaml)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
repo_root <- normalizePath(file.path(script_dir, ".."))
cfg <- yaml::read_yaml(file.path(repo_root, "config.yaml"))
resolve <- function(p) if (grepl("^([A-Za-z]:|/|\\\\\\\\)", p)) p else file.path(cfg$data_root, p)
out_dir <- file.path(repo_root, cfg$output_dir, "rnaseq_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c("S4-50-ngmL-IFNg-MSC", "S5-50-ngmL-IFNg-MSC", "S6-50-ngmL-IFNg-MSC",
             "S10-50-ngmL-IFNg-Trained", "S11-50-ngmL-IFNg-Trained", "S12-50-ngmL-IFNg-Trained")
files <- file.path(resolve(cfg$rnaseq$rsem_dir), paste0(samples, ".genes.results"))
stopifnot(all(file.exists(files)))

txi <- tximport(files, type = "rsem", txOut = FALSE)

conditions <- c("control", "control", "control", "treated", "treated", "treated")
colData <- data.frame(sampleName = samples, condition = factor(conditions))
rownames(colData) <- colData$sampleName

# Genes with zero or negative average transcript length break the offset DESeq2 builds
# from txi$length, so they go before any filtering on counts.
zero_length_genes <- which(rowSums(txi$length <= 0) > 0)

txi_filtered <- txi
txi_filtered$counts <- txi_filtered$counts[-zero_length_genes, ]
if (!is.null(txi_filtered$length)) {
  txi_filtered$length <- txi_filtered$length[-zero_length_genes, ]
}

rowSumsCounts <- rowSums(txi_filtered$counts)
filterThreshold <- 10
txi_filtered$counts <- txi_filtered$counts[rowSumsCounts >= filterThreshold, ]
if (!is.null(txi_filtered$length)) {
  txi_filtered$length <- txi_filtered$length[rowSumsCounts >= filterThreshold, ]
}

dds <- DESeqDataSetFromTximport(txi_filtered, colData = colData, design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)

write.csv(as.data.frame(res),
          file.path(out_dir, "DESeq2_results_Naive+IFNg_VS_Trained+IFNg.csv"))

selected_genes_down <- res[which(res$log2FoldChange < -0.5 & res$padj < 0.01), ]
selected_genes_up   <- res[which(res$log2FoldChange >  0.5 & res$padj < 0.01), ]
cat(sprintf("up: %d   down: %d\n", nrow(selected_genes_up), nrow(selected_genes_down)))

write.csv(selected_genes_down,
          file.path(out_dir, "DESeq2_results_Naive+IFNg_VS_Trained+IFNg_DOWN.csv"),
          row.names = TRUE)
write.csv(selected_genes_up,
          file.path(out_dir, "DESeq2_results_Naive+IFNg_VS_Trained+IFNg_UP.csv"),
          row.names = TRUE)
