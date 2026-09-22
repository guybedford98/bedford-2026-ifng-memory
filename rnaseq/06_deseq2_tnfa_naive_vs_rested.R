# TNFa Naive vs Rested DEGs, the mirror of stage 05 on the TNFa cross-training set.
# DESeq2 across all 4 TNFa conditions, then contrast Rested vs Naive (padj < 0.05).
# Four conditions, 2 biological replicates each:
#   Naive        = untreated NTC                  (S1, S2)
#   Stimulated   = first-time 50 ng/mL TNFa       (S3, S4)
#   Rested       = TNFa-trained, then rested      (S7, S8)
#   Rechallenged = TNFa-trained, rested, re-stim  (S9, S10)

library(DESeq2)
library(yaml)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
repo_root <- normalizePath(file.path(script_dir, ".."))
cfg <- yaml::read_yaml(file.path(repo_root, "config.yaml"))
resolve <- function(p) if (grepl("^([A-Za-z]:|/|\\\\\\\\)", p)) p else file.path(cfg$data_root, p)
out_dir <- file.path(repo_root, cfg$output_dir, "rnaseq_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "S1-MSC-NTC", "S2-MSC-NTC",
  "S3-MSC-50-ngml-TNFa", "S4-MSC-50-ngml-TNFa",
  "S7-TNFa-Trained", "S8-TNFa-Trained",
  "S9-TNFa-Trained-50-ngml-TNFa", "S10-TNFa-Trained-50-ngml-TNFa"
)
files <- file.path(resolve(cfg$rnaseq$tnfa_rsem_dir), paste0(samples, ".genes.results"))
stopifnot(all(file.exists(files)))

mats <- lapply(files, function(f) {
  d <- read.delim(f, stringsAsFactors = FALSE)
  setNames(d$expected_count, d$gene_id)
})
common     <- Reduce(intersect, lapply(mats, names))
counts_mat <- sapply(mats, function(x) x[common])
rownames(counts_mat) <- common
colnames(counts_mat) <- samples
counts_mat <- round(as.matrix(counts_mat))

coldata <- data.frame(
  condition = factor(
    c("Naive","Naive","Stimulated","Stimulated",
      "Rested","Rested","Rechallenged","Rechallenged"),
    levels = c("Naive","Stimulated","Rested","Rechallenged")),
  row.names = colnames(counts_mat))

counts_filt <- counts_mat[rowSums(counts_mat) >= 10, ]
dds <- DESeqDataSetFromMatrix(counts_filt, coldata, ~condition)
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition","Rested","Naive"), alpha = 0.05)
res_df <- as.data.frame(res)
res_df <- res_df[order(res_df$padj), ]
res_df$gene <- rownames(res_df)

n_sig  <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
n_sig_lfc <- sum(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, na.rm = TRUE)

cat(sprintf("TNFa Naive vs Rested DEGs (padj<0.05): %d total | %d up | %d down\n", n_sig, n_up, n_down))
cat(sprintf("  of those, |log2FC|>1: %d\n", n_sig_lfc))

write.csv(res_df, file.path(out_dir, "DEG_TNFa_Naive_vs_Rested.csv"), row.names = FALSE)
cat("Saved DEG_TNFa_Naive_vs_Rested.csv\n")
