# Naive vs Rested DEGs with all four IFNg conditions fitted in one ~condition model.
# Input: RSEM expected counts for the 12 IFNg-training libraries.
# This is the run behind the DEG count printed on the Fig S2D volcano; stage 04 fits
# the two conditions on their own and gives the count quoted in the Results text.

library(DESeq2)
library(readxl)
library(dplyr)
library(yaml)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
repo_root <- normalizePath(file.path(script_dir, ".."))
cfg <- yaml::read_yaml(file.path(repo_root, "config.yaml"))
resolve <- function(p) if (grepl("^([A-Za-z]:|/|\\\\\\\\)", p)) p else file.path(cfg$data_root, p)
out_dir <- file.path(repo_root, cfg$output_dir, "rnaseq_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

raw <- read_excel(resolve(cfg$rnaseq$counts))

# Round RSEM expected-counts to integers and use gene_id as rownames
counts_mat <- raw %>%
  select(-`transcript_id(s)`) %>%
  mutate(across(-gene_id, ~ round(.x))) %>%
  tibble::column_to_rownames("gene_id") %>%
  as.matrix()

coldata <- data.frame(
  sample    = colnames(counts_mat),
  condition = factor(
    c("Naive","Naive","Naive",
      "Stimulated","Stimulated","Stimulated",
      "Rested","Rested","Rested",
      "Rechallenged","Rechallenged","Rechallenged"),
    levels = c("Naive","Stimulated","Rested","Rechallenged")
  ),
  row.names = colnames(counts_mat)
)

keep        <- rowSums(counts_mat) >= 10
counts_filt <- counts_mat[keep, ]
cat(sprintf("Genes before filter: %d  |  after filter: %d\n",
            nrow(counts_mat), nrow(counts_filt)))

dds_all <- DESeqDataSetFromMatrix(
  countData = counts_filt,
  colData   = coldata,
  design    = ~ condition
)
dds_all <- DESeq(dds_all)

res <- results(dds_all,
               contrast  = c("condition", "Rested", "Naive"),
               alpha     = 0.05)
res_df <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  arrange(padj)

n_sig  <- sum(res_df$padj < 0.05, na.rm = TRUE)
n_up   <- sum(res_df$padj < 0.05 & res_df$log2FoldChange > 0, na.rm = TRUE)
n_down <- sum(res_df$padj < 0.05 & res_df$log2FoldChange < 0, na.rm = TRUE)
cat(sprintf("\nNaive vs. Rested DEGs (padj < 0.05): %d total  |  %d up  |  %d down\n",
            n_sig, n_up, n_down))

write.csv(res_df, file.path(out_dir, "DEG_Naive_vs_Rested_full_model.csv"),
          row.names = FALSE)
cat("DEG table saved.\n")
