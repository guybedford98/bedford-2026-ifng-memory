# Fig S3C: PCA of the four TNFa training conditions, VST counts, PC2 on x and PC1 on y.
# Mirrors make_Fig1G.R on the TNFa cross-training RNA-seq (RSEM/hg38), 2 replicates
# per condition:
#   Naive        = untreated NTC                  (S1, S2)
#   Stimulated   = first-time 50 ng/mL TNFa       (S3, S4)
#   Rested       = TNFa-trained, then rested      (S7, S8)
#   Rechallenged = TNFa-trained, rested, re-stim  (S9, S10)

library(DESeq2)
library(ggplot2)
library(yaml)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
repo_root <- normalizePath(file.path(script_dir, ".."))
cfg <- yaml::read_yaml(file.path(repo_root, "config.yaml"))
resolve <- function(p) if (grepl("^([A-Za-z]:|/|\\\\\\\\)", p)) p else file.path(cfg$data_root, p)
fig_dir <- file.path(repo_root, cfg$output_dir, "FigS3C")
tab_dir <- file.path(repo_root, cfg$output_dir, "rnaseq_tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

samples <- c(
  "S1-MSC-NTC", "S2-MSC-NTC",
  "S3-MSC-50-ngml-TNFa", "S4-MSC-50-ngml-TNFa",
  "S7-TNFa-Trained", "S8-TNFa-Trained",
  "S9-TNFa-Trained-50-ngml-TNFa", "S10-TNFa-Trained-50-ngml-TNFa"
)
files <- file.path(resolve(cfg$rnaseq$tnfa_rsem_dir), paste0(samples, ".genes.results"))
stopifnot(all(file.exists(files)))

# Build the RSEM expected-count matrix, aligned by gene_id
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
    c("Naive","Naive",
      "Stimulated","Stimulated",
      "Rested","Rested",
      "Rechallenged","Rechallenged"),
    levels = c("Naive","Stimulated","Rested","Rechallenged")),
  row.names = colnames(counts_mat))

counts_filt <- counts_mat[rowSums(counts_mat) >= 10, ]
cat(sprintf("Genes before filter: %d  |  after filter: %d\n",
            nrow(counts_mat), nrow(counts_filt)))

dds <- DESeqDataSetFromMatrix(counts_filt, coldata, ~condition)
dds <- estimateSizeFactors(dds)
vsd <- vst(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
cat(sprintf("PC1 %s%%  PC2 %s%%\n", pct_var[1], pct_var[2]))

condition_colors <- c(
  "Naive"        = "#909499",
  "Stimulated"   = "#4D9AD4",
  "Rested"       = "#4A4A4B",
  "Rechallenged" = "#0B519D"
)

p <- ggplot(pca_data, aes(x = PC2, y = PC1, color = condition, shape = condition)) +
  geom_point(size = 4, stroke = 0.8) +
  scale_color_manual(values = condition_colors, name = "Condition") +
  scale_shape_manual(values = c(16, 17, 15, 18), name = "Condition") +
  coord_fixed() +
  labs(
    x = sprintf("PC2  (%s%%)", pct_var[2]),
    y = sprintf("PC1  (%s%%)", pct_var[1])
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    legend.title    = element_text(face = "bold"),
    axis.line       = element_line(linewidth = 0.5),
    panel.grid      = element_blank()
  )

ggsave(file.path(fig_dir, "FigS3C.pdf"), p, width = 8, height = 5, useDingbats = FALSE)
ggsave(file.path(fig_dir, "FigS3C.png"), p, width = 8, height = 5, dpi = 600)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(fig_dir, "FigS3C.svg"), p, width = 8, height = 5)
} else {
  message("svglite not installed - skipping SVG export (the PDF is vector-editable).")
}

write.csv(pca_data, file.path(tab_dir, "TNFa_PCA_coordinates.csv"), row.names = FALSE)
