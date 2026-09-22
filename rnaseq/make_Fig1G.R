# Fig 1G: PCA of the four IFNg training conditions, VST counts, PC2 on x and PC1 on y.

library(DESeq2)
library(ggplot2)
library(readxl)
library(yaml)

args <- commandArgs(trailingOnly = FALSE)
script_dir <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
repo_root <- normalizePath(file.path(script_dir, ".."))
cfg <- yaml::read_yaml(file.path(repo_root, "config.yaml"))
resolve <- function(p) if (grepl("^([A-Za-z]:|/|\\\\\\\\)", p)) p else file.path(cfg$data_root, p)
fig_dir <- file.path(repo_root, cfg$output_dir, "Fig1G")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

raw        <- read_excel(resolve(cfg$rnaseq$counts))
counts_mat <- raw[, c(-1, -2)]
rownames(counts_mat) <- raw$gene_id
counts_mat <- round(as.matrix(counts_mat))

coldata <- data.frame(
  condition = factor(
    c("Naive","Naive","Naive",
      "Stimulated","Stimulated","Stimulated",
      "Rested","Rested","Rested",
      "Rechallenged","Rechallenged","Rechallenged"),
    levels = c("Naive","Stimulated","Rested","Rechallenged")),
  row.names = colnames(counts_mat))

counts_filt <- counts_mat[rowSums(counts_mat) >= 10, ]
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

ggsave(file.path(fig_dir, "Fig1G.pdf"), p, width = 8, height = 5, useDingbats = FALSE)
ggsave(file.path(fig_dir, "Fig1G.png"), p, width = 8, height = 5, dpi = 600)
if (requireNamespace("svglite", quietly = TRUE)) {
  ggsave(file.path(fig_dir, "Fig1G.svg"), p, width = 8, height = 5)
} else {
  message("svglite not installed - skipping SVG export (the PDF is vector-editable).")
}
