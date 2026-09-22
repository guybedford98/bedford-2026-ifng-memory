# Rscript 00_run_gsea.R <DESeq2 results csv> <output GSEA table csv>

args <- commandArgs(trailingOnly = TRUE)
deseq2_csv <- args[1]
out_csv    <- args[2]

# ------------------------------
#   LOAD PACKAGES
# ------------------------------
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)


# ------------------------------
#   LOAD DESEQ2 RESULTS CSV
# ------------------------------
res <- read.csv(deseq2_csv,
                header = TRUE,
                row.names = 1,
                check.names = FALSE)

# Remove rows with missing stat values
res <- res[!is.na(res$stat), ]


# ------------------------------
#   CREATE RANKED GENE LIST
# ------------------------------
# The DESeq2 design makes stat > 0 = UP IN TREATED
geneList <- res$stat
names(geneList) <- rownames(res)

# Sort decreasing as required by GSEA
geneList <- sort(geneList, decreasing = TRUE)


# ------------------------------
#   RUN GSEA (GO Biological Process)
# ------------------------------
ego_bp <- gseGO(
  geneList     = geneList,
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  ont          = "BP",
  minGSSize    = 10,
  maxGSSize    = 500,
  pvalueCutoff = 0.5,
  verbose      = TRUE
)

# Save table
write.csv(as.data.frame(ego_bp),
          out_csv,
          row.names = FALSE)
