# Rscript make_FigS2A.R <top5 csv> <output pdf>

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
top5_csv <- args[1]
out_pdf  <- args[2]

# 1. Load table ------------------------------------------------------------
# The two "pvalue" columns hold clusterProfiler's p.adjust, not the raw permutation p.
df <- read.csv(top5_csv, check.names = FALSE)

# 2. Build long format -----------------------------------------------------
df_long <- bind_rows(
  df %>%
    transmute(
      ID,
      Description,
      Condition = "STIM",
      NES      = `NES STIM`,
      pvalue   = `pvalue STIM`
    ),
  df %>%
    transmute(
      ID,
      Description,
      Condition = "Rested",
      NES      = `NES Rested`,
      pvalue   = `pvalue Rested`
    )
)

# Stim first
df_long$Condition <- factor(df_long$Condition, levels = c("STIM", "Rested"))

# Order by NES in Stim
df_long <- df_long %>%
  group_by(Description) %>%
  mutate(orderNES = max(NES)) %>%
  ungroup() %>%
  arrange(desc(orderNES))

# Wrap long pathway names
df_long$Description <- stringr::str_wrap(df_long$Description, width = 30)
df_long$Description <- factor(df_long$Description, levels = unique(df_long$Description))

# Compute -log10(p)
df_long$neglogp <- -log10(df_long$pvalue)

# 3. Plot ------------------------------------------------------------------
p <- ggplot(df_long, aes(x = neglogp, y = Description, fill = NES)) +
  geom_col(width = 0.6, color = "black") +

  scale_fill_gradient2(
    name = "NES",
    low  = "darkblue",
    mid  = "white",
    high = "darkred",
    midpoint = 0,
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +

  facet_wrap(~ Condition, ncol = 2) +

  theme_bw(base_size = 14) +
  theme(
    strip.background   = element_rect(fill = "white"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position    = "right",
    plot.margin        = margin(5.5, 20, 5.5, 200),
    axis.text.y        = element_text(hjust = 1)
  ) +

  labs(
    x = "-log10(p-value)",
    y = "GO Biological Process"
  )

# 4. Save PDF ---------------------------------------------------------------
ggsave(
  filename = out_pdf,
  plot     = p,
  width    = 12,
  height   = 5,
  units    = "in"
)
