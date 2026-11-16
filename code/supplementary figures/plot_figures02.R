# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3, scales 1.4.0
#
# Inputs: VIP-PAD MOFA diagnostics, factor stability, scores and clustering results.
#
# Outputs: ./figures-output/FigureS02A.pdf through FigureS02F.pdf
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})



# ==============================================================================
# FIGURE S02 (PANEL A) — Read model diagnostics
# ==============================================================================

info <- fread("./results/VIP-PAD/mofa/Model_Diagnostics.csv")
mt <- fread(
  "./results/VIP-PAD/clustering/Cluster_Stability_Metrics.csv"
)



# ==============================================================================
# FIGURE S02 (PANEL A) — Model selection
# ==============================================================================
# 01. Model selection across 100 random initializations.
blue <- "#326A91"
grey <- "#989FA6"
theme_set(
  theme_classic(base_family = "Arial", base_size = 9) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 8),
      axis.title = element_text(size = 9),
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(
        size = 10,
        face = "bold",
        margin = margin(b = 5)
      ),
      plot.subtitle = element_text(size = 8, margin = margin(b = 7)),
      legend.title = element_text(size = 8),
      legend.text = element_text(size = 8),
      legend.key.height = grid::unit(3, "mm"),
      plot.margin = margin(7, 8, 6, 7),
      legend.position = "bottom"
    )
)
info[,
  selection := factor(
    ifelse(selected, "Selected model", "Other runs"),
    levels = c("Other runs", "Selected model")
  )
]
p <- ggplot(info, aes(seed, delta_ELBO_from_best, fill = selection)) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    colour = "#777777"
  ) +
  geom_point(shape = 21, size = 2.5, stroke = 0.3, colour = "black") +
  scale_fill_manual(values = c(`Other runs` = grey, `Selected model` = blue)) +
  geom_text(
    data = info[selected == TRUE],
    aes(label = paste("Seed", seed)),
    vjust = -1,
    hjust = 0.5,
    family = "Arial",
    size = 2.7
  ) +
  scale_x_continuous(breaks = c(42, seq(52, 132, 10), 141)) +
  labs(
    x = "Random seed",
    y = "ELBO gap from best",
    fill = NULL,
    title = "Selection by model fit",
    subtitle = sprintf(
      "Highest ELBO among converged models: seed %d",
      info[selected == TRUE, seed]
    )
  )
p <- p +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(sigma = 100),
    breaks = c(0, 100, 1000, 10000, 50000, 1e+05),
    labels = scales::label_comma()
  )
p <- p +
  labs(title = NULL, subtitle = NULL, caption = NULL) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    plot.caption = element_blank(),
    panel.border = element_rect(colour = "#333333", fill = NA, linewidth = 0.35)
  )



# ==============================================================================
# FIGURE S02 (PANEL A) — Save plot
# ==============================================================================

dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
ggsave(
  "./figures-output/FigureS02A.pdf",
  p,
  device = grDevices::cairo_pdf,
  width = 125,
  height = 95,
  units = "mm",
  dpi = 300,
  bg = "white"
)



# ==============================================================================
# FIGURE S02 (PANEL D) — Cluster stability
# ==============================================================================
theme_set(
  theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 14),
      axis.title = element_text(size = 15),
      axis.line = element_line(linewidth = 0.55),
      axis.ticks = element_line(linewidth = 0.45),
      axis.ticks.length = grid::unit(2.2, "mm"),
      plot.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.key.height = grid::unit(4, "mm"),
      panel.spacing.x = grid::unit(9, "mm"),
      panel.spacing.y = grid::unit(7, "mm"),
      plot.margin = margin(8, 10, 7, 8),
      legend.position = "bottom"
    )
)
mt <- melt(
  mt,
  id.vars = "K",
  measure.vars = c("silhouette", "PAC", "ARI_median", "min_Jaccard")
)
lab <- c(
  silhouette = "Silhouette",
  PAC = "PAC",
  ARI_median = "Median adjusted Rand index",
  min_Jaccard = "Minimum mean Jaccard"
)
mt[, metric := factor(lab[as.character(variable)], levels = lab)]
ps <- ggplot(mt, aes(K, value)) +
  geom_line(colour = "black") +
  geom_point(
    aes(fill = factor(K == 2, levels = c(FALSE, TRUE))),
    shape = 21,
    size = 3,
    stroke = 0.3,
    colour = "black"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 2:5, expand = expansion(add = 0.2)) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.1))) +
  scale_fill_manual(
    values = c(`FALSE` = "#989FA6", `TRUE` = "#326A91"),
    labels = c("Other K", "Selected K = 2"),
    name = NULL
  ) +
  labs(x = "Number of groups", y = NULL) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )
ggsave(
  "./figures-output/FigureS02D.pdf",
  ps,
  width = 8.4,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S02 (PANEL B) — Factor stability
# ==============================================================================
runs <- fread(
  "./results/VIP-PAD/mofa/Factor_Stability_All_Runs.csv"
)
factor_levels <- paste0("Factor", 1:10)
plotted <- runs[matched == TRUE & is.finite(score_abs_r)]
plotted[, factor := factor(factor, levels = factor_levels)]

plot <- ggplot(plotted, aes(factor, score_abs_r)) +
  geom_boxplot(outlier.size = 0.9, fill = "#BBDEEC", linewidth = 0.4) +
  coord_cartesian(ylim = c(0, 1)) +
  scale_y_continuous(
    breaks = seq(0, 1, 0.25),
    labels = scales::label_number(accuracy = 0.01)
  ) +
  labs(x = "VIP-PAD MOFA factors", y = "Matched absolute score correlation") +
  theme_classic(base_size = 14, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.text.x = element_text(size = 13, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 15),
    axis.line = element_line(linewidth = 0.55),
    axis.ticks = element_line(linewidth = 0.45),
    plot.margin = margin(8, 12, 9, 9)
  )



# ==============================================================================
# FIGURE S02 (PANEL B) — Save plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS02B.pdf",
  plot,
  width = 195,
  height = 122,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S02 (PANEL C) — Factor correlations
# ==============================================================================
scores <- fread("./results/VIP-PAD/mofa/Factor_Scores_Raw.csv")
x <- as.matrix(scores[, ..factor_levels])
# Compare factor scores across participants.
cc <- cor(x, method = "pearson")

d <- CJ(col = 1:10, row = 1:10)[col >= row]
d[, r := cc[cbind(row, col)]]
labels <- data.table(
  col = (1:10) - 0.57,
  row = 1:10,
  label = paste("Factor", 1:10)
)
p <- ggplot(d, aes(col, 11 - row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", r)),
    family = "Arial",
    size = 2.4,
    colour = ifelse(abs(d$r) > 0.65, "white", "black")
  ) +
  geom_text(
    data = labels,
    aes(x = col, y = 11 - row, label = label),
    inherit.aes = FALSE,
    family = "Arial",
    size = 3,
    hjust = 1
  ) +
  scale_fill_gradient2(
    low = "#326A91",
    mid = "#F5F3F0",
    high = "#A94E39",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r",
    guide = guide_colourbar(
      barheight = grid::unit(33, "mm"),
      theme = theme(
        legend.frame = element_rect(colour = "#333333", linewidth = 0.35)
      )
    )
  ) +
  scale_x_continuous(
    breaks = 1:10,
    labels = paste("Factor", 1:10),
    position = "top",
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(xlim = c(0.5, 10.5), ylim = c(0.5, 10.5), clip = "off") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 9) +
  theme(
    text = element_text(colour = "black"),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = 0.5,
      size = 9,
      colour = "black",
      margin = margin(b = 6)
    ),
    axis.text.y = element_blank(),
    panel.border = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8),
    plot.margin = margin(9, 8, 7, 48)
  )



# ==============================================================================
# FIGURE S02 (PANEL C) — Save plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS02C.pdf",
  p,
  width = 150,
  height = 145,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S02 (PANEL E) — Consensus clustering
# ==============================================================================
x <- readRDS(
  "./results/VIP-PAD/clustering/Clustering_Stability.rds"
)
g <- x$assignments$K2
ord <- order(g, x$Z[, 1])
n <- length(g)
edge <- sum(g == 1) + 0.5
centers <- c((sum(g == 1) + 1) / 2, sum(g == 1) + (sum(g == 2) + 1) / 2)
cm <- x$consensus[["2"]][ord, ord]
dd <- data.table(
  x = rep(seq_len(n), each = n),
  y = rep(seq_len(n), n),
  value = as.vector(cm)
)
base <- theme_void(base_family = "Noto Sans", base_size = 15) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 15),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 12, 8, 10),
    legend.box.margin = margin(8, 0, 0, 0)
  )
p1 <- ggplot(dd, aes(x, y, fill = value)) +
  geom_raster() +
  geom_vline(xintercept = edge, colour = "white", linewidth = 0.35) +
  geom_hline(yintercept = edge, colour = "white", linewidth = 0.35) +
  scale_x_continuous(
    breaks = centers,
    labels = c("Cluster 1", "Cluster 2"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = centers,
    labels = c("Cluster 1", "Cluster 2"),
    expand = c(0, 0)
  ) +
  scale_fill_gradient(
    low = "white",
    high = "#315E86",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    name = "Co-clustering"
  ) +
  coord_equal() +
  base +
  guides(
    fill = guide_colourbar(
      barwidth = unit(48, "mm"),
      barheight = unit(3.5, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  )

ggsave(
  "./figures-output/FigureS02E.pdf",
  p1,
  width = 5.3,
  height = 5.5,
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S02 (PANEL F) — Factor profiles
# ==============================================================================
# Standardize factor scores for the heatmap.
zz <- as.data.table(as.table(t(x$Z_display[ord, , drop = FALSE])))
setnames(zz, c("Factor", "Patient", "value"))
zz[, position := rep(seq_len(n), each = 10)]
zz[, Factor := factor(Factor, levels = rev(colnames(x$Z_display)))]
p2 <- ggplot(zz, aes(position, Factor, fill = value)) +
  geom_raster() +
  geom_vline(xintercept = edge, colour = "white", linewidth = 0.4) +
  scale_x_continuous(
    breaks = centers,
    labels = c("Cluster 1", "Cluster 2"),
    expand = c(0, 0)
  ) +
  scale_y_discrete(expand = c(0, 0), labels = function(z) {
    sub("^(Factor|F0)", "F", z)
  }) +
  scale_fill_gradient2(
    low = "#315E86",
    mid = "white",
    high = "#A94743",
    limits = c(-2.5, 2.5),
    oob = scales::squish,
    breaks = c(-2, 0, 2),
    name = "Z score"
  ) +
  base +
  theme(axis.text.y = element_text(size = 15, margin = margin(r = 7))) +
  guides(
    fill = guide_colourbar(
      barwidth = unit(48, "mm"),
      barheight = unit(3.5, "mm"),
      title.position = "top",
      title.hjust = 0.5
    )
  )

ggsave(
  "./figures-output/FigureS02F.pdf",
  p2,
  width = 8,
  height = 4.4,
  device = cairo_pdf,
  bg = "white"
)
