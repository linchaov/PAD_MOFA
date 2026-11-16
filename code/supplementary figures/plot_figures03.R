# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggh4x 0.3.1, ggplot2 4.0.3, scales 1.4.0
#
# Inputs: UKB MOFA diagnostics, factor stability, variance, scores and clustering results.
#
# Outputs: ./figures-output/FigureS03A.pdf through FigureS03E.pdf
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
# FIGURE S03 (PANEL A) — Read model diagnostics
# ==============================================================================

info <- fread("./results/UKB/mofa/Model_Diagnostics.csv")



# ==============================================================================
# FIGURE S03 (PANEL A) — Model selection
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
    subtitle = paste(
      "Highest ELBO among converged models: seed",
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
# FIGURE S03 (PANEL A) — Save plot
# ==============================================================================

dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
ggsave(
  "./figures-output/FigureS03A.pdf",
  p,
  device = grDevices::cairo_pdf,
  width = 125,
  height = 95,
  units = "mm",
  dpi = 300,
  bg = "white"
)



# ==============================================================================
# FIGURE S03 (PANEL C) — Factor correlations
# ==============================================================================
scores <- fread("./results/UKB/mofa/Factor_Scores_Raw.csv")
factor_cols <- paste0("Factor", 1:11)
x <- as.matrix(scores[, ..factor_cols])
# Compare factor scores across participants.
cc <- cor(x, method = "pearson")

n_factors <- length(factor_cols)
d <- CJ(col = seq_len(n_factors), row = seq_len(n_factors))[col >= row]
d[, r := cc[cbind(row, col)]]
labels <- data.table(
  col = seq_len(n_factors) - 0.57,
  row = seq_len(n_factors),
  label = paste("Factor", seq_len(n_factors))
)
p <- ggplot(d, aes(col, n_factors + 1 - row, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(
    aes(label = sprintf("%.2f", r)),
    family = "Arial",
    size = 2.4,
    colour = ifelse(abs(d$r) > 0.65, "white", "black")
  ) +
  geom_text(
    data = labels,
    aes(x = col, y = n_factors + 1 - row, label = label),
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
    breaks = seq_len(n_factors),
    labels = paste("Factor", seq_len(n_factors)),
    position = "top",
    expand = c(0, 0)
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_fixed(
    xlim = c(0.5, n_factors + 0.5),
    ylim = c(0.5, n_factors + 0.5),
    clip = "off"
  ) +
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
# FIGURE S03 (PANEL C) — Save plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS03C.pdf",
  p,
  width = 165,
  height = 155,
  units = "mm",
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S03 (PANEL E) — Cluster stability
# ==============================================================================
metrics <- fread("./results/UKB/clustering/Cluster_Stability_Metrics.csv")
selected_k <- metrics[selected == TRUE, K]
theme_set(
  theme_classic(base_size = 12, base_family = "Arial") +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 14),
      axis.title = element_text(size = 15),
      axis.line = element_line(linewidth = 0.55),
      axis.ticks = element_line(linewidth = 0.45),
      plot.title = element_text(size = 13, face = "bold"),
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 13),
      legend.key.height = grid::unit(4, "mm"),
      plot.margin = margin(7, 8, 6, 7),
      legend.position = "bottom"
    )
)
mt <- melt(
  metrics,
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
    aes(fill = factor(K == selected_k, levels = c(FALSE, TRUE))),
    shape = 21,
    size = 3,
    stroke = 0.3,
    colour = "black"
  ) +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = 2:5) +
  scale_fill_manual(
    values = c(`FALSE` = "#989FA6", `TRUE` = "#326A91"),
    labels = c("Other K", paste("Selected K =", selected_k)),
    name = NULL
  ) +
  labs(x = "Number of groups", y = NULL) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14)
  )
ps <- ps +
  ggh4x::facetted_pos_scales(
    y = list(
      metric == "Minimum mean Jaccard" ~ scale_y_continuous(
        limits = c(0.7, 1),
        breaks = seq(0.7, 0.95, 0.05),
        labels = scales::label_number(accuracy = 0.01),
        expand = expansion(mult = c(0.02, 0.02))
      )
    )
  )



# ==============================================================================
# FIGURE S03 (PANEL E) — Save plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS03E.pdf",
  ps,
  width = 8.4,
  height = 6,
  device = cairo_pdf,
  bg = "white"
)



# ==============================================================================
# FIGURE S03 (PANELS D AND B) — Set plot theme
# ==============================================================================

theme_set(
  theme_classic(base_size = 9, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 10),
      legend.title = element_blank(),
      plot.margin = margin(6, 8, 6, 6)
    )
)



# ==============================================================================
# FIGURE S03 (PANEL D) — Variance explained
# ==============================================================================
v <- fread("./results/UKB/mofa/Variance_Explained.csv")
v <- melt(v, id.vars = "factor", variable.name = "view", value.name = "R2")
v[, factor := factor(factor, levels = unique(factor))]
p <- ggplot(v, aes(factor, R2, fill = view)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.75) +
  geom_text(
    aes(label = ifelse(R2 < 0.01, "<0.01", sprintf("%.2f", R2))),
    position = position_dodge(width = 0.75),
    vjust = -0.5,
    size = 3.3,
    colour = "black",
    show.legend = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_fill_manual(
    values = c(Protein = "#4575B4", Metabolite = "#F46D43"),
    labels = c(Protein = "Olink protein", Metabolite = "Nightingale NMR")
  ) +
  labs(x = "UKB MOFA factors", y = "Variance explained (%)", title = NULL) +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12, colour = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )



# ==============================================================================
# FIGURE S03 (PANEL D) — Save plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS03D.pdf",
  p,
  width = 210,
  height = 125,
  units = "mm",
  device = cairo_pdf
)



# ==============================================================================
# FIGURE S03 (PANEL B) — Factor stability
# ==============================================================================
s <- fread("./results/UKB/mofa/Factor_Stability_All_Runs.csv")
s[, factor := factor(factor, levels = unique(factor))]
s <- s[matched == TRUE & is.finite(score_abs_r)]
p <- ggplot(s, aes(factor, score_abs_r)) +
  geom_boxplot(outlier.size = 0.6, fill = "#BBDEEC", linewidth = 0.3) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(x = "UKB MOFA factors", y = "Matched absolute score correlation") +
  theme(
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 13)
  )
ggsave(
  "./figures-output/FigureS03B.pdf",
  p,
  width = 180,
  height = 110,
  units = "mm",
  device = cairo_pdf
)
