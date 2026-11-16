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
# Inputs: SingleFactor_OOF.rds and Results.rds from VIP-PAD 12-month prediction.
#
# Outputs: ./figures-output/ (FigureS25A, FigureS25B, FigureS25C_MACCE, FigureS25C_MALE)
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})



# ==============================================================================
# FIGURE S25 (PANEL C) — Risk reclassification
# ==============================================================================
# Function to draw event and non-event reclassification tables.
draw_nri <- function(z) {
  grid.newpage()
  # Function to position text in panel coordinates.
  txt <- function(s, x, y, sz = 12, b = FALSE, just = 'centre') {
    grid.text(
      s,
      x,
      y,
      just = just,
      gp = gpar(
        fontfamily = 'sans',
        fontsize = sz,
        fontface = if (b) 'bold' else 'plain',
        col = '#282828'
      )
    )
  }
  txt(z$endpoint, .09, .905, 20, TRUE, 'left')
  th <- z$threshold * 100
  for (k in 1:2) {
    v <- z$reclassification[status == c('Event', 'Non-event')[k]]
    m <- matrix(0, 2, 2)
    for (i in 1:2) {
      for (j in 1:2) {
        m[i, j] <- v[old == i & new == j]$IPCW_weighted_n
      }
    }
    top <- if (k == 1) .83 else .46
    txt(
      if (k == 1) 'Among events' else 'Among non-events',
      .39,
      top + .035,
      14,
      TRUE
    )
    txt('Clinical + Factors', .39, top - .012, 14, TRUE)
    txt('Clinical', .13, top - .075, 12, TRUE)
    txt(paste0('<', th, '%'), .30, top - .075, 12, TRUE)
    txt(paste0('\u2265', th, '%'), .46, top - .075, 12, TRUE)
    txt('Total', .59, top - .075, 12, TRUE)
    for (i in 1:2) {
      yy <- top - .13 - (i - 1) * .065
      txt(paste0(if (i == 1) '<' else '\u2265', th, '%'), .13, yy, 12, TRUE)
      for (j in 1:2) {
        fill <- if (i == j) {
          'white'
        } else if ((k == 1 && j > i) || (k == 2 && j < i)) {
          '#A6C8E0'
        } else {
          '#FBC892'
        }
        grid.rect(
          x = c(.30, .46)[j],
          y = yy,
          width = .16,
          height = .065,
          gp = gpar(fill = fill, col = '#62676B', lwd = .8)
        )
        txt(sprintf('%.1f', m[i, j]), c(.30, .46)[j], yy)
      }
      txt(sprintf('%.1f', sum(m[i, ])), .59, yy)
    }
    yy <- top - .26
    txt('Total', .13, yy, 12, TRUE)
    txt(sprintf('%.1f', sum(m[, 1])), .30, yy)
    txt(sprintf('%.1f', sum(m[, 2])), .46, yy)
    txt(sprintf('%.1f', sum(m)), .59, yy)
    a <- if (k == 1) m[1, 2] else m[2, 1]
    b <- if (k == 1) m[2, 1] else m[1, 2]
    txt(
      if (k == 1) 'Event NRI' else 'Non-event NRI',
      .67,
      top - .075,
      12,
      TRUE,
      'left'
    )
    txt(
      sprintf('= (%.1f - %.1f) / %.1f', a, b, sum(m)),
      .67,
      top - .14,
      12,
      FALSE,
      'left'
    )
    txt(sprintf('= %.3f', z$nri[k]), .67, top - .205, 12, FALSE, 'left')
  }
  txt(
    bquote(bold('Total NRI') ~ '=' ~ .(sprintf('%.3f', z$nri['Total']))),
    .67,
    .145,
    12,
    FALSE,
    'left'
  )
  pp <- z$bootstrap[metric == 'NRI_Total']$P
  txt(
    if (pp < .001) {
      expression(italic(P) < 0.001)
    } else {
      bquote(italic(P) == .(sprintf('%.3f', pp)))
    },
    .67,
    .095,
    12,
    FALSE,
    'left'
  )
}
# Function to separate risk-category changes for participants with and without events.
save_nri <- function(result, path) {

  cairo_pdf(path, width = 9, height = 5.3)
  on.exit(dev.off())
  draw_nri(result)
}



# ==============================================================================
# FIGURE S25 (PANELS A–B) — Single-factor discrimination
# ==============================================================================
# Function to draw each ROC curve from out-of-fold predictions.
plot_factor_roc <- function(z, ep) {
  factors <- z$metrics$factor
  g <- ggplot(z$roc[factor %in% factors], aes(FPR, TPR, colour = factor)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      linewidth = .4,
      colour = 'black'
    ) +
    geom_path(linewidth = .65) +
    scale_colour_manual(values = cols) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 1, .25),
      labels = scales::label_percent()
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, .25),
      labels = scales::label_percent()
    ) +
    labs(title = ep, x = '1 - Specificity', y = 'Sensitivity') +
    base_theme +
    theme(legend.position = 'none')
  ys <- .05 + .065 * rev(seq_along(factors) - 1)
  for (j in seq_along(factors)) {
    av <- z$metrics[factor == factors[j]]$AUC
    g <- g +
      annotate(
        'segment',
        x = .54,
        xend = .58,
        y = ys[j],
        yend = ys[j],
        colour = cols[factors[j]],
        linewidth = .8
      ) +
      annotate(
        'text',
        x = .60,
        y = ys[j],
        label = sprintf('%s  AUC %.3f', factors[j], av),
        hjust = 0,
        size = 6.2
      )
  }

  g
}



# ==============================================================================
# FIGURE S25 (PANELS A–C) — Prepare prediction results
# ==============================================================================

base_theme <- theme_classic(base_size = 22) +
  theme(
    panel.border = element_rect(fill = NA, colour = 'black', linewidth = .5),
    axis.line = element_blank(),
    axis.text = element_text(colour = 'black', size = 20),
    plot.title = element_text(size = 22),
    plot.margin = margin(10, 10, 10, 10)
  )

cols <- c(
  F1 = '#5878A6',
  F2 = '#D9808C',
  F3 = '#58A091',
  F4 = '#9B7AAF',
  F8 = '#CE913B'
)
single <- readRDS("./results/prediction/SingleFactor_OOF.rds")
performance <- readRDS("./results/prediction/Results.rds")



# ==============================================================================
# FIGURE S25 (PANELS A–C) — Draw and save panels
# ==============================================================================

dir.create("./figures-output", recursive = TRUE, showWarnings = FALSE)
# FIGURE S25 (PANEL A) — MALE discrimination
ggsave(
  "./figures-output/FigureS25A.pdf",
  plot_factor_roc(single$MALE, "MALE"),
  width = 7.5,
  height = 7,
  device = cairo_pdf
)
# FIGURE S25 (PANEL B) — MACCE discrimination
ggsave(
  "./figures-output/FigureS25B.pdf",
  plot_factor_roc(single$MACCE, "MACCE"),
  width = 7.5,
  height = 7,
  device = cairo_pdf
)
# FIGURE S25 (PANEL C) — MALE and MACCE reclassification
save_nri(performance$MALE, "./figures-output/FigureS25C_MALE.pdf")
save_nri(performance$MACCE, "./figures-output/FigureS25C_MACCE.pdf")
