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
# Inputs: Results.rds and Cluster_NRI.rds from VIP-PAD 12-month prediction.
#
# Outputs: ./figures-output/ (Figure6A_Calibration, Figure6A_DCA, Figure6A_ROC,
#          Figure6B_Calibration, Figure6B_DCA, Figure6B_ROC, Figure6C_MACCE, Figure6C_MALE)
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
# FUNCTIONS
# ==============================================================================

# Function to write the panel to its figure-specific PDF.
export_plot <- function(g, path) {

  ggsave(path, g, width = 7.5, height = 7, device = cairo_pdf)
}
# Function to draw endpoint-specific ROC curves.
plot_roc <- function(z) {
  g <- ggplot(z$roc, aes(FPR, TPR, colour = model)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      linewidth = .4,
      colour = 'black'
    ) +
    geom_line(linewidth = .65) +
    scale_colour_manual(values = pal) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(
      breaks = seq(0, 1, .25),
      labels = scales::label_percent()
    ) +
    scale_y_continuous(
      breaks = seq(0, 1, .25),
      labels = scales::label_percent()
    ) +
    labs(title = z$endpoint, x = '1 - Specificity', y = 'Sensitivity') +
    base_theme +
    theme(legend.position = 'none')
  ys <- c(.18, .115, .05)
  av <- z$metrics$AUC[match(mods, z$metrics$model)]
  for (j in 1:3) {
    g <- g +
      annotate(
        'segment',
        x = .195,
        xend = .235,
        y = ys[j],
        yend = ys[j],
        colour = pal[j],
        linewidth = .8
      ) +
      annotate(
        'text',
        x = .255,
        y = ys[j],
        label = sprintf('%s  AUC %.3f', mods[j], av[j]),
        hjust = 0,
        size = 6.2
      )
  }
  pv <- z$bootstrap$P[match(c('DeltaAUC_C', 'DeltaAUC_F'), z$bootstrap$metric)]
  stars <- ifelse(
    pv < .001,
    '***',
    ifelse(pv < .01, '**', ifelse(pv < .05, '*', 'ns'))
  )
  for (j in 1:2) {
    bx <- c(.89, .945)[j]
    bot <- ys[j + 1]
    g <- g +
      annotate(
        'segment',
        x = bx,
        xend = bx,
        y = bot,
        yend = ys[1],
        linewidth = .45
      ) +
      annotate(
        'segment',
        x = bx - .012,
        xend = bx,
        y = bot,
        yend = bot,
        linewidth = .45
      ) +
      annotate(
        'segment',
        x = bx - .012,
        xend = bx,
        y = ys[1],
        yend = ys[1],
        linewidth = .45
      ) +
      annotate(
        'text',
        x = bx + .026,
        y = mean(c(bot, ys[1])),
        label = stars[j],
        angle = 90,
        size = 8
      )
  }
  g
}
# Function to compare predicted and observed risks by risk group.
plot_calibration <- function(
  z,
  show_ci = TRUE,
  limit = NULL,
  metric = "brier"
) {
  d <- z$calibration
  m <- z$metrics[match(mods, model)]
  labels <- if (metric == 'brier_inline') {
    sprintf(
      '%s  Brier %s',
      m$model,
      formatC(m$Brier, format = 'fg', digits = 3, flag = '#')
    )
  } else if (metric == 'brier') {
    sprintf('%s\nBrier %.4f', m$model, m$Brier)
  } else {
    sprintf('%s\nIntercept %.3f   Slope %.3f', m$model, m$Intercept, m$Slope)
  }
  lim <- if (!is.null(limit)) {
    limit
  } else {
    max(.3, ceiling(max(d$upper, d$predicted) * 10) / 10)
  }
  ggplot(d, aes(predicted, observed, colour = model, group = model)) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = 2,
      linewidth = .4,
      colour = 'black'
    ) +
    {
      if (show_ci) {
        geom_errorbar(
          aes(ymin = lower, ymax = upper),
          width = .006,
          linewidth = .4,
          alpha = .65
        )
      }
    } +
    geom_line(linewidth = .65) +
    geom_point(size = 3) +
    scale_colour_manual(values = pal, breaks = mods, labels = labels) +
    scale_x_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, lim, .1)
    ) +
    scale_y_continuous(
      labels = scales::label_percent(accuracy = 1),
      breaks = seq(0, lim, .1)
    ) +
    coord_equal(xlim = c(0, lim), ylim = c(0, lim), expand = FALSE) +
    labs(
      title = z$endpoint,
      x = 'Predicted 12-month risk',
      y = 'Observed 12-month risk',
      colour = NULL
    ) +
    base_theme +
    theme(
      legend.position = 'inside',
      legend.position.inside = if (z$endpoint == 'MALE') {
        c(.995, .005)
      } else {
        c(.97, .04)
      },
      legend.justification = c(1, 0),
      legend.text = element_text(size = 16, lineheight = 1.15),
      legend.key.spacing.y = unit(8, 'pt'),
      legend.background = element_blank()
    )
}
# Function to plot net benefit across decision thresholds.
plot_dca <- function(z) {
  levels <- c(mods, 'Treat all', 'Treat none')
  colors <- c(pal, 'Treat all' = '#AAAAAA', 'Treat none' = '#222222')
  top <- ceiling(max(z$dca$net_benefit) * 20) / 20
  ggplot(z$dca, aes(threshold, net_benefit, colour = model, linetype = model)) +
    geom_line(linewidth = .65) +
    scale_colour_manual(values = colors, breaks = levels) +
    scale_linetype_manual(
      values = setNames(
        c('solid', 'solid', 'solid', 'dashed', 'dotted'),
        levels
      ),
      breaks = levels
    ) +
    scale_x_continuous(
      breaks = seq(0, .3, .05),
      labels = scales::label_percent(accuracy = 1)
    ) +
    scale_y_continuous(breaks = seq(0, top, .05)) +
    coord_cartesian(xlim = c(0, .30), ylim = c(-.015, top), expand = FALSE) +
    labs(
      title = z$endpoint,
      x = 'Threshold probability',
      y = 'Net benefit',
      colour = NULL,
      linetype = NULL
    ) +
    base_theme +
    theme(
      legend.position = 'inside',
      legend.position.inside = c(.97, .97),
      legend.justification = c(1, 1),
      plot.margin = margin(10, 24, 10, 10),
      legend.text = element_text(size = 17),
      legend.key.spacing.y = unit(4, 'pt'),
      legend.background = element_blank()
    )
}
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
    txt('Clinical + C1/C2', .39, top - .012, 14, TRUE)
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
# Function to export the reclassification panel.
save_nri <- function(result, path) {

  cairo_pdf(path, width = 9, height = 5.3)
  on.exit(dev.off())
  draw_nri(result)
}



# ==============================================================================
# FIGURE 6 (PANELS A–C) — Read prediction results
# ==============================================================================

mods <- c('Clinical', 'Clinical + C1/C2', 'Clinical + Factors')
pal <- setNames(c('#666666', '#D9808C', '#5878A6'), mods)
base_theme <- theme_classic(base_size = 22) +
  theme(
    panel.border = element_rect(fill = NA, colour = 'black', linewidth = .5),
    axis.line = element_blank(),
    axis.text = element_text(colour = 'black', size = 20),
    plot.title = element_text(size = 22),
    plot.margin = margin(10, 10, 10, 10)
  )

performance <- readRDS("./results/prediction/Results.rds")
cluster_nri <- readRDS("./results/prediction/Cluster_NRI.rds")



# ==============================================================================
# FIGURE 6 (PANEL A) — MALE discrimination, calibration and decision curve
# ==============================================================================
dir.create("./figures-output", recursive = TRUE, showWarnings = FALSE)
export_plot(plot_roc(performance$MALE), "./figures-output/Figure6A_ROC.pdf")
export_plot(
  plot_calibration(performance$MALE, metric = "brier_inline"),
  "./figures-output/Figure6A_Calibration.pdf"
)
export_plot(plot_dca(performance$MALE), "./figures-output/Figure6A_DCA.pdf")



# ==============================================================================
# FIGURE 6 (PANEL B) — MACCE discrimination, calibration and decision curve
# ==============================================================================
export_plot(plot_roc(performance$MACCE), "./figures-output/Figure6B_ROC.pdf")
export_plot(
  plot_calibration(performance$MACCE, metric = "brier_inline"),
  "./figures-output/Figure6B_Calibration.pdf"
)
export_plot(plot_dca(performance$MACCE), "./figures-output/Figure6B_DCA.pdf")



# ==============================================================================
# FIGURE 6 (PANEL C) — Reclassification after adding cluster membership
# ==============================================================================
save_nri(cluster_nri$MALE, "./figures-output/Figure6C_MALE.pdf")
save_nri(cluster_nri$MACCE, "./figures-output/Figure6C_MACCE.pdf")
