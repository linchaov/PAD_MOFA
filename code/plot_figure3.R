# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: circlize 0.4.18, data.table 1.18.6.1, ggplot2 4.0.3, patchwork 1.3.2,
#                  RColorBrewer 1.1-3, scales 1.4.0, stringi 1.8.9, survival 3.8-6
#
# Inputs: Cross-cohort pathway and loading comparisons; UKB Cox models, cluster characteristics
#         and
#          survival estimates.
#
# Outputs: ./figures-output/Figure3A.pdf through Figure3G.pdf
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(patchwork)
  library(circlize)
  library(RColorBrewer)
  library(survival)
})



# ==============================================================================
# FIGURE 3 (PANEL A) — Protein pathway enrichment across cohorts
# ==============================================================================
# Function to plot protein pathway enrichment across cohorts.
plot_3A <- function() {
  cohorts <- c("VIP-PAD", "UKB")
  A <- lapply(cohorts, function(co) {
    fread(sprintf(
      "./results/cross-cohort/enrichment/%s/Displayed_Terms.csv",
      co
    ))
  })
  names(A) <- cohorts
  S <- lapply(cohorts, function(co) {
    fread(sprintf(
      "./results/cross-cohort/enrichment/%s/Signed_Displayed_Terms.csv",
      co
    ))
  })
  names(S) <- cohorts
  mx <- max(vapply(
    c(S, A),
    function(x) max(-log10(pmax(x$padj, .Machine$double.xmin))),
    numeric(1)
  ))
  upper <- max(50, ceiling(mx / 10) * 10)
  colfun <- colorRamp2(
    c(0, 2, 5, 10, 20, 35, upper),
    c("white", brewer.pal(9, "YlOrRd")[4:9])
  )
  themes <- c(
    "Inflammation & immunity",
    "Complement",
    "Coagulation & platelets",
    "Lipid metabolism",
    "Endothelium & angiogenesis",
    "ECM & remodeling",
    "Metabolism & oxidative stress",
    "Cell injury & repair",
    "Other"
  )
  labels <- c(
    "Inflammation\n& immunity",
    "Complement",
    "Coagulation\n& platelets",
    "Lipid\nmetabolism",
    "Endothelium\n& angiogenesis",
    "ECM &\nremodeling",
    "Metabolism &\noxidative stress",
    "Cell injury\n& repair",
    "Other"
  )
  W <- 300
  H <- 248.5
  left <- 15
  pw <- 29
  gap <- 1.2
  panel_height <- 80
  # Function to draw pathway enrichment by cohort and factor.
  plot <- function() {
    grid.newpage()
    pushViewport(viewport(xscale = c(0, W), yscale = c(0, H)))
    # Function to position text in panel coordinates.
    tx <- function(t, x, y, size = 11, just = "centre", rot = 0, bold = FALSE) {
      grid.text(
        t,
        x = unit(x, "native"),
        y = unit(y, "native"),
        just = just,
        rot = rot,
        gp = gpar(
          fontfamily = "Arial",
          fontsize = size *
            1.25,
          fontface = if (bold) {
            "bold"
          } else {
            "plain"
          },
          lineheight = 1.1
        )
      )
    }
    for (co in cohorts) {
      K <- nrow(fread(file.path("./results", co, "mofa/Variance_Explained.csv"),
                      select = "factor"))
      rh <- panel_height / K
      top <- if (co == "VIP-PAD") {
        230.5
      } else {
        134.5
      }
      tx(paste(co, "Absolute loadings"), left, top + 7, 13, "left", bold = TRUE)
      for (k in 1:K) {
        tx(paste0("F", k), left - 2, top - (k - 0.5) * rh, 10.5, "right")
      }
      for (i in seq_along(themes)) {
        th <- themes[i]
        terms <- stringi::stri_sort(
          unique(c(
            S[[co]][theme == th, pathway],
            A[[co]][theme == th, pathway]
          )),
          locale = "en_US"
        )
        if (!length(terms)) {
          terms <- "NO_SELECTED_TERM"
        }
        m <- matrix(0, K, length(terms))
        d <- A[[co]][theme == th]
        m[cbind(d$factor, match(d$pathway, terms))] <- -log10(pmax(
          d$padj,
          .Machine$double.xmin
        ))
        x0 <- left + (i - 1) * (pw + gap)
        cw <- pw / ncol(m)
        for (k in 1:K) {
          grid.rect(
            x = unit(x0 + (seq_len(ncol(m)) - 0.5) * cw, "native"),
            y = unit(
              top -
                (k - 0.5) * rh,
              "native"
            ),
            width = unit(cw, "native"),
            height = unit(rh, "native"),
            gp = gpar(fill = colfun(m[k, ]), col = NA)
          )
        }
        grid.rect(
          x = unit(x0 + pw / 2, "native"),
          y = unit(top - K * rh / 2, "native"),
          width = unit(pw, "native"),
          height = unit(K * rh, "native"),
          gp = gpar(fill = NA, col = "#444444", lwd = 0.6)
        )
        if (co == "UKB") {
          tx(
            labels[i],
            x0 + pw / 2 + 6,
            top - K * rh - 4,
            12,
            just = "right",
            rot = 45
          )
        }
      }
    }
    # One shared legend above the panels; equal scale in both cohorts.
    lx <- 220
    ly <- 242.5
    lw <- 53
    for (i in 1:200) {
      grid.rect(
        x = unit(lx + (i - 0.5) * lw / 200, "native"),
        y = unit(ly, "native"),
        width = unit(lw / 200 + 0.01, "native"),
        height = unit(2.8, "native"),
        gp = gpar(fill = colfun((i - 1) / 199 * upper), col = NA)
      )
    }
    tx("-log10(adj. P)", lx - 3, ly, 10, "right")
    for (t in c(0, 10, 25, upper)) {
      tx(as.character(t), lx + t / upper * lw, ly - 4, 9)
    }
    popViewport()
  }

  cairo_pdf(
    "./figures-output/Figure3A.pdf",
    width = W / 25.4,
    height = H / 25.4
  )
  plot()
  dev.off()
}



# ==============================================================================
# FIGURE 3 (PANEL B) — Correlations of factor loadings
# ==============================================================================
# Function to plot correlations of factor loadings.
plot_3B <- function(x) {
  theme_set(
    theme_classic(base_family = "Arial", base_size = 9) +
      theme(
        text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8),
        legend.position = "bottom",
        plot.margin = margin(6, 8, 6, 6)
      )
  )
  d <- copy(x$correlations[view == "Protein"])
  d[, `:=`(
    ff = factor(
      fuwai_factor,
      levels = max(d$fuwai_factor):1,
      labels = paste0("F", max(d$fuwai_factor):1)
    ),
    uf = factor(ukb_factor, levels = 1:11, labels = paste0("F", 1:11))
  )]
  heat <- ggplot(d, aes(uf, ff, fill = pearson)) +
    geom_tile() +
    geom_text(
      aes(label = sprintf("%.2f", ifelse(abs(pearson) < 0.005, 0, pearson))),
      size = 3.5
    ) +
    scale_fill_gradient2(
      low = "#326A91",
      mid = "white",
      high = "#B84F50",
      midpoint = 0,
      limits = c(-1, 1),
      breaks = seq(-1, 1, 0.5),
      name = "Pearson r"
    ) +
    coord_equal() +
    labs(
      x = "UKB factors",
      y = "VIP-PAD factors",
      title = "a  Protein-loading correspondence",
      subtitle = "879 shared proteins; all 110 factor pairs"
    ) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      legend.key.width = grid::unit(12, "mm")
    )
  heat <- heat + labs(title = NULL, subtitle = NULL, y = "VIP-PAD factors")
  heat <- heat +
    theme(
      axis.text = element_text(size = 13, colour = "black"),
      axis.title = element_text(size = 15, colour = "black"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      plot.title = element_text(size = 17, face = "bold"),
      plot.margin = margin(2, 2, 2, 2)
    ) +
    guides(
      fill = guide_colourbar(
        barwidth = grid::unit(50, "mm"),
        barheight = grid::unit(5, "mm")
      )
    )

  ggsave(
    "./figures-output/Figure3B.pdf",
    heat,
    width = 200,
    height = 190,
    units = "mm",
    device = cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 3 (PANEL C) — Loading correlations across 879 shared proteins
# ==============================================================================
# Function to plot loading correlations across 879 shared proteins.
plot_3C <- function(comparison) {
  w <- comparison$loadings
  pairs <- comparison$candidates[, .(f = fuwai_factor, u = ukb_factor)]
  ds <- lapply(1:3, function(i) {
    data.table(
      x = w[[paste0("F", pairs$f[i])]],
      y = w[[paste0("U", pairs$u[i])]]
    )
  })
  xr <- range(unlist(lapply(ds, function(d) d$x)))
  yr <- range(unlist(lapply(ds, function(d) d$y)))
  xr <- xr + c(-1, 1) * diff(xr) * 0.07
  yr <- yr + c(-1, 1) * diff(yr) * 0.09
  xb <- seq(xr[1], xr[2], length.out = 49)
  yb <- seq(yr[1], yr[2], length.out = 49)
  # Count shared proteins within bins of paired factor loadings.
  bins <- lapply(ds, function(d) {
    b <- d[, .(
      ix = cut(x, xb, labels = FALSE, include.lowest = TRUE),
      iy = cut(y, yb, labels = FALSE, include.lowest = TRUE)
    )][, .N, by = .(ix, iy)]
    b[, `:=`(x = (xb[ix] + xb[ix + 1]) / 2, y = (yb[iy] + yb[iy + 1]) / 2)]
    b
  })
  maxn <- max(unlist(lapply(bins, function(b) b$N)))
  plots <- lapply(1:3, function(i) {
    d <- ds[[i]]
    f <- pairs$f[i]
    u <- pairs$u[i]
    ggplot(d, aes(x, y)) +
      geom_tile(
        data = bins[[i]],
        aes(fill = N),
        width = diff(xb)[1],
        height = diff(yb)[1]
      ) +
      scale_fill_gradientn(
        colours = c("#E3EFF5", "#8DBCD0", "#3D88A9", "#154D70"),
        trans = "log10",
        limits = c(1, maxn),
        breaks = c(1, 5, 10, 20, 100),
        name = "Proteins per bin"
      ) +
      geom_smooth(
        method = "lm",
        formula = y ~ x,
        se = FALSE,
        colour = "black",
        linewidth = 0.85
      ) +
      scale_x_continuous(
        limits = xr,
        breaks = scales::breaks_pretty(4),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        limits = yr,
        breaks = scales::breaks_pretty(4),
        expand = c(0, 0)
      ) +
      labs(
        title = paste0("VIP-PAD F", f, " vs UKB F", u),
        subtitle = sprintf(
          "r = %.3f   rho = %.3f",
          cor(d$x, d$y),
          cor(d$x, d$y, method = "spearman")
        ),
        x = paste0("VIP-PAD F", f, " loading"),
        y = paste0("UKB F", u, " loading")
      ) +
      theme_classic(base_size = 13) +
      theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 15, hjust = 0.5, colour = "black"),
        axis.title = element_text(size = 13),
        axis.text = element_text(size = 11, colour = "black"),
        axis.line = element_line(linewidth = 0.55),
        axis.ticks = element_line(linewidth = 0.5),
        plot.margin = margin(8, 10, 7, 7),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
      ) +
      guides(
        fill = guide_colourbar(
          barwidth = grid::unit(65, "mm"),
          barheight = grid::unit(4, "mm")
        )
      )
  })
  p <- wrap_plots(plots, ncol = 1, guides = "collect") &
    theme(legend.position = "bottom")

  ggsave(
    "./figures-output/Figure3C.pdf",
    p,
    width = 120,
    height = 285,
    units = "mm",
    device = cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 3 (PANEL D) — Correspondence between VIP-PAD and UKB factors
# ==============================================================================
# Function to plot correspondence between VIP-PAD and UKB factors.
plot_3D <- function(x) {
  theme_set(
    theme_classic(base_family = "Arial", base_size = 9) +
      theme(
        text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8),
        legend.position = "bottom",
        plot.margin = margin(6, 8, 6, 6)
      )
  )
  d <- copy(x$correlations[view == "Protein"])
  d[, `:=`(
    ff = factor(
      fuwai_factor,
      levels = max(d$fuwai_factor):1,
      labels = paste0("F", max(d$fuwai_factor):1)
    ),
    uf = factor(ukb_factor, levels = 1:11, labels = paste0("F", 1:11))
  )]
  ed <- d[display_link == TRUE]
  ed[, direction := ifelse(pearson < 0, "Negative", "Positive")]
  ed[, `:=`(y1 = 12 - fuwai_factor, y2 = 12 - ukb_factor)]
  nodes <- rbind(
    data.table(
      cohort = "VIP-PAD",
      xx = 0,
      yy = 12 - (seq_len(max(d$fuwai_factor))),
      label = paste0("F", seq_len(max(d$fuwai_factor)))
    ),
    data.table(
      cohort = "UKB",
      xx = 1,
      yy = 12 - (1:11),
      label = paste0("F", 1:11)
    )
  )
  net <- ggplot() +
    geom_segment(
      data = ed,
      aes(
        x = 0.12,
        xend = 0.88,
        y = y1,
        yend = y2,
        linewidth = abs_r,
        colour = direction
      ),
      alpha = 0.85
    ) +
    scale_colour_manual(
      values = c(Negative = "#326A91", Positive = "#B84F50"),
      name = "Correlation"
    ) +
    geom_point(data = nodes, aes(xx, yy), size = 5.3, colour = "#326A91") +
    geom_text(
      data = nodes[cohort == "VIP-PAD"],
      aes(xx - 0.08, yy, label = label),
      hjust = 1,
      size = 6.6
    ) +
    geom_text(
      data = nodes[cohort == "UKB"],
      aes(xx + 0.08, yy, label = label),
      hjust = 0,
      size = 6.6
    ) +
    scale_linewidth_continuous(range = c(1.2, 2.6), guide = "none") +
    scale_x_continuous(
      limits = c(-0.3, 1.3),
      breaks = c(0, 1),
      labels = c("VIP-PAD", "UKB"),
      position = "top"
    ) +
    scale_y_continuous(limits = c(0.85, 11.5), expand = c(0, 0)) +
    labs(
      title = "b  Descriptive links",
      subtitle = "Top 3 by magnitude; colours retain the sign",
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.line = element_blank(),
      axis.ticks = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank()
    )
  net <- net + labs(title = NULL, subtitle = NULL)
  net <- net +
    theme(
      axis.text.x = element_text(
        size = 20,
        colour = "black",
        face = "plain",
        hjust = 0.5,
        margin = margin(b = 10)
      ),
      plot.title = element_text(size = 18, face = "bold"),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 13),
      legend.key.width = grid::unit(14, "mm"),
      plot.margin = margin(2, 2, 2, 2)
    )
  net <- net + theme(legend.position = "none")

  ggsave(
    "./figures-output/Figure3D.pdf",
    net,
    width = 130,
    height = 205,
    units = "mm",
    device = cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 3 (PANEL E) — UKB factor associations with MALE and MACCE
# ==============================================================================
# Function to plot UKB factor associations with MALE and MACCE.
plot_3E <- function() {
  r <- fread("./results/UKB/outcomes/Cox/Associations.csv")
  r[,
    status := fifelse(
      FDR < 0.05,
      fifelse(HR < 1, "Lower", "Higher"),
      "Nonsignificant"
    )
  ]
  pal <- c(Lower = "#426F9D", Higher = "#BD5857", Nonsignificant = "#6B7380")
  theme_set(
    theme_classic(base_family = "Arial", base_size = 14) +
      theme(
        axis.text = element_text(colour = "black", size = 14),
        axis.text.y = element_text(size = 16),
        axis.title = element_text(size = 16),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 19),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "none",
        plot.margin = margin(8, 8, 8, 8)
      )
  )
  plots <- list()
  for (ep in c("MALE", "MACCE")) {
    dd <- r[endpoint == ep]
    dd[, row := 12 - factor]
    lim <- range(c(dd$CI_low, dd$CI_high, 1))
    lim <- exp(log(lim) + c(-1, 1) * 0.08 * diff(log(lim)))
    candidates <- c(
      0.1,
      0.2,
      0.25,
      0.33,
      0.5,
      0.75,
      1,
      1.25,
      1.5,
      2,
      2.5,
      3,
      4,
      5,
      7.5,
      10
    )
    br <- candidates[candidates >= lim[1] & candidates <= lim[2]]
    if (length(br) > 7) {
      br <- br[!br %in% c(0.33, 1.25, 2.5, 7.5)]
    }
    br <- br[br != 1.25]
    p <- ggplot(dd, aes(y = row, x = HR, colour = status)) +
      geom_vline(
        xintercept = 1,
        linetype = "dashed",
        colour = "#A8ADB7",
        linewidth = 0.35
      ) +
      geom_segment(
        aes(x = CI_low, xend = CI_high, yend = row),
        linewidth = 0.5
      ) +
      geom_segment(
        aes(x = CI_low, xend = CI_low, y = row - 0.11, yend = row + 0.11),
        linewidth = 0.5
      ) +
      geom_segment(
        aes(x = CI_high, xend = CI_high, y = row - 0.11, yend = row + 0.11),
        linewidth = 0.5
      ) +
      geom_point(aes(fill = status), shape = 21, size = 3.3, stroke = 0.65) +
      scale_colour_manual(values = pal) +
      scale_fill_manual(
        values = c(
          Lower = unname(pal["Lower"]),
          Higher = unname(pal["Higher"]),
          Nonsignificant = "white"
        )
      ) +
      scale_x_log10(limits = lim, breaks = br, labels = function(v) {
        vapply(
          v,
          function(z) format(z, trim = TRUE, scientific = FALSE),
          character(1)
        )
      }) +
      scale_y_continuous(
        breaks = 1:11,
        labels = paste0("F", 11:1),
        limits = c(0.4, 11.6),
        expand = c(0, 0)
      ) +
      labs(title = ep, subtitle = NULL, x = "Adjusted HR (95% CI)", y = NULL)
    if (ep == "MACCE") {
      p <- p + theme(axis.text.y = element_blank())
    }
    plots[[ep]] <- p
  }
  p <- plots$MALE | plots$MACCE

  ggsave(
    "./figures-output/Figure3E.pdf",
    p,
    width = 160,
    height = 155,
    units = "mm",
    device = cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 3 (PANEL F) — Clinical traits across UKB molecular clusters
# ==============================================================================
# Function to plot clinical traits across UKB molecular clusters.
plot_3F <- function() {
  d <- readRDS(
    "./results/UKB/clinical/Distributions.rds"
  )
  st <- fread(
    "./results/UKB/clinical/Descriptive_Statistics.csv"
  )
  st[, group := factor(group, levels = c("C1", "C2"))]
  variables <- c("Age", "BMI", "SBP", "LDL", "eGFR")
  tests <- fread(
    "./results/UKB/clinical/Wilcoxon_Tests.csv"
  )
  tests[,
    label := fcase(
      P < 0.001 , "***" , P < 0.01 , "**" , P < 0.05 , "*" ,
      default = ""
    )
  ]
  pal <- c(C1 = "#F4A7AE", C2 = "#5878A6")
  units <- c(
    Age = "years",
    BMI = "kg/m²",
    SBP = "mmHg",
    LDL = "mmol/L",
    eGFR = "mL/min/1.73 m²"
  )
  pp <- lapply(variables, function(v) {
    z <- d[variable == v]
    ss <- st[variable == v]
    rr <- range(z$value)
    pad <- diff(rr) * 0.08
    lims <- c(max(0, rr[1] - pad), rr[2] + diff(rr) * 0.22)
    p <- ggplot(z, aes(x = x, y = value, color = group, fill = group))
    # Estimate the displayed distribution separately for each cluster.
    den <- rbindlist(lapply(c("C1", "C2"), function(gg) {
      u <- z[group == gg]
      de <- density(u$value, from = min(u$value), to = max(u$value), n = 256)
      x0 <- as.numeric(u$group[1])
      data.table(
        group = gg,
        x = c(rep(x0 - 0.07, 256), rev(x0 - 0.07 - 0.34 * de$y / max(de$y))),
        value = c(de$x, rev(de$x))
      )
    }))
    p <- p +
      geom_polygon(
        data = den,
        aes(group = group),
        alpha = 0.25,
        linewidth = 0.45
      ) +
      geom_point(
        aes(
          x = x +
            0.15
        ),
        position = position_jitter(width = 0.09, seed = 42),
        size = 0.9,
        alpha = 0.27
      ) +
      geom_boxplot(
        aes(group = group),
        width = 0.11,
        outlier.shape = NA,
        fill = "white",
        linewidth = 0.55
      )
    by <- rr[2] + diff(rr) * 0.08
    star <- tests[variable == v, label]
    if (nzchar(star)) {
      p <- p +
        annotate(
          "text",
          x = 1.5,
          y = by,
          label = star,
          size = 7,
          family = "Arial",
          colour = "black"
        )
    }
    p +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      scale_x_continuous(
        breaks = c(1, 2),
        labels = paste0(
          c("C1", "C2"),
          "\n(n = ",
          ss$N[match(c("C1", "C2"), ss$group)],
          ")"
        ),
        limits = c(0.48, 2.52)
      ) +
      scale_y_continuous(expand = expansion(mult = 0)) +
      coord_cartesian(ylim = lims) +
      labs(x = NULL, y = paste0(v, " (", units[v], ")")) +
      theme_classic(base_size = 13, base_family = "Arial") +
      theme(
        legend.position = "none",
        axis.title.y = element_text(size = 14, margin = margin(r = 7)),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(
          size = 11,
          hjust = 0.5,
          margin = margin(b = 10)
        ),
        axis.text.x = element_text(
          size = 11.5,
          color = "#222222",
          lineheight = 1.15
        ),
        axis.text.y = element_text(size = 11, color = "#333333"),
        axis.line = element_line(linewidth = 0.5),
        axis.ticks = element_line(linewidth = 0.5),
        axis.ticks.length = unit(3, "pt"),
        plot.margin = margin(8, 12, 5, 5)
      )
  })
  pp <- lapply(pp, function(q) {
    q +
      theme(
        axis.title.y = element_text(size = 17),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13.5)
      )
  })
  p <- wrap_plots(pp, nrow = 1)

  ggsave(
    "./figures-output/Figure3F.pdf",
    p,
    width = 15.5,
    height = 3.8,
    device = cairo_pdf
  )
}



# ==============================================================================
# FIGURE 3 (PANEL G) — MACCE-free survival by UKB cluster
# ==============================================================================
# Function to plot MACCE-free survival by UKB cluster.
plot_3G <- function() {
  cv <- fread("./results/UKB/outcomes/KM_MACCE/curves.csv")
  cv[, time := time / 12]
  fit <- readRDS(
    "./results/UKB/outcomes/KM_MACCE/KM_Fits.rds"
  )[["MACCE"]]
  # Show follow-up in years; the fitted survival times are in months.
  times <- seq(0, 18, by = 3)
  # Extract the at-risk counts at the displayed follow-up times.
  rt <- rbindlist(lapply(names(fit$strata), function(stratum) {
    group <- sub("^Group=", "", stratum)
    group_fit <- fit[stratum]
    n <- vapply(
      times,
      function(year) {
        index <- which(group_fit$time >= year * 12)
        if (length(index)) {
          group_fit$n.risk[index[1L]]
        } else {
          0
        }
      },
      numeric(1)
    )
    data.table(endpoint = "MACCE", Group = group, time = times, n_risk = n)
  }))
  ct <- fread("./results/UKB/outcomes/KM_MACCE/counts.csv")
  tt <- fread("./results/UKB/outcomes/KM_MACCE/tests.csv")
  o <- list(curves = cv, risk = rt, counts = ct, tests = tt)
  pal <- c(C1 = "#F4A7AE", C2 = "#5878A6")
  # Function to draw the Kaplan–Meier curves and corresponding numbers at risk.
  draw_km <- function(e, o) {
    W <- 7.45
    H <- 4.9
    left <- 1.46
    right <- 7.2
    bottom <- 2.04
    top <- 4.36
    cv <- o$curves[endpoint == e]
    rt <- o$risk[endpoint == e]
    ct <- o$counts[endpoint == e]
    pv <- o$tests[endpoint == e, Logrank_P]
    ylo <- 0.4
    ymajor <- seq(ylo, 1, 0.1)
    yminor <- seq(ylo, 1, 0.05)
    # Function to map a value to its horizontal panel coordinate.
    xx <- function(t) left + (right - left) * t / 18
    # Function to map a value to its vertical panel coordinate.
    yy <- function(v) bottom + (top - bottom) * (v - ylo) / (1 - ylo)
    ptxt <- if (pv < 0.001) {
      expression(paste("Log-rank ", italic(P), " < 0.001"))
    } else {
      bquote(paste("Log-rank ", italic(P), " = ", .(sprintf("%.3f", pv))))
    }
    # Function to draw survival curves and the risk table.
    draw <- function() {
      grid.newpage()
      # Function to position text in panel coordinates.
      txt <- function(label, x, y, size = 12, just = "centre", rot = 0) {
        grid.text(
          label,
          x = unit(x, "in"),
          y = unit(y, "in"),
          just = just,
          rot = rot,
          gp = gpar(fontfamily = "Noto Sans", fontsize = size, col = "black")
        )
      }
      # Function to draw a line in panel coordinates.
      ln <- function(x, y, x2, y2, col = "#E8EBED", lwd = 0.7) {
        grid.lines(
          x = unit(c(x, x2), "in"),
          y = unit(c(y, y2), "in"),
          gp = gpar(col = col, lwd = lwd)
        )
      }
      # Function to draw a rectangle in panel coordinates.
      rect <- function(x, y, w, h, col = "#888888", fill = "white", lwd = 0.6) {
        grid.rect(
          x = unit(x, "in"),
          y = unit(y, "in"),
          width = unit(w, "in"),
          height = unit(h, "in"),
          just = c("left", "bottom"),
          gp = gpar(col = col, fill = fill, lwd = lwd)
        )
      }
      txt(e, (left + right) / 2, 4.64, 14)
      for (t in seq(0, 18, 1)) {
        ln(
          xx(t),
          bottom,
          xx(t),
          top,
          if (t %% 3 == 0) {
            "#E5E9EC"
          } else {
            "#F0F2F3"
          },
          if (t %% 3 == 0) {
            0.8
          } else {
            0.6
          }
        )
      }
      for (v in yminor) {
        major <- any(abs(ymajor - v) < 1e-08)
        ln(
          left,
          yy(v),
          right,
          yy(v),
          if (major) {
            "#E5E9EC"
          } else {
            "#F0F2F3"
          },
          if (major) {
            0.8
          } else {
            0.6
          }
        )
      }
      for (t in c(0, 3, 6, 9, 12, 15, 18)) {
        ln(xx(t), bottom, xx(t), bottom - 0.045, "#444444", 0.6)
        txt(as.character(t), xx(t), bottom - 0.16, 12)
      }
      for (v in ymajor) {
        ln(left - 0.04, yy(v), left, yy(v), "#444444", 0.6)
        txt(sprintf("%.1f", v), left - 0.09, yy(v), 12, "right")
      }
      txt("Event-free survival", left - 0.55, (top + bottom) / 2, 13, rot = 90)
      txt("Years", (left + right) / 2, 1.6, 13)
      # Draw stepwise 95% log-log confidence bands beneath both survival curves.
      for (g in names(pal)) {
        endg <- min(18, max(fit[paste0("Group=", g)]$time) / 12)
        d <- rbind(
          data.table(time = 0, lower = 1, upper = 1),
          cv[Group == g & time <= endg, .(time, lower, upper)]
        )
        setorder(d, time)
        last <- copy(tail(d, 1))
        last[, time := endg]
        d <- rbind(d, last)
        tx <- rep(d$time, each = 2)[-1]
        lo <- head(rep(d$lower, each = 2), -1)
        hi <- head(rep(d$upper, each = 2), -1)
        grid.polygon(
          x = unit(c(xx(tx), rev(xx(tx))), "in"),
          y = unit(c(yy(lo), rev(yy(hi))), "in"),
          gp = gpar(fill = adjustcolor(pal[g], alpha.f = 0.16), col = NA)
        )
      }
      for (g in names(pal)) {
        endg <- min(18, max(fit[paste0("Group=", g)]$time) / 12)
        d <- rbind(
          data.table(time = 0, S = 1),
          cv[Group == g & time <= endg, .(time, S)]
        )
        setorder(d, time)
        d <- rbind(d, data.table(time = endg, S = tail(d$S, 1)))
        tx <- rep(d$time, each = 2)[-1]
        sy <- rep(d$S, each = 2)
        sy <- sy[-length(sy)]
        grid.lines(
          x = unit(xx(tx), "in"),
          y = unit(yy(sy), "in"),
          gp = gpar(col = pal[g], lwd = 2.2, linejoin = "round")
        )
      }
      rect(
        left + 0.1,
        bottom + 0.1,
        0.73,
        0.53,
        col = "#777777",
        fill = "white"
      )
      for (i in seq_along(pal)) {
        gy <- bottom + 0.49 - (i - 1) * 0.23
        ln(left + 0.19, gy, left + 0.38, gy, pal[i], 2.2)
        txt(names(pal)[i], left + 0.47, gy, 12, "left")
      }
      txt(
        ptxt,
        xx(7.5),
        bottom +
          if (e == "MACCE") {
            0.28
          } else {
            0.4
          },
        11.5,
        "left"
      )
      riskleft <- left - 0.55
      riskright <- right
      txt("Number at risk", riskleft, 1.36, 12, "left")
      rect(
        riskleft,
        0.23,
        riskright - riskleft,
        1,
        col = "#888888",
        fill = "white",
        lwd = 0.65
      )
      riskx <- seq(riskleft + 2.68, riskright - 0.23, length.out = 7)
      txt("Years", riskleft + 0.11, 1.06, 10.5, "left")
      for (j in seq_along(riskx)) {
        txt(as.character(c(0, 3, 6, 9, 12, 15, 18)[j]), riskx[j], 1.06, 10.5)
      }
      for (i in seq_along(pal)) {
        g <- names(pal)[i]
        ry <- if (i == 1) {
          0.73
        } else {
          0.41
        }
        z <- ct[Group == g]
        txt(
          sprintf("%s (n = %d; events = %d)", g, z$N, z$Events),
          riskleft + 0.11,
          ry,
          10.5,
          "left"
        )
        for (j in seq_along(riskx)) {
          txt(
            as.character(rt[
              Group == g & time == c(0, 3, 6, 9, 12, 15, 18)[j],
              n_risk
            ]),
            riskx[j],
            ry,
            11
          )
        }
      }
    }

    cairo_pdf(
      "./figures-output/Figure3G.pdf",
      width = W,
      height = H,
      family = "Noto Sans"
    )
    draw()
    dev.off()
  }
  draw_km("MACCE", o)
}



# ==============================================================================
# FIGURE 3 (PANELS A–G) — Draw and save panels
# ==============================================================================

dir.create("./figures-output", recursive = TRUE, showWarnings = FALSE)
# FIGURE 3 (PANEL A)
plot_3A()
comparison <- readRDS("./results/cross-cohort/Comparison.rds")
# FIGURE 3 (PANEL B)
plot_3B(comparison)
# FIGURE 3 (PANEL C)
plot_3C(comparison)
# FIGURE 3 (PANEL D)
plot_3D(comparison)
# FIGURE 3 (PANEL E)
plot_3E()
# FIGURE 3 (PANEL F)
plot_3F()
# FIGURE 3 (PANEL G)
plot_3G()
