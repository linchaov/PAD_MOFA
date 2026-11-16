# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3, ggrepel 0.9.8, patchwork 1.3.2
#
# Inputs: Clinical traits, factor scores and Rutherford comparisons; Cox, abundance and cluster
#          survival results for VIP-PAD.
#
# Outputs: ./figures-output/Figure2A.pdf through Figure2I.pdf
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
  library(ggrepel)
})



# ==============================================================================
# FIGURE 2 (PANEL A) — Variance explained and clinical trait associations
# ==============================================================================
# Function to plot variance explained and clinical trait associations.
plot_2A <- function() {
  traits <- as.data.table(readRDS(
    "./results/VIP-PAD/clinical/Clinical_Traits_Updated.rds"
  ))
  # Display BMI immediately after sex within clinical characteristics.
  bmi_id <- traits[toupper(label) == "BMI", id]
  sex_id <- traits[grepl("^Sex", label), id]
  traits[id == bmi_id, group := traits[id == sex_id, group]]
  trait_order <- traits[id != bmi_id, id]
  trait_order <- append(trait_order, bmi_id, after = match(sex_id, trait_order))
  traits <- traits[match(trait_order, id)]
  r <- fread(
    "./results/VIP-PAD/clinical/Figure_2_Updated_Associations.csv"
  )
  fm <- fread(
    "./results/VIP-PAD/clinical/Variance_Explained.csv"
  )
  setnames(
    fm,
    c("Protein", "Metabolite"),
    c("Protein_R2_percent", "Metabolite_R2_percent")
  )
  colour_limit <- 1
  group_gap <- 0.3 / 2.54
  column_width <- 0.38
  figure_width <- 2.77 + nrow(traits) * column_width + 3 * group_gap
  figure_height <- 6.7
  # Function to lay out the clinical correlations and variance explained.
  draw_figure <- function() {
    grid.newpage()
    pushViewport(viewport(
      xscale = c(0, figure_width),
      yscale = c(0, figure_height)
    ))
    # Function to position text in panel coordinates.
    txt <- function(
      s,
      x,
      y,
      size = 10,
      just = "centre",
      rot = 0,
      bold = FALSE,
      col = "black"
    ) {
      grid.text(
        s,
        x = unit(x, "native"),
        y = unit(y, "native"),
        just = just,
        rot = rot,
        gp = gpar(
          fontfamily = "Noto Sans",
          fontsize = size * 1.12,
          col = col,
          fontface = if (bold) {
            "bold"
          } else {
            "plain"
          }
        )
      )
    }
    # Function to draw a rectangle in panel coordinates.
    rect <- function(x, y, w, h, fill, col = NA) {
      grid.rect(
        x = unit(x, "native"),
        y = unit(y, "native"),
        width = unit(w, "native"),
        height = unit(h, "native"),
        gp = gpar(fill = fill, col = col)
      )
    }
    # Function to draw a line segment in panel coordinates.
    seg <- function(x0, y0, x1, y1, col = "#DFE3E7", lwd = 0.4) {
      grid.lines(
        x = unit(c(x0, x1), "native"),
        y = unit(c(y0, y1), "native"),
        gp = gpar(col = col, lwd = lwd)
      )
    }
    # Function to map the adjusted P-value to a bubble diameter.
    bubble_diameter <- function(q) {
      3 + 5 * sqrt(pmin(6, -log10(pmax(q, 1e-300))))
    }
    ys <- 4.68 - (0:9) * 0.4
    left <- 2.25
    group_index <- match(traits$group, unique(traits$group)) - 1L
    column_left <- left +
      (seq_len(nrow(traits)) - 1) * column_width +
      group_index * group_gap
    purp <- colorRampPalette(c("#F8F7FB", "#D0C2E2", "#8468BA", "#42258C"))(301)
    pal <- colorRampPalette(c("darkblue", "white", "darkred"))(4001)
    for (k in 1:10) {
      txt(paste("Factor", k), 0.92, ys[k], 9, "right")
      for (v in 1:2) {
        a <- fm[[c("Protein_R2_percent", "Metabolite_R2_percent")[v]]][k]
        x <- 1.24 + (v - 1) * 0.57
        rect(x, ys[k], 0.55, 0.38, purp[1 + round(min(a, 30) * 10)])
        txt(
          sprintf("%.1f", a),
          x,
          ys[k],
          8,
          col = if (a > 17) {
            "white"
          } else {
            "black"
          }
        )
      }
    }
    txt("Proteomics", 1.24, 5.02, 9, "left", 55)
    txt("LC-MS metabolomics", 1.81, 5.02, 9, "left", 55)
    for (g in unique(traits$group)) {
      ii <- which(traits$group == g)
      lo <- column_left[min(ii)]
      hi <- column_left[max(ii)] + column_width
      txt(g, (lo + hi) / 2, 6.4, 10, bold = TRUE)
      seg(lo, 4.94, hi, 4.94, "#7F8A94", 0.6)
      for (boundary_index in 0:length(ii)) {
        tick_x <- lo + boundary_index * column_width
        seg(tick_x, 4.94, tick_x, 4.885, "#7F8A94", 0.5)
      }
    }
    for (j in 1:29) {
      x <- column_left[j] + column_width / 2
      txt(traits$label[j], x, 5.02, 9, "left", 55)
      for (k in 1:10) {
        a <- r[id == traits$id[j] & factor == fm$factor[k]]
        diam <- bubble_diameter(a$FDR)
        col <- pal[
          1 + round((a$estimate + colour_limit) / (2 * colour_limit) * 4000)
        ]
        grid.circle(
          x = unit(x, "native"),
          y = unit(ys[k], "native"),
          r = unit(diam / 2, "pt"),
          gp = gpar(fill = col, col = NA, lwd = 0.35)
        )
        if (a$stars != "") {
          txt(
            a$stars,
            x,
            ys[k] + 0.005,
            6.6,
            col = if (abs(a$estimate) > 0.65) {
              "white"
            } else {
              "black"
            }
          )
        }
      }
    }
    txt("Variance (%)", 1.55, 0.63, 8)
    for (i in 1:150) {
      rect(0.99 + (i - 0.5) * 0.0073, 0.41, 0.0075, 0.11, purp[round(i * 2)])
    }
    txt("0", 0.99, 0.19, 8)
    txt("15", 1.54, 0.19, 8)
    txt("30", 2.09, 0.19, 8)
    correlation_width <- 1.35
    txt("Correlation", left + correlation_width / 2, 0.63, 8)
    for (i in 1:200) {
      rect(
        left + (i - 0.5) * correlation_width / 200,
        0.41,
        correlation_width / 200 + 2e-04,
        0.11,
        pal[1 + round((i - 1) / 199 * 4000)]
      )
    }
    for (tick in c(
      -colour_limit,
      -colour_limit / 2,
      0,
      colour_limit / 2,
      colour_limit
    )) {
      xx <- left +
        (tick + colour_limit) / (2 * colour_limit) * correlation_width
      seg(xx, 0.35, xx, 0.3, "#666666", 0.4)
      txt(as.character(tick), xx, 0.19, 8)
    }
    popViewport()
  }

  cairo_pdf(
    "./figures-output/Figure2A.pdf",
    width = figure_width,
    height = figure_height,
    family = "Noto Sans"
  )
  draw_figure()
  invisible(dev.off())
}



# ==============================================================================
# FIGURE 2 (PANEL B) — Factor scores by Rutherford grade
# ==============================================================================
# Function to plot factor scores by Rutherford grade.
plot_2B <- function() {
  long <- fread(
    "./results/VIP-PAD/clinical/Rutherford_Factor_Scores.csv"
  )
  stats <- fread(
    "./results/VIP-PAD/clinical/Rutherford_Factor_Statistics.csv"
  )
  levels_factor <- sprintf("F%02d", 1:10)
  long[, factor := factor(factor, levels = levels_factor)]
  stats[, factor := factor(factor, levels = levels_factor)]
  sample_counts <- long[factor == levels_factor[1L], .N, by = severity]
  stats[, factor_num := as.integer(factor)]
  stats[,
    y := vapply(
      as.character(factor),
      function(f) max(long$score[long$factor == f]) + 0.25,
      numeric(1)
    )
  ]
  # Lower the F6 annotation to leave clearance below the upper border.
  stats[factor == "F06", y := y - 0.65]
  long[, x := as.integer(factor) + ifelse(severity == "0-3", -0.19, 0.19)]
  p <- ggplot(
    long,
    aes(
      x = x,
      y = score,
      group = interaction(factor, severity),
      fill = severity
    )
  ) +
    geom_hline(
      yintercept = 0,
      colour = "grey75",
      linetype = "dashed",
      linewidth = 0.25
    ) +
    geom_point(
      shape = 21,
      size = 0.55,
      stroke = 0.18,
      alpha = 0.25,
      position = position_jitter(width = 0.105, height = 0, seed = 42)
    ) +
    geom_boxplot(
      width = 0.25,
      outlier.shape = NA,
      linewidth = 0.35,
      alpha = 0.95
    ) +
    geom_text(
      data = stats,
      aes(x = factor_num, y = y, label = annotation),
      inherit.aes = FALSE,
      size = 6.5,
      lineheight = 1.15,
      vjust = 0
    ) +
    scale_fill_manual(
      values = c(`0-3` = "white", `4-6` = "#438CA5"),
      labels = sprintf(
        "Rutherford %s (n=%d)",
        c("0-3", "4-6"),
        sample_counts$N[match(c("0-3", "4-6"), sample_counts$severity)]
      )
    ) +
    scale_x_continuous(breaks = 1:10, labels = 1:10, limits = c(0.5, 10.5)) +
    expand_limits(y = max(long$score) + 0.25) +
    scale_y_continuous(expand = expansion(mult = c(0.025, 0.065))) +
    labs(x = "MOFA factors", y = "Factor score", fill = NULL) +
    theme_classic(base_size = 15, base_family = "Arial") +
    theme(
      legend.position = "top",
      legend.justification = "left",
      legend.key.size = grid::unit(3, "mm"),
      legend.text = element_text(size = 13),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = grid::unit(1, "mm"),
      plot.title = element_text(size = 11, face = "bold"),
      plot.subtitle = element_text(size = 9),
      plot.caption = element_text(size = 7, hjust = 0),
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.65
      ),
      axis.line = element_blank(),
      axis.ticks = element_line(linewidth = 0.55),
      axis.text = element_text(colour = "black", size = 13),
      axis.title = element_text(size = 16),
      plot.margin = margin(8, 8, 8, 8)
    )

  ggsave(
    "./figures-output/Figure2B.pdf",
    p,
    width = 9,
    height = 3.5,
    device = cairo_pdf
  )
}



# ==============================================================================
# FIGURE 2 (PANEL C) — Age- and sex-adjusted factor associations with MALE and MACCE
# ==============================================================================
# Function to plot age- and sex-adjusted factor associations with MALE and MACCE.
plot_2C <- function() {
  estimates <- read.csv(
    "./results/VIP-PAD/outcomes/Cox/Associations.csv",
    stringsAsFactors = FALSE
  )
  endpoints <- c("MALE", "MACCE")
  factors <- sprintf("F%02d", 1:10)
  ticks <- list(
    MALE = c(0.5, 0.75, 1, 1.25, 1.5),
    MACCE = c(0.5, 0.75, 1, 1.5, 2, 2.5)
  )
  colours <- c(
    Lower = "#426F9B",
    Higher = "#B75C59",
    Nonsignificant = "#626A73"
  )
  fills <- colours
  fills["Nonsignificant"] <- "white"
  panels <- list()
  for (endpoint in endpoints) {
    d <- estimates[estimates$endpoint == endpoint, ]
    d$row <- 11L - match(d$factor, factors)
    d$status <- ifelse(
      d$FDR < 0.05,
      ifelse(d$HR < 1, "Lower", "Higher"),
      "Nonsignificant"
    )
    limits <- range(ticks[[endpoint]])
    p <- ggplot(d, aes(x = HR, y = row, colour = status)) +
      geom_vline(
        xintercept = 1,
        linetype = "dashed",
        colour = "#ADB4BD",
        linewidth = 0.35
      ) +
      geom_segment(
        aes(x = CI_low, xend = CI_high, yend = row),
        linewidth = 0.45
      ) +
      geom_segment(
        aes(x = CI_low, xend = CI_low, y = row - 0.11, yend = row + 0.11),
        linewidth = 0.45
      ) +
      geom_segment(
        aes(
          x = CI_high,
          xend = CI_high,
          y = row - 0.11,
          yend = row +
            0.11
        ),
        linewidth = 0.45
      ) +
      geom_point(aes(fill = status), shape = 21, size = 2.5, stroke = 0.55) +
      scale_colour_manual(values = colours) +
      scale_fill_manual(values = fills) +
      scale_x_log10(
        limits = limits,
        breaks = ticks[[endpoint]],
        labels = as.character(ticks[[endpoint]]),
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        breaks = 1:10,
        labels = paste("Factor", 10:1),
        limits = c(0.4, 10.6),
        expand = c(0, 0)
      ) +
      labs(title = endpoint, x = "Adjusted HR (95% CI)", y = NULL) +
      theme_classic(base_family = "Arial", base_size = 10) +
      theme(
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "none",
        plot.margin = margin(8, 6, 8, 6)
      )
    if (endpoint == "MACCE") {
      p <- p + theme(axis.text.y = element_blank())
    }
    panels[[endpoint]] <- p
  }

  ggsave(
    "./figures-output/Figure2C.pdf",
    panels$MALE | panels$MACCE,
    width = 5.4,
    height = 4.7,
    units = "in",
    device = grDevices::cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 2 (PANEL D) — Numbers of differentially abundant proteins and metabolites
# ==============================================================================
# Function to plot numbers of differentially abundant proteins and metabolites.
plot_2D <- function() {
  s <- fread(
    "./results/VIP-PAD/abundance/Contrast_Summary.csv"
  )
  d <- melt(
    s,
    id.vars = "modality",
    measure.vars = c("higher_C2", "higher_C1"),
    variable.name = "direction",
    value.name = "count"
  )
  d[, signed_count := ifelse(direction == "higher_C2", -count, count)]
  d[,
    modality := factor(
      modality,
      levels = c("Metabolite", "Protein"),
      labels = c("Metabolites", "Proteins")
    )
  ]
  d[, label_x := signed_count + ifelse(signed_count < 0, -45, 45)]
  d[, label_hjust := ifelse(signed_count < 0, 1, 0)]
  # Shared linear scale allows direct comparison across assays; all labels are counts.
  p <- ggplot(d, aes(x = signed_count, y = modality, fill = direction)) +
    geom_vline(xintercept = 0, colour = "#555555", linewidth = 0.4) +
    geom_col(width = 0.5, position = "identity") +
    geom_text(
      aes(
        x = label_x,
        label = format(count, big.mark = ",", trim = TRUE),
        hjust = label_hjust
      ),
      family = "Noto Sans",
      size = 4.8,
      colour = "black"
    ) +
    scale_fill_manual(
      values = c(higher_C2 = "#5878A6", higher_C1 = "#F4A7AE"),
      breaks = c("higher_C2", "higher_C1"),
      labels = c("C2 higher", "C1 higher"),
      name = NULL
    ) +
    scale_x_continuous(
      limits = c(-1850, 1850),
      breaks = seq(-1500, 1500, 500),
      labels = function(x) format(abs(x), big.mark = ",", trim = TRUE),
      expand = c(0, 0)
    ) +
    scale_y_discrete(expand = expansion(add = 0.55)) +
    labs(x = "Number of differential molecules", y = NULL) +
    theme_classic(base_family = "Noto Sans", base_size = 14) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(size = 14, margin = margin(r = 10)),
      axis.text.x = element_text(size = 11),
      axis.title.x = element_text(size = 13, margin = margin(t = 8)),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "right",
      legend.box.spacing = grid::unit(10, "mm"),
      legend.text = element_text(size = 12),
      legend.key.size = grid::unit(4, "mm"),
      legend.key.height = grid::unit(8, "mm"),
      legend.key.spacing.y = grid::unit(3, "mm"),
      plot.margin = margin(6, 10, 7, 8)
    )

  ggsave(
    "./figures-output/Figure2D.pdf",
    p,
    width = 7.9,
    height = 2.5,
    device = cairo_pdf,
    family = "Noto Sans"
  )
}



# ==============================================================================
# FIGURE 2 (PANEL E) — Clinical distributions in the molecular clusters
# ==============================================================================
# Function to plot clinical distributions in the molecular clusters.
plot_2E <- function() {
  d <- readRDS(
    "./results/VIP-PAD/clinical/Distributions.rds"
  )
  st <- fread(
    "./results/VIP-PAD/clinical/Descriptive_Statistics.csv"
  )
  st[, group := factor(group, levels = c("C1", "C2"))]
  variables <- c("Age", "BMI", "SBP", "LDL", "eGFR")
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
    "./figures-output/Figure2E.pdf",
    p,
    width = 15.5,
    height = 3.8,
    device = cairo_pdf
  )
}



# ==============================================================================
# FIGURE 2 (PANELS F–G) — Protein and metabolite abundance contrasts between clusters
# ==============================================================================
# Function to plot protein and metabolite abundance contrasts between clusters.
plot_2FG <- function() {
  r <- fread(
    "./results/VIP-PAD/abundance/All_Molecular_Contrasts.csv"
  )
  counts <- fread(
    "./results/VIP-PAD/abundance/Contrast_Summary.csv"
  )
  # Function to wrap long labels to fit the panel.
  wrap_label <- function(s, width = 34) {
    s <- sub("\\s*\\[L[0-9]+\\]$", "", s)
    out <- character()
    while (nchar(s) > width) {
      pos <- gregexpr("[ -]", substr(s, 1, width))[[1]]
      cut <- if (max(pos) > width * 0.5) {
        max(pos)
      } else {
        width
      }
      out <- c(out, substr(s, 1, cut))
      s <- trimws(substr(s, cut + 1, nchar(s)))
    }
    paste(c(out, s), collapse = "\n")
  }
  for (view in c("Protein", "Metabolite")) {
    z <- copy(r[modality == view])
    nlab <- if (view == "Protein") {
      6L
    } else {
      4L
    }
    # Label leading features from each direction of the abundance contrast.
    labels <- rbindlist(lapply(c(-1, 1), function(s) {
      q <- z[FDR < 0.05 & sign(log2FC_C1_vs_C2) == s][order(
        FDR,
        -abs(log2FC_C1_vs_C2),
        feature
      )]
      q <- q[!duplicated(label)]
      head(q, nlab)
    }))
    labels[,
      display_label := vapply(
        label,
        wrap_label,
        character(1),
        width = if (view == "Protein") {
          34
        } else {
          29
        }
      )
    ]
    count <- counts[modality == view]
    legend_labels <- c(
      `Higher in C2` = sprintf("C2 higher (%d)", count$higher_C2),
      `Not significant` = sprintf("FDR >= 0.05 (%d)", count$not_significant),
      `Higher in C1` = sprintf("C1 higher (%d)", count$higher_C1)
    )
    xhalf <- max(
      if (view == "Protein") 4.2 else 14.8,
      ceiling(max(abs(z$log2FC_C1_vs_C2)) * 1.08)
    )
    xlim <- c(-xhalf, xhalf)
    ymax <- max(z$neglog10FDR)
    height <- 6.2
    ytop <- max(
      if (view == "Protein") 250 else 150,
      ceiling(ymax * 1.12 / 50) * 50
    )
    ybottom <- -ytop * 0.055
    p <- ggplot(z, aes(log2FC_C1_vs_C2, neglog10FDR)) +
      geom_hline(
        yintercept = -log10(0.05),
        linetype = "dashed",
        colour = "#777777",
        linewidth = 0.45
      ) +
      geom_vline(
        xintercept = 0,
        linetype = "dashed",
        colour = "#AAAAAA",
        linewidth = 0.4
      ) +
      geom_point(
        data = z[category == "Not significant"],
        aes(colour = category),
        size = 1.45,
        alpha = 0.55
      ) +
      geom_point(
        data = z[category != "Not significant"],
        aes(colour = category),
        size = 1.55,
        alpha = 0.72
      ) +
      geom_text_repel(
        data = labels,
        aes(label = display_label),
        family = "Noto Sans",
        colour = "black",
        size = 4.8,
        box.padding = 0.5,
        point.padding = 0.18,
        min.segment.length = 0,
        segment.colour = "#666666",
        segment.size = 0.3,
        nudge_y = if (view == "Protein") {
          6
        } else {
          4
        },
        nudge_x = ifelse(labels$log2FC_C1_vs_C2 > 0, 1, -1) *
          if (view == "Protein") {
            0.16
          } else {
            0.55
          },
        seed = 42,
        max.overlaps = Inf,
        max.iter = 1e+05,
        max.time = 10,
        force = 3,
        lineheight = 0.94
      ) +
      scale_colour_manual(
        values = c(
          `Higher in C2` = "#5878A6",
          `Not significant` = "#BEBEBE",
          `Higher in C1` = "#F4A7AE"
        ),
        breaks = names(legend_labels),
        labels = legend_labels,
        drop = FALSE,
        name = NULL
      ) +
      scale_x_continuous(
        limits = xlim,
        breaks = if (view == "Protein") {
          seq(-4, 4, 2)
        } else {
          seq(-10, 10, 5)
        },
        expand = c(0, 0)
      ) +
      scale_y_continuous(
        limits = c(ybottom, ytop),
        breaks = seq(0, ytop, 50),
        expand = c(0, 0)
      ) +
      labs(
        x = expression("Cluster 1 vs. Cluster 2 (log"[2] * " FC)"),
        y = expression(-log[10](FDR))
      ) +
      theme_classic(base_family = "Noto Sans", base_size = 14) +
      theme(
        text = element_text(colour = "black"),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(size = 17),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        legend.position = "none",
        legend.text = element_text(size = 10),
        legend.key.width = grid::unit(5, "mm"),
        legend.spacing.x = grid::unit(2, "mm"),
        plot.margin = margin(12, 14, 10, 10),
        panel.border = element_rect(
          fill = NA,
          colour = "#555555",
          linewidth = 0.45
        ),
        axis.line = element_blank()
      ) +
      guides(
        colour = guide_legend(
          nrow = 1,
          override.aes = list(alpha = 1, size = 2.5)
        )
      )

    ggsave(
      c(
        Protein = "./figures-output/Figure2F.pdf",
        Metabolite = "./figures-output/Figure2G.pdf"
      )[[view]],
      p,
      width = 7.7,
      height = height,
      device = cairo_pdf,
      family = "Noto Sans",
      limitsize = FALSE
    )
  }
}



# ==============================================================================
# FIGURE 2 (PANELS H–I) — MALE- and MACCE-free survival by molecular cluster
# ==============================================================================
# Function to plot MALE- and MACCE-free survival by molecular cluster.
plot_2HI <- function() {
  pal <- c(C1 = "#F4A7AE", C2 = "#5878A6")
  # Function to draw the Kaplan–Meier curves and corresponding numbers at risk.
  draw_km <- function(e, o) {
    W <- 6
    H <- 4.9
    left <- 1.46
    right <- 5.76
    bottom <- 2.04
    top <- 4.36
    cv <- o$curves[endpoint == e]
    rt <- o$risk[endpoint == e]
    ct <- o$counts[endpoint == e]
    pv <- o$tests[endpoint == e, Logrank_P]
    ylo <- if (e == "MALE") {
      0
    } else {
      min(0.7, floor(min(cv[time <= 36, lower], 1, na.rm = TRUE) * 10) / 10)
    }
    ymajor <- seq(
      ylo,
      1,
      if (e == "MALE") {
        0.2
      } else {
        0.1
      }
    )
    yminor <- seq(
      ylo,
      1,
      if (e == "MALE") {
        0.1
      } else {
        0.05
      }
    )
    # Function to map a value to its horizontal panel coordinate.
    xx <- function(t) left + (right - left) * t / 36
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
      for (t in seq(0, 36, 3)) {
        ln(
          xx(t),
          bottom,
          xx(t),
          top,
          if (t %% 12 == 0) {
            "#E5E9EC"
          } else {
            "#F0F2F3"
          },
          if (t %% 12 == 0) {
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
      for (t in c(0, 12, 24, 36)) {
        ln(xx(t), bottom, xx(t), bottom - 0.045, "#444444", 0.6)
        txt(as.character(t), xx(t), bottom - 0.16, 12)
      }
      for (v in ymajor) {
        ln(left - 0.04, yy(v), left, yy(v), "#444444", 0.6)
        txt(sprintf("%.1f", v), left - 0.09, yy(v), 12, "right")
      }
      txt("Event-free survival", left - 0.55, (top + bottom) / 2, 13, rot = 90)
      txt("Months", (left + right) / 2, 1.6, 13)
      # Draw stepwise 95% log-log confidence bands beneath both survival curves.
      for (g in names(pal)) {
        d <- rbind(
          data.table(time = 0, lower = 1, upper = 1),
          cv[Group == g & time <= 36, .(time, lower, upper)]
        )
        setorder(d, time)
        last <- copy(tail(d, 1))
        last[, time := 36]
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
        d <- rbind(
          data.table(time = 0, S = 1),
          cv[Group == g & time <= 36, .(time, S)]
        )
        setorder(d, time)
        d <- rbind(d, data.table(time = 36, S = tail(d$S, 1)))
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
        xx(18),
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
      riskx <- seq(riskleft + 2.85, riskright - 0.25, length.out = 4)
      txt("Months", riskleft + 0.11, 1.06, 10.5, "left")
      for (j in seq_along(riskx)) {
        txt(as.character(c(0, 12, 24, 36)[j]), riskx[j], 1.06, 10.5)
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
            as.character(rt[Group == g & time == c(0, 12, 24, 36)[j], n_risk]),
            riskx[j],
            ry,
            11
          )
        }
      }
    }

    cairo_pdf(
      c(
        MALE = "./figures-output/Figure2H.pdf",
        MACCE = "./figures-output/Figure2I.pdf"
      )[[e]],
      width = W,
      height = H,
      family = "Noto Sans"
    )
    draw()
    dev.off()
  }
  o <- list(
    curves = fread("./results/VIP-PAD/outcomes/KM/curves.csv"),
    risk = fread("./results/VIP-PAD/outcomes/KM/risk.csv"),
    counts = fread("./results/VIP-PAD/outcomes/KM/counts.csv"),
    tests = fread("./results/VIP-PAD/outcomes/KM/tests.csv")
  )
  for (e in c("MALE", "MACCE")) {
    draw_km(e, o)
  }
}



# ==============================================================================
# FIGURE 2 (PANELS A–I) — Draw and save panels
# ==============================================================================

dir.create("./figures-output", recursive = TRUE, showWarnings = FALSE)
# FIGURE 2 (PANEL A)
plot_2A()
# FIGURE 2 (PANEL B)
plot_2B()
# FIGURE 2 (PANEL C)
plot_2C()
# FIGURE 2 (PANEL D)
plot_2D()
# FIGURE 2 (PANEL E)
plot_2E()
# FIGURE 2 (PANELS F–G)
plot_2FG()
# FIGURE 2 (PANELS H–I)
plot_2HI()
