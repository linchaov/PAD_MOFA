# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggbeeswarm 0.7.3, ggplot2 4.0.3, ggrepel 0.9.8,
#                  scales 1.4.0
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds heatmap measurements, VIP-PAD Factor 2 loadings, feature
#          annotations, cluster assignments, contrasts, enrichment, survival and SHAP results.
#
# Outputs: ./figures-output/Figure5A.pdf through Figure5I.pdf
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
  library(ggrepel)
  library(ggbeeswarm)
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to wrap labels to the available width.
wrap_label <- function(z, n = 26) {
  vapply(
    z,
    function(s) {
      s <- gsub("(?<=[-/])", " ", s, perl = TRUE)
      paste(strwrap(s, width = n), collapse = "\n")
    },
    character(1)
  )
}

# Function to draw a rectangle in panel coordinates.
panel_rect <- function(x, y, w, hh, fill, col = NA) {
  grid.rect(
    x = unit(x, "in"),
    y = unit(y, "in"),
    width = unit(w, "in"),
    height = unit(hh, "in"),
    just = c("left", "top"),
    gp = gpar(fill = fill, col = col, lwd = 0.6)
  )
}

# Function to extract gene symbols from VIP-PAD protein identifiers.
vip_symbols <- function(features) {
  genes <- sub("_[^_]+$", "", sub("^PROT_", "", features))
  genes[genes == "Unmapped"] <- ""
  chartr("_", "-", genes)
}
# Function to join signed factor loadings to molecule labels.
factor_loadings <- function() {
  # Join signed loadings to protein and metabolite labels.
  loads <- rbindlist(lapply(c("Protein", "Metabolite"), function(v) {
    w <- fit$W[[v]]
    ids <- rownames(w)
    label <- if (v == "Protein") {
      vip_symbols(ids)
    } else {
      metabolites$label[match(ids, metabolites$feature)]
    }
    label[is.na(label) | label == ""] <- rownames(w)[is.na(label) | label == ""]
    data.table(
      feature = rownames(w),
      modality = v,
      label = label,
      loading_raw = w[, "Factor2"]
    )
  }))
  loads[,
    relative_weight := abs(loading_raw) / max(abs(loading_raw)),
    by = modality
  ]
}

# Function to rank metabolites by absolute loading (top 15) and retain proteins with relative weight > 0.8.
leading_molecules <- function(loads) {
  sel <- rbindlist(list(
    loads[modality == "Protein" & relative_weight > 0.8][order(
      -abs(loading_raw),
      feature
    )],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:15]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  sel
}



# ==============================================================================
# FIGURE 5 (PANEL A) — Leading proteins and metabolites
# ==============================================================================
# Function to plot leading proteins and metabolites.
plot_5A <- function() {
  sel <- leading_molecules(loads)
  font <- "Noto Sans"
  # Function to position text in panel coordinates.
  txt <- function(
    s,
    x,
    y,
    size = 12,
    just = "centre",
    rot = 0,
    col = "black",
    bold = FALSE
  ) {
    grid.text(
      s,
      x = unit(x, "in"),
      y = unit(y, "in"),
      just = just,
      rot = rot,
      gp = gpar(
        fontfamily = font,
        fontsize = size,
        col = col,
        fontface = if (bold) {
          "bold"
        } else {
          "plain"
        }
      )
    )
  }

  ca <- copy(sel)
  ca[, feature_key := factor(feature, levels = feature)]
  ca[, label_plot := ifelse(modality == "Metabolite", wrap_label(label, 18), label)]
  ca[
    grepl("Hydroxyphenylacetic", label, ignore.case = TRUE),
    label_plot := "ortho-\nHydroxyphenyl\nacetic acid"
  ]
  ca[
    grepl("Tetramethylchromano", label, ignore.case = TRUE),
    label_plot := "Tetramethyl\nchromanol\nglucoside"
  ]
  ca[, value2 := pmax(0, relative_weight - 0.45 / .N)]
  ca[,
    number_y := relative_weight - 0.05 - ifelse(seq_len(.N) %% 2 == 0, 0.055, 0)
  ]
  pa <- ggplot(ca, aes(feature_key, relative_weight, fill = modality)) +
    geom_hline(
      yintercept = c(0.25, 0.5, 0.75, 1),
      colour = "black",
      linewidth = 0.45
    ) +
    geom_col(width = 1, alpha = 0.6) +
    geom_col(aes(y = value2), width = 0.85) +
    geom_text(
      aes(y = number_y, label = sprintf("%.2f", relative_weight)),
      size = 5.5,
      family = font
    ) +
    geom_text(
      aes(y = 0.5, label = ifelse(direction == "Positive", "+", "−")),
      size = 8,
      family = font
    ) +
    scale_fill_manual(
      values = c(Protein = "#4575b4", Metabolite = "#f46d43"),
      labels = c("Metabolites", "Proteins"),
      name = NULL
    ) +
    scale_x_discrete(labels = setNames(ca$label_plot, ca$feature)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    coord_radial(inner.radius = 0.025, clip = "off") +
    guides(theta = guide_axis_theta(angle = 90)) +
    theme_void(base_family = font) +
    theme(
      axis.text.theta = element_text(size = 23, colour = "black"),
      legend.position = "bottom",
      legend.text = element_text(size = 13),
      plot.margin = margin(100, 110, 110, 170)
    ) +
    labs(caption = "Relative absolute weight") +
    theme(plot.caption = element_text(hjust = 0.5, size = 12))
  pa <- pa + theme(legend.position = "none", plot.caption = element_blank())

  # Function to draw the ranked molecular loadings in radial coordinates.
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(
      x = 0,
      y = unit(0, "in"),
      width = unit(14.5, "in"),
      height = unit(14.5, "in"),
      just = c("left", "bottom")
    ))
    print(pa, newpage = FALSE)
    txt("Factor 2 top features", 0.48, 7.8, 36, rot = 90)
    panel_rect(0.65, 2.85, 0.32, 0.32, "#4575b4")
    txt("Proteins", 1.1, 2.69, 28, "left")
    panel_rect(0.65, 2.28, 0.32, 0.32, "#f46d43")
    txt("Metabolites", 1.1, 2.12, 28, "left")
    popViewport()
  }

  cairo_pdf(
    "./figures-output/Figure5A.pdf",
    width = 14.5,
    height = 14.5,
    family = font
  )
  draw()
  dev.off()
}



# ==============================================================================
# FIGURE 5 (PANEL B) — Molecular abundance across the two clusters
# ==============================================================================
# Function to plot molecular abundance across the two clusters.
plot_5B <- function() {
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:4],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:8]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  cl <- readRDS(
    "./results/VIP-PAD/clustering/Clustering_Stability.rds"
  )
  ids <- rownames(fit$Z)
  ix <- match(ids, cohort$study_id)
  cx <- match(ids, cl$assignments$sample_id)
  groups <- cl$assignments$K2[cx]
  abundance <- as.matrix(cohort[ix, sel$feature, drop = FALSE])
  rownames(abundance) <- ids
  h <- list(
    groups = groups,
    patient_order = list(F02 = order(groups, fit$Z[, "Factor2"], ids)),
    molecular_matrix = abundance
  )
  font <- "Noto Sans"
  # Function to position text in panel coordinates.
  txt <- function(
    s,
    x,
    y,
    size = 12,
    just = "centre",
    rot = 0,
    col = "black",
    bold = FALSE
  ) {
    grid.text(
      s,
      x = unit(x, "in"),
      y = unit(y, "in"),
      just = just,
      rot = rot,
      gp = gpar(
        fontfamily = font,
        fontsize = size,
        lineheight = 0.85,
        col = col,
        fontface = if (bold) {
          "bold"
        } else {
          "plain"
        }
      )
    )
  }
  # Function to draw a line in panel coordinates.
  ln <- function(x, y, x2, y2, col = "#C7C7C7", width = 0.6, lty = 1) {
    grid.lines(
      x = unit(c(x, x2), "in"),
      y = unit(c(y, y2), "in"),
      gp = gpar(col = col, lwd = width, lty = lty)
    )
  }

  # Function to draw the cluster abundance heatmap.
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(
      x = unit(-0.65, "in"),
      y = 0,
      width = unit(15.7, "in"),
      height = unit(8.7, "in"),
      just = c("left", "bottom")
    ))
    ord <- h$patient_order$F02
    left <- 5.8
    hw <- 9.5
    g <- h$groups[ord]
    boundary <- sum(g == 1) / length(g) * hw
    rh <- 0.26
    pal <- colorRampPalette(c(
      "#4575B4",
      "#ABD9E9",
      "#FFFFBF",
      "#FDAE61",
      "#D73027"
    ))(257)
    # Function to map standardized abundance to heatmap colors.
    cols <- function(z) {
      idx <- floor((pmin(3, pmax(-3, z)) + 3) / 6 * 256) + 1
      {
        ans <- matrix(pal[idx], nrow = nrow(z), ncol = ncol(z))
        ans[is.na(ans)] <- "#888888"
        ans
      }
    }
    # Function to render the abundance matrix as a heatmap.
    raster <- function(mat, x, top, w, hh) {
      grid.raster(
        as.raster(mat),
        x = unit(x, "in"),
        y = unit(top, "in"),
        width = unit(w, "in"),
        height = unit(hh, "in"),
        just = c("left", "top"),
        interpolate = FALSE
      )
    }
    # Brackets mark the two molecular clusters.
    for (k in 1:2) {
      x0 <- if (k == 1) {
        left
      } else {
        left + boundary
      }
      x1 <- if (k == 1) {
        left + boundary
      } else {
        left + hw
      }
      gap <- 0.035
      ln(x0 + gap, 8.09, x1 - gap, 8.09, "black", 1)
      ln(x0 + gap, 8.09, x0 + gap, 7.97, "black", 1)
      ln(x1 - gap, 8.09, x1 - gap, 7.97, "black", 1)
      txt(paste("Cluster", k), (x0 + x1) / 2, 8.32, 15)
    }
    for (v in c("Metabolite", "Protein")) {
      q <- sel[modality == v]
      nr <- nrow(q)
      top <- if (v == "Metabolite") {
        7.8
      } else {
        3.2
      }
      rh <- 0.48
      wx <- 4.12
      ww <- 1.28
      txt(
        if (v == "Protein") {
          "Proteins"
        } else {
          "Metabolites"
        },
        wx + ww / 2,
        top + 0.29,
        15,
        bold = TRUE
      )
      z <- t(scale(h$molecular_matrix[, q$feature, drop = FALSE]))[,
        ord,
        drop = FALSE
      ]
      raster(cols(z), left, top, hw, nr * rh)
      ln(left + boundary, top, left + boundary, top - nr * rh, "white", 1)
      for (vv in c(0, 0.3, 0.6, 0.9)) {
        xx <- wx + ww * vv / 1.1
        ln(xx, top, xx, top - nr * rh, "#C3C3C3", 0.8)
        txt(sprintf("%.1f", vv), xx, top - nr * rh - 0.18, 11)
      }
      for (j in 0:nr) {
        ln(wx, top - j * rh, wx + ww, top - j * rh, "#DDDDDD", 0.6)
      }
      panel_rect(wx, top, ww, nr * rh, NA, "#444444")
      panel_rect(left, top, hw, nr * rh, NA, "#999999")
      for (i in seq_len(nr)) {
        yy <- top - (i - 0.5) * rh
        txt(
          wrap_label(
            q$label[i],
            if (v == "Protein") {
              37
            } else {
              24
            }
          ),
          3.94,
          yy,
          if (v == "Protein") {
            12.5
          } else {
            11.2
          },
          "right"
        )
        xe <- wx + ww * q$relative_weight[i] / 1.1
        ln(wx, yy, xe, yy, "black", 1.8)
        grid.points(
          x = unit(xe, "in"),
          y = unit(yy, "in"),
          pch = 16,
          size = unit(0.1, "in"),
          gp = gpar(col = "black")
        )
        txt(
          if (q$direction[i] == "Positive") {
            "+"
          } else {
            "−"
          },
          5.59,
          yy,
          22,
          col = "black"
        )
      }
    }
    raster(matrix(pal, nrow = 1), 6.25, 0.44, 2, 0.12)
    txt("Z score", 7.25, 0.57, 11)
    for (j in 0:6) {
      txt(as.character(j - 3), 6.25 + j / 3, 0.16, 10)
    }
    txt("+  Positive weight     −  Negative weight", 11, 0.38, 12)

    popViewport()
  }

  cairo_pdf(
    "./figures-output/Figure5B.pdf",
    width = 15.05,
    height = 8.7,
    family = font
  )
  draw()
  dev.off()
}



# ==============================================================================
# FIGURE 5 (PANEL C) — Lipid classes ranked by Factor 2 loading
# ==============================================================================
# Function to plot Factor 2 metabolite loadings by lipid class.
plot_5C <- function(frozen_lipid) {
  class_results <- copy(frozen_lipid$lipid_results)
  lipid_molecules <- copy(frozen_lipid$lipid_molecules)
  class_order <- c(
    "AcCar",
    "LPC",
    "LPC-P",
    "LPE",
    "PC",
    "MG",
    "MG-O",
    "DG",
    "Other"
  )
  lipid_molecules[, plot_class := factor(plot_class, levels = class_order)]
  stars <- class_results[, .(plot_class = lipid_class, FDR)]
  stars[,
    label := fifelse(
      FDR < 0.001,
      "***",
      fifelse(FDR < 0.01, "**", fifelse(FDR < 0.05, "*", ""))
    )
  ]
  stars[, y := max(lipid_molecules$loading_raw) + 0.06]
  palette <- c(
    AcCar = "#5F91B4",
    LPC = "#D69A51",
    `LPC-P` = "#B8A04C",
    LPE = "#72A587",
    PC = "#8BAE78",
    MG = "#AC8BB5",
    `MG-O` = "#858DB1",
    DG = "#D4777F",
    Other = "#B7B7B7"
  )
  lipid_plot <- ggplot(
    lipid_molecules,
    aes(plot_class, loading_raw, colour = plot_class)
  ) +
    geom_hline(
      yintercept = 0,
      linetype = 2,
      colour = "grey45",
      linewidth = 0.35
    ) +
    ggbeeswarm::geom_quasirandom(width = 0.23, size = 1.5, alpha = 0.8) +
    geom_text(
      data = stars,
      aes(x = plot_class, y = y, label = label),
      inherit.aes = FALSE,
      family = "Noto Sans",
      size = 6.4,
      colour = "black"
    ) +
    scale_colour_manual(values = palette) +
    scale_x_discrete(labels = class_order) +
    scale_y_continuous(
      limits = c(
        min(lipid_molecules$loading_raw) - 0.04,
        max(lipid_molecules$loading_raw) +
          0.13
      ),
      expand = expansion(mult = 0)
    ) +
    labs(x = NULL, y = "Factor 2 weight") +
    theme_classic(base_family = "Noto Sans", base_size = 15) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 16),
      axis.title = element_text(size = 18),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(8, 10, 5, 8)
    )

  ggsave(
    "./figures-output/Figure5C.pdf",
    lipid_plot,
    device = cairo_pdf,
    width = 180,
    height = 138,
    units = "mm"
  )
}



# ==============================================================================
# FIGURE 5 (PANEL D) — Distribution of signed molecular loadings
# ==============================================================================
# Function to plot distribution of signed molecular loadings.
plot_5D <- function(x) {
  d <- rbindlist(list(
    copy(x$protein)[, .(
      feature,
      label,
      weight = loading_raw,
      view = "Proteomics"
    )],
    copy(x$metabolite)[,
      .(feature, label, weight = loading_raw, view = "Metabolomics")
    ]
  ))
  d[, scaled_weight := weight / max(abs(weight))]
  d[, view := factor(view, levels = c("Proteomics", "Metabolomics"))]
  d[, direction := ifelse(weight > 0, "Positive", "Negative")]
  d[, selected := FALSE]
  ids <- d[
    order(-abs(weight), feature),
    head(.SD, 5),
    by = .(view, direction)
  ]$feature
  d[feature %in% ids, selected := TRUE]
  set.seed(52)
  d[, position := runif(.N, -0.18, 0.18)]
  # Function to wrap long molecule names for display.
  wrap_name <- function(z, width = 18L) {
    lines <- character()
    while (nchar(z) > width) {
      candidates <- gregexpr("[ -]", substr(z, 1, width))[[1]]
      candidates <- candidates[candidates >= floor(width * 0.55)]
      cut <- if (length(candidates)) {
        max(candidates)
      } else {
        width
      }
      lines <- c(lines, trimws(substr(z, 1, cut)))
      z <- trimws(substring(z, cut + 1L))
    }
    paste(c(lines, z), collapse = "\n")
  }
  d[, display_label := vapply(label, wrap_name, character(1))]
  trans <- scales::trans_new(
    "signed_power_1.8",
    function(z) sign(z) * abs(z)^1.8,
    function(z) {
      sign(z) *
        abs(z)^(1 / 1.8)
    }
  )
  # Place labels in nonoverlapping vertical slots; leaders retain the true weights.
  label_data <- copy(d[selected == TRUE])
  label_data[, label_lines := lengths(strsplit(display_label, "\n"))]
  setorder(label_data, view, -scaled_weight)
  label_data[,
    label_y := {
      slot <- label_lines + 0.75
      top <- if (direction[1] == "Positive") {
        0.98
      } else {
        -0.02
      }
      bottom <- if (direction[1] == "Positive") {
        0.02
      } else {
        -0.98
      }
      target <- top - (cumsum(slot) - slot / 2) / sum(slot) * (top - bottom)
      sign(target) * abs(target)^(1 / 1.8)
    },
    by = .(view, direction)
  ]
  label_data[, label_x := 0]
  p <- ggplot(d, aes(position, scaled_weight)) +
    geom_hline(
      yintercept = 0,
      colour = "grey65",
      linetype = 2,
      linewidth = 0.35
    ) +
    geom_point(
      data = d[selected == FALSE],
      colour = "#929292",
      alpha = 0.23,
      size = 1.6
    ) +
    geom_point(
      data = d[selected == TRUE],
      aes(colour = direction),
      size = 2.6
    ) +
    geom_segment(
      data = label_data,
      aes(x = position, y = scaled_weight, xend = label_x, yend = label_y),
      colour = "grey60",
      linewidth = 0.25
    ) +
    geom_text(
      data = label_data[view == "Proteomics"],
      aes(x = label_x, y = label_y, label = display_label, colour = direction),
      family = "Noto Sans",
      size = 4.8,
      lineheight = 0.9
    ) +
    geom_text(
      data = label_data[view == "Metabolomics"],
      aes(x = label_x, y = label_y, label = display_label, colour = direction),
      family = "Noto Sans",
      size = 3.6,
      lineheight = 0.9
    ) +
    facet_wrap(~view, nrow = 1) +
    scale_x_continuous(limits = c(-0.72, 0.72), breaks = NULL) +
    scale_y_continuous(
      transform = trans,
      limits = c(-1.04, 1.04),
      breaks = seq(-1, 1, 0.2),
      labels = function(z) format(round(z, 1), nsmall = 1)
    ) +
    scale_colour_manual(
      values = c(Positive = "#C85B63", Negative = "#526EA4")
    ) +
    labs(
      title = NULL,
      x = "View",
      y = "Scaled weight",
      colour = "Weight direction"
    ) +
    theme_classic(base_family = "Noto Sans", base_size = 17) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 17),
      plot.title = element_text(size = 17, hjust = 0.5),
      panel.spacing = grid::unit(5, "mm"),
      legend.position = "none",
      plot.margin = margin(8, 12, 5, 8)
    )

  cairo_pdf(
    "./figures-output/Figure5D.pdf",
    width = 4.15,
    height = 6.6,
    family = "Noto Sans"
  )
  panel_plot <- ggplotGrob(p)
  panel_columns <- unique(panel_plot$layout$l[grepl(
    "^panel",
    panel_plot$layout$name
  )])
  panel_plot$widths[panel_columns] <- grid::unit(c(0.65, 0.8), "null")
  grid::grid.draw(panel_plot)
  dev.off()
}



# ==============================================================================
# FIGURE 5 (PANEL E) — Molecular abundance contrasts by Rutherford grade and cluster
# ==============================================================================
# Function to plot molecular abundance contrasts by Rutherford grade and cluster.
plot_5E <- function() {
  d <- fread(
    "./results/VIP-PAD/abundance/Plot_Data.csv"
  )
  d[,
    View := factor(
      ifelse(modality == "Protein", "Proteins", "Metabolites"),
      levels = c("Proteins", "Metabolites")
    )
  ]
  selected <- leading_molecules(loads)
  d[, F2_display := feature %in% selected$feature]
  fg <- d[F2_display == TRUE]
  font <- "Noto Sans"
  # Function to wrap labels to the available width.
  wrap <- function(z, n = 26) {
    vapply(
      z,
      function(s) {
        paste(
          strwrap(gsub("(?<=[-/])", " ", s, perl = TRUE), width = n),
          collapse = "\n"
        )
      },
      character(1)
    )
  }
  annot <- copy(fg)
  annot[, label_plot := wrap(label, 24)]
  met_labels <- copy(annot[View == "Metabolites"])
  setorder(met_labels, y)
  met_labels[, label_x := rep(c(-2.8, 2.8), length.out = .N)]
  met_labels[, label_y := seq(-16, 8, length.out = .N), by = label_x]
  met_labels[, label_plot := wrap(label, 24)]
  # Compact observed ranges; y=x still passes through numerical (0,0).
  p <- ggplot(d, aes(x, y)) +
    scale_x_continuous(
      expand = expansion(mult = c(0.08, 0.08)),
      breaks = scales::breaks_pretty(n = 5)
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0.1, 0.12)),
      breaks = scales::breaks_pretty(n = 6)
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      linetype = "dotted",
      colour = "black",
      linewidth = 0.6
    ) +
    geom_point(colour = "#C8C8C8", size = 0.85, alpha = 0.55) +
    geom_point(data = fg, colour = "#D73027", size = 2.1) +
    geom_text_repel(
      data = annot[View == "Proteins"],
      aes(label = label_plot),
      family = font,
      size = 6.6,
      seed = 43,
      colour = "black",
      max.overlaps = Inf,
      min.segment.length = 0,
      segment.colour = "#555555",
      segment.size = 0.3,
      box.padding = 1.1,
      point.padding = 0.4,
      force = 12,
      force_pull = 0.1,
      max.time = 60,
      max.iter = 2e+05
    ) +
    geom_blank(
      data = data.table(
        View = factor("Metabolites", levels = levels(d$View)),
        x = c(-5.5, 5.5),
        y = c(-18, 10)
      )
    ) +
    geom_segment(
      data = met_labels,
      aes(x = x, y = y, xend = label_x, yend = label_y),
      colour = "#777777",
      linewidth = 0.25
    ) +
    geom_text(
      data = met_labels,
      aes(x = label_x, y = label_y, label = label_plot),
      family = font,
      size = 5,
      lineheight = 0.85
    ) +
    facet_wrap(~View, scales = "free", ncol = 1) +
    labs(
      title = NULL,
      x = "log₂ FC (Rutherford 4–6 vs. 0–3)",
      y = "log₂ FC (Cluster 1 vs. Cluster 2)"
    ) +
    theme_classic(base_family = font, base_size = 20) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        margin = margin(t = 8, b = 4)
      ),
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 20),
      axis.title = element_text(size = 23),
      axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.65),
      axis.ticks.length = unit(2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 27),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      plot.margin = margin(8, 12, 8, 14),
      panel.spacing = unit(1.2, "cm")
    )

  ggsave(
    "./figures-output/Figure5E.pdf",
    p,
    width = 6.2,
    height = 13.8,
    device = cairo_pdf,
    bg = "white"
  )
}



# ==============================================================================
# FIGURE 5 (PANEL F) — Enriched protein pathways
# ==============================================================================
# Function to plot enriched protein pathways.
plot_5F <- function(enrichment_input) {
  enrichment <- copy(enrichment_input$enrichment)
  membership <- readRDS("../inputs/vip-sets.rds")$B
  chosen <- list()
  for (db in c("GO BP", "KEGG", "Reactome")) {
    cand <- enrichment[
      database == db & direction == "all" & FDR < 0.05 & statistic > 0
    ][order(FDR, size, pathway)]
    keep <- character()
    for (set_identifier in cand$set_id) {
      test <- cand[set_id == set_identifier, test_id]
      candidate_members <- which(membership[test, ] == 1)
      redundant <- any(vapply(
        keep,
        # Compare pathway membership with terms already retained for display.
        function(previous_id) {
          previous_members <- which(
            membership[cand[set_id == previous_id, test_id], ] == 1
          )
          length(intersect(previous_members, candidate_members)) /
            length(union(previous_members, candidate_members)) >=
            0.6
        },
        logical(1)
      ))
      if (!redundant) {
        keep <- c(keep, set_identifier)
      }
      if (length(keep) >= 6L) {
        break
      }
    }
    chosen[[paste(db, "all")]] <- cand[set_id %in% keep]
  }
  sel <- rbindlist(chosen)
  ids <- unique(sel$set_id)
  display_data <- enrichment[set_id %in% ids & FDR < 0.05 & statistic > 0]
  rows <- unique(sel[, .(database, set_id, description)])
  rows <- merge(
    rows,
    enrichment[set_id %in% ids, .(best_FDR = min(FDR)), by = set_id],
    by = "set_id"
  )
  rows[, db_order := match(database, c("GO BP", "KEGG", "Reactome"))]
  rows[,
    gene_ratio := enrichment$size[match(set_id, enrichment$set_id)] /
      ncol(membership)
  ]
  setorderv(
    rows,
    c("db_order", "gene_ratio", "best_FDR", "set_id"),
    c(1, -1, 1, 1)
  )
  rows[,
    display := vapply(
      description,
      function(description_text) {
        if (
          grepl("^Regulation of Insulin-like Growth Factor", description_text)
        ) {
          description_text <- "IGF transport and uptake by IGFBPs"
        }
        paste(strwrap(description_text, width = 43), collapse = "\n")
      },
      character(1)
    )
  ]
  display_data[, y := factor(set_id, levels = rev(rows$set_id))]
  display_data[,
    database := factor(database, levels = c("GO BP", "KEGG", "Reactome"))
  ]
  display_data[, logFDR := -log10(pmax(FDR, .Machine$double.xmin))]
  labels <- setNames(rows$display, rows$set_id)
  display_data[, gene_ratio := size / ncol(membership)]
  pathway_plot <- ggplot(display_data, aes(x = gene_ratio, y = y)) +
    geom_point(
      aes(size = size, fill = logFDR),
      shape = 21,
      colour = "#333333",
      stroke = 0.3
    ) +
    facet_grid(database ~ ., scales = "free_y", space = "free_y", drop = TRUE) +
    scale_y_discrete(labels = labels, expand = expansion(add = 0.65)) +
    scale_x_continuous(
      limits = c(0, ceiling(max(display_data$gene_ratio) * 1.15 / 0.02) * 0.02),
      breaks = seq(
        0,
        ceiling(
          max(display_data$gene_ratio) *
            1.15 /
            0.02
        ) *
          0.02,
        0.02
      ),
      labels = function(v) sprintf("%.2f", v),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    scale_size_area(max_size = 6.3, name = "Count", breaks = c(20, 40, 80)) +
    scale_fill_gradient(
      low = "#F6C1C7",
      high = "#8D0A15",
      limits = c(0, 8),
      breaks = c(0, 2, 4, 6, 8),
      oob = scales::squish,
      name = expression(-log[10](FDR))
    ) +
    labs(x = "Gene ratio", y = NULL) +
    theme_bw(base_size = 12.5, base_family = "Noto Sans") +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.text.y = element_text(size = 12),
      axis.ticks.y = element_line(linewidth = 0.3),
      axis.ticks.length = grid::unit(1.2, "mm"),
      panel.grid.major.x = element_line(colour = "#E7E7E7", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(colour = "#E7E7E7", linewidth = 0.35),
      panel.border = element_rect(colour = "#4D4D4D", linewidth = 0.4),
      panel.spacing.x = grid::unit(5, "mm"),
      panel.spacing.y = grid::unit(2.5, "mm"),
      strip.background = element_rect(
        fill = "#E7E7E7",
        colour = "#4D4D4D",
        linewidth = 0.4
      ),
      strip.text.x = element_text(size = 12),
      strip.text.y = element_text(size = 12),
      legend.position = "right",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10.5),
      legend.key.height = grid::unit(6, "mm"),
      legend.spacing.y = grid::unit(4, "mm"),
      legend.box.spacing = grid::unit(3, "mm"),
      plot.margin = margin(8, 8, 8, 8)
    ) +
    guides(
      size = guide_legend(
        order = 2,
        override.aes = list(fill = "white", shape = 21)
      ),
      fill = guide_colourbar(
        order = 3,
        barheight = grid::unit(30, "mm"),
        barwidth = grid::unit(4, "mm")
      )
    )
  lines <- sum(vapply(strsplit(rows$display, "\n"), length, integer(1)))
  height <- max(170, 40 + nrow(rows) * 5.5 + (lines - nrow(rows)) * 3.6)
  pathway_plot <- pathway_plot +
    theme(
      legend.position = "bottom",
      legend.location = "plot",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.spacing.x = grid::unit(4, "mm"),
      legend.spacing.y = grid::unit(1, "mm"),
      legend.box.spacing = grid::unit(3, "mm")
    ) +
    guides(
      size = guide_legend(
        order = 2,
        nrow = 1,
        title.position = "top",
        override.aes = list(fill = "white", shape = 21)
      ),
      fill = guide_colourbar(
        order = 3,
        direction = "horizontal",
        title.position = "top",
        barwidth = grid::unit(35, "mm"),
        barheight = grid::unit(4, "mm")
      )
    )
  height <- height + 15

  ggsave(
    "./figures-output/Figure5F.pdf",
    pathway_plot,
    device = cairo_pdf,
    width = 185,
    height = height,
    units = "mm",
    family = "Noto Sans",
    limitsize = FALSE
  )
}



# ==============================================================================
# FIGURE 5 (PANEL G) — Event-free survival across factor tertiles
# ==============================================================================
# Function to plot event-free survival across factor tertiles.
plot_5G <- function() {
  a <- readRDS("./results/VIP-PAD/outcomes/Factor_2_MACCE/km_results.rds")
  lev <- c("Low", "Medium", "High")
  pal <- c(Low = "#42A5D9", Medium = "#A88F00", High = "#F25A24")
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

  W <- 5.95
  H <- 5.3
  left <- 1.55
  right <- 5.63
  bottom <- 2.26
  top <- 4.72
  # Function to map a value to its horizontal panel coordinate.
  xx <- function(t) left + (right - left) * t / 36
  ylo <- 0.4
  # Function to map a value to its vertical panel coordinate.
  yy <- function(v) bottom + (top - bottom) * (v - ylo) / (1 - ylo)

  cairo_pdf(
    "./figures-output/Figure5G.pdf",
    width = W,
    height = H,
    family = "Noto Sans"
  )
  grid.newpage()
  txt("MACCE", (left + right) / 2, 5.1, 14)
  ln(left, bottom, left, top, "#333333", 0.7)
  ln(left, bottom, right, bottom, "#333333", 0.7)
  yticks <- seq(0.4, 1, 0.1)
  for (v in yticks) {
    ln(left - 0.04, yy(v), left, yy(v), "#555555", 0.6)
    txt(sprintf("%.1f", v), left - 0.09, yy(v), 12, "right")
  }
  for (t in c(0, 12, 24, 36)) {
    ln(xx(t), bottom, xx(t), bottom - 0.04, "#555555", 0.6)
    txt(as.character(t), xx(t), bottom - 0.16, 12)
  }
  txt("Event-free survival", left - 0.55, (top + bottom) / 2, 13, rot = 90)
  txt("Months", (left + right) / 2, 1.82, 13)
  ds <- list()
  for (g in lev) {
    d <- rbind(
      data.table(time = 0, S = 1, lower = 1, upper = 1, n_censor = 0),
      a$curves[
        endpoint == "MACCE" & Tertile == g & time <= 36,
        .(time, S, lower, upper, n_censor)
      ]
    )
    setorder(d, time)
    last <- copy(tail(d, 1))
    last[, time := 36]
    d <- rbind(d, last)
    ds[[g]] <- d
    tx <- rep(d$time, each = 2)[-1]
    lower <- rep(d$lower, each = 2)
    lower <- lower[-length(lower)]
    upper <- rep(d$upper, each = 2)
    upper <- upper[-length(upper)]
    grid.polygon(
      x = unit(c(xx(tx), rev(xx(tx))), "in"),
      y = unit(c(yy(lower), rev(yy(upper))), "in"),
      gp = gpar(fill = adjustcolor(pal[g], alpha.f = 0.16), col = NA)
    )
  }
  for (g in lev) {
    d <- ds[[g]]
    tx <- rep(d$time, each = 2)[-1]
    sy <- rep(d$S, each = 2)
    sy <- sy[-length(sy)]
    grid.lines(
      x = unit(xx(tx), "in"),
      y = unit(yy(sy), "in"),
      gp = gpar(col = pal[g], lwd = 2, linejoin = "round")
    )
    cen <- d[n_censor > 0 & time < 36]
    if (nrow(cen)) {
      grid.segments(
        x0 = unit(xx(cen$time), "in"),
        x1 = unit(xx(cen$time), "in"),
        y0 = unit(
          yy(cen$S) -
            0.025,
          "in"
        ),
        y1 = unit(yy(cen$S) + 0.025, "in"),
        gp = gpar(col = pal[g], lwd = 0.6)
      )
    }
  }
  for (i in 1:3) {
    lx <- left + 0.45 + (i - 1) * 1.15
    ly <- top + 0.15
    rect(
      lx - 0.16,
      ly - 0.045,
      0.22,
      0.09,
      col = NA,
      fill = adjustcolor(pal[i], alpha.f = 0.3)
    )
    ln(lx - 0.16, ly, lx + 0.06, ly, pal[i], 2)
    txt(lev[i], lx + 0.12, ly, 11, "left")
  }
  te <- a$tests[endpoint == "MACCE"]
  pv <- if (te$P < 0.001) {
    expression(paste("Log-rank ", italic(P), " < 0.001"))
  } else {
    bquote(paste("Log-rank ", italic(P), " = ", .(sprintf("%.3f", te$P))))
  }
  txt(pv, left + 0.12, bottom + 0.5, 11, "left")
  rleft <- left -
    0.55 -
    convertHeight(
      grobHeight(textGrob(
        "Event-free survival",
        gp = gpar(fontfamily = "Noto Sans", fontsize = 13)
      )),
      "in",
      valueOnly = TRUE
    ) /
      2
  txt("Number at risk", rleft, 1.54, 12, "left")
  rect(rleft, 0.2, right - rleft, 1.17, col = "#888888")
  rx <- seq(3.75, right - 0.12, length.out = 4)
  txt("Months", rleft + 0.1, 1.21, 10.5, "left")
  for (j in 1:4) {
    txt(as.character(c(0, 12, 24, 36)[j]), rx[j], 1.21, 10.5)
  }
  for (i in 1:3) {
    g <- lev[i]
    ry <- 1.21 - i * 0.28
    ct <- a$counts[endpoint == "MACCE" & Tertile == g]
    txt(
      sprintf("%s (n = %d; events = %d)", g, ct$N, ct$Events),
      rleft + 0.1,
      ry,
      10.5,
      "left"
    )
    for (j in 1:4) {
      txt(
        as.character(a$risk[
          endpoint == "MACCE" & Tertile == g & time == c(0, 12, 24, 36)[j],
          n_risk
        ]),
        rx[j],
        ry,
        11
      )
    }
  }
  dev.off()
}



# ==============================================================================
# FIGURE 5 (PANELS H–I) — SHAP ranking in training and molecular contributions in held-out patients
# ==============================================================================
# Function to plot training-set SHAP ranks and holdout SHAP values.
plot_5HI <- function() {
  set.seed(42)
  # Function to position text in panel coordinates.
  txt <- function(
    label,
    x,
    y,
    size = 12,
    just = "centre",
    rot = 0,
    col = "black"
  ) {
    grid.text(
      label,
      x = unit(x, "in"),
      y = unit(y, "in"),
      just = just,
      rot = rot,
      gp = gpar(fontfamily = "Noto Sans", fontsize = size, col = col)
    )
  }
  # Function to draw a line in panel coordinates.
  ln <- function(x, y, x2, y2, col = "#E8E8E8", lwd = 0.7, lty = 1) {
    grid.lines(
      x = unit(c(x, x2), "in"),
      y = unit(c(y, y2), "in"),
      gp = gpar(col = col, lwd = lwd, lty = lty)
    )
  }
  # Function to wrap long labels to fit the panel.
  wrap_label <- function(x, width = 38) {
    # Wrap chemical names at hyphens.
    parts <- strsplit(x, "(?<=-)|(?<= )", perl = TRUE)[[1]]
    ans <- character()
    line <- ""
    for (p in parts) {
      if (nchar(paste0(line, p)) > width && nchar(line) > 0) {
        ans <- c(ans, trimws(line))
        line <- p
      } else {
        line <- paste0(line, p)
      }
    }
    paste(c(ans, trimws(line)), collapse = "\n")
  }
  a <- readRDS(
    "./results/VIP-PAD/outcomes/Factor_2_MACCE/MACCE_shap_results.rds"
  )
  imp <- a$importance
  n <- nrow(imp)

  # FIGURE 5 (PANEL H) — Training-set SHAP ranking
  cairo_pdf(
    "./figures-output/Figure5H.pdf",
    width = 4.2,
    height = 4.2,
    family = "Noto Sans"
  )
  grid.newpage()
  left <- 0.85
  right <- 3.9
  bottom <- 0.72
  top <- 3.62
  mx <- max(imp$Mean_abs_SHAP) * 1.08
  # Function to map a value to its horizontal panel coordinate.
  xx <- function(v) left + (right - left) * v / n
  # Function to map a value to its vertical panel coordinate.
  yy <- function(v) bottom + (top - bottom) * v / mx
  for (v in pretty(c(0, mx), 5)) {
    if (v >= 0 && v <= mx) {
      ln(left, yy(v), right, yy(v))
      txt(format(v, digits = 2), left - 0.12, yy(v), 11, "right")
    }
  }
  for (v in c(0, 50, 100, 150, 200, n)) {
    ln(xx(v), bottom, xx(v), top)
    txt(as.character(v), xx(v), bottom - 0.16, 11)
  }
  ln(left, bottom, right, bottom, "#444444")
  ln(left, bottom, left, top, "#444444")
  grid.lines(
    x = unit(xx(imp$rank), "in"),
    y = unit(yy(imp$Mean_abs_SHAP), "in"),
    gp = gpar(col = "black", lwd = 0.8)
  )
  grid.points(
    x = unit(xx(imp$rank), "in"),
    y = unit(yy(imp$Mean_abs_SHAP), "in"),
    pch = 16,
    size = unit(0.035, "in"),
    gp = gpar(col = "black")
  )
  ep <- imp[rank == a$elbow]
  grid.points(
    x = unit(xx(ep$rank), "in"),
    y = unit(yy(ep$Mean_abs_SHAP), "in"),
    pch = 16,
    size = unit(0.095, "in"),
    gp = gpar(col = "#E32620")
  )
  txt(
    sprintf("Elbow point = %d molecules", a$elbow),
    xx(ep$rank) + 0.18,
    yy(ep$Mean_abs_SHAP) + 0.17,
    12,
    "left"
  )
  txt("MACCE: SHAP elbow detection", (left + right) / 2, 3.94, 13)
  txt("Molecule rank", (left + right) / 2, 0.2, 12)
  txt("Mean |SHAP|", 0.22, (bottom + top) / 2, 12, rot = 90)
  dev.off()
  # Match molecular labels and values by identifier.
  feat <- a$panel
  k <- length(feat)
  labs <- imp[match(feat, feature), label]
  labs <- gsub("\\s*\\[L[0-9]+\\]$", "", labs)
  labs <- vapply(labs, wrap_label, character(1))
  rowheight <- max(0.44, max(lengths(strsplit(labs, "\n"))) * 0.17 + 0.1)
  H <- 1.65 + k * rowheight
  W <- 8.4
  left <- 3.85
  right <- 7.35
  bottom <- 0.95
  top <- H - 0.6
  # Function to position a row within the panel.
  y <- function(i) top - (i - 0.5) * rowheight
  sh <- a$holdout_SHAP[, feat, drop = FALSE]
  X <- a$holdout_X[, feat, drop = FALSE]
  span <- max(abs(sh)) * 1.1
  xr <- c(-span, span)
  # Function to map a value to its horizontal panel coordinate.
  xx <- function(v) left + (right - left) * (v + span) / (2 * span)
  cp <- colorRampPalette(c("#2315BE", "#B799E9", "#F4D7E6", "#FF2318"))(256)

  # FIGURE 5 (PANEL I) — Holdout molecular contributions
  cairo_pdf(
    "./figures-output/Figure5I.pdf",
    width = W,
    height = H,
    family = "Noto Sans"
  )
  grid.newpage()
  for (i in seq_len(k)) {
    ln(left, y(i), right, y(i), "#EDEDED", 0.6, 3)
  }
  ln(xx(0), bottom, xx(0), top, "#777777", 0.7)
  for (i in seq_len(k)) {
    vals <- sh[, i]
    v <- X[, i]
    coloridx <- 1 +
      round(255 * (rank(v, ties.method = "average") - 1) / (length(v) - 1))
    bins <- cut(
      vals,
      breaks = seq(-span, span, length.out = 55),
      include.lowest = TRUE,
      labels = FALSE
    )
    off <- numeric(length(vals))
    maxbin <- max(table(bins))
    for (b in unique(bins)) {
      ii <- which(bins == b)
      ii <- sample(ii, length(ii))
      m <- length(ii)
      oo <- seq_len(m) - (m + 1) / 2
      off[ii] <- oo / max(1, maxbin / 2) * rowheight * 0.35
    }
    # Draw points in a reproducible shuffled order.
    oo <- sample(seq_along(vals))
    grid.points(
      x = unit(xx(vals[oo]), "in"),
      y = unit(y(i) + off[oo], "in"),
      pch = 16,
      size = unit(0.041, "in"),
      gp = gpar(col = adjustcolor(cp[coloridx[oo]], alpha.f = 0.78))
    )
    txt(labs[i], left - 0.12, y(i), 12, "right")
  }
  ln(left, bottom, right, bottom, "#444444", 0.7)
  ticks <- pretty(xr, 5)
  ticks <- ticks[ticks >= -span & ticks <= span]
  for (v in ticks) {
    ln(xx(v), bottom, xx(v), bottom - 0.05, "#444444", 0.6)
    txt(sprintf("%.1f", v), xx(v), bottom - 0.19, 11)
  }
  txt("SHAP value (log event-free time)", (left + right) / 2, 0.55, 12)
  txt("Shorter", left + 0.1, 0.37, 10.5, "left")
  txt("Longer", right - 0.1, 0.37, 10.5, "right")
  grid.lines(
    x = unit(c(left + 0.38, left + 0.06), "in"),
    y = unit(c(0.55, 0.55), "in"),
    arrow = arrow(length = unit(0.07, "in"), ends = "last"),
    gp = gpar(col = "black", lwd = 0.7)
  )
  grid.lines(
    x = unit(c(right - 0.38, right - 0.06), "in"),
    y = unit(c(0.55, 0.55), "in"),
    arrow = arrow(length = unit(0.07, "in"), ends = "last"),
    gp = gpar(col = "black", lwd = 0.7)
  )
  txt("MACCE", (left + right) / 2, H - 0.22, 13)
  cbx <- 7.66
  cby <- (top + bottom) / 2 - 1
  cbh <- 2
  grid.raster(
    as.raster(matrix(rev(cp), ncol = 1)),
    x = unit(cbx, "in"),
    y = unit(cby, "in"),
    width = unit(0.13, "in"),
    height = unit(cbh, "in"),
    just = c("left", "bottom"),
    interpolate = FALSE
  )
  txt("High", cbx + 0.2, cby + cbh, 10, "left")
  txt("Low", cbx + 0.2, cby, 10, "left")
  txt("Molecular level", cbx + 0.55, cby + cbh / 2, 11, rot = 90)
  dev.off()
}



# ==============================================================================
# FIGURE 5 (PANELS A–I) — Prepare Factor 2 loadings
# ==============================================================================

fit <- readRDS(
  "./results/VIP-PAD/mofa/Selected_Fit.rds"
)
metabolites <- readRDS("../inputs/metabolites.rds")$metabolite

loads <- factor_loadings()



# ==============================================================================
# FIGURE 5 (PANELS A–I) — Draw and save panels
# ==============================================================================

dir.create("./figures-output", recursive = TRUE, showWarnings = FALSE)
# FIGURE 5 (PANEL A)
plot_5A()
# FIGURE 5 (PANEL B)
plot_5B()
factor2 <- readRDS("./results/VIP-PAD/enrichment/F2.rds")
# FIGURE 5 (PANEL C)
plot_5C(factor2)
# FIGURE 5 (PANEL D)
plot_5D(factor2)
# FIGURE 5 (PANEL E)
plot_5E()
# FIGURE 5 (PANEL F)
plot_5F(factor2)
# FIGURE 5 (PANEL G)
plot_5G()
# FIGURE 5 (PANELS H–I)
plot_5HI()
