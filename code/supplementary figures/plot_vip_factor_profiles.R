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
# Global: S12, S15 and S18 show the molecular profiles of VIP-PAD Factors 3, 4 and 8.
#
# Inputs: VIP-PAD measurements, factor loadings, clusters, abundance contrasts and enrichment results.
# Outputs: Individual panel PDFs under ./figures-output/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(grid)
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(ggbeeswarm)
})



# ==============================================================================
# SHARED FUNCTIONS
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

# Function to draw a line in panel coordinates.
panel_line <- function(x, y, x2, y2, col = "#C7C7C7", width = 0.6, lty = 1) {
  grid.lines(
    x = unit(c(x, x2), "in"),
    y = unit(c(y, y2), "in"),
    gp = gpar(col = col, lwd = width, lty = lty)
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

# Function to extract protein and metabolite loadings for one VIP-PAD factor.
factor_loadings <- function(fit, factor_name, metabolites) {
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
      loading_raw = w[, factor_name]
    )
  }))
  loads[,
    relative_weight := abs(loading_raw) / max(abs(loading_raw)),
    by = modality
  ]
  loads
}

# Function to plot enriched VIP-PAD protein pathways.
plot_protein_pathways <- function(enrichment_input, output, axis_padding, bottom_margin) {
  enrichment <- copy(enrichment_input$enrichment)
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
    scale_y_discrete(labels = labels, expand = expansion(add = axis_padding)) +
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
      plot.margin = margin(8, 8, bottom_margin, 8)
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
    output,
    pathway_plot,
    device = cairo_pdf,
    width = 205,
    height = height,
    units = "mm",
    family = "Noto Sans",
    limitsize = FALSE
  )
}

# Function to plot signed molecular loadings for VIP-PAD Factors 4 and 8.
plot_signed_loadings <- function(loads, output) {
  d <- loads[, .(
    feature,
    label,
    weight = loading_raw,
    view = ifelse(modality == "Protein", "Proteomics", "Metabolomics")
  )]
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
  wrap_name <- function(z, width = 26L) {
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
    function(z) sign(z) * abs(z)^(1 / 1.8)
  )
  p <- ggplot(d, aes(scaled_weight, position)) +
    geom_vline(
      xintercept = 0,
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
    geom_point(data = d[selected == TRUE], aes(colour = direction), size = 2.6) +
    geom_text_repel(
      data = d[selected == TRUE & view == "Proteomics"],
      aes(label = display_label, colour = direction),
      family = "Noto Sans",
      size = 4.1,
      seed = 52,
      max.overlaps = Inf,
      box.padding = 0.65,
      force = 3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.colour = "grey60",
      segment.size = 0.25,
      max.iter = 1e+05,
      max.time = 10
    ) +
    geom_text_repel(
      data = d[selected == TRUE & view == "Metabolomics"],
      aes(label = display_label, colour = direction),
      family = "Noto Sans",
      size = 3.4,
      seed = 52,
      direction = "both",
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.colour = "grey60",
      segment.size = 0.25,
      max.iter = 1e+05,
      max.time = 10
    ) +
    facet_wrap(~view, ncol = 1) +
    scale_y_continuous(limits = c(-0.72, 0.72), breaks = NULL) +
    scale_x_continuous(
      transform = trans,
      limits = c(-1.04, 1.04),
      breaks = seq(-1, 1, 0.2),
      labels = function(z) format(round(z, 1), nsmall = 1)
    ) +
    scale_colour_manual(values = c(Positive = "#C85B63", Negative = "#526EA4")) +
    labs(
      title = NULL,
      x = "Scaled weight",
      y = "View",
      colour = "Weight direction"
    ) +
    theme_classic(base_family = "Noto Sans", base_size = 15) +
    theme(
      axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      plot.title = element_text(size = 17, hjust = 0.5),
      panel.spacing = grid::unit(5, "mm"),
      legend.position = "none",
      plot.margin = margin(8, 12, 5, 8)
    )
  ggsave(
    output,
    p,
    width = 8.5,
    height = 5.8,
    device = cairo_pdf
  )
}

# Function to compare abundance by molecular cluster and Rutherford category.
plot_abundance_contrasts <- function(selected_features, output, label_range) {
  d <- copy(abundance_contrasts)
  d[,
    View := factor(
      ifelse(modality == "Protein", "Proteins", "Metabolites"),
      levels = c("Proteins", "Metabolites")
    )
  ]
  d[, display := feature %in% selected_features]
  fg <- d[display == TRUE]
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
  met_labels[, label_y := seq(label_range[1], label_range[2], length.out = .N), by = label_x]
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
      size = 5.4,
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
      size = 4.8,
      lineheight = 0.8
    ) +
    facet_wrap(~View, scales = "free", nrow = 1) +
    labs(
      title = NULL,
      x = "log₂ FC (Rutherford 4–6 vs. 0–3)",
      y = "log₂ FC (Cluster 1 vs. Cluster 2)"
    ) +
    theme_classic(base_family = font, base_size = 17) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        margin = margin(t = 8, b = 4)
      ),
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 17),
      axis.title = element_text(size = 20),
      axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.65),
      axis.ticks.length = unit(2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 23),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      plot.margin = margin(8, 12, 8, 8),
      panel.spacing = unit(1.2, "cm")
    )
  ggsave(
    output,
    p,
    width = 12.2,
    height = 6.5,
    device = cairo_pdf,
    bg = "white"
  )
}

# Function to plot metabolite loadings and lipid-class enrichment.
plot_lipid_classes <- function(enrichment_input, output, axis_label) {
  class_results <- copy(enrichment_input$lipid_results)
  lipid_molecules <- copy(enrichment_input$lipid_molecules)
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
  lipid_padding <- diff(range(lipid_molecules$loading_raw)) * 0.12
  stars[, y := max(lipid_molecules$loading_raw) + lipid_padding]
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
  # Match the compact Factor 1 lipid panel while enlarging text for assembly.
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
        min(lipid_molecules$loading_raw) - lipid_padding,
        max(lipid_molecules$loading_raw) +
          2 * lipid_padding
      ),
      expand = expansion(mult = 0)
    ) +
    labs(x = NULL, y = axis_label) +
    theme_classic(base_family = "Noto Sans", base_size = 15) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 16),
      axis.title = element_text(size = 18),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(12, 10, 5, 8)
    )
  ggsave(
    output,
    lipid_plot,
    device = cairo_pdf,
    width = 180,
    height = 115,
    units = "mm"
  )
}



# ==============================================================================
# FIGURES S12, S15 AND S18 (PANELS A–F) — Shared inputs
# ==============================================================================

fit <- readRDS("./results/VIP-PAD/mofa/Selected_Fit.rds")
cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
cl <- readRDS("./results/VIP-PAD/clustering/Clustering_Stability.rds")
metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
membership <- readRDS("../inputs/vip-sets.rds")$B
abundance_contrasts <- data.table::fread("./results/VIP-PAD/abundance/Plot_Data.csv")



# ==============================================================================
# FIGURE S12 (PANELS A–F) — Molecular profile of VIP-PAD Factor 3
# ==============================================================================

local({
  # FIGURE S12 (PANEL A) — Highest absolute factor loadings
  loads <- factor_loadings(fit, "Factor3", metabolites)
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:15],
    loads[modality == "Metabolite" & relative_weight > 0.8][order(
      -abs(loading_raw),
      feature
    )]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  ca <- copy(sel)
  ca[, feature_key := factor(feature, levels = feature)]
  ca[, label_plot := ifelse(modality == "Metabolite", wrap_label(label, 16), label)]
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
      y = unit(-1.1, "in"),
      width = unit(14.5, "in"),
      height = unit(14.5, "in"),
      just = c("left", "bottom")
    ))
    print(pa, newpage = FALSE)
    txt("Factor 3 top features", 0.48, 7.8, 36, rot = 90)
    panel_rect(0.65, 2.85, 0.32, 0.32, "#4575b4")
    txt("Proteins", 1.1, 2.69, 28, "left")
    panel_rect(0.65, 2.28, 0.32, 0.32, "#f46d43")
    txt("Metabolites", 1.1, 2.12, 28, "left")
    popViewport()
  }
  dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
  cairo_pdf(
    "./figures-output/FigureS12A.pdf",
    width = 14.5,
    height = 13.4,
    family = font
  )
  draw()
  dev.off()
  selected_features <- sel$feature

  # FIGURE S12 (PANEL B) — Molecular abundance by cluster
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:8],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:4]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  # Align molecular data and fixed cluster memberships by ID.
  ids <- rownames(fit$Z)
  ix <- match(ids, cohort$study_id)
  cx <- match(ids, cl$assignments$sample_id)
  groups <- cl$assignments$K2[cx]
  abundance <- as.matrix(cohort[ix, sel$feature, drop = FALSE])
  rownames(abundance) <- ids
  h <- list(
    groups = groups,
    patient_order = list(F03 = order(groups, fit$Z[, "Factor3"], ids)),
    molecular_matrix = abundance
  )
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  # Function to draw the cluster abundance heatmap.
  draw <- function() {
    grid.newpage()
    # Match the Factor 2 row height and font scale for 12 displayed features.
    pushViewport(viewport(
      x = unit(-1.65, "in"),
      y = 0,
      width = unit(15.7, "in"),
      height = unit(8.7, "in"),
      just = c("left", "bottom")
    ))
    ord <- h$patient_order$F03
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
      panel_line(x0 + gap, 8.09, x1 - gap, 8.09, "black", 1)
      panel_line(x0 + gap, 8.09, x0 + gap, 7.97, "black", 1)
      panel_line(x1 - gap, 8.09, x1 - gap, 7.97, "black", 1)
      txt(paste("Cluster", k), (x0 + x1) / 2, 8.32, 15)
    }
    for (v in c("Protein", "Metabolite")) {
      q <- sel[modality == v]
      nr <- nrow(q)
      top <- if (v == "Protein") {
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
      panel_line(left + boundary, top, left + boundary, top - nr * rh, "white", 1)
      for (vv in c(0, 0.3, 0.6, 0.9)) {
        xx <- wx + ww * vv / 1.1
        panel_line(xx, top, xx, top - nr * rh, "#C3C3C3", 0.8)
        txt(sprintf("%.1f", vv), xx, top - nr * rh - 0.18, 11)
      }
      for (j in 0:nr) {
        panel_line(wx, top - j * rh, wx + ww, top - j * rh, "#DDDDDD", 0.6)
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
              16
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
        panel_line(wx, yy, xe, yy, "black", 1.8)
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
    "./figures-output/FigureS12B.pdf",
    width = 14,
    height = 8.7,
    family = font
  )
  draw()
  dev.off()

  # FIGURE S12 (PANEL C) — Signed molecular loadings
  d <- loads[, .(
    feature,
    label,
    weight = loading_raw,
    view = ifelse(modality == "Protein", "Proteomics", "Metabolomics")
  )]
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
  d[,
    display_label := vapply(
      label,
      function(z) paste(strwrap(z, width = 20), collapse = "\n"),
      character(1)
    )
  ]
  trans <- scales::trans_new(
    "signed_power_1.8",
    function(z) sign(z) * abs(z)^1.8,
    function(z) sign(z) * abs(z)^(1 / 1.8)
  )
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
    geom_point(data = d[selected == TRUE], aes(colour = direction), size = 2.6) +
    geom_text_repel(
      data = d[selected == TRUE & view == "Proteomics"],
      aes(label = display_label, colour = direction),
      family = "Noto Sans",
      size = 4.1,
      seed = 52,
      max.overlaps = Inf,
      box.padding = 0.65,
      force = 3,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.colour = "grey60",
      segment.size = 0.25,
      max.iter = 1e+05,
      max.time = 10
    ) +
    geom_text_repel(
      data = d[selected == TRUE & view == "Metabolomics"],
      aes(label = display_label, colour = direction),
      family = "Noto Sans",
      size = 3.1,
      seed = 52,
      direction = "y",
      force = 15,
      force_pull = 0.05,
      nudge_x = -d[
        selected == TRUE &
          view == "Metabolomics"
      ]$position,
      max.overlaps = Inf,
      box.padding = 0.6,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.colour = "grey60",
      segment.size = 0.25,
      max.iter = 1e+05,
      max.time = 10
    ) +
    facet_wrap(~view, nrow = 1) +
    scale_x_continuous(limits = c(-0.72, 0.72), breaks = NULL) +
    scale_y_continuous(
      transform = trans,
      limits = c(-1.04, 1.04),
      breaks = seq(-1, 1, 0.2),
      labels = function(z) format(round(z, 1), nsmall = 1)
    ) +
    scale_colour_manual(values = c(Positive = "#C85B63", Negative = "#526EA4")) +
    labs(
      title = NULL,
      x = "View",
      y = "Scaled weight",
      colour = "Weight direction"
    ) +
    theme_classic(base_family = "Noto Sans", base_size = 15) +
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(size = 16),
      plot.title = element_text(size = 17, hjust = 0.5),
      panel.spacing = grid::unit(5, "mm"),
      legend.position = "none",
      plot.margin = margin(8, 12, 5, 8)
    )
  panel_plot <- ggplotGrob(p)
  panel_columns <- unique(panel_plot$layout$l[grepl(
    "^panel",
    panel_plot$layout$name
  )])
  panel_plot$widths[panel_columns] <- grid::unit(c(0.65, 0.8), "null")
  ggsave(
    "./figures-output/FigureS12C.pdf",
    panel_plot,
    width = 5.2,
    height = 6.6,
    device = cairo_pdf
  )

  # FIGURE S12 (PANEL D) — Cluster and Rutherford abundance contrasts
  d <- copy(abundance_contrasts)
  d[,
    View := factor(
      ifelse(modality == "Protein", "Proteins", "Metabolites"),
      levels = c("Proteins", "Metabolites")
    )
  ]
  d[, F3_display := feature %in% selected_features]
  fg <- d[F3_display == TRUE]
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
      size = 5.4,
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
      size = 4.8,
      lineheight = 0.8
    ) +
    facet_wrap(~View, scales = "free", nrow = 1) +
    labs(
      title = NULL,
      x = "log₂ FC (Rutherford 4–6 vs. 0–3)",
      y = "log₂ FC (Cluster 1 vs. Cluster 2)"
    ) +
    theme_classic(base_family = font, base_size = 17) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(
        size = 20,
        face = "bold",
        hjust = 0,
        margin = margin(t = 8, b = 4)
      ),
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 17),
      axis.title = element_text(size = 20),
      axis.line = element_line(linewidth = 0.65),
      axis.ticks = element_line(linewidth = 0.65),
      axis.ticks.length = unit(2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(size = 23),
      legend.position = "bottom",
      legend.text = element_text(size = 12),
      plot.margin = margin(8, 12, 8, 8),
      panel.spacing = unit(1.2, "cm")
    )
  ggsave(
    "./figures-output/FigureS12D.pdf",
    p,
    width = 12.2,
    height = 6.5,
    device = cairo_pdf,
    bg = "white"
  )

  # FIGURE S12 (PANEL F) — Enriched protein pathways
  enrichment_input <- readRDS("./results/VIP-PAD/enrichment/F3.rds")
  plot_protein_pathways(
    enrichment_input, "./figures-output/FigureS12F.pdf",
    axis_padding = 1, bottom_margin = 8
  )

  # FIGURE S12 (PANEL E) — Metabolite loadings by lipid class
  lipid_molecules <- copy(enrichment_input$lipid_molecules)
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
  lipid_padding <- diff(range(lipid_molecules$loading_raw)) * 0.12
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
    scale_colour_manual(values = palette) +
    scale_x_discrete(labels = class_order) +
    scale_y_continuous(
      limits = c(
        min(lipid_molecules$loading_raw) -
          lipid_padding,
        max(lipid_molecules$loading_raw) + 2 * lipid_padding
      ),
      expand = expansion(mult = 0)
    ) +
    labs(x = NULL, y = "Factor 3 weight") +
    theme_classic(base_family = "Noto Sans", base_size = 13) +
    theme(
      text = element_text(colour = "black"),
      axis.text = element_text(colour = "black", size = 13),
      axis.title = element_text(size = 15),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4),
      legend.position = "none",
      plot.margin = margin(8, 10, 5, 8)
    )
  ggsave(
    "./figures-output/FigureS12E.pdf",
    lipid_plot,
    device = cairo_pdf,
    width = 180,
    height = 115,
    units = "mm"
  )
})



# ==============================================================================
# FIGURE S15 (PANELS A–F) — Molecular profile of VIP-PAD Factor 4
# ==============================================================================

local({
  # FIGURE S15 (PANEL A) — Highest absolute factor loadings
  loads <- factor_loadings(fit, "Factor4", metabolites)
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:6],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:10]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  ca <- copy(sel)
  ca[, feature_key := factor(feature, levels = feature)]
  ca[, label_plot := ifelse(modality == "Metabolite", wrap_label(label, 16), label)]
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
      axis.text.theta = element_text(size = 21, colour = "black"),
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
      y = unit(-1.1, "in"),
      width = unit(14.5, "in"),
      height = unit(14.5, "in"),
      just = c("left", "bottom")
    ))
    print(pa, newpage = FALSE)
    txt("Factor 4 top features", 0.48, 7.8, 32, rot = 90)
    panel_rect(0.65, 2.85, 0.32, 0.32, "#4575b4")
    txt("Proteins", 1.1, 2.69, 24, "left")
    panel_rect(0.65, 2.28, 0.32, 0.32, "#f46d43")
    txt("Metabolites", 1.1, 2.12, 24, "left")
    popViewport()
  }
  dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
  cairo_pdf(
    "./figures-output/FigureS15A.pdf",
    width = 14.5,
    height = 13.4,
    family = font
  )
  draw()
  dev.off()
  selected_features <- sel$feature

  # FIGURE S15 (PANEL B) — Molecular abundance by cluster
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:4],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:8]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  # Align molecular data and fixed cluster memberships by ID.
  ids <- rownames(fit$Z)
  ix <- match(ids, cohort$study_id)
  cx <- match(ids, cl$assignments$sample_id)
  groups <- cl$assignments$K2[cx]
  abundance <- as.matrix(cohort[ix, sel$feature, drop = FALSE])
  rownames(abundance) <- ids
  h <- list(
    groups = groups,
    patient_order = list(F04 = order(groups, fit$Z[, "Factor4"], ids)),
    molecular_matrix = abundance
  )
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
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
    ord <- h$patient_order$F04
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
      panel_line(x0 + gap, 8.09, x1 - gap, 8.09, "black", 1)
      panel_line(x0 + gap, 8.09, x0 + gap, 7.97, "black", 1)
      panel_line(x1 - gap, 8.09, x1 - gap, 7.97, "black", 1)
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
      panel_line(left + boundary, top, left + boundary, top - nr * rh, "white", 1)
      for (vv in c(0, 0.3, 0.6, 0.9)) {
        xx <- wx + ww * vv / 1.1
        panel_line(xx, top, xx, top - nr * rh, "#C3C3C3", 0.8)
        txt(sprintf("%.1f", vv), xx, top - nr * rh - 0.18, 11)
      }
      for (j in 0:nr) {
        panel_line(wx, top - j * rh, wx + ww, top - j * rh, "#DDDDDD", 0.6)
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
              19
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
        panel_line(wx, yy, xe, yy, "black", 1.8)
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
    "./figures-output/FigureS15B.pdf",
    width = 15.05,
    height = 8.7,
    family = font
  )
  draw()
  dev.off()

  # FIGURE S15 (PANEL D) — Signed molecular loadings
  plot_signed_loadings(loads, "./figures-output/FigureS15D.pdf")

  # FIGURE S15 (PANEL E) — Cluster and Rutherford abundance contrasts
  plot_abundance_contrasts(selected_features, "./figures-output/FigureS15E.pdf", c(-10, 6))

  # FIGURE S15 (PANEL F) — Enriched protein pathways
  enrichment_input <- readRDS("./results/VIP-PAD/enrichment/F4.rds")
  plot_protein_pathways(
    enrichment_input, "./figures-output/FigureS15F.pdf",
    axis_padding = 0.65, bottom_margin = 18
  )

  # FIGURE S15 (PANEL C) — Metabolite loadings by lipid class
  plot_lipid_classes(
    enrichment_input, "./figures-output/FigureS15C.pdf", "Factor 4 weight"
  )
})



# ==============================================================================
# FIGURE S18 (PANELS A–F) — Molecular profile of VIP-PAD Factor 8
# ==============================================================================

local({
  # FIGURE S18 (PANEL A) — Highest absolute factor loadings
  loads <- factor_loadings(fit, "Factor8", metabolites)
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:15],
    loads[modality == "Metabolite" & relative_weight >= 0.75][order(
      -abs(loading_raw),
      feature
    )]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  ca <- copy(sel)
  ca[, feature_key := factor(feature, levels = feature)]
  ca[, label_plot := ifelse(modality == "Metabolite", wrap_label(label, 16), label)]
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
      y = unit(-1.55, "in"),
      width = unit(14.5, "in"),
      height = unit(14.5, "in"),
      just = c("left", "bottom")
    ))
    print(pa, newpage = FALSE)
    txt("Factor 8 top features", 0.48, 7.8, 32, rot = 90)
    panel_rect(0.65, 2.85, 0.32, 0.32, "#4575b4")
    txt("Proteins", 1.1, 2.69, 24, "left")
    panel_rect(0.65, 2.28, 0.32, 0.32, "#f46d43")
    txt("Metabolites", 1.1, 2.12, 24, "left")
    popViewport()
  }
  dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
  cairo_pdf(
    "./figures-output/FigureS18A.pdf",
    width = 14.5,
    height = 12.95,
    family = font
  )
  draw()
  dev.off()
  selected_features <- sel$feature

  # FIGURE S18 (PANEL B) — Molecular abundance by cluster
  sel <- rbindlist(list(
    loads[modality == "Protein"][order(-abs(loading_raw), feature)][1:8],
    loads[modality == "Metabolite"][order(-abs(loading_raw), feature)][1:2]
  ))
  sel[, direction := ifelse(loading_raw > 0, "Positive", "Negative")]
  sel[, vieworder := match(modality, c("Protein", "Metabolite"))]
  setorder(sel, vieworder, -relative_weight, feature)
  # Align molecular data and fixed cluster memberships by ID.
  ids <- rownames(fit$Z)
  ix <- match(ids, cohort$study_id)
  cx <- match(ids, cl$assignments$sample_id)
  groups <- cl$assignments$K2[cx]
  abundance <- as.matrix(cohort[ix, sel$feature, drop = FALSE])
  rownames(abundance) <- ids
  h <- list(
    groups = groups,
    patient_order = list(F08 = order(groups, fit$Z[, "Factor8"], ids)),
    molecular_matrix = abundance
  )
  font <- "Noto Sans"
  viewcols <- c(Protein = "#6196B7", Metabolite = "#E58A75")
  signcols <- c(Positive = "#B75C59", Negative = "#426F9B")
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
        fontface = if (bold) "bold" else "plain"
      )
    )
  }

  # Function to draw the cluster abundance heatmap.
  draw <- function() {
    grid.newpage()
    # Crop 1.65 inches from the left margin.
    pushViewport(viewport(
      x = unit(-1.65, "in"),
      y = 0,
      width = unit(15.7, "in"),
      height = unit(6.5, "in"),
      just = c("left", "bottom")
    ))
    ord <- h$patient_order$F08
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
      panel_line(x0 + gap, 5.89, x1 - gap, 5.89, "black", 1)
      panel_line(x0 + gap, 5.89, x0 + gap, 5.77, "black", 1)
      panel_line(x1 - gap, 5.89, x1 - gap, 5.77, "black", 1)
      txt(paste("Cluster", k), (x0 + x1) / 2, 6.12, 15)
    }
    for (v in c("Protein", "Metabolite")) {
      q <- sel[modality == v]
      nr <- nrow(q)
      top <- if (v == "Protein") {
        5.6
      } else {
        2.15
      }
      rh <- 0.34
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
      panel_line(left + boundary, top, left + boundary, top - nr * rh, "white", 1)
      for (vv in c(0, 0.3, 0.6, 0.9)) {
        xx <- wx + ww * vv / 1.1
        panel_line(xx, top, xx, top - nr * rh, "#C3C3C3", 0.8)
        txt(sprintf("%.1f", vv), xx, top - nr * rh - 0.18, 11)
      }
      for (j in 0:nr) {
        panel_line(wx, top - j * rh, wx + ww, top - j * rh, "#DDDDDD", 0.6)
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
              16
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
        panel_line(wx, yy, xe, yy, "black", 1.8)
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
    "./figures-output/FigureS18B.pdf",
    width = 14.05,
    height = 6.5,
    family = font
  )
  draw()
  dev.off()

  # FIGURE S18 (PANEL D) — Signed molecular loadings
  plot_signed_loadings(loads, "./figures-output/FigureS18D.pdf")

  # FIGURE S18 (PANEL E) — Cluster and Rutherford abundance contrasts
  plot_abundance_contrasts(selected_features, "./figures-output/FigureS18E.pdf", c(-3, 4))

  # FIGURE S18 (PANEL F) — Enriched protein pathways
  enrichment_input <- readRDS("./results/VIP-PAD/enrichment/F8.rds")
  plot_protein_pathways(
    enrichment_input, "./figures-output/FigureS18F.pdf",
    axis_padding = 0.65, bottom_margin = 18
  )

  # FIGURE S18 (PANEL C) — Metabolite loadings by lipid class
  plot_lipid_classes(
    enrichment_input, "./figures-output/FigureS18C.pdf", "Factor 8 weight"
  )
})
