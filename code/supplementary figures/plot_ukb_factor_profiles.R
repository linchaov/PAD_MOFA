# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3, ggrepel 0.9.8, scales 1.4.0
#
# Global: S09 and S22 show the molecular profiles of UKB Factors 5 and 2.
#
# Inputs: UKB measurements, factor loadings, clusters and protein and NMR enrichment results.
# Outputs: Individual panel PDFs under ./figures-output/.
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
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to wrap labels to the available width.
wrap_label <- function(z, n = 26) {
  vapply(
    z,
    function(s) {
      s <- gsub('(?<=[-/])', ' ', s, perl = TRUE)
      paste(strwrap(s, width = n), collapse = '\n')
    },
    character(1)
  )
}

# Function to draw a line in panel coordinates.
panel_line <- function(x, y, x2, y2, col = '#C7C7C7', width = .6, lty = 1) {
  grid.lines(
    x = unit(c(x, x2), 'in'),
    y = unit(c(y, y2), 'in'),
    gp = gpar(col = col, lwd = width, lty = lty)
  )
}

# Function to draw a rectangle in panel coordinates.
panel_rect <- function(x, y, w, hh, fill, col = NA) {
  grid.rect(
    x = unit(x, 'in'),
    y = unit(y, 'in'),
    width = unit(w, 'in'),
    height = unit(hh, 'in'),
    just = c('left', 'top'),
    gp = gpar(fill = fill, col = col, lwd = .6)
  )
}

# Function to wrap labels to the available width.
wrap <- function(x, n = 24) {
  vapply(
    x,
    function(s) paste(strwrap(s, width = n), collapse = '\n'),
    character(1)
  )
}
# Function to format NMR trait labels for display.
pretty_nmr <- function(s) {
  s <- gsub('([a-z0-9])([A-Z])', '\\1 \\2', s, perl = TRUE)
  s <- gsub('([A-Z])([A-Z][a-z])', '\\1 \\2', s, perl = TRUE)
  s <- gsub('Metabolomics', '(NMR)', s, fixed = TRUE)
  s <- gsub('Conc ', 'Particle concentration: ', s, fixed = TRUE)
  s <- gsub('Very Small', 'XS', s, fixed = TRUE)
  s <- gsub('Very Large', 'XL', s, fixed = TRUE)
  s <- gsub('Large', 'L', s, fixed = TRUE)
  s <- gsub('Medium', 'M', s, fixed = TRUE)
  s <- gsub('Small', 'S', s, fixed = TRUE)
  s <- gsub('Polyunsaturated', 'PUFA', s, fixed = TRUE)
  s <- gsub('Monounsaturated', 'MUFA', s, fixed = TRUE)
  s <- gsub('Docosahexaenoic Acid', 'DHA', s, fixed = TRUE)
  s <- sub('^Pct (.*)$', '\\1 (%)', s)
  s <- gsub(' To Total Lipids In ', ' / total lipids in ', s, fixed = TRUE)
  s
}
# Function to join molecular loadings to the feature annotations for one factor.
loads_for <- function(k) {
  r <- rbindlist(lapply(names(fit$W), function(v) {
    w <- fit$W[[v]][, k]
    data.table(feature = names(w), modality = v, loading_raw = as.numeric(w))
  }))
  r[, label := ann$display_label[match(feature, ann$feature)]]
  r[,
    relative_weight := abs(loading_raw) / max(abs(loading_raw)),
    by = modality
  ]
  r[, direction := ifelse(loading_raw > 0, 'Positive', 'Negative')]
  r
}



# ==============================================================================
# FIGURES S09 AND S22 (PANELS A–E) — Shared data and plot settings
# ==============================================================================

cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
fit <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
nmr_enrichment <- data.table::fread(
  "./results/UKB/factor-profiles/UKB_F2_F5_NMR_Category_Enrichment.csv"
)
ann <- as.data.table(attr(cohort, "feature_dictionary"))
cl <- fread("./results/UKB/clustering/Patient_Assignments.csv")
font <- 'Noto Sans'
colsign <- c(Positive = '#C85B63', Negative = '#526EA4')
ann[,
  display_label := ifelse(
    modality == 'Protein',
    toupper(label),
    pretty_nmr(label)
  )
]
theme_main <- theme_classic(base_family = font, base_size = 15) +
  theme(
    axis.text = element_text(colour = 'black'),
    plot.margin = margin(10, 12, 10, 10),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11)
  )

# Function to plot top molecular loadings.
plot_top_loadings <- function(sel, k, output) {
  font <- 'Noto Sans'
  # Function to position text in panel coordinates.
  txt <- function(
    s,
    x,
    y,
    size = 12,
    just = 'centre',
    rot = 0,
    col = 'black',
    bold = FALSE
  ) {
    grid.text(
      s,
      x = unit(x, 'in'),
      y = unit(y, 'in'),
      just = just,
      rot = rot,
      gp = gpar(
        fontfamily = font,
        fontsize = size,
        col = col,
        fontface = if (bold) 'bold' else 'plain'
      )
    )
  }

  ca <- copy(sel)
  ca[, feature_key := factor(feature, levels = feature)]
  ca[, label_plot := ifelse(modality == 'Metabolite', wrap_label(label, 18), label)]
  ca[
    grepl('Hydroxyphenylacetic', label, ignore.case = TRUE),
    label_plot := 'ortho-\nHydroxyphenyl\nacetic acid'
  ]
  ca[
    grepl('Tetramethylchromano', label, ignore.case = TRUE),
    label_plot := 'Tetramethyl\nchromanol\nglucoside'
  ]
  ca[, value2 := pmax(0, relative_weight - 0.45 / .N)]
  ca[,
    number_y := relative_weight - .05 - ifelse(seq_len(.N) %% 2 == 0, .055, 0)
  ]
  pa <- ggplot(ca, aes(feature_key, relative_weight, fill = modality)) +
    geom_hline(
      yintercept = c(.25, .5, .75, 1),
      colour = 'black',
      linewidth = .45
    ) +
    geom_col(width = 1, alpha = .6) +
    geom_col(aes(y = value2), width = .85) +
    geom_text(
      aes(y = number_y, label = sprintf('%.2f', relative_weight)),
      size = 5.5,
      family = font
    ) +
    geom_text(
      aes(y = .50, label = ifelse(direction == 'Positive', '+', '\u2212')),
      size = 8,
      family = font
    ) +
    scale_fill_manual(
      values = c(Protein = '#4575b4', Metabolite = '#f46d43'),
      name = NULL
    ) +
    scale_x_discrete(labels = setNames(ca$label_plot, ca$feature)) +
    scale_y_continuous(limits = c(0, 1.02), expand = c(0, 0)) +
    coord_radial(inner.radius = .025, clip = 'off') +
    guides(theta = guide_axis_theta(angle = 90)) +
    theme_void(base_family = font) +
    theme(
      axis.text.theta = element_text(size = 23, colour = 'black'),
      legend.position = 'bottom',
      legend.text = element_text(size = 13),
      plot.margin = margin(100, 110, 110, 170)
    ) +
    labs(caption = 'Relative absolute weight') +
    theme(plot.caption = element_text(hjust = .5, size = 12))
  pa <- pa + theme(legend.position = 'none', plot.caption = element_blank())

  # Function to draw the ranked molecular loadings in radial coordinates.
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(
      x = 0,
      y = unit(0, 'in'),
      width = unit(14.5, 'in'),
      height = unit(14.5, 'in'),
      just = c('left', 'bottom')
    ))
    print(pa, newpage = FALSE)
    txt(paste('UKB Factor', k, 'top features'), .48, 7.8, 36, rot = 90)
    panel_rect(.65, 2.85, .32, .32, '#4575b4')
    txt('Proteins', 1.10, 2.69, 28, 'left')
    panel_rect(.65, 2.28, .32, .32, '#f46d43')
    txt('Metabolites', 1.10, 2.12, 28, 'left')
    popViewport()
  }

  # Panel A
  cairo_pdf(
    output,
    width = 14.5,
    height = 14.5,
    family = font
  )
  draw()
  dev.off()
}

# Function to plot abundance by molecular cluster.
plot_cluster_abundance <- function(sel, k, top_view, output) {
  ids <- rownames(fit$Z)
  ix <- match(ids, cohort$study_id)
  cx <- match(ids, cl$sample_id)
  groups <- cl$K2[cx]
  # Match features to measured columns and standardize each heatmap row.
  fields <- ann$source_field[match(sel$feature, ann$feature)]
  abundance <- as.matrix(cohort[ix, fields, drop = FALSE])
  dimnames(abundance) <- list(ids, sel$feature)
  h <- list(
    groups = groups,
    patient_order = order(groups, fit$Z[, k], ids),
    molecular_matrix = abundance
  )
  font <- 'Noto Sans'
  # Function to position text in panel coordinates.
  txt <- function(
    s,
    x,
    y,
    size = 12,
    just = 'centre',
    rot = 0,
    col = 'black',
    bold = FALSE
  ) {
    grid.text(
      s,
      x = unit(x, 'in'),
      y = unit(y, 'in'),
      just = just,
      rot = rot,
      gp = gpar(
        fontfamily = font,
        fontsize = size * 1.16,
        lineheight = .85,
        col = col,
        fontface = if (bold) 'bold' else 'plain'
      )
    )
  }

  # Function to draw the cluster abundance heatmap.
  draw <- function() {
    grid.newpage()
    pushViewport(viewport(
      x = unit(-.65, 'in'),
      y = 0,
      width = unit(15.7, 'in'),
      height = unit(8.7, 'in'),
      just = c('left', 'bottom')
    ))
    ord <- h$patient_order
    left <- 5.8
    hw <- 9.5
    g <- h$groups[ord]
    boundary <- sum(g == 1) / length(g) * hw
    rh <- .26
    pal <- colorRampPalette(c(
      '#4575B4',
      '#ABD9E9',
      '#FFFFBF',
      '#FDAE61',
      '#D73027'
    ))(257)
    # Function to map standardized abundance to heatmap colors.
    cols <- function(z) {
      idx <- floor((pmin(3, pmax(-3, z)) + 3) / 6 * 256) + 1
      {
        ans <- matrix(pal[idx], nrow = nrow(z), ncol = ncol(z))
        ans[is.na(ans)] <- '#888888'
        ans
      }
    }
    # Function to render the abundance matrix as a heatmap.
    raster <- function(mat, x, top, w, hh) {
      grid.raster(
        as.raster(mat),
        x = unit(x, 'in'),
        y = unit(top, 'in'),
        width = unit(w, 'in'),
        height = unit(hh, 'in'),
        just = c('left', 'top'),
        interpolate = FALSE
      )
    }
    # Brackets mark the two molecular clusters.
    for (k in 1:2) {
      x0 <- if (k == 1) left else left + boundary
      x1 <- if (k == 1) left + boundary else left + hw
      gap <- .035
      panel_line(x0 + gap, 8.09, x1 - gap, 8.09, 'black', 1)
      panel_line(x0 + gap, 8.09, x0 + gap, 7.97, 'black', 1)
      panel_line(x1 - gap, 8.09, x1 - gap, 7.97, 'black', 1)
      txt(paste('Cluster', k), (x0 + x1) / 2, 8.32, 15)
    }
    for (v in c('Metabolite', 'Protein')) {
      q <- sel[modality == v]
      nr <- nrow(q)
      top <- if (v == top_view) 7.80 else 3.20
      rh <- .48
      wx <- 4.12
      ww <- 1.28
      txt(
        if (v == 'Protein') 'Proteins' else 'Metabolites',
        wx + ww / 2,
        top + .29,
        15,
        bold = TRUE
      )
      z <- t(scale(h$molecular_matrix[, q$feature, drop = FALSE]))[,
        ord,
        drop = FALSE
      ]
      raster(cols(z), left, top, hw, nr * rh)
      panel_line(left + boundary, top, left + boundary, top - nr * rh, 'white', 1)
      for (vv in c(0, .3, .6, .9)) {
        xx <- wx + ww * vv / 1.1
        panel_line(xx, top, xx, top - nr * rh, '#C3C3C3', .8)
        txt(sprintf('%.1f', vv), xx, top - nr * rh - .18, 11)
      }
      for (j in 0:nr) {
        panel_line(wx, top - j * rh, wx + ww, top - j * rh, '#DDDDDD', .6)
      }
      panel_rect(wx, top, ww, nr * rh, NA, '#444444')
      panel_rect(left, top, hw, nr * rh, NA, '#999999')
      for (i in seq_len(nr)) {
        yy <- top - (i - .5) * rh
        txt(
          wrap_label(q$label[i], if (v == 'Protein') 37 else 24),
          3.94,
          yy,
          if (v == 'Protein') 12.5 else 11.2,
          'right'
        )
        xe <- wx + ww * q$relative_weight[i] / 1.1
        panel_line(wx, yy, xe, yy, 'black', 1.8)
        grid.points(
          x = unit(xe, 'in'),
          y = unit(yy, 'in'),
          pch = 16,
          size = unit(.10, 'in'),
          gp = gpar(col = 'black')
        )
        txt(
          if (q$direction[i] == 'Positive') '+' else '\u2212',
          5.59,
          yy,
          22,
          col = 'black'
        )
      }
    }
    raster(matrix(pal, nrow = 1), 6.25, .44, 2, .12)
    txt('Z score', 7.25, .57, 11)
    for (j in 0:6) {
      txt(as.character(j - 3), 6.25 + j / 3, .16, 10)
    }
    txt('+  Positive weight     \u2212  Negative weight', 11, .38, 12)

    popViewport()
  }

  # Panel B
  cairo_pdf(
    output,
    width = 15.05,
    height = 8.7,
    family = font
  )
  draw()
  dev.off()
}

# Function to plot NMR loadings by category.
plot_nmr_categories <- function(k, input, output) {
  classes <- c("VLDL", "IDL", "LDL", "HDL", "Fatty acids")
  d <- fread(
    input
  )
  d[, category := factor(category, levels = c(classes, "Other"))]
  q <- copy(nmr_enrichment[factor == k])
  q[, category := factor(pathway, levels = c(classes, 'Other'))]
  span <- diff(range(d$loading_raw))
  q[, y := max(d$loading_raw) + span * .09]
  cp <- c(
    VLDL = '#7295B3',
    IDL = '#D19B5E',
    LDL = '#B4A15C',
    HDL = '#87AA98',
    'Fatty acids' = '#A68EB6',
    Other = '#B6B6B6'
  )
  p <- ggplot(d, aes(category, loading_raw, colour = category)) +
    geom_hline(yintercept = 0, linetype = 2, colour = 'grey65') +
    geom_point(
      position = position_jitter(width = .15, height = 0, seed = 51),
      size = 1.8,
      alpha = .8
    ) +
    geom_text(
      data = q[stars != ''],
      aes(x = category, y = y, label = stars),
      inherit.aes = FALSE,
      colour = 'black',
      size = 6.4
    ) +
    scale_colour_manual(values = cp) +
    scale_y_continuous(
      limits = c(
        min(d$loading_raw) - span * .04,
        max(d$loading_raw) + span * .17
      )
    ) +
    labs(x = NULL, y = paste('Factor', k, 'weight')) +
    theme_main +
    theme(legend.position = 'none', axis.text.x = element_text(size = 12))
  p <- p +
    theme(
      axis.text = element_text(size = 16),
      axis.text.x = element_text(size = 16),
      axis.title.y = element_text(size = 19)
    )
  p <- p +
    theme(
      axis.text = element_text(size = 19),
      axis.text.x = element_text(size = 19),
      axis.title.y = element_text(size = 23)
    )

  # Panel C
  ggsave(
    output,
    p,
    width = 8,
    height = 5,
    device = cairo_pdf
  )
}

# Function to plot enriched UKB protein pathways.
plot_pathways <- function(input, output, stack_guides) {
  q <- fread(
    input
  )
  q[,
    collection := factor(
      collection,
      levels = c("GO:BP", "KEGG", "Reactome", "REACTOME")
    )
  ]
  setorder(q, collection, -ratio, padj)
  q[, display := factor(pathway, levels = rev(unique(pathway)))]
  lab <- setNames(wrap(q$display_label, 45), q$pathway)
  q[, logFDR := -log10(pmax(padj, .Machine$double.xmin))]
  p <- ggplot(q, aes(ratio, display)) +
    geom_point(
      aes(size = size, fill = logFDR, shape = direction),
      colour = '#333333',
      stroke = .3,
      position = position_dodge(width = .35, orientation = 'y')
    ) +
    facet_grid(collection ~ ., scales = 'free_y', space = 'free_y') +
    scale_y_discrete(labels = lab, expand = expansion(add = .7)) +
    scale_x_continuous(
      limits = c(0, max(.02, ceiling(max(q$ratio) * 100) / 100)),
      breaks = scales::breaks_pretty(n = 3),
      labels = function(x) sprintf('%.2f', x),
      expand = expansion(mult = c(.03, .08))
    ) +
    scale_size_area(
      max_size = 6,
      name = 'Count',
      breaks = if (max(q$size) > 300) {
        c(20, 100, 500)
      } else if (max(q$size) > 100) {
        c(20, 100, 200)
      } else {
        c(10, 20, 50)
      }
    ) +
    scale_shape_manual(
      values = c(negative = 25, positive = 24, all = 21),
      labels = c(
        negative = 'Negative\nweights',
        positive = 'Positive\nweights',
        all = 'Absolute\nweights'
      ),
      name = NULL
    ) +
    scale_fill_gradient(
      low = '#F6C1C7',
      high = '#8D0A15',
      breaks = scales::breaks_pretty(n = 3),
      name = expression(-log[10](FDR))
    ) +
    labs(x = 'Gene ratio', y = NULL) +
    theme_bw(base_size = 13, base_family = font) +
    theme(
      axis.text = element_text(colour = 'black'),
      panel.grid.minor = element_blank(),
      legend.position = 'bottom',
      legend.box = 'horizontal',
      legend.spacing.x = unit(2, 'mm'),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 11),
      legend.key.width = unit(4, 'mm'),
      strip.text = element_text(size = 12),
      strip.text.y = element_text(angle = 0, size = 11),
      plot.margin = margin(8, 10, 8, 8)
    ) +
    guides(
      shape = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(size = 3, fill = 'grey80')
      ),
      size = guide_legend(
        order = 2,
        nrow = 1,
        title.position = 'top',
        override.aes = list(fill = 'white', shape = 21)
      ),
      fill = guide_colourbar(
        order = 3,
        direction = 'horizontal',
        title.position = 'top',
        barwidth = unit(25, 'mm'),
        barheight = unit(3, 'mm')
      )
    )
  p <- p +
    theme(
      axis.text = element_text(size = 15),
      axis.title = element_text(size = 16),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 13),
      strip.text.y = element_text(size = 13)
    )
  if (stack_guides) {
    p <- p +
      theme(
        legend.position = 'right',
        legend.box = 'vertical',
        legend.spacing.y = unit(3, 'mm')
      ) +
      guides(
        shape = guide_legend(
          order = 1,
          ncol = 1,
          override.aes = list(size = 3, fill = 'grey80')
        ),
        size = guide_legend(
          order = 2,
          ncol = 1,
          title.position = 'top',
          override.aes = list(fill = 'white', shape = 21)
        ),
        fill = guide_colourbar(
          order = 3,
          direction = 'vertical',
          title.position = 'top',
          barwidth = unit(4, 'mm'),
          barheight = unit(28, 'mm')
        )
      )
  }
  p <- p + guides(shape = 'none')
  lines <- sum(vapply(strsplit(unname(lab), '\n'), length, integer(1)))
  ggsave(
    output,
    p,
    width = 9.2,
    height = max(6, 2 + .47 * nrow(q) + .23 * (lines - nrow(q))),
    device = cairo_pdf
  )
}


dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)



# ==============================================================================
# FIGURE S09 (PANELS A–E) — Molecular profile of UKB Factor 5
# ==============================================================================

local({
  k <- 5

  # FIGURE S09 (PANEL A) — Top molecular loadings
  loads <- loads_for(k)
  sel <- rbind(
    loads[modality == 'Protein'][order(-relative_weight, feature)][1:15],
    loads[modality == 'Metabolite' & relative_weight > .8][order(
      -relative_weight,
      feature
    )]
  )
  plot_top_loadings(sel, k, "./figures-output/FigureS09A.pdf")

  # FIGURE S09 (PANEL B) — Abundance by molecular cluster
  sel <- rbind(
    loads[modality == 'Protein'][order(-relative_weight, feature)][1:8],
    loads[modality == 'Metabolite'][order(-relative_weight, feature)][1:4]
  )
  plot_cluster_abundance(sel, k, "Protein", "./figures-output/FigureS09B.pdf")

  # FIGURE S09 (PANEL C) — NMR loadings by category
  plot_nmr_categories(k, "./results/UKB/factor-profiles/F5_NMR_Categories.csv", "./figures-output/FigureS09C.pdf")

  # FIGURE S09 (PANEL D) — Enriched protein pathways
  plot_pathways(
    "./results/UKB/factor-profiles/F5_Absolute_Displayed_Pathways.csv", "./figures-output/FigureS09D.pdf",
    stack_guides = TRUE
  )

  # FIGURE S09 (PANEL E) — Scaled protein and metabolite loadings
  local({
    d <- copy(loads)
    # Cross-view plot: same global scaling and signed-power display as VIP-PAD.
    d[, scaled_weight := loading_raw / max(abs(loading_raw))]
    d[,
      view := factor(
        ifelse(modality == 'Protein', 'Proteomics', 'Metabolomics'),
        levels = c('Proteomics', 'Metabolomics')
      )
    ]
    ids <- d[
      order(-abs(loading_raw), feature),
      head(.SD, 5),
      by = .(view, direction)
    ]$feature
    d[, selected := feature %in% ids]
    set.seed(51)
    d[, position := runif(.N, -.18, .18)]
    d[, display_label := wrap(label, 21)]
    tr <- scales::trans_new(
      'signed_power_1.8',
      function(z) sign(z) * abs(z)^1.8,
      function(z) sign(z) * abs(z)^(1 / 1.8)
    )
    labd <- copy(d[selected == TRUE])
    setorder(labd, view, -scaled_weight)
    labd[,
      label_y := sign(scaled_weight) * seq(.94, .13, length.out = .N)^(1 / 1.8),
      by = .(view, direction)
    ]
    labd[
      direction == 'Negative',
      label_y := -seq(.13, .94, length.out = .N)^(1 / 1.8),
      by = view
    ]
    labd[, label_x := rep(c(-.42, .42), length.out = .N), by = .(view, direction)]
    labd[
      modality == 'Metabolite' & abs(label_y) < .2^(1 / 1.8),
      label_y := sign(label_y) * .2^(1 / 1.8)
    ]
    labd[
      modality == 'Metabolite' &
        direction == 'Negative' &
        abs(label_y) <= .2^(1 / 1.8),
      label_y := -.14^(1 / 1.8)
    ]
    p <- ggplot(d, aes(position, scaled_weight)) +
      geom_hline(
        yintercept = 0,
        colour = 'grey65',
        linetype = 2,
        linewidth = .35
      ) +
      geom_point(
        data = d[selected == FALSE],
        colour = '#929292',
        alpha = .23,
        size = 1.4
      ) +
      geom_point(
        data = d[selected == TRUE],
        aes(colour = direction),
        size = 2.6
      ) +
      geom_text_repel(
        data = d[selected == TRUE & modality == 'Protein'],
        aes(label = display_label, colour = direction),
        size = 5.2,
        family = font,
        seed = 51,
        max.overlaps = Inf,
        box.padding = .65,
        force = 5,
        min.segment.length = 0,
        segment.colour = 'grey60',
        max.iter = 100000,
        max.time = 5
      ) +
      geom_segment(
        data = labd[modality == 'Metabolite'],
        aes(x = position, y = scaled_weight, xend = label_x, yend = label_y),
        colour = 'grey65',
        linewidth = .25
      ) +
      geom_label(
        data = labd[modality == 'Metabolite'],
        aes(x = label_x, y = label_y, label = display_label, colour = direction),
        size = 4.2,
        family = font,
        lineheight = .92,
        fill = 'white',
        linewidth = 0,
        label.padding = unit(.05, 'lines')
      ) +
      facet_wrap(~view, nrow = 1) +
      scale_x_continuous(limits = c(-.75, .75), breaks = NULL) +
      scale_y_continuous(
        transform = tr,
        limits = c(-1.06, 1.06),
        breaks = seq(-1, 1, .2)
      ) +
      scale_colour_manual(values = colsign) +
      labs(x = 'View', y = 'Scaled weight') +
      theme_main +
      theme(
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        panel.border = element_rect(fill = NA, colour = 'black', linewidth = .5),
        axis.line = element_blank()
      )
    p <- p +
      coord_flip() +
      facet_wrap(~view, ncol = 1) +
      labs(x = NULL) +
      scale_y_continuous(
        transform = tr,
        limits = c(-1.06, 1.06),
        breaks = c(-1, -.8, -.6, -.4, 0, .4, .6, .8, 1)
      )
    p <- p +
      theme(
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 17)
      )
    p <- p +
      theme(
        strip.text = element_text(size = 19),
        axis.text = element_text(size = 15),
        axis.title = element_text(size = 18)
      )

    # Panel E
    ggsave(
      "./figures-output/FigureS09E.pdf",
      p,
      width = 9.6,
      height = 5.8,
      device = cairo_pdf
    )
  })
})



# ==============================================================================
# FIGURE S22 (PANELS A–E) — Molecular profile of UKB Factor 2
# ==============================================================================

local({
  k <- 2

  # FIGURE S22 (PANEL A) — Top molecular loadings
  loads <- loads_for(k)
  sel <- rbind(
    loads[modality == 'Protein' & relative_weight > .8][order(
      -relative_weight,
      feature
    )],
    loads[modality == 'Metabolite'][order(-relative_weight, feature)][1:15]
  )
  plot_top_loadings(sel, k, "./figures-output/FigureS22A.pdf")

  # FIGURE S22 (PANEL B) — Abundance by molecular cluster
  sel <- rbind(
    loads[modality == 'Protein'][order(-relative_weight, feature)][1:4],
    loads[modality == 'Metabolite'][order(-relative_weight, feature)][1:8]
  )
  plot_cluster_abundance(sel, k, "Metabolite", "./figures-output/FigureS22B.pdf")

  # FIGURE S22 (PANEL C) — NMR loadings by category
  plot_nmr_categories(k, "./results/UKB/factor-profiles/F2_NMR_Categories.csv", "./figures-output/FigureS22C.pdf")

  # FIGURE S22 (PANEL D) — Scaled protein and metabolite loadings
  local({
    d <- copy(loads)
    # Cross-view plot: same global scaling and signed-power display as VIP-PAD.
    d[, scaled_weight := loading_raw / max(abs(loading_raw))]
    d[,
      view := factor(
        ifelse(modality == 'Protein', 'Proteomics', 'Metabolomics'),
        levels = c('Proteomics', 'Metabolomics')
      )
    ]
    ids <- d[
      order(-abs(loading_raw), feature),
      head(.SD, 5),
      by = .(view, direction)
    ]$feature
    d[, selected := feature %in% ids]
    set.seed(51)
    d[, position := runif(.N, -.18, .18)]
    d[, display_label := wrap(label, 21)]
    d[
      display_label == 'Phospholipids In IDL',
      display_label := 'Phospholipids\nIn IDL'
    ]
    tr <- scales::trans_new(
      'signed_power_1.8',
      function(z) sign(z) * abs(z)^1.8,
      function(z) sign(z) * abs(z)^(1 / 1.8)
    )
    labd <- copy(d[selected == TRUE])
    setorder(labd, view, -scaled_weight)
    labd[,
      label_y := sign(scaled_weight) * seq(.94, .13, length.out = .N)^(1 / 1.8),
      by = .(view, direction)
    ]
    labd[
      direction == 'Negative',
      label_y := -seq(.13, .94, length.out = .N)^(1 / 1.8),
      by = view
    ]
    labd[, label_x := rep(c(-.42, .42), length.out = .N), by = .(view, direction)]
    labd[
      modality == 'Metabolite' & abs(label_y) < .2^(1 / 1.8),
      label_y := sign(label_y) * .2^(1 / 1.8)
    ]
    p <- ggplot(d, aes(position, scaled_weight)) +
      geom_hline(
        yintercept = 0,
        colour = 'grey65',
        linetype = 2,
        linewidth = .35
      ) +
      geom_point(
        data = d[selected == FALSE],
        colour = '#929292',
        alpha = .23,
        size = 1.4
      ) +
      geom_point(
        data = d[selected == TRUE],
        aes(colour = direction),
        size = 2.6
      ) +
      geom_text_repel(
        data = d[selected == TRUE & modality == 'Protein'],
        aes(label = display_label, colour = direction),
        size = 4.8,
        family = font,
        seed = 51,
        max.overlaps = Inf,
        box.padding = .65,
        force = 5,
        min.segment.length = 0,
        segment.colour = 'grey60',
        max.iter = 100000,
        max.time = 5
      ) +
      geom_segment(
        data = labd[modality == 'Metabolite'],
        aes(x = position, y = scaled_weight, xend = label_x, yend = label_y),
        colour = 'grey65',
        linewidth = .25
      ) +
      geom_label(
        data = labd[modality == 'Metabolite'],
        aes(x = label_x, y = label_y, label = display_label, colour = direction),
        size = 3.9,
        family = font,
        lineheight = .92,
        fill = 'white',
        linewidth = 0,
        label.padding = unit(.05, 'lines')
      ) +
      facet_wrap(~view, nrow = 1) +
      scale_x_continuous(limits = c(-.75, .75), breaks = NULL) +
      scale_y_continuous(
        transform = tr,
        limits = c(-1.06, 1.06),
        breaks = seq(-1, 1, .2)
      ) +
      scale_colour_manual(values = colsign) +
      labs(x = 'View', y = 'Scaled weight') +
      theme_main +
      theme(
        legend.position = 'none',
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        panel.border = element_rect(fill = NA, colour = 'black', linewidth = .5),
        axis.line = element_blank()
      )
    p <- p +
      coord_flip() +
      facet_wrap(~view, ncol = 1) +
      labs(x = NULL) +
      scale_y_continuous(
        transform = tr,
        limits = c(-1.06, 1.06),
        breaks = c(-1, -.8, -.6, -.4, 0, .4, .6, .8, 1)
      )
    p <- p +
      theme(
        strip.text = element_text(size = 18),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 17)
      )

    # Panel D
    ggsave(
      "./figures-output/FigureS22D.pdf",
      p,
      width = 10,
      height = 5.8,
      device = cairo_pdf
    )
  })

  # FIGURE S22 (PANEL E) — Enriched protein pathways
  plot_pathways(
    "./results/UKB/factor-profiles/F2_Absolute_Displayed_Pathways.csv", "./figures-output/FigureS22E.pdf",
    stack_guides = FALSE
  )
})
