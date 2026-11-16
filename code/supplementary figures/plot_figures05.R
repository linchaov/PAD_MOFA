# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3, ggrepel 0.9.8, patchwork 1.3.2, scales 1.4.0
#
# Inputs: Prepared muscle and plaque atlas embeddings, ranked genes, loadings and expression
#          summaries.
#
# Outputs: ./figures-output/FigureS05A.pdf through FigureS05C.pdf
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})



# ==============================================================================
# FIGURE S05 (PANELS A–B) — Tissue expression maps
# ==============================================================================
# Function to blend tissue colors with white.
soften <- function(palette, fraction) {
  setNames(
    grDevices::rgb(
      t(grDevices::col2rgb(palette) * (1 - fraction) + 255 * fraction),
      maxColorValue = 255
    ),
    names(palette)
  )
}
# Function to draw the shared cell map used in the tissue panels.
base_map <- function(cells, embedding) {
  xr <- range(cells$x)
  yr <- range(cells$y)
  dx <- diff(xr)
  dy <- diff(yr)
  ox <- xr[1] - 0.01 * dx
  oy <- yr[1] - 0.09 * dy
  al <- 0.13 * min(dx, dy)
  ggplot(cells, aes(x, y)) +
    coord_fixed(
      xlim = xr + c(-0.06, 0.06) * dx,
      ylim = yr + c(-0.15, 0.06) * dy,
      clip = "off"
    ) +
    annotate(
      "segment",
      x = ox,
      xend = ox + al,
      y = oy,
      yend = oy,
      linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
    ) +
    annotate(
      "segment",
      x = ox,
      xend = ox,
      y = oy,
      yend = oy + al,
      linewidth = 0.35,
      arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed")
    ) +
    annotate(
      "text",
      x = ox +
        al / 2,
      y = oy - 0.025 * dy,
      label = paste(embedding, "1"),
      family = "Arial",
      size = 2.7,
      vjust = 1
    ) +
    annotate(
      "text",
      x = ox - 0.025 * dx,
      y = oy + al / 2,
      label = paste(embedding, "2"),
      family = "Arial",
      size = 2.7,
      angle = 90,
      vjust = 0.5
    )
}



# ==============================================================================
# FIGURE S05 (PANEL C) — Expression by cell population
# ==============================================================================
# Function to plot gene expression by cell population.
make_dot <- function(key, title) {
  dt <- copy(expression_data[tissue == key])
  cts <- unique(as.character(dt$celltype))
  dt[, gene := factor(gene, levels = lev)]
  dt[, celltype := factor(celltype, levels = cts)]
  labs_ct <- c(
    `Endothelial cells` = "Endothelial",
    `Smooth muscle/Pericytes` = "SMC / pericytes",
    `Macrophages/Monocytes` = "Macro / mono",
    `Satellite cells` = "Satellite",
    `Mast cells` = "Mast",
    `NK/T/B cells` = "NK / T / B",
    Neutrophils = "Neutrophils",
    Myonuclei = "Myonuclei",
    FAPs = "FAPs"
  )
  # Function to format cell population labels for the horizontal axis.
  xlab <- function(x) {
    v <- labs_ct[x]
    v[is.na(v)] <- x[is.na(v)]
    unname(v)
  }
  ggplot(dt, aes(celltype, gene)) +
    geom_point(
      data = dt[display_dot == TRUE],
      aes(size = detection_pct, fill = z),
      shape = 21,
      stroke = 0.2,
      colour = "#555555"
    ) +
    geom_point(
      data = dt[present == FALSE],
      shape = 4,
      size = 1.7,
      colour = "grey65"
    ) +
    scale_fill_gradient2(
      low = "#2155A0",
      mid = "#FFFFFF",
      high = "#B92725",
      limits = c(-3, 3),
      breaks = c(-3, 0, 3),
      oob = scales::squish,
      name = "Relative RNA (Z-score)"
    ) +
    scale_size_area(
      max_size = 5.5,
      limits = c(0, 100),
      breaks = c(10, 25, 50, 100),
      name = "Detected cells (%)"
    ) +
    scale_y_discrete(
      limits = lev,
      drop = FALSE,
      expand = expansion(add = 0.5)
    ) +
    scale_x_discrete(labels = xlab, drop = FALSE) +
    labs(title = title) +
    common +
    theme(axis.text.y = element_blank(), legend.position = "bottom")
}



# ==============================================================================
# FIGURE S05 (PANELS A–B) — Prepare tissue maps
# ==============================================================================

theme_set(
  theme_void(base_family = "Arial", base_size = 10) +
    theme(
      text = element_text(colour = "black"),
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.margin = margin(7, 7, 7, 7),
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 10)
    )
)
major_pal <- c(
  Vascular = "#477DAA",
  Stromal = "#B88B48",
  Myogenic = "#8A6CAA",
  Immune = "#4B9A85",
  Lymphoid = "#6676B5",
  Myeloid = "#C47A68"
)
cell_pal <- c(
  `Endothelial cells` = "#4B83B4",
  `Smooth muscle/Pericytes` = "#8793C6",
  FAPs = "#D5AA57",
  `Satellite cells` = "#B084BE",
  Myonuclei = "#8265A7",
  `Macrophages/Monocytes` = "#D57965",
  `NK/T/B cells` = "#59A391",
  Neutrophils = "#B49D83",
  `Mast cells` = "#C77DAD",
  `T cell` = "#58A997",
  `NK cell` = "#5B84B9",
  `B cell` = "#B58ABD",
  Monocyte = "#E39A64",
  Macrophage = "#C76C63",
  `Dendritic cell` = "#C7B15A",
  Neutrophil = "#9EAEB1",
  `Mast cell` = "#C77DAD",
  `Plasma cell` = "#996EAA",
  EC = "#427EA1",
  `Smooth Muscle Cell` = "#858AC0",
  Fibroblast = "#B99A65",
  Fibromyocyte = "#927653"
)
major_pal_original <- major_pal
cell_pal_original <- cell_pal
short <- c(
  `Endothelial cells` = "Endothelial",
  `Smooth muscle/Pericytes` = "SMC / pericyte",
  FAPs = "FAPs",
  `Satellite cells` = "Satellite",
  Myonuclei = "Myonuclei",
  `Macrophages/Monocytes` = "Macrophage / monocyte",
  `NK/T/B cells` = "NK / T / B",
  Neutrophils = "Neutrophil",
  `Mast cells` = "Mast",
  `T cell` = "T cell",
  `NK cell` = "NK cell",
  `B cell` = "B cell",
  Monocyte = "Monocyte",
  Macrophage = "Macrophage",
  `Dendritic cell` = "Dendritic",
  Neutrophil = "Neutrophil",
  `Mast cell` = "Mast",
  `Plasma cell` = "Plasma cell",
  EC = "Endothelial",
  `Smooth Muscle Cell` = "SMC",
  Fibroblast = "Fibroblast",
  Fibromyocyte = "Fibromyocyte"
)



# ==============================================================================
# FIGURE S05 (PANELS A–B) — Draw and save tissue maps
# ==============================================================================

dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
for (tissue in c("Muscle", "Plaque")) {
  fraction <- if (tissue == "Muscle") {
    0.35
  } else {
    0.25
  }
  major_pal <- soften(major_pal_original, fraction)
  cell_pal <- soften(cell_pal_original, fraction)
  z <- readRDS(file.path(
    "./results/VIP-PAD/cell-inputs/Common_Atlas",
    paste0(tissue, "_Input.rds")
  ))
  cells <- copy(z$cells)
  if (tissue == "Muscle") {
    mapping <- c(
      `Endothelial cells` = "Vascular",
      `Smooth muscle/Pericytes` = "Vascular",
      FAPs = "Stromal",
      `Satellite cells` = "Myogenic",
      Myonuclei = "Myogenic",
      `Macrophages/Monocytes` = "Immune",
      `NK/T/B cells` = "Immune",
      Neutrophils = "Immune",
      `Mast cells` = "Immune"
    )
    ctorder <- c(
      "Endothelial cells",
      "Smooth muscle/Pericytes",
      "FAPs",
      "Satellite cells",
      "Myonuclei",
      "Macrophages/Monocytes",
      "NK/T/B cells",
      "Neutrophils",
      "Mast cells"
    )
  } else {
    mapping <- c(
      `T cell` = "Lymphoid",
      `NK cell` = "Lymphoid",
      `B cell` = "Lymphoid",
      `Plasma cell` = "Lymphoid",
      Monocyte = "Myeloid",
      Macrophage = "Myeloid",
      `Dendritic cell` = "Myeloid",
      Neutrophil = "Myeloid",
      `Mast cell` = "Myeloid",
      EC = "Vascular",
      `Smooth Muscle Cell` = "Vascular",
      Fibroblast = "Stromal",
      Fibromyocyte = "Stromal"
    )
    ctorder <- c(
      "T cell",
      "NK cell",
      "B cell",
      "Plasma cell",
      "Monocyte",
      "Macrophage",
      "Dendritic cell",
      "Neutrophil",
      "Mast cell",
      "EC",
      "Smooth Muscle Cell",
      "Fibroblast",
      "Fibromyocyte"
    )
  }
  cells[, compartment := unname(mapping[celltype])]
  cells[, celltype := factor(celltype, levels = ctorder)]
  set.seed(20261018)
  draw <- cells[sample(.N)]
  lab <- cells[, .(x = median(x), y = median(y), n = .N), by = celltype]
  lab[, label := unname(short[as.character(celltype)])]
  # Plaque includes small residual non-immune populations. Keep them in all panels.
  pB <- base_map(draw, z$embedding) +
    geom_point(aes(colour = celltype), size = 0.16, stroke = 0, alpha = 0.84) +
    scale_colour_manual(values = cell_pal, guide = "none") +
    geom_text_repel(
      data = lab,
      aes(label = label),
      colour = "black",
      family = "Arial",
      size = 3.3,
      seed = 42,
      box.padding = 0.4,
      point.padding = 0.2,
      min.segment.length = 0,
      segment.color = "#666666",
      segment.size = 0.25,
      max.overlaps = Inf,
      max.time = 5
    )
  title_text <- if (tissue == "Muscle") {
    "PAD skeletal muscle | GSE235143"
  } else {
    "Femoral plaque | GSE234077"
  }
  pB <- pB +
    labs(title = title_text) +
    theme(
      plot.title = element_text(
        family = "Arial",
        size = 13,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 10)
      )
    )
  panel <- if (tissue == "Muscle") {
    "FigureS05A.pdf"
  } else {
    "FigureS05B.pdf"
  }

  ggsave(
    file.path("./figures-output", panel),
    pB,
    width = 5.2,
    height = 4.8,
    device = cairo_pdf,
    family = "Arial",
    bg = "white"
  )
}



# ==============================================================================
# FIGURE S05 (PANEL C) — Prepare factor loadings and cell expression
# ==============================================================================

set.seed(42)
fs <- sprintf("F%02d", c(1:4, 8))
sel <- fread(
  "./results/VIP-PAD/cell-expression/Gene_selection.csv"
)
genes <- unique(sel$gene)
lev <- rev(genes)
ld <- fread(
  "./results/VIP-PAD/cell-expression/Loading_plot_data.csv"
)
ld[, gene := factor(gene, levels = lev)]
ld[, factor := factor(factor, levels = fs)]
expression_data <- fread(
  "./results/VIP-PAD/cell-expression/Cell_expression_plot_data.csv"
)
theme_set(theme_bw(base_size = 12, base_family = "Arial"))
common <- theme(
  panel.grid = element_blank(),
  axis.title = element_blank(),
  axis.ticks = element_blank(),
  axis.line = element_blank(),
  axis.text.y = element_text(size = 11, colour = "black"),
  axis.text.x = element_text(
    size = 14,
    colour = "black",
    angle = 45,
    hjust = 1
  ),
  panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.45),
  plot.title = element_text(size = 13, face = "bold"),
  legend.title = element_text(size = 12),
  legend.text = element_text(size = 11),
  plot.margin = margin(6, 7, 6, 7)
)
pl <- ggplot(ld, aes(factor, gene, fill = loading_per_score_SD)) +
  geom_tile(colour = "white", linewidth = 0.2) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    limits = c(-1, 1),
    breaks = c(-1, 0, 1),
    name = "Protein loading (per 1 SD)"
  ) +
  scale_y_discrete(limits = lev, drop = FALSE, expand = expansion(add = 0.5)) +
  scale_x_discrete(
    limits = fs,
    breaks = fs,
    labels = paste0("F", c(1:4, 8)),
    drop = FALSE
  ) +
  labs(title = "MOFA loadings") +
  common +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 1),
    legend.position = "bottom"
  )
pm <- make_dot("Muscle", "PAD muscle (4 samples)  GSE235143")
pp <- make_dot("Plaque_main", "PAD plaque (6 samples)  GSE234077")
pl <- pl +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(38, "mm")
    )
  )
pm <- pm +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(38, "mm")
    ),
    size = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
  )
pp <- pp +
  guides(
    fill = guide_colourbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = grid::unit(38, "mm")
    ),
    size = guide_legend(title.position = "top", title.hjust = 0.5, nrow = 1)
  )
fig <- pl +
  pm +
  pp +
  plot_layout(widths = c(2.6, 4.5, 4), guides = "collect") &
  theme(legend.position = "bottom")
h <- max(8, 0.23 * length(genes) + 2) + 0.35



# ==============================================================================
# FIGURE S05 (PANEL C) — Save expression plot
# ==============================================================================

ggsave(
  "./figures-output/FigureS05C.pdf",
  fig,
  width = 11.5,
  height = h,
  device = cairo_pdf,
  dpi = 300,
  bg = "white"
)
