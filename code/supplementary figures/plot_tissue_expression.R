# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3
#
# Global: S06, S07, S11, S14, S17, S21 and S24 show tissue expression of factor-associated genes.
#
# Inputs: Cohort-specific gene lists, tissue embeddings and cell-expression summaries.
# Outputs: Individual panel PDFs under ./figures-output/.
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
# FUNCTIONS
# ==============================================================================

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
# PLOT SETTINGS
# ==============================================================================


# Function to plot expression in skeletal muscle and femoral plaque.
plot_tissue_expression <- function(cohort, k, output_prefix) {
  mapping_dir <- file.path("./results", cohort, "mofa", sprintf("Factor_%d_Cell_Mapping", k))
  cell_dir <- file.path("./results", cohort, "cell-inputs", paste0("Factor_", k))
  S <- fread(
    file.path(mapping_dir, sprintf("F%d_Features.csv", k))
  )

  theme_set(
    theme_void(base_family = "Arial", base_size = 10) +
      theme(
        text = element_text(colour = "black"),
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        plot.margin = margin(7, 7, 7, 7),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8.5)
      )
  )
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

  dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)

  # FIGURES S06, S07, S11, S14, S17, S21 AND S24 (PANELS A, C) — Tissue expression density
  # Overlay joint expression density on each tissue embedding.
  for (tissue in c("Muscle", "Plaque")) {
    panel_title <- paste(tissue, "| Factor", k, "top features")
    z <- readRDS(file.path(
      cell_dir,
      paste0(tissue, "_Input.rds")
    ))
    cells <- as.data.table(z$cells)
    den <- fread(file.path(
      mapping_dir,
      paste0(tissue, "_Top15_Joint_Density.csv.gz")
    ))
    d <- merge(
      cells[, .(cell_id, x, y)],
      den[, .(cell_id, joint_density_relative)],
      by = "cell_id"
    )
    setorder(d, joint_density_relative)
    pp <- base_map(d, z$embedding) +
      geom_point(aes(colour = joint_density_relative), size = 0.19, stroke = 0) +
      scale_colour_gradientn(
        colours = c(
          "#D9D9DC",
          "#453276",
          "#345F8D",
          "#239A8C",
          "#73C46B",
          "#F1E448"
        ),
        values = c(0, 0.015, 0.1, 0.3, 0.6, 1),
        limits = c(0, 1),
        name = "Relative density",
        breaks = c(0, 0.25, 0.5, 0.75, 1)
      ) +
      annotate(
        "text",
        x = mean(range(cells$x)),
        y = max(cells$y) +
          0.045 *
            diff(range(cells$y)),
        label = panel_title,
        family = "Arial",
        size = 3.3
      ) +
      theme(legend.position = "bottom") +
      guides(
        colour = guide_colourbar(
          nbin = 1024,
          display = "raster",
          direction = "horizontal",
          barwidth = grid::unit(3.6, "cm"),
          barheight = grid::unit(0.23, "cm"),
          title.position = "top",
          ticks = FALSE
        )
      )

    # Panels A, C
    ggsave(
      c(
        Muscle = paste0(output_prefix, "A.pdf"),
        Plaque = paste0(output_prefix, "C.pdf")
      )[[tissue]],
      pp,
      width = 5.2,
      height = 4.8,
      device = cairo_pdf,
      family = "Arial",
      bg = "white"
    )
    # Summarize expression and detection by cell type.

    ########################################################################
    # FIGURES S06, S07, S11, S14, S17, S21 AND S24 (PANELS B, D) — Expression by cell population
    ########################################################################
    ds <- fread(file.path(
      mapping_dir,
      paste0(tissue, "_Dotplot_Values.csv")
    ))

    if (tissue == "Muscle") {
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
    ds[, celltype := factor(celltype, levels = rev(ctorder))]
    ds[, gene := factor(gene, levels = S$gene)]
    pD <- ggplot(ds, aes(gene, celltype)) +
      geom_blank() +
      geom_point(
        data = ds[display_dot == TRUE],
        aes(size = pct, fill = mean_expression),
        shape = 21,
        stroke = 0.18,
        colour = "#444444"
      ) +
      geom_point(
        data = ds[is.na(pct)],
        shape = 4,
        size = 1.2,
        colour = "#999999",
        stroke = 0.35
      ) +
      scale_size_area(
        max_size = 4.1,
        limits = c(0, 100),
        breaks = c(10, 25, 50, 75, 100),
        name = "Cells expressing (%)"
      ) +
      scale_fill_viridis_c(
        option = "D",
        begin = 0.04,
        end = 0.96,
        name = "Average expression",
        limits = c(0, ceiling(max(ds$mean_expression, na.rm = TRUE)))
      ) +
      scale_x_discrete(limits = function(x) S$gene[S$gene %in% x]) +
      scale_y_discrete(labels = short) +
      labs(x = NULL, y = NULL, title = panel_title) +
      theme_minimal(base_family = "Arial", base_size = 10) +
      theme(
        text = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        plot.title = element_text(
          size = 12,
          hjust = 0.5,
          margin = margin(b = 10)
        ),
        axis.text.x = element_text(angle = 55, hjust = 1, vjust = 1, size = 9),
        axis.text.y = element_text(size = 9),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_line(colour = "#ECECEF", linewidth = 0.3),
        strip.text = element_text(size = 9),
        legend.position = "bottom",
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        panel.spacing.x = grid::unit(5, "mm"),
        plot.margin = margin(8, 8, 8, 8)
      ) +
      guides(
        fill = guide_colourbar(
          order = 1,
          title.position = "top",
          barwidth = grid::unit(3.2, "cm"),
          barheight = grid::unit(0.25, "cm"),
          ticks = FALSE
        ),
        size = guide_legend(
          order = 2,
          title.position = "top",
          nrow = 1,
          override.aes = list(fill = "#666F7F")
        )
      )

    # Panels B, D
    ggsave(
      c(
        Muscle = paste0(output_prefix, "B.pdf"),
        Plaque = paste0(output_prefix, "D.pdf")
      )[[tissue]],
      pD,
      width = 8,
      height = if (tissue == "Muscle") {
        4.7
      } else {
        5.4
      },
      device = cairo_pdf,
      family = "Arial",
      bg = "white"
    )
  }
}



# ==============================================================================
# FIGURE S06 (PANELS A–D) — Tissue expression of genes encoding VIP-PAD Factor 1 proteins
# ==============================================================================

# FIGURE S06 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S06 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("VIP-PAD", 1, "./figures-output/FigureS06")



# ==============================================================================
# FIGURE S07 (PANELS A–D) — Tissue expression of genes encoding VIP-PAD Factor 2 proteins
# ==============================================================================

# FIGURE S07 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S07 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("VIP-PAD", 2, "./figures-output/FigureS07")



# ==============================================================================
# FIGURE S11 (PANELS A–D) — Tissue expression of genes encoding UKB Factor 5 proteins
# ==============================================================================

# FIGURE S11 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S11 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("UKB", 5, "./figures-output/FigureS11")



# ==============================================================================
# FIGURE S14 (PANELS A–D) — Tissue expression of genes encoding VIP-PAD Factor 3 proteins
# ==============================================================================

# FIGURE S14 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S14 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("VIP-PAD", 3, "./figures-output/FigureS14")



# ==============================================================================
# FIGURE S17 (PANELS A–D) — Tissue expression of genes encoding VIP-PAD Factor 4 proteins
# ==============================================================================

# FIGURE S17 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S17 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("VIP-PAD", 4, "./figures-output/FigureS17")



# ==============================================================================
# FIGURE S21 (PANELS A–D) — Tissue expression of genes encoding VIP-PAD Factor 8 proteins
# ==============================================================================

# FIGURE S21 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S21 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("VIP-PAD", 8, "./figures-output/FigureS21")



# ==============================================================================
# FIGURE S24 (PANELS A–D) — Tissue expression of genes encoding UKB Factor 2 proteins
# ==============================================================================

# FIGURE S24 — Muscle density (PANEL A) and cell populations (PANEL B)
# FIGURE S24 — Plaque density (PANEL C) and cell populations (PANEL D)
plot_tissue_expression("UKB", 2, "./figures-output/FigureS24")
