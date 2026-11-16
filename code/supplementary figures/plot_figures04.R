# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, ggplot2 4.0.3, gtable 0.3.6
#
# Inputs: VIP-PAD factor-comorbidity model estimates in Figure_3_Plot_Data.rds.
#
# Outputs: ./figures-output/FigureS04.pdf
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
# FIGURE S04 — Prepare comorbidity associations
# ==============================================================================

x <- readRDS("./results/VIP-PAD/clinical/Figure_3_Plot_Data.rds")
d <- copy(x$statistics)
d[, disease := factor(id, levels = rev(x$traits$id))]
d[,
  panel := factor(
    paste("Factor", as.integer(sub("F", "", factor))),
    levels = paste("Factor", 1:10)
  )
]
d[,
  significance := factor(
    ifelse(P < 0.05, "P < 0.05", "P >= 0.05"),
    levels = c("P < 0.05", "P >= 0.05")
  )
]
labels <- setNames(x$traits$label, x$traits$id)
bands <- data.table(y = seq(2, 10, 2))
# 02. Draw estimates and confidence intervals.
p <- ggplot(d, aes(x = OR, y = disease, colour = significance)) +
  geom_vline(
    xintercept = 1,
    colour = "#B7BCC2",
    linetype = "dashed",
    linewidth = 0.35
  ) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    orientation = "y",
    width = 0.2,
    linewidth = 0.65
  ) +
  geom_point(size = 2.2) +
  facet_wrap(~panel, ncol = 5, axes = "all_x", axis.labels = "all_x") +
  scale_x_log10(
    limits = c(
      min(0.4, min(d$CI_low) * 0.97),
      max(
        3.2,
        max(d$CI_high) *
          1.03
      )
    ),
    breaks = c(0.5, 1, 2, 3),
    labels = c("0.5", "1", "2", "3")
  ) +
  scale_y_discrete(labels = labels) +
  scale_colour_manual(
    values = c(`P < 0.05` = "#000000", `P >= 0.05` = "#B3B3B3"),
    name = NULL
  ) +
  labs(x = "Adjusted OR (95% CI)", y = NULL) +
  theme_classic(base_family = "Noto Sans", base_size = 12) +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black", size = 11),
    axis.title = element_text(size = 12),
    axis.line.x = element_line(linewidth = 0.55),
    axis.ticks.x = element_line(colour = "black", linewidth = 0.55),
    axis.ticks.length = grid::unit(2, "mm"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 14),
    panel.spacing.x = grid::unit(0.55, "cm"),
    panel.spacing.y = grid::unit(0.9, "cm"),
    panel.background = element_blank(),
    legend.position = "none",
    plot.margin = margin(10, 12, 10, 10)
  )
# 03. Continuous row bands span panel gaps; axes remain separate.
# Construct the layout with the export font metrics to prevent label clipping.
layout_device <- tempfile(fileext = ".pdf")



# ==============================================================================
# DRAW AND SAVE FIGURE S04
# ==============================================================================

dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
grDevices::cairo_pdf(
  layout_device,
  width = 12,
  height = 10.5,
  family = "Noto Sans"
)
build <- ggplot_build(p)
y_range <- build$layout$panel_params[[1]]$y.range
grob <- ggplotGrob(p)
panels <- grob$layout[grepl("^panel", grob$layout$name), ]
for (top in unique(panels$t)) {
  cells <- panels[panels$t == top, ]
  band_grob <- grid::rectGrob(
    x = 0.5,
    y = (bands$y - y_range[1]) / diff(y_range),
    width = 1,
    height = 1 / diff(y_range),
    default.units = "npc",
    gp = grid::gpar(fill = "#F3F5F7", col = NA)
  )
  grob <- gtable::gtable_add_grob(
    grob,
    band_grob,
    t = top,
    b = max(cells$b),
    l = min(cells$l),
    r = max(cells$r),
    z = 0,
    clip = "on",
    name = paste0("row_bands_", top)
  )
}
invisible(dev.off())
unlink(layout_device)
ggsave(
  "./figures-output/FigureS04.pdf",
  grob,
  width = 12,
  height = 10.5,
  device = cairo_pdf,
  dpi = 180,
  bg = "white"
)
