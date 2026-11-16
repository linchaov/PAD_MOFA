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
# Inputs: UKB clinical and renal associations in clinical/Factor_Associations.rds.
#
# Outputs: ./figures-output/FigureS08A.pdf through FigureS08B.pdf
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
# FIGURE S08 (PANEL A) — Factor–trait correlations
# ==============================================================================
# Function to draw the factor–trait heatmap.
draw_heat <- function() {
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(xscale = c(0, 11.15), yscale = c(0, 6.15)))
  # Function to position text in panel coordinates.
  txt <- function(
    t,
    x,
    y,
    size = 10,
    just = "centre",
    rot = 0,
    bold = FALSE,
    col = "black"
  ) {
    grid::grid.text(
      t,
      grid::unit(x, "native"),
      grid::unit(y, "native"),
      just = just,
      rot = rot,
      gp = grid::gpar(
        fontfamily = "Arial",
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
  # Function to draw a line in panel coordinates.
  line <- function(x0, y0, x1, y1, col = "#7F8A94", lwd = 0.6) {
    grid::grid.lines(
      grid::unit(c(x0, x1), "native"),
      grid::unit(c(y0, y1), "native"),
      gp = grid::gpar(col = col, lwd = lwd)
    )
  }
  # Function to draw a rectangle in panel coordinates.
  rect <- function(x, y, w, h, fill) {
    grid::grid.rect(
      grid::unit(x, "native"),
      grid::unit(y, "native"),
      grid::unit(w, "native"),
      grid::unit(h, "native"),
      gp = grid::gpar(fill = fill, col = NA)
    )
  }
  pal <- colorRampPalette(c("darkblue", "white", "darkred"))(4001)
  cw <- 0.37
  left <- 1.02
  gap <- 0.16
  xl <- left +
    (seq_len(nrow(v)) - 1) * cw +
    (match(v$group, unique(v$group)) - 1) * gap
  yy <- 4.3 - (0:10) * 0.34
  for (k in 1:11) {
    txt(paste("Factor", k), 0.9, yy[k], 10, "right")
  }
  for (g in unique(v$group)) {
    ii <- which(v$group == g)
    lo <- xl[min(ii)]
    hi <- xl[max(ii)] + cw
    txt(g, (lo + hi) / 2, 5.98, 11, bold = TRUE)
    line(lo, 4.54, hi, 4.54)
    for (j in 0:length(ii)) {
      line(lo + j * cw, 4.54, lo + j * cw, 4.485, lwd = 0.5)
    }
  }
  for (j in seq_len(nrow(v))) {
    x <- xl[j] + cw / 2
    txt(v$label[j], x, 4.62, 10, "left", 55)
    for (k in 1:11) {
      z <- r[factor_order == k & variable_order == j]
      diameter <- 3 + 5 * sqrt(pmin(6, -log10(pmax(z$FDR, 1e-300))))
      colour <- pal[1 + round((z$r + 1) / 2 * 4000)]
      grid::grid.circle(
        grid::unit(x, "native"),
        grid::unit(yy[k], "native"),
        r = grid::unit(diameter / 2, "pt"),
        gp = grid::gpar(fill = colour, col = NA)
      )
      if (z$stars != "") {
        txt(
          z$stars,
          x,
          yy[k] + 0.005,
          7.4,
          col = if (abs(z$r) > 0.65) {
            "white"
          } else {
            "black"
          }
        )
      }
    }
  }
  lx <- left
  lw <- 1.35
  txt("Correlation", lx + lw / 2, 0.59, 9)
  for (i in 1:200) {
    rect(
      lx + (i - 0.5) * lw / 200,
      0.39,
      lw / 200 + 2e-04,
      0.1,
      pal[1 + round((i - 1) / 199 * 4000)]
    )
  }
  for (tick in c(-1, -0.5, 0, 0.5, 1)) {
    x <- lx + (tick + 1) / 2 * lw
    line(x, 0.34, x, 0.29, lwd = 0.4)
    txt(as.character(tick), x, 0.17, 9)
  }
  grid::popViewport()
}

# Function to format a value in scientific notation.
fmt <- function(x) formatC(x, format = "e", digits = 2)



# ==============================================================================
# FIGURE S08 (PANEL B) — Factor scores and clinical measurements
# ==============================================================================
# Function to plot the clinical measurement against factor scores and annotate their correlation.
scatter <- function(f, vname, ylab, colour) {
  z <- data.frame(score = d[[f]], value = d[[vname]])
  z <- z[is.finite(z$score) & is.finite(z$value), ]
  index <- which(r[["factor"]] == f & r[["variable"]] == vname)
  s <- r[index]
  # The line shows the unadjusted linear trend.
  ggplot(z, aes(score, value)) +
    geom_point(colour = colour, alpha = 0.58, size = 1.45) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = FALSE,
      colour = colour,
      linewidth = 0.65
    ) +
    labs(
      subtitle = sprintf("N = %d   r = %.3f   FDR = %s", s$N, s$r, fmt(s$FDR)),
      x = paste("UKB Factor", sub("F", "", f), "score"),
      y = ylab
    ) +
    base +
    theme(
      axis.title = element_text(size = 13),
      axis.text = element_text(size = 12),
      plot.subtitle = element_text(size = 12)
    )
}



# ==============================================================================
# FIGURE S08 (PANELS A–B) — Prepare clinical associations
# ==============================================================================

a <- readRDS(
  "./results/UKB/clinical/Factor_Associations.rds"
)
r <- copy(a$associations)
v <- a$variables
d <- a$data
base <- theme_classic(base_size = 11, base_family = "Arial") +
  theme(
    text = element_text(colour = "black"),
    axis.text = element_text(colour = "black"),
    plot.subtitle = element_text(size = 12),
    plot.margin = margin(8, 10, 8, 8)
  )
b <- scatter("F5", "eGFR_calculated", "eGFR (mL/min/1.73 m²)", "#597DA1")
c <- scatter("F5", "CystatinC_i0_Blood", "Cystatin C (mg/L)", "#BE7783")



# ==============================================================================
# FIGURE S08 (PANELS A–B) — Draw and save panels
# ==============================================================================

# All 275 factor–phenotype tests share one BH-FDR family.
dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)
# FIGURE S08 (PANEL A) — Factor–trait correlations
cairo_pdf(
  "./figures-output/FigureS08A.pdf",
  width = 11.15,
  height = 6.15,
  family = "Arial"
)
draw_heat()
invisible(dev.off())
# FIGURE S08 (PANEL B) — Factor 5 and kidney markers
cairo_pdf(
  "./figures-output/FigureS08B.pdf",
  width = 8.6,
  height = 3.6,
  family = "Arial"
)
gb <- ggplotGrob(b)
gc <- ggplotGrob(c)
gb$widths <- gc$widths <- grid::unit.pmax(gb$widths, gc$widths)
gb$heights <- gc$heights <- grid::unit.pmax(gb$heights, gc$heights)
for (j in 1:2) {
  grid::pushViewport(grid::viewport(
    x = (j - 0.5) / 2,
    y = 0.5,
    width = 0.5,
    height = 1
  ))
  grid::grid.draw(
    if (j == 1) {
      gb
    } else {
      gc
    }
  )
  grid::popViewport()
}
invisible(dev.off())
