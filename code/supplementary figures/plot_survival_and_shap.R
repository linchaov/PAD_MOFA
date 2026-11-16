# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1
#
# Global: S10, S13, S16, S19, S20 and S23 show event-free survival and SHAP contributions.
#
# Inputs: Factor-specific Kaplan–Meier estimates and SHAP results in VIP-PAD and UKB.
# Outputs: Individual panel PDFs under ./figures-output/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(grid)
})



# ==============================================================================
# FIGURES S13, S16, S19 AND S20 (PANEL A) — VIP-PAD event-free survival
# ==============================================================================

# Function to plot event-free survival by factor-score tertile.
draw_km_vip <- function(result_dir, output_prefix, k, outcome, high_colour) {

  a <- readRDS(file.path(result_dir, "km_results.rds"))
  lev <- c("Low", "Medium", "High")
  pal <- c(Low = "#42A5D9", Medium = "#A88F00", High = high_colour)
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
  H <- 5.6
  left <- 1.55
  right <- 5.63
  bottom <- 2.26
  top <- 4.72
  # Function to map a value to its horizontal panel coordinate.
  xx <- function(t) left + (right - left) * t / 36
  ylo <- max(
    0,
    floor(
      min(a$curves[endpoint == outcome & time <= 36, lower], na.rm = TRUE) * 10
    ) /
      10
  )
  # Function to map a value to its vertical panel coordinate.
  yy <- function(v) bottom + (top - bottom) * (v - ylo) / (1 - ylo)

  cairo_pdf(
    paste0(output_prefix, "A.pdf"),
    width = W,
    height = H,
    family = "Noto Sans"
  )
  grid.newpage()
  txt(outcome, (left + right) / 2, 5.4, 14)
  txt(paste("Factor", k, "tertiles"), left - 0.15, 5.02, 10.5, "left")
  ln(left, bottom, left, top, "#333333", 0.7)
  ln(left, bottom, right, bottom, "#333333", 0.7)
  yticks <- seq(ylo, 1, 0.1)
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
        endpoint == outcome & Tertile == g & time <= 36,
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
    lx <- left + 1.28 + (i - 1) * 1.02
    ly <- 5.02
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
  te <- a$tests[endpoint == outcome]
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
    ct <- a$counts[endpoint == outcome & Tertile == g]
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
          endpoint == outcome & Tertile == g & time == c(0, 12, 24, 36)[j],
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
# FIGURES S13, S16, S19 AND S20 (PANELS B–C) — VIP-PAD SHAP contributions
# ==============================================================================
# Function to plot training-set SHAP ranks and holdout SHAP values.
draw_shap_vip <- function(result_dir, output_prefix, outcome, high_colour) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
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
    file.path(result_dir, paste0(outcome, "_shap_results.rds"))
  )
  imp <- a$importance
  n <- nrow(imp)

  cairo_pdf(
    paste0(output_prefix, "B.pdf"),
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
  txt(paste0(outcome, ": SHAP elbow detection"), (left + right) / 2, 3.94, 13)
  txt("Molecule rank", (left + right) / 2, 0.2, 12)
  txt("Mean |SHAP|", 0.22, (bottom + top) / 2, 12, rot = 90)
  dev.off()
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
  cp <- colorRampPalette(c("#2315BE", "#B799E9", "#F4D7E6", high_colour))(256)

  cairo_pdf(
    paste0(output_prefix, "C.pdf"),
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
  txt(outcome, (left + right) / 2, H - 0.22, 13)
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
# FIGURES S10 AND S23 (PANEL A) — UKB event-free survival
# ==============================================================================

# Function to plot event-free survival by factor-score tertile.
draw_km_ukb <- function(result_dir, output_prefix, k, outcome, high_colour) {

  a <- readRDS(
    file.path(result_dir, "km_results.rds")
  )
  lev <- c("Low", "Medium", "High")
  pal <- c(Low = "#42A5D9", Medium = "#A88F00", High = high_colour)
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
  H <- 5.6
  left <- 1.55
  right <- 5.63
  bottom <- 2.26
  top <- 4.72
  # Function to map a value to its horizontal panel coordinate.
  xx <- function(t) left + (right - left) * t / 18
  ylo <- max(
    0,
    floor(
      min(a$curves[endpoint == outcome & time <= 18, lower], na.rm = TRUE) * 10
    ) /
      10
  )
  # Function to map a value to its vertical panel coordinate.
  yy <- function(v) bottom + (top - bottom) * (v - ylo) / (1 - ylo)

  cairo_pdf(
    paste0(output_prefix, "A.pdf"),
    width = W,
    height = H,
    family = "Noto Sans"
  )
  grid.newpage()
  txt(paste0("UKB: ", outcome), (left + right) / 2, 5.4, 14)
  txt(paste("Factor", k, "tertiles"), left - 0.15, 5.02, 10.5, "left")
  ln(left, bottom, left, top, "#333333", 0.7)
  ln(left, bottom, right, bottom, "#333333", 0.7)
  yticks <- seq(ylo, 1, 0.1)
  for (v in yticks) {
    ln(left - 0.04, yy(v), left, yy(v), "#555555", 0.6)
    txt(sprintf("%.1f", v), left - 0.09, yy(v), 12, "right")
  }
  for (t in c(0, 6, 12, 18)) {
    ln(xx(t), bottom, xx(t), bottom - 0.04, "#555555", 0.6)
    txt(as.character(t), xx(t), bottom - 0.16, 12)
  }
  txt("Event-free survival", left - 0.55, (top + bottom) / 2, 13, rot = 90)
  txt("Years", (left + right) / 2, 1.82, 13)
  ds <- list()
  for (g in lev) {
    d <- rbind(
      data.table(time = 0, S = 1, lower = 1, upper = 1, n_censor = 0),
      a$curves[
        endpoint == outcome & Tertile == g & time <= 18,
        .(time, S, lower, upper, n_censor)
      ]
    )
    setorder(d, time)
    last <- copy(tail(d, 1))
    last[, time := min(18, max(a$curves[endpoint == outcome & Tertile == g, time]))]
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
    cen <- d[n_censor > 0 & time < 18]
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
    lx <- left + 1.28 + (i - 1) * 1.02
    ly <- 5.02
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
  te <- a$tests[endpoint == outcome]
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
  txt("Years", rleft + 0.1, 1.21, 10.5, "left")
  for (j in 1:4) {
    txt(as.character(c(0, 6, 12, 18)[j]), rx[j], 1.21, 10.5)
  }
  for (i in 1:3) {
    g <- lev[i]
    ry <- 1.21 - i * 0.28
    ct <- a$counts[endpoint == outcome & Tertile == g]
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
          endpoint == outcome & Tertile == g & time == c(0, 6, 12, 18)[j],
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
# FIGURES S10 AND S23 (PANELS B–C) — UKB SHAP contributions
# ==============================================================================
# Function to plot training-set SHAP ranks and holdout SHAP values.
draw_shap_ukb <- function(result_dir, output_prefix, outcome, high_colour) {
  invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
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
    file.path(result_dir, paste0(outcome, "_shap_results.rds"))
  )
  imp <- a$importance
  n <- nrow(imp)

  cairo_pdf(
    paste0(output_prefix, "B.pdf"),
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
  for (v in c(0, 50, 100, 150)) {
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
    sprintf("Elbow = %d molecules", a$elbow),
    right - 0.05,
    top - 0.2,
    11,
    "right"
  )
  txt(
    paste0("UKB: ", outcome, ": SHAP elbow detection"),
    (left + right) / 2,
    3.94,
    13
  )
  txt("Molecule rank", (left + right) / 2, 0.2, 12)
  txt("Mean |SHAP|", 0.22, (bottom + top) / 2, 12, rot = 90)
  dev.off()
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
  cp <- colorRampPalette(c("#2315BE", "#B799E9", "#F4D7E6", high_colour))(256)

  cairo_pdf(
    paste0(output_prefix, "C.pdf"),
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
    pointcols <- cp[coloridx]
    pointcols[is.na(v)] <- "#A0A0A0"
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
    oo <- sample(seq_along(vals))
    grid.points(
      x = unit(xx(vals[oo]), "in"),
      y = unit(y(i) + off[oo], "in"),
      pch = 16,
      size = unit(0.041, "in"),
      gp = gpar(col = adjustcolor(pointcols[oo], alpha.f = 0.78))
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
  txt(paste0("UKB: ", outcome), (left + right) / 2, H - 0.22, 13)
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


dir.create("./figures-output", showWarnings = FALSE, recursive = TRUE)



# ==============================================================================
# FIGURE S10 (PANELS A–C) — UKB Factor 5, MACCE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/UKB/outcomes/Factor_5_MACCE"
output_prefix <- "./figures-output/FigureS10"
# FIGURE S10 (PANEL A) — Event-free survival
draw_km_ukb(result_dir, output_prefix, 5, "MACCE", "#F35A24")
# FIGURE S10 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_ukb(result_dir, output_prefix, "MACCE", "#FF3318")



# ==============================================================================
# FIGURE S13 (PANELS A–C) — VIP-PAD Factor 3, MACCE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/VIP-PAD/outcomes/Factor_3_MACCE"
output_prefix <- "./figures-output/FigureS13"
# FIGURE S13 (PANEL A) — Event-free survival
draw_km_vip(result_dir, output_prefix, 3, "MACCE", "#F35A24")
# FIGURE S13 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_vip(result_dir, output_prefix, "MACCE", "#FF3318")



# ==============================================================================
# FIGURE S16 (PANELS A–C) — VIP-PAD Factor 4, MACCE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/VIP-PAD/outcomes/Factor_4_MACCE"
output_prefix <- "./figures-output/FigureS16"
# FIGURE S16 (PANEL A) — Event-free survival
draw_km_vip(result_dir, output_prefix, 4, "MACCE", "#F45A24")
# FIGURE S16 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_vip(result_dir, output_prefix, "MACCE", "#FF2318")



# ==============================================================================
# FIGURE S19 (PANELS A–C) — VIP-PAD Factor 8, MALE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/VIP-PAD/outcomes/Factor_8_MALE"
output_prefix <- "./figures-output/FigureS19"
# FIGURE S19 (PANEL A) — Event-free survival
draw_km_vip(result_dir, output_prefix, 8, "MALE", "#F85A24")
# FIGURE S19 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_vip(result_dir, output_prefix, "MALE", "#FF2318")



# ==============================================================================
# FIGURE S20 (PANELS A–C) — VIP-PAD Factor 8, MACCE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/VIP-PAD/outcomes/Factor_8_MACCE"
output_prefix <- "./figures-output/FigureS20"
# FIGURE S20 (PANEL A) — Event-free survival
draw_km_vip(result_dir, output_prefix, 8, "MACCE", "#F85A24")
# FIGURE S20 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_vip(result_dir, output_prefix, "MACCE", "#FF2318")



# ==============================================================================
# FIGURE S23 (PANELS A–C) — UKB Factor 2, MACCE-free survival and molecular contributions
# ==============================================================================

result_dir <- "./results/UKB/outcomes/Factor_2_MACCE"
output_prefix <- "./figures-output/FigureS23"
# FIGURE S23 (PANEL A) — Event-free survival
draw_km_ukb(result_dir, output_prefix, 2, "MACCE", "#F35A24")
# FIGURE S23 — SHAP elbow (PANEL B) and molecular contributions (PANEL C)
draw_shap_ukb(result_dir, output_prefix, "MACCE", "#FF3318")
