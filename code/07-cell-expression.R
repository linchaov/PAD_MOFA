# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, Nebulosa 1.22.0
#
# Global: Localize genes encoding high-loading plasma proteins in muscle and femoral plaque.
#
# Inputs: Factor loadings, protein-to-gene annotations, tissue atlas objects and
#         donor-level expression summaries.
#
# Outputs: Tissue expression, density and donor-balanced summaries under ./results/VIP-PAD/ and
#         ./results/UKB/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
})



# ==============================================================================
# EXPRESSION OF FACTOR-ASSOCIATED GENES | FIGURE S05
# ==============================================================================

# Function to summarize cell expression of genes encoding high-loading proteins.
vip_build_loading_cell_dot_f1_f4_f8 <- function() {
  # VIP-PAD loading and cell-expression summaries for Factors 1–4 and 8.
  # Inputs include the fitted MOFA model and PAD single-cell donor aggregates.
  # Map and deduplicate the top 5% of proteins by absolute loading;
  # retain up to four genes per loading sign and factor.
  # RNA gate: >=3 independent donor-cell groups each >=10% nonzero cells;
  # each eligible donor-cell group >=20 cells. Z-score before filtering.
  invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
  dir.create(
    "./results/VIP-PAD/cell-expression",
    recursive = TRUE,
    showWarnings = FALSE
  )
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  meta <- readRDS("./results/VIP-PAD/Variable_Metadata.rds")
  ann <- as.data.table(meta$annotations$protein)
  mapping <- fread("./results/Tissue/results/tables/Protein_Gene_Mapping.csv")
  fs <- sprintf('F%02d', c(1:4, 8))
  # Match factor-specific protein weights to the tissue gene annotations.
  p <- rbindlist(lapply(c(1:4, 8), function(k) {
    w <- fit$W$Protein[, paste0('Factor', k)]
    ii <- match(names(w), ann$variable_name)
    x <- data.table(
      feature = names(w),
      label = ann$GeneSymbol[ii],
      loading = as.numeric(w)
    )
    x[,
      top5pct := rank(-abs(loading), ties.method = 'first') <= ceiling(.05 * .N)
    ]
    x <- merge(x, mapping[, .(label, gene, method)], by = 'label', all.x = TRUE)
    x <- x[!is.na(gene) & gene != ''][order(-abs(loading), feature)][
      !duplicated(gene)
    ]
    x[, `:=`(
      factor = sprintf('F%02d', k),
      loading_per_score_SD = loading * sd(fit$Z[, paste0('Factor', k)])
    )]
    x
  }))
  sel <- copy(p[top5pct == TRUE])
  sel[, direction := ifelse(loading_per_score_SD > 0, 'Positive', 'Negative')]
  sel <- sel[
    order(factor, direction, -abs(loading_per_score_SD), gene),
    head(.SD, 4),
    by = .(factor, direction)
  ]
  genes <- unique(sel$gene)
  ld <- p[factor %in% fs & gene %in% genes]

  fwrite(
    sel,
    "./results/VIP-PAD/cell-expression/Gene_selection.csv"
  )
  fwrite(
    ld,
    "./results/VIP-PAD/cell-expression/Loading_plot_data.csv"
  )
  o <- readRDS("./results/Tissue/data/cell_localization_objects.rds")
  all_dt <- list()
  for (key in c('Muscle', 'Plaque_main')) {
    obj <- o[[key]]
    cts <- colnames(obj$mean)
    pg <- as.data.table(obj$primary_groups)
    if (key == 'Muscle') {
      detmat <- readRDS("./results/Tissue/data/muscle_pseudobulk.rds")$detection
    } else {
      raw <- fread("./results/Tissue/data/plaque_main_detection.csv.gz")
      detmat <- as.matrix(raw[, -1])
      rownames(detmat) <- raw[[1]]
    }
    dt <- CJ(gene = genes, celltype = cts, sorted = FALSE)
    dt[, present := gene %in% rownames(obj$mean)]
    dt[, `:=`(
      mean_logCPM = NA_real_,
      detection_pct = NA_real_,
      z = NA_real_,
      n_samples_detected_10pct = 0L
    )]
    for (g in genes) {
      if (g %in% rownames(obj$mean)) {
        v <- obj$mean[g, cts]
        zv <- if (sd(v) > 0) as.numeric(scale(v)) else rep(0, length(v))
        dt[
          gene == g,
          `:=`(
            mean_logCPM = as.numeric(v),
            detection_pct = 100 * as.numeric(obj$detection[g, cts]),
            z = zv
          )
        ]
      }
    }
    for (ct in cts) {
      groups <- pg[celltype == ct, group]
      for (g in genes) {
        if (g %in% rownames(detmat)) {
          vals <- as.numeric(detmat[g, groups])
          dt[
            gene == g & celltype == ct,
            n_samples_detected_10pct := sum(vals >= .1)
          ]
        }
      }
    }
    dt[, `:=`(
      display_dot = present & n_samples_detected_10pct >= 3L,
      tissue = key
    )]
    all_dt[[key]] <- dt
  }

  fwrite(
    rbindlist(all_dt),
    "./results/VIP-PAD/cell-expression/Cell_expression_plot_data.csv"
  )
  cat(
    'UNIQUE GENES',
    length(genes),
    'MAX LOADING',
    max(abs(ld$loading_per_score_SD)),
    '\n'
  )
}



# ==============================================================================
# TISSUE EXPRESSION BY COHORT AND FACTOR
# ==============================================================================

# Function to estimate tissue-specific gene expression and cell density.
tissue_density_and_expression <- function(cohort, k) {
  data_dir <- file.path(
    "./results",
    cohort,
    "cell-inputs",
    paste0("Factor_", k)
  )
  table_dir <- file.path(
    "./results",
    cohort,
    "mofa",
    paste0("Factor_", k, "_Cell_Mapping")
  )
  S <- fread(file.path(table_dir, paste0("F", k, "_Features.csv")))
  for (tissue in c("Muscle", "Plaque")) {
    z <- readRDS(file.path(data_dir, paste0(tissue, '_Input.rds')))
    cells <- copy(z$cells)
    stopifnot(
      all(cells$condition == 'PAD'),
      all(!cells$pooled),
      identical(colnames(z$expression), cells$cell_id)
    )
    if (tissue == 'Muscle') {
      mapping <- c(
        'Endothelial cells' = 'Vascular',
        'Smooth muscle/Pericytes' = 'Vascular',
        'FAPs' = 'Stromal',
        'Satellite cells' = 'Myogenic',
        'Myonuclei' = 'Myogenic',
        'Macrophages/Monocytes' = 'Immune',
        'NK/T/B cells' = 'Immune',
        'Neutrophils' = 'Immune',
        'Mast cells' = 'Immune'
      )
      ctorder <- c(
        'Endothelial cells',
        'Smooth muscle/Pericytes',
        'FAPs',
        'Satellite cells',
        'Myonuclei',
        'Macrophages/Monocytes',
        'NK/T/B cells',
        'Neutrophils',
        'Mast cells'
      )
    } else {
      mapping <- c(
        'T cell' = 'Lymphoid',
        'NK cell' = 'Lymphoid',
        'B cell' = 'Lymphoid',
        'Plasma cell' = 'Lymphoid',
        'Monocyte' = 'Myeloid',
        'Macrophage' = 'Myeloid',
        'Dendritic cell' = 'Myeloid',
        'Neutrophil' = 'Myeloid',
        'Mast cell' = 'Myeloid',
        'EC' = 'Vascular',
        'Smooth Muscle Cell' = 'Vascular',
        'Fibroblast' = 'Stromal',
        'Fibromyocyte' = 'Stromal'
      )
      ctorder <- c(
        'T cell',
        'NK cell',
        'B cell',
        'Plasma cell',
        'Monocyte',
        'Macrophage',
        'Dendritic cell',
        'Neutrophil',
        'Mast cell',
        'EC',
        'Smooth Muscle Cell',
        'Fibroblast',
        'Fibromyocyte'
      )
    }
    cells[, compartment := unname(mapping[celltype])]
    cells[, celltype := factor(celltype, levels = ctorder)]
    cache <- file.path(data_dir, paste0(tissue, '_Gene_Densities.rds'))

    if (file.exists(cache)) {
      kd <- readRDS(cache)
      stopifnot(
        identical(kd$cell_id, cells$cell_id),
        identical(kd$genes, rownames(z$expression))
      )
      den <- kd$density
    } else {
      den <- matrix(
        NA_real_,
        nrow = nrow(cells),
        ncol = nrow(z$expression),
        dimnames = list(cells$cell_id, rownames(z$expression))
      )
      for (g in colnames(den)) {
        w <- z$expression[g, ]
        cat('DENSITY', tissue, g, '\n')
        flush.console()
        if (sum(w) > 0) {
          den[, g] <- Nebulosa:::calculate_density(
            w,
            as.matrix(cells[, .(x, y)]),
            method = 'ks'
          )
        }
      }

      saveRDS(
        list(
          cell_id = cells$cell_id,
          genes = rownames(z$expression),
          density = den
        ),
        cache
      )
    }
    requested <- S$gene
    use <- intersect(requested, colnames(den))
    use <- use[vapply(
      use,
      function(g) all(is.finite(den[, g])) && any(den[, g] > 0),
      logical(1)
    )]
    # Product is mathematically identical to Nebulosa's joint density; computed
    # on log scale to avoid floating-point underflow; then rescaled for display.
    ld <- rowSums(log(pmax(den[, use, drop = FALSE], .Machine$double.xmin)))
    rel <- exp(ld - max(ld))
    out <- cells[, .(cell_id, donor, celltype)]
    out[, `:=`(
      direction = "Top15",
      joint_density_relative = rel,
      log_joint_density = ld
    )]
    dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
    fwrite(
      out,
      file.path(table_dir, paste0(tissue, "_Top15_Joint_Density.csv.gz"))
    )
    dots <- list()
    donors <- list()
    for (ct in ctorder) {
      ii <- which(cells$celltype == ct)
      if (!length(ii)) {
        next
      }
      genes <- rownames(z$expression)
      dots[[ct]] <- data.table(
        gene = genes,
        celltype = ct,
        n_cells = length(ii),
        mean_expression = log1p(rowMeans(expm1(z$expression[,
          ii,
          drop = FALSE
        ]))),
        pct = 100 * rowMeans(z$counts[, ii, drop = FALSE] > 0)
      )
      for (donor in unique(cells$donor[ii])) {
        jj <- which(cells$celltype == ct & cells$donor == donor)
        donors[[paste(ct, donor)]] <- data.table(
          gene = genes,
          celltype = ct,
          donor = donor,
          n_cells = length(jj),
          mean_expression = log1p(rowMeans(expm1(z$expression[,
            jj,
            drop = FALSE
          ]))),
          pct = 100 * rowMeans(z$counts[, jj, drop = FALSE] > 0)
        )
      }
    }
    ds <- rbindlist(dots)
    don <- rbindlist(donors)
    balanced <- don[
      n_cells >= 20,
      .(
        mean_expression = mean(mean_expression),
        pct = mean(pct),
        n_donors = .N
      ),
      by = .(gene, celltype)
    ]

    fwrite(
      balanced,
      file.path(table_dir, paste0(tissue, '_Donor_Balanced_Expression.csv'))
    )
    ds <- merge(
      CJ(gene = S$gene, celltype = ctorder),
      ds,
      by = c('gene', 'celltype'),
      all.x = TRUE
    )
    support <- don[,
      .(
        n_samples_with_celltype = uniqueN(donor),
        n_samples_detected_10pct = uniqueN(donor[pct >= 10])
      ),
      by = .(gene, celltype)
    ]
    ds <- merge(ds, support, by = c('gene', 'celltype'), all.x = TRUE)
    ds[,
      display_dot := !is.na(pct) &
        !is.na(n_samples_detected_10pct) &
        n_samples_detected_10pct >= 3
    ]
    ds <- merge(
      ds,
      S[, .(gene, direction, absolute_loading_rank)],
      by = 'gene',
      all.x = TRUE
    )
    ds[, gene := factor(gene, levels = S$gene)]
    ds[, celltype := factor(celltype, levels = rev(ctorder))]
    ds[,
      direction := factor(
        direction,
        levels = c('Positive', 'Negative'),
        labels = c('Positive loadings', 'Negative loadings')
      )
    ]

    fwrite(
      ds,
      file.path(table_dir, paste0(tissue, '_Dotplot_Values.csv'))
    )
  }
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  vip_build_loading_cell_dot_f1_f4_f8()
  for (k in c(1, 2, 3, 4, 8)) {
    tissue_density_and_expression("VIP-PAD", k)
  }
  for (k in c(2, 5)) {
    tissue_density_and_expression("UKB", k)
  }
}
