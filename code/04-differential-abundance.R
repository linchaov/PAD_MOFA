# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, limma 3.68.5, statmod 1.5.2
#
# Global: Differential abundance by molecular cluster and Rutherford category.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds; fitted MOFA model, cluster assignments and feature
#         annotations.
#
# Outputs: Feature-level fold changes and molecule contrasts under
#         ./results/VIP-PAD/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to extract gene symbols from VIP-PAD protein identifiers.
vip_symbols <- function(features) {
  genes <- sub("_[^_]+$", "", sub("^PROT_", "", features))
  genes[genes == "Unmapped"] <- ""
  chartr("_", "-", genes)
}

# Function to restore the original measurement scales for abundance contrasts.
vip_abundance <- function() {
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  attrs <- attr(cohort, "omics_matrix_attributes")
  ids <- cohort$study_id
  rownames(cohort) <- ids
  assignments <- readRDS(
    "./results/VIP-PAD/clustering/Clustering_Stability.rds"
  )$assignments
  group <- assignments$K2[match(ids, assignments$sample_id)]
  stopifnot(!anyDuplicated(assignments$sample_id), setequal(group, 1:2))
  # Function to reverse feature-wise standardization with the saved centers and scales.
  original_scale <- function(view) {
    center <- attrs[[view]][["scaled:center"]]
    spread <- attrs[[view]][["scaled:scale"]]
    fields <- names(spread)
    x <- as.matrix(cohort[, fields, drop = FALSE])
    values <- sweep(sweep(x, 2, spread[fields], "*"), 2, center[fields], "+")
    stopifnot(all(is.finite(values)))
    values
  }
  clinical <- cohort[, c("Age", "Sex", "RutherfordCategory"), drop = FALSE]
  clinical$SampleID <- ids
  clinical$cluster <- group
  list(
    sample_id = ids,
    clinical = clinical,
    proteomics = original_scale("proteomics"),
    metabolomics = original_scale("metabolomics")
  )
}

# Function to identify leading Factor 1 proteins and metabolites by loading magnitude.
leading_f1_features <- function() {
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
  # Join signed loadings to protein and metabolite labels.
  loads <- rbindlist(lapply(c('Protein', 'Metabolite'), function(v) {
    w <- fit$W[[v]]
    ids <- rownames(w)
    label <- if (v == "Protein") {
      vip_symbols(ids)
    } else {
      metabolites$label[match(ids, metabolites$feature)]
    }
    label[is.na(label) | label == ''] <- rownames(w)[is.na(label) | label == '']
    data.table(
      feature = rownames(w),
      modality = v,
      label = label,
      loading_raw = w[, 'Factor1']
    )
  }))
  loads[,
    relative_weight := abs(loading_raw) / max(abs(loading_raw)),
    by = modality
  ]
  sel <- rbindlist(list(
    loads[modality == 'Protein'][order(-abs(loading_raw), feature)][1:15],
    loads[modality == 'Metabolite' & relative_weight > 0.8][order(
      -abs(loading_raw),
      feature
    )]
  ))
  sel[, direction := ifelse(loading_raw > 0, 'Positive', 'Negative')]
  sel[, vieworder := match(modality, c('Protein', 'Metabolite'))]
  setorder(sel, vieworder, -relative_weight, feature)
  sel$feature
}



# ==============================================================================
# ABUNDANCE DIFFERENCES BETWEEN VIP-PAD CLUSTERS | FIGURE 2 (PANELS D, F AND G)
# ==============================================================================

# Function to estimate adjusted abundance contrasts between clusters.
vip_compute_c1_vs_c2_log2fc <- function() {
  # Estimate age- and sex-adjusted log2 fold changes between clusters.
  a <- vip_abundance()
  metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
  w <- a$clinical
  X <- cbind(a$proteomics, a$metabolomics)
  grp <- w$cluster
  cov <- data.frame(
    group = factor(ifelse(grp == 1, 'C1', 'C2'), levels = c('C2', 'C1')),
    age = as.numeric(w$Age),
    sex = factor(w$Sex, levels = c('Female', 'Male'))
  )
  # Combine assay annotations in measurement order.
  features <- rbindlist(lapply(c('Protein', 'Metabolite'), function(v) {
    ids <- if (v == "Protein") colnames(a$proteomics) else colnames(a$metabolomics)
    labels <- if (v == "Protein") {
      vip_symbols(ids)
    } else {
      metabolites$label[match(ids, metabolites$feature)]
    }
    data.table(feature = ids, modality = v, label = labels)
  }))
  features <- features[match(colnames(X), feature)]
  features[is.na(label) | label == '', label := feature]
  stopifnot(
    identical(colnames(X), features$feature),
    all(apply(X, 2, sd) > 0)
  )
  design <- model.matrix(~ group + age + sex, cov)
  stopifnot(qr(design)$rank == ncol(design))
  results <- list()
  for (view in c('Protein', 'Metabolite')) {
    idx <- which(features$modality == view)
    Y <- t(X[, idx, drop = FALSE])
    fit <- eBayes(lmFit(Y, design), trend = FALSE, robust = TRUE)
    b <- fit$coefficients[, 'groupC1']
    se <- fit$stdev.unscaled[, 'groupC1'] * sqrt(fit$s2.post)
    critical <- qt(.975, fit$df.total)
    P <- fit$p.value[, 'groupC1']
    r <- copy(features[idx])
    r[, `:=`(
      N = nrow(cov),
      N_C1 = sum(cov$group == 'C1'),
      N_C2 = sum(cov$group == 'C2'),
      log2FC_C1_vs_C2 = b,
      SE_moderated = se,
      CI95_lower = b - critical * se,
      CI95_upper = b + critical * se,
      t_moderated = fit$t[, 'groupC1'],
      P = P,
      FDR = p.adjust(P, 'BH'),
      mean_log2_C1 = colMeans(X[cov$group == 'C1', idx, drop = FALSE]),
      mean_log2_C2 = colMeans(X[cov$group == 'C2', idx, drop = FALSE])
    )]
    results[[view]] <- r
  }
  r <- rbindlist(results)
  r[,
    category := fifelse(
      FDR >= .05,
      'Not significant',
      fifelse(log2FC_C1_vs_C2 > 0, 'Higher in C1', 'Higher in C2')
    )
  ]
  r[,
    category := factor(
      category,
      levels = c('Higher in C2', 'Not significant', 'Higher in C1')
    )
  ]
  # Bound plotting coordinates and export the original FDR.
  r[, neglog10FDR := -log10(pmax(FDR, 1e-300))]
  dir.create(
    "./results/VIP-PAD/abundance",
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(
    r,
    "./results/VIP-PAD/abundance/All_Molecular_Contrasts.csv",
    bom = TRUE
  )
  counts <- r[,
    .(
      tested = .N,
      FDR_below_005 = sum(FDR < .05),
      higher_C1 = sum(FDR < .05 & log2FC_C1_vs_C2 > 0),
      higher_C2 = sum(FDR < .05 & log2FC_C1_vs_C2 < 0),
      not_significant = sum(FDR >= .05),
      min_FDR = min(FDR),
      minimum_beta = min(log2FC_C1_vs_C2),
      maximum_beta = max(log2FC_C1_vs_C2),
      maximum_logFDR = max(neglog10FDR)
    ),
    by = modality
  ]
  fwrite(
    counts,
    "./results/VIP-PAD/abundance/Contrast_Summary.csv"
  )

}



# ==============================================================================
# FACTOR 1 ABUNDANCE CONTRASTS | FIGURE 4 (PANEL D)
# ==============================================================================

# Function to compare Factor 1 molecules by cluster and Rutherford grade.
vip_compute_f1_abundance_contrasts <- function() {
  # Fit separate limma contrasts for Rutherford grade and cluster. Adjust P-values within each
  # omics layer.
  if (.Platform$OS.type == 'windows') {
    invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
  }

  dir.create(
    "./results/VIP-PAD/abundance",
    recursive = TRUE,
    showWarnings = FALSE
  )
  a <- vip_abundance()
  metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
  leading_features <- leading_f1_features()
  w <- a$clinical
  X <- cbind(a$proteomics, a$metabolomics)
  ruth <- w$RutherfordCategory
  grp <- w$cluster
  cov <- data.frame(
    group = factor(ifelse(grp == 1, 'C1', 'C2'), levels = c('C2', 'C1')),
    severity = factor(
      ifelse(is.na(ruth), NA, ifelse(ruth >= 4, 'R4-6', 'R0-3')),
      levels = c('R0-3', 'R4-6')
    ),
    age = as.numeric(w$Age),
    sex = factor(w$Sex, levels = c('Female', 'Male'))
  )
  keep <- complete.cases(cov) & is.finite(cov$age)
  cov <- droplevels(cov[keep, ])
  xx <- X[keep, , drop = FALSE]
  # Combine assay annotations in measurement order.
  features <- rbindlist(lapply(c('Protein', 'Metabolite'), function(v) {
    ids <- if (v == "Protein") colnames(a$proteomics) else colnames(a$metabolomics)
    labels <- if (v == "Protein") {
      vip_symbols(ids)
    } else {
      metabolites$label[match(ids, metabolites$feature)]
    }
    data.table(feature = ids, modality = v, label = labels)
  }))
  features <- features[match(colnames(X), feature)]
  res <- list()
  for (v in c('Protein', 'Metabolite')) {
    for (contrast in c('Rutherford', 'Cluster')) {
      idx <- which(features$modality == v)
      des <- if (contrast == 'Rutherford') {
        model.matrix(~ severity + age + sex, cov)
      } else {
        model.matrix(~ group + age + sex, cov)
      }
      stopifnot(qr(des)$rank == ncol(des))
      cn <- if (contrast == 'Rutherford') 'severityR4-6' else 'groupC1'
      fit <- eBayes(
        lmFit(t(xx[, idx, drop = FALSE]), des),
        trend = FALSE,
        robust = TRUE
      )
      tab <- copy(features[idx])
      tab[, `:=`(
        comparison = contrast,
        N = nrow(cov),
        log2FC = fit$coefficients[, cn],
        P = fit$p.value[, cn]
      )]
      tab[, FDR := p.adjust(P, 'BH')]
      ss <- fit$stdev.unscaled[, cn] * sqrt(fit$s2.post)
      tab[, `:=`(
        CI_lower = log2FC - qt(.975, fit$df.total) * ss,
        CI_upper = log2FC + qt(.975, fit$df.total) * ss
      )]
      res[[paste(v, contrast)]] <- tab
    }
  }
  all <- rbindlist(res)
  d <- merge(
    all[
      comparison == 'Rutherford',
      .(feature, modality, label, x = log2FC, P_x = P, FDR_x = FDR)
    ],
    all[comparison == 'Cluster', .(feature, y = log2FC, P_y = P, FDR_y = FDR)],
    by = 'feature'
  )
  d[, F1_display := feature %in% leading_features]

  fwrite(
    d,
    "./results/VIP-PAD/abundance/Plot_Data.csv"
  )

}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  vip_compute_c1_vs_c2_log2fc()
  vip_compute_f1_abundance_contrasts()
}
