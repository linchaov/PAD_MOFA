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
# Global: Cross-cohort alignment and correlations of protein loadings.
#
# Inputs: Fitted MOFA models, VIP-PAD assay annotations and the feature dictionary in
#         ./Datasets_PAD_MOFA/UKB.rds.
#
# Outputs: Shared-protein mappings and factor-pair correlations under
#         ./results/cross-cohort/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
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



# ==============================================================================
# CROSS-COHORT LOADING CORRELATIONS | FIGURE 3 (PANELS B–D)
# ==============================================================================

# Function to compare cohort-specific loadings for shared proteins.
cross_compare_factors <- function() {
  A <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  B <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
  ids <- rownames(A$W$Protein)
  pa <- data.table(
    variable_name = ids,
    Accession = sub(".*_", "", ids),
    GeneSymbol = vip_symbols(ids)
  )
  ukb <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  pb <- as.data.table(attr(ukb, "feature_dictionary"))[modality == "Protein"]
  # Match proteins by unique UniProt accession in both cohorts.
  pa[,
    unique_id := !duplicated(Accession) &
      !duplicated(Accession, fromLast = TRUE) &
      grepl('^[A-Z0-9]+$', Accession)
  ]
  pb[,
    unique_id := !duplicated(uniprot_id) &
      !duplicated(uniprot_id, fromLast = TRUE) &
      grepl('^[A-Z0-9]+$', uniprot_id)
  ]
  P <- merge(
    pa[
      unique_id == TRUE,
      .(id = Accession, fuwai = variable_name, label = GeneSymbol)
    ],
    pb[unique_id == TRUE, .(id = uniprot_id, ukb = feature, ukb_label = label)],
    by = 'id'
  )
  dir.create(
    "./results/cross-cohort",
    recursive = TRUE,
    showWarnings = FALSE
  )
  # Compare the same UniProt accessions across all pairs of MOFA factors.
  vip_weights <- A$W$Protein[P$fuwai, , drop = FALSE]
  ukb_weights <- B$W$Protein[P$ukb, , drop = FALSE]
  pearson_r <- cor(vip_weights, ukb_weights)
  spearman_r <- cor(vip_weights, ukb_weights, method = "spearman")
  D <- CJ(
    fuwai_factor = seq_len(ncol(vip_weights)),
    ukb_factor = seq_len(ncol(ukb_weights))
  )
  pair_index <- as.matrix(D)
  D[, `:=`(
    view = "Protein",
    n_shared = nrow(P),
    rho = spearman_r[pair_index],
    pearson = pearson_r[pair_index]
  )]
  D[, `:=`(abs_rho = abs(rho), abs_r = abs(pearson))]
  D[,
    reciprocal_best := ukb_factor ==
      max.col(abs(pearson_r), ties.method = "first")[fuwai_factor] &
      fuwai_factor ==
        max.col(t(abs(pearson_r)), ties.method = "first")[ukb_factor]
  ]
  D[, display_link := FALSE]

  # Display the three factor pairs with the largest absolute Pearson correlations.
  candidates <- D[order(-abs_r)][1:3]
  candidates[, display_link := TRUE]
  D[,
    display_link := view == 'Protein' &
      paste(fuwai_factor, ukb_factor) %in%
        paste(candidates$fuwai_factor, candidates$ukb_factor)
  ]
  # Function to extract matched protein loadings and label the columns by cohort.
  extract <- function(fit, features, prefix) {
    m <- fit$W$Protein[features, , drop = FALSE]
    colnames(m) <- paste0(prefix, seq_len(ncol(m)))
    as.data.table(m)
  }
  shared <- cbind(
    P[, .(feature_a = fuwai, feature_b = ukb, label)],
    extract(A, P$fuwai, 'F'),
    extract(B, P$ukb, 'U')
  )

  saveRDS(
    list(correlations = D, candidates = candidates, loadings = shared),
    "./results/cross-cohort/Comparison.rds"
  )
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) cross_compare_factors()
