# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, fgsea 1.38.0, lipidr 2.26.0, MOFA2 1.22.1
#
# Global: Protein pathway and metabolite-class enrichment in VIP-PAD and UKB.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds and UKB.rds, fitted MOFA models, assay and lipid
#         annotations, and study gene sets.
#
# Outputs: Enrichment results and displayed pathways under ./results/VIP-PAD/, ./results/UKB/
#         and ./results/cross-cohort/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(MOFA2)
  library(lipidr)
  library(fgsea)
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to assign overlapping biological themes to pathways in either cohort.
pathway_themes <- function(info) {
  themes <- list(
    'Inflammation & immunity' = c(
      'IMMUN',
      'INFLAM',
      'LEUKOCYTE',
      'LYMPHOCYTE',
      'CYTOKINE',
      'INTERLEUKIN',
      'INTERFERON',
      'CHEMOKINE',
      'NEUTROPHIL',
      'MACROPHAGE',
      'TOLL_LIKE',
      'ANTIGEN',
      '_T_CELL',
      '_B_CELL',
      'DEFENSE',
      'TUMOR_NECROSIS'
    ),
    'Complement' = c('COMPLEMENT'),
    'Coagulation & platelets' = c(
      'COAGULATION',
      'PLATELET',
      'HEMOSTASIS',
      'HAEMOSTASIS',
      'FIBRIN',
      'THROMB',
      'PLASMIN'
    ),
    'Lipid metabolism' = c(
      'LIPID',
      'LIPOPROTEIN',
      'CHOLESTEROL',
      'FATTY_ACID',
      'TRIGLYCERIDE',
      'TRIACYLGLYCEROL',
      'PHOSPHOLIPID',
      'SPHINGO',
      'STEROL',
      'APOLIPOPROTEIN'
    ),
    'Endothelium & angiogenesis' = c(
      'ENDOTHEL',
      'ANGIOGEN',
      'VASCUL',
      'BLOOD_VESSEL',
      'VEGF',
      'VASODIL',
      'VASOCONSTR',
      'NITRIC_OXIDE'
    ),
    'ECM & remodeling' = c(
      'EXTRACELLULAR_MATRIX',
      'COLLAGEN',
      'INTEGRIN',
      'CELL_ADHESION',
      'CELL_SUBSTRATE_ADHESION',
      'BASEMENT_MEMBRANE',
      'ELASTIN',
      'MATRIX_METALLO',
      'CONNECTIVE_TISSUE',
      'FIBROSIS'
    ),
    'Metabolism & oxidative stress' = c(
      'MITOCHON',
      'OXIDATIVE',
      'REACTIVE_OXYGEN',
      'RESPIRATORY',
      'ELECTRON_TRANSPORT',
      'GLYCOL',
      'GLUCOSE',
      'TRICARBOXYLIC',
      'CITRIC_ACID',
      'ATP_SYNTH',
      'ANTIOXIDANT',
      'GLUTATHIONE'
    ),
    'Cell injury & repair' = c(
      'APOPT',
      'AUTOPHAG',
      'CELL_DEATH',
      'WOUND',
      'REPAIR',
      'UNFOLDED_PROTEIN',
      'HEAT_SHOCK',
      'CHAPERONE',
      'PROTEASOM',
      'PROTEOSTASIS',
      'PROTEIN_FOLDING',
      'HYPOXIA',
      'ISCHEM',
      'SENESCENCE'
    )
  )
  # Assign pathways to the displayed biological themes.
  map <- rbindlist(lapply(names(themes), function(th) {
    info[
      grepl(paste(themes[[th]], collapse = '|'), pathway),
      .(pathway, theme = th)
    ]
  }))
  map <- map[
    !(theme == 'Inflammation & immunity' & grepl('COMPLEMENT', pathway))
  ]
  rbind(
    map,
    info[!pathway %in% map$pathway, .(pathway, theme = 'Other')]
  )
}

# Function to extract gene symbols from VIP-PAD protein identifiers.
vip_symbols <- function(features) {
  genes <- sub("_[^_]+$", "", sub("^PROT_", "", features))
  genes[genes == "Unmapped"] <- ""
  chartr("_", "-", genes)
}



# ==============================================================================
# VIP-PAD PROTEIN AND LIPID ENRICHMENT | FIGURE 4 (PANEL E) AND FIGURE 5 (PANEL F)
# ==============================================================================

# Function to test VIP-PAD factor loadings for protein pathways and lipid classes.
vip_factor_enrichment <- function(k) {
  factor_name <- paste0("Factor", k)
  factor_label <- sprintf("F%02d", k)
  # Use two factors in the MOFA2 enrichment object and extract the target factor.
  pair <- c(1L, max(2L, k))
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  mol <- readRDS("../inputs/metabolites.rds")
  measured <- names(attr(cohort, "omics_matrix_attributes")$proteomics[["scaled:scale"]])
  mol$protein <- data.table(feature = measured, label = vip_symbols(measured))
  sets <- readRDS("../inputs/vip-sets.rds")
  ids <- rownames(fit$Z)
  stopifnot(!anyDuplicated(ids), setequal(ids, cohort$study_id))
  rownames(cohort) <- cohort$study_id

  # Match factor loadings to the annotated proteins and metabolites.
  for (view in c('protein', 'metabolite')) {
    w <- fit$W[[if (view == 'protein') 'Protein' else 'Metabolite']][,
      factor_name
    ]
    tab <- as.data.table(mol[[view]])
    tab[, loading_raw := unname(w[feature])]
    mol[[view]] <- tab
  }
  lip <- as.data.table(mol$lipid_molecules)
  lip[,
    loading_raw := mol$metabolite$loading_raw[match(
      feature,
      mol$metabolite$feature
    )]
  ]
  mol$lipid_molecules <- lip

  # Test eight lipid classes against the measured lipid background.
  classes <- c('AcCar', 'LPC', 'LPC-P', 'LPE', 'PC', 'MG', 'MG-O', 'DG')
  class_vec <- ifelse(
    lip$plot_class %in% classes,
    lip$plot_class,
    NA_character_
  )
  li <- data.frame(
    Molecule = lip$feature,
    Class = class_vec,
    total_cl = NA_real_,
    total_cs = NA_real_,
    istd = FALSE,
    contrast = paste0('Factor_', k),
    logFC = lip$loading_raw
  )
  set.seed(20261020)
  lr <- as.data.table(lipidr::lsea(
    li,
    rank.by = 'logFC',
    min_size = 5,
    nproc = 1,
    eps = 0,
    sampleSize = 1001,
    nPermSimple = 100000
  ))
  lr[, `:=`(
    FDR = p.adjust(pval, 'BH'),
    lipid_class = sub('^Class_', '', set),
    direction = ifelse(NES > 0, 'Positive', 'Negative')
  )]
  lr[,
    leadingEdge_features := vapply(
      leadingEdge,
      paste,
      collapse = '; ',
      FUN.VALUE = character(1)
    )
  ]
  lr[,
    leadingEdge_names := vapply(
      leadingEdge,
      function(ids) {
        paste(lip$unified_label[match(ids, lip$feature)], collapse = '; ')
      },
      character(1)
    )
  ]
  lr[, c('leadingEdge', 'padj') := NULL]
  mol$lipid_results <- lr
  # Test protein pathways against the 3,030 mapped genes using the MOFA2 rank-sum test.
  B <- sets$B
  map <- sets$map
  stopifnot(identical(colnames(B), map$feature), !anyDuplicated(map$feature))
  W <- fit$W$Protein[map$feature, paste0('Factor', pair), drop = FALSE]
  Z <- fit$Z[ids, paste0('Factor', pair), drop = FALSE]
  colnames(W) <- colnames(Z) <- sprintf('F%02d', pair)
  X <- t(as.matrix(cohort[, map$feature, drop = FALSE])[ids, , drop = FALSE])
  obj <- new(
    'MOFA',
    data = list(Protein = list(PAD = X)),
    expectations = list(W = list(Protein = W), Z = list(PAD = Z)),
    data_options = list(views = 'Protein', groups = 'PAD'),
    dimensions = list(
      M = 1L,
      G = 1L,
      N = c(PAD = nrow(Z)),
      D = c(Protein = nrow(W)),
      K = 2L
    ),
    model_options = list(likelihoods = c(Protein = 'gaussian')),
    intercepts = list(
      Protein = list(PAD = setNames(rep(0, nrow(W)), rownames(W)))
    ),
    status = 'trained',
    on_disk = FALSE
  )
  obj@features_metadata <- data.frame(feature = map$feature, view = 'Protein')
  obj@samples_metadata <- data.frame(sample = ids, group = 'PAD')
  a <- MOFA2::run_enrichment(
    obj,
    view = 'Protein',
    feature.sets = B,
    factors = sprintf('F%02d', pair),
    set.statistic = 'rank.sum',
    statistical.test = 'parametric',
    sign = 'all',
    min.size = 10,
    p.adj.method = 'BH',
    alpha = .05,
    verbose = FALSE
  )
  ua <- data.table(
    test_id = rownames(a$pval),
    direction = "all",
    P = a$pval[, factor_label],
    statistic = a$set.statistics[, factor_label]
  )
  ua[, FDR := p.adjust(P, 'BH')]
  ea <- merge(
    sets$aliases,
    ua,
    by = 'test_id',
    allow.cartesian = TRUE,
    sort = FALSE
  )
  mol$enrichment <- ea
  dir.create("./results/VIP-PAD/enrichment", recursive = TRUE, showWarnings = FALSE)

  saveRDS(mol, sprintf("./results/VIP-PAD/enrichment/F%d.rds", k))
}



# ==============================================================================
# VIP-PAD ABSOLUTE LOADING ENRICHMENT | FIGURE 3 (PANEL A)
# ==============================================================================

# Function to test VIP-PAD protein pathways using absolute factor loadings.
cross_absolute_enrichment <- function() {
  dir.create(
    "./results/cross-cohort/enrichment/VIP-PAD",
    recursive = TRUE,
    showWarnings = FALSE
  )
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  ids <- names(attr(cohort, "omics_matrix_attributes")$proteomics[["scaled:scale"]])
  annotation <- data.table(variable_name = ids, GeneSymbol = vip_symbols(ids))
  loads <- annotation[
    variable_name %in% rownames(fit$W$Protein),
    .(feature = variable_name, gene = toupper(GeneSymbol))
  ]
  mem <- fread(
    "../inputs/gene-sets.csv"
  )
  loads[,
    retained := !is.na(gene) &
      grepl('^[A-Z0-9-]+$', gene) &
      !duplicated(gene) &
      !duplicated(gene, fromLast = TRUE)
  ]
  loads <- loads[retained == TRUE]

  fwrite(loads[, .(feature, gene)], "./results/cross-cohort/enrichment/VIP-PAD/Gene_Mapping.csv")
  info <- mem[
    gene %in% loads$gene,
    .(measured_size = uniqueN(gene)),
    by = .(pathway, collection)
  ]
  info[, eligible := measured_size >= 10L]
  fwrite(
    info,
    "./results/cross-cohort/enrichment/VIP-PAD/Gene_Set_Coverage.csv"
  )
  info <- info[eligible == TRUE]
  W <- fit$W$Protein[loads$feature, , drop = FALSE]
  Z <- fit$Z
  cat('Mapped proteins:', nrow(W), 'eligible sets:', nrow(info), '\n')
  flush.console()
  rownames(cohort) <- cohort$study_id
  protein <- as.matrix(cohort[, loads$feature, drop = FALSE])
  X <- t(protein[rownames(Z), , drop = FALSE])
  B <- matrix(
    0,
    nrow(info),
    nrow(loads),
    dimnames = list(info$pathway, loads$feature)
  )
  mm <- mem[pathway %in% info$pathway & gene %in% loads$gene]
  B[cbind(match(mm$pathway, rownames(B)), match(mm$gene, loads$gene))] <- 1
  obj <- new(
    'MOFA',
    data = list(Protein = list(PAD = X)),
    expectations = list(W = list(Protein = W), Z = list(PAD = Z)),
    data_options = list(views = 'Protein', groups = 'PAD'),
    dimensions = list(
      M = 1L,
      G = 1L,
      N = c(PAD = nrow(Z)),
      D = c(Protein = nrow(W)),
      K = ncol(W)
    ),
    model_options = list(likelihoods = c(Protein = 'gaussian')),
    intercepts = list(
      Protein = list(PAD = setNames(rep(0, nrow(W)), rownames(W)))
    ),
    status = 'trained',
    on_disk = FALSE
  )
  obj@features_metadata <- data.frame(feature = rownames(W), view = 'Protein')
  obj@samples_metadata <- data.frame(sample = rownames(Z), group = 'PAD')
  r <- run_enrichment(
    obj,
    view = 'Protein',
    feature.sets = B,
    factors = colnames(W),
    set.statistic = 'rank.sum',
    statistical.test = 'parametric',
    sign = 'all',
    min.size = 10,
    p.adj.method = 'BH',
    alpha = .05,
    verbose = FALSE
  )
  # Collect the pathway statistics for each factor.
  G <- rbindlist(lapply(seq_len(ncol(W)), function(k) {
    data.table(
      pathway = rownames(r$pval),
      factor = k,
      pval = r$pval[, k],
      padj = r$pval.adj[, k],
      statistic = r$set.statistics[, k]
    )
  }))
  G <- merge(G, info[, .(pathway, collection)], by = 'pathway')

  fwrite(
    G,
    "./results/cross-cohort/enrichment/VIP-PAD/All_Tests.csv"
  )
  map <- pathway_themes(info)

  fwrite(
    map,
    "./results/cross-cohort/enrichment/VIP-PAD/Theme_Mapping.csv"
  )
  sg <- merge(
    G[padj < .05 & statistic > 0],
    map,
    by = 'pathway',
    allow.cartesian = TRUE
  )
  sg[, display_rank := rank(padj, ties.method = 'min'), by = .(factor, theme)]
  sel <- sg[display_rank <= 30]
  fwrite(
    sel,
    "./results/cross-cohort/enrichment/VIP-PAD/Displayed_Terms.csv"
  )
}



# ==============================================================================
# SIGNED LOADING ENRICHMENT | FIGURE 3 (PANEL A)
# ==============================================================================

# Function to test pathways using signed protein loadings in both cohorts.
cross_signed_enrichment <- function() {
  # Directional terms determine the published heatmap column order, including blank cells.
  base <- "./results/cross-cohort/enrichment"
  mem <- fread(
    "../inputs/gene-sets.csv"
  )
  for (cohort in c("VIP-PAD", "UKB")) {
    input <- file.path(base, cohort)
    w <- fread(file.path(input, "Gene_Mapping.csv"), select = c("feature", "gene"))
    fit <- readRDS(file.path(
      "./results",
      cohort,
      "mofa/Selected_Fit.rds"
    ))
    W <- fit$W$Protein[w$feature, , drop = FALSE]
    K <- ncol(W)
    N <- nrow(W)
    A <- fread(file.path(input, "All_Tests.csv"))
    info <- unique(A[, .(pathway, collection)])
    map <- fread(file.path(input, "Theme_Mapping.csv"))
    B <- matrix(0, nrow(info), N, dimnames = list(info$pathway, w$feature))
    mm <- mem[pathway %in% info$pathway & gene %in% w$gene]
    B[cbind(match(mm$pathway, rownames(B)), match(mm$gene, w$gene))] <- 1
    m <- rowSums(B)
    n <- N - m
    # Use mean ranks for ties and the MOFA2 parametric rank-sum variance.
    # Function to test pathway enrichment using ranked loadings.
    engine <- function(T, direction) {
      ranks <- apply(T, 2, rank, ties.method = 'average')
      u <- B %*% ranks - m * (m + 1) / 2
      z <- (u - m * n / 2) / sqrt(m * n * (N + 1) / 12)
      p <- 2 * pnorm(-abs(z))
      if (direction != 'all') {
        p[z < 0] <- 1
      }
      # Collect the pathway statistics for each factor.
      rbindlist(lapply(seq_len(K), function(k) {
        data.table(
          pathway = rownames(B),
          factor = k,
          direction = direction,
          pval = p[, k],
          statistic = z[, k]
        )
      }))
    }
    G <- rbind(engine(pmax(W, 0), 'positive'), engine(pmax(-W, 0), 'negative'))
    # Joint BH over both directions and every eligible set within a factor.
    G[, padj := p.adjust(pval, 'BH'), by = factor]
    G <- merge(G, info, by = 'pathway')
    sg <- merge(
      G[padj < .05 & statistic > 0],
      map,
      by = 'pathway',
      allow.cartesian = TRUE
    )
    sg[,
      display_rank := rank(padj, ties.method = 'min'),
      by = .(factor, direction, theme)
    ]
    sel <- sg[display_rank <= 30]

    fwrite(sel, file.path(input, "Signed_Displayed_Terms.csv"))
  }
}



# ==============================================================================
# UKB ABSOLUTE LOADING ENRICHMENT | FIGURE 3 (PANEL A)
# ==============================================================================

# Function to test UKB protein pathways using absolute factor loadings.
cross_ukb_absolute_enrichment <- function() {
  dir.create(
    "./results/cross-cohort/enrichment/UKB",
    recursive = TRUE,
    showWarnings = FALSE
  )
  fit <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
  cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  annotation <- as.data.table(attr(cohort, "feature_dictionary"))
  loads <- annotation[
    modality == "Protein" & feature %in% rownames(fit$W$Protein),
    .(feature, source_field, gene = toupper(label))
  ]
  loads <- loads[
    !is.na(gene) &
      grepl("^[A-Z0-9-]+$", gene) &
      !duplicated(gene) &
      !duplicated(gene, fromLast = TRUE)
  ]
  fwrite(
    loads[, .(feature, gene)],
    "./results/cross-cohort/enrichment/UKB/Gene_Mapping.csv"
  )
  mem <- fread(
    "../inputs/gene-sets.csv"
  )
  info <- mem[
    gene %in% loads$gene,
    .(measured_size = uniqueN(gene)),
    by = .(pathway, collection)
  ]
  info[, eligible := measured_size >= 10L]

  fwrite(
    info,
    "./results/cross-cohort/enrichment/UKB/Gene_Set_Coverage.csv"
  )
  info <- info[eligible == TRUE]
  W <- fit$W$Protein[loads$feature, , drop = FALSE]
  Z <- fit$Z
  rownames(cohort) <- cohort$study_id
  # Match UKB fitting: standardize each source field over all observed participants.
  protein <- as.matrix(cohort[, loads$source_field, drop = FALSE])
  for (j in seq_len(ncol(protein))) {
    z <- protein[, j]
    s <- sd(z, na.rm = TRUE)
    protein[, j] <- (z - mean(z, na.rm = TRUE)) / s
  }
  colnames(protein) <- loads$feature
  X <- t(protein[rownames(Z), , drop = FALSE])
  B <- matrix(
    0,
    nrow(info),
    nrow(loads),
    dimnames = list(info$pathway, loads$feature)
  )
  mm <- mem[pathway %in% info$pathway & gene %in% loads$gene]
  B[cbind(match(mm$pathway, rownames(B)), match(mm$gene, loads$gene))] <- 1
  obj <- new(
    'MOFA',
    data = list(Protein = list(PAD = X)),
    expectations = list(W = list(Protein = W), Z = list(PAD = Z)),
    data_options = list(views = 'Protein', groups = 'PAD'),
    dimensions = list(
      M = 1L,
      G = 1L,
      N = c(PAD = nrow(Z)),
      D = c(Protein = nrow(W)),
      K = ncol(W)
    ),
    model_options = list(likelihoods = c(Protein = 'gaussian')),
    intercepts = list(
      Protein = list(PAD = setNames(rep(0, nrow(W)), rownames(W)))
    ),
    status = 'trained',
    on_disk = FALSE
  )
  obj@features_metadata <- data.frame(feature = rownames(W), view = 'Protein')
  obj@samples_metadata <- data.frame(sample = rownames(Z), group = 'PAD')
  r <- run_enrichment(
    obj,
    view = 'Protein',
    feature.sets = B,
    factors = colnames(W),
    set.statistic = 'rank.sum',
    statistical.test = 'parametric',
    sign = 'all',
    min.size = 10,
    p.adj.method = 'BH',
    alpha = .05,
    verbose = FALSE
  )
  # Collect the pathway statistics for each factor.
  G <- rbindlist(lapply(seq_len(ncol(W)), function(k) {
    data.table(
      pathway = rownames(r$pval),
      factor = k,
      pval = r$pval[, k],
      padj = r$pval.adj[, k],
      statistic = r$set.statistics[, k]
    )
  }))
  G <- merge(G, info[, .(pathway, collection)], by = 'pathway')

  fwrite(
    G,
    "./results/cross-cohort/enrichment/UKB/All_Tests.csv"
  )
  map <- pathway_themes(info)

  fwrite(
    map,
    "./results/cross-cohort/enrichment/UKB/Theme_Mapping.csv"
  )
  sg <- merge(
    G[padj < .05 & statistic > 0],
    map,
    by = 'pathway',
    allow.cartesian = TRUE
  )
  sg[, display_rank := rank(padj, ties.method = 'min'), by = .(factor, theme)]
  sel <- sg[display_rank <= 30]
  fwrite(
    sel,
    "./results/cross-cohort/enrichment/UKB/Displayed_Terms.csv"
  )
}



# ==============================================================================
# UKB NMR CLASS ENRICHMENT
# ==============================================================================

# Function to test enrichment of NMR classes in UKB factors.
ukb_nmr_enrichment <- function() {
  cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  ann <- as.data.table(attr(cohort, "feature_dictionary"))
  fit <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
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
  ann[,
    display_label := ifelse(
      modality == 'Protein',
      toupper(label),
      pretty_nmr(label)
    )
  ]
  dir.create(
    "./results/UKB/factor-profiles",
    recursive = TRUE,
    showWarnings = FALSE
  )
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

  classes <- c('VLDL', 'IDL', 'LDL', 'HDL', 'Fatty acids')
  res <- list()
  for (k in c(2, 5)) {
    # Group UKB NMR traits by measure name.
    nmr <- loads_for(k)[modality == 'Metabolite']
    nmr[, original := ann$label[match(feature, ann$feature)]]
    nmr[,
      category := fcase(
        grepl('VLDL', original) , 'VLDL' , grepl('IDL', original) , 'IDL' , grepl('LDL', original) , 'LDL' , grepl('HDL', original) , 'HDL' , grepl('Fatty|Unsaturat|Omega|Docosa|Linoleic', original) , 'Fatty acids' ,
        default = 'Other'
      )
    ]
    nmr[,
      category := factor(
        category,
        levels = c('VLDL', 'IDL', 'LDL', 'HDL', 'Fatty acids', 'Other')
      )
    ]

    fwrite(nmr, sprintf("./results/UKB/factor-profiles/F%d_NMR_Categories.csv", k))
    sets <- lapply(classes, function(g) nmr[category == g, feature])
    names(sets) <- classes
    setorder(nmr, -loading_raw, feature)
    stats <- setNames(nmr$loading_raw, nmr$feature)
    set.seed(20261024 + k)
    z <- as.data.table(fgseaMultilevel(
      pathways = sets,
      stats = stats,
      minSize = 5,
      maxSize = 250,
      scoreType = 'std',
      gseaParam = 1,
      sampleSize = 1001,
      nPermSimple = 100000,
      eps = 0,
      nproc = 1
    ))
    z[, factor := k]
    z[, FDR_within_factor := padj]
    z[, padj := NULL]
    z[,
      leading_edge_features := vapply(
        leadingEdge,
        paste,
        collapse = '; ',
        FUN.VALUE = character(1)
      )
    ]
    z[,
      leading_edge_labels := vapply(
        leadingEdge,
        function(ids) paste(nmr$label[match(ids, nmr$feature)], collapse = '; '),
        character(1)
      )
    ]
    z[, leadingEdge := NULL]
    res[[as.character(k)]] <- z
  }
  r <- rbindlist(res)
  r[, FDR := p.adjust(pval, method = 'BH')]
  r[, direction := ifelse(NES > 0, 'Positive', 'Negative')]
  r[,
    stars := fcase(
      FDR < .001 , '***' , FDR < .01 , '**' , FDR < .05 , '*' ,
      default = ''
    )
  ]

  fwrite(
    r,
    "./results/UKB/factor-profiles/UKB_F2_F5_NMR_Category_Enrichment.csv"
  )
}



# ==============================================================================
# UKB PROTEIN PATHWAYS FOR DISPLAY
# ==============================================================================

# Function to export pathways displayed for the UKB factors.
export_ukb_displayed_pathways <- function() {
  proteins <- fread(
    "./results/cross-cohort/enrichment/UKB/Gene_Mapping.csv",
    select = c("feature", "gene")
  )
  members <- fread(
    "../inputs/gene-sets.csv"
  )[gene %in% proteins$gene]
  sets <- lapply(split(members$gene, members$pathway), unique)
  tests <- fread(
    "./results/cross-cohort/enrichment/UKB/All_Tests.csv"
  )
  tests[,
    collection := fcase(
      collection == "CP:REACTOME"     , "Reactome" ,
      collection == "CP:KEGG_MEDICUS" , "KEGG"     ,
      default = collection
    )
  ]
  for (k in c(2, 5)) {
    eligible <- tests[
      factor == k &
        padj < .05 &
        statistic > 0 &
        collection %in% c("GO:BP", "KEGG", "Reactome", "REACTOME")
    ]
    setorder(eligible, collection, padj, pathway)
    selected <- eligible[,
      {
        keep <- character()
        for (id in pathway) {
          overlap <- vapply(
            keep,
            function(other) {
              length(intersect(sets[[id]], sets[[other]])) /
                length(union(sets[[id]], sets[[other]]))
            },
            numeric(1)
          )
          if (!any(overlap >= .6)) {
            keep <- c(keep, id)
          }
          if (length(keep) == 6L) break
        }
        .SD[pathway %in% keep]
      },
      by = collection
    ]
    selected[, size := lengths(sets[pathway])]
    selected[, ratio := size / nrow(proteins)]
    selected[, direction := "all"]
    label <- tolower(gsub(
      "_",
      " ",
      sub("^(GOBP_|REACTOME_|KEGG_)", "", selected$pathway)
    ))
    label <- paste0(toupper(substr(label, 1, 1)), substring(label, 2))
    label[grepl(
      "Regulation of insulin like growth factor",
      label
    )] <- "IGF transport and uptake by IGFBPs"
    selected[,
      display_label := sub("Medicus reference ", "", label, fixed = TRUE)
    ]
    dir.create(
      "./results/UKB/factor-profiles",
      recursive = TRUE,
      showWarnings = FALSE
    )

    fwrite(
      selected,
      sprintf(
        "./results/UKB/factor-profiles/F%d_Absolute_Displayed_Pathways.csv",
        k
      )
    )
  }
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  cross_absolute_enrichment()
  cross_ukb_absolute_enrichment()
  export_ukb_displayed_pathways()
  cross_signed_enrichment()
  for (k in c(1L, 2L, 3L, 4L, 8L)) {
    vip_factor_enrichment(k)
  }
  ukb_nmr_enrichment()
}
