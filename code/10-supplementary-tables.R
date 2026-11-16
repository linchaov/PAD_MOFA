# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, fgsea 1.38.0, MOFA2 1.22.1
#
# Global: Prepare Supplementary Tables S15–S22 from factor loadings and enrichment results.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds and UKB.rds, fitted MOFA models, enrichment outputs from 06-pathway-enrichment.R, pathway
#         annotations and VIP-PAD enrichment results.
#
# Outputs: Supplementary Tables S15-S22 are printed to the console.
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

# Function to print the table number and all rows and columns.
print_table <- function(id, rows) {
  cat("\n", id, "\n", sep = "")
  # Excel displays P-values in scientific notation and uses General for S18 FDR.
  for (column in intersect(c("P-value", if (id != "S18") "FDR"), names(rows))) {
    set(rows, j = column, value = sprintf("%.2E", rows[[column]]))
  }
  print(rows, nrows = Inf, trunc.cols = FALSE, digits = 15)
  invisible(NULL)
}

# Function to extract gene symbols from VIP-PAD protein identifiers.
vip_symbols <- function(features) {
  genes <- sub("_[^_]+$", "", sub("^PROT_", "", features))
  genes[genes == "Unmapped"] <- ""
  chartr("_", "-", genes)
}



# ==============================================================================
# COMBINE CROSS-COHORT PATHWAY RESULTS | TABLE S15
# ==============================================================================

# Function to assemble the cross-cohort pathway table.
pathways_s15 <- function() {
  # Format the displayed pathways exactly as in the Figure 3A supplementary table.
  # Use measured gene counts from the current fit and pathway identifiers from the study
  # annotations.
  themes <- c(
    'Inflammation & immunity',
    'Complement',
    'Coagulation & platelets',
    'Lipid metabolism',
    'Endothelium & angiogenesis',
    'ECM & remodeling',
    'Metabolism & oxidative stress',
    'Cell injury & repair',
    'Other'
  )

  meta <- fread("../inputs/pathways.csv")
  rows <- list()
  for (co in c("VIP-PAD", "UKB")) {
    folder <- file.path(
      './results/cross-cohort/enrichment',
      co
    )
    a <- fread(file.path(folder, "Displayed_Terms.csv"))
    ix <- match(a$pathway, meta$pathway)
    stopifnot(!anyNA(ix))
    coverage <- fread(file.path(folder, "Gene_Set_Coverage.csv"))
    url <- meta$gs_url[ix]
    pid <- ifelse(
      grepl('GO:[0-9]+', url),
      sub('.*(GO:[0-9]+).*', '\\1', url),
      ifelse(
        grepl('R-HSA-[0-9]+', url),
        sub('.*(R-HSA-[0-9]+).*', '\\1', url),
        ifelse(
          grepl('/entry/N[0-9]+', url),
          sub('.*/entry/(N[0-9]+).*', '\\1', url),
          NA_character_
        )
      )
    )
    b <- data.table(
      Cohort = co,
      Factor = paste('Factor', a$factor),
      Biological_theme = a$theme,
      Database = sub('^CP:', '', a$collection),
      Pathway_ID = pid,
      Pathway_name = gsub(
        '_',
        ' ',
        sub('^(GOBP_|GOCC_|GOMF_|REACTOME_|KEGG_MEDICUS_)', '', a$pathway)
      ),
      Gene_count = as.integer(coverage$measured_size[match(
        a$pathway,
        coverage$pathway
      )]),
      P_value = a$pval,
      FDR = a$padj,
      MSigDB_gene_set = a$pathway,
      factor_order = a$factor,
      theme_order = match(a$theme, themes)
    )
    setorder(b, factor_order, theme_order, FDR, MSigDB_gene_set)
    rows[[co]] <- b
  }
  print_table("S15", rbindlist(rows)[, .(
    Cohort, Factor, `Biological theme` = Biological_theme, Database,
    `Pathway ID` = Pathway_ID, `Pathway name` = Pathway_name,
    `Gene count` = Gene_count, `P-value` = P_value, FDR
  )])
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



# ==============================================================================
# UKB NMR CLASS ENRICHMENT | TABLE S22
# ==============================================================================

# Function to test enrichment of NMR classes in UKB factors.
ukb_nmr_enrichment <- function() {
  cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  ann <- as.data.table(attr(cohort, "feature_dictionary"))
  fit <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
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
    z <- as.data.table(fgsea::fgseaMultilevel(
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
# UKB PROTEIN PATHWAYS | TABLE S21
# ==============================================================================

# Function to prepare UKB protein pathway enrichment for Table S21.
protein_s21 <- function() {
  # UKB protein enrichment for Factors 2 and 5.

  coverage <- fread(
    "./results/cross-cohort/enrichment/UKB/Gene_Set_Coverage.csv"
  )
  features <- fread(
    "./results/cross-cohort/enrichment/UKB/Gene_Mapping.csv"
  )
  tests <- fread(
    "./results/cross-cohort/enrichment/UKB/All_Tests.csv"
  )
  ukb_protein <- tests[factor %in% c(2, 5) & padj < .05 & statistic > 0]
  ukb_protein[,
    count := coverage$measured_size[match(pathway, coverage$pathway)]
  ]
  ukb_protein[, gene_ratio := count / nrow(features)]
  name <- sub("^(GOBP_|GOCC_|GOMF_|REACTOME_|KEGG_)", "", ukb_protein$pathway)
  name <- tolower(gsub("_", " ", name, fixed = TRUE))
  ukb_protein[, description := paste0(toupper(substr(name, 1, 1)), substring(name, 2))]
  table <- ukb_protein[, .(
    Factor = paste0("F", factor),
    Database = fcase(
      collection == "CP:REACTOME", "Reactome",
      collection == "CP:KEGG_MEDICUS", "KEGG",
      default = sub("^GO:", "GO ", collection)
    ),
    `Pathway ID` = pathway, Description = description,
    Count = count, `Gene ratio` = gene_ratio,
    `P-value` = pval,
    FDR = padj
  )]
  setorderv(table, c("Factor", "Database", "FDR", "Pathway ID"))
  print_table("S21", table)
}



# ==============================================================================
# UKB NMR CATEGORIES | TABLE S22
# ==============================================================================

# Function to prepare UKB NMR enrichment for Table S22.
lipids_s22 <- function() {
  result_file <- "./results/UKB/factor-profiles/UKB_F2_F5_NMR_Category_Enrichment.csv"
  if (!file.exists(result_file)) ukb_nmr_enrichment()
  x <- fread(result_file)
  print_table(
    "S22",
    x[
      FDR < .05,
      .(
        Factor = paste0("F", factor),
        Category = pathway, Direction = direction, `Set size` = size,
        NES,
        `P-value` = pval,
        FDR,
        `Leading-edge measures` = leading_edge_labels
      )
    ]
  )
}



# ==============================================================================
# TOP PROTEIN AND METABOLITE LOADINGS | TABLES S16, S17 AND S20
# ==============================================================================

# Function to prepare the top molecular loadings by factor.
loading_tables <- function() {
  # Supplementary Tables S16, S17 and S20. Use the original within-view
  # loading ranks for the factors shown in the tables.

  for (cohort in c("VIP-PAD", "UKB")) {
    vip <- cohort == "VIP-PAD"
    base <- file.path("./results", cohort)
    fit <- readRDS(file.path(base, "mofa/Selected_Fit.rds"))
    if (vip) {
      metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
    } else {
      cohort_data <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
      ukb_labels <- as.data.table(attr(cohort_data, "feature_dictionary"))
      ukb_labels[,
        display_label := ifelse(
          modality == "Protein",
          toupper(label),
          pretty_nmr(label)
        )
      ]
    }
    factors <- if (vip) c(1, 2, 3, 4, 8) else c(2, 5)
    for (view in c("Protein", "Metabolite")) {
      # Rank molecular features by their loading within each factor.
      ranks <- rbindlist(lapply(factors, function(k) {
        weights <- fit$W[[view]][, k]
        rows <- data.table(
          feature = names(weights),
          loading = as.numeric(weights)
        )
        if (vip) {
          if (view == "Protein") {
            rows[, label := vip_symbols(feature)]
          } else {
            rows[, label := metabolites$label[match(feature, metabolites$feature)]]
          }
        } else {
          rows[,
            label := ukb_labels$display_label[match(
              feature,
              ukb_labels$feature
            )]
          ]
          if (view == "Protein") rows[, label := toupper(label)]
        }
        rows[is.na(label) | label == "", label := feature]
        rows[, abs_weight := abs(loading)]
        setorder(rows, -abs_weight, feature)
        rows[
          1:50,
          .(
            factor = paste0("F", k),
            view,
            feature,
            label,
            loading,
            absolute_rank = 1:50
          )
        ]
      }))
      if (vip) {
        table_id <- if (view == "Protein") "S16" else "S17"
        table <- ranks[, .(Factor = factor, label, Loading = loading, `Absolute rank` = absolute_rank)]
        setnames(table, "label", if (view == "Protein") "Protein name" else "Metabolite name")
        print_table(table_id, table)
      } else if (view == "Protein") {
        proteins <- ranks
      } else {
        ukb <- rbind(proteins, ranks)
        ukb[,
          view := factor(
            view,
            levels = c("Protein", "Metabolite"),
            labels = c("Protein", "NMR")
          )
        ]
        setorder(ukb, factor, view, absolute_rank)
        print_table("S20", ukb[, .(
          Factor = factor, `Omics layer` = view, `Molecule or measure` = label,
          Loading = loading, `Absolute rank` = absolute_rank
        )])
      }
    }
  }
}



# ==============================================================================
# VIP-PAD PROTEIN AND LIPID ENRICHMENT | TABLES S18 AND S19
# ==============================================================================

# Function to prepare VIP-PAD protein and lipid enrichment tables.
vip_enrichment_tables <- function() {
  # Tables S18 and S19 report enriched protein pathways and lipid classes in VIP-PAD.

  sets <- readRDS("../inputs/vip-sets.rds")
  protein <- list()
  lipid <- list()
  for (k in c(1, 2, 3, 4, 8)) {
    result <- readRDS(sprintf("./results/VIP-PAD/enrichment/F%d.rds", k))
    x <- result$enrichment
    # Several pathway names can describe the same measured gene set.
    # Keep one canonical pathway per enriched test.
    protein[[as.character(k)]] <- x[
      FDR < .05 & statistic > 0 & set_id == test_id,
      .(
        Factor = paste0("F", k),
        Database = database, `Pathway ID` = pathway, Description = description,
        Count = size, `Gene ratio` = size / ncol(sets$B),
        `P-value` = P,
        FDR
      )
    ]
    x <- result$lipid_results
    lipid[[as.character(k)]] <- x[
      FDR < .05,
      .(
        Factor = paste0("F", k),
        `Lipid class` = lipid_class, Direction = direction, `Set size` = size,
        NES,
        `P-value` = pval,
        FDR,
        `Leading-edge metabolites` = leadingEdge_names
      )
    ]
  }
  pathways <- rbindlist(protein)
  setorderv(pathways, c("Factor", "Database", "FDR", "Pathway ID"))
  print_table("S18", pathways)
  print_table("S19", rbindlist(lipid))
}



# ==============================================================================
# ASSEMBLE SUPPLEMENTARY TABLES | TABLES S15–S22
# ==============================================================================

# Function to prepare the supplementary tables from the corresponding analyses.
all_tables <- function() {
  loading_tables()
  pathways_s15()
  protein_s21()
  vip_enrichment_tables()
  lipids_s22()
}



# ==============================================================================
# RUN AND PRINT SUPPLEMENTARY TABLES
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  # RStudio Source runs all tables; Rscript also accepts a single table ID.
  args <- if (sys.nframe() == 0L) commandArgs(trailingOnly = TRUE) else character()
  tables <- list(ALL = all_tables, S15 = pathways_s15, S21 = protein_s21, S22 = lipids_s22)
  table <- if (length(args)) toupper(args) else "ALL"
  if (length(table) != 1L || !table %in% names(tables)) {
    stop("Choose one of all, S15, S21 or S22.")
  }
  tables[[table]]()
}
