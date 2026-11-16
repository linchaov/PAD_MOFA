# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, logistf 1.26.1, lpSolve 5.6.23
#
# Global: Baseline characteristics, factor–trait correlations, comorbidity models and cluster
#         comparisons.
#
# Inputs: VIP.rds and UKB.rds in ./Datasets_PAD_MOFA/; fitted factors and cluster assignments
#         from 01-mofa.R.
#
# Outputs: Baseline tables S10 and S11 are printed to the console. Association tables,
#          clinical plot data and Rutherford factor comparisons are saved
#         under ./results/.
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

# Function to compare clinical measurements between clusters using rank tests.
cluster_wilcoxon <- function(d) {
  res <- d[,
    .(
      P = wilcox.test(
        value ~ group,
        alternative = 'two.sided',
        exact = FALSE,
        correct = TRUE
      )$p.value
    ),
    by = variable
  ]
  res[,
    label := fcase(
      P < .001, '***', P < .01, '**', P < .05, '*',
      default = 'ns'
    )
  ]
  res
}

# Function to set out-of-range measurements to missing.
clean_range <- function(x, lo, hi) {
  x[!is.finite(x) | x < lo | x > hi] <- NA_real_
  x
}

# Function to code recorded conditions and treatments from the VIP-PAD fields.
source_binary <- function(value, pattern = NULL) {
  if (is.null(pattern)) return(as.numeric(value != 'No'))
  ifelse(is.na(value), NA_real_, as.numeric(grepl(pattern, value)))
}

# Function to extract clinical fields and harmonize VIP-PAD identifiers and outcome columns.
vip_clinical_data <- function() {
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  attrs <- attr(cohort, "omics_matrix_attributes")
  omics <- unlist(
    lapply(attrs, function(a) names(a[["scaled:scale"]])),
    use.names = FALSE
  )
  clinical <- cohort[, setdiff(names(cohort), omics), drop = FALSE]
  stopifnot(!anyNA(clinical$study_id), !anyDuplicated(clinical$study_id))
  clinical$SampleID <- clinical$study_id
  clinical$FollowUp_Months <- 12 * clinical$MALE_time_years
  clinical$MALE_binary <- clinical$MALE_event
  clinical$MACCE_binary <- clinical$MACCE_event
  clinical
}

# Function to extract UKB clinical fields using the omics feature dictionary.
ukb_clinical_data <- function() {
  cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  omics <- attr(cohort, "feature_dictionary")$source_field
  clinical <- cohort[, setdiff(names(cohort), omics), drop = FALSE]
  stopifnot(!anyNA(clinical$study_id), !anyDuplicated(clinical$study_id))
  clinical
}

# Function to define labels, units and groups for VIP-PAD clinical traits.
vip_trait_dictionary <- function() {
  data.table::fread(
    text = "id|label|group|type|definition
age|Age|Clinical characteristics|Pearson|Age in years
male|Sex (male)|Clinical characteristics|Point-biserial|Male=1, Female=0
rutherford|Rutherford grade|Clinical characteristics|Spearman|Ordered severity 0-6; unknown is missing
barthel|Barthel index|Clinical characteristics|Pearson|Preoperative recorded total, 0-100
coronary|Coronary disease|Comorbidities|Point-biserial|Any explicit CAD/MI/angina mention; other recorded heart disease counts as no coronary history
stroke|Stroke|Comorbidities|Point-biserial|Ischaemic or haemorrhagic stroke explicitly recorded
diabetes|Diabetes|Comorbidities|Point-biserial|Recorded diabetes history
ckd|CKD|Comorbidities|Point-biserial|Any recorded CKD stage versus none
htn|Hypertension|Comorbidities|Point-biserial|Recorded history: present=1, absent=0
copd|COPD|Comorbidities|Point-biserial|Recorded COPD : any severity/use=1; explicit none=0; missing excluded pairwise
anemia|Anemia|Comorbidities|Point-biserial|Recorded Anemia : any severity/use=1; explicit none=0; missing excluded pairwise
dyslipidemia|Dyslipidemia|Comorbidities|Point-biserial|Recorded Dyslipidemia : any severity/use=1; explicit none=0; missing excluded pairwise
hyperhomocysteinemia|Hyperhomocysteinemia|Comorbidities|Point-biserial|Recorded Hyperhomocysteinemia : any severity/use=1; explicit none=0; missing excluded pairwise
hyperuricemia|Hyperuricemia|Comorbidities|Point-biserial|Recorded Hyperuricemia : any severity/use=1; explicit none=0; missing excluded pairwise
antihypertensive|Antihypertensive|Medications|Point-biserial|Recorded preoperative class use (yes=1, no=0)
lipidlowering|Lipid-lowering|Medications|Point-biserial|Recorded LipidLowering : any severity/use=1; explicit none=0; missing excluded pairwise
antithrombotic|Antithrombotic|Medications|Point-biserial|Recorded Antithrombotic : any severity/use=1; explicit none=0; missing excluded pairwise
vasodilator|Vasodilator|Medications|Point-biserial|Recorded Vasodilator : any severity/use=1; explicit none=0; missing excluded pairwise
sbp|SBP|Biomarkers|Pearson|mmHg; inherited clinical QC
dbp|DBP|Biomarkers|Pearson|mmHg; inherited clinical QC
bmi|BMI|Biomarkers|Pearson|Recalculated weight/height^2, inherited QC
glucose|Fasting glucose|Biomarkers|Pearson|mmol/L; original scale, inherited QC
tg|Triglycerides|Biomarkers|Pearson|mmol/L; original scale, inherited QC
hdl|HDL-C|Biomarkers|Pearson|mmol/L; inherited QC
ldl|LDL-C|Biomarkers|Pearson|mmol/L; inherited QC
egfr|eGFR|Biomarkers|Pearson|eGFR calculated using the 2021 CKD-EPI equation
crp|CRP|Biomarkers|Pearson|mg/L; original scale, inherited QC
hemoglobin|Hemoglobin|Biomarkers|Pearson|Recorded g/L
albumin|Albumin|Biomarkers|Pearson|Albumin g/L; baseline table QC 10-70; Pearson",
    sep = "|"
  )
}

# Function to align the two molecular clusters with the clinical rows by participant ID.
cluster_groups <- function(ids, assignments) {
  group <- assignments$K2[match(ids, assignments$sample_id)]
  stopifnot(!anyDuplicated(assignments$sample_id), setequal(group, 1:2))
  group
}



# ==============================================================================
# UKB CLINICAL MEASURES BY CLUSTER | FIGURE 3 (PANEL F)
# ==============================================================================

# Function to summarize clinical measurements by UKB molecular cluster.
ukb_clinical_distributions <- function() {
  w <- as.data.table(ukb_clinical_data())
  a <- fread("./results/UKB/clustering/Patient_Assignments.csv")
  ids <- w$study_id
  group <- cluster_groups(ids, a)
  # Preserve the definitions and native units in the UKB baseline table.
  vars <- c(
    Age = 'Age_i0',
    BMI = 'BMI_i0',
    SBP = 'SBP_i0_PMBP',
    LDL = 'LDLDirect_i0_Blood',
    eGFR = 'eGFR_calculated'
  )
  # Reshape each clinical trait for the cluster comparison.
  d <- rbindlist(lapply(names(vars), function(v) {
    data.table(
      sample_id = ids,
      group = factor(paste0('C', group), levels = c('C1', 'C2')),
      variable = v,
      value = w[[vars[v]]]
    )
  }))
  d <- d[!is.na(value)]
  d[, x := as.numeric(group)]
  st <- d[,
    .(
      N = .N,
      mean = mean(value),
      SD = sd(value),
      median = median(value),
      Q1 = unname(quantile(value, .25)),
      Q3 = unname(quantile(value, .75)),
      min = min(value),
      max = max(value)
    ),
    by = .(variable, group)
  ]
  dir.create(
    "./results/UKB/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    d,
    "./results/UKB/clinical/Distributions.rds"
  )
  fwrite(
    st,
    "./results/UKB/clinical/Descriptive_Statistics.csv"
  )
  fwrite(cluster_wilcoxon(d), "./results/UKB/clinical/Wilcoxon_Tests.csv")
}



# ==============================================================================
# UKB FACTOR–TRAIT ASSOCIATIONS | FIGURE S08 (PANELS A–B)
# ==============================================================================

# Function to correlate UKB factor scores with clinical traits.
ukb_factor_clinical_associations <- function() {
  # Pearson correlations, pairwise complete observations, BH across 275 tests.
  clinical <- as.data.table(ukb_clinical_data())
  fit <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
  ids <- rownames(fit$Z)
  stopifnot(
    !anyDuplicated(ids),
    setequal(ids, clinical$study_id),
    all(is.finite(fit$Z))
  )
  d <- clinical[match(ids, study_id)]
  # Match clinical records to factor scores by participant ID.
  for (k in seq_len(ncol(fit$Z))) {
    set(d, j = paste0('F', k), value = fit$Z[, k])
  }

  # Clinical phenotypes: demographics, eight comorbidities and biomarkers.
  spec <- data.table(
    variable = c(
      'Age_i0',
      'BMI_i0',
      'SBP_i0_PMBP',
      'DBP_i0_PMBP',
      'Glucose_i0_Blood',
      'HbA1c_i0_Blood',
      'Cholesterol_i0_Blood',
      'LDLDirect_i0_Blood',
      'HDLCholesterol_i0_Blood',
      'Triglycerides_i0_Blood',
      'eGFR_calculated',
      'Creatinine_i0_Blood',
      'CystatinC_i0_Blood',
      'CRP_i0_Blood',
      'Hemoglobin_gL',
      'Albumin_i0_Blood'
    ),
    label = c(
      'Age',
      'BMI',
      'SBP',
      'DBP',
      'Glucose',
      'HbA1c',
      'Total cholesterol',
      'LDL-C',
      'HDL-C',
      'Triglycerides',
      'eGFR',
      'Creatinine',
      'Cystatin C',
      'CRP',
      'Hemoglobin',
      'Albumin'
    ),
    unit = c(
      'years',
      'kg/m²',
      'mmHg',
      'mmHg',
      'mmol/L',
      '%',
      'mmol/L',
      'mmol/L',
      'mmol/L',
      'mmol/L',
      'mL/min/1.73 m²',
      'µmol/L',
      'mg/L',
      'mg/L',
      'g/L',
      'g/L'
    ),
    group = c(
      rep('Demographics', 2),
      rep('Blood pressure', 2),
      rep('Glycemia', 2),
      rep('Lipids', 4),
      rep('Renal function', 3),
      rep('Other biomarkers', 3)
    )
  )
  # Missing binary traits remain missing.
  comorb <- c(
    'Hypertension',
    'Dyslipidemia',
    'Diabetes',
    'CoronaryHeartDisease',
    'Stroke',
    'ChronicKidneyDisease',
    'HeartFailure',
    'COPD'
  )
  binary_vars <- paste0('Baseline_', comorb, '_i0')
  d[, Sex_male := as.numeric(as.character(Sex) == 'Male')]
  for (i in seq_along(binary_vars)) {
    value <- as.character(d[[binary_vars[i]]])
    set(d, j = binary_vars[i], value = as.numeric(value == 'Yes'))
  }
  sex <- data.table(
    variable = 'Sex_male',
    label = 'Sex (male)',
    unit = 'Male=1; Female=0',
    group = 'Clinical'
  )
  com <- data.table(
    variable = binary_vars,
    label = c(
      'Hypertension',
      'Dyslipidemia',
      'Diabetes',
      'CHD',
      'Stroke',
      'CKD',
      'Heart failure',
      'COPD'
    ),
    unit = 'Yes=1; No=0',
    group = 'Comorbidities'
  )
  spec <- rbindlist(list(spec[1], sex, spec[2], com, spec[3:16]))
  spec[1:3, group := 'Clinical']
  spec[12:25, group := 'Biomarkers']
  spec[,
    type := ifelse(
      variable %in% c('Sex_male', binary_vars),
      'Point-biserial',
      'Pearson'
    )
  ]
  spec[, source_column := ifelse(variable == 'Sex_male', 'Sex', variable)]
  # Summarize missing and outlying clinical measurements.
  qc <- rbindlist(lapply(spec$variable, function(v) {
    x <- d[[v]]
    a <- x[is.finite(x)]
    q <- quantile(a, c(.25, .75))
    iq <- q[2] - q[1]
    data.table(
      variable = v,
      N = length(a),
      missing_n = sum(is.na(x)),
      missing_percent = 100 * mean(is.na(x)),
      min = min(a),
      Q1 = q[1],
      median = median(a),
      Q3 = q[2],
      max = max(a),
      outside_3IQR = sum(a < q[1] - 3 * iq | a > q[2] + 3 * iq)
    )
  }))
  spec <- merge(spec, qc, by = 'variable', sort = FALSE)[match(
    spec$variable,
    variable
  )]
  spec[type == 'Point-biserial', outside_3IQR := NA_integer_]
  spec[,
    positive_n := vapply(
      variable,
      function(x) {
        if (x %in% c('Sex_male', binary_vars)) {
          sum(d[[x]] == 1, na.rm = TRUE)
        } else {
          NA_integer_
        }
      },
      integer(1)
    )
  ]
  spec[, source_file := "./Datasets_PAD_MOFA/UKB.rds"]
  spec[,
    notes := 'Baseline clinical data; original scale; all finite observed values retained'
  ]
  spec[
    type == 'Point-biserial',
    notes := paste('Baseline binary coding:', unit)
  ]
  spec[
    variable == 'eGFR_calculated',
    notes := 'CKD-EPI 2021 creatinine equation; verified against baseline age and sex'
  ]
  spec[
    variable == 'Hemoglobin_gL',
    notes := 'Clinical hemoglobin supplied in g/L in the prepared cohort'
  ]
  spec[
    variable == 'HbA1c_i0_Blood',
    notes := 'Baseline value expressed as a percentage'
  ]
  spec[
    variable == 'CystatinC_i0_Blood',
    notes := 'Clinical assay, mg/L: UKB field 30720; https://biobank.ndph.ox.ac.uk/showcase/field.cgi?id=30720'
  ]

  # Test each factor against the clinical traits with available measurements.
  results <- rbindlist(lapply(seq_len(ncol(fit$Z)), function(k) {
    rbindlist(lapply(spec$variable, function(v) {
      ok <- is.finite(d[[v]]) & is.finite(d[[paste0('F', k)]])
      z <- cor.test(d[[paste0('F', k)]][ok], d[[v]][ok], method = 'pearson')
      data.table(
        factor = paste0('F', k),
        variable = v,
        N = sum(ok),
        r = unname(z$estimate),
        CI_lower = z$conf.int[1],
        CI_upper = z$conf.int[2],
        P = z$p.value
      )
    }))
  }))
  results[, FDR := p.adjust(P, 'BH')]
  results[,
    stars := fifelse(
      FDR < .001,
      '***',
      fifelse(FDR < .01, '**', fifelse(FDR < .05, '*', ''))
    )
  ]
  results <- merge(
    results,
    spec[, .(variable, label, unit, type)],
    by = 'variable',
    sort = FALSE
  )
  results[, factor_order := as.integer(sub('F', '', factor))]
  results[, variable_order := match(variable, spec$variable)]
  setorder(results, factor_order, variable_order)
  # Save clinical values and association estimates for Figure S08.
  dir.create(
    "./results/UKB/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    list(
      data = as.data.frame(d[,
        c(paste0('F', seq_len(ncol(fit$Z))), spec$variable),
        with = FALSE
      ]),
      associations = results,
      variables = spec
    ),
    "./results/UKB/clinical/Factor_Associations.rds"
  )
}

# Function to use Monte Carlo sampling when the exact RxC Fisher test fails.
fisher_rxc <- function(tt) {
  tryCatch(
    fisher.test(tt, workspace = 2e7),
    error = function(e) fisher.test(tt, simulate.p.value = TRUE, B = 100000)
  )
}



# ==============================================================================
# BASELINE CHARACTERISTICS | TABLES S10 AND S11
# ==============================================================================

# Function to print baseline characteristics by molecular cluster in each cohort.
baseline_table <- function(cohort) {
  set.seed(42)
  invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
  w <- if (cohort == 'UKB') ukb_clinical_data() else vip_clinical_data()
  assignment <- fread(file.path('./results', cohort, 'clustering/Patient_Assignments.csv'))
  group <- cluster_groups(w$study_id, assignment)
  rows <- list()
  # Function to count or format missing observations for the baseline table.
  miss <- function(x) {
    sprintf('%d (%.1f%%)', sum(is.na(x)), 100 * mean(is.na(x)))
  }
  # Function to format P-values for the baseline table.
  fmtp <- function(p) {
    if (is.na(p)) {
      'NA'
    } else if (p < .001) {
      '<0.001'
    } else {
      sprintf('%.3f', p)
    }
  }
  # Function to append a formatted characteristic to the baseline table.
  add <- function(
    label,
    value = '',
    missing = '',
    c1 = '',
    c2 = '',
    p = ''
  ) {
    if (label == 'Sex (male)' && c1 == '') {
      x <- w$Sex == 'Male'
      c1 <- pct(x[group == 1])
      c2 <- pct(x[group == 2])
      p <- test(x)
    }
    rows[[length(rows) + 1]] <<- data.table(
      Characteristic = label,
      Overall = value,
      Missing = missing,
      Cluster1 = c1,
      Cluster2 = c2,
      P = p
    )
  }
  # Function to format counts and percentages for categorical characteristics.
  pct <- function(x) {
    if (!sum(!is.na(x))) {
      return('NA')
    }
    sprintf('%d (%.1f%%)', sum(x, na.rm = TRUE), 100 * mean(x, na.rm = TRUE))
  }
  # Function to compare a baseline characteristic between clusters.
  test <- function(x, continuous = FALSE) {
    ok <- !is.na(x)
    p <- NA_real_
    a <- x[ok]
    g <- factor(group[ok], levels = 1:2)
    if (continuous) {
      if (length(unique(a)) > 1 && all(table(g) > 1)) {
        p <- wilcox.test(
          a ~ g,
          alternative = 'two.sided',
          exact = FALSE,
          correct = TRUE
        )$p.value
      }
    } else {
      tt <- table(a, g)
      tt <- tt[rowSums(tt) > 0, , drop = FALSE]
      if (nrow(tt) > 1 && all(colSums(tt) > 0)) {
        if (nrow(tt) == 2L) {
          p <- fisher.test(tt)$p.value
        } else {
          p <- fisher_rxc(tt)$p.value
        }
      }
    }
    fmtp(p)
  }
  # Function to summarize a continuous characteristic and compare clusters.
  cont <- function(
    v,
    label,
    d = 1,
    force_median = FALSE,
    force_mean = FALSE
  ) {
    x <- w[[v]]
    z <- x[is.finite(x)]
    sk <- mean((z - mean(z))^3) / mean((z - mean(z))^2)^1.5
    med <- !force_mean && (force_median || abs(sk) > 1)
    # Function to format the within-group distribution for this characteristic.
    summarize <- function(a) {
      a <- a[is.finite(a)]
      if (!length(a)) {
        return('NA')
      }
      # Function to format a value to the specified decimal precision.
      f <- function(y) formatC(y, digits = d, format = 'f')
      if (med) {
        paste0(
          f(median(a)),
          ' [',
          f(quantile(a, .25)),
          ', ',
          f(quantile(a, .75)),
          ']'
        )
      } else {
        paste0(f(mean(a)), ' ± ', f(sd(a)))
      }
    }
    add(
      label,
      summarize(x),
      miss(x),
      c1 = summarize(x[group == 1]),
      c2 = summarize(x[group == 2]),
      p = test(x, continuous = TRUE)
    )
  }
  # Function to tabulate all levels of a categorical characteristic.
  catvar <- function(
    v,
    label,
    levels,
    labels = levels
  ) {
    x <- trimws(as.character(w[[v]]))
    add(label, '', miss(x), p = test(x))
    for (i in seq_along(levels)) {
      y <- x == levels[i]
      add(
        paste0('    ', labels[i]),
        pct(y),
        c1 = pct(y[group == 1]),
        c2 = pct(y[group == 2])
      )
    }
  }
  if (cohort == 'UKB') {
    # Table S11. UKB baseline characteristics.
    # Function to convert recorded binary values to analysis indicators.
    bin <- function(v, label, positive = 'Yes') {
      x <- as.character(w[[v]])
      y <- x == positive
      add(
        label,
        pct(y),
        miss(y),
        c1 = pct(y[group == 1]),
        c2 = pct(y[group == 2]),
        p = test(y)
      )
    }
    add('Demographics and lifestyle')
    cont('Age_i0', 'Age, years', force_mean = TRUE)
    bin('Sex', 'Sex (male)', positive = 'Male')
    bin('Ethnic_White_i0', 'White ethnicity')
    cont('BMI_i0', 'BMI, kg/m²', 1, TRUE)
    bin('Education_College_i0', 'Education (College / university degree or above)')
    cont('TownsendDeprivationIndex', 'Townsend deprivation index', 2, force_mean = TRUE)
    catvar('SmokingStatus_i0', 'Smoking status',
           c('Never', 'Previous', 'Current'), c('Never', 'Former', 'Current'))
    bin('Alcohol_Heavy_i0', 'Heavy alcohol consumption')
    bin('PAshort_MV_i0', 'Healthy physical activity')
    add('PAD characteristics')
    cont('PAD_duration_years', 'PAD duration, years', 1, TRUE)
    bin('Baseline_PriorPADProcedure_i0', 'Prior PAD-related procedure')
    add('Prevalent comorbidities')
    labels <- c(
      Hypertension = 'Hypertension',
      Dyslipidemia = 'Dyslipidemia',
      Diabetes = 'Diabetes mellitus',
      CoronaryHeartDisease = 'Coronary heart disease',
      Stroke = 'Stroke',
      ChronicKidneyDisease = 'Chronic kidney disease',
      HeartFailure = 'Heart failure',
      COPD = 'COPD'
    )
    for (v in names(labels)) {
      bin(paste0('Baseline_', v, '_i0'), labels[[v]])
    }
    add('Medications')
    labels <- c(
      AntiHypertensive = 'Antihypertensive medication',
      AntiHypertensive_Beta = '    Beta-blocker',
      AntiHypertensive_CCB = '    Calcium channel blocker',
      AntiHypertensive_ACEI = '    ACE inhibitor',
      AntiHypertensive_ARB = '    Angiotensin receptor blocker',
      AntiHypertensive_Diuretic = '    Diuretic',
      LipidLowering = 'Lipid-lowering medication',
      LipidLowering_Statin = '    Statin use'
    )
    for (v in names(labels)) {
      bin(paste0('Baseline_', v, '_i0'), labels[[v]])
    }
    bin('Other_lipid_lowering', '    Other lipid-lowering (without statin)')
    labels <- c(
      AntiPlatelet = '    Antiplatelet',
      AntiCoagulation = '    Anticoagulant',
      Vasodilator = 'Vasodilator'
    )
    for (v in names(labels)) {
      bin(paste0('Baseline_', v, '_i0'), labels[[v]])
    }
    add('Blood pressure and biomarkers')
    labels <- c(
      SBP_i0_PMBP = 'SBP, mmHg',
      DBP_i0_PMBP = 'DBP, mmHg',
      Glucose_i0_Blood = 'Glucose, mmol/L',
      HbA1c_i0_Blood = 'HbA1c, %',
      Cholesterol_i0_Blood = 'Total cholesterol, mmol/L',
      Triglycerides_i0_Blood = 'Triglycerides, mmol/L',
      HDLCholesterol_i0_Blood = 'HDL-C, mmol/L',
      LDLDirect_i0_Blood = 'LDL-C, mmol/L',
      Creatinine_i0_Blood = 'Creatinine, µmol/L',
      eGFR_calculated = 'eGFR, mL/min/1.73 m²',
      CRP_i0_Blood = 'CRP, mg/L',
      Hemoglobin_gL = 'Hemoglobin, g/L',
      Albumin_i0_Blood = 'Albumin, g/L'
    )
    for (v in names(labels)) {
      cont(
        v,
        labels[[v]],
        if (grepl('Cholesterol|LDL|Glucose|Triglycerides|CRP', v)) 2 else 1,
        force_mean = v %in% c('SBP_i0_PMBP', 'DBP_i0_PMBP')
      )
    }

    cat("\nS11. UKB baseline characteristics\n")
  } else {
    # Table S10. VIP-PAD baseline characteristics.
    # Function to summarize a binary characteristic by cluster.
    binary <- function(v, label, pattern = NULL) {
      x <- trimws(as.character(w[[v]]))
      y <- if (is.null(pattern)) x != 'No' else grepl(pattern, x)
      y[is.na(x)] <- NA
      add(
        label,
        pct(y),
        miss(y),
        c1 = pct(y[group == 1]),
        c2 = pct(y[group == 2]),
        p = test(y)
      )
    }
    add('Demographics')
    cont('Age', 'Age, years')
    add(
      'Sex (male)',
      sprintf('%d (%.1f%%)', sum(w$Sex == 'Male'), 100 * mean(w$Sex == 'Male')),
      miss(w$Sex)
    )
    cont('BMI', 'BMI, kg/m²')
    catvar(
      'SmokingStatus',
      'Smoking status',
      c('Non-smoker', 'Former smoker', 'Current smoker')
    )
    add('PAD characteristics')
    cont(
      'DiseaseDuration_Years_Calc',
      'Self-reported disease duration, years',
      1,
      TRUE
    )
    catvar(
      'PriorPADProcedure',
      'Prior PAD-related procedure',
      c(
        'No documented prior PAD procedure',
        'Prior PAD revascularization',
        'Prior PAD-related amputation',
        'Prior PAD arterial surgery, type unspecified'
      ),
      c('None', 'Revascularization', 'Amputation', 'Other PAD arterial surgery')
    )
    catvar(
      'RutherfordCategory',
      'Rutherford category',
      0:6,
      c(
        '0: Asymptomatic',
        '1: Mild claudication',
        '2: Moderate claudication',
        '3: Severe claudication',
        '4: Rest pain',
        '5: Minor tissue loss',
        '6: Major tissue loss'
      )
    )
    # Group disease course at admission for Table S10.
    w$DiseaseCourse_Grouped <- w$DiseaseCourse
    w$DiseaseCourse_Grouped[w$DiseaseCourse %in% c('Acute', 'Acute-on-chronic')] <-
      'Acute onset or exacerbation'
    cont('ClaudicationDistance', 'Preoperative claudication distance, m', 0, TRUE)
    cont('BarthelTotal', 'Barthel index', 0, TRUE)
    cont('ABI_Left', 'Left ABI', 2)
    cont('ABI_Right', 'Right ABI', 2)
    cont('LVEF', 'LVEF, %', 1)
    catvar(
      'DiseaseCourse_Grouped',
      'Disease course at this admission',
      c('Acute onset or exacerbation', 'Subacute', 'Chronic')
    )
    catvar(
      'IndexProcedureCategory',
      'Procedure during this admission',
      c(
        'No procedure documented',
        'Endovascular',
        'Open surgery',
        'Hybrid',
        'Procedure, approach unrecorded'
      ),
      c(
        'No procedure',
        'Endovascular',
        'Open surgery',
        'Hybrid',
        'Procedure, approach unrecorded'
      )
    )
    add('Prevalent comorbidities')
    binary('Hypertension', 'Hypertension')
    binary('DiabetesCategory', 'Diabetes mellitus')
    binary('HeartDisease', 'Coronary artery disease', 'Coronary artery disease|Myocardial infarction|Angina')
    binary('CerebrovascularDisease', 'Stroke', 'stroke')
    binary('CKDStage_Recorded', 'Chronic kidney disease')
    binary('Dyslipidemia', 'Dyslipidemia')
    add('Medications')
    for (v in c(
      'Antihypertensive',
      'LipidLowering',
      'Antithrombotic',
      'Vasodilator'
    )) {
      binary(
        v,
        switch(
          v,
          Antihypertensive = 'Antihypertensive',
          LipidLowering = 'Lipid-lowering',
          Antithrombotic = 'Antithrombotic',
          Vasodilator = 'Vasodilator'
        )
      )
    }
    add('Blood pressure and biomarkers')
    for (v in c(
      'SBP',
      'DBP',
      'FastingGlucose',
      'Triglycerides',
      'HDLCholesterol',
      'LDLCholesterol',
      'eGFR_CKD_EPI_2021',
      'CRP',
      'Hemoglobin',
      'Albumin'
    )) {
      lab <- switch(
        v,
        SBP = 'SBP, mmHg',
        DBP = 'DBP, mmHg',
        FastingGlucose = 'Fasting glucose, mmol/L',
        Triglycerides = 'Triglycerides, mmol/L',
        HDLCholesterol = 'HDL-C, mmol/L',
        LDLCholesterol = 'LDL-C, mmol/L',
        eGFR_CKD_EPI_2021 = 'eGFR, mL/min/1.73 m²',
        CRP = 'CRP, mg/L',
        Hemoglobin = 'Hemoglobin, g/L',
        Albumin = 'Albumin, g/L'
      )
      cont(
        v,
        lab,
        if (
          v %in%
            c(
              'FastingGlucose',
              'Triglycerides',
              'HDLCholesterol',
              'LDLCholesterol',
              'CRP'
            )
        ) {
          2
        } else {
          1
        }
      )
    }

    cat("\nS10. VIP-PAD baseline characteristics\n")
  }
  tab <- rbindlist(rows)
  tab[, Characteristic := trimws(Characteristic)]
  if (cohort == "UKB") {
    # Report a common missing count once per section, as in Table S11.
    for (heading in c("Prevalent comorbidities", "Medications")) {
      first <- match(heading, tab$Characteristic)
      last <- match(if (heading == "Medications") "Vasodilator" else "COPD", tab$Characteristic)
      values <- tab$Missing[(first + 1L):last]
      if (length(unique(values)) == 1L) {
        tab[first, Missing := values[1L]]
        tab[(first + 1L):last, Missing := ""]
      }
    }
  }
  setnames(tab, c(
    "Characteristic", sprintf("Overall (N = %d)", nrow(w)), "Missing, n (%)",
    sprintf("Cluster 1 (N = %d)", sum(group == 1, na.rm = TRUE)),
    sprintf("Cluster 2 (N = %d)", sum(group == 2, na.rm = TRUE)), "P-value"
  ))
  print(tab, nrows = Inf, trunc.cols = FALSE)
}



# ==============================================================================
# VIP-PAD FACTOR–TRAIT ASSOCIATIONS | FIGURE 2 (PANEL A)
# ==============================================================================

# Function to relate VIP-PAD factor scores to clinical traits.
vip_clinical <- function() {
  set.seed(42)
  # Pearson/point-biserial; Spearman for Rutherford; BH across nrow(traits)*K tests.
  traits <- vip_trait_dictionary()
  # Function to assemble clinical trait values for factor associations.
  clinical_values <- function(w) {
    sbp <- clean_range(w$SBP, 50, 300)
    dbp <- clean_range(w$DBP, 30, 180)
    bad <- which(!is.na(sbp) & !is.na(dbp) & sbp <= dbp)
    sbp[bad] <- NA_real_
    dbp[bad] <- NA_real_
    x <- list(
      age = w$Age,
      male = as.numeric(w$Sex == 'Male'),
      rutherford = w$RutherfordCategory,
      barthel = w$BarthelTotal,
      coronary = source_binary(w$HeartDisease, 'Coronary artery disease|Myocardial infarction|Angina'),
      stroke = source_binary(
        w$CerebrovascularDisease,
        'stroke'
      ),
      diabetes = source_binary(w$DiabetesCategory),
      ckd = source_binary(w$CKDStage_Recorded),
      htn = source_binary(w$Hypertension),
      antihypertensive = source_binary(w$Antihypertensive),
      sbp = sbp,
      dbp = dbp,
      bmi = clean_range(w$BMI, 10, 70),
      glucose = clean_range(w$FastingGlucose, 1, 50),
      tg = clean_range(w$Triglycerides, .1, 30),
      hdl = clean_range(w$HDLCholesterol, .1, 5),
      ldl = clean_range(w$LDLCholesterol, .1, 15),
      egfr = w$eGFR_CKD_EPI_2021,
      crp = clean_range(w$CRP, .01, 500),
      hemoglobin = w$Hemoglobin,
      albumin = clean_range(as.numeric(w$Albumin), 10, 70)
    )
    for (n in c(
      'COPD',
      'Anemia',
      'Dyslipidemia',
      'Hyperhomocysteinemia',
      'Hyperuricemia',
      'LipidLowering',
      'Antithrombotic',
      'Vasodilator'
    )) {
      x[[tolower(n)]] <- source_binary(w[[n]])
    }
    x[traits$id]
  }
  dir.create(
    "./results/VIP-PAD/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )
  clinical <- vip_clinical_data()
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  Z <- fit$Z[as.character(clinical$study_id), , drop = FALSE]
  K <- ncol(Z)
  colnames(Z) <- sprintf('F%02d', seq_len(K))
  variance <- fit$variance$r2_per_factor$PAD
  fm <- data.table(
    factor = colnames(Z),
    Protein = as.numeric(variance[, "Protein"]),
    Metabolite = as.numeric(variance[, "Metabolite"])
  )
  values <- clinical_values(clinical)
  stopifnot(all(lengths(values) == nrow(clinical)))
  trait_N <- vapply(values, function(x) sum(is.finite(x)), integer(1))
  traits[, N := as.integer(trait_N[id])]
  traits[,
    n_positive := vapply(
      id,
      function(key) {
        if (type[match(key, id)] == "Point-biserial") {
          sum(values[[key]] == 1, na.rm = TRUE)
        } else {
          NA_integer_
        }
      },
      integer(1)
    )
  ]
  # Match each trait to its prespecified correlation test.
  r <- rbindlist(lapply(seq_len(nrow(traits)), function(i) {
    x <- values[[traits$id[i]]]
    # Correlate this clinical trait with each factor score.
    rbindlist(lapply(seq_len(K), function(k) {
      y <- Z[, k]
      ok <- is.finite(x) & is.finite(y)
      method <- if (traits$id[i] == 'rutherford') 'spearman' else 'pearson'
      test <- cor.test(x[ok], y[ok], method = method, exact = FALSE)
      data.table(
        id = traits$id[i],
        factor = fm$factor[k],
        N = sum(ok),
        estimate = unname(test$estimate),
        P = test$p.value
      )
    }))
  }))
  r[, FDR := p.adjust(P, 'BH')]
  r[,
    stars := ifelse(
      FDR < .001,
      '***',
      ifelse(FDR < .01, '**', ifelse(FDR < .05, '*', ''))
    )
  ]

  fwrite(
    fm,
    "./results/VIP-PAD/clinical/Variance_Explained.csv"
  )
  associations <- merge(r, traits[, .(id, label, group, type)], by = 'id', sort = FALSE)
  fwrite(associations, "./results/VIP-PAD/clinical/Figure_2_Updated_Associations.csv")
  saveRDS(traits, "./results/VIP-PAD/clinical/Clinical_Traits_Updated.rds")
}



# ==============================================================================
# VIP-PAD CLINICAL MEASURES BY CLUSTER | FIGURE 2 (PANEL E)
# ==============================================================================

# Function to prepare VIP-PAD clinical measurements by cluster.
vip_prepare_clinical_data <- function() {
  # Clinical distributions by cluster; two-sided Wilcoxon rank-sum tests.
  # Summarize available measurements separately for each trait.
  w <- as.data.table(vip_clinical_data())
  cl <- fread(
    "./results/VIP-PAD/clustering/Patient_Assignments.csv"
  )
  g <- factor(paste0('C', cluster_groups(w$SampleID, cl)), levels = c('C1', 'C2'))
  sbp <- clean_range(w$SBP, 50, 300)
  dbp <- clean_range(w$DBP, 30, 180)
  bad <- which(!is.na(sbp) & !is.na(dbp) & sbp <= dbp)
  sbp[bad] <- NA_real_
  values <- list(
    Age = w$Age,
    BMI = clean_range(w$BMI, 10, 70),
    SBP = sbp,
    LDL = clean_range(w$LDLCholesterol, .1, 15),
    eGFR = w$eGFR_CKD_EPI_2021
  )
  d <- rbindlist(lapply(names(values), function(v) {
    data.table(variable = v, group = g, value = values[[v]])
  }))
  d <- d[is.finite(value)]
  d[, x := as.numeric(group)]
  st <- d[,
    .(
      N = .N,
      median = median(value),
      Q1 = quantile(value, .25),
      Q3 = quantile(value, .75),
      minimum = min(value),
      maximum = max(value)
    ),
    by = .(variable, group)
  ]
  dir.create(
    "./results/VIP-PAD/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(
    st,
    "./results/VIP-PAD/clinical/Descriptive_Statistics.csv"
  )
  saveRDS(
    d,
    "./results/VIP-PAD/clinical/Distributions.rds"
  )
  print(cluster_wilcoxon(d))
}



# ==============================================================================
# VIP-PAD FACTOR–COMORBIDITY MODELS | FIGURE S04
# ==============================================================================

# Function to estimate factor associations with comorbidities.
vip_comorbidity_models <- function() {
  # Age- and sex-adjusted logistic associations between factors and comorbidities.
  set.seed(42)
  trait_layout <- vip_trait_dictionary()
  wide_data <- vip_clinical_data()
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  required <- c(
    "HeartDisease",
    "CerebrovascularDisease",
    "DiabetesCategory",
    "CKDStage_Recorded",
    "Hypertension",
    "COPD",
    "Anemia",
    "Dyslipidemia",
    "Hyperhomocysteinemia",
    "Hyperuricemia"
  )
  wide_data <- wide_data[, c("study_id", "Age", "Sex", required), drop = FALSE]
  factor_metadata <- data.table(
    factor = sprintf('F%02d', seq_len(ncol(fit$Z))),
    score_variable = sprintf('F%02d', seq_len(ncol(fit$Z)))
  )
  wide_data[factor_metadata$score_variable] <- as.data.frame(
    fit$Z[as.character(wide_data$study_id), , drop = FALSE]
  )
  traits <- copy(trait_layout[group == 'Comorbidities'])
  disease_values <- list(
    coronary = source_binary(wide_data$HeartDisease, 'Coronary artery disease|Myocardial infarction|Angina'),
    stroke = source_binary(
      wide_data$CerebrovascularDisease,
      'stroke'
    ),
    diabetes = source_binary(wide_data$DiabetesCategory),
    ckd = source_binary(wide_data$CKDStage_Recorded),
    htn = source_binary(wide_data$Hypertension)
  )
  for (variable in c(
    'COPD',
    'Anemia',
    'Dyslipidemia',
    'Hyperhomocysteinemia',
    'Hyperuricemia'
  )) {
    disease_values[[tolower(variable)]] <- source_binary(wide_data[[variable]])
  }
  results <- list()
  for (factor_index in seq_len(ncol(fit$Z))) {
    factor_variable <- factor_metadata$score_variable[factor_index]
    factor_values <- wide_data[[factor_variable]]
    standardized_score <- as.numeric(scale(factor_values))
    for (trait_id in traits$id) {
      model_data <- data.frame(
        disease = disease_values[[trait_id]],
        score = standardized_score,
        Age = wide_data$Age,
        Sex = factor(wide_data$Sex, levels = c('Female', 'Male'))
      )
      model_data <- droplevels(model_data[complete.cases(model_data), ])
      design_matrix <- model.matrix(~ score + Age + Sex, data = model_data)
      stopifnot(qr(design_matrix)$rank == ncol(design_matrix))
      normalized_design <- sweep(
        design_matrix,
        2,
        apply(abs(design_matrix), 2, max),
        '/'
      )
      signed_design <- normalized_design * (2 * model_data$disease - 1)
      parameter_count <- ncol(signed_design)
      split_design <- cbind(signed_design, -signed_design)
      separation_fit <- lpSolve::lp(
        'max',
        colSums(split_design),
        rbind(split_design, rep(1, 2 * parameter_count)),
        c(rep('>=', nrow(signed_design)), '<='),
        c(rep(0, nrow(signed_design)), 1)
      )
      separated <- separation_fit$objval > 1e-7
      method <- 'MLE logistic; profile-likelihood CI'
      if (separated) {
        model_fit <- logistf::logistf(
          disease ~ score + Age + Sex,
          data = model_data,
          pl = TRUE
        )
        coefficient <- unname(coef(model_fit)['score'])
        interval <- unname(c(
          model_fit$ci.lower['score'],
          model_fit$ci.upper['score']
        ))
        p_value <- unname(model_fit$prob['score'])
        method <- 'Firth logistic; penalized profile-likelihood CI'
      } else {
        model_fit <- glm(
          disease ~ score + Age + Sex,
          data = model_data,
          family = binomial(),
          control = glm.control(maxit = 100)
        )
        stopifnot(model_fit$converged, all(is.finite(coef(model_fit))))
        coefficient <- unname(coef(model_fit)['score'])
        interval <- unname(suppressMessages(stats::confint(
          model_fit,
          parm = 'score'
        )))
        p_value <- coef(summary(model_fit))['score', 'Pr(>|z|)']
      }
      result_key <- paste(trait_id, factor_index)
      results[[result_key]] <- data.table(
        id = trait_id,
        factor = factor_metadata$factor[factor_index],
        N = nrow(model_data),
        cases = sum(model_data$disease),
        log_OR = coefficient,
        OR = exp(coefficient),
        CI_low = exp(interval[1]),
        CI_high = exp(interval[2]),
        P = p_value,
        method = method,
        adjustment = 'Age + Sex'
      )
    }
  }
  statistics <- merge(
    rbindlist(results),
    traits[, .(id, label)],
    by = 'id',
    sort = FALSE
  )

  dir.create(
    "./results/VIP-PAD/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    list(
      traits = traits,
      statistics = statistics
    ),
    "./results/VIP-PAD/clinical/Figure_3_Plot_Data.rds"
  )
}



# ==============================================================================
# FACTOR SCORES BY RUTHERFORD GRADE | FIGURE 2 (PANEL B)
# ==============================================================================

# Function to compare factor scores across Rutherford categories.
rutherford_associations <- function() {
  clinical <- as.data.table(vip_clinical_data())
  scores <- fread("./results/VIP-PAD/mofa/Factor_Scores_Raw.csv")
  joined <- merge(
    scores,
    clinical[, .(SampleID, RutherfordCategory)],
    by = 'SampleID',
    all.x = TRUE,
    sort = FALSE
  )
  joined[, rutherford := RutherfordCategory]
  joined <- joined[!is.na(rutherford)]
  joined[,
    severity := factor(
      ifelse(rutherford <= 3, '0-3', '4-6'),
      levels = c('0-3', '4-6')
    )
  ]
  cols <- setdiff(names(scores), 'SampleID')
  long <- melt(
    joined,
    id.vars = c('SampleID', 'severity'),
    measure.vars = cols,
    variable.name = 'factor',
    value.name = 'score'
  )
  long[,
    factor := factor(
      sprintf('F%02d', match(as.character(factor), cols)),
      levels = sprintf('F%02d', seq_along(cols))
    )
  ]
  stats <- long[,
    {
      a <- score[severity == '0-3']
      b <- score[severity == '4-6']
      test <- wilcox.test(b, a, exact = FALSE, correct = TRUE)
      list(
        n_0_3 = length(a),
        n_4_6 = length(b),
        median_0_3 = median(a),
        median_4_6 = median(b),
        P = test$p.value
      )
    },
    by = factor
  ]
  stats[,
    annotation := ifelse(
      P < .001,
      '***',
      ifelse(P < .01, '**', ifelse(P < .05, '*', ''))
    )
  ]
  dir.create(
    "./results/VIP-PAD/clinical",
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(
    stats,
    "./results/VIP-PAD/clinical/Rutherford_Factor_Statistics.csv"
  )
  fwrite(
    long[, .(severity, factor, score)],
    "./results/VIP-PAD/clinical/Rutherford_Factor_Scores.csv"
  )
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  # Baseline characteristics (Tables S10 and S11)
  baseline_table('VIP-PAD')
  baseline_table('UKB')

  # Factor associations with clinical traits
  vip_clinical()
  rutherford_associations()
  vip_comorbidity_models()
  ukb_factor_clinical_associations()

  # Clinical measures by molecular cluster
  vip_prepare_clinical_data()
  ukb_clinical_distributions()
}
