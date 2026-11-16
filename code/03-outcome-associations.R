# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, survival 3.8-6
#
# Global: Age- and sex-adjusted Cox models and event-free survival by molecular cluster.
#
# Inputs: VIP.rds and UKB.rds in ./Datasets_PAD_MOFA/; fitted MOFA models and cluster
#         assignments.
#
# Outputs: Cox and Kaplan–Meier plot data under ./results/. Tables S13 and S14 print to the console.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to convert binary values to analysis indicators.
ukb_binary <- function(v) {
  v <- as.character(v)
  ifelse(is.na(v), NA_integer_, as.integer(v %in% c('Yes', '1')))
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
  stopifnot(
    !anyNA(clinical$study_id),
    !anyDuplicated(clinical$study_id),
    isTRUE(all.equal(clinical$MALE_time_years, clinical$MACCE_time_years))
  )
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



# ==============================================================================
# UKB FACTOR–OUTCOME ASSOCIATIONS | FIGURE 3 (PANEL E)
# ==============================================================================

# Function to estimate factor-specific hazard ratios in UKB.
ukb_factor_outcome_cox <- function() {
  # Fit age- and sex-adjusted Cox models for each UKB factor and endpoint.
  set.seed(42)
  x <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
  d <- as.data.table(ukb_clinical_data())
  d <- d[match(rownames(x$Z), study_id)]
  z <- scale(x$Z)
  d[, age := as.numeric(Age_i0)]
  d[, sex := factor(Sex, levels = c('Female', 'Male'))]
  results <- flow <- checks <- list()
  for (ep in c('MALE', 'MACCE')) {
    e <- if (ep == 'MALE') 'MajorAdverseLimbEvent' else 'MACCE'
    base <- ukb_binary(d[[paste0('Baseline_', e, '_i0')]])
    event <- ukb_binary(d[[paste0('Incident_', e, '_i0')]])
    days <- as.numeric(d[[paste0('Time_', e, '_i0')]])
    eligible <- !is.na(base) &
      base == 0 &
      !is.na(event) &
      is.finite(days) &
      days > 0
    ok <- eligible & complete.cases(d$age, d$sex)
    dc <- data.frame(
      time = days[ok] / 365.25,
      event = event[ok],
      age = d$age[ok],
      sex = d$sex[ok]
    )
    flow[[ep]] <- data.table(
      endpoint = ep,
      MOFA_N = nrow(d),
      prior_endpoint = sum(base == 1, na.rm = TRUE),
      unknown_history = sum(is.na(base)),
      valid_endpoint_N = sum(eligible),
      covariate_exclusions = sum(eligible & !ok),
      N = nrow(dc),
      Events = sum(dc$event)
    )
    print(table(dc$sex, dc$event))
    if (any(table(dc$sex, dc$event) == 0)) {
      stop(
        'Sparse zero sex-event cell: assess separation and estimator before plotting.'
      )
    }
    for (k in seq_len(ncol(z))) {
      dc$score <- z[ok, k]
      fit <- coxph(
        Surv(time, event) ~ score + age + sex,
        data = dc,
        ties = 'efron',
        x = TRUE,
        y = TRUE
      )
      s <- summary(fit)
      key <- paste(ep, k)
      results[[key]] <- data.table(
        endpoint = ep,
        factor = k,
        N = nrow(dc),
        Events = sum(dc$event),
        HR = s$conf.int['score', 'exp(coef)'],
        CI_low = s$conf.int['score', 'lower .95'],
        CI_high = s$conf.int['score', 'upper .95'],
        P = s$coefficients['score', 'Pr(>|z|)']
      )
      ph <- cox.zph(fit)
      checks[[key]] <- data.table(
        endpoint = ep,
        factor = k,
        PH_factor_P = ph$table['score', 'p'],
        PH_global_P = ph$table['GLOBAL', 'p']
      )
    }
  }
  r <- rbindlist(results)
  r[, FDR := p.adjust(P, 'BH'), by = endpoint]
  fl <- rbindlist(flow)
  ph <- rbindlist(checks)
  ph[, PH_factor_FDR := p.adjust(PH_factor_P, 'BH'), by = endpoint]
  dir.create(
    "./results/UKB/outcomes/Cox",
    recursive = TRUE,
    showWarnings = FALSE
  )

  fwrite(r, "./results/UKB/outcomes/Cox/Associations.csv")
  print(fl)
  print(ph)

}



# ==============================================================================
# UKB MACCE-FREE SURVIVAL BY CLUSTER | FIGURE 3 (PANEL G)
# ==============================================================================

# Function to estimate MACCE-free survival by UKB cluster.
ukb_macce_cluster_km <- function() {
  # Estimate MACCE-free survival by UKB molecular cluster.
  d <- as.data.table(ukb_clinical_data())
  assignment <- fread(
    "./results/UKB/clustering/Patient_Assignments.csv"
  )
  d[,
    Group := factor(
      paste0('C', assignment$K2[match(study_id, assignment$sample_id)]),
      levels = c('C1', 'C2')
    )
  ]
  base <- ukb_binary(d[['Baseline_MACCE_i0']])
  event <- ukb_binary(d[['Incident_MACCE_i0']])
  days <- as.numeric(d[['Time_MACCE_i0']])
  ok <- !is.na(base) &
    base == 0 &
    !is.na(event) &
    is.finite(days) &
    days > 0 &
    !is.na(d$Group)
  z <- data.table(
    Group = d$Group[ok],
    time = days[ok] / (365.25 / 12),
    event = event[ok]
  )
  fit <- survfit(
    Surv(time, event) ~ Group,
    data = z,
    conf.type = 'log-log',
    conf.int = .95
  )
  sk <- summary(fit, censored = TRUE)
  sr <- summary(fit, times = c(0, 12, 24, 36), extend = FALSE)
  curves <- data.table(
    endpoint = "MACCE",
    Group = sub('Group=', '', as.character(sk$strata)),
    time = sk$time,
    S = sk$surv,
    lower = sk$lower,
    upper = sk$upper
  )
  risk <- data.table(
    endpoint = "MACCE",
    Group = sub('Group=', '', as.character(sr$strata)),
    time = sr$time,
    n_risk = sr$n.risk
  )
  counts <- z[,
    .(
      N = .N,
      Events = sum(event),
      Events_36months = sum(event == 1 & time <= 36),
      Max_months = max(time)
    ),
    by = Group
  ][, endpoint := "MACCE"]
  lr <- survdiff(Surv(time, event) ~ Group, data = z, rho = 0)
  tests <- data.table(
    endpoint = "MACCE",
    Logrank_P = pchisq(lr$chisq, df = 1, lower.tail = FALSE),
    chi_square = lr$chisq
  )
  flow <- data.table(
    endpoint = "MACCE",
    MOFA_N = nrow(d),
    prior_endpoint = sum(base == 1, na.rm = TRUE),
    unknown_history = sum(is.na(base)),
    N = nrow(z),
    Events = sum(z$event)
  )
  o <- list(curves = curves, risk = risk, counts = counts, tests = tests)
  dir.create("./results/UKB/outcomes/KM_MACCE", recursive = TRUE, showWarnings = FALSE)
  for (n in names(o)) {
    fwrite(o[[n]], file.path("./results/UKB/outcomes/KM_MACCE", paste0(n, '.csv')))
  }
  print(flow)
  saveRDS(list(MACCE = fit), "./results/UKB/outcomes/KM_MACCE/KM_Fits.rds")

}



# ==============================================================================
# VIP-PAD FACTOR–OUTCOME ASSOCIATIONS | FIGURE 2 (PANEL C) AND TABLES S13 AND S14
# ==============================================================================

# Function to tabulate VIP-PAD factor associations with outcomes.
vip_factor_outcomes <- function() {
  # VIP-PAD: full-follow-up factor associations, adjusted for age and sex only.
  # Estimate associations per standard deviation and across score tertiles.
  d <- as.data.table(vip_clinical_data())
  f <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  f$Z <- f$Z[as.character(d$SampleID), , drop = FALSE]
  d[, Sex := factor(Sex, levels = c('Female', 'Male'))]

  eligible <- complete.cases(d$Age, d$Sex) &
    is.finite(d$FollowUp_Months) &
    d$FollowUp_Months > 0
  res <- list()
  flows <- diagnostics <- list()
  for (ep in c('MALE', 'MACCE')) {
    event <- d[[paste0(ep, '_binary')]]
    included <- eligible & !is.na(event)
    flows[[ep]] <- data.table(
      endpoint = ep, N_initial = nrow(d), N = sum(included),
      Events = sum(event[included]), known_status = sum(!is.na(event)),
      base_eligible = sum(eligible)
    )
    for (k in seq_len(ncol(f$Z))) {
      score <- as.numeric(scale(f$Z[, k]))
      cp <- quantile(f$Z[, k], c(1 / 3, 2 / 3), type = 7)
      group <- cut(f$Z[, k], c(-Inf, cp, Inf), labels = c('T1', 'T2', 'T3'))
      q <- data.table(
        months = d$FollowUp_Months,
        Age = d$Age,
        Sex = d$Sex,
        time = d$FollowUp_Months / 12,
        event = event,
        score = score,
        group = group
      )[included]
      mc <- coxph(
        Surv(months, event) ~ score + Age + Sex,
        data = q, ties = 'efron', x = TRUE, y = TRUE
      )
      mt <- coxph(
        Surv(time, event) ~ group + Age + Sex,
        data = q, ties = 'efron', x = TRUE, y = TRUE
      )
      for (row in c('Per 1 SD', 'T1', 'T2', 'T3')) {
        sub <- if (row == 'Per 1 SD') q else q[group == row]
        ne <- sum(sub$event)
        py <- sum(sub$time)
        ci <- poisson.test(ne, T = py)$conf.int * 1000
        h <- lo <- hi <- pv <- NA_real_
        if (row != 'T1') {
          m <- if (row == 'Per 1 SD') mc else mt
          nm <- if (row == 'Per 1 SD') 'score' else paste0('group', row)
          beta <- coef(m)[nm]
          se <- sqrt(vcov(m)[nm, nm])
          h <- exp(beta)
          lo <- exp(beta - 1.96 * se)
          hi <- exp(beta + 1.96 * se)
          pv <- summary(m)$coefficients[nm, 'Pr(>|z|)']
        }
        res[[length(res) + 1]] <- data.table(
          endpoint = ep,
          factor = k,
          category = row,
          N = nrow(sub),
          events = ne,
          person_years = py,
          rate = ne / py * 1000,
          rate_low = ci[1],
          rate_high = ci[2],
          HR = h,
          HR_low = lo,
          HR_high = hi,
          P = pv
        )
      }
      for (nm in c('continuous', 'tertiles')) {
        ph <- cox.zph(if (nm == 'continuous') mc else mt)
        diagnostics[[length(diagnostics) + 1]] <- data.table(
          endpoint = ep,
          factor = k,
          model = nm,
          PH_factor_P = ph$table[1, 'p'],
          PH_global_P = ph$table['GLOBAL', 'p']
        )
      }
    }
  }
  r <- rbindlist(res)
  associations <- r[category == 'Per 1 SD', .(
    endpoint, factor = sprintf('F%02d', factor), N, Events = events,
    HR, CI_low = HR_low, CI_high = HR_high, P
  )]
  associations[, FDR := p.adjust(P, 'BH'), by = endpoint]
  dir.create("./results/VIP-PAD/outcomes/Cox", recursive = TRUE, showWarnings = FALSE)
  fwrite(associations, "./results/VIP-PAD/outcomes/Cox/Associations.csv")
  print(rbindlist(flows))
  print(rbindlist(diagnostics))
  headers <- c(
    'Factor / category',
    'Total',
    'Events',
    'Person-years',
    'Incidence rate (95% CI)\nper 1,000 person-years',
    'Adjusted HR (95% CI)',
    'P-value'
  )
  for (ep in c('MALE', 'MACCE')) {
    rows <- list()
    for (k in seq_len(ncol(f$Z))) {
      rows[[length(rows) + 1]] <- c(paste('Factor', k), rep('', 6))
      for (i in which(r$endpoint == ep & r$factor == k)) {
        x <- r[i]
        rows[[length(rows) + 1]] <- c(
          if (x$category == 'Per 1 SD') 'Continuous, per 1 SD' else x$category,
          as.character(x$N),
          as.character(x$events),
          sprintf('%.2f', x$person_years),
          sprintf('%.2f (%.2f–%.2f)', x$rate, x$rate_low, x$rate_high),
          if (x$category == 'T1') {
            'Reference'
          } else {
            sprintf('%.2f (%.2f–%.2f)', x$HR, x$HR_low, x$HR_high)
          },
          if (x$category == 'T1') {
            '—'
          } else if (x$P < .001) {
            '<0.001'
          } else {
            sprintf('%.3f', x$P)
          }
        )
      }
    }
    mat <- do.call(rbind, rows)
    tab <- as.data.table(mat)
    setnames(tab, headers)
    cat(if (ep == 'MALE') '\nS13. MALE\n' else '\nS14. MACCE\n')
    print(tab, nrows = Inf, trunc.cols = FALSE)
  }

}



# ==============================================================================
# VIP-PAD CLUSTER SURVIVAL | FIGURE 2 (PANELS H AND I)
# ==============================================================================

# Function to estimate VIP-PAD event-free survival by molecular cluster.
vip_cluster_survival <- function() {
  set.seed(42)
  # Compare event-free survival between the two clusters using the log-rank test.
  dir.create(
    "./results/VIP-PAD/outcomes/KM",
    recursive = TRUE,
    showWarnings = FALSE
  )
  d <- as.data.table(vip_clinical_data())
  assignment <- fread(
    "./results/VIP-PAD/clustering/Patient_Assignments.csv"
  )
  stopifnot(
    setequal(assignment$sample_id, d$SampleID),
    !anyDuplicated(assignment$sample_id)
  )
  d[,
    Group := factor(
      paste0('C', assignment$K2[match(SampleID, assignment$sample_id)]),
      levels = c('C1', 'C2')
    )
  ]

  d[, valid_time := is.finite(FollowUp_Months) & FollowUp_Months > 0]
  counts <- curves <- risk <- tests <- follow <- list()
  for (e in c('MALE', 'MACCE')) {
    z <- copy(d[valid_time & !is.na(get(paste0(e, '_binary')))])
    z[, event := get(paste0(e, '_binary'))]
    counts[[e]] <- z[, .(N = .N, Events = sum(event)), by = Group][,
      endpoint := e
    ]
    km <- survfit(
      Surv(FollowUp_Months, event) ~ Group,
      z,
      conf.type = 'log-log'
    )
    sk <- summary(km, censored = TRUE)
    curves[[e]] <- data.table(
      endpoint = e,
      Group = sub('Group=', '', as.character(sk$strata)),
      time = sk$time,
      S = sk$surv,
      lower = sk$lower,
      upper = sk$upper,
      n_censor = sk$n.censor
    )
    sr <- summary(km, times = c(0, 12, 24, 36), extend = TRUE)
    risk[[e]] <- data.table(
      endpoint = e,
      Group = sub('Group=', '', as.character(sr$strata)),
      time = sr$time,
      n_risk = sr$n.risk
    )
    lr <- survdiff(Surv(FollowUp_Months, event) ~ Group, z)
    tests[[e]] <- data.table(
      endpoint = e,
      Logrank_P = pchisq(lr$chisq, 1, lower.tail = FALSE)
    )
    follow[[e]] <- z[,
      .(N = .N, median_followup_months = median(FollowUp_Months)),
      by = Group
    ][, endpoint := e]
  }
  o <- list(
    counts = rbindlist(counts),
    curves = rbindlist(curves),
    risk = rbindlist(risk),
    tests = rbindlist(tests),
    followup = rbindlist(follow)
  )

  fwrite(
    o$counts,
    "./results/VIP-PAD/outcomes/KM/counts.csv"
  )
  fwrite(
    o$curves,
    "./results/VIP-PAD/outcomes/KM/curves.csv"
  )
  fwrite(
    o$risk,
    "./results/VIP-PAD/outcomes/KM/risk.csv"
  )
  fwrite(o$tests, "./results/VIP-PAD/outcomes/KM/tests.csv")
  print(o$followup)
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  vip_factor_outcomes()
  vip_cluster_survival()
  ukb_factor_outcome_cox()
  ukb_macce_cluster_km()
}
