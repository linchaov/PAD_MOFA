# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: data.table 1.18.6.1, survival 3.8-6, xgboost 3.2.1.1
#
# Global: Factor tertiles, event-free survival and molecular contributions to predicted
#         event-free time.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds, ./Datasets_PAD_MOFA/UKB.rds, fitted MOFA models
#         from 01-mofa.R. VIP-PAD metabolite labels are read from ../inputs/metabolites.rds.
#
# Outputs: Kaplan-Meier estimates, AFT performance estimates, SHAP rankings and holdout contributions
#          in ./results/VIP-PAD/ and ./results/UKB/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
  library(survival)
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

# Function to encode event times and censoring bounds for XGBoost.
make_dm <- function(X, time, event) {
  z <- xgb.DMatrix(X)
  setinfo(z, 'label_lower_bound', time)
  setinfo(z, 'label_upper_bound', ifelse(event == 1, time, Inf))
  z
}
# Function to assign event-stratified cross-validation folds.
fold_assign <- function(event, seed) {
  set.seed(seed)
  f <- integer(length(event))
  for (v in 0:1) {
    j <- which(event == v)
    f[j] <- sample(rep(1:5, length.out = length(j)))
  }
  f
}
# Function to estimate boosting rounds and molecular importance in training folds.
cv_train <- function(
  X,
  time,
  event,
  folds,
  seed,
  prefix,
  shap,
  params
) {
  records <- list()
  oof <- matrix(NA_real_, nrow(X), ncol(X), dimnames = list(NULL, colnames(X)))
  op <- rep(NA_real_, nrow(X))
  for (k in 1:5) {
    it <- which(folds != k)
    iv <- which(folds == k)
    par <- params
    par$seed <- seed + k
    dt <- make_dm(X[it, , drop = FALSE], time[it], event[it])
    dv <- make_dm(X[iv, , drop = FALSE], time[iv], event[iv])
    m <- xgb.train(
      par,
      dt,
      nrounds = 1500,
      evals = list(valid = dv),
      early_stopping_rounds = 50,
      verbose = 0
    )
    el <- as.data.table(attr(m, 'evaluation_log'))
    metric <- grep('aft', names(el), value = TRUE)[1]
    best <- which.min(el[[metric]])
    # Refit using the round with the lowest validation loss.
    m <- xgb.train(par, dt, nrounds = best, verbose = 0)
    op[iv] <- predict(m, dv, outputmargin = TRUE)
    if (shap) {
      sc <- predict(m, dv, predcontrib = TRUE, approxcontrib = FALSE)
      oof[iv, ] <- sc[, seq_len(ncol(X)), drop = FALSE]
    }
    records[[k]] <- data.table(
      Fold = k,
      Train_N = length(it),
      Validation_N = length(iv),
      Train_events = sum(event[it]),
      Validation_events = sum(event[iv]),
      Best_rounds = best,
      Best_validation_nloglik = min(el[[metric]]),
      Reached_max_rounds = nrow(el) == 1500
    )
    cat(prefix, 'fold', k, 'best rounds', best, '\n')
    flush.console()
  }
  records <- rbindlist(records)
  nr <- max(1L, as.integer(round(median(records$Best_rounds))))
  list(
    rounds = nr,
    records = records,
    oof_SHAP = if (shap) oof else NULL,
    oof_prediction = op
  )
}

# Fit endpoint-specific AFT models within each cohort.



# ==============================================================================
# SURVIVAL MODELLING AND SHAP ATTRIBUTION
# ==============================================================================

# Function to fit and evaluate the cohort-specific accelerated failure time model.
fit_aft <- function(cohort, factor_number, endpoint, a) {
  model_id <- paste("Factor", factor_number, endpoint, sep = "_")
  invisible(Sys.setlocale('LC_CTYPE', 'English_United States.utf8'))
  vip <- cohort == "VIP-PAD"
  score_name <- if (vip) paste0("Factor_", factor_number) else "Factor_score"
  set.seed(42)
  input_dir <- file.path("./results", cohort, "outcomes", model_id)
  params <- list(
    objective = 'survival:aft',
    eval_metric = 'aft-nloglik',
    aft_loss_distribution = 'normal',
    aft_loss_distribution_scale = 1.2,
    tree_method = 'hist',
    eta = .03,
    max_depth = 2,
    min_child_weight = 10,
    lambda = 5,
    alpha = 0,
    subsample = .8,
    colsample_bytree = .8,
    nthread = 4
  )
  # Function to summarize survival prediction performance for the fitted models.
  metric_table <- function(q, pred, tag, endpoint) {
    d <- data.frame(time = q$time, event = q$event, pred = pred)
    hc <- concordance(
      Surv(time, event) ~ pred,
      data = d,
      reverse = FALSE,
      timewt = 'n'
    )
    uc <- concordance(
      Surv(time, event) ~ pred,
      data = d,
      reverse = FALSE,
      timewt = 'n/G2',
      ymax = if (vip) 36 else 18
    )
    metrics <- data.table(
      endpoint = endpoint,
      model = tag,
      holdout_N = nrow(d),
      holdout_events = sum(d$event),
      Harrell_C = hc$concordance,
      Harrell_SE = sqrt(hc$var),
      Harrell_lower = max(0, hc$concordance - 1.96 * sqrt(hc$var)),
      Harrell_upper = min(1, hc$concordance + 1.96 * sqrt(hc$var)),
      Uno_C = uc$concordance,
      Uno_SE = sqrt(uc$var)
    )
    setnames(metrics, "Uno_C", if (vip) "Uno_C_36months" else "Uno_C_18years")
    metrics
  }
  for (e in endpoint) {
    seed <- if (e == 'MALE') 42L else 43L
    set.seed(seed)
    q <- copy(a$endpoint_data[[e]])
    X <- a$X[match(q$sample_id, a$sample_id), , drop = FALSE]
    it <- sort(unlist(lapply(0:1, function(v) {
      j <- which(q$event == v)
      sample(j, floor(.7 * length(j)))
    })))
    iv <- setdiff(seq_len(nrow(q)), it)
    folds <- fold_assign(q$event[it], seed + 100)
    split <- q[, .(sample_id, time, event)]
    split[, subset := ifelse(seq_len(.N) %in% it, 'Training', 'Holdout')]
    cat(
      e,
      'training',
      length(it),
      'events',
      sum(q$event[it]),
      '; holdout',
      length(iv),
      'events',
      sum(q$event[iv]),
      '\n'
    )
    flush.console()
    cvfull <- cv_train(
      X[it, , drop = FALSE],
      q$time[it],
      q$event[it],
      folds,
      seed,
      paste0(e, '_full'),
      TRUE,
      params
    )
    imp <- data.table(
      feature = colnames(X),
      Mean_abs_SHAP = colMeans(abs(cvfull$oof_SHAP))
    )
    imp <- merge(
      imp,
      a$features[, .(feature, label, modality, weight, within_view_rank)],
      by = 'feature',
      sort = FALSE
    )
    imp <- imp[order(-Mean_abs_SHAP, feature)]
    imp[, rank := seq_len(.N)]
    active <- imp[Mean_abs_SHAP > 1e-12]
    active[, rank_normalized := (rank - 1) / (.N - 1)]
    active[,
      importance_normalized := (Mean_abs_SHAP - min(Mean_abs_SHAP)) /
        (max(Mean_abs_SHAP) - min(Mean_abs_SHAP))
    ]
    active[,
      elbow_distance := (1 - rank_normalized - importance_normalized) / sqrt(2)
    ]
    elbow <- active$rank[which.max(active$elbow_distance)]
    panel <- imp[rank <= elbow, feature]
    cvpanel <- cv_train(
      X[it, panel, drop = FALSE],
      q$time[it],
      q$event[it],
      folds,
      seed + 200,
      paste0(e, '_panel'),
      FALSE,
      params
    )
    par <- params
    par$seed <- seed
    full <- xgb.train(
      par,
      make_dm(X[it, , drop = FALSE], q$time[it], q$event[it]),
      nrounds = cvfull$rounds,
      verbose = 0
    )
    model <- xgb.train(
      par,
      make_dm(X[it, panel, drop = FALSE], q$time[it], q$event[it]),
      nrounds = cvpanel$rounds,
      verbose = 0
    )
    # Evaluate the fitted models in held-out patients.
    dvfull <- make_dm(X[iv, , drop = FALSE], q$time[iv], q$event[iv])
    dvpanel <- make_dm(X[iv, panel, drop = FALSE], q$time[iv], q$event[iv])
    pf <- predict(full, dvfull, outputmargin = TRUE)
    pp <- predict(model, dvpanel, outputmargin = TRUE)
    sc <- predict(model, dvpanel, predcontrib = TRUE, approxcontrib = FALSE)
    residual <- max(abs(rowSums(sc) - pp))
    colnames(sc) <- c(panel, 'BIAS')
    # Training-fitted factor-only AFT model provides a compact descriptive comparator.
    sf <- survreg(
      reformulate(score_name, response = "Surv(time,event)"),
      data = as.data.frame(q[it]),
      dist = 'lognormal'
    )
    p1 <- predict(sf, newdata = as.data.frame(q[iv]), type = 'lp')
    met <- rbindlist(list(
      metric_table(q[iv], pf, paste0('Full_', ncol(X), '_molecules'), e),
      metric_table(q[iv], pp, paste0('Elbow_panel_', elbow), e),
      metric_table(
        q[iv],
        p1,
        paste0('F', factor_number, '_only_lognormal_AFT'),
        e
      )
    ))
    result <- list(
      endpoint = e,
      seed = seed,
      params = params,
      train_index = it,
      holdout_index = iv,
      folds = folds,
      importance = imp,
      elbow_curve = active,
      elbow = elbow,
      panel = panel,
      cv_full = cvfull,
      cv_panel = cvpanel,
      holdout_SHAP = sc,
      holdout_X = X[iv, panel, drop = FALSE],
      holdout_sample_id = q$sample_id[iv],
      metrics = met,
      additivity_max_error = residual,
      selected_feature_labels = a$features,
      split_summary = split[, .(N = .N, Events = sum(event)), by = subset]
    )

    saveRDS(
      result,
      file.path(input_dir, paste0(e, '_shap_results.rds'))
    )
    cat(e, 'elbow=', elbow, '\n')
    flush.console()
  }

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
# UKB FACTOR SURVIVAL AND ATTRIBUTION | FIGURES S10 AND S23
# ==============================================================================

# Function to prepare survival and attribution inputs for UKB Factors 2 and 5.
ukb_prepare_f2_f5_outcome_characterization <- function() {
  f <- readRDS("./results/UKB/mofa/Selected_Fit.rds")
  cohort <- readRDS("./Datasets_PAD_MOFA/UKB.rds")
  ids <- cohort$study_id
  stopifnot(!anyNA(ids), !anyDuplicated(ids), setequal(ids, rownames(f$Z)))
  meta <- as.data.table(attr(cohort, "feature_dictionary"))
  meta[,
    display_label := ifelse(
      modality == "Protein",
      toupper(label),
      pretty_nmr(label)
    )
  ]
  # Function to extract measurements for one omics view.
  raw_view <- function(view) {
    dictionary <- meta[modality == view]
    x <- as.matrix(cohort[, dictionary$source_field, drop = FALSE])
    dimnames(x) <- list(ids, dictionary$feature)
    x
  }
  raw <- list(
    Protein = raw_view("Protein"),
    Metabolite = raw_view("Metabolite"),
    sample_id = ids
  )
  d <- cohort
  # Function to convert recorded binary values to analysis indicators.
  bin <- function(x) {
    z <- as.character(x)
    ifelse(is.na(z), NA_integer_, as.integer(z %in% c('Yes', '1')))
  }
  for (k in c(2, 5)) {
    for (e in if (k == 2) 'MACCE' else c('MALE', 'MACCE')) {
      key <- paste0('Factor_', k, '_', e)
      score <- f$Z[ids, k]
      cp <- quantile(score, c(1 / 3, 2 / 3), type = 7)
      tt <- cut(score, c(-Inf, cp, Inf), labels = c('Low', 'Medium', 'High'))
      nm <- if (e == 'MALE') 'MajorAdverseLimbEvent' else 'MACCE'
      base <- bin(d[[paste0('Baseline_', nm, '_i0')]])
      event <- bin(d[[paste0('Incident_', nm, '_i0')]])
      days <- as.numeric(d[[paste0('Time_', nm, '_i0')]])
      ok <- !is.na(base) &
        base == 0 &
        !is.na(event) &
        is.finite(days) &
        days > 0
      q <- data.table(
        sample_id = ids[ok],
        Factor_score = score[ok],
        Tertile = tt[ok],
        time = days[ok] / 365.25,
        event = event[ok]
      )
      km <- survfit(Surv(time, event) ~ Tertile, q, conf.type = 'log-log')
      sk <- summary(km, censored = TRUE)
      curves <- data.table(
        endpoint = e,
        Tertile = sub('Tertile=', '', as.character(sk$strata)),
        time = sk$time,
        S = sk$surv,
        lower = sk$lower,
        upper = sk$upper,
        n_risk = sk$n.risk,
        n_event = sk$n.event,
        n_censor = sk$n.censor
      )
      sr <- summary(km, times = c(0, 6, 12, 18), extend = TRUE)
      risk <- data.table(
        endpoint = e,
        Tertile = sub('Tertile=', '', as.character(sr$strata)),
        time = sr$time,
        n_risk = sr$n.risk
      )
      counts <- q[,
        .(
          N = .N,
          Events = sum(event),
          Events_by_18 = sum(event == 1 & time <= 18)
        ),
        by = Tertile
      ][, endpoint := e]
      lr <- survdiff(Surv(time, event) ~ Tertile, q)
      tests <- data.table(
        endpoint = e,
        P = pchisq(lr$chisq, 2, lower.tail = FALSE)
      )
      endpoints <- setNames(list(q), e)
      dir.create(
        file.path("./results/UKB/outcomes", key),
        recursive = TRUE,
        showWarnings = FALSE
      )

      saveRDS(
        list(
          cutpoints = cp,
          counts = counts,
          tests = tests,
          curves = curves,
          risk = risk,
          endpoint_data = endpoints
        ),
        file.path("./results/UKB/outcomes", key, "km_results.rds")
      )
      # Collect the molecular weights for the target factor across views.
      weights <- rbindlist(lapply(names(f$W), function(v) {
        w <- f$W[[v]][, k]
        r <- data.table(
          feature = names(w),
          weight = as.numeric(w),
          modality = v
        )
        r[, label := meta$display_label[match(feature, meta$feature)]]
        r[v == 'Protein', label := toupper(label)]
        r[, abs_weight := abs(weight)]
        setorder(r, -abs_weight, feature)
        r[, within_view_rank := seq_len(.N)]
        r[, selected := within_view_rank <= ceiling(.05 * .N)]
        r
      }))
      features <- weights[selected == TRUE]
      X <- cbind(raw$Protein, raw$Metabolite)[
        match(ids, raw$sample_id),
        features$feature,
        drop = FALSE
      ]
      rownames(X) <- ids
      fit_aft("UKB", k, e, list(
        X = X, features = features, sample_id = ids, endpoint_data = endpoints
      ))
    }
  }
}

# Function to look up display labels for molecular features.
vip_feature_labels <- function(features, view, metabolites) {
  labels <- if (view == "Protein") {
    vip_symbols(features)
  } else {
    metabolites$label[match(features, metabolites$feature)]
  }
  missing <- is.na(labels) | labels == ""
  labels[missing] <- features[missing]
  labels
}

# VIP-PAD stores standardized molecular measurements. Restore the vendor log2 scale
# using the feature-specific centers and standard deviations supplied with the cohort.
# Function to prepare VIP-PAD abundance and event-time data.
vip_aft_data <- function(cohort, outcome) {
  ids <- cohort$study_id
  rownames(cohort) <- ids
  scale_attributes <- attr(cohort, "omics_matrix_attributes")
  # Function to assemble feature abundance values for the model input.
  abundance <- function(view) {
    pars <- scale_attributes[[view]]
    features <- names(pars[["scaled:scale"]])
    x <- as.matrix(cohort[, features, drop = FALSE])
    sweep(
      sweep(x, 2, pars[["scaled:scale"]][features], "*"),
      2,
      pars[["scaled:center"]][features],
      "+"
    )
  }
  clinical <- data.table::data.table(
    SampleID = ids,
    FollowUp_Months = 12 * cohort[[paste0(outcome, "_time_years")]],
    MALE_binary = cohort$MALE_event,
    MACCE_binary = cohort$MACCE_event
  )
  list(
    sample_id = ids,
    clinical = clinical,
    proteomics = abundance("proteomics"),
    metabolomics = abundance("metabolomics")
  )
}



# ==============================================================================
# VIP-PAD FACTOR SURVIVAL AND ATTRIBUTION
# ==============================================================================

# Function to prepare tertile survival and molecular attribution inputs for one VIP-PAD factor.
vip_factor_outcome_inputs <- function(k, outcome) {
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  a <- vip_aft_data(readRDS("./Datasets_PAD_MOFA/VIP.rds"), outcome)
  metabolites <- readRDS("../inputs/metabolites.rds")$metabolite
  factor_name <- paste0("Factor", k)
  score_column <- paste0("Factor_", k)
  event_column <- paste0(outcome, "_binary")
  # Define score tertiles across VIP-PAD before restricting to the endpoint cohort.
  ids <- a$sample_id
  score <- fit$Z[as.character(ids), factor_name]
  d <- as.data.table(a$clinical)[match(ids, SampleID)]
  d[, sample_id := ids]
  d[, (score_column) := score]
  cutpoints <- quantile(score, c(1 / 3, 2 / 3), type = 7)
  d[,
    Tertile := cut(
      score,
      c(-Inf, cutpoints, Inf),
      labels = c('Low', 'Medium', 'High')
    )
  ]

  q <- d[
    is.finite(FollowUp_Months) &
      FollowUp_Months > 0 &
      !is.na(get(event_column)),
    .(
      sample_id,
      factor_score = get(score_column),
      Tertile,
      time = FollowUp_Months,
      event = get(event_column)
    )
  ]
  setnames(q, "factor_score", score_column)
  km <- survfit(Surv(time, event) ~ Tertile, q, conf.type = 'log-log')
  sk <- summary(km, censored = TRUE)
  curves <- data.table(
    endpoint = outcome,
    Tertile = sub('Tertile=', '', as.character(sk$strata)),
    time = sk$time,
    S = sk$surv,
    lower = sk$lower,
    upper = sk$upper,
    n_risk = sk$n.risk,
    n_event = sk$n.event,
    n_censor = sk$n.censor
  )
  sr <- summary(km, times = c(0, 12, 24, 36), extend = TRUE)
  risk <- data.table(
    endpoint = outcome,
    Tertile = sub('Tertile=', '', as.character(sr$strata)),
    time = sr$time,
    n_risk = sr$n.risk
  )
  counts <- q[,
    .(
      N = .N,
      Events = sum(event),
      Events_by_36 = sum(event == 1 & time <= 36),
      Median_observed_time = median(time)
    ),
    by = Tertile
  ][, endpoint := outcome]
  lr <- survdiff(Surv(time, event) ~ Tertile, q)
  tests <- data.table(
    endpoint = outcome,
    P = pchisq(lr$chisq, 2, lower.tail = FALSE)
  )
  endpoint_data <- setNames(list(q), outcome)
  dir.create(
    sprintf("./results/VIP-PAD/outcomes/Factor_%d_%s", k, outcome),
    recursive = TRUE,
    showWarnings = FALSE
  )

  saveRDS(
    list(
      cutpoints = cutpoints,
      counts = counts,
      tests = tests,
      curves = curves,
      risk = risk,
      endpoint_data = endpoint_data
    ),
    sprintf(
      "./results/VIP-PAD/outcomes/Factor_%d_%s/km_results.rds",
      k,
      outcome
    )
  )
  # Rank features by absolute loading within each omics layer and retain the top 5%.
  weights <- rbindlist(lapply(c('Protein', 'Metabolite'), function(v) {
    w <- fit$W[[v]][, factor_name]
    r <- data.table(
      feature = names(w),
      weight = as.numeric(w),
      abs_weight = abs(w),
      modality = v
    )
    r[, label := vip_feature_labels(feature, v, metabolites)]
    setorder(r, -abs_weight, feature)
    r[, within_view_rank := seq_len(.N)]
    r[, selected := within_view_rank <= ceiling(.05 * .N)]
    r
  }))
  features <- weights[selected == TRUE]
  X <- cbind(a$proteomics, a$metabolomics)[, features$feature, drop = FALSE]
  stopifnot(all(is.finite(X)))
  fit_aft("VIP-PAD", k, outcome, list(
    X = X, features = features, sample_id = ids, endpoint_data = endpoint_data
  ))
}



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

if (sys.nframe() == 0L || interactive()) {
  ukb_prepare_f2_f5_outcome_characterization()
  # Factor 8 was associated with both outcomes.
  analyses <- data.frame(
    factor = c(1L, 2L, 3L, 4L, 8L, 8L),
    outcome = c("MALE", "MACCE", "MACCE", "MACCE", "MACCE", "MALE")
  )
  for (i in seq_len(nrow(analyses))) {
    vip_factor_outcome_inputs(analyses$factor[i], analyses$outcome[i])
  }
}
