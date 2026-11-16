# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: coxphf 1.13.4, data.table 1.18.6.1, glmnet 5.0,
#                  mice 3.19.0, pec 2025.06.24, survival 3.8-6
#
# Global: Nested cross-validation of 12-month MALE and MACCE risk in VIP-PAD.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds and VIP-PAD factor scores and cluster assignments
#         from 01-mofa.R. Candidate names and coding for S23 and S24 are listed below.
#
# Outputs: Fitted models, held-out predictions, performance measures and reclassification
#          under ./results/prediction/; supplementary tables S23 and S24 are printed to the console.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(mice)
  library(survival)
  library(glmnet)
})



# ==============================================================================
# FUNCTIONS
# ==============================================================================

# Function to align clinical predictors, factor scores and cluster assignments by study ID.
prediction_inputs <- function(cfg) {
  cohort <- readRDS("./Datasets_PAD_MOFA/VIP.rds")
  fit <- readRDS(
    "./results/VIP-PAD/mofa/Selected_Fit.rds"
  )
  assignments <- fread(
    "./results/VIP-PAD/clustering/Patient_Assignments.csv"
  )
  ids <- cohort$study_id
  stopifnot(
    !anyNA(ids),
    !anyDuplicated(ids),
    !anyDuplicated(assignments$sample_id),
    setequal(ids, rownames(fit$Z)),
    setequal(ids, assignments$sample_id)
  )
  groups <- assignments$K2[match(ids, assignments$sample_id)]
  cols <- grep("^prediction_", names(cohort), value = TRUE)
  # Function to assemble eligible participants, event times and indicators for one endpoint.
  build_endpoint <- function(endpoint) {
    time <- cohort[[paste0(endpoint, "_time_years")]]
    event <- cohort[[paste0(endpoint, "_event")]]
    ok <- is.finite(time) & time > 0 & event %in% 0:1
    z <- cohort[ok, cols, drop = FALSE]
    names(z) <- sub("^prediction_", "", cols)
    for (v in names(z)) {
      if (is.character(z[[v]])) z[[v]] <- factor(z[[v]])
    }
    factors <- fit$Z[
      ids[ok],
      sub("^F", "Factor", cfg$factors[[endpoint]]),
      drop = FALSE
    ]
    colnames(factors) <- cfg$factors[[endpoint]]
    cluster <- matrix(
      as.numeric(groups[ok] == 2L),
      ncol = 1L,
      dimnames = list(NULL, "clusterC2")
    )
    list(
      z = z,
      sample_id = ids[ok],
      event = as.integer(event[ok]),
      time = time[ok],
      factors = factors,
      cluster = cluster,
      source = "Datasets_PAD_MOFA/VIP.rds"
    )
  }
  setNames(lapply(names(cfg$factors), build_endpoint), names(cfg$factors))
}

# Function to assign participants to event-stratified cross-validation folds.
make_folds <- function(time, event, k, seed) {
  set.seed(seed)
  strata <- interaction(event, time <= 1, drop = TRUE)
  fold <- integer(length(time))
  for (s in levels(strata)) {
    ix <- which(strata == s)
    fold[ix] <- sample(rep(seq_len(k), length.out = length(ix)))
  }
  fold
}

# Fit imputation on training-fold clinical predictors and apply it to the holdout.
# Function to impute missing predictors using the training fold.
impute_split <- function(train, test, cfg, seed) {
  if (is.null(test)) {
    test <- train[FALSE, , drop = FALSE]
  }
  dat <- rbind(train, test)
  init <- dat
  set.seed(seed)
  for (nm in names(dat)) {
    donors <- train[[nm]][!is.na(train[[nm]])]
    missing <- which(is.na(dat[[nm]]))
    if (length(missing)) {
      init[[nm]][missing] <- sample(donors, length(missing), replace = TRUE)
    }
  }
  method <- setNames(ifelse(colSums(is.na(dat)) > 0, "rf", ""), names(dat))
  pm <- matrix(
    1L,
    ncol(dat),
    ncol(dat),
    dimnames = list(names(dat), names(dat))
  )
  diag(pm) <- 0L
  # Determine unusable predictor columns from training data only.
  constant <- vapply(
    train,
    function(v) length(unique(v[!is.na(v)])) < 2L,
    logical(1)
  )
  pm[, constant] <- 0L
  pm[method == "", ] <- 0L
  mi <- mice(
    dat,
    m = cfg$m,
    maxit = cfg$iterations,
    method = method,
    predictorMatrix = pm,
    data.init = init,
    ignore = c(rep(FALSE, nrow(train)), rep(TRUE, nrow(test))),
    seed = seed,
    ntree = cfg$trees,
    rfPackage = "ranger",
    num.threads = 1L,
    remove.constant = FALSE,
    remove.collinear = FALSE,
    printFlag = FALSE
  )
  completed <- lapply(seq_len(cfg$m), function(k) complete(mi, k))
  list(data = completed, final_imputer = mi)
}

# Function to construct comparable predictor matrices for each imputation.
matrices <- function(completed) {
  designs <- lapply(completed, function(d) model.matrix(~ ., d))
  first <- designs[[1L]]
  predictors <- colnames(first)[-1L]
  mapping <- setNames(names(completed[[1L]])[attr(first, "assign")[-1L]], predictors)
  x <- lapply(designs, function(d) d[, predictors, drop = FALSE])
  list(x = x, mapping = mapping)
}

# Function to fit a penalized Cox model across the tuning path.
fit_path <- function(x, time, event, alpha, cfg) {
  first <- glmnet(
    x,
    Surv(time, event),
    family = "cox",
    cox.ties = "breslow",
    alpha = alpha,
    nlambda = 30,
    lambda.min.ratio = 0.1,
    standardize = TRUE,
    control = list(thresh = 1e-07)
  )
  lm <- first$lambda[1]
  fit <- glmnet(
    x,
    Surv(time, event),
    family = "cox",
    cox.ties = "breslow",
    alpha = alpha,
    lambda = lm * cfg$ratios,
    standardize = TRUE,
    control = list(thresh = 1e-07)
  )
  list(
    beta = as.matrix(coef(fit, s = lm * cfg$ratios)),
    lambda = lm * cfg$ratios
  )
}

# Function to calculate the validation loss per event.
cox_loss <- function(time, event, lp) {
  lp <- as.matrix(lp)
  ord <- order(time, decreasing = TRUE)
  time <- time[ord]
  event <- event[ord]
  lp <- lp[ord, , drop = FALSE]
  lp <- sweep(lp, 2, apply(lp, 2, max))
  risk <- apply(exp(lp), 2, cumsum)
  if (is.null(dim(risk))) {
    risk <- matrix(risk, ncol = 1)
  }
  last <- length(time) - match(time, rev(time)) + 1L
  -2 *
    colSums(
      lp[event == 1, , drop = FALSE] -
        log(risk[last[event == 1], , drop = FALSE])
    ) /
    sum(event)
}

# Function to tune the clinical model within training folds.
tune_clinical <- function(z, time, event, cfg, seed) {
  folds <- make_folds(time, event, cfg$inner_k, seed)
  losses <- vector("list", cfg$inner_k)
  for (f in seq_len(cfg$inner_k)) {
    tr <- which(folds != f)
    va <- which(folds == f)
    im <- impute_split(
      z[tr, , drop = FALSE],
      z[va, , drop = FALSE],
      cfg,
      seed + f
    )
    mx <- matrices(im$data)$x
    # Evaluate penalty settings on the held-out inner fold.
    losses[[f]] <- rbindlist(lapply(cfg$alphas, function(a) {
      loss <- sapply(mx, function(x) {
        fit <- fit_path(
          x[seq_along(tr), , drop = FALSE],
          time[tr],
          event[tr],
          a,
          cfg
        )
        cox_loss(
          time[va],
          event[va],
          x[length(tr) + seq_along(va), , drop = FALSE] %*% fit$beta
        )
      })
      data.table(
        alpha = a,
        index = seq_along(cfg$ratios),
        fold = f,
        loss = rowMeans(loss)
      )
    }))
  }
  cv <- rbindlist(losses)[,
    .(loss = mean(loss), se = sd(loss) / sqrt(.N)),
    by = .(alpha, index)
  ]
  best <- cv[which.min(loss)]
  eligible <- cv[alpha == best$alpha & loss <= best$loss + best$se]
  one <- eligible[which.min(index)]
  list(
    alpha = best$alpha,
    min_index = best$index,
    one_index = one$index,
    cv = cv
  )
}

# Function to identify predictors with nonzero penalized coefficients across training imputations.
select_clinical <- function(mx, time, event, tuning, cfg) {
  paths <- lapply(
    mx$x,
    fit_path,
    time = time,
    event = event,
    alpha = tuning$alpha,
    cfg = cfg
  )
  # Function to calculate selection frequencies across imputed training sets.
  frequency <- function(index) {
    selected <- lapply(paths, function(p) {
      unique(unname(mx$mapping[rownames(p$beta)[abs(p$beta[, index]) > 1e-08]]))
    })
    vars <- unique(unname(mx$mapping))
    setNames(
      vapply(
        vars,
        function(v) mean(vapply(selected, function(s) v %in% s, logical(1))),
        numeric(1)
      ),
      vars
    )
  }
  freq <- frequency(tuning$one_index)
  rule <- "lambda.1se"
  chosen <- names(freq)[freq >= cfg$selection_fraction]
  if (!length(chosen)) {
    freq <- frequency(tuning$min_index)
    rule <- "lambda.min"
    chosen <- names(freq)[freq >= cfg$selection_fraction]
  }
  list(variables = chosen, frequency = freq, rule = rule)
}

# Function to fit a Cox risk model with optional Firth correction.
fit_risk <- function(
  train,
  test,
  time,
  event,
  keep_model = FALSE,
  firth = FALSE
) {
  vary <- if (ncol(train)) {
    vapply(
      seq_len(ncol(train)),
      function(j) length(unique(train[, j])) > 1,
      logical(1)
    )
  } else {
    logical()
  }
  train <- train[, vary, drop = FALSE]
  test <- test[, vary, drop = FALSE]
  audit <- NULL
  if (ncol(train)) {
    warnings <- character()
    if (firth) {
      d <- data.frame(time = time, event = event, train, check.names = TRUE)
      fit <- withCallingHandlers(
        coxphf::coxphf(
          Surv(time, event) ~ .,
          data = d,
          pl = FALSE,
          maxit = 1000,
          maxstep = 0.25
        ),
        warning = function(w) stop(conditionMessage(w), call. = FALSE)
      )
    } else {
      fit <- withCallingHandlers(
        coxph(
          Surv(time, event) ~ train,
          ties = "breslow",
          singular.ok = FALSE,
          control = coxph.control(iter.max = 100)
        ),
        warning = function(w) {
          warnings <<- c(warnings, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      if (length(warnings)) {
        if (
          any(!grepl("coefficient may be infinite", warnings, fixed = TRUE))
        ) {
          stop(
            "Unexpected Cox warning: ",
            paste(unique(warnings), collapse = "; ")
          )
        }
        # Identify predictor columns implicated in the Cox convergence warnings.
        audit <- rbindlist(lapply(unique(warnings), function(w) {
          indices <- as.integer(strsplit(
            trimws(sub(".*variable +([0-9 ,]+).*", "\\1", w)),
            "[ ,]+"
          )[[1]])
          rbindlist(lapply(indices, function(j) {
            v <- train[, j]
            binary <- all(v %in% c(0, 1))
            data.table(
              variable = colnames(train)[j],
              coefficient = unname(coef(fit)[j]),
              n_train = length(time),
              events = sum(event),
              binary = binary,
              n_zero = if (binary) {
                sum(v == 0)
              } else {
                NA_integer_
              },
              events_zero = if (binary) {
                sum(event[v == 0])
              } else {
                NA_integer_
              },
              n_one = if (binary) {
                sum(v == 1)
              } else {
                NA_integer_
              },
              events_one = if (binary) {
                sum(event[v == 1])
              } else {
                NA_integer_
              },
              warning = w
            )
          }))
        }))
        return(list(
          risk = rep(NA_real_, nrow(test)),
          model = NULL,
          audit = audit
        ))
      }
    }
    beta <- coef(fit)
    if (any(!is.finite(beta))) {
      stop("Non-finite Cox coefficients")
    }
    lp <- drop(train %*% beta)
    newlp <- drop(test %*% beta)
  } else {
    beta <- numeric()
    lp <- rep(0, length(time))
    newlp <- rep(0, nrow(test))
  }
  center <- mean(lp)
  r <- exp(lp - center)
  times <- sort(unique(time[event == 1 & time <= 1]))
  h <- sum(vapply(
    times,
    function(tt) sum(event == 1 & time == tt) / sum(r[time >= tt]),
    numeric(1)
  ))
  p <- -expm1(-h * exp(newlp - center))
  list(
    risk = p,
    model = if (keep_model) {
      list(
        columns = colnames(train),
        beta = beta,
        center = center,
        H0_12m = h,
        estimator = if (firth) "Firth" else "ordinary Cox"
      )
    } else {
      NULL
    },
    audit = audit
  )
}



# ==============================================================================
# FIT AND PREDICT WITHIN A FOLD
# ==============================================================================

# Function to fit models on the training fold and predict held-out risks.
run_fold <- function(job, cfg) {
  x <- inputs[[job$endpoint]]
  n <- nrow(x$z)
  seed <- cfg$seed +
    if (job$endpoint == "MACCE") {
      100000L
    } else {
      0L
    }
  seed <- seed + job$repeat_id * 1000L + job$fold * 10L
  folds <- make_folds(x$time, x$event, cfg$outer_k, cfg$seed + job$repeat_id)
  full <- job$fold == 0L
  tr <- if (full) {
    seq_len(n)
  } else {
    which(folds != job$fold)
  }
  te <- if (full) {
    integer()
  } else {
    which(folds == job$fold)
  }
  tuning <- tune_clinical(
    x$z[tr, , drop = FALSE],
    x$time[tr],
    x$event[tr],
    cfg,
    seed
  )
  im <- impute_split(
    x$z[tr, , drop = FALSE],
    x$z[te, , drop = FALSE],
    cfg,
    seed + 8L
  )
  mx <- matrices(im$data)
  train_mx <- list(
    x = lapply(mx$x, function(a) a[seq_along(tr), , drop = FALSE]),
    mapping = mx$mapping
  )
  sel <- select_clinical(train_mx, x$time[tr], x$event[tr], tuning, cfg)
  cols <- names(mx$mapping)[mx$mapping %in% sel$variables]
  risks <- array(NA_real_, c(length(te), 3, cfg$m))
  models <- vector("list", cfg$m)
  audit <- list()
  use_firth <- FALSE
  rowix <- c(tr, te)
  extras <- list(
    matrix(numeric(), nrow = length(rowix), ncol = 0),
    x$cluster[rowix, , drop = FALSE],
    x$factors[rowix, , drop = FALSE]
  )
  for (pass in 1:2) {
    for (k in seq_len(cfg$m)) {
      a <- mx$x[[k]][, cols, drop = FALSE]
      models[[k]] <- vector("list", 3)
      for (j in 1:3) {
        design <- cbind(a, extras[[j]])
        result <- fit_risk(
          design[seq_along(tr), , drop = FALSE],
          design[length(tr) + seq_along(te), , drop = FALSE],
          x$time[tr],
          x$event[tr],
          full,
          firth = use_firth
        )
        if (!is.null(result$audit)) {
          v <- result$audit
          v[, `:=`(
            endpoint = job$endpoint,
            repeat_id = job$repeat_id,
            fold = job$fold,
            imputation = k,
            model = c("Clinical", "Clinical + C1/C2", "Clinical + Factors")[j]
          )]
          audit[[length(audit) + 1L]] <- v
        }
        risks[, j, k] <- result$risk
        models[[k]][[j]] <- result$model
      }
    }
    if (!length(audit) || use_firth) {
      break
    }
    use_firth <- TRUE
    # Record convergence diagnostics before the Firth fit.

    fwrite(
      rbindlist(audit),
      file.path(
        "./results/prediction",
        paste0(
          "Convergence_",
          job$endpoint,
          "_",
          job$repeat_id,
          "_",
          job$fold,
          ".csv"
        )
      )
    )
  }
  list(
    endpoint = job$endpoint,
    repeat_id = job$repeat_id,
    fold = job$fold,
    rows = te,
    predictions = if (length(te)) apply(risks, c(1, 2), mean) else NULL,
    selection = sel,
    tuning = tuning,
    estimator = if (use_firth) "Firth" else "ordinary Cox",
    convergence_audit = rbindlist(audit),
    final_model = if (full) {
      list(models = models, imputer = im$final_imputer)
    } else {
      NULL
    }
  )
}

# Function to calculate inverse probability of censoring weights at 12 months.
ipcw <- function(time, event) {
  ord <- order(time)
  g <- pec::ipcw(
    Surv(time, status) ~ 1,
    data = data.frame(time = time[ord], status = event[ord]),
    method = "marginal",
    times = 1,
    subjectTimes = time[ord],
    subjectTimesLag = 1
  )
  gt <- numeric(length(time))
  gt[ord] <- g$IPCW.subjectTimes
  y <- as.numeric(event == 1 & time <= 1)
  w <- numeric(length(time))
  w[y == 1] <- 1 / gt[y == 1]
  w[time > 1] <- 1 / as.numeric(g$IPCW.times)
  list(y = y, w = w)
}
# Function to calculate the sensitivity and specificity at each risk threshold.
roc_points <- function(p, q) {
  d <- data.table(p = p, yes = q$w * q$y, no = q$w * (1 - q$y))
  d <- d[, .(yes = sum(yes), no = sum(no)), by = p][order(-p)]
  data.table(
    FPR = c(0, cumsum(d$no) / sum(d$no)),
    TPR = c(0, cumsum(d$yes) / sum(d$yes))
  )
}
# Function to calculate the area under the ROC curve.
auc_value <- function(p, q) {
  r <- roc_points(p, q)
  sum(diff(r$FPR) * (head(r$TPR, -1) + tail(r$TPR, -1)) / 2)
}
# Function to count weighted risk-category changes by event status.
reclassification_table <- function(old_risk, new_risk, q, threshold) {
  rbindlist(lapply(0:1, function(s) {
    d <- data.table(
      old = as.integer(old_risk >= threshold) + 1L,
      new = as.integer(new_risk >= threshold) + 1L,
      weight = q$w * (q$y == s),
      observed = as.integer(q$w > 0 & q$y == s)
    )
    v <- merge(
      CJ(old = 1:2, new = 1:2),
      d[,
        .(IPCW_weighted_n = sum(weight), observed_n = sum(observed)),
        by = .(old, new)
      ],
      all.x = TRUE
    )
    v[is.na(IPCW_weighted_n), `:=`(IPCW_weighted_n = 0, observed_n = 0)]
    v[,
      status := if (s == 1) {
        "Event"
      } else {
        "Non-event"
      }
    ]
    v
  }))
}

# Function to calculate prediction performance at the specified time horizon.
metrics <- function(p, q) {
  lp <- qlogis(pmin(pmax(p, 1e-07), 1 - 1e-07))
  d <- data.frame(y = q$y, w = q$w, lp = lp)
  h <- glm(y ~ offset(lp), data = d, weights = w, family = quasibinomial())
  g <- glm(y ~ lp, data = d, weights = w, family = quasibinomial())
  c(
    AUC = auc_value(p, q),
    Brier = sum(q$w * (q$y - p)^2) / length(p),
    Intercept = unname(coef(h)[1]),
    Slope = unname(coef(g)[2])
  )
}
# Function to calculate event and non-event reclassification at the risk threshold.
nri_value <- function(old, new, q, threshold) {
  move <- as.numeric(new >= threshold) - as.numeric(old >= threshold)
  a <- sum(q$w * q$y * move) / sum(q$w * q$y)
  b <- -sum(q$w * (1 - q$y) * move) / sum(q$w * (1 - q$y))
  c(Event = a, Non_event = b, Total = a + b)
}

# Function to estimate the 12-month event rate.
event_rate_12m <- function(time, event) {
  unname(
    1 - summary(survfit(Surv(time, event) ~ 1), times = 1, extend = FALSE)$surv
  )
}

# Function to resample participants to quantify uncertainty in prediction metrics.
bootstrap_metrics <- function(time, event, predictions, threshold, B, seed) {
  n <- length(time)
  # Function to compute the statistic for one resample.
  stat <- function(ix) {
    q <- ipcw(time[ix], event[ix])
    p <- predictions[ix, , drop = FALSE]
    aa <- unname(apply(p, 2, auc_value, q = q))
    c(
      setNames(aa, paste0("AUC_", 1:3)),
      DeltaAUC_C = aa[2] - aa[1],
      DeltaAUC_F = aa[3] - aa[1],
      setNames(
        nri_value(
          p[,
            1
          ],
          p[, 3],
          q,
          threshold
        ),
        c("NRI_Event", "NRI_Non_event", "NRI_Total")
      )
    )
  }
  estimate <- stat(seq_len(n))
  set.seed(seed)
  boot <- t(replicate(B, stat(sample.int(n, n, replace = TRUE))))
  data.table(
    metric = names(estimate),
    estimate = unname(estimate),
    lower = apply(boot, 2, quantile, 0.025),
    upper = apply(boot, 2, quantile, 0.975),
    # Compute bootstrap P-values for the discrimination and reclassification contrasts.
    P = sapply(seq_along(estimate), function(j) {
      if (j <= 3) {
        NA_real_
      } else {
        (1 + sum(abs(boot[, j] - estimate[j]) >= abs(estimate[j]))) / (B + 1)
      }
    })
  )
}



# ==============================================================================
# EVALUATE ENDPOINT PREDICTIONS | FIGURE 6
# ==============================================================================

# Function to evaluate the cross-validated risk predictions for one endpoint.
evaluate_endpoint <- function(endpoint, state, cfg) {
  x <- inputs[[endpoint]]
  n <- nrow(x$z)
  mods <- c("Clinical", "Clinical + C1/C2", "Clinical + Factors")
  runs <- Filter(function(v) v$endpoint == endpoint && v$fold > 0L, state)
  predictions <- array(NA_real_, c(n, 3, cfg$repeats))
  # Collect held-out predictions and selected predictors across folds.
  selection <- rbindlist(lapply(runs, function(v) {
    predictions[v$rows, , v$repeat_id] <<- v$predictions
    data.table(
      repeat_id = v$repeat_id,
      fold = v$fold,
      rule = v$selection$rule,
      variable = names(v$selection$frequency),
      mi_fraction = unname(v$selection$frequency),
      selected = names(v$selection$frequency) %in% v$selection$variables
    )
  }))
  p <- apply(predictions, c(1, 2), mean)
  colnames(p) <- mods
  q <- ipcw(x$time, x$event)
  # Calculate performance within each cross-validation repeat.
  metric_table <- rbindlist(lapply(seq_len(cfg$repeats), function(r) {
    rbindlist(lapply(1:3, function(j) {
      cbind(
        data.table(repeat_id = r, model = mods[j]),
        as.data.table(as.list(metrics(predictions[, j, r], q)))
      )
    }))
  }))
  # Summarize performance from the repeat-averaged predictions.
  aggregate <- rbindlist(lapply(1:3, function(j) {
    cbind(
      data.table(model = mods[j]),
      as.data.table(as.list(metrics(
        p[,
          j
        ],
        q
      )))
    )
  }))
  roc <- rbindlist(lapply(1:3, function(j) {
    v <- roc_points(p[, j], q)
    v[, model := mods[j]]
    v
  }))
  # Group predicted risks to estimate observed event rates for calibration.
  cal <- rbindlist(lapply(1:3, function(j) {
    breaks <- unique(quantile(p[, j], seq(0, 1, 0.2)))
    grp <- if (length(breaks) > 1) {
      cut(p[, j], breaks, include.lowest = TRUE, labels = FALSE)
    } else {
      rep(1L, n)
    }
    rbindlist(lapply(sort(unique(grp)), function(k) {
      ix <- which(grp == k)
      sf <- summary(
        survfit(Surv(x$time[ix], x$event[ix]) ~ 1, conf.type = "log-log"),
        times = 1,
        extend = FALSE
      )
      data.table(
        model = mods[j],
        group = k,
        n = length(ix),
        predicted = mean(p[ix, j]),
        observed = 1 - sf$surv[[1L]],
        lower = 1 - sf$upper[[1L]],
        upper = 1 - sf$lower[[1L]]
      )
    }))
  }))
  # Calculate net benefit across risk thresholds.
  dca <- rbindlist(lapply(seq(0.01, 0.3, 0.0025), function(th) {
    nb <- sapply(1:3, function(j) {
      sum(q$w * (p[, j] >= th) * (q$y - (1 - q$y) * th / (1 - th))) / n
    })
    data.table(
      threshold = th,
      model = c(mods, "Treat all", "Treat none"),
      net_benefit = c(
        nb,
        sum(
          q$w *
            (q$y -
              (1 - q$y) * th / (1 - th))
        ) /
          n,
        0
      )
    )
  }))
  th <- cfg$thresholds[[endpoint]]
  nr <- nri_value(p[, 1], p[, 3], q, th)
  table <- reclassification_table(p[, 1], p[, 3], q, th)
  boot <- bootstrap_metrics(
    x$time,
    x$event,
    p,
    th,
    cfg$bootstrap,
    cfg$seed +
      if (endpoint == "MALE") {
        700000L
      } else {
        800000L
      }
  )
  rate <- event_rate_12m(x$time, x$event)
  final <- state[[paste(endpoint, 0, 0, sep = "_")]]
  list(
    endpoint = endpoint,
    models = mods,
    predictions = p,
    time = x$time,
    event = x$event,
    predictions_by_repeat = predictions,
    metrics = aggregate,
    metrics_by_repeat = metric_table,
    roc = roc,
    calibration = cal,
    dca = dca,
    nri = nr,
    reclassification = table,
    bootstrap = boot,
    threshold = th,
    event_rate = rate,
    selection = selection,
    final_selection = final$selection,
    source = x$source
  )
}



# ==============================================================================
# EXPORT PREDICTION RESULTS | FIGURE 6
# ==============================================================================

# Function to export the prediction results for the manuscript.
save_results <- function(state, cfg) {
  results <- setNames(
    lapply(c("MALE", "MACCE"), evaluate_endpoint, state = state, cfg = cfg),
    c("MALE", "MACCE")
  )

  saveRDS(results, "./results/prediction/Results.rds")
  for (ep in names(results)) {
    z <- results[[ep]]
    cat(ep, "FINAL:", paste(z$final_selection$variables, collapse = ", "), "\n")
  }
  invisible(results)
}
# Function to run a prediction job and save its fold results.
work <- function(i) {
  j <- as.list(jobs[i, ])
  key <- paste(j$endpoint, j$repeat_id, j$fold, sep = "_")
  f <- file.path(
    "./results/prediction",
    paste0("Job_", key, ".rds")
  )
  if (!file.exists(f)) {
    a <- run_fold(j, cfg)

    saveRDS(a, paste0(f, ".partial"))
    stopifnot(file.rename(paste0(f, ".partial"), f))
  }
  cat(format(Sys.time()), key, "complete\n")
  key
}
# Function to coordinate the prediction folds and aggregate their predictions.
run <- function() {
  keys <- with(jobs, paste(endpoint, repeat_id, fold, sep = "_"))
  todo <- which(
    !file.exists(file.path(
      "./results/prediction",
      paste0("Job_", keys, ".rds")
    ))
  )
  if (length(todo)) {
    cl <- parallel::makePSOCKcluster(
      min(length(todo), cfg$workers),
      outfile = "./results/prediction/Workers.log"
    )
    on.exit(try(parallel::stopCluster(cl), silent = TRUE))
    parallel::clusterEvalQ(cl, {
      library(data.table)
      library(mice)
      library(survival)
      library(glmnet)
      setDTthreads(1L)
      NULL
    })
    parallel::clusterExport(
      cl,
      c(
        ls(envir = .GlobalEnv)[vapply(
          mget(ls(envir = .GlobalEnv), .GlobalEnv),
          is.function,
          logical(1)
        )],
        "cfg",
        "jobs",
        "inputs"
      )
    )
    parallel::parLapplyLB(cl, todo, work, chunk.size = 1)
    parallel::stopCluster(cl)
    on.exit(NULL)
  }
  state <- setNames(
    lapply(
      file.path(
        "./results/prediction",
        paste0("Job_", keys, ".rds")
      ),
      readRDS
    ),
    keys
  )
  for (ep in c("MALE", "MACCE")) {
    k <- paste(ep, 0, 0, sep = "_")

    saveRDS(
      state[[k]]$final_model,
      file.path(
        "./results/prediction",
        paste0(ep, "_FinalModel.rds")
      )
    )
    state[[k]]$final_model <- NULL
  }
  save_results(state, cfg)
}



# ==============================================================================
# MODEL SETTINGS
# ==============================================================================

config <- list(
  seed = 20261023L,
  outer_k = 5L,
  repeats = 10L,
  inner_k = 5L,
  m = 20L,
  iterations = 10L,
  trees = 10L,
  workers = as.integer(Sys.getenv("VIP_PRED_WORKERS", "8")),
  alphas = c(0.1, 0.5, 1),
  ratios = exp(seq(0, log(0.001), length.out = 40)),
  selection_fraction = 0.5,
  bootstrap = 2000L,
  factors = list(MALE = c("F1", "F8"), MACCE = c("F2", "F3", "F4", "F8")),
  thresholds = c(MALE = 0.16, MACCE = 0.1)
)



# ==============================================================================
# RUN ANALYSIS
# ==============================================================================

cfg <- config
inputs <- prediction_inputs(cfg)
dir.create("./results/prediction", recursive = TRUE, showWarnings = FALSE)
jobs <- expand.grid(
  fold = 1:5,
  endpoint = c("MALE", "MACCE"),
  repeat_id = 1:10,
  stringsAsFactors = FALSE
)
jobs <- rbind(
  jobs,
  data.frame(fold = 0L, endpoint = c("MALE", "MACCE"), repeat_id = 0L)
)
run()



# ==============================================================================
# SINGLE-FACTOR DISCRIMINATION | FIGURE S25 (PANELS A–B)
# ==============================================================================
data.table::setDTthreads(1L)
allout <- list()
for (ep in c("MALE", "MACCE")) {
  x <- inputs[[ep]]
  fac <- config$factors[[ep]]
  n <- length(x$time)
  pred <- array(
    NA_real_,
    c(n, length(fac), 10),
    dimnames = list(NULL, fac, NULL)
  )
  for (repeat_id in 1:10) {
    folds <- make_folds(x$time, x$event, 5, config$seed + repeat_id)
    for (f in 1:5) {
      tr <- which(folds != f)
      te <- which(folds == f)
      for (j in seq_along(fac)) {
        z <- fit_risk(
          x$factors[tr, fac[j], drop = FALSE],
          x$factors[te, fac[j], drop = FALSE],
          x$time[tr],
          x$event[tr]
        )
        pred[te, j, repeat_id] <- z$risk
      }
    }
  }
  p <- apply(pred, c(1, 2), mean)
  q <- ipcw(x$time, x$event)
  # Summarize discrimination for each individual factor.
  mm <- rbindlist(lapply(seq_along(fac), function(j) {
    data.table(
      factor = fac[j],
      AUC = auc_value(p[, j], q)
    )
  }))
  rr <- rbindlist(lapply(seq_along(fac), function(j) {
    z <- roc_points(p[, j], q)
    z[, factor := fac[j]]
    z
  }))
  allout[[ep]] <- list(roc = rr, metrics = mm)
}



# ==============================================================================
# SAVE OUTPUT
# ==============================================================================
saveRDS(allout, "./results/prediction/SingleFactor_OOF.rds")



# ==============================================================================
# RECLASSIFICATION WITH CLUSTER MEMBERSHIP | FIGURE 6 (PANEL C)
# ==============================================================================

results <- readRDS("./results/prediction/Results.rds")
cluster_nri <- results
for (ep in names(cluster_nri)) {
  z <- cluster_nri[[ep]]
  p <- z$predictions
  th <- config$thresholds[[ep]]
  q <- ipcw(z$time, z$event)
  z$nri <- nri_value(p[, 1], p[, 2], q, th)
  z$reclassification <- reclassification_table(p[, 1], p[, 2], q, th)
  # Function to compute the statistic for one resample.
  stat <- function(ix) {
    nri_value(p[ix, 1], p[ix, 2], ipcw(z$time[ix], z$event[ix]), th)
  }
  set.seed(
    config$seed +
      if (ep == "MALE") {
        700000L
      } else {
        800000L
      }
  )
  n <- length(z$time)
  B <- 2000L
  boot <- t(replicate(B, stat(sample.int(n, n, replace = TRUE))))
  z$bootstrap <- data.table(
    metric = c("NRI_Event", "NRI_Non_event", "NRI_Total"),
    estimate = unname(z$nri),
    lower = apply(boot, 2, quantile, 0.025),
    upper = apply(boot, 2, quantile, 0.975),
    P = vapply(
      1:3,
      function(j) {
        (1 +
          sum(abs(boot[, j] - z$nri[j]) >= abs(z$nri[j]))) /
          (B + 1)
      },
      numeric(1)
    )
  )
  z$threshold <- th
  cluster_nri[[ep]] <- z

  cat(ep, "NRI", z$nri, "P", z$bootstrap$P[3], "\n")
}
saveRDS(cluster_nri, "./results/prediction/Cluster_NRI.rds")



# ==============================================================================
# PRINT SUPPLEMENTARY TABLES S23 AND S24
# ==============================================================================
labels <- rbindlist(list(
  list("Age, years", "Continuous", "Age"),
  list("Sex", "Female / Male (ref: Female)", "Sex"),
  list("Body mass index, kg/m²", "Continuous", "BMI"),
  list(
    "Current smoking",
    "Non-current / Current (ref: Non-current)",
    "SmokingStatus"
  ),
  list("Disease duration, years", "Continuous", "DiseaseDuration_Years_Calc"),
  list(
    "Prior PAD-related procedure",
    "No / Yes (ref: No)",
    "PriorPADProcedure"
  ),
  list("Rutherford category", "0-3 / 4-6 (ref: 0-3)", "RutherfordCategory"),
  list("Barthel index, points", "Continuous", "BarthelTotal"),
  list(
    "Disease course at admission",
    "Acute onset or exacerbation / Chronic / Subacute (ref: Acute onset or exacerbation)",
    "DiseaseCourse_Grouped"
  ),
  list("History of hypertension", "No / Yes (ref: No)", "Hypertension"),
  list("History of diabetes", "No / Yes (ref: No)", "DiabetesCategory"),
  list(
    "History of coronary artery disease",
    "No / Yes (ref: No)",
    "HeartDisease"
  ),
  list("History of stroke", "No / Yes (ref: No)", "CerebrovascularDisease"),
  list("Recorded history of CKD", "No / Yes (ref: No)", "CKDStage_Recorded"),
  list("History of dyslipidemia", "No / Yes (ref: No)", "Dyslipidemia"),
  list(
    "Preoperative antihypertensive therapy",
    "No / Yes (ref: No)",
    "Antihypertensive"
  ),
  list(
    "Preoperative lipid-lowering therapy",
    "No / Yes (ref: No)",
    "LipidLowering"
  ),
  list(
    "Preoperative antithrombotic therapy",
    "No / Yes (ref: No)",
    "Antithrombotic"
  ),
  list("Preoperative vasodilator therapy", "No / Yes (ref: No)", "Vasodilator"),
  list("Systolic blood pressure, mmHg", "Continuous", "SBP"),
  list("Diastolic blood pressure, mmHg", "Continuous", "DBP"),
  list("Fasting glucose, mmol/L", "Continuous", "FastingGlucose"),
  list("Triglycerides, mmol/L", "Continuous", "Triglycerides"),
  list("HDL cholesterol, mmol/L", "Continuous", "HDLCholesterol"),
  list("LDL cholesterol, mmol/L", "Continuous", "LDLCholesterol"),
  list(
    "eGFR (CKD-EPI 2021), mL/min/1.73 m²",
    "Continuous",
    "eGFR_CKD_EPI_2021"
  ),
  list("Hemoglobin, g/L", "Continuous", "Hemoglobin"),
  list("Albumin, g/L", "Continuous", "Albumin"),
  list("Mean corpuscular volume (MCV), fL", "Continuous", "MCV"),
  list(
    "Mean corpuscular hemoglobin concentration (MCHC), g/L",
    "Continuous",
    "MCHC"
  ),
  list("Red cell distribution width (RDW-CV), %", "Continuous", "RDW_CV"),
  list("Platelet count, 10⁹/L", "Continuous", "PlateletCount"),
  list("Neutrophil count, 10⁹/L", "Continuous", "NeutrophilCount"),
  list("Lymphocyte count, 10⁹/L", "Continuous", "LymphocyteCount"),
  list("Monocyte count, 10⁹/L", "Continuous", "MonocyteCount"),
  list("Total protein, g/L", "Continuous", "TotalProtein"),
  list("Globulin, g/L", "Continuous", "Globulin"),
  list("Alanine aminotransferase (ALT), U/L", "Continuous", "ALT"),
  list("Aspartate aminotransferase (AST), U/L", "Continuous", "AST"),
  list("Urea, mmol/L", "Continuous", "Urea"),
  list("Uric acid, µmol/L", "Continuous", "Urate"),
  list("Total cholesterol, mmol/L", "Continuous", "Cholesterol")
))
setnames(labels, c("Variable", "Coding", "variable"))
ids <- c(MALE = "S23", MACCE = "S24")
for (ep in names(ids)) {
  x <- inputs[[ep]]
  z <- results[[ep]]
  vars <- names(x$z)
  tab <- labels[match(vars, variable)]
  tab[, missing_n := as.integer(colSums(is.na(x$z)))]
  tab[,
    selected_folds := z$selection[, .(n = sum(selected)), by = variable][
      match(vars, variable),
      n
    ]
  ]
  tab[, final_selected := vars %in% z$final_selection$variables]
  tab[, candidate_order := match(variable, labels$variable)]
  setorder(tab, -selected_folds, candidate_order)
  folds <- cfg$outer_k * cfg$repeats
  tab <- tab[, .(
    Variable, Coding,
    `Missing before imputation, n (%)` = sprintf(
      "%d (%.1f%%)", missing_n, 100 * missing_n / nrow(x$z)
    ),
    Selected = sprintf("%d (%.0f%%)", selected_folds, 100 * selected_folds / folds),
    `Final model` = ifelse(final_selected, "Yes", "No")
  )]
  setnames(tab, "Selected", sprintf("Selected /%d (%%)", folds))
  cat("\n", ids[[ep]], " (", ep, ")\n", sep = "")
  print(tab, nrows = Inf, trunc.cols = FALSE)
}
