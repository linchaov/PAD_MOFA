# ==============================================================================
#  METADATA
# ==============================================================================
#
# Author: Chaoyang Lin, MD
# Contact: chaoyanglint@163.com
# License: MIT
#
# Author's environment: R 4.6.1, Windows 11 x64.
# R packages used: clue 0.3-68, cluster 2.1.8.2, data.table 1.18.6.1, digest 0.6.39,
#                  filelock 1.0.3, MOFA2 1.22.1, reticulate 1.47.0
# Python backend: Python 3.12.14, mofapy2 0.7.5.
#
# Global: Cohort-specific MOFA fitting, factor stability and molecular clustering.
#
# Inputs: ./Datasets_PAD_MOFA/VIP.rds and ./Datasets_PAD_MOFA/UKB.rds.
#
# Outputs: Models, factor scores, loadings and cluster assignments under
#          ./results/VIP-PAD/ and ./results/UKB/.
#
# ==============================================================================



# ==============================================================================
# SETUP
# ==============================================================================

# Sys.setenv(RETICULATE_PYTHON = "/your/path/to/python")

suppressPackageStartupMessages({
  library(MOFA2)
  library(data.table)
  library(digest)
  library(clue)
  library(cluster)
})



# ==============================================================================
# COHORT MOFA FITTING AND CLUSTERING | FIGURES S02 AND S03
# ==============================================================================

# Function to fit cohort-specific MOFA models and assess cluster stability.
fit_cohort <- function(cohort_name, seed = NULL) {
  base <- file.path("./results", cohort_name)
  out <- file.path(base, "mofa")
  vip <- cohort_name == "VIP-PAD"
  set.seed(42)

  dir.create(
    file.path(base, "mofa/runs"),
    recursive = TRUE,
    showWarnings = FALSE
  )

  # Function to fit one initialization or dispatch the parallel seed runs.
  train_models <- function(seed = NULL) {
    threads <- if (vip) 1L else 2L
    Sys.setenv(OMP_NUM_THREADS = threads, OPENBLAS_NUM_THREADS = threads,
               MKL_NUM_THREADS = threads)
    if (vip) Sys.setenv(NUMEXPR_NUM_THREADS = 1L)
    input <- if (vip) "./Datasets_PAD_MOFA/VIP.rds" else "./Datasets_PAD_MOFA/UKB.rds"
    seeds <- 42:141
    contract <- list(
      input_sha256 = digest(file = input, algo = 'sha256'),
      seeds = seeds,
      data = list(scale_views = TRUE, use_float32 = FALSE),
      model = list(num_factors = 25L, spikeslab_weights = TRUE),
      training = list(maxiter = 10000L, convergence_mode = "slow",
                      drop_factor_threshold = .02),
      MOFA2_version = as.character(packageVersion('MOFA2')),
      workers = if (vip) 12L else 8L,
      threads_per_worker = threads
    )
    target <- file.path(base, "mofa/Training_Contract.rds")
    if (is.null(seed) && !file.exists(target)) {

      saveRDS(contract, target)
    } else {
      stopifnot(identical(contract, readRDS(target)))
    }

    if (is.null(seed)) {
      cluster <- parallel::makePSOCKcluster(contract$workers)
      on.exit(parallel::stopCluster(cluster), add = TRUE)
      parallel::clusterExport(
        cluster,
        c('out', 'cohort_name'),
        envir = environment()
      )
      parallel::parLapplyLB(cluster, seeds, function(seed) {
        # Resume completed initializations for the same input and model settings.
        if (
          file.exists(file.path(out, sprintf('runs/seed_%03d_fit.rds', seed)))
        ) {
          return(0L)
        }
        system2(
          file.path(R.home('bin'), 'Rscript'),
          c("code/01-mofa.R", "--seed", cohort_name, seed),
          stdout = file.path(
            out, 'runs',
            sprintf('seed_%03d.log', seed)
          ),
          stderr = ''
        )
      })
    } else {
      prefix <- file.path(out, sprintf('runs/seed_%03d', seed))
      lock <- filelock::lock(paste0(prefix, "_training.lock"))
      on.exit(filelock::unlock(lock), add = TRUE)
      complete <- paste0(prefix, "_fit.rds")
      if (file.exists(complete)) {
        stopifnot(identical(readRDS(complete)$contract, contract))
        quit(status = 0)
      }
      cohort <- readRDS(input)
      ids <- cohort$study_id
      rownames(cohort) <- ids
      views <- if (vip) {
        # Use standardized VIP-PAD values in the feature order stored for each assay.
        attrs <- attr(cohort, "omics_matrix_attributes")
        # Function to transpose the standardized VIP-PAD measurements into a MOFA view.
        vip_view <- function(view) {
          fields <- names(attrs[[view]][["scaled:scale"]])
          t(as.matrix(cohort[, fields, drop = FALSE]))
        }
        list(
          Protein = vip_view("proteomics"),
          Metabolite = vip_view("metabolomics")
        )
      } else {
        # The UKB dictionary maps source columns to the MOFA feature names.
        # Standardize each feature using its observed values and preserve missingness.
        dictionary <- attr(cohort, "feature_dictionary")
        # Function to build a UKB MOFA view using fields in the feature dictionary.
        ukb_view <- function(view) {
          ann <- dictionary[dictionary$modality == view, , drop = FALSE]
          fields <- as.character(ann$source_field)
          x <- as.matrix(cohort[, fields, drop = FALSE])
          for (j in seq_len(ncol(x))) {
            z <- x[, j]
            s <- sd(z, na.rm = TRUE)
            x[, j] <- (z - mean(z, na.rm = TRUE)) / s
          }
          colnames(x) <- as.character(ann$feature)
          t(x)
        }
        list(
          Protein = ukb_view("Protein"),
          Metabolite = ukb_view("Metabolite")
        )
      }
      model <- create_mofa(views, groups = rep("PAD", length(ids)))
      d <- get_default_data_options(model)
      d[names(contract$data)] <- contract$data
      m <- get_default_model_options(model)
      m[names(contract$model)] <- contract$model
      t <- get_default_training_options(model)
      t[names(contract$training)] <- contract$training
      t$seed <- seed
      t$outfile <- paste0(prefix, "_model.hdf5")
      model <- prepare_mofa(
        model,
        data_options = d,
        model_options = m,
        training_options = t
      )
      model <- run_mofa(model, use_basilisk = FALSE)

      saveRDS(model, paste0(prefix, "_model.rds"))
      fit <- list(
        contract = contract,
        seed = seed,
        Z = get_factors(model)[[1]],
        elbo = tail(na.omit(as.numeric(model@training_stats$elbo)), 1)
      )
      saveRDS(fit, paste0(complete, '.tmp'))
      stopifnot(file.rename(paste0(complete, '.tmp'), complete))
      cat('SEED_COMPLETE', seed, '\n')
    }
  }

  # Function to select the converged model with the highest ELBO.
  select_model <- function() {
    contract <- readRDS(
      file.path(base, "mofa/Training_Contract.rds")
    )
    fits <- lapply(contract$seeds, function(s) {
      readRDS(file.path(out, sprintf('runs/seed_%03d_fit.rds', s)))
    })
    stopifnot(
      all(vapply(fits, function(f) identical(f$contract, contract), logical(1)))
    )
    # Read seed-specific convergence and ELBO before comparing fits.
    info <- rbindlist(lapply(fits, function(f) {
      log <- readLines(
        file.path(
          out, 'runs',
          sprintf('seed_%03d.log', f$seed)
        ),
        warn = FALSE
      )
      elbo <- if (vip) {
        lines <- grep('^Iteration [0-9]+:.* ELBO=', log, value = TRUE)
        as.numeric(sub('.* ELBO=([^,]+),.*', '\\1', tail(lines, 1)))
      } else {
        f$elbo
      }
      data.table(seed = f$seed, elbo = elbo,
                 converged = any(grepl('Converged!', log, fixed = TRUE)))
    }))
    best <- info[converged == TRUE & is.finite(elbo)][order(-elbo, seed)]
    seed <- best$seed[1]
    ref <- fits[[match(seed, contract$seeds)]]
    info[, selected := seed == best$seed[1]]
    info[, delta_ELBO_from_best := best$elbo[1] - elbo]

    fwrite(info, file.path(base, "mofa/Model_Diagnostics.csv"))
    model <- readRDS(file.path(out, sprintf('runs/seed_%03d_model.rds', seed)))
    saveRDS(
      model,
      file.path(base, "mofa/Selected_Model.rds")
    )
    ref <- list(seed = seed, Z = ref$Z, W = MOFA2::get_weights(model),
                variance = MOFA2::get_variance_explained(model))
    variance_pf <- ref$variance$r2_per_factor$PAD
    fwrite(
      data.table(
        factor = if (vip) colnames(ref$Z) else rownames(variance_pf),
        Protein = variance_pf[, 'Protein'],
        Metabolite = variance_pf[, 'Metabolite']
      ),
      file.path(base, "mofa/Variance_Explained.csv")
    )

    saveRDS(
      ref,
      file.path(base, "mofa/Selected_Fit.rds")
    )
    fwrite(
      data.table(SampleID = rownames(ref$Z), ref$Z),
      file.path(base, "mofa/Factor_Scores_Raw.csv")
    )
    # Match the converged runs to the selected model for stability estimates.
    rows <- lapply(fits, function(f) {
      if (f$seed == seed || !info[seed == f$seed, converged]) {
        return(NULL)
      }
      cc <- cor(ref$Z, f$Z)
      k <- ncol(ref$Z)
      kc <- ncol(f$Z)
      cost <- matrix(0, max(k, kc), max(k, kc))
      cost[1:k, 1:kc] <- abs(cc)
      assignment <- as.integer(solve_LSAP(cost, maximum = TRUE))[1:k]
      # Record the correlation for each matched factor pair.
      rbindlist(lapply(1:k, function(i) {
        j <- assignment[i]
        data.table(
          seed = f$seed,
          factor = colnames(ref$Z)[i],
          matched = j <= kc,
          score_abs_r = if (j <= kc) abs(cc[i, j]) else NA_real_
        )
      }))
    })

    fwrite(
      rbindlist(rows),
      file.path(base, "mofa/Factor_Stability_All_Runs.csv")
    )
  }

  # Function to assess cluster stability across subsamples.
  cluster_stability <- function() {
    model <- readRDS(
      file.path(base, "mofa/Selected_Model.rds")
    )
    Z <- get_factors(model)[[1]]
    # Cluster patients using all retained MOFA factors.
    fac <- colnames(Z)
    ids <- rownames(Z)
    # Function to count unordered pairs.
    choose2 <- function(x) x * (x - 1) / 2
    # Function to calculate the adjusted Rand index between two cluster assignments.
    ari <- function(a, b) {
      tab <- table(a, b)
      s <- sum(choose2(tab))
      sa <- sum(choose2(rowSums(tab)))
      sb <- sum(choose2(colSums(tab)))
      ex <- sa * sb / choose2(length(a))
      den <- (sa + sb) / 2 - ex
      if (abs(den) < 1e-12) {
        return(1)
      }
      (s - ex) / den
    }
    # Function to label existing clusters by mean Factor 1 score (ascending; ties use Factor 2).
    relabel <- function(g) {
      centers <- do.call(
        rbind,
        lapply(sort(unique(g)), function(j) colMeans(Z[g == j, , drop = FALSE]))
      )
      label_order <- if (ncol(centers) > 1L) {
        order(centers[, 1], centers[, 2])
      } else {
        order(centers[, 1])
      }
      match(g, label_order)
    }
    set.seed(20261015)
    subsets <- replicate(
      100,
      sort(sample.int(nrow(Z), floor(nrow(Z) * .8))),
      simplify = FALSE
    )
    D <- dist(Z)
    metrics <- list()
    assignments <- data.table(sample_id = ids)
    for (k in 2:5) {
      set.seed(20261015 + k)
      f <- cluster_samples(
        model,
        k = k,
        factors = fac,
        nstart = 100,
        iter.max = 200
      )
      stopifnot(f$ifault == 0)
      g <- relabel(f$cluster)
      if (k == 2L) assignments[, K2 := g]
      den <- num <- matrix(0L, nrow(Z), nrow(Z))
      jac <- matrix(NA_real_, 100, k)
      ar <- numeric(100)
      for (b in 1:100) {
        ii <- subsets[[b]]
        submodel <- model
        submodel@expectations$Z <- list(PAD = Z[ii, , drop = FALSE])
        submodel@dimensions$N <- c(PAD = length(ii))
        set.seed(20261015 + k * 1000 + b)
        ff <- cluster_samples(
          submodel,
          k = k,
          factors = fac,
          nstart = 50,
          iter.max = 200
        )
        stopifnot(ff$ifault == 0)
        den[ii, ii] <- den[ii, ii] + 1L
        num[ii, ii] <- num[ii, ii] + outer(ff$cluster, ff$cluster, '==')
        ar[b] <- ari(g[ii], ff$cluster)
        for (j in 1:k) {
          a <- g[ii] == j
          jac[b, j] <- max(vapply(
            1:k,
            function(h) {
              z <- ff$cluster == h
              sum(a & z) / sum(a | z)
            },
            numeric(1)
          ))
        }
      }
      C <- num / den
      s <- silhouette(g, D)[, 'sil_width']
      sizes <- as.integer(table(factor(g, levels = 1:k)))
      metrics[[k]] <- data.table(
        K = k,
        min_size = min(sizes),
        silhouette = mean(s),
        PAC = mean(C[upper.tri(C)] > .1 & C[upper.tri(C)] < .9),
        ARI_median = median(ar),
        min_Jaccard = min(colMeans(jac))
      )
      if (vip && k == 2L) consensus_k2 <- C
    }
    mt <- rbindlist(metrics)
    eligible <- mt[min_size >= 30 & min_Jaccard >= .75][order(
      PAC,
      -silhouette,
      K
    )]
    mt[, selected := K == eligible$K[1]]
    if (vip) {
      dir.create(file.path(base, "clustering"), recursive = TRUE, showWarnings = FALSE)

      saveRDS(
        list(Z = Z, Z_display = scale(Z), assignments = assignments,
             consensus = list("2" = consensus_k2)),
        file.path(base, "clustering/Clustering_Stability.rds")
      )
    }
    dir.create(
      file.path(base, "clustering"),
      recursive = TRUE,
      showWarnings = FALSE
    )
    fwrite(
      mt[, .(K, silhouette, PAC, ARI_median, min_Jaccard, selected)],
      file.path(base, "clustering/Cluster_Stability_Metrics.csv")
    )
    fwrite(
      assignments,
      file.path(base, "clustering/Patient_Assignments.csv")
    )
  }

  train_models(seed)
  if (is.null(seed)) {
    select_model()
    cluster_stability()
  }
}



# ==============================================================================
# RUN BOTH COHORTS
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (sys.nframe() == 0L && identical(args[1L], "--seed")) {
  fit_cohort(args[2L], as.integer(args[3L]))
} else {
  fit_cohort("VIP-PAD")
  fit_cohort("UKB")
}
