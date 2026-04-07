#!/usr/bin/env Rscript

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OMP_MAX_ACTIVE_LEVELS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  KMP_SHM_DISABLE = "1",
  KMP_SHM_DIR = "/tmp",
  TMPDIR = "/tmp"
)

suppressPackageStartupMessages({
  library(dplyr)
  library(PMA)
  library(pracma)
})

find_project_root <- function(start = getwd()) {
  current <- normalizePath(start, winslash = "/", mustWork = TRUE)

  repeat {
    if (file.exists(file.path(current, "CCAR3.Rproj")) || dir.exists(file.path(current, ".git"))) {
      return(current)
    }

    parent <- dirname(current)
    if (identical(parent, current)) {
      stop("Could not locate the CCAR3 project root from ", start, call. = FALSE)
    }
    current <- parent
  }
}

project_root <- find_project_root()
setwd(project_root)

source("~/Documents/ccar3/R/helpers.r")
source("~/Documents/ccar3/R/globals.R")
source("~/Documents/ccar3/R/utils.R")
source("~/Documents/ccar3/R/reduced_rank_regression.R")

source(file.path(project_root, "experiments", "evaluation.R"))
source(file.path(project_root, "experiments", "alternative_methods", "SAR.R"))
source(file.path(project_root, "experiments", "alternative_methods", "Parkhomenko.R"))
source(file.path(project_root, "experiments", "alternative_methods", "Witten_CrossValidation.R"))
source(file.path(project_root, "experiments", "alternative_methods", "Waaijenborg.R"))
source(file.path(project_root, "experiments", "alternative_methods", "scca_chao.R"))
source(file.path(project_root, "experiments", "alternative_methods", "GCA", "utils.R"))
source(file.path(project_root, "experiments", "alternative_methods", "GCA", "gca_to_cca.R"))
source(file.path(project_root, "experiments", "alternative_methods", "GCA", "init_process.R"))
source(file.path(project_root, "experiments", "alternative_methods", "GCA", "sgca_init.R"))
source(file.path(project_root, "experiments", "alternative_methods", "GCA", "sgca_tgd.R"))

params <- list(
  out_dir = file.path(
    project_root,
    "experiments",
    "real-data-experiments",
    "nutrimouse",
    "nutrimouse_benchmark_all_methods"
  ),
  rank = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_RANK", unset = "5")),
  repeats = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_REPEATS", unset = "5")),
  outer_folds = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_OUTER_FOLDS", unset = "10")),
  inner_folds = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_INNER_FOLDS", unset = "8")),
  full_model_cv_folds = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_FULL_MODEL_CV_FOLDS", unset = "8")),
  seed = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_SEED", unset = "123")),
  method_timeout_sec = as.numeric(Sys.getenv("CCAR3_NUTRIMOUSE_TIMEOUT_SEC", unset = "900")),
  rrr_lambdas = 10^seq(-3, 3, length.out = 30),
  witten_penalties = 10^seq(-3, 0, length.out = 30),
  sar_lambdas = c(1, 0.7, 0.5, 0.3, 0.1, 0.05, 0.01),
  parkhomenko_lambdas = c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2),
  waaijenborg_lambdas = c(1, 0.7, 0.5, 0.3, 0.1, 0.05, 0.01),
  pma_nperms = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_PMA_NPERMS", unset = "25")),
  fantope_rho_multipliers = c(0.5, 1, 2, 4),
  gca_lambda_grid = c(0.0025, 0.005, 0.01, 0.02),
  gca_k_grid = c(10L, 20L, 30L),
  gca_eta = as.numeric(Sys.getenv("CCAR3_NUTRIMOUSE_GCA_ETA", unset = "0.00025")),
  gca_maxiter = as.integer(Sys.getenv("CCAR3_NUTRIMOUSE_GCA_MAXITER", unset = "1500")),
  gca_convergence = 1e-6,
  gca_init_maxiter = 1000,
  chao_rho = 1,
  chao_lambda_max = 0.1,
  chao_num_lambda = 6,
  chao_niter = 500,
  chao_thresh = 0.01,
  #methods = c("cca_rrr", "witten_cv", "witten_pma", "sar", "parkhomenko", "waaijenborg")
  methods = c("cca_rrr", "witten_cv", "witten_pma", "sar", "parkhomenko", "waaijenborg", "scca_gao", "gca", "fantope")
)

scale_matrix <- function(mat) {
  mu <- colMeans(mat)
  sdv <- apply(mat, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1
  # Center and scale each feature column.
  out <- sweep(mat, 2, mu, "-")
  out <- sweep(out, 2, sdv, "/")
  out[!is.finite(out)] <- 0
  out
}

scale_train_test <- function(train, test) {
  mu <- colMeans(train)
  sdv <- apply(train, 2, sd)
  sdv[!is.finite(sdv) | sdv == 0] <- 1

  train_s <- sweep(train, 2, mu, "-")
  train_s <- sweep(train_s, 2, sdv, "/")
  test_s <- sweep(test, 2, mu, "-")
  test_s <- sweep(test_s, 2, sdv, "/")

  train_s[!is.finite(train_s)] <- 0
  test_s[!is.finite(test_s)] <- 0

  list(train = train_s, test = test_s)
}

make_folds <- function(indices, k, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  shuffled <- sample(indices, length(indices), replace = FALSE)
  split(shuffled, cut(seq_along(shuffled), breaks = k, labels = FALSE))
}

safe_mean <- function(x) {
  out <- mean(x, na.rm = TRUE)
  if (is.nan(out)) {
    return(NA_real_)
  }
  out
}

safe_sd <- function(x) {
  out <- stats::sd(x, na.rm = TRUE)
  if (is.nan(out)) {
    return(NA_real_)
  }
  out
}

cache_key <- function(...) {
  values <- list(...)
  pieces <- vapply(values, function(value) {
    if (length(value) == 1 && is.numeric(value)) {
      return(sprintf("%.12g", value))
    }
    paste(as.character(value), collapse = "-")
  }, character(1))
  paste(pieces, collapse = "__")
}

get_cached_value <- function(cache_env, key, builder) {
  if (!exists(key, envir = cache_env, inherits = FALSE)) {
    assign(key, builder(), envir = cache_env)
  }
  get(key, envir = cache_env, inherits = FALSE)
}

extract_param_value <- function(optimal_row, name) {
  if (is.null(optimal_row) || !(name %in% names(optimal_row))) {
    return(NA_real_)
  }

  value <- optimal_row[[name]][[1]]
  if (is.null(value) || length(value) == 0) {
    return(NA_real_)
  }

  as.numeric(value)
}

last_non_na <- function(x) {
  values <- as.numeric(x)
  values <- values[is.finite(values)]
  if (length(values) == 0) {
    return(NA_real_)
  }
  tail(values, 1)
}

extract_parkhomenko_uv <- function(fit) {
  if (!is.null(fit$uhat) && !is.null(fit$vhat)) {
    return(list(U = as.matrix(fit$uhat), V = as.matrix(fit$vhat)))
  }

  if (!is.null(fit$a) && !is.null(fit$b)) {
    return(list(U = matrix(fit$a, ncol = 1), V = matrix(fit$b, ncol = 1)))
  }

  stop("Unexpected Parkhomenko output structure.", call. = FALSE)
}

diag_cor_values <- function(A, B) {
  cor_mat <- suppressWarnings(stats::cor(A, B))
  if (is.null(dim(cor_mat))) {
    return(abs(as.numeric(cor_mat)))
  }
  diag(cor_mat)
}

calc_metrics <- function(U, V, X, Y, prefix) {
  x_scores <- X %*% U
  y_scores <- Y %*% V
  component_cor <- diag_cor_values(x_scores, y_scores)
  out <- list(
    mean_mse = mean((x_scores - y_scores)^2),
    mean_cor = safe_mean(component_cor),
    distance = subdistance(x_scores, y_scores),
    component_cor = component_cor,
    x_scores = x_scores,
    y_scores = y_scores
  )
  names(out)[1:3] <- paste0(prefix, c("_mse", "_cor", "_dist"))
  out
}

normalize_uv_by_cov <- function(U, V, X_train, Y_train) {
  Sx <- crossprod(X_train) / nrow(X_train)
  Sy <- crossprod(Y_train) / nrow(Y_train)
  U2 <- U
  V2 <- V

  if (!is.null(U2) && ncol(U2) > 0) {
    Mx <- t(U2) %*% Sx %*% U2
    inv_sqrt_x <- tryCatch(pracma::sqrtm(Mx)$Binv, error = function(e) NULL)
    if (!is.null(inv_sqrt_x)) {
      U2 <- U2 %*% inv_sqrt_x
    }
  }

  if (!is.null(V2) && ncol(V2) > 0) {
    My <- t(V2) %*% Sy %*% V2
    inv_sqrt_y <- tryCatch(pracma::sqrtm(My)$Binv, error = function(e) NULL)
    if (!is.null(inv_sqrt_y)) {
      V2 <- V2 %*% inv_sqrt_y
    }
  }

  list(U = U2, V = V2)
}

cca_normalization_error <- function(W, Sigma) {
  if (is.null(W) || ncol(W) == 0) {
    return(NA_real_)
  }
  gram <- t(W) %*% Sigma %*% W
  max(abs(gram - diag(ncol(W))), na.rm = TRUE)
}

normalize_and_check_cca <- function(U, V, X, Y, method_name, tol = 1e-3) {
  normalized <- normalize_uv_by_cov(as.matrix(U), as.matrix(V), X, Y)
  Sx <- crossprod(X) / nrow(X)
  Sy <- crossprod(Y) / nrow(Y)
  x_err <- cca_normalization_error(normalized$U, Sx)
  y_err <- cca_normalization_error(normalized$V, Sy)

  if ((is.finite(x_err) && x_err > tol) || (is.finite(y_err) && y_err > tol)) {
    stop(
      sprintf(
        "CCA normalization failed for %s (x_err=%.4g, y_err=%.4g).",
        method_name, x_err, y_err
      ),
      call. = FALSE
    )
  }

  list(
    U = normalized$U,
    V = normalized$V,
    x_error = x_err,
    y_error = y_err
  )
}

all_permutations <- function(values) {
  if (length(values) <= 1) {
    return(list(values))
  }

  out <- list()
  idx <- 0L
  for (i in seq_along(values)) {
    remainder <- values[-i]
    tails <- all_permutations(remainder)
    for (tail in tails) {
      idx <- idx + 1L
      out[[idx]] <- c(values[[i]], tail)
    }
  }
  out
}

best_cluster_accuracy <- function(predicted, truth) {
  predicted <- as.character(predicted)
  truth <- as.character(truth)
  pred_levels <- sort(unique(predicted))
  truth_levels <- sort(unique(truth))

  contingency <- table(
    factor(predicted, levels = pred_levels),
    factor(truth, levels = truth_levels)
  )

  if (nrow(contingency) != ncol(contingency)) {
    stop("Predicted clusters and truth labels must have the same number of classes.", call. = FALSE)
  }

  perms <- all_permutations(seq_along(truth_levels))
  scores <- vapply(perms, function(perm) {
    sum(contingency[cbind(seq_len(nrow(contingency)), perm)])
  }, numeric(1))
  best_perm <- perms[[which.max(scores)]]
  label_map <- setNames(truth_levels[best_perm], pred_levels)
  mapped_pred <- factor(label_map[predicted], levels = truth_levels)
  truth_factor <- factor(truth, levels = truth_levels)
  confusion <- table(predicted = mapped_pred, truth = truth_factor)

  list(
    accuracy = sum(diag(confusion)) / sum(confusion),
    confusion = confusion,
    mapping = label_map
  )
}

stratified_train_indices <- function(labels, prop = 0.75, seed = 107) {
  set.seed(seed)
  groups <- split(seq_along(labels), labels)
  train_idx <- lapply(groups, function(idx) {
    n_train <- floor(length(idx) * prop)
    n_train <- max(1L, min(length(idx) - 1L, n_train))
    sample(idx, n_train)
  })
  sort(unlist(train_idx, use.names = FALSE))
}

compute_full_fit_embedding <- function(U, X, n_components = min(params$rank, ncol(U))) {
  embedding <- 1 / sqrt(nrow(X)) * X %*% U[, seq_len(n_components), drop = FALSE]
  colnames(embedding) <- paste0("comp_", seq_len(n_components))
  embedding
}

run_lda_task <- function(embedding, labels, train_idx) {
  labels <- factor(labels)
  test_idx <- setdiff(seq_len(nrow(embedding)), train_idx)
  scaled <- scale_train_test(embedding[train_idx, , drop = FALSE], embedding[test_idx, , drop = FALSE])
  train_df <- data.frame(label = labels[train_idx], scaled$train, check.names = FALSE)
  test_df <- data.frame(scaled$test, check.names = FALSE)

  fit <- MASS::lda(label ~ ., data = train_df)
  predictions <- predict(fit, newdata = test_df)$class
  truth <- labels[test_idx]
  confusion <- table(predicted = predictions, truth = truth)

  list(
    accuracy = mean(predictions == truth),
    confusion = confusion,
    predictions = predictions,
    truth = truth
  )
}

run_mclust_task <- function(embedding, labels, G) {
  if (!requireNamespace("mclust", quietly = TRUE)) {
    stop("Package 'mclust' is required for the clustering comparisons.", call. = FALSE)
  }

  if (!("package:mclust" %in% search())) {
    suppressPackageStartupMessages(library(mclust))
  }

  model <- Mclust(embedding, G = G:G, verbose = FALSE)
  match <- best_cluster_accuracy(model$classification, labels)

  list(
    accuracy = match$accuracy,
    ari = mclust::adjustedRandIndex(as.numeric(model$classification), as.numeric(factor(labels))),
    confusion = match$confusion,
    mapping = match$mapping,
    classification = model$classification,
    model = model
  )
}

evaluate_embedding_cluster_task <- function(embedding, labels, label_name) {
  cluster_task <- run_mclust_task(embedding, labels, length(unique(labels)))
  summary <- data.frame(
    cluster_accuracy = cluster_task$accuracy,
    cluster_ari = cluster_task$ari,
    stringsAsFactors = FALSE
  )
  names(summary) <- paste0(label_name, c("_cluster_accuracy", "_cluster_ari"))

  list(
    summary = summary,
    artifact = list(
      embedding = embedding,
      cluster = cluster_task
    )
  )
}

evaluate_embedding_label_task <- function(embedding, labels, label_name, train_prop = 0.75, seed = 107) {
  train_idx <- stratified_train_indices(labels, prop = train_prop, seed = seed)
  lda_task <- run_lda_task(embedding, labels, train_idx)
  cluster_task <- evaluate_embedding_cluster_task(embedding, labels, label_name)

  summary <- data.frame(
    lda_accuracy = lda_task$accuracy,
    cluster_accuracy = cluster_task$summary[[paste0(label_name, "_cluster_accuracy")]],
    cluster_ari = cluster_task$summary[[paste0(label_name, "_cluster_ari")]],
    stringsAsFactors = FALSE
  )
  names(summary) <- paste0(label_name, c("_lda_accuracy", "_cluster_accuracy", "_cluster_ari"))

  list(
    summary = summary,
    artifact = list(
      embedding = embedding,
      train_idx = train_idx,
      lda = lda_task,
      cluster = cluster_task$artifact$cluster
    )
  )
}

evaluate_diet_prediction_task <- function(U, X, metadata, train_prop = 0.75, seed = 107, embedding = NULL) {
  if (is.null(embedding)) {
    embedding <- compute_full_fit_embedding(U, X)
  }
  evaluate_embedding_label_task(embedding, metadata$diet, "diet", train_prop = train_prop, seed = seed)
}

evaluate_genotype_prediction_task <- function(U, X, metadata, train_prop = 0.75, seed = 108, embedding = NULL) {
  if (is.null(embedding)) {
    embedding <- compute_full_fit_embedding(U, X)
  }
  evaluate_embedding_label_task(embedding, metadata$genotype, "genotype", train_prop = train_prop, seed = seed)
}

evaluate_full_fit_tasks <- function(U, X, metadata) {
  embedding <- compute_full_fit_embedding(U, X)
  diet_task <- evaluate_diet_prediction_task(U, X, metadata, seed = 107, embedding = embedding)
  genotype_task <- evaluate_genotype_prediction_task(U, X, metadata, seed = 108, embedding = embedding)

  list(
    summary = bind_cols(diet_task$summary, genotype_task$summary),
    artifact = list(
      embedding = embedding,
      diet_train_idx = diet_task$artifact$train_idx,
      genotype_train_idx = genotype_task$artifact$train_idx,
      diet_lda = diet_task$artifact$lda,
      genotype_lda = genotype_task$artifact$lda,
      diet_cluster = diet_task$artifact$cluster,
      genotype_cluster = genotype_task$artifact$cluster
    )
  )
}

run_method_with_timeout <- function(fun, timeout_sec) {
  timed_out <- FALSE
  error_message <- NA_character_

  result <- tryCatch({
    if (requireNamespace("R.utils", quietly = TRUE)) {
      R.utils::withTimeout(fun(), timeout = timeout_sec, onTimeout = "error")
    } else {
      setTimeLimit(elapsed = timeout_sec, transient = TRUE)
      on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE), add = TRUE)
      fun()
    }
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (inherits(e, "TimeoutException") || grepl("reached elapsed time limit", msg, fixed = TRUE)) {
      timed_out <<- TRUE
    }
    error_message <<- msg
    NULL
  })

  list(result = result, timed_out = timed_out, error_message = error_message)
}

save_benchmark_artifact <- function(out_dir, repeat_id, fold_id, method, artifact) {
  details_dir <- file.path(out_dir, "benchmark_details")
  dir.create(details_dir, recursive = TRUE, showWarnings = FALSE)
  safe_method <- gsub("[^A-Za-z0-9_]+", "_", method)
  path <- file.path(
    details_dir,
    sprintf("repeat_%02d_fold_%02d_%s.rds", repeat_id, fold_id, safe_method)
  )
  saveRDS(artifact, path)
  path
}

write_benchmark_outputs <- function(results_df, out_dir) {
  if (nrow(results_df) == 0) {
    return(invisible(NULL))
  }

  utils::write.csv(
    results_df,
    file.path(out_dir, "benchmark_runs.csv"),
    row.names = FALSE
  )

  summary_df <- results_df %>%
    group_by(method) %>%
    summarise(
      mean_train_mse = safe_mean(train_mse),
      sd_train_mse = safe_sd(train_mse),
      mean_test_mse = safe_mean(test_mse),
      sd_test_mse = safe_sd(test_mse),
      mean_train_cor = safe_mean(train_cor),
      mean_test_cor = safe_mean(test_cor),
      mean_test_dist = safe_mean(test_dist),
      mean_test_genotype_cluster_accuracy = safe_mean(test_genotype_cluster_accuracy),
      mean_test_genotype_cluster_ari = safe_mean(test_genotype_cluster_ari),
      mean_runtime_sec = safe_mean(runtime_sec),
      n_timeouts = sum(timed_out),
      n_completed = sum(!is.na(test_mse)),
      .groups = "drop"
    ) %>%
    arrange(mean_test_mse)

  utils::write.csv(
    summary_df,
    file.path(out_dir, "benchmark_summary.csv"),
    row.names = FALSE
  )

  invisible(summary_df)
}

load_nutrimouse_data <- function() {
  if (!requireNamespace("mixOmics", quietly = TRUE)) {
    stop("Package 'mixOmics' is required for the nutrimouse dataset.", call. = FALSE)
  }

  env <- new.env(parent = emptyenv())
  data("nutrimouse", package = "mixOmics", envir = env)
  nutrimouse <- get("nutrimouse", envir = env)

  X_raw <- as.matrix(nutrimouse$gene)
  Y_raw <- as.matrix(nutrimouse$lipid)

  metadata <- data.frame(
    sample_id = seq_len(nrow(X_raw)),
    genotype = nutrimouse$genotype,
    diet = nutrimouse$diet,
    stringsAsFactors = FALSE
  )

  list(
    X_raw = X_raw,
    Y_raw = Y_raw,
    metadata = metadata
  )
}

prepare_full_nutrimouse_data <- function(data_obj) {
  X_raw <- scale(data_obj$X_raw)
  Y_raw <- scale(data_obj$Y_raw)

  list(
    X_raw = X_raw,
    Y_raw = Y_raw,
    X_full = scale_matrix(X_raw),
    Y_full = scale_matrix(Y_raw),
    metadata = data_obj$metadata
  )
}

fantope_default_rho <- function(X, Y) {
  0.5 * sqrt(log(ncol(X) + ncol(Y)) / nrow(X))
}

rho_grid_from_default <- function(X, Y, multipliers = params$fantope_rho_multipliers) {
  sort(unique(as.numeric(fantope_default_rho(X, Y) * multipliers)))
}

fit_witten_with_penalties <- function(X, Y, penalty_x, penalty_y) {
  fit <- PMA::CCA(
    x = X,
    z = Y,
    typex = "standard",
    typez = "standard",
    K = params$rank,
    penaltyx = penalty_x,
    penaltyz = penalty_y,
    trace = FALSE,
    standardize = FALSE
  )
  normalized <- normalize_and_check_cca(fit$u, fit$v, X, Y, "Witten")
  list(
    U = normalized$U,
    V = normalized$V,
    fit = fit,
    normalized = normalized
  )
}

fit_cca_rrr_with_lambda <- function(X, Y, lambda) {
  fit <- cca_rrr(
    X = X,
    Y = Y,
    Sx = NULL,
    Sy = NULL,
    lambda = lambda,
    r = params$rank,
    highdim = TRUE,
    rho = 1,
    niter = 2 * 1e4,
    thresh = 1e-6,
    solver = "ADMM",
    thresh_0 = 0,
    standardize = FALSE,
    preprocess = FALSE,
    LW_Sy = TRUE,
    verbose = FALSE
  )
  list(U = fit$U, V = fit$V, fit = fit)
}

build_mask <- function(p1, p2) {
  p <- p1 + p2
  mask <- matrix(0, p, p)
  mask[seq_len(p1), seq_len(p1)] <- 1
  idx_y <- (p1 + 1):p
  mask[idx_y, idx_y] <- 1
  mask
}

build_gca_init <- function(X, Y, rho) {
  p1 <- ncol(X)
  p2 <- ncol(Y)
  S <- stats::cov(cbind(X, Y))
  mask <- build_mask(p1, p2)
  sigma0hat <- S * mask

  init_solver <- sgca_init(
    A = S,
    B = sigma0hat,
    rho = rho,
    K = params$rank,
    maxiter = params$gca_init_maxiter,
    trace = FALSE
  )

  list(
    S = S,
    sigma0hat = sigma0hat,
    p1 = p1,
    p2 = p2,
    ainit = init_process(init_solver$Pi, params$rank),
    init_solver = init_solver,
    rho = rho
  )
}

fit_fantope_with_rho <- function(X, Y, rho) {
  fit <- Fantope(X, Y, r = params$rank, rho = rho)
  normalized <- normalize_and_check_cca(fit$u, fit$v, X, Y, "Fantope")
  list(U = normalized$U, V = normalized$V, fit = fit, normalized = normalized, rho = rho)
}

cv_fantope_with_info <- function(X, Y, kfolds = params$inner_folds) {
  folds <- make_folds(seq_len(nrow(X)), k = kfolds)
  rho_grid <- rho_grid_from_default(X, Y)

  if (length(rho_grid) == 1) {
    final <- fit_fantope_with_rho(X, Y, rho_grid[[1]])
    return(list(
      U = final$U,
      V = final$V,
      rho = rho_grid[[1]],
      artifact = list(
        fit = final$fit,
        cv_summary = data.frame(rho = rho_grid[[1]], mse = NA_real_, se = NA_real_)
      )
    ))
  }

  cv_summary <- lapply(rho_grid, function(rho) {
    fold_mse <- vapply(folds, function(val_idx) {
      fit <- tryCatch(
        fit_fantope_with_rho(X[-val_idx, , drop = FALSE], Y[-val_idx, , drop = FALSE], rho),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NA_real_)
      }
      mean((X[val_idx, , drop = FALSE] %*% fit$U - Y[val_idx, , drop = FALSE] %*% fit$V)^2)
    }, numeric(1))

    data.frame(
      rho = rho,
      mse = safe_mean(fold_mse),
      se = safe_sd(fold_mse) / sqrt(sum(!is.na(fold_mse))),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  best_rho <- cv_summary$rho[which.min(cv_summary$mse)]
  final <- fit_fantope_with_rho(X, Y, best_rho)

  list(
    U = final$U,
    V = final$V,
    rho = best_rho,
    artifact = list(fit = final$fit, cv_summary = cv_summary)
  )
}

fit_gca_with_params <- function(X, Y, lambda, k, rho = fantope_default_rho(X, Y), init_fit = NULL) {
  if (is.null(init_fit)) {
    init_fit <- build_gca_init(X, Y, rho)
  }

  a_est <- sgca_tgd(
    A = init_fit$S,
    B = init_fit$sigma0hat,
    r = params$rank,
    init = init_fit$ainit,
    k = k,
    lambda = lambda,
    eta = params$gca_eta,
    convergence = params$gca_convergence,
    maxiter = params$gca_maxiter,
    plot = FALSE
  )
  converted <- gca_to_cca(a_est, init_fit$S, c(init_fit$p1, init_fit$p2))
  normalized <- normalize_and_check_cca(converted$u, converted$v, X, Y, "GCA")

  list(
    U = normalized$U,
    V = normalized$V,
    fit = a_est,
    init = init_fit$ainit,
    init_fit = init_fit$init_solver,
    normalized = normalized,
    rho = rho
  )
}

cv_gca_with_info <- function(X, Y, kfolds = params$inner_folds) {
  folds <- make_folds(seq_len(nrow(X)), k = kfolds)
  grid <- expand.grid(
    rho = rho_grid_from_default(X, Y),
    lambda = params$gca_lambda_grid,
    k = params$gca_k_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(grid) == 1) {
    final <- fit_gca_with_params(
      X,
      Y,
      lambda = grid$lambda[[1]],
      k = grid$k[[1]],
      rho = grid$rho[[1]]
    )
    return(list(
      U = final$U,
      V = final$V,
      lambda_x = grid$lambda[[1]],
      k = grid$k[[1]],
      rho = grid$rho[[1]],
      artifact = list(
        fit = final$fit,
        init = final$init,
        init_fit = final$init_fit,
        cv_summary = data.frame(
          rho = grid$rho[[1]],
          lambda_x = grid$lambda[[1]],
          k = grid$k[[1]],
          mse = NA_real_,
          se = NA_real_
        )
      )
    ))
  }

  init_cache <- new.env(parent = emptyenv())
  get_fold_init <- function(fold_id, rho) {
    key <- cache_key("gca", fold_id, rho)
    get_cached_value(init_cache, key, function() {
      val_idx <- folds[[fold_id]]
      tryCatch(
        build_gca_init(
          X[-val_idx, , drop = FALSE],
          Y[-val_idx, , drop = FALSE],
          rho = rho
        ),
        error = function(e) NULL
      )
    })
  }

  cv_summary <- lapply(seq_len(nrow(grid)), function(i) {
    rho <- grid$rho[i]
    lambda <- grid$lambda[i]
    k <- grid$k[i]
    fold_mse <- vapply(seq_along(folds), function(fold_id) {
      val_idx <- folds[[fold_id]]
      init_fit <- get_fold_init(fold_id, rho)
      if (is.null(init_fit)) {
        return(NA_real_)
      }
      fit <- tryCatch(
        fit_gca_with_params(
          X[-val_idx, , drop = FALSE],
          Y[-val_idx, , drop = FALSE],
          lambda = lambda,
          k = k,
          rho = rho,
          init_fit = init_fit
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NA_real_)
      }
      mean((X[val_idx, , drop = FALSE] %*% fit$U - Y[val_idx, , drop = FALSE] %*% fit$V)^2)
    }, numeric(1))

    data.frame(
      rho = rho,
      lambda_x = lambda,
      k = k,
      mse = safe_mean(fold_mse),
      se = safe_sd(fold_mse) / sqrt(sum(!is.na(fold_mse))),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  best_row <- cv_summary %>% arrange(mse, se) %>% slice(1)
  full_init <- build_gca_init(X, Y, best_row$rho)
  final <- fit_gca_with_params(
    X,
    Y,
    lambda = best_row$lambda_x,
    k = best_row$k,
    rho = best_row$rho,
    init_fit = full_init
  )

  list(
    U = final$U,
    V = final$V,
    lambda_x = best_row$lambda_x,
    k = best_row$k,
    rho = best_row$rho,
    artifact = list(
      fit = final$fit,
      init = final$init,
      init_fit = final$init_fit,
      cv_summary = cv_summary
    )
  )
}

group_lasso_lambda_grid <- function() {
  params$chao_rho * seq(0, 1, length.out = params$chao_num_lambda + 1) * params$chao_lambda_max
}

cv_group_lasso_with_info <- function(X, Y_tilde, lambda_values, kfolds = params$inner_folds) {
  folds <- make_folds(seq_len(nrow(X)), k = kfolds)

  cv_summary <- lapply(lambda_values, function(lambda) {
    fold_mse <- vapply(folds, function(val_idx) {
      train_idx <- setdiff(seq_len(nrow(X)), val_idx)
      fit <- tryCatch(
        group_lasso(
          X[train_idx, , drop = FALSE],
          Y_tilde[train_idx, , drop = FALSE],
          lambda = lambda,
          rho = params$chao_rho,
          niter = params$chao_niter,
          thresh = params$chao_thresh
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NA_real_)
      }
      mean((X[val_idx, , drop = FALSE] %*% fit - Y_tilde[val_idx, , drop = FALSE])^2)
    }, numeric(1))

    data.frame(
      lambda = lambda,
      mse = safe_mean(fold_mse),
      se = safe_sd(fold_mse) / sqrt(sum(!is.na(fold_mse))),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  best_lambda <- cv_summary$lambda[which.min(cv_summary$mse)]
  final <- group_lasso(
    X,
    Y_tilde,
    lambda = best_lambda,
    rho = params$chao_rho,
    niter = params$chao_niter,
    thresh = params$chao_thresh
  )

  list(B = final, lambda = best_lambda, cv_summary = cv_summary)
}

fit_scca_gao_with_lambdas <- function(
  X,
  Y,
  lambda_x,
  lambda_y,
  rho = fantope_default_rho(X, Y),
  init_fit = NULL
) {
  if (is.null(init_fit)) {
    init_fit <- fit_fantope_with_rho(X, Y, rho)
  }

  u_raw <- group_lasso(
    X,
    Y %*% init_fit$V,
    lambda = lambda_x,
    rho = params$chao_rho,
    niter = params$chao_niter,
    thresh = params$chao_thresh
  )
  v_raw <- group_lasso(
    Y,
    X %*% init_fit$U,
    lambda = lambda_y,
    rho = params$chao_rho,
    niter = params$chao_niter,
    thresh = params$chao_thresh
  )

  normalized <- normalize_and_check_cca(u_raw, v_raw, X, Y, "scca_gao")

  list(
    U = normalized$U,
    V = normalized$V,
    rho = rho,
    init_fit = init_fit,
    raw_u = u_raw,
    raw_v = v_raw
  )
}

fit_scca_gao_with_info <- function(X, Y, kfolds = params$inner_folds) {
  folds <- make_folds(seq_len(nrow(X)), k = kfolds)
  rho_grid <- rho_grid_from_default(X, Y)
  lambda_grid <- group_lasso_lambda_grid()
  grid <- expand.grid(
    rho = rho_grid,
    lambda_x = lambda_grid,
    lambda_y = lambda_grid,
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  if (nrow(grid) == 1) {
    init_fit <- fit_fantope_with_rho(X, Y, grid$rho[[1]])
    final <- fit_scca_gao_with_lambdas(
      X,
      Y,
      lambda_x = grid$lambda_x[[1]],
      lambda_y = grid$lambda_y[[1]],
      rho = grid$rho[[1]],
      init_fit = init_fit
    )
    return(list(
      U = final$U,
      V = final$V,
      lambda_x = grid$lambda_x[[1]],
      lambda_y = grid$lambda_y[[1]],
      rho = grid$rho[[1]],
      artifact = list(
        init_fit = init_fit$fit,
        cv_summary = data.frame(
          rho = grid$rho[[1]],
          lambda_x = grid$lambda_x[[1]],
          lambda_y = grid$lambda_y[[1]],
          mse = NA_real_,
          se = NA_real_
        ),
        raw_u = final$raw_u,
        raw_v = final$raw_v
      )
    ))
  }

  init_cache <- new.env(parent = emptyenv())
  get_fold_init <- function(fold_id, rho) {
    key <- cache_key("scca_gao", fold_id, rho)
    get_cached_value(init_cache, key, function() {
      val_idx <- folds[[fold_id]]
      tryCatch(
        fit_fantope_with_rho(
          X[-val_idx, , drop = FALSE],
          Y[-val_idx, , drop = FALSE],
          rho = rho
        ),
        error = function(e) NULL
      )
    })
  }

  cv_summary <- lapply(seq_len(nrow(grid)), function(i) {
    rho <- grid$rho[i]
    lambda_x <- grid$lambda_x[i]
    lambda_y <- grid$lambda_y[i]

    fold_mse <- vapply(seq_along(folds), function(fold_id) {
      val_idx <- folds[[fold_id]]
      init_fit <- get_fold_init(fold_id, rho)
      if (is.null(init_fit)) {
        return(NA_real_)
      }
      fit <- tryCatch(
        fit_scca_gao_with_lambdas(
          X[-val_idx, , drop = FALSE],
          Y[-val_idx, , drop = FALSE],
          lambda_x = lambda_x,
          lambda_y = lambda_y,
          rho = rho,
          init_fit = init_fit
        ),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NA_real_)
      }
      mean((X[val_idx, , drop = FALSE] %*% fit$U - Y[val_idx, , drop = FALSE] %*% fit$V)^2)
    }, numeric(1))

    data.frame(
      rho = rho,
      lambda_x = lambda_x,
      lambda_y = lambda_y,
      mse = safe_mean(fold_mse),
      se = safe_sd(fold_mse) / sqrt(sum(!is.na(fold_mse))),
      stringsAsFactors = FALSE
    )
  }) %>%
    bind_rows()

  best_row <- cv_summary %>% arrange(mse, se) %>% slice(1)
  init_fit <- fit_fantope_with_rho(X, Y, best_row$rho)
  final <- fit_scca_gao_with_lambdas(
    X,
    Y,
    lambda_x = best_row$lambda_x,
    lambda_y = best_row$lambda_y,
    rho = best_row$rho,
    init_fit = init_fit
  )

  list(
    U = final$U,
    V = final$V,
    lambda_x = best_row$lambda_x,
    lambda_y = best_row$lambda_y,
    rho = best_row$rho,
    artifact = list(
      init_fit = init_fit$fit,
      cv_summary = cv_summary,
      raw_u = final$raw_u,
      raw_v = final$raw_v
    )
  )
}

fit_witten_cv_with_info <- function(X, Y, kfolds = params$inner_folds) {
  cv <- Witten.CV(
    X = X,
    Y = Y,
    n.cv = kfolds,
    rank = params$rank,
    lambdax = params$witten_penalties,
    lambday = c(0),
    standardize = FALSE
  )
  fit <- fit_witten_with_penalties(X, Y, cv$lambdax.opt, cv$lambday.opt)

  list(
    U = fit$U,
    V = fit$V,
    lambda_x = cv$lambdax.opt,
    lambda_y = cv$lambday.opt,
    rho = NA_real_,
    k = NA_real_,
    artifact = list(cv = cv, fit = fit$fit, normalized = fit$normalized)
  )
}

fit_cca_rrr_cv_with_info <- function(X, Y, kfolds = params$inner_folds) {
  fit <- cca_rrr_cv(
    X = X,
    Y = Y,
    r = params$rank,
    lambdas = params$rrr_lambdas,
    kfolds = kfolds,
    solver = "ADMM",
    parallelize = FALSE,
    LW_Sy = TRUE,
    standardize = FALSE,
    rho = 1,
    thresh_0 = 0,
    niter = 2 * 1e4,
    thresh = 1e-6,
    verbose = FALSE
  )

  list(
    U = fit$U,
    V = fit$V,
    lambda_x = fit$lambda_x,
    lambda_y = NA_real_,
    rho = NA_real_,
    k = NA_real_,
    artifact = fit
  )
}

fit_sar_with_info <- function(X, Y, kfolds = params$inner_folds) {
  fit <- SparseCCA(
    X = X,
    Y = Y,
    lambdaAseq = params$rrr_lambdas,
    lambdaBseq = c(0),
    rank = params$rank,
    selection.criterion = 2,
    n.cv = kfolds,
    max.iter = 20,
    conv = 1e-2
  )
  normalized <- normalize_and_check_cca(fit$uhat, fit$vhat, X, Y, "SAR")

  list(
    U = normalized$U,
    V = normalized$V,
    lambda_x = fit$lambdaA,
    lambda_y = fit$lambdaB,
    rho = NA_real_,
    k = NA_real_,
    artifact = list(fit = fit, normalized = normalized)
  )
}

fit_parkhomenko_with_info <- function(X, Y, kfolds = params$inner_folds) {
  fit <- SCCA_Parkhomenko(
    x.data = X,
    y.data = Y,
    n.cv = kfolds,
    lambda.v.seq = params$witten_penalties,
    lambda.u.seq = c(0),
    Krank = params$rank
  )
  uv <- extract_parkhomenko_uv(fit)
  normalized <- normalize_and_check_cca(uv$U, uv$V, X, Y, "Parkhomenko")

  list(
    U = normalized$U,
    V = normalized$V,
    lambda_x = fit$lambda.uopt,
    lambda_y = fit$lambda.vopt,
    rho = NA_real_,
    k = NA_real_,
    artifact = list(fit = fit, normalized = normalized)
  )
}

fit_waaijenborg_with_info <- function(X, Y, kfolds = params$inner_folds) {
  fit <- Waaijenborg(
    X = X,
    Y = Y,
    lambdaxseq = params$rrr_lambdas,
    lambdayseq = c(0),
    rank = params$rank,
    selection.criterion = 2,
    n.cv = kfolds,
    max.iter = 20,
    conv = 1e-3
  )
  normalized <- normalize_and_check_cca(fit$vhat, fit$uhat, X, Y, "Waaijenborg")

  list(
    U = normalized$U,
    V = normalized$V,
    lambda_x = last_non_na(fit$lambdax_FINAL),
    lambda_y = last_non_na(fit$lambday_FINAL),
    rho = NA_real_,
    k = NA_real_,
    artifact = list(fit = fit, normalized = normalized)
  )
}

benchmark_fitters <- list(
  cca_rrr = function(X, Y) {
    fit_cca_rrr_cv_with_info(X, Y, kfolds = params$inner_folds)
  },
  witten_cv = function(X, Y) {
    fit_witten_cv_with_info(X, Y, kfolds = params$inner_folds)
  },
  witten_pma = function(X, Y) {
    perm <- PMA::CCA.permute(
      x = X,
      z = Y,
      typex = "standard",
      typez = "standard",
      penaltyxs = params$witten_penalties,
      penaltyzs = c(0),
      standardize = FALSE,
      nperms = params$pma_nperms,
      trace = FALSE
    )
    fit <- fit_witten_with_penalties(X, Y, perm$bestpenaltyx, perm$bestpenaltyz)

    list(
      U = fit$U,
      V = fit$V,
      lambda_x = perm$bestpenaltyx,
      lambda_y = perm$bestpenaltyz,
      rho = NA_real_,
      k = NA_real_,
      artifact = list(perm = perm, fit = fit$fit, normalized = fit$normalized)
    )
  },
  sar = function(X, Y) {
    fit_sar_with_info(X, Y, kfolds = params$inner_folds)
  },
  parkhomenko = function(X, Y) {
    fit_parkhomenko_with_info(X, Y, kfolds = params$inner_folds)
  },
  waaijenborg = function(X, Y) {
    fit_waaijenborg_with_info(X, Y, kfolds = params$inner_folds)
  },
  scca_gao = function(X, Y) {
    fit_scca_gao_with_info(X, Y, kfolds = params$inner_folds)
  },
  gca = function(X, Y) {
    cv_gca_with_info(X, Y, kfolds = params$inner_folds)
  },
  fantope = function(X, Y) {
    cv_fantope_with_info(X, Y, kfolds = params$inner_folds)
  }
)

full_model_cv_fitters <- list(
  cca_rrr = function(X, Y) {
    fit_cca_rrr_cv_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  witten_cv = function(X, Y) {
    fit_witten_cv_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  witten_pma = function(X, Y) {
    fit_witten_cv_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  sar = function(X, Y) {
    fit_sar_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  parkhomenko = function(X, Y) {
    fit_parkhomenko_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  waaijenborg = function(X, Y) {
    fit_waaijenborg_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  scca_gao = function(X, Y) {
    fit_scca_gao_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  gca = function(X, Y) {
    cv_gca_with_info(X, Y, kfolds = params$full_model_cv_folds)
  },
  fantope = function(X, Y) {
    cv_fantope_with_info(X, Y, kfolds = params$full_model_cv_folds)
  }
)

choose_optimal_params <- function(results_df) {
  successful <- results_df %>%
    filter(!timed_out, is.finite(test_mse))

  if (nrow(successful) == 0) {
    return(data.frame())
  }

  method_names <- unique(successful$method)
  method <- method_names[[1]]
  param_cols <- c("lambda_x", "lambda_y", "rho", "k")
  param_cols <- param_cols[vapply(successful[param_cols], function(x) any(!is.na(x)), logical(1))]

  if (length(param_cols) == 0) {
    return(data.frame(
      method = method,
      mean_test_mse = safe_mean(successful$test_mse),
      n = nrow(successful),
      stringsAsFactors = FALSE
    ))
  }

  successful %>%
    group_by(across(all_of(param_cols))) %>%
    summarise(
      mean_test_mse = safe_mean(test_mse),
      sd_test_mse = safe_sd(test_mse),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(mean_test_mse) %>%
    slice(1) %>%
    mutate(method = method, .before = 1)
}

build_full_fit_row <- function(method, fit, X, Y, metadata, optimal_row) {
  metrics <- calc_metrics(fit$U, fit$V, X, Y, "full")
  tasks <- evaluate_full_fit_tasks(fit$U, X, metadata)
  component_cor <- metrics$component_cor
  component_names <- paste0("cor_", seq_along(component_cor))
  component_df <- as.data.frame(as.list(stats::setNames(component_cor, component_names)))

  row <- data.frame(
    method = method,
    lambda_x = extract_param_value(optimal_row, "lambda_x"),
    lambda_y = extract_param_value(optimal_row, "lambda_y"),
    rho = extract_param_value(optimal_row, "rho"),
    k = extract_param_value(optimal_row, "k"),
    cv_folds = extract_param_value(optimal_row, "cv_folds"),
    selection_runtime_sec = extract_param_value(optimal_row, "selection_runtime_sec"),
    full_mse = metrics$full_mse,
    full_cor = metrics$full_cor,
    full_dist = metrics$full_dist,
    stringsAsFactors = FALSE
  )

  list(
    summary_row = bind_cols(row, component_df, tasks$summary),
    artifact = list(
      method = method,
      params = optimal_row,
      U = fit$U,
      V = fit$V,
      x_scores = metrics$x_scores,
      y_scores = metrics$y_scores,
      component_cor = component_cor,
      metadata = metadata,
      cv_folds = extract_param_value(optimal_row, "cv_folds"),
      selection_runtime_sec = extract_param_value(optimal_row, "selection_runtime_sec"),
      tasks = tasks$artifact,
      fit = fit$artifact
    )
  )
}

run_full_model_diet_tasks <- function(full_fit_artifacts, X, metadata, out_dir = NULL) {
  diet_rows <- list()
  diet_artifacts <- list()

  for (method in names(full_fit_artifacts)) {
    fit_artifact <- full_fit_artifacts[[method]]
    if (is.null(fit_artifact$U)) {
      next
    }

    diet_task <- if (!is.null(fit_artifact$tasks) &&
      !is.null(fit_artifact$tasks$diet_lda) &&
      !is.null(fit_artifact$tasks$diet_cluster) &&
      !is.null(fit_artifact$tasks$embedding)) {
      list(
        summary = data.frame(
          diet_lda_accuracy = fit_artifact$tasks$diet_lda$accuracy,
          diet_cluster_accuracy = fit_artifact$tasks$diet_cluster$accuracy,
          diet_cluster_ari = fit_artifact$tasks$diet_cluster$ari,
          stringsAsFactors = FALSE
        ),
        artifact = list(
          embedding = fit_artifact$tasks$embedding,
          train_idx = fit_artifact$tasks$diet_train_idx,
          lda = fit_artifact$tasks$diet_lda,
          cluster = fit_artifact$tasks$diet_cluster
        )
      )
    } else {
      evaluate_diet_prediction_task(fit_artifact$U, X, metadata)
    }

    optimal_row <- fit_artifact$params
    diet_rows[[method]] <- data.frame(
      method = method,
      lambda_x = extract_param_value(optimal_row, "lambda_x"),
      lambda_y = extract_param_value(optimal_row, "lambda_y"),
      rho = extract_param_value(optimal_row, "rho"),
      k = extract_param_value(optimal_row, "k"),
      cv_folds = extract_param_value(optimal_row, "cv_folds"),
      selection_runtime_sec = extract_param_value(optimal_row, "selection_runtime_sec"),
      diet_lda_accuracy = diet_task$summary$diet_lda_accuracy,
      diet_cluster_accuracy = diet_task$summary$diet_cluster_accuracy,
      diet_cluster_ari = diet_task$summary$diet_cluster_ari,
      stringsAsFactors = FALSE
    )

    diet_artifacts[[method]] <- list(
      method = method,
      params = optimal_row,
      U = fit_artifact$U,
      diet = diet_task$artifact
    )
  }

  diet_summary <- bind_rows(diet_rows) %>%
    arrange(desc(diet_lda_accuracy), desc(diet_cluster_accuracy), desc(diet_cluster_ari))

  if (!is.null(out_dir)) {
    utils::write.csv(
      diet_summary,
      file.path(out_dir, "full_model_diet_summary.csv"),
      row.names = FALSE
    )
    saveRDS(
      list(summary = diet_summary, artifacts = diet_artifacts),
      file.path(out_dir, "full_model_diet_tasks.rds")
    )
  }

  list(summary = diet_summary, artifacts = diet_artifacts)
}

main <- function() {
  dir.create(params$out_dir, recursive = TRUE, showWarnings = FALSE)

  data_obj <- load_nutrimouse_data()
  full_data <- prepare_full_nutrimouse_data(data_obj)
  X_raw <- full_data$X_raw
  Y_raw <- full_data$Y_raw
  metadata <- full_data$metadata

  message(
    "Benchmarking nutrimouse with methods: ",
    paste(params$methods, collapse = ", ")
  )

  benchmark_rows <- list()
  row_idx <- 0L

  # Cross-validation benchmark loop:
  # 1. Re-split the nutrimouse samples into repeated outer folds.
  # 2. Tune each method only on the outer-training split, including rho for Fantope-based initializations.
  # 3. Score the fitted canonical variates on the held-out outer fold.
  # 4. Track held-out genotype clustering on the fitted XU embedding.
  # 5. Aggregate the selected parameters, then refit each method on all samples.
  # for (repeat_id in seq_len(params$repeats)) {
  #   set.seed(params$seed + repeat_id)
  #   folds <- make_folds(seq_len(nrow(X_raw)), k = params$outer_folds)
  # 
  #   for (fold_id in seq_along(folds)) {
  #     test_idx <- folds[[fold_id]]
  #     test_metadata <- metadata[test_idx, , drop = FALSE]
  # 
  #     #X_scaled <- scale_train_test(X_raw[-test_idx, , drop = FALSE], X_raw[test_idx, , drop = FALSE])
  #     #Y_scaled <- scale_train_test(Y_raw[-test_idx, , drop = FALSE], Y_raw[test_idx, , drop = FALSE])
  # 
  #     #X_train <- X_scaled$train
  #     #X_test <- X_scaled$test
  #     #Y_train <- Y_scaled$train
  #     #Y_test <- Y_scaled$test
  #     X_train <- X_raw[-test_idx, , drop = FALSE]
  #     X_test <- X_raw[test_idx, , drop = FALSE]
  #     Y_train <- Y_raw[-test_idx, , drop = FALSE]
  #     Y_test <- Y_raw[test_idx, , drop = FALSE]
  # 
  #     for (method in params$methods) {
  #       message(sprintf("Repeat %d/%d | fold %d/%d | %s", repeat_id, params$repeats, fold_id, params$outer_folds, method))
  #       start_time <- proc.time()[["elapsed"]]
  # 
  #       timed_fit <- run_method_with_timeout(
  #         function() benchmark_fitters[[method]](X_train, Y_train),
  #         params$method_timeout_sec
  #       )
  # 
  #       runtime_sec <- proc.time()[["elapsed"]] - start_time
  #       lambda_x <- NA_real_
  #       lambda_y <- NA_real_
  #       rho <- NA_real_
  #       k <- NA_real_
  #       train_mse <- NA_real_
  #       train_cor <- NA_real_
  #       train_dist <- NA_real_
  #       test_mse <- NA_real_
  #       test_cor <- NA_real_
  #       test_dist <- NA_real_
  #       test_genotype_cluster_accuracy <- NA_real_
  #       test_genotype_cluster_ari <- NA_real_
  #       artifact_payload <- list(status = "error", error_message = timed_fit$error_message)
  #       benchmark_task_payload <- list()
  # 
  #       if (!is.null(timed_fit$result)) {
  #         fit <- timed_fit$result
  #         lambda_x <- if (!is.null(fit$lambda_x)) fit$lambda_x else NA_real_
  #         lambda_y <- if (!is.null(fit$lambda_y)) fit$lambda_y else NA_real_
  #         rho <- if (!is.null(fit$rho)) fit$rho else NA_real_
  #         k <- if (!is.null(fit$k)) fit$k else NA_real_
  # 
  #         train_metrics <- calc_metrics(fit$U, fit$V, X_train, Y_train, "train")
  #         test_metrics <- calc_metrics(fit$U, fit$V, X_test, Y_test, "test")
  # 
  #         train_mse <- train_metrics$train_mse
  #         train_cor <- train_metrics$train_cor
  #         train_dist <- train_metrics$train_dist
  #         test_mse <- test_metrics$test_mse
  #         test_cor <- test_metrics$test_cor
  #         test_dist <- test_metrics$test_dist
  # 
  #         test_embedding <- compute_full_fit_embedding(fit$U, X_test)
  #         genotype_cluster <- tryCatch(
  #           evaluate_embedding_cluster_task(test_embedding, test_metadata$genotype, "genotype"),
  #           error = function(e) NULL
  #         )
  #         if (!is.null(genotype_cluster)) {
  #           test_genotype_cluster_accuracy <- genotype_cluster$summary$genotype_cluster_accuracy
  #           test_genotype_cluster_ari <- genotype_cluster$summary$genotype_cluster_ari
  #           benchmark_task_payload <- list(
  #             test_embedding = test_embedding,
  #             test_genotype_cluster = genotype_cluster$artifact$cluster
  #           )
  #         }
  # 
  #         artifact_payload <- if (!is.null(fit$artifact)) fit$artifact else fit
  #       }
  # 
  #       artifact_path <- save_benchmark_artifact(
  #         out_dir = params$out_dir,
  #         repeat_id = repeat_id,
  #         fold_id = fold_id,
  #         method = method,
  #         artifact = list(
  #           method = method,
  #           repeat_id = repeat_id,
  #           fold_id = fold_id,
  #           lambda_x = lambda_x,
  #           lambda_y = lambda_y,
  #           rho = rho,
  #           k = k,
  #           runtime_sec = runtime_sec,
  #           timed_out = timed_fit$timed_out,
  #           error_message = timed_fit$error_message,
  #           result = artifact_payload,
  #           benchmark_tasks = benchmark_task_payload
  #         )
  #       )
  # 
  #       row_idx <- row_idx + 1L
  #       benchmark_rows[[row_idx]] <- data.frame(
  #         repeat_id = repeat_id,
  #         fold_id = fold_id,
  #         train_n = nrow(X_train),
  #         test_n = nrow(X_test),
  #         method = method,
  #         train_mse = train_mse,
  #         train_cor = train_cor,
  #         train_dist = train_dist,
  #         test_mse = test_mse,
  #         test_cor = test_cor,
  #         test_dist = test_dist,
  #         test_genotype_cluster_accuracy = test_genotype_cluster_accuracy,
  #         test_genotype_cluster_ari = test_genotype_cluster_ari,
  #         lambda_x = lambda_x,
  #         lambda_y = lambda_y,
  #         rho = rho,
  #         k = k,
  #         runtime_sec = runtime_sec,
  #         timed_out = timed_fit$timed_out,
  #         error_message = ifelse(is.na(timed_fit$error_message), "", timed_fit$error_message),
  #         artifact_path = artifact_path,
  #         stringsAsFactors = FALSE
  #       )
  # 
  #       write_benchmark_outputs(bind_rows(benchmark_rows), params$out_dir)
  #     }
  #   }
  # }
  # 
  # benchmark_df <- bind_rows(benchmark_rows)
  # benchmark_summary <- write_benchmark_outputs(benchmark_df, params$out_dir)
  # 
  # optimal_params <- lapply(params$methods, function(method) {
  #   choose_optimal_params(benchmark_df %>% filter(method == !!method))
  # }) %>%
  #   bind_rows()
  # 
  # utils::write.csv(
  #   optimal_params,
  #   file.path(params$out_dir, "optimal_parameters.csv"),
  #   row.names = FALSE
  # )

  X_full <- full_data$X_full
  Y_full <- full_data$Y_full
  full_model_cv_rows <- list()
  full_fit_rows <- list()
  full_fit_artifacts <- list()

  for (method in params$methods) {
    message(
      sprintf(
        "Full-data fit with %d-fold CV: %s",
        params$full_model_cv_folds,
        method
      )
    )
    start_time <- proc.time()[["elapsed"]]
    timed_fit <- run_method_with_timeout(
      function() full_model_cv_fitters[[method]](X_full, Y_full),
      params$method_timeout_sec
    )
    runtime_sec <- proc.time()[["elapsed"]] - start_time

    fit <- timed_fit$result
    full_model_cv_row <- data.frame(
      method = method,
      lambda_x = if (!is.null(fit$lambda_x)) fit$lambda_x else NA_real_,
      lambda_y = if (!is.null(fit$lambda_y)) fit$lambda_y else NA_real_,
      rho = if (!is.null(fit$rho)) fit$rho else NA_real_,
      k = if (!is.null(fit$k)) fit$k else NA_real_,
      cv_folds = as.numeric(params$full_model_cv_folds),
      selection_runtime_sec = runtime_sec,
      timed_out = timed_fit$timed_out,
      error_message = ifelse(is.na(timed_fit$error_message), "", timed_fit$error_message),
      stringsAsFactors = FALSE
    )
    full_model_cv_rows[[method]] <- full_model_cv_row

    if (is.null(fit)) {
      warning(
        "Skipping full-data fit for ",
        method,
        ": full-model CV failed",
        if (!is.na(timed_fit$error_message)) paste0(" (", timed_fit$error_message, ")") else "."
      )
      next
    }

    full_result <- build_full_fit_row(
      method,
      fit,
      X_full,
      Y_full,
      metadata,
      full_model_cv_row[1, , drop = FALSE]
    )
    full_fit_rows[[method]] <- full_result$summary_row
    full_fit_artifacts[[method]] <- full_result$artifact
  }

  full_model_cv_summary <- bind_rows(full_model_cv_rows)
  utils::write.csv(
    full_model_cv_summary,
    file.path(params$out_dir, "full_model_cv_summary.csv"),
    row.names = FALSE
  )

  full_fit_summary <- bind_rows(full_fit_rows) %>% arrange(full_mse)
  diet_prediction_results <- run_full_model_diet_tasks(
    full_fit_artifacts,
    X_full,
    metadata,
    out_dir = params$out_dir
  )

  utils::write.csv(
    full_fit_summary,
    file.path(params$out_dir, "full_data_fit_summary.csv"),
    row.names = FALSE
  )

  saveRDS(
    list(
      params = params,
      #benchmark_runs = benchmark_df,
      #benchmark_summary = benchmark_summary,
      #optimal_parameters = optimal_params,
      full_model_cv = full_model_cv_summary,
      full_data_fits = full_fit_artifacts,
      full_model_diet = diet_prediction_results
    ),
    file.path(params$out_dir, "full_data_fits.rds")
  )

  message("Benchmark and full-data fits stored in: ", params$out_dir)
}
 
if (sys.nframe() == 0) {
  main()
}
