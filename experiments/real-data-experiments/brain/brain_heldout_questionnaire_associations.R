raw_args <- commandArgs(trailingOnly = FALSE)
args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!startsWith(arg, "--")) {
      next
    }
    pieces <- strsplit(sub("^--", "", arg), "=", fixed = TRUE)[[1]]
    key <- pieces[1]
    value <- if (length(pieces) > 1) {
      paste(pieces[-1], collapse = "=")
    } else {
      "TRUE"
    }
    out[[key]] <- value
  }
  out
}

opts <- parse_args(args)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || identical(x, "")) y else x
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

as_num <- function(x, default) {
  if (is.null(x) || identical(x, "")) {
    return(default)
  }
  as.numeric(x)
}

script_file <- grep("^--file=", raw_args, value = TRUE)
script_dir <- if (length(script_file) > 0) {
  normalizePath(dirname(sub("^--file=", "", script_file[1])),
                winslash = "/", mustWork = TRUE)
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

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

project_root <- find_project_root(script_dir)
setwd(project_root)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(readxl)
  library(tidyr)
  library(magrittr)
  library(Matrix)
  library(tibble)
})

source_existing_files <- function(paths) {
  for (path in paths) {
    source(path)
  }
}

source_filtered_r_file <- function(path, replacements = list()) {
  lines <- readLines(path, warn = FALSE)
  if (length(replacements) > 0) {
    for (replacement in replacements) {
      lines <- gsub(
        replacement$pattern,
        replacement$replacement,
        lines,
        perl = TRUE
      )
    }
  }
  eval(parse(text = lines), envir = .GlobalEnv)
}

resolve_ccar3_source_dir <- function(candidate) {
  if (is.null(candidate) || !nzchar(candidate)) {
    return(NULL)
  }

  candidate <- path.expand(candidate)
  if (!dir.exists(candidate)) {
    return(NULL)
  }

  if (file.exists(file.path(candidate, "helpers.r")) &&
      file.exists(file.path(candidate, "reduced_rank_regression.R"))) {
    return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
  }

  r_dir <- file.path(candidate, "R")
  if (file.exists(file.path(r_dir, "helpers.r")) &&
      file.exists(file.path(r_dir, "reduced_rank_regression.R"))) {
    return(normalizePath(r_dir, winslash = "/", mustWork = TRUE))
  }

  NULL
}

resolve_ccar3_repo_root <- function(candidate) {
  source_dir <- resolve_ccar3_source_dir(candidate)
  if (is.null(source_dir)) {
    return(NULL)
  }

  if (file.exists(file.path(source_dir, "DESCRIPTION"))) {
    return(source_dir)
  }

  repo_root <- dirname(source_dir)
  if (file.exists(file.path(repo_root, "DESCRIPTION"))) {
    return(normalizePath(repo_root, winslash = "/", mustWork = TRUE))
  }

  NULL
}

patch_ccar3_parallel_worker_loader <- function(repo_root) {
  if (!exists("initialize_parallel_workers", mode = "function", inherits = TRUE)) {
    return(invisible(NULL))
  }

  original_fun <- get("initialize_parallel_workers", mode = "function", inherits = TRUE)
  patched_fun <- function(cl, pkg = "ccar3", pkg_path = NULL, verbose = FALSE) {
    if (is.null(pkg_path) || !nzchar(pkg_path)) {
      env_pkg_path <- Sys.getenv("CCAR3_PKG_PATH", unset = "")
      if (nzchar(env_pkg_path)) {
        pkg_path <- env_pkg_path
      } else {
        pkg_path <- repo_root
      }
    }

    original_fun(cl = cl, pkg = pkg, pkg_path = pkg_path, verbose = verbose)
  }

  patch_env <- new.env(parent = environment(original_fun))
  patch_env$original_fun <- original_fun
  patch_env$repo_root <- repo_root
  environment(patched_fun) <- patch_env
  assign("initialize_parallel_workers", patched_fun, envir = .GlobalEnv)
  invisible(NULL)
}

create_cv_folds <- function(indices, k, list = TRUE, returnTrain = FALSE) {
  indices <- as.integer(indices)
  n <- length(indices)
  if (n < 2) {
    stop("Need at least 2 observations to create cross-validation folds.", call. = FALSE)
  }

  k <- max(2L, min(as.integer(k), n))
  shuffled <- sample(indices, n)
  fold_id <- rep(seq_len(k), length.out = n)
  folds <- split(shuffled, fold_id)

  if (isTRUE(returnTrain)) {
    train_folds <- lapply(folds, function(test_idx) setdiff(indices, test_idx))
    if (!isTRUE(list)) {
      return(train_folds)
    }
    return(train_folds)
  }

  if (!isTRUE(list)) {
    fold_lookup <- integer(max(indices))
    for (i in seq_along(folds)) {
      fold_lookup[folds[[i]]] <- i
    }
    return(fold_lookup[indices])
  }

  folds
}

source_ccar3_package_source <- function(ccar3_dir) {
  ccar3_source_dir <- resolve_ccar3_source_dir(ccar3_dir)
  ccar3_repo_root <- resolve_ccar3_repo_root(ccar3_dir)
  if (is.null(ccar3_source_dir)) {
    stop(
      "Could not locate CCAR3 repo source under ",
      ccar3_dir,
      ". Expected either the repo root or its R/ directory.",
      call. = FALSE
    )
  }

  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' must be installed to use CCAR3 methods.", call. = FALSE)
  }
  suppressPackageStartupMessages(library(foreach))

  has_caret <- requireNamespace("caret", quietly = TRUE)
  replacements <- list()
  if (!has_caret) {
    replacements <- list(
      list(pattern = "^library\\((caret)\\).*", replacement = ""),
      list(pattern = "caret::createFolds", replacement = "create_cv_folds"),
      list(pattern = "\\.packages = c\\(\"Matrix\", \"caret\"\\)", replacement = ".packages = c(\"Matrix\")")
    )
  }

  source_filtered_r_file(file.path(ccar3_source_dir, "helpers.r"))
  source_filtered_r_file(file.path(ccar3_source_dir, "utils.R"))
  source_filtered_r_file(
    file.path(ccar3_source_dir, "reduced_rank_regression.R"),
    replacements = replacements
  )
  source_filtered_r_file(
    file.path(ccar3_source_dir, "group_reduced_rank_regression.R"),
    replacements = replacements
  )
  source_filtered_r_file(
    file.path(ccar3_source_dir, "graph_reduced_rank_regression.R"),
    replacements = replacements
  )
  patch_ccar3_parallel_worker_loader(ccar3_repo_root)
  invisible(ccar3_source_dir)
}

source_ccar3_methods <- function(ccar3_dir = "~/Documents/ccar3/R") {
  package_candidates <- c(
    ccar3_dir,
    Sys.getenv("CCAR3_PKG_PATH", unset = ""),
    file.path(path.expand("~"), "Documents", "ccar3", "R"),
    file.path(path.expand("~"), "Documents", "ccar3")
  )
  package_candidates <- unique(package_candidates[nzchar(package_candidates)])

  for (candidate in package_candidates) {
    source_dir <- resolve_ccar3_source_dir(candidate)
    repo_root <- resolve_ccar3_repo_root(candidate)
    if (is.null(source_dir)) {
      next
    }

    options(ccar3_pkg_path = repo_root %||% source_dir)
    Sys.setenv(CCAR3_PKG_PATH = repo_root %||% source_dir)
    message("Loading CCAR3 methods from repo source: ", source_dir)
    source_ccar3_package_source(source_dir)
    return(invisible(source_dir))
  }

  stop(
    "Could not load CCAR3 core methods from repo source. Checked:\n- ",
    paste(package_candidates, collapse = "\n- "),
    call. = FALSE
  )
}

source_sar_method <- function(project_root) {
  source_existing_files(file.path(
    project_root,
    c(
      "experiments/alternative_methods/SAR.R",
      "experiments/alternative_methods/GCA/gca_to_cca.R"
    )
  ))
}

resolve_existing_file <- function(candidates, label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    sprintf(
      "Could not find %s. Checked:\n- %s",
      label,
      paste(candidates, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

resolve_existing_dir <- function(candidates, label) {
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  for (candidate in candidates) {
    if (dir.exists(candidate)) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    sprintf(
      "Could not find %s. Checked:\n- %s",
      label,
      paste(candidates, collapse = "\n- ")
    ),
    call. = FALSE
  )
}

drop_non_numeric_columns <- function(df) {
  numeric_df <- df[, vapply(df, is.numeric, logical(1)), drop = FALSE]
  as.matrix(numeric_df)
}

read_matrix_csv <- function(path) {
  df <- if (requireNamespace("data.table", quietly = TRUE)) {
    data.table::fread(path, data.table = FALSE)
  } else {
    suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  }
  mat <- drop_non_numeric_columns(df)
  if (ncol(mat) == 0) {
    stop("No numeric columns found in ", path, call. = FALSE)
  }
  mat
}

prepare_train_test_matrix <- function(train_mat, test_mat = NULL, scale = FALSE) {
  train_mat <- as.matrix(train_mat)
  test_mat <- if (is.null(test_mat)) NULL else as.matrix(test_mat)

  train_means <- colMeans(train_mat, na.rm = TRUE)

  impute_with_reference <- function(mat, reference_means) {
    if (is.null(mat)) {
      return(NULL)
    }

    mat <- as.matrix(mat)
    for (j in seq_len(ncol(mat))) {
      missing <- is.na(mat[, j])
      if (any(missing)) {
        mat[missing, j] <- reference_means[j]
      }
    }
    mat
  }

  train_mat <- impute_with_reference(train_mat, train_means)
  test_mat <- impute_with_reference(test_mat, train_means)

  centered_train <- sweep(train_mat, 2, train_means, FUN = "-")
  centered_test <- if (is.null(test_mat)) NULL else sweep(test_mat, 2, train_means, FUN = "-")

  if (!isTRUE(scale)) {
    return(list(
      train = centered_train,
      test = centered_test,
      center = train_means,
      scale = rep(1, length(train_means))
    ))
  }

  train_sds <- apply(centered_train, 2, stats::sd)
  train_sds[!is.finite(train_sds) | train_sds == 0] <- 1

  list(
    train = sweep(centered_train, 2, train_sds, FUN = "/"),
    test = if (is.null(centered_test)) NULL else sweep(centered_test, 2, train_sds, FUN = "/"),
    center = train_means,
    scale = train_sds
  )
}

infer_positions_name <- function(positions) {
  if ("name" %in% names(positions)) {
    return(positions$name)
  }
  if ("Region" %in% names(positions)) {
    return(vapply(
      strsplit(positions$Region, "_", fixed = TRUE),
      function(parts) {
        if (length(parts) >= 3 &&
            parts[1] %in% c("L", "R") &&
            parts[length(parts)] %in% c("L", "R")) {
          return(paste(parts[2:(length(parts) - 1)], collapse = "_"))
        }
        if (length(parts) >= 2) {
          return(parts[2])
        }
        parts[1]
      },
      character(1)
    ))
  }
  rep(NA_character_, nrow(positions))
}

build_positions_table <- function(coordinates_path, label_path = NULL) {
  positions <- suppressMessages(readr::read_csv(coordinates_path, show_col_types = FALSE))
  required_cols <- c("x", "y", "z")
  if (!all(required_cols %in% names(positions))) {
    stop("Coordinates file must contain columns x, y, z.", call. = FALSE)
  }

  positions$name <- infer_positions_name(positions)

  if ("Gyrus" %in% names(positions)) {
    return(positions)
  }

  if (is.null(label_path) || !file.exists(label_path)) {
    return(positions)
  }

  labels <- readxl::read_xlsx(label_path)
  if (!all(c("Label ID.L", "...6") %in% names(labels))) {
    return(positions)
  }

  gyrus_col <- if ("Gyrus...10" %in% names(labels)) {
    "Gyrus...10"
  } else if ("Gyrus" %in% names(labels)) {
    "Gyrus"
  } else {
    NULL
  }

  labels <- labels %>%
    dplyr::filter(`Label ID.L` < 211) %>%
    dplyr::mutate(
      name = vapply(strsplit(`...6`, ","), function(x) x[1], character(1)),
      Gyrus = vapply(
        strsplit(`...6`, ","),
        function(x) ifelse(length(x) >= 2, trimws(x[2]), NA_character_),
        character(1)
      )
    )

  if (!is.null(gyrus_col)) {
    labels <- labels %>%
      tidyr::fill(dplyr::all_of(gyrus_col)) %>%
      dplyr::mutate(Gyrus = dplyr::coalesce(.data[[gyrus_col]], Gyrus))
  }

  labels <- labels %>%
    dplyr::select(name, Gyrus)

  positions <- positions %>%
    dplyr::left_join(labels, by = "name")

  if (nrow(positions) >= 19) {
    positions$Gyrus[1:19] <- positions$name[1:19]
  }
  if (nrow(positions) >= 15) {
    positions$Gyrus[15] <- "Brain Stem"
    positions$name[15] <- "Brain Stem"
  }

  positions
}

read_roi_group_metadata <- function(group_path) {
  if (is.null(group_path) || !file.exists(group_path)) {
    return(NULL)
  }

  group_df <- if (grepl("\\.xlsx$", group_path, ignore.case = TRUE)) {
    readxl::read_xlsx(group_path)
  } else {
    suppressMessages(readr::read_csv(group_path, show_col_types = FALSE))
  }

  if (!"Gyrus" %in% names(group_df)) {
    return(NULL)
  }

  group_df
}

apply_roi_group_metadata <- function(positions, group_path = NULL) {
  group_df <- read_roi_group_metadata(group_path)
  if (is.null(group_df)) {
    return(positions)
  }

  join_keys <- intersect(c("Index", "Region", "name"), intersect(names(positions), names(group_df)))
  if (length(join_keys) == 0) {
    return(positions)
  }

  group_df <- group_df %>%
    dplyr::select(dplyr::all_of(join_keys), Gyrus) %>%
    dplyr::distinct()

  positions %>%
    dplyr::left_join(group_df, by = join_keys, suffix = c("", ".from_file")) %>%
    dplyr::mutate(Gyrus = dplyr::coalesce(.data[["Gyrus.from_file"]], Gyrus)) %>%
    dplyr::select(-dplyr::any_of("Gyrus.from_file"))
}

read_group_assignments <- function(group_path) {
  if (is.null(group_path) || !file.exists(group_path)) {
    return(NULL)
  }

  if (grepl("\\.xlsx$", group_path, ignore.case = TRUE)) {
    group_df <- readxl::read_xlsx(group_path, col_names = FALSE)
  } else {
    group_df <- suppressMessages(
      readr::read_csv(group_path, show_col_types = FALSE, col_names = FALSE)
    )
  }

  suppressWarnings(as.integer(group_df[[1]]))
}

normalize_group_assignments <- function(group_ids, expected_n_features) {
  if (is.null(group_ids)) {
    return(NULL)
  }

  group_ids <- as.integer(group_ids)
  if (length(group_ids) == expected_n_features) {
    return(group_ids)
  }

  filtered_ids <- group_ids[!is.na(group_ids) & group_ids != 0]
  if (length(filtered_ids) == expected_n_features) {
    return(filtered_ids)
  }

  NULL
}

aggregate_features_by_group <- function(X, group_ids) {
  keep <- !is.na(group_ids) & group_ids != 0
  if (!all(keep)) {
    X <- X[, keep, drop = FALSE]
    group_ids <- group_ids[keep]
  }

  grouped_sums <- rowsum(t(X), group = as.integer(group_ids), reorder = TRUE)
  group_sizes <- as.numeric(table(group_ids)[rownames(grouped_sums)])
  t(grouped_sums / group_sizes)
}

align_brain_matrix_to_positions <- function(X, positions, group_path = NULL) {
  if (ncol(X) == nrow(positions)) {
    return(list(
      X = X,
      positions = positions,
      aggregated = FALSE,
      message = NULL
    ))
  }

  group_ids_raw <- read_group_assignments(group_path)
  normalized_group_ids <- normalize_group_assignments(group_ids_raw, ncol(X))
  if (is.null(normalized_group_ids)) {
    stop(
      sprintf(
        "X has %d features but ROI metadata has %d rows, and %s could not be aligned to the X columns.",
        ncol(X), nrow(positions), group_path %||% "the group assignment file"
      ),
      call. = FALSE
    )
  }

  aggregated_X <- aggregate_features_by_group(X, normalized_group_ids)
  colnames(aggregated_X) <- positions$name

  list(
    X = aggregated_X,
    positions = positions,
    aggregated = TRUE,
    message = sprintf(
      "Aggregated voxel-level X from %d features to %d parcels using %s.",
      ncol(X), ncol(aggregated_X), basename(group_path)
    )
  )
}

build_groups <- function(positions) {
  if (!"Gyrus" %in% names(positions)) {
    stop("Could not build ROI groups: no Gyrus labels available.", call. = FALSE)
  }
  split(seq_len(nrow(positions)), positions$Gyrus)
}

build_incidence_matrix <- function(edges, n_nodes, weight = 1) {
  Gamma <- matrix(0, nrow(edges), n_nodes)
  for (edge_idx in seq_len(nrow(edges))) {
    Gamma[edge_idx, edges$Source[edge_idx]] <- weight
    Gamma[edge_idx, edges$Target[edge_idx]] <- -weight
  }
  Gamma
}

build_knn_graph_incidence <- function(positions, k = 4) {
  coords <- as.matrix(positions[, c("x", "y", "z")])
  dist_mat <- as.matrix(dist(coords))
  edge_pairs <- vector("list", nrow(dist_mat))

  for (i in seq_len(nrow(dist_mat))) {
    nn_idx <- order(dist_mat[i, ])[seq_len(k + 1)]
    edge_pairs[[i]] <- cbind(rep(i, k), nn_idx[-1])
  }

  edges <- do.call(rbind, edge_pairs) %>%
    as.data.frame()
  names(edges) <- c("Source", "Target")
  edges <- edges %>%
    dplyr::mutate(
      source_min = pmin(Source, Target),
      target_max = pmax(Source, Target)
    ) %>%
    dplyr::transmute(Source = source_min, Target = target_max) %>%
    dplyr::distinct()

  build_incidence_matrix(edges, nrow(positions), weight = 1)
}

load_benchmark_lambdas <- function(path, methods) {
  df <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  required_cols <- c("method", "fold", "split", "lambda")
  if (!all(required_cols %in% names(df))) {
    stop("Benchmark results CSV is missing required columns.", call. = FALSE)
  }

  df %>%
    dplyr::filter(split == "train", method %in% methods) %>%
    dplyr::group_by(method, fold) %>%
    dplyr::summarise(lambda = dplyr::first(lambda), .groups = "drop")
}

fit_ccar3_final <- function(method, X_train, Y_train, lambda, r, groups, Gamma,
                            niter, rho, thresh) {
  if (method == "cca_rrr_cv") {
    return(cca_rrr(
      X = X_train,
      Y = Y_train,
      lambda = lambda,
      r = r,
      standardize = FALSE,
      preprocess = "none",
      LW_Sy = TRUE,
      rho = rho,
      niter = niter,
      thresh = thresh
    ))
  }

  if (method == "cca_group_rrr_cv") {
    return(cca_group_rrr(
      X = X_train,
      Y = Y_train,
      groups = groups,
      lambda = lambda,
      r = r,
      standardize = FALSE,
      preprocess = "none",
      LW_Sy = TRUE,
      rho = rho,
      niter = niter,
      thresh = thresh
    ))
  }

  if (method == "cca_graph_rrr_cv") {
    return(cca_graph_rrr(
      X = X_train,
      Y = Y_train,
      Gamma = Gamma,
      lambda = lambda,
      r = r,
      standardize = FALSE,
      preprocess = "none",
      LW_Sy = TRUE,
      rho = rho,
      niter = niter,
      thresh = thresh
    ))
  }

  stop("Unknown CCAR3 method: ", method, call. = FALSE)
}

fit_sar_cv_train <- function(X_train, Y_train, r, lambda_grid) {
  method <- SparseCCA(
    X = X_train,
    Y = Y_train,
    rank = r,
    lambdaAseq = lambda_grid,
    lambdaBseq = c(0),
    standardize = FALSE,
    selection.criterion = 2,
    n.cv = 5,
    max.iter = 100,
    conv = 1e-2
  )

  S <- cov(cbind(X_train, Y_train))
  gca_to_cca(rbind(method$uhat, method$vhat), S, c(ncol(X_train), ncol(Y_train)))
}

extract_method_uv <- function(fit) {
  U <- fit$U %||% fit$u
  V <- fit$V %||% fit$v
  if (is.null(U) || is.null(V)) {
    stop("Could not extract U/V from fit object.", call. = FALSE)
  }
  list(U = as.matrix(U), V = as.matrix(V))
}

compute_fit_metrics <- function(X, Y, U, V, rank_r) {
  X_scores <- as.matrix(X) %*% as.matrix(U)
  Y_scores <- as.matrix(Y) %*% as.matrix(V)
  n <- nrow(X_scores)
  correlations <- as.numeric(diag(crossprod(X_scores, Y_scores) / n))

  if (length(correlations) < rank_r) {
    correlations <- c(correlations, rep(NA_real_, rank_r - length(correlations)))
  } else if (length(correlations) > rank_r) {
    correlations <- correlations[seq_len(rank_r)]
  }

  list(
    X_scores = X_scores,
    Y_scores = Y_scores,
    global_mse = mean((X_scores - Y_scores)^2),
    correlations = correlations
  )
}

outcome_labels <- c(
  bas_drive = "BAS Drive",
  bas_funseeking = "BAS Fun Seeking",
  bas_rewardrespons = "BAS Reward Responsiveness",
  bis_total = "BIS Total",
  masq_distress = "MASQ Distress",
  masq_anhedonia = "MASQ Anhedonia",
  masq_anxarousal = "MASQ Anxious Arousal",
  panas_positive = "PANAS Positive Affect",
  panas_negative = "PANAS Negative Affect"
)

evaluate_questionnaire_associations <- function(method, fold_id, split_name, scores, Y,
                                                fit_time_sec, selected_lambda = NA_real_) {
  score_df <- as.data.frame(scores)
  names(score_df) <- paste0("score_", seq_len(ncol(score_df)))

  outcome_tables <- lapply(seq_len(ncol(Y)), function(j) {
    y <- Y[, j]
    outcome_name <- colnames(Y)[j]
    fit_df <- data.frame(y = y, score_df, check.names = FALSE)
    fit_df <- fit_df[stats::complete.cases(fit_df), , drop = FALSE]

    fit <- stats::lm(y ~ ., data = fit_df)
    fit_summary <- summary(fit)
    fstat <- fit_summary$fstatistic
    overall_p_value <- if (is.null(fstat) || any(!is.finite(fstat))) {
      NA_real_
    } else {
      stats::pf(fstat["value"], fstat["numdf"], fstat["dendf"], lower.tail = FALSE)
    }

    tibble::tibble(
      fold = fold_id,
      split = split_name,
      method = method,
      outcome = outcome_name,
      outcome_label = outcome_labels[[outcome_name]] %||% outcome_name,
      adjusted_r2 = unname(fit_summary$adj.r.squared),
      r2 = unname(fit_summary$r.squared),
      overall_p_value = unname(overall_p_value),
      significant_0_05 = is.finite(overall_p_value) && overall_p_value < 0.05,
      n = nrow(fit_df),
      n_predictors = ncol(score_df),
      lambda = selected_lambda,
      fit_time_sec = fit_time_sec
    )
  })

  dplyr::bind_rows(outcome_tables)
}

build_fit_metadata_row <- function(fold_id, method, fit_time_sec, selected_lambda,
                                   global_mse, correlations, rank_r) {
  correlation_values <- correlations
  if (length(correlation_values) < rank_r) {
    correlation_values <- c(correlation_values, rep(NA_real_, rank_r - length(correlation_values)))
  } else if (length(correlation_values) > rank_r) {
    correlation_values <- correlation_values[seq_len(rank_r)]
  }

  tibble::as_tibble_row(c(
    list(
      fold = fold_id,
      method = method,
      fit_time_sec = fit_time_sec,
      lambda = selected_lambda,
      global_mse = global_mse
    ),
    stats::setNames(as.list(correlation_values), paste0("correlation_", seq_len(rank_r)))
  ))
}

method_order <- c("cca_rrr_cv", "cca_group_rrr_cv", "cca_graph_rrr_cv", "SAR_CV")
ccar3_methods <- c("cca_rrr_cv", "cca_group_rrr_cv", "cca_graph_rrr_cv")

data_dir <- resolve_existing_dir(
  c(
    opts$data_dir,
    Sys.getenv("BRAIN_DATA_DIR", unset = ""),
    file.path(script_dir, "data"),
    file.path(project_root, "experiments", "data", "fMRI-data"),
    file.path(project_root, "data", "fMRI-data"),
    file.path(project_root, "experiments", "data"),
    file.path(project_root, "data")
  ),
  "brain data directory"
)

output_dir <- normalizePath(
  file.path(script_dir, "results", "heldout_questionnaire_associations"),
  winslash = "/",
  mustWork = FALSE
)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

x_path <- resolve_existing_file(
  c(
    opts$x_path,
    Sys.getenv("BRAIN_X_PATH", unset = ""),
    file.path(data_dir, "activations_X_preprocessed.csv"),
    file.path(project_root, "data", "activations_X_preprocessed.csv")
  ),
  "brain X matrix"
)

y_path <- resolve_existing_file(
  c(
    opts$y_path,
    Sys.getenv("BRAIN_Y_PATH", unset = ""),
    file.path(data_dir, "activations_Y_preprocessed.csv"),
    file.path(project_root, "data", "activations_Y_preprocessed.csv")
  ),
  "brain Y matrix"
)

coordinates_path <- resolve_existing_file(
  c(
    opts$coordinates_path,
    Sys.getenv("BRAIN_COORDINATES_PATH", unset = ""),
    file.path(data_dir, "parcellation_coordinates.csv"),
    file.path(project_root, "data", "parcellation_coordinates.csv")
  ),
  "ROI coordinate metadata"
)

label_candidates <- c(
  opts$label_path,
  Sys.getenv("BRAIN_LABEL_PATH", unset = ""),
  file.path(data_dir, "BNA_subregions.xlsx"),
  file.path(data_dir, "data", "atlases", "BNA_subregions.xlsx"),
  file.path(project_root, "data", "BNA_subregions.xlsx")
)
label_path <- label_candidates[file.exists(label_candidates)][1]
if (is.na(label_path)) {
  label_path <- NULL
}

roi_group_candidates <- c(
  opts$roi_group_path,
  opts$group_path,
  Sys.getenv("BRAIN_GROUP_PATH", unset = ""),
  file.path(data_dir, "groups.csv"),
  file.path(project_root, "data", "groups.csv")
)
roi_group_path <- roi_group_candidates[file.exists(roi_group_candidates)][1]
if (is.na(roi_group_path)) {
  roi_group_path <- NULL
}

aggregation_group_candidates <- c(
  opts$aggregation_group_path,
  Sys.getenv("BRAIN_AGGREGATION_GROUP_PATH", unset = ""),
  file.path(data_dir, "activation_groups.xlsx"),
  file.path(data_dir, "activation_groups.csv"),
  file.path(project_root, "data", "activation_groups.xlsx"),
  file.path(project_root, "data", "activation_groups.csv")
)
aggregation_group_path <- aggregation_group_candidates[file.exists(aggregation_group_candidates)][1]
if (is.na(aggregation_group_path)) {
  aggregation_group_path <- NULL
}

benchmark_results_dir <- resolve_existing_dir(
  c(
    opts$benchmark_results_dir,
    file.path(project_root, "results_20_folds"),
    file.path(project_root, "results_old"),
    file.path(project_root, "results")
  ),
  "benchmark results directory"
)
benchmark_results_path <- resolve_existing_file(
  c(file.path(benchmark_results_dir, "brain_method_comparison_results.csv")),
  "benchmark results CSV"
)

rank_r <- as.integer(as_num(opts$r, 2))
outer_folds <- as.integer(as_num(opts$outer_folds, 20))
seed <- as.integer(as_num(opts$seed, 123))
lambda_grid <- 10 ^ seq(
  from = as_num(opts$lambda_min_log10, -3),
  to = as_num(opts$lambda_max_log10, 2),
  length.out = as.integer(as_num(opts$nlambda, 20))
)
niter <- as.integer(as_num(opts$niter, 10000))
rho <- as_num(opts$rho, 1)
thresh <- as_num(opts$thresh, 1e-5)
scale_inputs <- as_bool(opts$standardize_inputs, default = FALSE)

methods <- c("cca_rrr_cv", "cca_group_rrr_cv", "cca_graph_rrr_cv", "SAR_CV")
if (!is.null(opts$methods) && nzchar(opts$methods)) {
  requested <- trimws(strsplit(opts$methods, ",", fixed = TRUE)[[1]])
  methods <- intersect(method_order, requested)
}
if (length(methods) == 0) {
  stop("No valid methods selected.", call. = FALSE)
}

message("Reading X from: ", x_path)
message("Reading Y from: ", y_path)
message("Reading coordinates from: ", coordinates_path)
message("Using archived benchmark lambdas from: ", benchmark_results_path)

X_raw <- read_matrix_csv(x_path)
Y_raw <- read_matrix_csv(y_path)
Y_names <- names(readr::read_csv(y_path, show_col_types = FALSE))
positions <- build_positions_table(coordinates_path, label_path)
positions <- apply_roi_group_metadata(positions, roi_group_path)

aligned_X <- align_brain_matrix_to_positions(
  X_raw,
  positions,
  group_path = aggregation_group_path
)
X_raw <- aligned_X$X
positions <- aligned_X$positions
if (!is.null(aligned_X$message)) {
  message(aligned_X$message)
}

groups <- build_groups(positions)
Gamma <- build_knn_graph_incidence(positions, k = as.integer(as_num(opts$graph_k, 4)))

source_ccar3_methods(opts$ccar3_dir %||% "~/Documents/ccar3/R")
if ("SAR_CV" %in% methods) {
  source_sar_method(project_root)
}

lambda_df <- load_benchmark_lambdas(benchmark_results_path, ccar3_methods)

set.seed(seed)
folds <- create_cv_folds(seq_len(nrow(X_raw)), k = outer_folds, list = TRUE, returnTrain = FALSE)

fit_metadata <- list()
association_results <- list()

for (fold_id in seq_along(folds)) {
  test_idx <- sort(as.integer(folds[[fold_id]]))
  train_idx <- setdiff(seq_len(nrow(X_raw)), test_idx)

  X_train_raw <- X_raw[train_idx, , drop = FALSE]
  Y_train_raw <- Y_raw[train_idx, , drop = FALSE]
  X_test_raw <- X_raw[test_idx, , drop = FALSE]
  Y_test_raw <- Y_raw[test_idx, , drop = FALSE]

  x_fold <- prepare_train_test_matrix(X_train_raw, X_test_raw, scale = scale_inputs)
  y_fold <- prepare_train_test_matrix(Y_train_raw, Y_test_raw, scale = scale_inputs)

  X_train <- x_fold$train
  X_test <- x_fold$test
  Y_train <- y_fold$train
  Y_test <- y_fold$test

  colnames(Y_train) <- Y_names
  colnames(Y_test) <- Y_names

  message(sprintf("Held-out fold %d/%d", fold_id, length(folds)))

  for (method in methods) {
    message("  Evaluating ", method)
    start_time <- proc.time()[["elapsed"]]

    if (method == "SAR_CV") {
      fit <- fit_sar_cv_train(X_train, Y_train, r = rank_r, lambda_grid = lambda_grid)
      selected_lambda <- NA_real_
    } else {
      selected_lambda <- lambda_df %>%
        dplyr::filter(method == !!method, fold == !!fold_id) %>%
        dplyr::pull(lambda)
      if (length(selected_lambda) != 1 || !is.finite(selected_lambda)) {
        stop(
          "Could not recover a single archived lambda for method ",
          method,
          " and fold ",
          fold_id,
          ".",
          call. = FALSE
        )
      }

      fit <- fit_ccar3_final(
        method = method,
        X_train = X_train,
        Y_train = Y_train,
        lambda = selected_lambda,
        r = rank_r,
        groups = groups,
        Gamma = Gamma,
        niter = niter,
        rho = rho,
        thresh = thresh
      )
    }

    fit_time_sec <- proc.time()[["elapsed"]] - start_time
    uv <- extract_method_uv(fit)
    test_metrics <- compute_fit_metrics(X_test, Y_test, uv$U, uv$V, rank_r = rank_r)

    fit_metadata[[paste(method, fold_id, sep = "_")]] <- build_fit_metadata_row(
      fold_id = fold_id,
      method = method,
      fit_time_sec = fit_time_sec,
      selected_lambda = selected_lambda,
      global_mse = test_metrics$global_mse,
      correlations = test_metrics$correlations,
      rank_r = rank_r
    )

    association_results[[paste(method, fold_id, sep = "_")]] <- evaluate_questionnaire_associations(
      method = method,
      fold_id = fold_id,
      split_name = "test",
      scores = test_metrics$X_scores,
      Y = Y_test,
      fit_time_sec = fit_time_sec,
      selected_lambda = selected_lambda
    )
  }
}

fit_metadata_df <- dplyr::bind_rows(fit_metadata) %>%
  dplyr::mutate(method = factor(method, levels = method_order)) %>%
  dplyr::arrange(fold, method)

association_long_df <- dplyr::bind_rows(association_results) %>%
  dplyr::mutate(method = factor(method, levels = method_order)) %>%
  dplyr::arrange(outcome_label, method, fold)

association_summary_df <- association_long_df %>%
  dplyr::group_by(outcome, outcome_label, method) %>%
  dplyr::summarise(
    mean_adjusted_r2 = mean(adjusted_r2, na.rm = TRUE),
    sd_adjusted_r2 = stats::sd(adjusted_r2, na.rm = TRUE),
    mean_r2 = mean(r2, na.rm = TRUE),
    prop_significant_0_05 = mean(significant_0_05, na.rm = TRUE),
    mean_fit_time_sec = mean(fit_time_sec, na.rm = TRUE),
    mean_lambda = mean(lambda, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(outcome_label, method)

association_wide_df <- association_summary_df %>%
  dplyr::select(outcome, outcome_label, method, mean_adjusted_r2) %>%
  tidyr::pivot_wider(names_from = method, values_from = mean_adjusted_r2) %>%
  dplyr::arrange(match(outcome, Y_names))

fit_metadata_path <- file.path(output_dir, "brain_heldout_method_fit_metrics.csv")
association_long_path <- file.path(output_dir, "brain_heldout_questionnaire_associations_long.csv")
association_summary_path <- file.path(output_dir, "brain_heldout_questionnaire_associations_summary.csv")
association_wide_path <- file.path(output_dir, "brain_heldout_questionnaire_associations_adjusted_r2.csv")

readr::write_csv(fit_metadata_df, fit_metadata_path)
readr::write_csv(association_long_df, association_long_path)
readr::write_csv(association_summary_df, association_summary_path)
readr::write_csv(association_wide_df, association_wide_path)

message("Saved held-out fit metrics to: ", fit_metadata_path)
message("Saved held-out long association table to: ", association_long_path)
message("Saved held-out summary association table to: ", association_summary_path)
message("Saved held-out adjusted R^2 table to: ", association_wide_path)
print(association_wide_df)
