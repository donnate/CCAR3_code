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
  library(magrittr)
  library(Matrix)
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

  has_caret <- requireNamespace("caret", quietly = TRUE)
  replacements <- list(
    list(pattern = "^library\\((pracma)\\).*", replacement = ""),
    list(pattern = "corpcor::cov\\.shrink\\(Y(?:,\\s*verbose\\s*=\\s*verbose)?\\)", replacement = "stats::cov(Y)")
  )
  if (!has_caret) {
    replacements <- c(
      list(
        list(pattern = "^library\\((caret)\\).*", replacement = ""),
        list(pattern = "caret::createFolds", replacement = "create_cv_folds"),
        list(pattern = "\\.packages = c\\(\"Matrix\", \"caret\"\\)", replacement = ".packages = c(\"Matrix\")")
      ),
      replacements
    )
  }

  if (!requireNamespace("foreach", quietly = TRUE)) {
    stop("Package 'foreach' must be installed to use CCAR3 cross-validation methods.", call. = FALSE)
  }
  suppressPackageStartupMessages(library(foreach))

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


source_ccar3_methods <- function(project_root, ccar3_dir = "~/Documents/ccar3/R") {
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
    return(invisible("repo-source"))
  }

  stop(
    "Could not load CCAR3 core methods from repo source. Checked:\n- ",
    paste(package_candidates, collapse = "\n- "),
    call. = FALSE
  )
}

source_benchmark_method_stack <- function(project_root) {
  benchmark_paths <- file.path(
    project_root,
    c(
      "experiments/evaluation.R",
      "experiments/alternative_methods/SAR.R",
      "experiments/alternative_methods/Parkhomenko.R",
      "experiments/alternative_methods/Witten_CrossValidation.R",
      "experiments/alternative_methods/Waaijenborg.R",
      "experiments/alternative_methods/scca_chao.R",
      "experiments/alternative_methods/GCA/utils.R",
      "experiments/alternative_methods/GCA/gca_to_cca.R",
      "experiments/alternative_methods/GCA/init_process.R",
      "experiments/alternative_methods/GCA/sgca_init.R",
      "experiments/alternative_methods/GCA/sgca_tgd.R",
      "experiments/experiment_functions.R"
    )
  )
  source_existing_files(benchmark_paths)
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
  positions <- suppressMessages(readr::read_csv(coordinates_path, 
  show_col_types = FALSE))
  required_cols <- c("x", "y", "z")
  if (!all(required_cols %in% names(positions))) {
    stop("Coordinates file must contain columns x, y, z.", call. = FALSE)
  }

  positions$name <- infer_positions_name(positions)

  if ("Gyrus" %in% names(positions)) {
    return(positions)
  }

  if (is.null(label_path) || !file.exists(label_path)) {
    stop(
      "Need either a Gyrus column in the coordinates file or a valid BNA label file.",
      call. = FALSE
    )
  }

  labels <- readxl::read_xlsx(label_path)
  if (!all(c("Label ID.L", "...6") %in% names(labels))) {
    stop(
      "BNA label file is missing one of: `Label ID.L` or `...6`.",
      call. = FALSE
    )
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
    stop(
      "ROI group metadata must share one of these columns with the positions table: Index, Region, name.",
      call. = FALSE
    )
  }

  group_df <- group_df %>%
    dplyr::select(dplyr::all_of(join_keys), Gyrus) %>%
    dplyr::distinct()

  positions <- positions %>%
    dplyr::left_join(group_df, by = join_keys, suffix = c("", ".from_file")) %>%
    dplyr::mutate(Gyrus = dplyr::coalesce(.data[["Gyrus.from_file"]], Gyrus)) %>%
    dplyr::select(-dplyr::any_of("Gyrus.from_file"))

  positions
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

  group_ids <- as.integer(group_ids)
  grouped_sums <- rowsum(t(X), group = group_ids, reorder = TRUE)
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

  nonzero_ids <- normalized_group_ids[normalized_group_ids != 0 & !is.na(normalized_group_ids)]
  unique_ids <- sort(unique(nonzero_ids))
  if (length(unique_ids) != nrow(positions)) {
    stop(
      sprintf(
        "After filtering group assignments, found %d nonzero ROI ids but coordinates describe %d ROIs.",
        length(unique_ids), nrow(positions)
      ),
      call. = FALSE
    )
  }

  if (!identical(unique_ids, seq_len(nrow(positions)))) {
    stop(
      "Expected nonzero group ids to match parcel indices 1..nrow(positions).",
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

build_groups <- function(positions, group_path = NULL, expected_n_features = NULL) {
  groups_from_file <- NULL
  if (!is.null(group_path) && file.exists(group_path) && !is.null(expected_n_features)) {
    group_ids <- normalize_group_assignments(read_group_assignments(group_path), expected_n_features)
    if (!is.null(group_ids)) {
      keep <- which(group_ids != 0 & !is.na(group_ids))
      groups_from_file <- split(keep, group_ids[keep])
    }
  }

  if (!is.null(groups_from_file)) {
    return(groups_from_file)
  }

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

compute_inverse_sqrt_transform <- function(gram, jitter = 1e-8, max_tries = 6) {
  gram <- (as.matrix(gram) + t(as.matrix(gram))) / 2
  p <- nrow(gram)

  for (k in 0:max_tries) {
    jitter_value <- if (k == 0) 0 else jitter * (10 ^ (k - 1))
    chol_fit <- tryCatch(
      chol(gram + diag(jitter_value, p)),
      error = function(e) NULL
    )
    if (!is.null(chol_fit)) {
      return(list(
        transform = backsolve(chol_fit, diag(p)),
        jitter = jitter_value
      ))
    }
  }

  eig <- eigen(gram, symmetric = TRUE)
  vals <- pmax(eig$values, jitter)
  list(
    transform = eig$vectors %*% diag(1 / sqrt(vals), nrow = length(vals)) %*% t(eig$vectors),
    jitter = NA_real_
  )
}

standardize_train_scores <- function(X_train, Y_train, U, V, jitter = 1e-8) {
  U <- as.matrix(U)
  V <- as.matrix(V)

  X_scores <- as.matrix(X_train) %*% U
  Y_scores <- as.matrix(Y_train) %*% V

  x_transform <- compute_inverse_sqrt_transform(crossprod(X_scores) / nrow(X_train), jitter = jitter)
  y_transform <- compute_inverse_sqrt_transform(crossprod(Y_scores) / nrow(Y_train), jitter = jitter)

  U_std <- U %*% x_transform$transform
  V_std <- V %*% y_transform$transform

  x_crossprod_std <- crossprod(as.matrix(X_train) %*% U_std) / nrow(X_train)
  y_crossprod_std <- crossprod(as.matrix(Y_train) %*% V_std) / nrow(Y_train)

  list(
    U = U_std,
    V = V_std,
    x_crossprod = x_crossprod_std,
    y_crossprod = y_crossprod_std,
    x_jitter = x_transform$jitter,
    y_jitter = y_transform$jitter
  )
}

inspect_train_scores <- function(X_train, Y_train, U, V) {
  U <- as.matrix(U)
  V <- as.matrix(V)

  list(
    U = U,
    V = V,
    x_crossprod = crossprod(as.matrix(X_train) %*% U) / nrow(X_train),
    y_crossprod = crossprod(as.matrix(Y_train) %*% V) / nrow(Y_train),
    x_jitter = NA_real_,
    y_jitter = NA_real_
  )
}

score_pair <- function(X, Y, U, V, split, fold_id, method, lambda, fit_time_sec,
                       n_components = NULL) {
  U <- as.matrix(U)
  V <- as.matrix(V)
  XU <- as.matrix(X) %*% U
  YV <- as.matrix(Y) %*% V
  n <- nrow(X)
  correlation_values <- as.numeric(diag(crossprod(XU, YV)/n))
  if (!is.null(n_components)) {
    if (length(correlation_values) < n_components) {
      correlation_values <- c(
        correlation_values,
        rep(NA_real_, n_components - length(correlation_values))
      )
    } else if (length(correlation_values) > n_components) {
      correlation_values <- correlation_values[seq_len(n_components)]
    }
  }
  correlation_cols <- stats::setNames(
    as.list(correlation_values),
    paste0("correlation_", seq_along(correlation_values))
  )

  tibble::as_tibble_row(c(
    list(
      method = method,
      lambda = lambda,
      fit_time_sec = fit_time_sec,
      fold = fold_id,
      split = split,
      global_mse = mean((XU - YV) ^ 2)
    ),
    correlation_cols
  ))
}

sanitize_filename_component <- function(x) {
  cleaned <- gsub("[^A-Za-z0-9]+", "_", x)
  cleaned <- gsub("^_+|_+$", "", cleaned)
  if (!nzchar(cleaned)) "value" else cleaned
}

append_csv_rows <- function(df, path) {
  if (nrow(df) == 0) {
    return(invisible(FALSE))
  }
  append_mode <- file.exists(path)
  readr::write_csv(df, path, append = append_mode, col_names = !append_mode)
  invisible(TRUE)
}

checkpoint_result_path <- function(checkpoint_dir, fold_id, method, suffix = "results") {
  file.path(
    checkpoint_dir,
    sprintf(
      "brain_method_comparison_fold_%02d_%s_%s.csv",
      fold_id,
      sanitize_filename_component(method),
      suffix
    )
  )
}

remove_existing_checkpoints <- function(checkpoint_dir) {
  checkpoint_files <- list.files(
    checkpoint_dir,
    pattern = "^brain_method_comparison_fold_.*\\.csv$",
    full.names = TRUE
  )
  if (length(checkpoint_files) == 0) {
    return(invisible(0L))
  }

  removed <- file.remove(checkpoint_files)
  invisible(sum(removed, na.rm = TRUE))
}

remove_if_exists <- function(paths) {
  existing_paths <- paths[file.exists(paths)]
  if (length(existing_paths) == 0) {
    return(invisible(logical()))
  }
  invisible(file.remove(existing_paths))
}

build_summary_df <- function(results_df) {
  if (nrow(results_df) == 0) {
    return(tibble::tibble())
  }

  correlation_cols <- grep("^correlation_", names(results_df), value = TRUE)
  results_with_mean_corr <- results_df
  results_with_mean_corr$mean_correlation <- if (length(correlation_cols) > 0) {
    rowMeans(as.matrix(results_df[, correlation_cols, drop = FALSE]), na.rm = TRUE)
  } else {
    NA_real_
  }

  fit_summary_df <- results_df %>%
    dplyr::group_by(method, split, fold) %>%
    dplyr::summarise(
      lambda = dplyr::first(lambda),
      fit_time_sec = dplyr::first(fit_time_sec),
      .groups = "drop"
    ) %>%
    dplyr::group_by(method, split) %>%
    dplyr::summarise(
      median_lambda = stats::median(lambda, na.rm = TRUE),
      mean_fit_time_sec = mean(fit_time_sec, na.rm = TRUE),
      .groups = "drop"
    )

  results_with_mean_corr %>%
    dplyr::group_by(method, split) %>%
    dplyr::summarise(
      mean_global_mse = mean(global_mse, na.rm = TRUE),
      mean_correlation = mean(mean_correlation, na.rm = TRUE),
      dplyr::across(
        dplyr::all_of(correlation_cols),
        ~ mean(.x, na.rm = TRUE),
        .names = "mean_{.col}"
      ),
      .groups = "drop"
    ) %>%
    dplyr::left_join(fit_summary_df, by = c("method", "split")) %>%
    dplyr::arrange(split, dplyr::desc(mean_correlation), mean_global_mse)
}

cv_fold_rmse_all_na <- function(fit) {
  !is.null(fit$cv_folds) &&
    "rmse" %in% names(fit$cv_folds) &&
    nrow(fit$cv_folds) > 0 &&
    all(is.na(fit$cv_folds$rmse))
}

fit_ccar3_method <- function(method, X_train, Y_train, r, lambda_grid,
                             groups = NULL, Gamma = NULL,
                             inner_kfolds = 20, niter = 10000,
                             rho = 1, thresh = 1e-5,
                             nb_cores = NULL) {
  if (method == "cca_rrr_cv") {
    return(cca_rrr_cv(
      X = X_train,
      Y = Y_train,
      r = r,
      lambdas = lambda_grid,
      kfolds = inner_kfolds,
      standardize = FALSE,
      preprocess="none",
      parallelize = TRUE,
      nb_cores = nb_cores,
      niter = niter,
      rho = rho,
      thresh = thresh,
      LW_Sy = TRUE
    ))
  }
  
  if (method == "cca_graph_rrr_cv") {
    fit <- cca_graph_rrr_cv(
      X = X_train,
      Y = Y_train,
      Gamma = Gamma,
      r = r,
      lambdas = lambda_grid,
      kfolds = inner_kfolds,
      standardize = FALSE,
      preprocess = "none",
      parallelize = TRUE,
      nb_cores = nb_cores,
      niter = niter,
      rho = rho,
      thresh = thresh,
      LW_Sy = TRUE
    )

    if (cv_fold_rmse_all_na(fit)) {
      message("    Parallel cca_graph_rrr_cv produced all-NA CV folds; retrying serially.")
      fit <- cca_graph_rrr_cv(
        X = X_train,
        Y = Y_train,
        Gamma = Gamma,
        r = r,
        lambdas = lambda_grid,
        kfolds = inner_kfolds,
        standardize = FALSE,
        preprocess = "none",
        parallelize = FALSE,
        nb_cores = nb_cores,
        niter = niter,
        rho = rho,
        thresh = thresh,
        LW_Sy = TRUE
      )
    }

    return(fit)
  }

  if (method == "cca_group_rrr_cv") {
    return(cca_group_rrr_cv(
      X = X_train,
      Y = Y_train,
      groups = groups,
      r = r,
      lambdas = lambda_grid,
      kfolds = inner_kfolds,
      standardize = FALSE,
      preprocess = "none",
      parallelize = TRUE,
      nb_cores = nb_cores,
      niter = niter,
      rho = rho,
      thresh = thresh,
      LW_Sy = TRUE
    ))
  }

  stop("Unknown CCAR3 method: ", method, call. = FALSE)
}

fit_baseline_method <- function(method, X_train, Y_train, r, lambda_grid, inner_kfolds) {
  additional_checks(
    X_train = X_train,
    Y_train = Y_train,
    S = NULL,
    rank = r,
    kfolds = inner_kfolds,
    method.type = method,
    lambdax = lambda_grid,
    lambday = c(0)
  )
}

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

output_dir <- normalizePath(file.path(script_dir, "results"),
                            winslash = "/", mustWork = FALSE)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
}

checkpoint_dir <- file.path(output_dir, "fold_method_checkpoints")
if (!dir.exists(checkpoint_dir)) {
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
}
removed_checkpoints <- remove_existing_checkpoints(checkpoint_dir)

results_path <- file.path(output_dir, "brain_method_comparison_results.csv")
summary_path <- file.path(output_dir, "brain_method_comparison_summary.csv")
errors_path <- file.path(output_dir, "brain_method_comparison_errors.csv")
rank_r <- as.integer(as_num(opts$r, 2))
correlation_cols <- paste0("correlation_", seq_len(rank_r))

for (path in c(results_path, summary_path, errors_path)) {
  if (file.exists(path)) {
    file.remove(path)
  }
}

readr::write_csv(
  tibble::as_tibble(c(
    list(
      method = character(),
      lambda = numeric(),
      fit_time_sec = numeric(),
      fold = integer(),
      split = character(),
      global_mse = numeric()
    ),
    stats::setNames(rep(list(numeric()), length(correlation_cols)), correlation_cols)
  )),
  results_path
)
readr::write_csv(
  tibble::as_tibble(c(
    list(
      method = character(),
      split = character(),
      mean_global_mse = numeric(),
      mean_correlation = numeric(),
      median_lambda = numeric(),
      mean_fit_time_sec = numeric()
    ),
    stats::setNames(rep(list(numeric()), length(correlation_cols)), paste0("mean_", correlation_cols))
  )),
  summary_path
)
readr::write_csv(
  tibble::tibble(
    fold = integer(),
    method = character(),
    fit_time_sec = numeric(),
    error = character()
  ),
  errors_path
)

message("Writing detailed results incrementally to: ", results_path)
message("Writing summary checkpoints to: ", summary_path)
message("Writing per-fit checkpoints to: ", checkpoint_dir)
if (removed_checkpoints > 0) {
  message("Cleared ", removed_checkpoints, " old checkpoint CSVs from: ", checkpoint_dir)
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
  opts$group_path,
  Sys.getenv("BRAIN_GROUP_PATH", unset = ""),
  file.path(data_dir, "groups.xlsx"),
  file.path(data_dir, "groups.csv"),
  file.path(project_root, "data", "groups.xlsx"),
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

X <- read_matrix_csv(x_path)
Y <- read_matrix_csv(y_path)

if (nrow(X) != nrow(Y)) {
  stop("X and Y must have the same number of rows.", call. = FALSE)
}

positions <- build_positions_table(coordinates_path, label_path)
if (!is.null(label_path)) {
  message("Using atlas label metadata from: ", label_path)
}
positions <- apply_roi_group_metadata(positions, roi_group_path)
if (!is.null(roi_group_path)) {
  message("Using ROI group metadata from: ", roi_group_path)
}
aligned_X <- align_brain_matrix_to_positions(X, positions, group_path = aggregation_group_path)
X <- aligned_X$X
positions <- aligned_X$positions
if (!is.null(aligned_X$message)) {
  message(aligned_X$message)
}

Gamma <- build_knn_graph_incidence(positions, k = as_num(opts$graph_k, 4))
groups <- build_groups(
  positions
)
group_count <- length(groups)
message("Constructed ", group_count, " anatomical groups for group-regularized CCAR3.")

outer_folds <- as.integer(as_num(opts$outer_folds, 20))
inner_kfolds <- as.integer(as_num(opts$inner_folds, 15))
seed <- as.integer(as_num(opts$seed, 123))
lambda_grid <- 10 ^ seq(
  from = as_num(opts$lambda_min_log10, -3),
  to = as_num(opts$lambda_max_log10, 2),
  length.out = as.integer(as_num(opts$nlambda, 20))
)
niter <- as.integer(as_num(opts$niter, 10000))
rho <- as_num(opts$rho, 1)
thresh <- as_num(opts$thresh, 1e-5)
nb_cores <- if (is.null(opts$nb_cores) || identical(opts$nb_cores, "")) NULL else as.integer(opts$nb_cores)

all_methods <- c(
  "cca_rrr_cv",
  # "cca_graph_rrr_cv",
  # "cca_group_rrr_cv",
  # "FIT_SAR_CV",
  # "FIT_SAR_BIC",
  # "Witten_Perm",
  # "Witten.CV",
  # "Waaijenborg-Author",
  # "Waaijenborg-CV",
  # "SCCA_Parkhomenko",
  # "Chao",
  # "Fantope",
  "SGCA"
)

requested_methods <- opts$methods %||% paste(all_methods, collapse = ",")
methods <- intersect(all_methods, strsplit(requested_methods, ",", fixed = TRUE)[[1]] %>% trimws())
if (length(methods) == 0) {
  stop("No valid methods selected.", call. = FALSE)
}

ccar3_methods <- methods[startsWith(methods, "cca_")]
baseline_methods <- setdiff(methods, ccar3_methods)
available_methods <- character()

if (length(ccar3_methods) > 0) {
  ccar3_dir <- opts$ccar3_dir %||% Sys.getenv(
    "CCAR3_PKG_PATH",
    unset = file.path(path.expand("~"), "Documents", "ccar3", "R")
  )
  ccar3_load_error <- tryCatch(
    {
      source_ccar3_methods(project_root, ccar3_dir = ccar3_dir)
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  if (is.null(ccar3_load_error)) {
    available_methods <- c(available_methods, ccar3_methods)
  } else {
    message("Skipping CCAR3 methods because they could not be loaded: ", ccar3_load_error)
  }
}

if (length(baseline_methods) > 0) {
  benchmark_load_error <- tryCatch(
    {
      source_benchmark_method_stack(project_root)
      NULL
    },
    error = function(e) conditionMessage(e)
  )

  if (is.null(benchmark_load_error)) {
    available_methods <- c(available_methods, baseline_methods)
  } else {
    message("Skipping benchmark methods because they could not be loaded: ", benchmark_load_error)
  }
}

methods <- unique(available_methods)
if (length(methods) == 0) {
  stop(
    "No selected methods are available in this R environment after loading local method code.",
    call. = FALSE
  )
}

if ("cca_group_rrr_cv" %in% methods && group_count != 30) {
  stop(
    "Expected 30 anatomical groups for `cca_group_rrr_cv`, but built ",
    group_count,
    ". Check ROI name parsing and atlas gyrus labels.",
    call. = FALSE
  )
}

set.seed(seed)
folds <- create_cv_folds(seq_len(nrow(X)), k = outer_folds, list = TRUE, returnTrain = FALSE)

results <- list()
errors <- list()
scale_inputs_with_train <- as_bool(opts$standardize_inputs, default = FALSE)

message("Running brain comparison on ", nrow(X), " subjects, ", ncol(X), " ROIs, ", ncol(Y), " outcomes.")
message("Methods: ", paste(methods, collapse = ", "))
if (scale_inputs_with_train) {
  message("Applying fold-wise train centering and scaling to X/Y, with test splits transformed using train statistics.")
} else {
  message("Applying fold-wise train centering to X/Y, with test splits transformed using train means.")
}
X = scale(X)
Y = scale(Y)
for (fold_id in seq_along(folds)) {
  test_idx <- sort(as.integer(folds[[fold_id]]))
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)

  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx, , drop = FALSE]
  X_test <- X[test_idx, , drop = FALSE]
  Y_test <- Y[test_idx, , drop = FALSE]

  #x_fold <- prepare_train_test_matrix(X_train, X_test, scale = scale_inputs_with_train)
  #y_fold <- prepare_train_test_matrix(Y_train, Y_test, scale = scale_inputs_with_train)
  #X_train <- x_fold$train
  #X_test <- x_fold$test
  #Y_train <- y_fold$train
  #Y_test <- y_fold$test
  # X_train <- scale(X_train, scale=FALSE)
  # Y_train <- scale(Y_train, scale=FALSE)
  # X_test <- scale(X_test, scale=FALSE)
  # Y_test <- scale(Y_test, scale=FALSE)

  message(sprintf("Outer fold %d/%d", fold_id, length(folds)))

  for (method in methods) {
    message("  Fitting ", method)

    start_time <- proc.time()[["elapsed"]]
    error_message <- NULL
    fit <- tryCatch(
      {
        if (startsWith(method, "cca_")) {
          fit_ccar3_method(
            method = method,
            X_train = X_train,
            Y_train = Y_train,
            r = rank_r,
            lambda_grid = lambda_grid,
            groups = groups,
            Gamma = Gamma,
            inner_kfolds = inner_kfolds,
            niter = niter,
            rho = rho,
            thresh = thresh,
            nb_cores = nb_cores
          )
        } else {
          fit_baseline_method(
            method = method,
            X_train = X_train,
            Y_train = Y_train,
            r = rank_r,
            lambda_grid = lambda_grid,
            inner_kfolds = inner_kfolds
          )
        }
      },
      error = function(e) {
        error_message <<- conditionMessage(e)
        NULL
      }
    )
    fit_time_sec <- proc.time()[["elapsed"]] - start_time

    if (is.null(fit)) {
      error_row <- tibble::tibble(
        fold = fold_id,
        method = method,
        fit_time_sec = fit_time_sec,
        error = error_message %||% "Unknown fitting error."
      )
      errors[[length(errors) + 1]] <- error_row
      append_csv_rows(error_row, errors_path)
      error_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "error")
      result_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "results")
      cv_summary_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "cv_summary")
      cv_folds_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "cv_folds")
      readr::write_csv(
        error_row,
        error_checkpoint_path
      )
      remove_if_exists(c(result_checkpoint_path, cv_summary_checkpoint_path, cv_folds_checkpoint_path))
      message("    Wrote error checkpoint: ", error_checkpoint_path)
      next
    }

    U <- fit$U %||% fit$u
    V <- fit$V %||% fit$v
    lambda_value <- fit$lambda %||% NA_real_
    standardized_fit <- if (startsWith(method, "cca_")) {
      inspect_train_scores(X_train, Y_train, U, V)
    } else {
      standardize_train_scores(X_train, Y_train, U, V)
    }
    U <- standardized_fit$U
    V <- standardized_fit$V

    if (!is.na(standardized_fit$x_jitter) && standardized_fit$x_jitter > 0) {
      message("    Added X-score whitening jitter: ", format(standardized_fit$x_jitter, scientific = TRUE))
    }
    if (!is.na(standardized_fit$y_jitter) && standardized_fit$y_jitter > 0) {
      message("    Added Y-score whitening jitter: ", format(standardized_fit$y_jitter, scientific = TRUE))
    }
    message("    crossprod(X_train %*% U) / nrow(X_train):")
    print(signif(standardized_fit$x_crossprod, 6))
    message("    crossprod(Y_train %*% V) / nrow(Y_train):")
    print(signif(standardized_fit$y_crossprod, 6))

    method_results <- dplyr::bind_rows(
      score_pair(X_train, Y_train, U, V, "train", fold_id, method, lambda_value, fit_time_sec,
                 n_components = rank_r),
      score_pair(X_test, Y_test, U, V, "test", fold_id, method, lambda_value, fit_time_sec,
                 n_components = rank_r)
    )
    results[[length(results) + 1]] <- method_results

    append_csv_rows(method_results, results_path)
    error_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "error")
    result_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "results")
    readr::write_csv(
      method_results,
      result_checkpoint_path
    )
    remove_if_exists(error_checkpoint_path)

    if (!is.null(fit$cv_summary)) {
      cv_summary_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "cv_summary")
      readr::write_csv(as.data.frame(fit$cv_summary), cv_summary_checkpoint_path)
      message("    Wrote CV summary checkpoint: ", cv_summary_checkpoint_path)
    }
    if (!is.null(fit$cv_folds)) {
      cv_folds_checkpoint_path <- checkpoint_result_path(checkpoint_dir, fold_id, method, suffix = "cv_folds")
      readr::write_csv(as.data.frame(fit$cv_folds), cv_folds_checkpoint_path)
      message("    Wrote CV folds checkpoint: ", cv_folds_checkpoint_path)
    }

    current_summary_df <- build_summary_df(dplyr::bind_rows(results))
    readr::write_csv(current_summary_df, summary_path)
    message("    Wrote results checkpoint: ", result_checkpoint_path)
    message("    Updated aggregate CSVs: ", results_path, " and ", summary_path)
  }
}

results_df <- dplyr::bind_rows(results)
errors_df <- dplyr::bind_rows(errors)

if (nrow(results_df) == 0) {
  stop("All method fits failed; no results were produced.", call. = FALSE)
}

summary_df <- build_summary_df(results_df)

readr::write_csv(results_df, results_path)
readr::write_csv(summary_df, summary_path)
if (nrow(errors_df) > 0) {
  readr::write_csv(errors_df, errors_path)
}

message("Saved detailed results to: ", results_path)
message("Saved summary results to: ", summary_path)
message("Saved per-fit checkpoints to: ", checkpoint_dir)
if (nrow(errors_df) > 0) {
  message("Saved fitting errors to: ", errors_path)
}
