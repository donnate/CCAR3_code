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

source_ccar3_package_source <- function(ccar3_dir) {
  replacements <- list(
    list(pattern = "^library\\((foreach|caret|pracma)\\).*", replacement = ""),
    list(pattern = "caret::createFolds", replacement = "create_cv_folds"),
    list(pattern = "corpcor::cov\\.shrink\\(Y(?:,\\s*verbose\\s*=\\s*verbose)?\\)", replacement = "stats::cov(Y)")
  )

  source_filtered_r_file(file.path(ccar3_dir, "R", "helpers.r"))
  source_filtered_r_file(file.path(ccar3_dir, "R", "utils.R"))
  source_filtered_r_file(
    file.path(ccar3_dir, "R", "reduced_rank_regression.R"),
    replacements = replacements
  )
  source_filtered_r_file(
    file.path(ccar3_dir, "R", "group_reduced_rank_regression.R"),
    replacements = replacements
  )
  source_filtered_r_file(
    file.path(ccar3_dir, "R", "graph_reduced_rank_regression.R"),
    replacements = replacements
  )
}

source_ccar3_methods <- function(project_root, ccar3_dir = NULL) {
  local_paths <- file.path(
    project_root,
    "src",
    c(
      "reduced_rank_regression.R",
      "group_reduced_rank_regression.R",
      "graph_reduced_rank_regression.R"
    )
  )

  if (all(file.exists(local_paths)) && requireNamespace("ccar3", quietly = TRUE)) {
    source_existing_files(local_paths)
    return(invisible("installed-package"))
  }

  package_candidates <- c(
    ccar3_dir,
    Sys.getenv("CCAR3_PKG_PATH", unset = ""),
    file.path(dirname(project_root), "ccar3"),
    file.path(dirname(project_root), "CCAR3")
  )
  package_candidates <- unique(package_candidates[nzchar(package_candidates)])

  for (candidate in package_candidates) {
    if (!dir.exists(candidate) || !file.exists(file.path(candidate, "DESCRIPTION"))) {
      next
    }

    candidate <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
    options(ccar3_pkg_path = candidate)
    Sys.setenv(CCAR3_PKG_PATH = candidate)
    source_ccar3_package_source(candidate)
    return(invisible("local-package-source"))
  }

  stop(
    "Could not load CCAR3 core methods. Expected either an installed `ccar3` package plus local src/ wrappers, or a local ccar3 package checkout via --ccar3_dir / CCAR3_PKG_PATH.",
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

standardize_with_mean_impute <- function(mat) {
  mat <- as.matrix(mat)
  col_means <- colMeans(mat, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) {
    missing <- is.na(mat[, j])
    if (any(missing)) {
      mat[missing, j] <- col_means[j]
    }
  }
  scale(mat)
}

infer_positions_name <- function(positions) {
  if ("name" %in% names(positions)) {
    return(positions$name)
  }
  if ("Region" %in% names(positions)) {
    return(vapply(strsplit(positions$Region, "_"), function(x) tail(x, 1), character(1)))
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

score_pair <- function(X, Y, U, V, split, fold_id, method, lambda) {
  XU <- as.matrix(X) %*% U
  YV <- as.matrix(Y) %*% V
  tibble::tibble(
    method = method,
    lambda = lambda,
    fold = fold_id,
    split = split,
    component = seq_len(ncol(U)),
    covariance = diag(stats::cov(XU, YV)),
    mse = colMeans((XU - YV) ^ 2),
    correlation = diag(stats::cor(XU, YV))
  )
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
      parallelize = FALSE,
      nb_cores = nb_cores,
      niter = niter,
      rho = rho,
      thresh = thresh,
      LW_Sy = TRUE
    ))
  }

  if (method == "cca_graph_rrr_cv") {
    return(cca_graph_rrr_cv(
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
    ))
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
      parallelize = FALSE,
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

group_candidates <- c(
  opts$group_path,
  Sys.getenv("BRAIN_GROUP_PATH", unset = ""),
  file.path(data_dir, "activation_groups.xlsx"),
  file.path(data_dir, "activation_groups.csv"),
  file.path(project_root, "data", "activation_groups.xlsx"),
  file.path(project_root, "data", "activation_groups.csv")
)
group_path <- group_candidates[file.exists(group_candidates)][1]
if (is.na(group_path)) {
  group_path <- NULL
}

X <- read_matrix_csv(x_path)
Y <- read_matrix_csv(y_path)

if (nrow(X) != nrow(Y)) {
  stop("X and Y must have the same number of rows.", call. = FALSE)
}

if (as_bool(opts$standardize_inputs, default = FALSE)) {
  X <- standardize_with_mean_impute(X)
  Y <- standardize_with_mean_impute(Y)
}

positions <- build_positions_table(coordinates_path, label_path)
aligned_X <- align_brain_matrix_to_positions(X, positions, group_path = group_path)
X <- aligned_X$X
positions <- aligned_X$positions
if (!is.null(aligned_X$message)) {
  message(aligned_X$message)
}

Gamma <- build_knn_graph_incidence(positions, k = as_num(opts$graph_k, 4))
groups <- build_groups(
  positions,
  group_path = if (isTRUE(aligned_X$aggregated)) NULL else group_path,
  expected_n_features = ncol(X)
)

rank_r <- as.integer(as_num(opts$r, 3))
outer_folds <- as.integer(as_num(opts$outer_folds, 20))
inner_kfolds <- as.integer(as_num(opts$inner_folds, outer_folds - 1))
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
  "cca_graph_rrr_cv",
  "cca_group_rrr_cv",
  "FIT_SAR_CV",
  "FIT_SAR_BIC",
  "Witten_Perm",
  "Witten.CV",
  "Waaijenborg-Author",
  "Waaijenborg-CV",
  "SCCA_Parkhomenko",
  "Chao",
  "Fantope",
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
  ccar3_dir <- opts$ccar3_dir %||% Sys.getenv("CCAR3_PKG_PATH", unset = "")
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

set.seed(seed)
folds <- create_cv_folds(seq_len(nrow(X)), k = outer_folds, list = TRUE, returnTrain = FALSE)

results <- list()
errors <- list()

message("Running brain comparison on ", nrow(X), " subjects, ", ncol(X), " ROIs, ", ncol(Y), " outcomes.")
message("Methods: ", paste(methods, collapse = ", "))

for (fold_id in seq_along(folds)) {
  test_idx <- sort(as.integer(folds[[fold_id]]))
  train_idx <- setdiff(seq_len(nrow(X)), test_idx)

  X_train <- X[train_idx, , drop = FALSE]
  Y_train <- Y[train_idx, , drop = FALSE]
  X_test <- X[test_idx, , drop = FALSE]
  Y_test <- Y[test_idx, , drop = FALSE]

  message(sprintf("Outer fold %d/%d", fold_id, length(folds)))

  for (method in methods) {
    message("  Fitting ", method)

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
        errors[[length(errors) + 1]] <<- tibble::tibble(
          fold = fold_id,
          method = method,
          error = conditionMessage(e)
        )
        NULL
      }
    )

    if (is.null(fit)) {
      next
    }

    U <- fit$U %||% fit$u
    V <- fit$V %||% fit$v
    lambda_value <- fit$lambda %||% NA_real_

    results[[length(results) + 1]] <- dplyr::bind_rows(
      score_pair(X_train, Y_train, U, V, "train", fold_id, method, lambda_value),
      score_pair(X_test, Y_test, U, V, "test", fold_id, method, lambda_value)
    )
  }
}

results_df <- dplyr::bind_rows(results)
errors_df <- dplyr::bind_rows(errors)

if (nrow(results_df) == 0) {
  stop("All method fits failed; no results were produced.", call. = FALSE)
}

summary_df <- results_df %>%
  dplyr::group_by(method, split) %>%
  dplyr::summarise(
    mean_covariance = mean(covariance, na.rm = TRUE),
    mean_correlation = mean(correlation, na.rm = TRUE),
    mean_mse = mean(mse, na.rm = TRUE),
    median_lambda = median(lambda, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(split, dplyr::desc(mean_correlation), mean_mse)

results_path <- file.path(output_dir, "brain_method_comparison_results.csv")
summary_path <- file.path(output_dir, "brain_method_comparison_summary.csv")
errors_path <- file.path(output_dir, "brain_method_comparison_errors.csv")

readr::write_csv(results_df, results_path)
readr::write_csv(summary_df, summary_path)
if (nrow(errors_df) > 0) {
  readr::write_csv(errors_df, errors_path)
}

message("Saved detailed results to: ", results_path)
message("Saved summary results to: ", summary_path)
if (nrow(errors_df) > 0) {
  message("Saved fitting errors to: ", errors_path)
}
