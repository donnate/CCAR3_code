#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

`%||%` <- function(x, y) {
  if (is.null(x)) {
    return(y)
  }
  x
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

env_or_default <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (nzchar(value)) {
    return(value)
  }
  default
}

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

safe_mean <- function(x) {
  out <- mean(x, na.rm = TRUE)
  if (is.nan(out)) {
    return(NA_real_)
  }
  out
}

make_folds <- function(indices, k, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  shuffled <- sample(indices, length(indices), replace = FALSE)
  split(shuffled, cut(seq_along(shuffled), breaks = k, labels = FALSE))
}

diag_cor_values <- function(A, B) {
  cor_mat <- suppressWarnings(stats::cor(A, B))
  if (is.null(dim(cor_mat))) {
    return(abs(as.numeric(cor_mat)))
  }
  diag(cor_mat)
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
    metadata = metadata,
    gene_names = colnames(X_raw),
    lipid_names = colnames(Y_raw)
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

resolve_test_correlation_dir <- function(benchmark_dir) {
  override_dir <- env_or_default("CCAR3_NUTRIMOUSE_TEST_COR_DIR", "")
  candidate_dirs <- c(
    if (nzchar(override_dir)) override_dir else character(0),
    benchmark_dir,
    paste0(benchmark_dir, "_cv"),
    paste0(benchmark_dir, "_1"),
    paste0(benchmark_dir, "_scaled_with_fantope")
  )

  for (candidate in unique(candidate_dirs)) {
    if (!nzchar(candidate)) {
      next
    }
    if (file.exists(file.path(candidate, "benchmark_runs.csv")) &&
        dir.exists(file.path(candidate, "benchmark_details"))) {
      return(candidate)
    }
  }

  NULL
}

load_test_component_correlations <- function(benchmark_dir, X_raw, Y_raw, methods = NULL, seed = 123L) {
  if (is.null(benchmark_dir)) {
    return(NULL)
  }

  benchmark_runs_path <- file.path(benchmark_dir, "benchmark_runs.csv")
  benchmark_details_dir <- file.path(benchmark_dir, "benchmark_details")
  if (!file.exists(benchmark_runs_path) || !dir.exists(benchmark_details_dir)) {
    return(NULL)
  }

  benchmark_runs <- utils::read.csv(benchmark_runs_path, stringsAsFactors = FALSE)
  if (!is.null(methods)) {
    benchmark_runs <- benchmark_runs[benchmark_runs$method %in% methods, , drop = FALSE]
  }
  if ("timed_out" %in% names(benchmark_runs)) {
    timed_out <- benchmark_runs$timed_out
    if (is.character(timed_out)) {
      timed_out <- tolower(trimws(timed_out)) %in% c("true", "t", "1")
    }
    benchmark_runs <- benchmark_runs[!timed_out, , drop = FALSE]
  }

  outer_folds <- length(unique(stats::na.omit(benchmark_runs$fold_id)))
  if (outer_folds == 0 || nrow(benchmark_runs) == 0) {
    return(NULL)
  }

  fold_cache <- new.env(parent = emptyenv())
  get_repeat_folds <- function(repeat_id) {
    key <- as.character(repeat_id)
    if (!exists(key, envir = fold_cache, inherits = FALSE)) {
      assign(
        key,
        make_folds(seq_len(nrow(X_raw)), k = outer_folds, seed = seed + repeat_id),
        envir = fold_cache
      )
    }
    get(key, envir = fold_cache, inherits = FALSE)
  }

  component_rows <- lapply(seq_len(nrow(benchmark_runs)), function(i) {
    run_row <- benchmark_runs[i, , drop = FALSE]
    detail_path <- file.path(
      benchmark_details_dir,
      sprintf("repeat_%02d_fold_%02d_%s.rds", run_row$repeat_id, run_row$fold_id, run_row$method)
    )
    if (!file.exists(detail_path)) {
      return(NULL)
    }

    artifact <- readRDS(detail_path)
    U <- artifact$result$U %||% NULL
    V <- artifact$result$V %||% NULL
    if (is.null(U) || is.null(V)) {
      return(NULL)
    }

    folds <- get_repeat_folds(artifact$repeat_id)
    test_idx <- folds[[artifact$fold_id]]
    if (length(test_idx) == 0) {
      return(NULL)
    }

    component_cor <- diag_cor_values(
      X_raw[test_idx, , drop = FALSE] %*% as.matrix(U),
      Y_raw[test_idx, , drop = FALSE] %*% as.matrix(V)
    )
    component_cor <- component_cor[seq_len(min(5L, length(component_cor)))]
    names(component_cor) <- paste0("cor_", seq_along(component_cor))

    data.frame(
      method = artifact$method,
      as.list(component_cor),
      stringsAsFactors = FALSE
    )
  })

  component_df <- dplyr::bind_rows(component_rows)
  if (nrow(component_df) == 0) {
    return(NULL)
  }

  component_names <- intersect(paste0("cor_", 1:5), names(component_df))
  if (length(component_names) == 0) {
    return(NULL)
  }

  summary_df <- component_df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(component_names), safe_mean),
      .groups = "drop"
    )

  summary_df
}

sanitize_file_stub <- function(x) {
  gsub("_+", "_", gsub("[^[:alnum:]]+", "_", tolower(x)))
}

method_label <- function(method) {
  labels <- c(
    cca_rrr = "CCA-RRR",
    witten_cv = "Witten CV",
    witten_pma = "Witten PMA",
    sar = "SAR",
    parkhomenko = "Parkhomenko",
    waaijenborg = "Waaijenborg",
    scca_gao = "SCCA Gao",
    gca = "GCA",
    fantope = "Fantope"
  )

  labels[[method]] %||% method
}

load_saved_full_fits <- function(out_dir) {
  full_fit_path <- file.path(out_dir, "full_data_fits.rds")
  diet_task_path <- file.path(out_dir, "full_model_diet_tasks.rds")
  summary_path <- file.path(out_dir, "full_data_fit_summary.csv")

  summary_df <- NULL
  if (file.exists(summary_path)) {
    summary_df <- utils::read.csv(summary_path, stringsAsFactors = FALSE)
  }

  if (file.exists(diet_task_path)) {
    bundle <- readRDS(diet_task_path)
    artifacts <- bundle$artifacts
  } else if (file.exists(full_fit_path)) {
    bundle <- readRDS(full_fit_path)
    artifacts <- bundle$full_data_fits %||% bundle$artifacts
  } else {
    stop(
      "Could not find saved full-data fits in ", out_dir,
      ". Expected full_model_diet_tasks.rds or full_data_fits.rds.",
      call. = FALSE
    )
  }

  # If the richer full-data artifact bundle exists, use it to fill in any
  # missing fields such as V loadings while keeping full_model_diet_tasks.rds
  # as the canonical source of saved methods.
  if (file.exists(full_fit_path)) {
    full_fit_bundle <- readRDS(full_fit_path)
    full_fit_artifacts <- full_fit_bundle$full_data_fits %||% full_fit_bundle$artifacts

    if (!is.null(full_fit_artifacts)) {
      for (method in intersect(names(artifacts), names(full_fit_artifacts))) {
        for (field in c("U", "V", "x_scores", "y_scores", "component_cor", "metadata", "fit")) {
          if (is.null(artifacts[[method]][[field]]) && !is.null(full_fit_artifacts[[method]][[field]])) {
            artifacts[[method]][[field]] <- full_fit_artifacts[[method]][[field]]
          }
        }
      }
    }
  }

  if (is.null(artifacts) || length(artifacts) == 0) {
    stop("Saved fit artifacts were found, but no methods were stored.", call. = FALSE)
  }

  artifact_names <- names(artifacts)
  for (i in seq_along(artifacts)) {
    artifacts[[i]]$method <- artifacts[[i]]$method %||% artifact_names[[i]]
  }

  if (!is.null(summary_df)) {
    ordered_methods <- summary_df$method[summary_df$method %in% names(artifacts)]
    remaining_methods <- setdiff(names(artifacts), ordered_methods)
    artifacts <- c(artifacts[ordered_methods], artifacts[remaining_methods])
  }

  list(
    artifacts = artifacts,
    summary = summary_df
  )
}

prepare_score_matrix <- function(mat, prefix, max_components = 5L) {
  if (is.null(mat)) {
    return(NULL)
  }

  mat <- as.matrix(mat)
  keep_cols <- seq_len(min(max_components, ncol(mat)))
  mat <- mat[, keep_cols, drop = FALSE]
  colnames(mat) <- paste0(prefix, keep_cols)
  mat
}

build_projection_df <- function(U, V, X, Y, metadata, x_scores = NULL, y_scores = NULL, max_components = 5L) {
  if (is.null(x_scores)) {
    U <- as.matrix(U)
    keep_u <- seq_len(min(max_components, ncol(U)))
    x_scores <- as.matrix(X %*% U[, keep_u, drop = FALSE])
  }
  x_scores <- prepare_score_matrix(x_scores, "XU", max_components = max_components)

  df <- data.frame(
    x_scores,
    genotype = metadata$genotype,
    diet = metadata$diet,
    stringsAsFactors = FALSE
  )

  if (is.null(y_scores) && !is.null(V)) {
    V <- as.matrix(V)
    keep_v <- seq_len(min(max_components, ncol(V)))
    y_scores <- as.matrix(Y %*% V[, keep_v, drop = FALSE])
  }
  y_scores <- prepare_score_matrix(y_scores, "YV", max_components = max_components)

  if (!is.null(y_scores)) {
    df <- cbind(df, as.data.frame(y_scores, stringsAsFactors = FALSE))
  }

  df
}

canonical_axis_label <- function(score_prefix, component_idx, cor_value = NULL, cor_source = NULL) {
  axis_name <- as.name(score_prefix)
  if (is.null(cor_value) || length(cor_value) == 0 || !is.finite(cor_value)) {
    return(as.expression(bquote(.(axis_name)[.(component_idx)])))
  }
  rounded_cor <- round(cor_value, 2)

  if (identical(cor_source, "test")) {
    return(as.expression(bquote(atop(.(axis_name)[.(component_idx)], hat(rho)[test] == .(rounded_cor)))))
  }

  if (identical(cor_source, "full")) {
    return(as.expression(bquote(atop(.(axis_name)[.(component_idx)], hat(rho)[full] == .(rounded_cor)))))
  }

  as.expression(bquote(atop(.(axis_name)[.(component_idx)], hat(rho) == .(rounded_cor))))
}

extract_method_component_correlations <- function(method, context) {
  component_names <- paste0("cor_", 1:5)

  if (!is.null(context$test_component_cor)) {
    test_row <- context$test_component_cor[context$test_component_cor$method == method, , drop = FALSE]
    available_names <- intersect(component_names, names(test_row))
    if (nrow(test_row) > 0 && length(available_names) > 0) {
      values <- rep(NA_real_, length(component_names))
      names(values) <- component_names
      values[available_names] <- as.numeric(test_row[1, available_names, drop = TRUE])
      return(list(values = values, source = "test"))
    }
  }

  if (!is.null(context$full_fit_summary)) {
    full_row <- context$full_fit_summary[context$full_fit_summary$method == method, , drop = FALSE]
    available_names <- intersect(component_names, names(full_row))
    if (nrow(full_row) > 0 && length(available_names) > 0) {
      values <- rep(NA_real_, length(component_names))
      names(values) <- component_names
      values[available_names] <- as.numeric(full_row[1, available_names, drop = TRUE])
      return(list(values = values, source = "full"))
    }
  }

  list(values = stats::setNames(rep(NA_real_, length(component_names)), component_names), source = NULL)
}

plot_embedding_pair <- function(df, x_col, y_col, out_file, method_name, legend_order, colors, labels_n,
                                component_cor = NULL, cor_source = NULL) {
  if (!all(c(x_col, y_col) %in% names(df))) {
    return(invisible(NULL))
  }

  x_idx <- as.integer(sub("^XU", "", x_col))
  y_idx <- as.integer(sub("^XU", "", y_col))
  x_cor <- component_cor[[paste0("cor_", x_idx)]] %||% NA_real_
  y_cor <- component_cor[[paste0("cor_", y_idx)]] %||% NA_real_

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = diet)) +
    geom_point(aes(shape = genotype), size = 4) +
    scale_color_manual(values = colors, breaks = legend_order, labels = labels_n) +
    scale_fill_manual(values = colors, breaks = legend_order, labels = labels_n) +
    stat_ellipse(level = 0.95) +
    xlab(canonical_axis_label("Xu", x_idx, cor_value = x_cor, cor_source = cor_source)) +
    ylab(canonical_axis_label("Xu", y_idx, cor_value = y_cor, cor_source = cor_source)) +
    labs(title = method_name, colour = "Diet", shape = "Genotype") +
    theme(legend.position = "none")

  ggsave(out_file, plot = p, width = 4, height = 4)
  invisible(p)
}

plot_arrow_overlay <- function(df, out_file, method_name) {
  required_cols <- c("XU1", "XU2", "YV1", "YV2")
  if (!all(required_cols %in% names(df))) {
    return(invisible(NULL))
  }

  p <- ggplot(df) +
    geom_point(aes(x = XU1, y = XU2), shape = 1, size = 2) +
    geom_point(aes(x = YV1, y = YV2), shape = 1, size = 2) +
    geom_segment(
      aes(x = XU1, y = XU2, xend = YV1, yend = YV2),
      arrow = arrow(length = grid::unit(0.3, "cm"))
    ) +
    xlab("Canonical variate 1") +
    ylab("Canonical variate 2") +
    labs(title = method_name) +
    theme_minimal(base_size = 11)

  ggsave(out_file, plot = p, width = 4, height = 4)
  invisible(p)
}

plot_loading_bars <- function(loadings, feature_names, out_file, facet_symbol, y_label, width, height, threshold = 0.1) {
  if (is.null(loadings)) {
    return(invisible(NULL))
  }

  loadings <- as.matrix(loadings)
  keep_cols <- seq_len(min(5L, ncol(loadings)))
  loadings <- loadings[, keep_cols, drop = FALSE]
  colnames(loadings) <- paste0(facet_symbol, keep_cols)

  plot_df <- as.data.frame(loadings, stringsAsFactors = FALSE)
  plot_df$name <- feature_names
  plot_df <- plot_df %>%
    pivot_longer(cols = -name, names_to = "component", values_to = "value") %>%
    filter(abs(value) > threshold)

  if (nrow(plot_df) == 0) {
    return(invisible(NULL))
  }

  plot_df <- plot_df %>%
    mutate(
      component_label = factor(
        component,
        levels = paste0(facet_symbol, keep_cols),
        labels = paste0(tolower(facet_symbol), "[", keep_cols, "]")
      )
    )

  p <- ggplot(plot_df) +
    geom_col(aes(x = value, y = name, fill = component_label), show.legend = FALSE) +
    facet_grid(cols = vars(component_label), labeller = label_parsed, scales = "free_x", space = "free_x") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    labs(y = y_label, x = "Value") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x = element_text(size = 11),
      axis.text.y = element_text(size = 11)
    )

  ggsave(out_file, plot = p, width = width, height = height)
  invisible(p)
}

plot_saved_method <- function(method, context = plot_context) {
  if (missing(context) || is.null(context)) {
    stop("No plotting context has been loaded yet.", call. = FALSE)
  }

  if (!(method %in% names(context$full_fit_artifacts))) {
    stop(
      "Unknown method '", method, "'. Available methods: ",
      paste(names(context$full_fit_artifacts), collapse = ", "),
      call. = FALSE
    )
  }

  fit <- context$full_fit_artifacts[[method]]
  if (is.null(fit$U)) {
    warning("Skipping ", method, ": missing saved U loadings.", call. = FALSE)
    return(invisible(NULL))
  }

  file_stub <- sanitize_file_stub(method)
  display_name <- method_label(method)
  method_component_cor <- extract_method_component_correlations(method, context)
  saved_x_scores <- fit$x_scores %||% NULL
  if (is.null(saved_x_scores) && !is.null(fit$diet$embedding)) {
    saved_x_scores <- sqrt(nrow(context$X)) * fit$diet$embedding
  }
  plot_df <- build_projection_df(
    U = fit$U,
    V = fit$V %||% NULL,
    X = context$X,
    Y = context$Y,
    metadata = context$metadata,
    x_scores = saved_x_scores,
    y_scores = fit$y_scores %||% NULL
  )

  legend_order <- c("lin", "sun", "fish", "ref", "coc")
  colors <- c("red", "orange", "dodgerblue", "black", "brown")
  labels_n <- c("LIN", "SUN", "FISH", "REF", "COC")

  if (!all(c("YV1", "YV2") %in% names(plot_df))) {
    message("Skipping XU/YV overlay for ", method, ": no saved Y-side scores or V loadings were found.")
  }

  plot_arrow_overlay(
    plot_df,
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-arrow.pdf")),
    display_name
  )
  plot_embedding_pair(
    plot_df,
    "XU1",
    "XU2",
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-var1.pdf")),
    display_name,
    legend_order,
    colors,
    labels_n,
    component_cor = method_component_cor$values,
    cor_source = method_component_cor$source
  )
  plot_embedding_pair(
    plot_df,
    "XU3",
    "XU4",
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-var2.pdf")),
    display_name,
    legend_order,
    colors,
    labels_n,
    component_cor = method_component_cor$values,
    cor_source = method_component_cor$source
  )
  plot_embedding_pair(
    plot_df,
    "XU1",
    "XU5",
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-var3.pdf")),
    display_name,
    legend_order,
    colors,
    labels_n,
    component_cor = method_component_cor$values,
    cor_source = method_component_cor$source
  )
  plot_loading_bars(
    fit$U,
    context$gene_names,
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-u.pdf")),
    facet_symbol = "U",
    y_label = "Genes",
    width = 10,
    height = 6.5
  )
  plot_loading_bars(
    fit$V %||% NULL,
    context$lipid_names,
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-v.pdf")),
    facet_symbol = "V",
    y_label = "Hepatic Fatty Acids",
    width = 10,
    height = 6
  )

  invisible(plot_df)
}

plot_saved_methods <- function(methods = names(plot_context$full_fit_artifacts), context = plot_context) {
  invisible(lapply(methods, plot_saved_method, context = context))
}

resolve_requested_methods <- function(available_methods) {
  requested <- env_or_default("CCAR3_NUTRIMOUSE_PLOT_METHODS", "")
  if (!nzchar(requested)) {
    return(available_methods)
  }

  methods <- trimws(strsplit(requested, ",", fixed = TRUE)[[1]])
  methods <- methods[nzchar(methods)]
  missing_methods <- setdiff(methods, available_methods)
  if (length(missing_methods) > 0) {
    warning(
      "Skipping unknown methods from CCAR3_NUTRIMOUSE_PLOT_METHODS: ",
      paste(missing_methods, collapse = ", "),
      call. = FALSE
    )
  }

  intersect(methods, available_methods)
}

main <- function() {
  project_root <- find_project_root()
  setwd(project_root)

  data_obj <- load_nutrimouse_data()
  full_data <- prepare_full_nutrimouse_data(data_obj)
  benchmark_dir <- env_or_default(
    "CCAR3_NUTRIMOUSE_BENCHMARK_DIR",
    file.path(
      project_root,
      "experiments",
      "real-data-experiments",
      "nutrimouse",
      "nutrimouse_benchmark_all_methods"
    )
  )
  plot_dir <- env_or_default(
    "CCAR3_NUTRIMOUSE_PLOT_DIR",
    file.path(benchmark_dir, "plots")
  )
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

  saved <- load_saved_full_fits(benchmark_dir)
  requested_methods <- resolve_requested_methods(names(saved$artifacts))
  test_cor_dir <- resolve_test_correlation_dir(benchmark_dir)
  test_component_cor <- load_test_component_correlations(
    benchmark_dir = test_cor_dir,
    X_raw = full_data$X_raw,
    Y_raw = full_data$Y_raw,
    methods = requested_methods,
    seed = as.integer(env_or_default("CCAR3_NUTRIMOUSE_SEED", "123"))
  )

  context <- list(
    project_root = project_root,
    benchmark_dir = benchmark_dir,
    plot_dir = plot_dir,
    X = full_data$X_full,
    Y = full_data$Y_full,
    metadata = full_data$metadata,
    gene_names = data_obj$gene_names,
    lipid_names = data_obj$lipid_names,
    full_fit_artifacts = saved$artifacts,
    full_fit_summary = saved$summary,
    test_component_cor = test_component_cor,
    test_cor_dir = test_cor_dir
  )

  for (method in requested_methods) {
    plot_saved_method(method, context = context)
  }

  message("Loaded saved nutrimouse full-data fits for: ", paste(names(context$full_fit_artifacts), collapse = ", "))
  message("Plots written to: ", context$plot_dir)

  context
}

plot_context <- main()
full_fit_artifacts <- plot_context$full_fit_artifacts
full_fit_summary <- plot_context$full_fit_summary
X <- plot_context$X
Y <- plot_context$Y
metadata <- plot_context$metadata
