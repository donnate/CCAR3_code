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

canonical_axis_label <- function(score_prefix, component_idx) {
  as.expression(bquote(.(as.name(score_prefix))[.(component_idx)]))
}

plot_embedding_pair <- function(df, x_col, y_col, out_file, method_name, legend_order, colors, labels_n) {
  if (!all(c(x_col, y_col) %in% names(df))) {
    return(invisible(NULL))
  }

  x_idx <- as.integer(sub("^XU", "", x_col))
  y_idx <- as.integer(sub("^XU", "", y_col))

  p <- ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]], colour = diet)) +
    geom_point(aes(shape = genotype), size = 4) +
    scale_color_manual(values = colors, breaks = legend_order, labels = labels_n) +
    scale_fill_manual(values = colors, breaks = legend_order, labels = labels_n) +
    stat_ellipse(level = 0.95) +
    xlab(canonical_axis_label("Xu", x_idx)) +
    ylab(canonical_axis_label("Xu", y_idx)) +
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
    labels_n
  )
  plot_embedding_pair(
    plot_df,
    "XU3",
    "XU4",
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-var2.pdf")),
    display_name,
    legend_order,
    colors,
    labels_n
  )
  plot_embedding_pair(
    plot_df,
    "XU1",
    "XU5",
    file.path(context$plot_dir, paste0("nutrimouse-", file_stub, "-var3.pdf")),
    display_name,
    legend_order,
    colors,
    labels_n
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
    full_fit_summary = saved$summary
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
