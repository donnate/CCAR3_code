get_script_path <- function() {
  file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(file_arg) == 1L) {
    return(normalizePath(sub("^--file=", "", file_arg), winslash = "/", mustWork = TRUE))
  }

  if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile, winslash = "/", mustWork = TRUE))
  }

  stop("Could not determine the script path.", call. = FALSE)
}

repo_root <- "~/Documents/CCAR3_code"

write_results_csv <- function(df, relative_path) {
  out_path <- file.path(repo_root, relative_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, out_path, row.names = FALSE)
  out_path
}

append_and_flush_results <- function(results, row, relative_path) {
  updated_results <- rbind(results, row)
  write_results_csv(updated_results, relative_path)
  updated_results
}

write_plot_file <- function(plot_obj, relative_path, width = 11, height = 7) {
  out_path <- file.path(repo_root, relative_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(out_path, plot = plot_obj, width = width, height = height, dpi = 300)
  out_path
}

elapsed_seconds <- function(timing) {
  as.numeric(timing["elapsed"])
}

parse_int_list <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    return(default)
  }

  as.integer(strsplit(value, ",", fixed = TRUE)[[1]])
}

theta_matrix <- function(strength_theta, r) {
  if (strength_theta == "high") {
    return(diag(seq(0.9, 0.75, length.out = r)))
  }

  if (strength_theta == "medium") {
    return(diag(seq(0.7, 0.55, length.out = r)))
  }

  diag(seq(0.5, 0.35, length.out = r))
}

source_local_ccar3_methods <- function(ccar3_dir = Sys.getenv("CCAR3_PKG_PATH", unset = "~/Documents/ccar3/"),
                                       ccar3_code_dir = repo_root) {
  ccar3_dir <- normalizePath(ccar3_dir, winslash = "/", mustWork = TRUE)
  ccar3_code_dir <- normalizePath(ccar3_code_dir, winslash = "/", mustWork = TRUE)
  options(ccar3_pkg_path = ccar3_dir)
  Sys.setenv(CCAR3_PKG_PATH = ccar3_dir)

  ccar3_env <- new.env(parent = globalenv())

  sys.source(file.path(ccar3_dir, "R", "helpers.r"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "utils.R"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "reduced_rank_regression.R"), envir = ccar3_env)
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "SAR.R"))

  gca_to_cca_local <- function(a_estimate, S, pp) {
    p1 <- pp[1]
    p2 <- pp[2]
    p <- p1 + p2
    nnz_indices <- which(apply(a_estimate, 1, function(x) sqrt(sum(x^2))) > 0)
    nnz_indices_x <- nnz_indices[nnz_indices < (p1 + 1)]
    nnz_indices_y <- nnz_indices[nnz_indices > p1]

    if (length(nnz_indices_x) > 0) {
      sigmaxhat <- S[nnz_indices_x, nnz_indices_x, drop = FALSE]
      gram_x <- t(a_estimate[nnz_indices_x, , drop = FALSE]) %*% sigmaxhat %*% a_estimate[nnz_indices_x, , drop = FALSE]
      a_estimate[nnz_indices_x, ] <- a_estimate[nnz_indices_x, , drop = FALSE] %*% pracma::sqrtm(gram_x)$Binv
    }

    if (length(nnz_indices_y) > 0) {
      sigmayhat <- S[nnz_indices_y, nnz_indices_y, drop = FALSE]
      gram_y <- t(a_estimate[nnz_indices_y, , drop = FALSE]) %*% sigmayhat %*% a_estimate[nnz_indices_y, , drop = FALSE]
      a_estimate[nnz_indices_y, ] <- a_estimate[nnz_indices_y, , drop = FALSE] %*% pracma::sqrtm(gram_y)$Binv
    }

    list(
      U = a_estimate[1:p1, , drop = FALSE],
      V = a_estimate[(p1 + 1):p, , drop = FALSE]
    )
  }

  ccar3_env$fit_sparsecca_cv <- function(X_train, Y_train, rank,
                                         lambdax = 10^seq(from = -3, to = 1, length.out = 30),
                                         lambday = c(0),
                                         standardize = TRUE) {
    X_train <- as.matrix(X_train)
    Y_train <- as.matrix(Y_train)
    S <- stats::cov(cbind(X_train, Y_train))
    method <- SparseCCA(
      X = X_train,
      Y = Y_train,
      rank = rank,
      lambdaAseq = lambdax,
      lambdaBseq = lambday,
      max.iter = 100,
      conv = 10^-2,
      selection.criterion = 2,
      n.cv = 5,
      standardize = standardize
    )

    gca_to_cca_local(rbind(method$uhat, method$vhat), S, c(ncol(X_train), ncol(Y_train)))
  }

  ccar3_env
}

normalize_projected_loadings <- function(M, data_matrix) {
  M <- as.matrix(M)

  if (ncol(M) == 0L || all(M == 0)) {
    return(M)
  }

  gram <- crossprod(data_matrix %*% M) / nrow(data_matrix)
  sqrt_inv <- tryCatch(pracma::sqrtm(gram)$Binv, error = function(e) NULL)

  if (!is.null(sqrt_inv) && all(is.finite(sqrt_inv))) {
    return(M %*% sqrt_inv)
  }

  col_scales <- sqrt(pmax(diag(gram), 0))
  valid <- which(is.finite(col_scales) & col_scales > 0)
  if (length(valid) == 0L) {
    return(M)
  }

  M[, valid, drop = FALSE] / rep(col_scales[valid], each = nrow(M))
}

build_result_row <- function(metrics, method, theta_strength, nnzeros, seed_id,
                             r_pca, n, lambda_opt, time, status = "ok",
                             error_message = NA_character_) {
  data.frame(
    metrics,
    method = method,
    theta_strength = theta_strength,
    nnzeros = nnzeros,
    seed_id = seed_id,
    r_pca = r_pca,
    n = n,
    lambda_opt = lambda_opt,
    time = time,
    status = status,
    error_message = error_message,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

empty_metrics <- function(p, q, r) {
  metrics <- evaluate(
    matrix(0, nrow = 10, ncol = p),
    matrix(0, nrow = 10, ncol = q),
    matrix(0, nrow = p, ncol = r),
    matrix(0, nrow = q, ncol = r),
    matrix(0, nrow = p, ncol = r),
    matrix(0, nrow = q, ncol = r),
    diag(p + q),
    diag(p + q)
  )
  metrics[1, setdiff(names(metrics), c("n_new", "p1", "p2", "r"))] <- NA_real_
  metrics
}

load_required_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0("Required package '", pkg, "' is not installed."),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

load_required_package("ggplot2")
load_required_package("dplyr")
load_required_package("pracma")
load_required_package("tidyr")

setwd(repo_root)
source("experiments/simulations/generate_example_rrr.R")
source("experiments/evaluation.R")
ccar3_methods <- source_local_ccar3_methods()

args <- commandArgs(trailingOnly = TRUE)
name_exp <- if (length(args) >= 1L) args[1] else paste0("nnzeros_compare_", format(Sys.time(), "%Y%m%d_%H%M%S"))

n <- as.integer(Sys.getenv("CCAR3_NNZ_N", unset = "500"))
p <- as.integer(Sys.getenv("CCAR3_NNZ_P", unset = "500"))
r <- as.integer(Sys.getenv("CCAR3_NNZ_R", unset = "3"))
r_pca <- as.integer(Sys.getenv("CCAR3_NNZ_R_PCA", unset = "3"))
repeats <- as.integer(Sys.getenv("CCAR3_NNZ_REPEATS", unset = "20"))
strength_theta <- Sys.getenv("CCAR3_NNZ_THETA", unset = "high")
q_vals <- parse_int_list(Sys.getenv("CCAR3_NNZ_QS", unset = "5,10,20,30"), c(5L, 10L, 20L, 30L))
nnzero_values <- parse_int_list(Sys.getenv("CCAR3_NNZ_VALUES", unset = "5,7,8,10,12,15,17,20,25,30,50"), c(5L, 7L, 8L, 10L, 12L, 15L, 17L, 20L, 25L, 30L, 50L))
parallelize_cv <- identical(tolower(Sys.getenv("CCAR3_NNZ_PARALLEL", unset = "false")), "true")
cv_kfolds <- as.integer(Sys.getenv("CCAR3_NNZ_KFOLDS", unset = "5"))
lambda_points <- as.integer(Sys.getenv("CCAR3_NNZ_LAMBDA_POINTS", unset = "30"))
rrr_niter <- as.integer(Sys.getenv("CCAR3_NNZ_RRR_NITER", unset = "20000"))
rrr_thresh <- as.numeric(Sys.getenv("CCAR3_NNZ_RRR_THRESH", unset = "1e-6"))
n_new <- as.integer(Sys.getenv("CCAR3_NNZ_N_NEW", unset = "5000"))

thetas <- theta_matrix(strength_theta, r)
results_rel_path <- file.path(
  "experiments",
  "simulations",
  "results",
  paste0("rrr_compare_nnzeros_new_posprocess_mse_", name_exp, ".csv")
)
results <- data.frame()

for (seed_n in seq_len(repeats)) {
  set.seed(1000 + seed_n)

  for (q in q_vals) {
    for (nnzeros in nnzero_values) {
      if (!(max(r_pca, r, nnzeros) < p && nnzeros > max(r_pca, r))) {
        next
      }

      generation_time <- system.time({
        gen <- generate_example_sparse_U(
          n = n,
          p1 = p,
          p2 = q,
          r_pca = r_pca,
          nnzeros = nnzeros,
          theta = thetas,
          lambda_pca = 1,
          r = r,
          overlapping_amount = 1,
          normalize_diagonal = TRUE,
          n_new = n_new
        )
      })

      Sigma0_sqrt <- pracma::sqrtm(gen$Sigma)$B
      Sigma_hat_sqrt <- pracma::sqrtm(gen$S)$B

      ccar3_time <- system.time({
        ccar3_fit <- tryCatch(
          ccar3_methods$cca_rrr_cv(
            gen$X,
            gen$Y,
            r = r,
            lambdas = 10^seq(-5, 1, length.out = lambda_points),
            kfolds = cv_kfolds,
            solver = "ADMM",
            parallelize = parallelize_cv,
            LW_Sy = TRUE,
            standardize = FALSE,
            preprocess = "none",
            cv_metric = "correlation",
            rho = 1,
            niter = rrr_niter,
            thresh = rrr_thresh
          ),
          error = function(e) e
        )
      })

	      if (inherits(ccar3_fit, "error")) {
	        results <- append_and_flush_results(
	          results,
	          build_result_row(
	            empty_metrics(p, q, r),
	            method = "ccar3",
            theta_strength = strength_theta,
            nnzeros = nnzeros,
            seed_id = seed_n,
            r_pca = r_pca,
            n = n,
            lambda_opt = NA_real_,
            time = NA_real_,
            status = "error",
            error_message = conditionMessage(ccar3_fit)
	          ),
	          results_rel_path
	        )
	      } else {
	        results <- append_and_flush_results(
	          results,
	          build_result_row(
	            evaluate(
	              gen$Xnew,
              gen$Ynew,
              ccar3_fit$U[, seq_len(r), drop = FALSE],
              ccar3_fit$V[, seq_len(r), drop = FALSE],
              gen$u,
              gen$v,
              Sigma_hat_sqrt = Sigma_hat_sqrt,
              Sigma0_sqrt = Sigma0_sqrt
            ),
            method = "ccar3",
            theta_strength = strength_theta,
            nnzeros = nnzeros,
            seed_id = seed_n,
            r_pca = r_pca,
            n = n,
            lambda_opt = ccar3_fit$lambda,
            time = elapsed_seconds(ccar3_time)
	          ),
	          results_rel_path
	        )
	      }

      sar_time <- system.time({
        sar_fit <- tryCatch(
          ccar3_methods$fit_sparsecca_cv(
            gen$X,
            gen$Y,
            rank = r,
            lambdax = 10^seq(-3, 1, length.out = lambda_points),
            lambday = c(0),
            standardize = TRUE
          ),
          error = function(e) e
        )
      })

	      if (inherits(sar_fit, "error")) {
	        results <- append_and_flush_results(
	          results,
	          build_result_row(
	            empty_metrics(p, q, r),
	            method = "FIT_SAR_CV",
            theta_strength = strength_theta,
            nnzeros = nnzeros,
            seed_id = seed_n,
            r_pca = r_pca,
            n = n,
            lambda_opt = NA_real_,
            time = NA_real_,
            status = "error",
            error_message = conditionMessage(sar_fit)
	          ),
	          results_rel_path
	        )
	      } else {
	        U_sar <- normalize_projected_loadings(sar_fit$U[, seq_len(r), drop = FALSE], gen$X)
	        V_sar <- normalize_projected_loadings(sar_fit$V[, seq_len(r), drop = FALSE], gen$Y)

	        results <- append_and_flush_results(
	          results,
	          build_result_row(
	            evaluate(
	              gen$Xnew,
              gen$Ynew,
              U_sar,
              V_sar,
              gen$u,
              gen$v,
              Sigma_hat_sqrt = Sigma_hat_sqrt,
              Sigma0_sqrt = Sigma0_sqrt
            ),
            method = "FIT_SAR_CV",
            theta_strength = strength_theta,
            nnzeros = nnzeros,
            seed_id = seed_n,
            r_pca = r_pca,
            n = n,
            lambda_opt = NA_real_,
            time = elapsed_seconds(sar_time)
	          ),
	          results_rel_path
	        )
	      }
	    }
	  }
	}

results_path <- write_results_csv(results, results_rel_path)

summary_results <- results %>%
  dplyr::filter(status == "ok", method %in% c("ccar3", "FIT_SAR_CV")) %>%
  dplyr::group_by(method, p2, nnzeros) %>%
  dplyr::summarise(
    distance_tot_mean = mean(distance_tot, na.rm = TRUE),
    distance_tot_sd = stats::sd(distance_tot, na.rm = TRUE),
    distance_tot_se = distance_tot_sd / sqrt(dplyr::n()),
    distance_tot_ci_lo = distance_tot_mean - 1.96 * distance_tot_se,
    distance_tot_ci_hi = distance_tot_mean + 1.96 * distance_tot_se,
    counts = dplyr::n(),
    .groups = "drop"
  )

if (nrow(summary_results) == 0L) {
  stop("No successful fits were produced, so no plot could be generated.", call. = FALSE)
}

summary_rel_path <- file.path(
  "experiments",
  "simulations",
  "results",
  paste0("rrr_compare_nnzeros_summary_", name_exp, ".csv")
)
summary_path <- write_results_csv(summary_results, summary_rel_path)

method_colors <- c("ccar3" = "#111111", "FIT_SAR_CV" = "#0B8F5A")
method_labels <- c("ccar3" = "ccar3", "FIT_SAR_CV" = "FIT_SAR_CV")

plot_obj <- ggplot2::ggplot(
  summary_results,
  ggplot2::aes(x = nnzeros, y = distance_tot_mean, color = method, group = method)
) +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = distance_tot_ci_lo, ymax = distance_tot_ci_hi),
    width = 0.8,
    alpha = 0.8
  ) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_point(size = 2) +
  ggplot2::facet_wrap(
    ~p2,
    nrow = 2,
    scales = "free_y",
    labeller = ggplot2::as_labeller(c(`5` = "q = 5", `10` = "q = 10", `20` = "q = 20", `30` = "q = 30"))
  ) +
  ggplot2::scale_color_manual(values = method_colors, labels = method_labels) +
  ggplot2::labs(
    title = paste0("Mean Subspace Distance vs Number of Nonzeros (theta = ", strength_theta, ", error bars = 95% CI)"),
    x = "Number of nonzero rows in U",
    y = "Mean subspace distance",
    color = "Method"
  ) +
  ggplot2::theme_bw(base_size = 13) +
  ggplot2::theme(
    legend.position = "top",
    panel.grid.minor = ggplot2::element_blank()
  )

plot_png_rel_path <- file.path(
  "experiments",
  "simulations",
  "results",
  paste0("rrr_compare_nnzeros_", name_exp, ".png")
)
plot_pdf_rel_path <- file.path(
  "experiments",
  "simulations",
  "results",
  paste0("rrr_compare_nnzeros_", name_exp, ".pdf")
)

plot_png_path <- write_plot_file(plot_obj, plot_png_rel_path)
plot_pdf_path <- write_plot_file(plot_obj, plot_pdf_rel_path)

cat("Wrote results to: ", results_path, "\n", sep = "")
cat("Wrote summary to: ", summary_path, "\n", sep = "")
cat("Wrote plot PNG to: ", plot_png_path, "\n", sep = "")
cat("Wrote plot PDF to: ", plot_pdf_path, "\n", sep = "")
