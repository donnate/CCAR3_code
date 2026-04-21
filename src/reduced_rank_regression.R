if (!requireNamespace("ccar3", quietly = TRUE)) {
  stop("Package 'ccar3' must be installed to use the reduced-rank regression helpers.")
}

.ccar3_prepare_xy_rrr <- function(X, Y, do.scale) {
  if (isTRUE(do.scale)) {
    return(list(X = scale(X), Y = scale(Y)))
  }

  list(
    X = scale(X, scale = FALSE),
    Y = scale(Y, scale = FALSE)
  )
}

.ccar3_warn_rrr_penalty <- function(Kx, lambda_Kx) {
  if (!is.null(Kx) || (!is.null(lambda_Kx) && !isTRUE(all.equal(lambda_Kx, 0)))) {
    warning(
      "Ignoring Kx/lambda_Kx: the ccar3 package interface does not expose this penalty.",
      call. = FALSE
    )
  }
}

.ccar3_match_cv_metric <- function(cv_metric) {
  if (missing(cv_metric) || is.null(cv_metric) || length(cv_metric) == 0L || all(is.na(cv_metric))) {
    metric <- "mse"
  } else {
    metric <- tolower(as.character(cv_metric[[1]]))
  }

  if (metric %in% c("mse", "rmse", "prediction")) {
    return("mse")
  }

  if (metric %in% c("correlation", "cor", "corr", "association")) {
    return("correlation")
  }

  stop(
    "Unsupported `cv_metric`. Use one of: 'mse' or 'correlation'.",
    call. = FALSE
  )
}

.ccar3_make_cv_folds <- function(n, kfolds) {
  create_cv_folds <- get0(
    ".create_cv_folds",
    envir = asNamespace("ccar3"),
    mode = "function",
    inherits = FALSE
  )

  if (is.function(create_cv_folds)) {
    return(create_cv_folds(n, kfolds))
  }

  fold_ids <- sample(rep(seq_len(kfolds), length.out = n))
  split(seq_len(n), fold_ids)
}

.ccar3_correlation_score <- function(X_scores, Y_scores) {
  X_scores <- as.matrix(X_scores)
  Y_scores <- as.matrix(Y_scores)

  if (nrow(X_scores) == 0L || nrow(Y_scores) == 0L ||
      ncol(X_scores) == 0L || ncol(Y_scores) == 0L) {
    return(0)
  }

  x_sd <- apply(X_scores, 2, stats::sd)
  y_sd <- apply(Y_scores, 2, stats::sd)

  keep_x <- which(is.finite(x_sd) & x_sd > 0)
  keep_y <- which(is.finite(y_sd) & y_sd > 0)

  if (length(keep_x) == 0L || length(keep_y) == 0L) {
    return(0)
  }

  C <- stats::cor(
    X_scores[, keep_x, drop = FALSE],
    Y_scores[, keep_y, drop = FALSE]
  )
  C[!is.finite(C)] <- 0

  sum(svd(C, nu = 0, nv = 0)$d)
}

.ccar3_rrr_cv_values <- function(X, Y, folds, r, lambda, solver, LW_Sy,
                                 rho, niter, thresh, thresh_0, cv_metric) {
  vapply(
    seq_along(folds),
    function(i) {
      idx_val <- folds[[i]]
      X_train <- X[-idx_val, , drop = FALSE]
      Y_train <- Y[-idx_val, , drop = FALSE]
      X_val <- X[idx_val, , drop = FALSE]
      Y_val <- Y[idx_val, , drop = FALSE]

      fit <- tryCatch(
        ccar3::cca_rrr(
          X = X_train,
          Y = Y_train,
          Sx = NULL,
          Sy = NULL,
          lambda = lambda,
          r = r,
          highdim = TRUE,
          solver = solver,
          LW_Sy = LW_Sy,
          standardize = FALSE,
          rho = rho,
          niter = niter,
          thresh = thresh,
          thresh_0 = thresh_0,
          verbose = FALSE
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        return(NA_real_)
      }

      if (cv_metric == "mse") {
        return(mean((X_val %*% fit$U - Y_val %*% fit$V)^2))
      }

      .ccar3_correlation_score(X_val %*% fit$U, Y_val %*% fit$V)
    },
    numeric(1)
  )
}

.ccar3_finalize_rrr_cv <- function(fit, resultsx, cv_metric) {
  list(
    ufinal = fit$U,
    vfinal = fit$V,
    U = fit$U,
    V = fit$V,
    lambda = fit$lambda,
    resultsx = resultsx,
    rmse = resultsx$rmse,
    cv_score = resultsx$cv_score,
    cv_metric = cv_metric,
    cor = fit$cor
  )
}

CCA_rrr <- function(X, Y, Sx = NULL, Sy = NULL,
                    lambda = 0, Kx = NULL, r, highdim = TRUE,
                    lambda_Kx = 0, solver = "ADMM",
                    LW_Sy = FALSE,
                    do.scale = TRUE,
                    scale = NULL,
                    rho = 1,
                    niter = 1e4,
                    thresh = 1e-4,
                    thresh_0 = 1e-6,
                    verbose = FALSE) {
  if (!is.null(scale)) {
    do.scale <- scale
  }

  .ccar3_warn_rrr_penalty(Kx, lambda_Kx)
  xy <- .ccar3_prepare_xy_rrr(X, Y, do.scale)

  fit <- ccar3::cca_rrr(
    X = xy$X,
    Y = xy$Y,
    Sx = Sx,
    Sy = Sy,
    lambda = lambda,
    r = r,
    highdim = highdim,
    solver = solver,
    LW_Sy = LW_Sy,
    standardize = FALSE,
    rho = rho,
    niter = niter,
    thresh = thresh,
    thresh_0 = thresh_0,
    verbose = verbose
  )

  fit$total_corr <- stats::cor(xy$X %*% fit$U, xy$Y %*% fit$V)
  fit
}


CCA_rrr.CV <- function(X, Y,
                       r = 2, Kx = NULL, lambda_Kx = 0,
                       param_lambda = 10^seq(-3, 1.5, length.out = 100),
                       kfolds = 14,
                       solver = "rrr",
                       parallelize = FALSE,
                       LW_Sy = FALSE,
                       do.scale = TRUE,
                       scale = NULL,
                       rho = 1,
                       niter = 1e4,
                       thresh = 1e-4,
                       thresh_0 = 1e-6,
                       cv_metric = "mse",
                       verbose = FALSE,
                       nb_cores = NULL) {
  if (!is.null(scale)) {
    do.scale <- scale
  }

  .ccar3_warn_rrr_penalty(Kx, lambda_Kx)
  xy <- .ccar3_prepare_xy_rrr(X, Y, do.scale)
  cv_metric <- .ccar3_match_cv_metric(cv_metric)

  if (cv_metric == "mse") {
    cv_folds <- get("cca_rrr_cv_folds", envir = asNamespace("ccar3"))

    resultsx <- data.frame(
      lambda = param_lambda,
      rmse = vapply(
        param_lambda,
        function(lambda) {
          cv_folds(
            X = xy$X,
            Y = xy$Y,
            Sx = NULL,
            Sy = NULL,
            kfolds = kfolds,
            lambda = lambda,
            r = r,
            solver = solver,
            LW_Sy = LW_Sy,
            standardize = FALSE,
            rho = rho,
            niter = niter,
            thresh = thresh,
            thresh_0 = thresh_0
          )
        },
        numeric(1)
      )
    )

    resultsx$cv_score <- resultsx$rmse
    resultsx$cv_metric <- cv_metric
    resultsx$rmse[is.na(resultsx$rmse) | resultsx$rmse == 0] <- 1e8
    resultsx <- resultsx[resultsx$rmse > 1e-5, , drop = FALSE]

    if (nrow(resultsx) == 0L) {
      stop("No valid lambda values were found during MSE cross-validation.", call. = FALSE)
    }

    fit <- ccar3::cca_rrr_cv(
      X = xy$X,
      Y = xy$Y,
      r = r,
      lambdas = param_lambda,
      kfolds = kfolds,
      solver = solver,
      parallelize = parallelize,
      LW_Sy = LW_Sy,
      standardize = FALSE,
      rho = rho,
      thresh_0 = thresh_0,
      niter = niter,
      thresh = thresh,
      verbose = verbose,
      nb_cores = nb_cores
    )
  } else {
    folds <- .ccar3_make_cv_folds(nrow(xy$X), kfolds)
    evaluate_lambda <- function(lambda) {
      .ccar3_rrr_cv_values(
        X = xy$X,
        Y = xy$Y,
        folds = folds,
        r = r,
        lambda = lambda,
        solver = solver,
        LW_Sy = LW_Sy,
        rho = rho,
        niter = niter,
        thresh = thresh,
        thresh_0 = thresh_0,
        cv_metric = cv_metric
      )
    }

    if (isTRUE(parallelize) && .Platform$OS.type != "windows") {
      detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
      if (!is.finite(detected_cores) || detected_cores < 1L) {
        detected_cores <- 1L
      }
      max_cores <- max(1L, detected_cores - 1L)
      nb_cores_effective <- if (is.null(nb_cores)) {
        min(length(param_lambda), max_cores)
      } else {
        min(as.integer(nb_cores), length(param_lambda), max_cores)
      }
      nb_cores_effective <- max(1L, nb_cores_effective)

      fold_scores <- parallel::mclapply(
        param_lambda,
        evaluate_lambda,
        mc.cores = nb_cores_effective
      )
    } else {
      fold_scores <- lapply(param_lambda, evaluate_lambda)
    }

    cv_score <- vapply(
      fold_scores,
      function(scores) {
        if (all(is.na(scores))) {
          return(NA_real_)
        }

        mean(scores, na.rm = TRUE)
      },
      numeric(1)
    )

    resultsx <- data.frame(
      lambda = param_lambda,
      rmse = -cv_score,
      cv_score = cv_score,
      cv_metric = cv_metric
    )
    resultsx <- resultsx[is.finite(resultsx$cv_score), , drop = FALSE]

    if (nrow(resultsx) == 0L) {
      stop("No valid lambda values were found during correlation cross-validation.", call. = FALSE)
    }

    opt_lambda <- resultsx$lambda[which.max(resultsx$cv_score)]
    fit <- ccar3::cca_rrr(
      X = xy$X,
      Y = xy$Y,
      Sx = NULL,
      Sy = NULL,
      lambda = opt_lambda,
      r = r,
      highdim = TRUE,
      solver = solver,
      LW_Sy = LW_Sy,
      standardize = FALSE,
      rho = rho,
      niter = niter,
      thresh = thresh,
      thresh_0 = thresh_0,
      verbose = verbose
    )
    fit$lambda <- opt_lambda
  }

  .ccar3_finalize_rrr_cv(fit, resultsx, cv_metric)
}
