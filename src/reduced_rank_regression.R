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

.ccar3_finalize_rrr_cv <- function(fit, resultsx) {
  list(
    ufinal = fit$U,
    vfinal = fit$V,
    U = fit$U,
    V = fit$V,
    lambda = fit$lambda,
    resultsx = resultsx,
    rmse = resultsx$rmse,
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
                       verbose = FALSE,
                       nb_cores = NULL) {
  if (!is.null(scale)) {
    do.scale <- scale
  }

  .ccar3_warn_rrr_penalty(Kx, lambda_Kx)
  xy <- .ccar3_prepare_xy_rrr(X, Y, do.scale)
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

  resultsx$rmse[is.na(resultsx$rmse) | resultsx$rmse == 0] <- 1e8
  resultsx <- resultsx[resultsx$rmse > 1e-5, , drop = FALSE]

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

  .ccar3_finalize_rrr_cv(fit, resultsx)
}
