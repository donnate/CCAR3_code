if (!requireNamespace("ccar3", quietly = TRUE)) {
  stop("Package 'ccar3' must be installed to use the group reduced-rank regression helpers.")
}

.ccar3_prepare_xy_group_rrr <- function(X, Y, do.scale) {
  if (isTRUE(do.scale)) {
    return(list(X = scale(X), Y = scale(Y)))
  }

  list(
    X = scale(X, scale = FALSE),
    Y = scale(Y, scale = FALSE)
  )
}

.ccar3_warn_group_rrr_penalty <- function(Kx, lambda_Kx) {
  if (!is.null(Kx) || (!is.null(lambda_Kx) && !isTRUE(all.equal(lambda_Kx, 0)))) {
    warning(
      "Ignoring Kx/lambda_Kx: the ccar3 package interface does not expose this penalty.",
      call. = FALSE
    )
  }
}

.ccar3_finalize_group_rrr_cv <- function(fit, resultsx) {
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

CCA_group_rrr <- function(X, Y,
                          groups,
                          Sx = NULL, Sy = NULL,
                          Sxy = NULL,
                          lambda = 0, Kx = NULL, r,
                          do.scale = FALSE, scale = NULL, lambda_Kx = 0,
                          LW_Sy = FALSE,
                          solver = "ADMM",
                          rho = 1,
                          niter = 1e4,
                          thresh = 1e-4,
                          thresh_0 = 1e-6,
                          verbose = FALSE) {
  if (!is.null(scale)) {
    do.scale <- scale
  }

  .ccar3_warn_group_rrr_penalty(Kx, lambda_Kx)
  xy <- .ccar3_prepare_xy_group_rrr(X, Y, do.scale)

  ccar3::cca_group_rrr(
    X = xy$X,
    Y = xy$Y,
    groups = groups,
    Sx = Sx,
    Sy = Sy,
    Sxy = Sxy,
    lambda = lambda,
    r = r,
    standardize = FALSE,
    LW_Sy = LW_Sy,
    solver = solver,
    rho = rho,
    niter = niter,
    thresh = thresh,
    thresh_0 = thresh_0,
    verbose = verbose
  )
}


CCA_group_rrr.CV <- function(X, Y,
                             groups,
                             r = 2, Kx = NULL, lambda_Kx = 0,
                             param_lambda = 10^seq(-3, 1.5, length.out = 10),
                             kfolds = 5,
                             parallelize = FALSE,
                             do.scale = FALSE,
                             scale = NULL,
                             LW_Sy = FALSE,
                             solver = "ADMM",
                             rho = 1,
                             niter = 1e4,
                             thresh = 1e-4,
                             thresh_0 = 1e-6,
                             verbose = FALSE,
                             nb_cores = NULL) {
  if (!is.null(scale)) {
    do.scale <- scale
  }

  .ccar3_warn_group_rrr_penalty(Kx, lambda_Kx)
  xy <- .ccar3_prepare_xy_group_rrr(X, Y, do.scale)
  cv_folds <- get("cca_group_rrr_cv_folds", envir = asNamespace("ccar3"))

  resultsx <- data.frame(
    lambda = param_lambda,
    rmse = vapply(
      param_lambda,
      function(lambda) {
        cv_folds(
          X = xy$X,
          Y = xy$Y,
          groups = groups,
          Sx = NULL,
          Sy = NULL,
          kfolds = kfolds,
          lambda = lambda,
          r = r,
          standardize = FALSE,
          LW_Sy = LW_Sy,
          solver = solver,
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

  fit <- ccar3::cca_group_rrr_cv(
    X = xy$X,
    Y = xy$Y,
    groups = groups,
    r = r,
    lambdas = param_lambda,
    kfolds = kfolds,
    parallelize = parallelize,
    standardize = FALSE,
    LW_Sy = LW_Sy,
    solver = solver,
    rho = rho,
    thresh_0 = thresh_0,
    niter = niter,
    thresh = thresh,
    verbose = verbose,
    nb_cores = nb_cores
  )

  .ccar3_finalize_group_rrr_cv(fit, resultsx)
}
