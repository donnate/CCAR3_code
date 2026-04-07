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

script_path <- get_script_path()
repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = TRUE)
setwd(repo_root)

write_results_csv <- function(df, relative_path) {
  out_path <- file.path(repo_root, relative_path)
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(df, out_path, row.names = FALSE)
}

elapsed_seconds <- function(timing) {
  as.numeric(timing["elapsed"])
}

source_local_ccar3_methods <- function(ccar3_dir = Sys.getenv("CCAR3_PKG_PATH", unset = "/scratch/midway3/cdonnat/ccar3"),
                                       ccar3_code_dir = repo_root) {
  ccar3_dir <- normalizePath(ccar3_dir, winslash = "/", mustWork = TRUE)
  ccar3_code_dir <- normalizePath(ccar3_code_dir, winslash = "/", mustWork = TRUE)
  options(ccar3_pkg_path = ccar3_dir)
  Sys.setenv(CCAR3_PKG_PATH = ccar3_dir)

  ccar3_env <- new.env(parent = globalenv())

  sys.source(file.path(ccar3_dir, "R", "helpers.r"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "utils.R"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "reduced_rank_regression.R"), envir = ccar3_env)
  source(file.path(ccar3_code_dir,  "experiments", "alternative_methods",  "SAR.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods",  "Parkhomenko.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods",  "Waaijenborg.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "utils.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "gca_to_cca.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "init_process.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "sgca_init.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "sgca_tgd.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "subdistance.R"))
  source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "scca_chao.R"))

  ccar3_env$sparse_CCA_benchmarks <- local({
    env <- ccar3_env

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

    function(X_train, Y_train, S = NULL,
             rank = 2, kfolds = 5, method.type = "FIT_SAR_CV",
             lambdax = 10^seq(from = -4, to = 2, length = 30),
             lambday = c(0),
             standardize = TRUE) {
      X_train <- as.matrix(data.frame(X_train) %>% dplyr::mutate_all(~tidyr::replace_na(., mean(., na.rm = TRUE))))
      Y_train <- as.matrix(data.frame(Y_train) %>% dplyr::mutate_all(~tidyr::replace_na(., mean(., na.rm = TRUE))))
      p1 <- ncol(X_train)
      p2 <- ncol(Y_train)
      pp <- c(p1, p2)

      if (is.null(S)) {
        S <- stats::cov(cbind(X_train, Y_train))
      }

      if (method.type == "FIT_SAR_BIC") {
        method <- env$SparseCCA(
          X = X_train, Y = Y_train, rank = rank,
          lambdaAseq = lambdax, lambdaBseq = lambday,
          max.iter = 100, conv = 10^-2,
          selection.criterion = 1, n.cv = kfolds,
          standardize = standardize
        )
        a_estimate <- rbind(method$uhat, method$vhat)
      } else if (method.type == "FIT_SAR_CV") {
        method <- env$SparseCCA(
          X = X_train, Y = Y_train, rank = rank,
          lambdaAseq = lambdax, lambdaBseq = lambday,
          max.iter = 100, conv = 10^-2,
          selection.criterion = 2, n.cv = kfolds,
          standardize = standardize
        )
        a_estimate <- rbind(method$uhat, method$vhat)
      } else if (method.type == "Witten_Perm") {
        method <- PMA::CCA.permute(
          x = X_train, z = Y_train,
          typex = "standard", typez = "standard",
          penaltyxs = lambdax[which(lambdax < 1)],
          penaltyzs = lambday[which(lambday < 1)],
          standardize = standardize,
          nperms = 50
        )
        fit <- PMA::CCA(
          x = X_train, z = Y_train,
          typex = "standard", typez = "standard", K = rank,
          penaltyx = method$bestpenaltyx,
          penaltyz = method$bestpenaltyz,
          trace = FALSE,
          standardize = standardize
        )
        a_estimate <- rbind(fit$u, fit$v)
      } else if (method.type == "Witten.CV") {
        method <- env$Witten.CV(
          X = X_train, Y = Y_train, n.cv = 5, rank = rank,
          lambdax = lambdax[which(lambdax < 1)],
          lambday = lambday[which(lambday < 1)],
          standardize = standardize
        )
        fit <- PMA::CCA(
          x = X_train, z = Y_train,
          typex = "standard", typez = "standard", K = rank,
          penaltyx = method$lambdax.opt,
          penaltyz = method$lambday.opt,
          trace = FALSE,
          standardize = standardize
        )
        a_estimate <- rbind(fit$u, fit$v)
      } else if (method.type == "SCCA_Parkhomenko") {
        method <- env$SCCA_Parkhomenko(
          x.data = X_train, y.data = Y_train, Krank = rank,
          lambda.v.seq = lambdax[which(lambdax < 2)],
          lambda.u.seq = lambday[which(lambday < 2)],
          standardize = standardize
        )
        a_estimate <- rbind(method$uhat, method$vhat)
      } else if (method.type == "SCCA_Waaijenborg") {
        method <- env$SCCA_Waaijenborg(
          x.data = X_train, y.data = Y_train, Krank = rank,
          lambda.v.seq = lambdax[which(lambdax < 2)],
          lambda.u.seq = lambday[which(lambday < 2)],
          standardize = standardize
        )
        a_estimate <- rbind(method$uhat, method$vhat)
      } else {
        stop("Unsupported method.type: ", method.type, call. = FALSE)
      }

      gca_to_cca_local(a_estimate, S, pp)
    }
  })

  ccar3_env
}

normalize_cv_fit <- function(fit) {
  if (!is.null(fit$U) && is.null(fit$ufinal)) {
    fit$ufinal <- fit$U
  }
  if (!is.null(fit$V) && is.null(fit$vfinal)) {
    fit$vfinal <- fit$V
  }
  fit
}

empty_evaluation_row <- function(X, Y, rank) {
  metrics <- evaluate(
    X,
    Y,
    matrix(0, ncol(X), rank),
    matrix(0, ncol(Y), rank),
    matrix(0, ncol(X), rank),
    matrix(0, ncol(Y), rank),
    diag(ncol(X) + ncol(Y)),
    diag(ncol(X) + ncol(Y))
  )
  metrics[1, setdiff(names(metrics), c("n_new", "p1", "p2", "r"))] <- NA_real_
  metrics
}

build_result_row <- function(metrics, noise, method, prop_missing,
                             overlapping_amount, nnzeros, theta_strength,
                             r_pca, n, exp, normalize_diagonal,
                             lambda_opt, time, status = "ok",
                             error_message = NA_character_) {
  data.frame(
    metrics,
    "noise" = noise,
    method = method,
    "prop_missing" = prop_missing,
    "overlapping_amount" = overlapping_amount,
    "nnzeros" = nnzeros,
    "theta_strength" = theta_strength,
    "r_pca" = r_pca,
    "n" = n,
    "exp" = exp,
    "normalize_diagonal" = normalize_diagonal,
    "lambda_opt" = lambda_opt,
    "time" = time,
    status = status,
    error_message = error_message,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

ccar3_benchmark_methods <- c(
  "FIT_SAR_CV",
  "FIT_SAR_BIC",
  "Witten_Perm",
  "Witten.CV",
  "SCCA_Parkhomenko"
)

ensure_loading_matrix <- function(M, n_rows, rank) {
  M <- as.matrix(M)

  if (nrow(M) != n_rows && ncol(M) == n_rows) {
    M <- t(M)
  }

  if (nrow(M) != n_rows) {
    stop(
      paste0(
        "Unexpected loading matrix shape: got ",
        nrow(M),
        "x",
        ncol(M),
        ", expected ",
        n_rows,
        " rows."
      ),
      call. = FALSE
    )
  }

  if (ncol(M) < rank) {
    M <- cbind(M, matrix(0, nrow = n_rows, ncol = rank - ncol(M)))
  }

  M[, seq_len(rank), drop = FALSE]
}

run_benchmark_method <- function(ccar3_methods, X, Y, rank, method, lambdax, lambday) {
  #### assumes everything is centered
  fit <- ccar3_methods$sparse_CCA_benchmarks(
    X,
    Y,
    S = NULL,
    rank = rank,
    kfolds = 5,
    method.type = method,
    lambdax = lambdax,
    lambday = lambday
  )
  list(
    u = ensure_loading_matrix(fit$U, ncol(X), rank),
    v = ensure_loading_matrix(fit$V, ncol(Y), rank)
  )
}

load_required_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Required package '",
        pkg,
        "' is not installed. Install it with: install.packages('",
        pkg,
        "')"
      ),
      call. = FALSE
    )
  }

  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
}

load_required_package("ggplot2")
load_required_package("dplyr")
load_required_package("tidyr")
load_required_package("pracma")

source("experiments/simulations/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source("experiments/evaluation.R")
ccar3_methods <- source_local_ccar3_methods()


#Simulation for missing values in both X and Y

args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 7) {
  seed <- as.numeric(args[1])
  name_exp <- args[2]
  n <- as.numeric(args[3])
  strength_theta <- args[4]
  p <- as.numeric(args[5])
  rs <- c(as.numeric(args[6]))
  q_vals <- c(as.numeric(args[7]))
} else if (length(args) == 5) {
  seed <- 1
  name_exp <- paste0("manual_", format(Sys.time(), "%Y%m%d_%H%M%S"))
  n <- as.numeric(args[1])
  strength_theta <- args[2]
  p <- as.numeric(args[3])
  rs <- c(as.numeric(args[4]))
  q_vals <- c(as.numeric(args[5]))
  cat("Detected 5 arguments. Using seed=1 and auto-generated name_exp=", name_exp, "\n", sep = "")
} else {
  stop(
    "Expected either 7 args: seed name_exp n strength_theta p r q, or 5 args: n strength_theta p r q",
    call. = FALSE
  )
}

print(seed)
set.seed(seed)
overlaps <- c(1)
#props <- c(0, 0.1, 0.2)
props <- c(0)
noise = 1
seeds = 1:50
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(10)
result = data.frame()
ran_any_configuration <- FALSE
for (seed_n in seeds){
  #for (n in c(100, 300, 500, 1000, 10000)){
  set.seed(seed * 100 + seed_n)
  for (nnzeros in nnzero_values){
    #for(p in c(100,  200, 300,  500, 800, 80, 20)){
    #for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in q_vals){
      #for(nnzeros in c(5, 10, 15, 20, 50)){
      for (r in rs){
        if ( strength_theta == "high"){
          thetas <- diag(seq(0.9, 0.75, length.out = r))
        }else{
          if ( strength_theta == "medium"){
            thetas <- diag(seq(0.7, 0.55, length.out = r))
          }
          else{
            thetas <- diag(seq(0.5, 0.35, length.out = r))
          }
        }
        for (r_pca in c(5)){
          nnzeros_effective <- min(nnzeros, p - 1)

          if (nnzeros_effective != nnzeros) {
            cat(
              "Adjusting nnzeros from ",
              nnzeros,
              " to ",
              nnzeros_effective,
              " so that nnzeros < p for this run.\n",
              sep = ""
            )
          }

          if ( (max(r_pca, r, nnzeros_effective) < p) & (nnzeros_effective > max(r_pca, r)) ) {
            ran_any_configuration <- TRUE
            for (overlapping_amount in overlaps){
              for(prop_missing in props){
                cat("seed:")
                cat(seed, " ")
                print(c(n, r, r_pca, strength_theta))
                #gen = generate(n, p, q, s, prop_missing)
                start_time_creation <- system.time({
                gen = generate_example_sparse_U(n, p, q,
                                                r_pca = r_pca,
                                                nnzeros = nnzeros_effective,
                                                theta = thetas,
                                                lambda_pca = 1,
                                                r = r,
                                                overlapping_amount = overlapping_amount,
                                                normalize_diagonal = normalize_diagonal,
                                                n_new = 5000) 
                })
                print(elapsed_seconds(start_time_creation))
                X = gen$X
                Y = gen$Y
                #Xna = gen$Xna
                #Yna = gen$Yna
                Sigma0_sqrt = sqrtm(gen$Sigma)$B
                Sigma_hat_sqrt = sqrtm(gen$S)$B
                
                

                #### Oracle
                if (requireNamespace("CCA", quietly = TRUE)) {
                  print("beginning oracle")
                  set_u =  which(apply(gen$u,1, norm)>0)
                  set_v =  which(apply(gen$v,1, norm)>0)
                  t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                  Uhat = matrix(0, p, r)
                  Vhat = matrix(0, q, r)
                  Uhat[set_u, ] <-  t$xcoef[, 1:r]
                  Vhat[set_v, ] <-  t$ycoef[, 1:r]
                  result <- rbind(result, build_result_row(
                    metrics = evaluate(gen$Xnew, gen$Ynew, Uhat,
                                       Vhat,
                                       gen$u, gen$v,
                                       Sigma_hat_sqrt = Sigma_hat_sqrt,
                                       Sigma0_sqrt = Sigma0_sqrt),
                    noise = noise,
                    method = "Oracle",
                    prop_missing = prop_missing,
                    overlapping_amount = overlapping_amount,
                    nnzeros = nnzeros_effective,
                    theta_strength = strength_theta,
                    r_pca = r_pca,
                    n = n,
                    exp = seed * 100 + seed_n,
                    normalize_diagonal = normalize_diagonal,
                    lambda_opt = 0,
                    time = 0
                  ))
                } else {
                  cat("Skipping Oracle: package 'CCA' is not available.\n")
                }
                
                print(paste0("Starting ", "Alt opt") )
                tryCatch({
                  start_time_alt3 <- system.time({
                    res_alt = ccar3_methods$cca_rrr_cv(X, Y,
                                                       r = r,
                                                       lambdas = c(10^seq(-3, 1, length.out = 30)),
                                                       kfolds = 5, solver = "ADMM", LW_Sy = LW_Sy,
                                                       parallelize = TRUE,
                                                       standardize = FALSE,
                                                       preprocess = "none",
                                                       rho = 1, niter = 2 * 1e4, thresh = 1e-6)
                  })
                  res_alt <- normalize_cv_fit(res_alt)
                  res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                  res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                  Uhat <- res_alt$ufinal[, 1:r]
                  Vhat <- res_alt$vfinal[, 1:r]
                  lambda_chosen = res_alt$lambda
                  if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                    #### Choose another lambda
                    while(sum(apply(Uhat, 1, function(x){sum(x!=0)}) >0) <r){
                      lambda_chosen = lambda_chosen / 2
                      #start_time_alt <- system.time({
                      res_alt <- ccar3_methods$cca_rrr(X, Y, Sx = NULL, Sy = NULL,
                                                       lambda = lambda_chosen,
                                                       r = r,
                                                       highdim = TRUE,
                                                       solver = "ADMM",
                                                       LW_Sy = LW_Sy, standardize = FALSE,
                                                       preprocess = "none",
                                                       thresh = 1e-6)
                      res_alt$U[which(is.na(res_alt$U))] <- 0
                      res_alt$V[which(is.na(res_alt$V))] <- 0
                      Uhat <- res_alt$U[, 1:r]
                      Vhat <- res_alt$V[, 1:r]
                      
                      #})
                      
                    }
                  }
                result <- rbind(result, build_result_row(
                  metrics = evaluate(gen$Xnew, gen$Ynew,
                                     Uhat,
                                     Vhat[, 1:r],
                                     gen$u, gen$v,
                                     Sigma_hat_sqrt = Sigma_hat_sqrt,
                                     Sigma0_sqrt = Sigma0_sqrt),
                  noise = noise,
                  method = "ccar3",
                  prop_missing = prop_missing,
                  overlapping_amount = overlapping_amount,
                  nnzeros = nnzeros_effective,
                  theta_strength = strength_theta,
                  r_pca = r_pca,
                  n = n,
                  exp = seed * 100 + seed_n,
                  normalize_diagonal = normalize_diagonal,
                  lambda_opt = res_alt$lambda,
                  time = elapsed_seconds(start_time_alt3)
                ))
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in crrr:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
                  
                  
		
		
		             print(paste0("Starting ", "Alt one it") )
                tryCatch({
                  start_time_alt4 <- system.time({
                                        res_alt1 <- ccar3_methods$cca_rrr(X, Y, Sx = NULL, Sy = NULL,
                                                                          lambda = 0.01,
                                                                          r = r,
                                                                           
                                                                          solver = "ADMM",
                                                                          LW_Sy = LW_Sy,        standardize = FALSE,
                                                                         preprocess = "none",
                                                                          thresh = 1e-6)
		  })
                  res_alt1$U[which(is.na( res_alt1$U))] <- 0
                  res_alt1$V[which(is.na( res_alt1$V))] <- 0
                  Uhat <- res_alt1$U[, 1:r]
                  Vhat <- res_alt1$V[, 1:r]
                  result <- rbind(result, build_result_row(
                    metrics = evaluate(gen$Xnew, gen$Ynew,
                                       Uhat,
                                       Vhat,
                                       gen$u, gen$v,
                                       Sigma_hat_sqrt = Sigma_hat_sqrt,
                                       Sigma0_sqrt = Sigma0_sqrt),
                    noise = noise,
                    method = "cca_rrr_-ADMM-one-iteration",
                    prop_missing = prop_missing,
                    overlapping_amount = overlapping_amount,
                    nnzeros = nnzeros_effective,
                    theta_strength = strength_theta,
                    r_pca = r_pca,
                    n = n,
                    exp = seed * 100 + seed_n,
                    normalize_diagonal = normalize_diagonal,
                    lambda_opt = 0.01,
                      time = elapsed_seconds(start_time_alt4)
                  ))
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
                write_results_csv(result, paste0("experiments/simulations/results/2026_newest_RRR_efficient_results", name_exp, ".csv"))





		
		for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                                 "SCCA_Parkhomenko", "Fantope", "Chao", "SGCA")){
                  
                  print(paste0("Starting ", method))
                  
                  
                  tryCatch({
                    start_time_additional_method <- system.time({
                      if (method %in% ccar3_benchmark_methods) {
                        test1 <- run_benchmark_method(
                          ccar3_methods = ccar3_methods,
                          X = gen$X,
                          Y = gen$Y,
                          rank = r,
                          method = method,
                          lambdax = 10^seq(-3, 1, length.out = 30),
                          lambday = c(0)
                        )
                      } else {
                        test1 <- additional_checks(gen$X,
                                                   gen$Y, S=NULL, 
                                                   rank=r, kfolds=5, 
                                                   method.type = method,
                                                   lambdax= 10^seq(-3,1, length.out = 30),
                                                   lambday = c(0))
                      }
                      test1$u = test1$u %*% pracma::sqrtm(crossprod(gen$X %*% test1$u)/ nrow(gen$X))$Binv
                      test1$v = test1$v %*% pracma::sqrtm(crossprod(gen$Y %*% test1$v)/ nrow(gen$Y))$Binv
                    })
                    
                    result <- rbind(result, build_result_row(
                      metrics = evaluate(gen$Xnew, gen$Ynew,
                                         test1$u[, 1:r],
                                         test1$v[, 1:r],
                                         gen$u, gen$v,
                                         Sigma_hat_sqrt = Sigma_hat_sqrt,
                                         Sigma0_sqrt = Sigma0_sqrt),
                      noise = noise,
                      method = method,
                      prop_missing = prop_missing,
                      overlapping_amount = overlapping_amount,
                      nnzeros = nnzeros_effective,
                      theta_strength = strength_theta,
                      r_pca = r_pca,
                      n = n,
                      exp = seed * 100 + seed_n,
                      normalize_diagonal = normalize_diagonal,
                      lambda_opt = 0,
                      time = elapsed_seconds(start_time_additional_method)
                    ))
                  }, error = function(e) {
                    result <- rbind(result, build_result_row(
                      metrics = empty_evaluation_row(gen$Xnew, gen$Ynew, r),
                      noise = noise,
                      method = method,
                      prop_missing = prop_missing,
                      overlapping_amount = overlapping_amount,
                      nnzeros = nnzeros_effective,
                      theta_strength = strength_theta,
                      r_pca = r_pca,
                      n = n,
                      exp = seed * 100 + seed_n,
                      normalize_diagonal = normalize_diagonal,
                      lambda_opt = NA_real_,
                      time = NA_real_,
                      status = "error",
                      error_message = conditionMessage(e)
                    ))

                    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                  })
                }
                
                write_results_csv(result, paste0("experiments/simulations/results/2026_newest_RRR_efficient_results", name_exp, ".csv"))
              
              }
            }
          } else {
            cat(
              "Skipping configuration with p=", p,
              ", r=", r,
              ", r_pca=", r_pca,
              ", nnzeros=", nnzeros_effective,
              ": require max(r_pca, r, nnzeros) < p and nnzeros > max(r_pca, r).\n",
              sep = ""
            )
          }
       # }
      }
      #}
    }
  }
#}
}
}

if (!ran_any_configuration) {
  stop(
    "No valid simulation configuration was executed. Check that p > max(r_pca, r) and choose nnzeros in (max(r_pca, r), p).",
    call. = FALSE
  )
}
