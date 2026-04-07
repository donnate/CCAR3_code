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

source_local_ccar3_methods <- function(ccar3_dir = Sys.getenv("CCAR3_PKG_PATH", unset = "/scratch/midway3/cdonnat/ccar3")) {
  ccar3_dir <- normalizePath(ccar3_dir, winslash = "/", mustWork = TRUE)
  options(ccar3_pkg_path = ccar3_dir)
  Sys.setenv(CCAR3_PKG_PATH = ccar3_dir)

  ccar3_env <- new.env(parent = globalenv())

  sys.source(file.path(ccar3_dir, "R", "helpers.r"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "utils.R"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "reduced_rank_regression.R"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "group_reduced_rank_regression.R"), envir = ccar3_env)
  sys.source(file.path(ccar3_dir, "R", "graph_reduced_rank_regression.R"), envir = ccar3_env)

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

library(ggplot2)
library(dplyr)
library(CCA)
library(tidyr)
library(pracma)

source("experiments/simulations/generate_example_rrr.R")
source('experiments/experiment_functions.R')
source('experiments/alternative_methods/SAR.R')
source('experiments/alternative_methods/Parkhomenko.R')
source('experiments/alternative_methods/Witten_CrossValidation.R')
source('experiments/alternative_methods/Waaijenborg.R')
source('experiments/alternative_methods/scca_chao.R')
source("experiments/evaluation.R")
ccar3_methods <- source_local_ccar3_methods()

args <- commandArgs(trailingOnly=TRUE)
seed <- as.numeric(args[1])
print(seed)
name_exp <- args[2]
set.seed(seed)
n <-500#as.numeric(args[3])
strength_theta <- args[3]
# <- as.numeric(300, 500, 700)
rs <- c(as.numeric(args[4]))
ps <- c(as.numeric(args[5]))
qvalues <- c(as.numeric(args[6]))
overlaps <- c(1)
#props <- c(0, 0.1, 0.2)
props <- c(0)
noise = 1
seeds = 1:100
normalize_diagonal = TRUE
LW_Sy = TRUE
nnzero_values = c(5)

result = c()
print(ps)
print(qvalues)
print(rs)
for(seed_n in seeds){
  #for (n in c(100, 300, 500, 1000, 10000)){
  set.seed(seed * 100 + seed_n)
  print("Start loop")
  for(nnzeros in nnzero_values){
    for(p in c(ps)){
    #for (p in c(20, 50, 80, 100, 200, 500, 1000)){
      for (q in qvalues){
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
		print("here")
		print(p)
		print(max(r_pca, r, nnzeros) )
		print( (max(r_pca, r, nnzeros) < p) )
          if ( (max(r_pca, r, nnzeros) < p) ) {
            for (overlapping_amount in overlaps){
              for(prop_missing in props){
                cat("seed:")
                cat(seed, " ")
                print(c(n, r, r_pca, strength_theta))
                start_time_creation <- system.time({
                gen = generate_example_group(n, p, q, 
                                  r_pca = r_pca,
                                  nnzeros = nnzeros,
                                  theta = thetas,
                                  lambda_pca = 1,
                                  r = r, overlapping_amount = overlapping_amount,
                                  normalize_diagonal = normalize_diagonal,
                                  n_new = 5000) 
              })
                print("Generation time")
              print(start_time_creation[[1]])
                print("We're done generating")
                X = gen$X
                Y = gen$Y
                groups = gen$groups
                Sigma0_sqrt = sqrtm(gen$Sigma)$B
                Sigma_hat_sqrt = sqrtm(gen$S)$B
                
                
                
                #### Oracle
                print("beginning oracle")
                set_u =  which(apply(gen$u,1, norm)>0)
                set_v =  which(apply(gen$v,1, norm)>0)
                t=CCA::cc(as.matrix(gen$X[,set_u]), as.matrix(gen$Y[, set_v]))
                Uhat = matrix(0, p, r)
                Vhat = matrix(0, q, r)
                Uhat[set_u, ] <-  t$xcoef[, 1:r]
                Vhat[set_v, ] <-  t$ycoef[, 1:r]
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, Uhat, 
                                                            Vhat, 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  method = "Oracle",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "exp" = seed * 100 + seed_n,
                                                   "lambda_opt" = 0,
                                                   "time" = 0
                                                   
                )
                )
                
                print(paste0("Starting ", "Alt opt") )
                
                tryCatch({
                start_time_alt2 <- system.time({
                  res_alt = ccar3_methods$cca_rrr_cv(X, Y,
                                                     r = r,
                                                     standardize = TRUE,
                                                     lambdas = c(10^seq(-3, 1, length.out = 30)),
                                                     kfolds = 3, solver = "rrr", LW_Sy = LW_Sy)
                })
                res_alt <- normalize_cv_fit(res_alt)
                res_alt$ufinal[which(is.na( res_alt$ufinal))] <- 0
                res_alt$vfinal[which(is.na( res_alt$vfinal))] <- 0
                Uhat <- res_alt$ufinal[, 1:r]
                Vhat <- res_alt$vfinal[, 1:r]
                if (sum(apply(res_alt$ufinal, 1, function(x){sum(x!=0)}) >0) <r){
                  #### Choose another lambda
                  lambda_chosen = max(res_alt$resultsx$lambda[which(res_alt$resultsx$rmse > 1.05 * min(res_alt$resultsx$rmse))])
                  start_time_alt <- system.time({
                    res_alt <- ccar3_methods$cca_rrr(X, Y, Sx = NULL, Sy = NULL,
                                                     lambda = lambda_chosen, r = r, highdim = TRUE,
                                                     solver = "rrr",
                                                     LW_Sy = LW_Sy, standardize = TRUE)
                    res_alt$U[which(is.na(res_alt$U))] <- 0
                    res_alt$V[which(is.na(res_alt$V))] <- 0
                    Uhat <- res_alt$U[, 1:r]
                    Vhat <- res_alt$V[, 1:r]
                    
                  })
                }
                print(res_alt$ufinal)
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                            Uhat, 
                                                            res_alt$vfinal[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  method = "RRR-opt",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = res_alt$lambda,
                                                   "time" = start_time_alt2[[1]]
                )
                )
                }, error = function(e) {
                  # Print the error message
                  cat("Error occurred in RRR CV:", conditionMessage(e), "\n")
                  # Skip to the next iteration
                })
              

                print(paste0("Starting ", "Alt opt2") )
                tryCatch({
                start_time_alt3 <- system.time({
                  res_alt = ccar3_methods$cca_rrr_cv(X, Y,
                                                     r = r,
                                                     lambdas = c(10^seq(-3, 1, length.out = 30)),
                                                     kfolds = 3, solver = "ADMM", LW_Sy = LW_Sy, standardize = TRUE,
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
                                                       LW_Sy = LW_Sy, standardize = TRUE,
                                                       thresh = 1e-6)
                      res_alt$U[which(is.na(res_alt$U))] <- 0
                      res_alt$V[which(is.na(res_alt$V))] <- 0
                      Uhat <- res_alt$U[, 1:r]
                      Vhat <- res_alt$V[, 1:r]
                      
                    #})
                    
                  }
                }
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew,
                                                            Uhat,
                                                            Vhat[, 1:r],
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt,
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,
                                                   method = "RRR-ADMM-opt",
                                                   "prop_missing" = prop_missing,
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = lambda_chosen,
                                                   "time" = start_time_alt3[[1]]
                )
                )
            }, error = function(e) {
              # Print the error message
              cat("Error occurred in CVXR CV:", conditionMessage(e), "\n")
              # Skip to the next iteration
            })
          
                
          print(paste0("Starting ", "CV opt-group") )

              tryCatch({
                start_time_alt4<- system.time({
                  res_alt = ccar3_methods$cca_group_rrr_cv(X, Y,
                                                           groups = gen$groups,
                                                           r = r,
                                                           lambdas = c(10^seq(-3, 1, length.out = 30)),
                                                           kfolds = 3, solver = "ADMM", LW_Sy = LW_Sy, standardize = TRUE,
                                                           rho = 1, niter = 2 * 1e4,
                                                           thresh = 1e-6)
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
                    res_alt <- ccar3_methods$cca_group_rrr(X, Y,
                                                           groups = gen$groups,
                                                           Sx = NULL, Sy = NULL,
                                                           lambda = lambda_chosen,
                                                           r = r,
                                                           solver = "ADMM",
                                                           LW_Sy = LW_Sy, standardize = TRUE,
                                                           rho = 1, niter = 2 * 1e4,
                                                           thresh = 1e-6)
                    res_alt$U[which(is.na(res_alt$U))] <- 0
                    res_alt$V[which(is.na(res_alt$V))] <- 0
                    Uhat <- res_alt$U[, 1:r]
                    Vhat <- res_alt$V[, 1:r]
                    
                    #})
                    
                  }
                }

                print(res_alt$ufinal)
                result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                            Uhat, 
                                                            Vhat[, 1:r], 
                                                            gen$u, gen$v,
                                                            Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                            Sigma0_sqrt = Sigma0_sqrt),
                                                   "noise" = noise,  
                                                   method = "CVX-opt-group",  
                                                   "prop_missing" = prop_missing, 
                                                   "nnzeros" = nnzeros,
                                                   "theta_strength" = strength_theta,
                                                   "overlapping_amount" = overlapping_amount,
                                                   "r_pca" = r_pca,
                                                   "n" = n,
                                                   "exp" = seed * 100 + seed_n,
                                                   "normalize_diagonal" = normalize_diagonal,
                                                   "lambda_opt" = lambda_chosen,
                                                   "time" = start_time_alt4[[1]]
                )
                )
          }, error = function(e) {
            # Print the error message
            cat("Error occurred in group CV", ":", conditionMessage(e), "\n")
            # Skip to the next iteration
          })
        

          for (method in c("FIT_SAR_CV", "FIT_SAR_BIC", "Witten_Perm",
                                 "Witten.CV", "Waaijenborg-Author", "Waaijenborg-CV",
                                 "SCCA_Parkhomenko", "Fantope", "Chao", "SGCA")){
                  
                  print(paste0("Starting ", method))
                  
                  
                  tryCatch({
                    start_time_additional_method <- system.time({
                      test1<-additional_checks(gen$X,
                                               gen$Y, S=NULL, 
                                               rank=r, kfolds=3, 
                                               method.type = method,
                                               lambdax= 10^seq(-3,1, length.out = 30),
                                               lambday = c(0))
                    })
                    result <- rbind(result, data.frame(evaluate(gen$Xnew, gen$Ynew, 
                                                                test1$u[, 1:r], 
                                                                test1$v[, 1:r], 
                                                                gen$u, gen$v,
                                                                Sigma_hat_sqrt = Sigma_hat_sqrt, 
                                                                Sigma0_sqrt = Sigma0_sqrt),
                                                       "noise" = noise,  method = method,  
                                                       "prop_missing" = prop_missing, 
                                                       "overlapping_amount" = overlapping_amount,
                                                       "nnzeros" = nnzeros,
                                                       "theta_strength" = strength_theta,
                                                       "r_pca" = r_pca,
                                                       "n" = n,
                                                       "exp" = seed * 100 + seed_n,
                                                       "normalize_diagonal" = normalize_diagonal,
                                                       "lambda_opt" = 0,
                                                       "time" = start_time_additional_method[[1]]
                    )
                    )
                  }, error = function(e) {
                    # Print the error message
                    cat("Error occurred in method", method, ":", conditionMessage(e), "\n")
                    # Skip to the next iteration
                  })
                }
                write_results_csv(result, paste0("experiments/simulations/results/group/2024-group-newest_RRR_efficient_results", name_exp, ".csv"))
                print("Done loop")
              
              }
              print("Done inner most loop")
            }
          }
        }
      }
    }
  }

  }
}
