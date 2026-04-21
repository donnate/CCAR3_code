library(MASS)
library(stats)
library(geigen)
library(pracma)
library(tidyverse)
library(CCA)
library(PMA)
library(mvtnorm)
library(glmnet)
library(caret)

load_optional_package <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Optional package '", pkg, "' is not installed; methods that depend on it will be skipped.")
    return(FALSE)
  }

  suppressPackageStartupMessages(
    library(pkg, character.only = TRUE)
  )
  TRUE
}

safe_source <- function(path) {
  tryCatch(
    source(path),
    error = function(e) {
      message("Skipping optional source '", path, "': ", conditionMessage(e))
      invisible(NULL)
    }
  )
}

load_optional_package("CVXR")
load_optional_package("VGAM")
load_optional_package("matlib")

source_local_methods <- function(ccar3_dir = Sys.getenv("CCAR3_PKG_PATH", unset = "~/Documents/ccar3"),
                                 ccar3_code_dir = getwd()) {
  ccar3_dir <- normalizePath(ccar3_dir, winslash = "/", mustWork = TRUE)
  ccar3_code_dir <- normalizePath(ccar3_code_dir, winslash = "/", mustWork = TRUE)

  options(ccar3_pkg_path = ccar3_dir)
  Sys.setenv(CCAR3_PKG_PATH = ccar3_dir)

  safe_source(file.path(ccar3_dir, "R", "alt_SAR.R"))
  safe_source(file.path(ccar3_dir, "R", "alt_Witten_CrossValidation.R"))
  safe_source(file.path(ccar3_dir, "R", "alt_Parkhomenko.R"))

  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "utils.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "gca_to_cca.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "init_process.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "sgca_init.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "sgca_tgd.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "GCA", "subdistance.R"))

  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "Waaijenborg.R"))
  safe_source(file.path(ccar3_code_dir, "experiments", "alternative_methods", "scca_chao.R"))
}

source_local_methods()

cv_function <- function(X, Y, 
                        kfolds=10, initu, initv,
                        lambdax, adaptive=TRUE, normalize=FALSE,
                        criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  nnz  <- numeric(length = kfolds)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    model <- adaptive_lasso(X_train, Y_train %*% initv, initu, adaptive=adaptive, 
                         lambdax, 
                         max.iter=5000, 
                         max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    # make predictions on validation data
    # compute RMSE on validation data
    ##### Normalize Uhat
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% model$Uhat - Y_val%*% initv)^2)
        if (norm(model$Uhat) == 0){ ####prevents selecting values that would make everything 0
          rmse[i] <- 1e8
        }
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% model$Uhat, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(model$Uhat^2, 1, sum) >1e-4)
    }else{
      sol <- gca_to_cca(rbind(model$Uhat, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
      rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
      nnz[i] <- sum(apply(sol$u^2, 1, sum) >1e-4)
     }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}


cv_function_tgd <- function(X, Y, Mask, kfolds=5, ainit,
                        lambda, r=2, k=20,  maxiter=1000, eta=0.001,
                        convergence=1e-3, normalize=FALSE, criterion="prediction") {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- numeric(length = kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]
  p <- p1 + p2;
  n <- nrow(X)
  pp <- c(p1,p2);
  S0 = cov(cbind(X, Y))
  
  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  for (i in seq_along(folds)) {
    # split data into training and validation sets
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    
    S = cov(cbind(X_train, Y_train))
    sigma0hat = S * Mask
    
    # fit model on training data with hyperparameters
    tryCatch(
    {
    final = sgca_tgd(A=S, B=sigma0hat,
             r=r,ainit, k, lambda = lambda, eta=eta,
             convergence=convergence,
             maxiter=maxiter, plot = FALSE, 
             scale=NULL)
    final <- gca_to_cca(final, S, pp)
    
    # make predictions on validation data
    # compute RMSE on validation data
    if (normalize == FALSE){
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% final$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*% final$u, Y_val%*% initv)))
      }
    }else{
      sol <- gca_to_cca(rbind(final$u, initv), 
                        cov(rbind(X_val, Y_val)), pp)
      if (criterion=="prediction"){
        rmse[i] <- sum((X_val %*% sol$u - Y_val%*% initv)^2)
      }else{
        rmse[i] <- sum(abs(cor(X_val %*%sol$u, Y_val%*% initv)))
      }
    }
    },
    error = function(e) {
      # If an error occurs, assign NA to the result
      rmse[i] <- NA
    })
  }
  
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
      return(1e8)
   }else{
  return(median(rmse, na.rm=TRUE))
   }
}

preselection <-function(Data, CorrelationMat, p1, r, alpha){
  p = ncol(Data)
  n=  nrow(Data)
  p2 =  p-p1
  t = apply(CorrelationMat -diag(diag(CorrelationMat)), 1, 
            function(x){max(x^2)})
  J = order(-t)[1: ceiling(alpha  * n/(r))]
  set_u = J[which(J <= p1)]
  set_v = J[which(J > p1)]
  t=CCA::cc(as.matrix(Data[,set_u]), as.matrix(Data[, set_v]))
  Uhat = matrix(0, p, r)
  Uhat[set_u, ] =  t$xcoef[,1:r]
  Uhat[set_v, ] =  t$ycoef[,1:r]
  return(Uhat)
}


pipeline_adaptive_lasso <- function(Data, Mask, sigma0hat, r, nu=1, Sigmax, 
                                    Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                    adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                    create_folds=TRUE, init ="Fantope", normalize=FALSE, alpha=0.5,
                                    criterion="prediction",
                                    fantope_solution=NULL){
  
  ### data splitting procedure 3 folds
  #maxiter=100;  lambdax=NULL;
 #adaptive=TRUE; kfolds=5;  param1=10^(seq(-5, 1, by = 0.5));
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  if(create_folds){
    folds <- createFolds(1:nrow(Data), k = 2, list = TRUE, returnTrain = FALSE)
    S1 <- cov(Data[folds[[1]],])
    S3 <- cov(Data[folds[[1]],])
    X = Data[folds[[2]],1:p1]
    Y = Data[folds[[2]],(p2+1):p]
    X1 = Data[folds[[2]],1:p1]
    Y1= Data[folds[[2]],(p2+1):p]
    sigma0hat1 <- S1 * Mask
  }else{
    S1 <- cov(as.matrix(Data))
    S3 <- S1
    X = Data[,1:p1]
    Y = Data[,(p2+1):p]
    X1 = Data[,1:p1]
    Y1 = Data[,(p2+1):p]
    sigma0hat1 <- sigma0hat
  }
 
  if (init == "Fantope"){
    if (is.null(fantope_solution)){
      ag <- sgca_init(A=S1, B=sigma0hat1, rho=0.5 * sqrt(log(p)/dim(X)[1]),
                      K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
      ainit <- init_process(ag$Pi, r) 
    }else{
      ainit <- fantope_solution
    }
      
  
  }else{
     ainit= preselection(Data, CorrelationMatrix, p1, r, alpha)
  }

  init <- gca_to_cca(ainit, S3, pp)
  print("Init done")
  initu<- init$u
  initv <- init$v


  if (is.null(lambdax)){
    ### do CV
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(X, Y, kfolds, initu, initv,
                                                  lambdax = .x, adaptive=adaptive,
                                                   normalize=normalize,
                                                  criterion=criterion)))

    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
    resultsx =NULL
  }
  if (is.null(lambday)){
    ### do CV
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(Y, X, kfolds, initv, initu,
                                                  lambdax = .x, adaptive=adaptive,
                                                  criterion=criterion)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
    results =NULL
  }
  
  ufinal = adaptive_lasso(X, Y %*% initv, initu, adaptive=adaptive, lambdax, 
                          max.iter=5000, 
                       max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
  
  vfinal = adaptive_lasso(Y, X %*% initu, initv, adaptive=adaptive, lambday, max.iter=5000, 
                       max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
  a_estimate = rbind(ufinal$Uhat, vfinal$Uhat)
  a_estimate <- gca_to_cca(a_estimate, S3, pp)
  return(list( ufinal = a_estimate$u, vfinal = a_estimate$v,
               initu=initu, initv=initv,
               Uhat= ufinal$Uhat, 
               Vhat=vfinal$Uhat,
               lambdax=lambdax, 
               lambday=lambday,
               resultsx=resultsx,
               resultsy=results
         )) #### Not too bad
}

pipeline_alternating_lasso <- function(Data, Mask, sigma0hat, r, nu=1, Sigmax, 
                                       Sigmay, maxiter=30, lambdax=NULL, lambday=NULL,
                                       adaptive=TRUE, kfolds=5, param1=10^(seq(-4, 2, by = 0.25)),
                                       create_folds=TRUE, init ="Fantope", normalize=FALSE, alpha=0.5,
                                       criterion="prediction",
                                       fantope_solution=NULL){
  
  ### Replaces the initialization step by an rcc
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  RCC_cv<-estim.regul_crossvalidation(Data[,1:p1], Data[,(p2+1):p],
                                      n.cv=5)
  method<-rcc(Data[,1:p1], Data[,(p2+1):p], 
              RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
  U0 = method$xcoef[,1:r]; ### Initial values
  V0 = method$ycoef[,1:r];
  X = Data[,1:p1]
  Y = Data[,(p2+1):p]
  Uinit = U0
  Vinit = V0
  converged = FALSE
  it = 0
  while(converged==FALSE & it <10){
    resultsx <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(X, Y, kfolds, U0, V0,
                                                  lambdax = .x, adaptive=adaptive,
                                                  normalize=normalize,
                                                  criterion=criterion)))
    
    
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6 + min(resultsx$rmse)) <0.05)
    lambdax = max(resultsx$param1[which_lambdax])
    Unew = adaptive_lasso(X, Y %*% V0, U0, adaptive=adaptive, lambdax, 
                          max.iter=5000, 
                          max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    results <- expand.grid(param1 = param1) %>%
      mutate(rmse = map_dbl(param1, ~ cv_function(Y, X, kfolds,  V0, Unew$Uhat,
                                                  lambdax = .x, adaptive=adaptive,
                                                  criterion=criterion)))
    best_hyperparams <- results[which.min(results$rmse), ]
    which_lambday = which(abs(results$rmse-min(results$rmse))/(1e-6 + min(results$rmse)) <0.05)
    lambday = max(results$param1[which_lambday])
    
    Vnew = adaptive_lasso(Y, X %*% Unew$Uhat, V0, adaptive=adaptive, lambday, max.iter=5000, 
                            max_k = 10, verbose = FALSE, THRESHOLD=1e-5)
    
    converged = (mean((Unew$Uhat-U0)^2) <1e-4) & (mean((Vnew$Uhat-V0)^2) <1e-4)
    it = it + 1
    U0  = Unew$Uhat
    V0 = Vnew$Uhat
    print(it)
    
  }
  
  return(list( ufinal = Unew, vfinal = Vnew,
               initu=Uinit, initv=Vinit,
               lambda.x=lambdax, 
               lambda.y=lambday
  ))
  
  
  
}

## Running initialization using convex relaxation

#(Sigma,Sigma0, lambda, rho, eta=0.001, nu=1,epsilon=5e-3,maxiter=1000,trace=FALSE)
pipeline_thresholded_gradient <- function(Data, Mask, sigma0hat, r=2, nu=1, Sigmax, 
                                          Sigmay, maxiter.init=30, 
                                          lambda=NULL, k=NULL, kfolds=5,
                                          maxiter=2000, convergence=1e-3, eta=1e-3,
                                          param1=10^(seq(-4, 1, by = 1)),
                                          param2=c(20, 1000), init="Fantope",
                                          normalize=FALSE,criterion="prediction",
                                          fantope_solution=NULL){
  p1 <- dim(Sigmax)[1]
  p2 <- dim(Sigmay)[1]
  p <- p1 + p2;
  n <- nrow(Data)
  pp <- c(p1,p2);
  S = cov(Data)

  if (init == "Fantope"){
    if (is.null(fantope_solution)){
      ag <- sgca_init(A=S, B=sigma0hat, rho=0.5 * sqrt(log(p)/n),
                      K=r ,nu=nu,trace=FALSE, maxiter = maxiter) ###needs to be changed to be a little more scalable
      ainit <- init_process(ag$Pi, r) 
    }else{
      ainit <- fantope_solution
    }
  }else{
      RCC_cv<-estim.regul_crossvalidation(Data[,1:p1], Data[,(p2+1):p],
                                          n.cv=5)
      method<-rcc(Data[,1:p1], Data[,(p2+1):p], 
                    RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
      ainit= rbind(method$xcoef[,1:r], method$ycoef[,1:r])
  }

  init <- gca_to_cca(ainit, S, pp)
  
  if (is.null(lambda) | is.null(k)){
    resultsx <- expand.grid(lambda = param1, k = param2) %>%
      mutate(rmse = map2_dbl(lambda, k, ~ cv_function_tgd(Data[, 1:p1], Data[, (p1+1):p], 
                                                          Mask, kfolds=5, ainit,
                                                          lambda = .x,
                                                          k = .y, r=r,
                                                          maxiter=maxiter, eta=eta, convergence=convergence,
                                                          normalize=normalize,
                                                          criterion=criterion)))
                                                          #X, Y, Mask, kfolds=5, ainit,lambda, r=2, k=20,  
                                                          #maxiter=1000, eta=0.001, convergence=1e-

    #print(resultsx)
    ###### (X, Y, Mask, kfolds=5, ainit, lambda, k=20)
  
    # print best hyperparameters and corresponding RMSE
    best_hyperparams <- resultsx[which.min(resultsx$rmse), ]
    which_lambdax = which(abs(resultsx$rmse-min(resultsx$rmse))/(1e-6  + min(resultsx$rmse)) <0.05)
    lambda = max(resultsx$lambda[which_lambdax])
    k = max(resultsx$k[which_lambdax])
    #print(c("selected", k, lambda))
  }
  final <- sgca_tgd(A=S, B=sigma0hat,
                    r=r, ainit,k=k, lambda = lambda, eta=eta,convergence=convergence,
                    maxiter=maxiter, plot=FALSE)
  a_estimate <- gca_to_cca(final, S, pp)
  return(list( ufinal = a_estimate$u, vfinal = a_estimate$v,
               initu=init$u, initv=init$v,
               final=final,
               lambda=lambda, 
               k=k,
               resultsx=resultsx
  ))
  
}

additional_checks <- function(X_train, Y_train, S=NULL, 
                              rank=2, kfolds=5, method.type="FIT_SAR_BIC",
                              lambdax = 10^seq(from=-3,to=2,length=100),
                              lambday = c(0, 1e-7, 1e-6, 1e-5)){

  force_rank_columns <- function(M, target_rows, target_rank, label) {
    M <- as.matrix(M)

    if (nrow(M) != target_rows && ncol(M) == target_rows) {
      M <- t(M)
    }

    if (nrow(M) != target_rows) {
      stop(
        paste0(
          "Unexpected shape for ",
          label,
          ": got ",
          nrow(M),
          "x",
          ncol(M),
          ", expected ",
          target_rows,
          " rows"
        )
      )
    }

    if (ncol(M) < target_rank) {
      M <- cbind(M, matrix(0, nrow = nrow(M), ncol = target_rank - ncol(M)))
    }

    if (ncol(M) > target_rank) {
      M <- M[, 1:target_rank, drop = FALSE]
    }

    M
  }

  normalize_loading_matrix <- function(M, Sigma) {
    M <- as.matrix(M)
    if (ncol(M) == 0 || all(M == 0)) {
      return(M)
    }

    active_cols <- which(colSums(M^2) > 0)
    if (length(active_cols) == 0) {
      return(M)
    }

    M_active <- M[, active_cols, drop = FALSE]
    gram <- t(M_active) %*% Sigma %*% M_active

    if (ncol(M_active) == 1) {
      scale_value <- as.numeric(gram)
      if (is.finite(scale_value) && scale_value > 0) {
        M_active <- M_active / sqrt(scale_value)
      }
    } else {
      sqrt_inv <- tryCatch(
        pracma::sqrtm(gram)$Binv,
        error = function(e) NULL
      )

      if (!is.null(sqrt_inv) && all(is.finite(sqrt_inv))) {
        M_active <- M_active %*% sqrt_inv
      } else {
        col_scales <- sqrt(pmax(diag(gram), 0))
        valid <- which(col_scales > 0)
        if (length(valid) > 0) {
          M_active[, valid] <- sweep(M_active[, valid, drop = FALSE], 2, col_scales[valid], "/")
        }
      }
    }

    M[, active_cols] <- M_active
    M
  }

  X_train = as.matrix(data.frame(X_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  Y_train = as.matrix(data.frame(Y_train) %>% mutate_all(~replace_na(., mean(., na.rm = TRUE))))
  p1 <- dim(X_train)[2]
  p2 <- dim(Y_train)[2]
  p <- p1 + p2;
  n <- nrow(X_train)
  pp <- c(p1,p2)
  if (is.null(S)) {
    S <- cov(cbind(X_train, Y_train))
  }

  
  if (method.type=="FIT_SAR_BIC"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                           lambdaAseq=lambdax,
                           lambdaBseq=lambday,
                           max.iter=100,conv=10^-2,
                           selection.criterion=1,n.cv=kfolds)
    Uhat <- force_rank_columns(method$uhat, p1, rank, "method$uhat")
    Vhat <- force_rank_columns(method$vhat, p2, rank, "method$vhat")
    a_estimate = rbind(Uhat, Vhat)
    
  }
  if(method.type=="FIT_SAR_CV"){
    method<-SparseCCA(X=X_train,Y=Y_train,rank=rank,
                          lambdaAseq=lambdax,
                          lambdaBseq=lambday,
                          max.iter=100,conv=10^-2, selection.criterion=2, n.cv=kfolds)
    Uhat <- force_rank_columns(method$uhat, p1, rank, "method$uhat")
    Vhat <- force_rank_columns(method$vhat, p2, rank, "method$vhat")
    a_estimate = rbind(Uhat, Vhat)
    
  }
  if (method.type=="Witten_Perm"){
    Witten_Perm <- PMA::CCA.permute(x=X_train,z=Y_train,
                               typex="standard",typez="standard", 
                               penaltyxs =lambdax[which(lambdax < 1)],
                               penaltyzs = lambday[which(lambday < 1)],
                               nperms=50)
    method<-PMA::CCA(x=X_train, z=Y_train, typex="standard",typez="standard",K=rank,
                         penaltyx=Witten_Perm$bestpenaltyx,penaltyz=Witten_Perm$bestpenaltyz,trace=FALSE)
    Uhat <- force_rank_columns(method$u, p1, rank, "method$u")
    Vhat <- force_rank_columns(method$v, p2, rank, "method$v")
    a_estimate = rbind(Uhat, Vhat)
  }
  if(method.type=="Witten.CV"){
    Witten_CV<-Witten.CV(X=X_train,Y=Y_train,
                        rank=rank,
                        n.cv=kfolds,lambdax=lambdax[which(lambdax < 1)],
                         lambday=c(lambday[which(lambday < 1)]))
    method <-PMA::CCA(x=X_train,z=Y_train,typex="standard",typez="standard",
                 K=rank,penaltyx=Witten_CV$lambdax.opt,
                 penaltyz=Witten_CV$lambday.opt,trace=FALSE)
    Uhat <- force_rank_columns(method$u, p1, rank, "method$u")
    Vhat <- force_rank_columns(method$v, p2, rank, "method$v")
    a_estimate = rbind(Uhat, Vhat)
    
  }
  if(method.type=="Waaijenborg-Author"){
    method<-Waaijenborg(X=X_train,Y=Y_train,
                        lambdaxseq=lambdax,
                        lambdayseq=lambday,
                        rank=rank,selection.criterion=1)
    Vhat <- force_rank_columns(method$vhat, p1, rank, "method$vhat")
    Uhat <- force_rank_columns(method$uhat, p2, rank, "method$uhat")
    a_estimate = rbind(Vhat, Uhat)
    
  }
  if(method.type=="Waaijenborg-CV"){
    method<-Waaijenborg(X=X_train,
                        Y=Y_train,lambdaxseq=lambdax,
                        lambdayseq=lambday,
                        rank=rank, selection.criterion=2)
    Vhat <- force_rank_columns(method$vhat, p1, rank, "method$vhat")
    Uhat <- force_rank_columns(method$uhat, p2, rank, "method$uhat")
    a_estimate = rbind(Vhat, Uhat)
    
  }
  if(method.type=="SCCA_Parkhomenko"){
    method<- SCCA_Parkhomenko(x.data=X_train, y.data=Y_train, Krank=rank,
                              lambda.v.seq = lambdax[which(lambdax < 2)],
                              lambda.u.seq = lambday[which(lambday < 2)])
    Uhat <- force_rank_columns(method$uhat, p1, rank, "method$uhat")
    Vhat <- force_rank_columns(method$vhat, p2, rank, "method$vhat")
    a_estimate = rbind(Uhat, Vhat)
    
  }
  if(method.type=="Canonical Ridge-Author"){
    RCC_cv<-estim.regul_crossvalidation(X_train,Y_train,n.cv=5, 
                                        lambda1grid=lambdax[which(lambdax < 1)],
                                        lambda2grid=lambday[which(lambdax < 1)])
    method<-rcc(X_train,Y_train, RCC_cv$lambda1.optim, RCC_cv$lambda2.optim)
    Uhat <- force_rank_columns(method$xcoef, p1, rank, "method$xcoef")
    Vhat <- force_rank_columns(method$ycoef, p2, rank, "method$ycoef")
    a_estimate = rbind(Uhat, Vhat)
    
    
  }
  if(method.type=="Fantope"){
    afant <- Fantope(X_train, Y_train, r = rank)
    Uhat <- force_rank_columns(afant$u, p1, rank, "afant$u")
    Vhat <- force_rank_columns(afant$v, p2, rank, "afant$v")
    a_estimate = rbind(Uhat, Vhat)
  }
  if(method.type=="Chao"){
    a_estimate1 = scca_chao(X_train, Y_train, 
                                rho = 1,
                                lambda_max = .1, num_lambda = 20, r = rank, niter = 500,
                                nfold = 8, thresh = .01)
    
    Uhat <- force_rank_columns(a_estimate1$u, p1, rank, "a_estimate1$u")
    Vhat <- force_rank_columns(a_estimate1$v, p2, rank, "a_estimate1$v")
    a_estimate = rbind(Uhat, Vhat)
  }
  if(method.type=="SGCA"){
    idx1 <- 1:p1
    idx2 <- (p1+ 1):(p1 + p2)
    Mask <- matrix(0, (p1 + p2), (p1 + p2))
    Mask[idx1, idx1] <- matrix(1,p1, p1)
    Mask[idx2, idx2] <- matrix(1, p2, p2)
    sigma0hat <- S * Mask
    
    ag <- sgca_init(A=S, B=sigma0hat, rho=0.5 * sqrt(log(p + p2)/n),
                    K=rank,  maxiter=1000, trace=FALSE)
    ainit <- init_process(ag$Pi, rank) 
    a_estimate <- sgca_tgd(A=S, B=sigma0hat,rank,ainit,k=20,lambda = 0.01, eta=0.00025,convergence=1e-6,maxiter=12000, plot = TRUE)

    
    
  }
  if (!(exists("Uhat") && exists("Vhat"))) {
    a_estimate <- gca_to_cca(a_estimate, S, pp)
    Uhat <- force_rank_columns(a_estimate$u, p1, rank, "a_estimate$u")
    Vhat <- force_rank_columns(a_estimate$v, p2, rank, "a_estimate$v")
  }

  Sigmax <- cov(X_train)
  Sigmay <- cov(Y_train)
  Uhat <- normalize_loading_matrix(Uhat, Sigmax)
  Vhat <- normalize_loading_matrix(Vhat, Sigmay)

  return(list(u = Uhat, v = Vhat))
}
