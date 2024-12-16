library(dplyr)
library(tidyr)
library(Matrix)
library(glmnet)
library(gglasso)
library(rrpack)
library(foreach)
library(doParallel)

CCA_rrr = function(X, Y, Sx=NULL, Sy=NULL,
                   lambda = 0, Kx=NULL, r, highdim=TRUE, 
                   lambda_Kx=0, solver="ADMM",
                   LW_Sy = FALSE,
                   do.scale = TRUE,
                   rho=1,
                   niter=1e4,
                   thresh = 1e-4, verbose=FALSE){
  # solve RRR: ||Y-XB|| + tr(Bt K B)
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  if (do.scale){
      X <- scale(X)
      Y <- scale(Y)
  }else{
      X <- scale(X, scale=FALSE)
      Y <- scale(Y, scale=FALSE)
  }
  
  if ( n <  min(q, p)) {
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (is.null(Sx)){
    Sx = t(X) %*% X /n
  }
  if (is.null(Sy)){
    Sy = t(Y) %*% Y /n
    if (LW_Sy){
      lw_cov <- corpcor::cov.shrink(Y)
      Sy <- as.matrix(lw_cov)
    }
  }

  svd_Sy = svd(Sy)
  sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sy$v)
  tilde_Y = Y %*% sqrt_inv_Sy
  
  if(!is.null(Kx)){
    Sx_tot = Sx + lambda_Kx * as.matrix(Kx) 
  }else{
    Sx_tot = Sx
  }
  Sxy = t(X) %*% tilde_Y/ n 
  if(!highdim){
    B_OLS = solve(Sx_tot) %*% Sxy
    svd_Sx = svd(Sx)
    sqrt_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sqrt_inv_Sx = svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sx$u)
    sol = svd(sqrt_Sx %*% B_OLS)
    V <- sqrt_inv_Sy %*% sol$v[, 1:r]
    U <- sqrt_inv_Sx %*% sol$u[, 1:r]
    if(verbose){
      print(t(U) %*% Sx %*% U)
      print(t(V) %*% Sy %*% V)
    }
  } else {
      if (solver == "CVX"){
        print("Using CVXR")
        ### Use CVXR
        B <- Variable(p, q)
        objective <- Minimize(1/n * sum_squares(tilde_Y - X %*% B) + lambda * sum(norm2(B, axis=1)))
        problem <- Problem(objective)
        result <- solve(problem)
        B_opt <- result$getValue(B)
      }else{
        if (solver == "ADMM") {
          print("Using ADMM")
          U <- matrix(0, p, q)
          Z <- matrix(0, p, q)
          prod_xy = t(X) %*% tilde_Y/ n
          invSx = solve(Sx_tot + rho * diag(rep(1, p)))
          for (i in 1:niter){
            #print(paste0("iter is ", i))
            Uold = U
            Zold = Z
            B = invSx %*% (prod_xy  + rho * (Z - U))
            Bold = B
            Z = B + U
            norm_col = sapply(1:nrow(Z), function(i){sqrt(sum(Z[i,]^2))})
            index_0 = which(norm_col < lambda/rho)
            index_pos = which(norm_col  > lambda/rho)
            if(length(index_0)>0){
              Z[index_0,] = 0
            }
            if (length(index_pos) > 1){
              Z[index_pos,] = diag(sapply(index_pos, function(x){ 1- (lambda/rho)/norm_col[x]})) %*% Z[index_pos,] 
            }else{
	      if (length(index_pos) == 1){
                Z[index_pos,] = ( 1- (lambda/rho)/norm_col[index_pos]) * Z[index_pos,]
              }
	    }
            U = U + B - Z
            if(verbose){
               print(c("ADMM iter", i, norm(Z - B), norm(Zold - Z), norm(Uold - U)))
            }
            if (max(c(norm(Z - B)/ sqrt(p), norm(Zold - Z)/ sqrt(p))) <thresh){
              break
            }
          }
          B_opt = B
        }else{
          print("Using Solver")
          test <- cv.srrr( tilde_Y, X, nrank = r,
                          method ="glasso",
                          nfold = 2, norder = NULL,
                          A0 = NULL,   V0 = NULL,
                          modstr = list("lamA" = rep(lambda, 10),
                                        "nlam" = 10))
          B_opt <- test$coef # test$U  %*% test$D %*% t(test$V)
        }
      }
      B_opt[which(abs(B_opt)<1e-5)] = 0
      I = which(apply(B_opt^2, 1, sum) > 0)
      if (length(I) >r-1){
        svd_Sx <- svd(Sx[I, I])
        sqrt_Sx <- svd_Sx$u %*% diag(sapply(svd_Sx$d, function(x){ifelse(x > 1e-4, sqrt(x), 0)}))  %*% t(svd_Sx$u)
        sol <- svd(sqrt_Sx %*% B_opt[I, ] )
        V <- sqrt_inv_Sy %*% sol$v[, 1:r]
        inv_D <- diag(sapply(1:r, FUN=function(x){ifelse(sol$d[x]<1e-4, 0, 1/sol$d[x])}))
        U <- B_opt %*% sol$v[, 1:r] %*% inv_D ### = U\lambda
      }else{
        U <- matrix(0, p, r)
        V <- matrix(0, q, r)
      }
      if(verbose){
        print(t(U) %*% Sx %*% U)
        print(t(V) %*% Sy %*% V)
      }
    }
    
  
  loss <- mean((Y %*% V - X %*% U)^2)
  correlation <- cor(Y %*% V, X %*% U)
  return(list(U = U, V = V, loss = loss,
              total_corr = correlation,
              cor = sapply(1:r, function(i){cov(X %*% U[,i], Y %*% V[,i])})))
}


CCA_rrr.CV<- function(X, Y, 
                      r=2, Kx = NULL, lambda_Kx = 0,
                      param_lambda=10^seq(-3, 1.5, length.out = 100),
                      kfolds=14,
                      solver="rrr",
                      parallelize = FALSE,
                      LW_Sy = FALSE,
                      do.scale=TRUE,
                      rho=1,
                      niter=1e4,
                      thresh = 1e-4
){
  print(paste0("Solver is ", solver))
  n = nrow(X)
  p = ncol(X)
  q = ncol(Y)
  if(q >n){
    X_temp <- X
    X <- Y
    Y <- X_temp
  }
  if (do.scale){
     X <- scale(X)
     Y <- scale(Y)
  }else{
    X <- scale(X, scale = FALSE)
    Y <- scale(Y, scale = FALSE)
  }
  Sx = t(X) %*% X /n
  Sy = t(Y) %*% Y /n
  if (LW_Sy){
    lw_cov <- corpcor::cov.shrink(Y)
    Sy <- as.matrix(lw_cov)
  }

  if ( n <  min(q, p)){
    print("Warning!!!! Both X and Y are high dimensional, method may fail")
  }
  if (solver %in% c("CVX", "CVXR", "ADMM")){
    if (parallelize){
      no_cores <- detectCores() - 2  # Save three cores for system processes
      registerDoParallel(cores=no_cores)
      resultsx <- foreach(lambda=param_lambda, .combine=rbind, .packages=c('CVXR',
                                                                           'Matrix')) %dopar% {
                                                                             rmse <- CCA_rrr.folds(X, Y, Sx=NULL, Sy=NULL,
                                                                                                   kfolds=kfolds, init,
                                                                                                   lambda=lambda, r=r, solver=solver,
                                                                                                    Kx = Kx, lambda_Kx =lambda_Kx,
                                                                                                   do.scale=do.scale,
                                                                                                    rho=rho,
                                                                                                    niter=niter,
                                                                                                   thresh = thresh)
                                                                             data.frame(lambda=lambda, rmse=rmse)
                                                                           }
    }else {
      resultsx <- expand.grid(lambda = param_lambda) %>%
        dplyr::mutate(rmse = purrr::map_dbl(lambda, 
                              ~CCA_rrr.folds(X, Y, Sx=NULL, Sy=NULL,
                              kfolds=kfolds,
                              lambda=.x, r=r, 
                              solver=solver,
                              Kx = Kx, lambda_Kx =lambda_Kx,
                              do.scale=do.scale, rho=rho,
                              niter=niter,
                              thresh = thresh)))
    }
    
    resultsx$rmse[which(is.na(resultsx$rmse))] = 1e8
    resultsx$rmse[which((resultsx$rmse) ==0)] = 1e8
    resultsx = resultsx %>% filter(rmse > 1e-5) 
    opt_lambda <- resultsx$lambda[which.min(resultsx$rmse)]
    #print(c("selected", opt_lambda))
    rmse <- resultsx$rmse
    if(is.na(opt_lambda) | is.null(opt_lambda)){
      opt_lambda = 0.1
    }
    final <- CCA_rrr(X, Y, Sx = NULL, Sy = NULL,
                     highdim = TRUE,
                     lambda = opt_lambda, Kx = Kx,
                     r = r, solver = solver,
                     lambda_Kx=lambda_Kx,
                     do.scale = do.scale,
                     LW_Sy =  LW_Sy, 
                     rho=rho,
                     niter=niter, thresh=thresh)
    #print(resultsx)
  }else{
    print("Using rrr CV")
    svd_Sy = svd(Sy)
    sqrt_inv_Sy = svd_Sy$u %*% diag(sapply(svd_Sy$d, function(x){ifelse(x > 1e-4, 1/sqrt(x), 0)}))  %*% t(svd_Sy$u)
    tilde_Y = Y %*% sqrt_inv_Sy
    print("Using rrr solver")
    test <- cv.srrr(tilde_Y, X, nrank = r,
                    method ="glasso",
                    nfold = kfolds,
                    norder = NULL,
                    A0 = NULL,
                    V0 = NULL,
                    modstr = list("lamA" = param_lambda, 
                                  "nlam" = length(param_lambda)))
    opt_lambda <- test$lambda[test$minid]
    if(is.na(opt_lambda) | is.null(opt_lambda)){
      opt_lambda = 0.1
    }
    rmse = test$cv.path
    final <-CCA_rrr(X, Y, Sx = NULL, Sy=NULL, 
                    lambda =opt_lambda, Kx=Kx,
                    r = r, highdim=TRUE, 
                    solver = solver,
                    lambda_Kx=lambda_Kx,
                    do.scale=do.scale,
                    rho=rho,
                    LW_Sy =  LW_Sy,
                    niter=niter,
                    thresh = thresh)
    resultsx = data.frame("lambda" = test$lambda,
                          "rmse" = test$cv.path[1:length(param_lambda)])
    final <- list(U = final$U,
                  V = sqrt_inv_Sy %*% final$V)
  }
  return(list( ufinal = final$U, 
               vfinal = final$V,
               lambda=opt_lambda,
               resultsx=resultsx,
               rmse = rmse,
               cor = sapply(1:r, function(i){cov(X %*% final$U[,i], Y %*% final$V[,i])})
  ))
}


CCA_rrr.folds<- function(X, Y, Sx, Sy, kfolds=5,
                         lambda=0.01,
                         r=2, Kx = NULL,
                         lambda_Kx = 0,
                         do.scale=TRUE,
                         solver = "CVX",
                         rho=1,
                         LW_Sy = TRUE,
                         niter=1e4,
                         thresh = 1e-4) {
  # define empty vector to store results
  folds <- createFolds(1:nrow(Y), k = kfolds, list = TRUE, returnTrain = FALSE)
  rmse <- rep(1e8, kfolds)
  p1 <- dim(X)[2]
  p2 <- dim(Y)[2]

  #init <- gca_to_cca(ainit, S0, pp)
  # loop over folds
  no_cores <- detectCores() - 2  # Save three cores for system processes
  registerDoParallel(cores=no_cores) 
  
  rmse <- foreach(i=seq_along(folds), .combine=c, .packages=c('CVXR', 'Matrix')) %dopar% {
    # ... same code inside your for loop ...
    X_train <- X[-folds[[i]], ]
    Y_train <- Y[-folds[[i]], ]
    X_val <- X[folds[[i]], ]
    Y_val <- Y[folds[[i]], ]
    tryCatch(
      {
        final = CCA_rrr(X_train, Y_train, Sx=NULL,
                        Sy =NULL, highdim = TRUE,
                        lambda = lambda, Kx=Kx, r=r, 
                        solver = solver,
                        lambda_Kx=lambda_Kx,
                        LW_Sy = LW_Sy,
                        do.scale=do.scale,
                        rho=rho,
                        niter=niter,
                        thresh = thresh)
        # make predictions on validation data
        # compute RMSE on validation data
        return(mean((X_val %*% final$U - Y_val%*% final$V)^2))
        #print(rmse)
      },
      error = function(e) {
        #    # If an error occurs, assign NA to the result
        print("An error has occured")
        return(NA)
      })
  }
  #print(c(lambda, rmse))
  # return mean RMSE across folds
  if (mean(is.na(rmse)) == 1){
    return(1e8)
  }else{
    return(mean(rmse, na.rm=TRUE))
  }
}





