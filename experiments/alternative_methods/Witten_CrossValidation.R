#' Sparse CCA by Witten and Tibshirani (2009)
#'
#' @param X Matrix of predictors (n x p)
#' @param Y Matrix of responses (n x q)
#' @param n.cv Number of cross-validation folds (default is 5)
#' @param lambdax Vector of sparsity parameters for X (default is a sequence from 0 to 1 with step 0.1)
#' @param lambday Vector of sparsity parameters for Y (default is a sequence from 0 to 1 with step 0.1)
#' @param standardize Standardize (center and scale) the data matrices X and Y (default is TRUE) before analysis
#' @param rank Number of canonical components to extract
#'
#' @return the appropriate levels of regularisation
#' @export

Witten.CV<-function(X,Y,n.cv=5,
                    rank=1,
                    lambdax=matrix(seq(from=0,to=1,by=0.1),nrow=1),
                    lambday=matrix(seq(from=0,to=1,by=0.1),nrow=1),
                    standardize = FALSE){


  n = nrow(X)
  whole.sample <- seq_len(n)
  fold_id <- sample(rep(seq_len(n.cv), length.out = n))
  lambdax=matrix(lambdax,nrow=1)
  lambday=matrix(lambday,nrow=1)

  cvscore<-array(NA,c(length(lambday),length(lambdax),n.cv)) #lambdax in columns, lambday in rows


  for (i in seq_len(n.cv)){
    testing.sample <- whole.sample[fold_id == i]
    training.sample <- whole.sample[fold_id != i]
    Xcv = as.matrix(X[training.sample, , drop = FALSE])
    Ycv = as.matrix(Y[training.sample, , drop = FALSE])
    Xtest = as.matrix(X[testing.sample, , drop = FALSE])
    Ytest = as.matrix(Y[testing.sample, , drop = FALSE])

    if (standardize){
      Xcv_scaled <- scale(Xcv, center = TRUE, scale = TRUE)
      Ycv_scaled <- scale(Ycv, center = TRUE, scale = TRUE)

      Xcv_center <- attr(Xcv_scaled, "scaled:center")
      Xcv_scale <- attr(Xcv_scaled, "scaled:scale")
      Ycv_center <- attr(Ycv_scaled, "scaled:center")
      Ycv_scale <- attr(Ycv_scaled, "scaled:scale")

      Xcv_scale[!is.finite(Xcv_scale) | Xcv_scale == 0] <- 1
      Ycv_scale[!is.finite(Ycv_scale) | Ycv_scale == 0] <- 1

      Xcv <- Xcv_scaled
      Ycv <- Ycv_scaled
      Xtest <- sweep(sweep(Xtest, 2, Xcv_center, "-"), 2, Xcv_scale, "/")
      Ytest <- sweep(sweep(Ytest, 2, Ycv_center, "-"), 2, Ycv_scale, "/")
    }

    cvscore[,,i]<-apply(lambdax,2, Witten.cv.lambdax, Xtrain=Xcv,Ytrain=Ycv,Xtest=Xtest,Ytest=Ytest,lambday=lambday, r=rank)
  }
  cvscore.vec<-c(apply(cvscore,c(1,2),mean, na.rm = TRUE))


  LAMBDAX<-c(matrix(rep(lambdax,length(lambday)),byrow=T,nrow=length(lambday),ncol=length(lambdax)))
  LAMBDAY<-c(matrix(rep(lambday,length(lambdax)),byrow=F,nrow=length(lambday),ncol=length(lambdax)))

  lambdax.opt<-LAMBDAX[which.max(cvscore.vec)]
  lambday.opt<-LAMBDAY[which.max(cvscore.vec)]

  X_fit <- as.matrix(X)
  Y_fit <- as.matrix(Y)
  if (standardize){
    X_fit <- scale(X_fit, center = TRUE, scale = TRUE)
    Y_fit <- scale(Y_fit, center = TRUE, scale = TRUE)
  }

  fit.witten <- PMA::CCA(
    x = X_fit,
    z = Y_fit,
    typex = "standard",
    typez = "standard",
    K = rank,
    penaltyx = lambdax.opt,
    penaltyz = lambday.opt,
    trace = FALSE
  )

  list(
    lambdax.opt = lambdax.opt,
    lambday.opt = lambday.opt,
    fit = fit.witten,
    cvscore = cvscore
  )
}


Witten.cv.lambdax<-function(U,Xtrain,Ytrain,Xtest,Ytest,lambday, r=1){ #AUXILIARY FUNCTION
  testcorrelations<-apply(lambday,2,Witten.cv.lambday,lambdaxfixed=U,Xtrain=Xtrain,Ytrain=Ytrain,Xtest=Xtest,Ytest=Ytest, r=r)
  return(testcorrelations)
}

Witten.cv.lambday<-function(V,Xtrain,Ytrain,Xtest,Ytest,
                            lambdaxfixed, r=1){ #AUXILIARY FUNCTION
  #print(lambdaxfixed)
  Fit.Witten<-PMA::CCA(x=Xtrain,z=Ytrain,typex="standard",typez="standard",K=r,penaltyx=lambdaxfixed,penaltyz=V,trace=F)
  test_cor <- suppressWarnings(cor(Xtest %*% Fit.Witten$u, Ytest %*% Fit.Witten$v))
  if (is.null(dim(test_cor))) {
    return(abs(test_cor[[1]]))
  }
  return(sum(diag(abs(test_cor))))
}
