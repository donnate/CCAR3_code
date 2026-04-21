# LIBRARIES
library(glmnet)
library(MASS)
library(pls)

force_column_matrix <- function(x) {
  x <- as.matrix(x)
  if (nrow(x) == 1 && ncol(x) > 1) {
    x <- t(x)
  }
  x
}

force_candidate_matrix <- function(x, expected_rows) {
  x <- as.matrix(x)
  if (nrow(x) != expected_rows && ncol(x) == expected_rows) {
    x <- t(x)
  }
  if (nrow(x) != expected_rows) {
    x <- matrix(as.numeric(x), nrow = expected_rows)
  }
  x
}

Waaijenborg<-function(X,Y,lambdaxseq=matrix(seq(from=1,to=0,by=-0.01),nrow=1),lambdayseq=matrix(seq(from=1,to=0,by=-0.01),nrow=1),rank,selection.criterion=1,n.cv=5,max.iter=20,conv=10^-3){
  ### Function to perform Sparse CCA based on Waaijenborg et al. (2008)
  # REFERENCE Waaijenborg et al. (2008), "Quantifying the Association between Gene Expressions and DNA-Markers by Penalized Canonical Correlation Analysis" in Statistical Applications in Genetics and Molecular Biology, Volume 7, Issue 1, Article 3
  
  ### INPUT
  # X                   : (nxp) data matrix
  # Y                   : (nxq) data matrix
  # lambdaxseq          : grid of sparsity parameters for lambdaA
  # lambdayseq          : grid of sparsity parameters for lambdaB
  # rank                : number of canonical vector pairs to extract
  # selection.criterion : 1 for cross-validation minimize difference between test and training sample correlation; 2 for maximizing test sample correlation
  # n.cv                : n.cv-fold cross-validation
  # max.iter            : maximum number of iterations
  # conv                : tolerance value for convergence
  
  ### OUTPUT
  # vhat                : (qxr) estimated canonical vectors correpsonding to the first data set
  # uhat                : (pxr) estimated canonical vectors correpsonding to the second data set
  # ksihat              : (nxr) estimated canonical variates corresponding to the first data set            
  # omegahat            : (nxr) estimated canonical variates corresponding to the second data set 
  # lambdax_FINAL       : value of the sparsity parameter lambdaA 
  # lambday_FINAL       : value of the sparsity parameter lambdaB
  # it                  : number of iterations         
  
  
  ### STORE RESULTS
  lambdaxseq=matrix(lambdaxseq,nrow=1)
  lambdayseq=matrix(lambdayseq,nrow=1)
  
  u_ALL<-matrix(NA,ncol=rank,nrow=ncol(Y))
  v_ALL<-matrix(NA,ncol=rank,nrow=ncol(X))
  ksi_ALL<-matrix(NA,ncol=rank,nrow=nrow(X))
  omega_ALL<-matrix(NA,ncol=rank,nrow=nrow(Y))
  cancors<-matrix(NA,ncol=rank,nrow=1)
  # Store sparsity parameters
  lambdax_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  lambday_ALL<-matrix(NA,ncol=rank,nrow=max.iter)
  iterations<-matrix(NA,ncol=rank,nrow=1)
  
  
  ### START CODE
  
  # Starting Values
  p<-ncol(Y)
  u.starting<-rep(1/p,p)
  q<-ncol(X)
  v.starting<-rep(1/q,q)
  

  # Sequentially extract the r canonical vector pairs
  for (i.r in 1:rank){
    
    if (i.r==1){ # for r=1: start from original data sets
      X_data<-X 
      Y_data<-Y
    } 
    
    # STARTING VALUES
    u.initial<-u.starting
    v.initial<-v.starting
    
    # STANDARDIZE DATA
    X_data_st <- as.matrix(stdize(X_data))
    if (ncol(X_data_st) != ncol(X_data)) {
      X_data_st <- matrix(X_data_st, ncol = ncol(X_data))
    }
    for (i.variable in 1:ncol(X_data)){
      if (is.na(apply(X_data_st,2,sum)[i.variable])==T) {
        X_data_st[,i.variable]<-0}    
    }
    
    Y_data_st <- as.matrix(stdize(Y_data))
    if (ncol(Y_data_st) != ncol(Y_data)) {
      Y_data_st <- matrix(Y_data_st, ncol = ncol(Y_data))
    }
    for (i.variable in 1:ncol(Y_data)){
      if (is.na(apply(Y_data_st,2,sum)[i.variable])==T) {
        Y_data_st[,i.variable]<-0}    
    }
    
    
    # CANONICAL VARIATES
    ksi <- force_column_matrix(X_data_st %*% v.initial)
    omega <- force_column_matrix(Y_data_st %*% u.initial)
    
    # INITIALIZE CONVERGENCE PARAMETERS
    it<-1
    diff.u<-conv*10
    diff.v<-conv*10
    
    
    while( (it<max.iter) & ((diff.u>conv) & (diff.v>conv))) {
      
      # Estimate canonical vector u
      U.hat <- force_candidate_matrix(apply(lambdayseq,2,UST,a=t(ksi)%*%Y_data_st), ncol(Y_data_st))
      keep_u <- which(colSums(abs(U.hat)) > 0)
      U.hat.reduced <- U.hat[, keep_u, drop = FALSE]
      lambday.reduced <- matrix(lambdayseq[keep_u], nrow=1)
      
      # Select tuning parameter
      if (selection.criterion==1){ #Minimize difference between test and training sample correlation 
            n = nrow(Y_data_st)
            n.cv.sample<-trunc(n/n.cv)
            whole.sample<-seq(1,n)
            Y.data<-Y_data_st
            ksi.data<-ksi
            cvscore<-matrix(NA,ncol=length(lambday.reduced),nrow=n.cv)
            for (i in 1:n.cv){
              testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
              training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
              Ycv = Y.data[training.sample, , drop = FALSE]
              ksicv = ksi[training.sample, , drop = FALSE]
              U.hat <- force_candidate_matrix(apply(lambday.reduced,2,UST,a=t(ksicv)%*%Ycv), ncol(Ycv))
              
              cvscore[i,]<-apply(U.hat,2,delta.correlation,Xtest=Y.data[testing.sample, , drop = FALSE],Ytest=ksi[testing.sample, , drop = FALSE],Xtrain=Y.data[training.sample, , drop = FALSE],Ytrain=ksi[training.sample, , drop = FALSE])  
            
            } 
            CVscore.mean<-apply(cvscore,2,mean)
            U.hat.final <- force_column_matrix(U.hat.reduced[,min(which.min(CVscore.mean)), drop = FALSE]) # OPTIMAL SOLUTION
            U.hat.final <- force_column_matrix(apply(U.hat.final,2,NORMALIZATION_UNIT))
            lambday_ALL[it,i.r]<-lambday.reduced[min(which.min(CVscore.mean))]
        
      } else {
            if (selection.criterion==2){ #Maximize test sample correlation
                  n = nrow(Y_data_st)
                  n.cv.sample<-trunc(n/n.cv)
                  whole.sample<-seq(1,n)
                  Y.data<-Y_data_st
                  ksi.data<-ksi
                  cvscore<-matrix(NA,ncol=length(lambday.reduced),nrow=n.cv)
                  for (i in 1:n.cv){
                    testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                    training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                    Ycv = Y.data[training.sample, , drop = FALSE]
                    ksicv = ksi[training.sample, , drop = FALSE]
                    U.hat <- force_candidate_matrix(apply(lambday.reduced,2,UST,a=t(ksicv)%*%Ycv), ncol(Ycv))
                    cvscore[i,]<-apply(U.hat,2,testsample.correlation,Xdata=Y.data[testing.sample, , drop = FALSE],yscore=ksi[testing.sample, , drop = FALSE] )  
                  } 
                  
                  CVscore.mean<-apply(cvscore,2,mean)
                  U.hat.final <- force_column_matrix(U.hat.reduced[,min(which.max(CVscore.mean)), drop = FALSE]) # OPTIMAL SOLUTION
                  U.hat.final <- force_column_matrix(apply(U.hat.final,2,NORMALIZATION_UNIT))
                  lambday_ALL[it,i.r]<-lambday.reduced[min(which.max(CVscore.mean))]
            } else {
              stop("selection.criterion needs to be equal to 1 (Difference test and training correlation) or 2 (Test correlation)")
            }
            
      }
            
      
      # Estimate canonical vector v
      V.hat <- force_candidate_matrix(apply(lambdaxseq,2,UST,a=t(omega)%*%X_data_st), ncol(X_data_st))
      keep_v <- which(colSums(abs(V.hat)) > 0)
      V.hat.reduced <- V.hat[, keep_v, drop = FALSE]
      lambdax.reduced <- matrix(lambdaxseq[keep_v],nrow=1)
      
      # Select tuning parameter
      if (selection.criterion==1){
              n = nrow(X_data_st)
              n.cv.sample<-trunc(n/n.cv)
              whole.sample<-seq(1,n)
              X.data<-X_data_st
              omega.data<-omega
              cvscore<-matrix(NA,ncol=length(lambdax.reduced),nrow=n.cv)
              for (i in 1:n.cv){
                testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                Xcv = X.data[training.sample, , drop = FALSE]
                omegacv = omega[training.sample, , drop = FALSE]
                V.hat <- force_candidate_matrix(apply(lambdax.reduced,2,UST,a=t(omegacv)%*%Xcv), ncol(Xcv))
                
                cvscore[i,]<-apply(V.hat,2,delta.correlation,Xtest=X.data[testing.sample, , drop = FALSE],Ytest=omega[testing.sample, , drop = FALSE],Xtrain=X.data[training.sample, , drop = FALSE],Ytrain=omega[training.sample, , drop = FALSE] )  
                
              } 
              CVscore.mean<-apply(cvscore,2,mean)
              V.hat.final <- force_column_matrix(V.hat.reduced[,min(which.min(CVscore.mean)), drop = FALSE]) # OPTIMAL SOLUTION
              V.hat.final <- force_column_matrix(apply(V.hat.final,2,NORMALIZATION_UNIT))
              lambdax_ALL[it,i.r]<-lambdax.reduced[min(which.min(CVscore.mean))]
        } else {
              if (selection.criterion==2){
                n = nrow(X_data_st)
                n.cv.sample<-trunc(n/n.cv)
                whole.sample<-seq(1,n)
                X.data<-X_data_st
                omega.data<-omega
                cvscore<-matrix(NA,ncol=length(lambdax.reduced),nrow=n.cv)
                for (i in 1:n.cv){
                  testing.sample<-whole.sample[((i-1)*n.cv.sample+1):(i*n.cv.sample)]
                  training.sample<-whole.sample[-(((i-1)*n.cv.sample+1):(i*n.cv.sample))]
                  Xcv = X.data[training.sample, , drop = FALSE]
                  omegacv = omega[training.sample, , drop = FALSE]
                  V.hat <- force_candidate_matrix(apply(lambdax.reduced,2,UST,a=t(omegacv)%*%Xcv), ncol(Xcv))
                  
                  cvscore[i,]<-apply(V.hat,2,testsample.correlation,Xdata=X.data[testing.sample, , drop = FALSE],yscore=omega[testing.sample, , drop = FALSE] )  
                } 
                
                CVscore.mean<-apply(cvscore,2,mean)
                V.hat.final <- force_column_matrix(V.hat.reduced[,min(which.max(CVscore.mean)), drop = FALSE]) # OPTIMAL SOLUTION
                V.hat.final <- force_column_matrix(apply(V.hat.final,2,NORMALIZATION_UNIT))
                lambdax_ALL[it,i.r]<-lambdax.reduced[min(which.min(CVscore.mean))]
              } else {
                stop("selection.criterion needs to be equal to 1 (Difference test and training correlation) or 2 (Test correlation)")
              }
        }
              
        # Convergence measures
        diff.u<-max(abs(u.initial-U.hat.final))
        diff.v<-max(abs(v.initial-V.hat.final))
             
        # Updated starting values 
        ksi <- force_column_matrix(X_data_st %*% V.hat.final)
        omega <- force_column_matrix(Y_data_st %*% U.hat.final)
        u.initial<-U.hat.final
        v.initial<-V.hat.final
      
        it<-it+1
      }
      
      # Number of ITERATIONS
      iterations[1,i.r]<-it
    
    # Final estimates of canonical vectors, variates and canonical correlation
    u_ALL[,i.r]<-U.hat.final
    v_ALL[,i.r]<-V.hat.final
    ksi_ALL[,i.r]<-ksi
    omega_ALL[,i.r]<-omega
    cancors[,i.r]<-abs(cor(ksi,omega))

    # Residual matrices
    gammahat <- ginv(crossprod(ksi)) %*% crossprod(ksi, X_data)
    thetahat <- ginv(crossprod(omega)) %*% crossprod(omega, Y_data)

    X_data <- X_data - ksi %*% gammahat
    Y_data <- Y_data - omega %*% thetahat

  }  # END FOR-LOOP
    
    
  ##OUTPUT 
  out<-list(uhat=u_ALL,vhat=v_ALL,ksihat=ksi_ALL,omegahat=omega_ALL,cancors=cancors,lambdax_FINAL=lambdax_ALL,lambday_FINAL=lambday_ALL,it=iterations)
  
}

UST<-function(a,U){ # Univariate Soft Thresholding
  matrix((abs(a)-U/2 + abs(abs(a)-U/2))/2*sign(a),ncol=1)
}

delta.correlation<-function(U,Xtest,Ytest,Xtrain,Ytrain){ # Cross-validation criterium
  if (all(U==0)) {
    return(0)
  } else {
    traincor=abs(cor(Xtrain%*%U,Ytrain))
    testcor=abs(cor(Xtest%*%U,Ytest))
    return(abs(traincor-testcor))
  }
  
}
