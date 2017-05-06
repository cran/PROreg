EffectsEst.multiroot <- function(y,m,beta,u,p,phi,D.,X,Z,maxiter){
  
  #number of observations
  nObs <- length(y)
  
  #number of fixed effects
  q <- length(beta)
  
  #number of random effects
  nRand <- length(u)
  
  #Initial values
  oldbeta <- beta
  oldu <- u
  
  #Defining the score equation functions
  EffectsEst.function <- function(x){
    beta <- x[1:q]
    u <- x[(q+1):length(x)]
    p <- 1/(1+exp(-(X%*%beta+Z%*%u)))
    S <- diag(c(p*(1-p)))
    
    t <- NULL
    for (j in 1:nObs){
      t1 <- 0
      t2 <- 0
      if (y[j]==0){
        t1 <- 0
      }else{
        for (k in 0:(y[j]-1)){
          t1 <- t1+1/(p[j]+k*phi)
        }
      }
      if (y[j]==m[j]){
        t2 <- 0
      }else{
        for (k in 0:(m[j]-y[j]-1)){
          t2 <- t2+1/(1-p[j]+k*phi)
        }
      }
      t <- c(t,t1-t2)
    }
    c(F1=t(X)%*%S%*%t,F2=t(Z)%*%S%*%t-D.%*%u)
  }
  
  #Geting the mle
  mle <- multiroot(f=EffectsEst.function,start=c(oldbeta,oldu),ctol=0.001,useFortran = TRUE,maxiter=maxiter)
  
  #Looking if convergence was reached
  iter <- mle$iter
  if (iter==maxiter){
    print("The maximum number of iterations was reached without convergence")
    conv <- "no"
    out <- list(conv=conv)
    return(out)
  }
  
  #Calculing the values
  beta <- mle$root[1:q]
  names(beta) <- colnames(X)
  u <- mle$root[(q+1):(nRand+q)]

  #Convergence
  conv <- "yes"
  
  #Random components names
  names(u) <- seq(1:length(u))
  
  #Variances
  eta <- X%*%beta+Z%*%u
  p <- 1/(1+exp(-eta))
  S <- diag(c(p*(1-p)))
  t <- NULL
  v <- NULL
  for (j in 1:nObs){
    t1 <- 0
    t2 <- 0
    v1 <- 0
    v2 <- 0
    if (y[j]==0){
      t1 <- 0
      v1 <- 0
    }else{
      for (k in 0:(y[j]-1)){
        t1 <- t1+1/(p[j]+k*phi)
        v1 <- v1+1/((p[j]+k*phi)^2)
      }
    }
    if (y[j]==m[j]){
      t2 <- 0
      v2 <- 0
    }else{
      for (k in 0:(m[j]-y[j]-1)){
        t2 <- t2+1/(1-p[j]+k*phi)
        v2 <- v2+1/((1-p[j]+k*phi)^2)
      }
    }
    t <- c(t,t1-t2)
    v <- c(v,v1+v2)
  }
  
  V <- diag(c(v))
  vcov.fixed <- solve(t(X)%*%S%*%V%*%S%*%X)
  var.random <- diag(t(Z)%*%S%*%V%*%S%*%Z+D.)
  names(var.random) <- seq(1:length(u))
  
  
  #return
  out <- list(fixed.est=beta,random.est=u,vcov.fixed=vcov.fixed,var.random=var.random,iter.fixrand=iter,conv.fixrand=conv)
  
}
