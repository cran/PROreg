VarEst <- function(y,m,p,X,Z,u,nRand,nComp,nRandComp,OLDall.sigma,OLDphi,q,nDim){
  
  #Number of observations
  nObs <- length(y)/nDim
  
  #The wiegth matrix S
  s <- c(p*(1-p))
  
  # Defining the variance of the random effects
  d <- d. <- NULL
  for (i in 1:nComp){
    d <- c(d,rep(OLDall.sigma[i]^2,nRandComp[i]))
    d. <- c(d.,rep(1/(OLDall.sigma[i])^2,nRandComp[i]))
  }
  D <- diag(d)
  D. <- diag(d.)
  
  ### ESTIMATION OF PSI
  function.psi <- function(psi){
    
    # Compute j and v
    v <- NULL
    hpsi <- NULL
    for (i in 1:nDim){
      for (l in ((i-1)*nObs+1):(i*nObs)){
        v1 <- 0
        v2 <- 0
        hpsi1 <- 0
        hpsi2 <- 0
        hpsi3 <- 0
        
        if (y[l]==0){
          v1 <- 0
          hpsi1 <- 0
        }else{
          for (k in 0:(y[l]-1)){
            v1 <- v1+1/(p[l]+k*exp(psi[i]))^2
            hpsi1 <- hpsi1+log(p[l]+k*exp(psi[i]))
          }
        }
        if (y[l]==m[l]){
          v2 <- 0
          hpsi2 <- 0
        }else{
          for (k in 0:(m[l]-y[l]-1)){
            v2 <- v2+1/(1-p[l]+k*exp(psi[i]))^2
            hpsi2 <- hpsi2+log(1-p[l]+k*exp(psi[i]))
          }
        }
        for (k in 0:(m[l]-1)){
          hpsi3 <- hpsi3+log(1+k*exp(psi[i]))
        }
        v <- c(v,v1+v2)
        hpsi <- c(hpsi,hpsi1+hpsi2-hpsi3)
      }
    }
    
    svs <- s*v*s
    W <- diag(c(svs))
    
    W[which(W<10^(-6))] <- 10^(-6)
    
    #Compute H (Expected Hessian Matrix)
    H1 <- cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z)
    H2 <- cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+D.)
    H <- rbind(H1,H2)
    eig <- svd(H)$d
    H.. <- sum(log(eig))
    out <- -(sum(hpsi)-(1/2)*H..)
    
    return(out)
    
  }
  
   OLDpsi <- log(OLDphi)
  mle.psi <- optim(OLDpsi,function.psi,hessian=TRUE,method="L-BFGS-B",lower=rep(-10,nDim),upper=rep(4,nDim),control=list(pgtol=1e-3))
  psi <- mle.psi$par
  
  if (min(psi)==-10){
    psi[order(psi)[1]] <- runif(1,-10,-8)
    cat("-Phi estimation has not converged-\n")
  }
  phi <- exp(psi)
  
  ##########
  
  #ESTIMATION OF SIGMA
  
  # Compute v because it is the same for all sigma score equation
  v <- NULL
  for (i in 1:nDim){
    for (l in ((i-1)*nObs+1):(i*nObs)){
      v1 <- 0
      v2 <- 0
      if (y[l]==0){
        v1 <- 0
      }else{
        for (k in 0:(y[l]-1)){
          v1 <- v1+1/(p[l]+k*phi[i])^2
        }
      }
      if (y[l]==m[l]){
        v2 <- 0
      }else{
        for (k in 0:(m[l]-y[l]-1)){
          v2 <- v2+1/(1-p[l]+k*phi[i])^2
        }
      }
      v <- c(v,v1+v2)
    }
  }
  
  
  #We get the score equation terms where are the same for all sigma
  SVS <- diag(s*v*s)
  #SVS[which(SVS<10^(-6))] <- 10^(-6) Cambian considerablemente las estimaciones! X
  P. <- t(Z)%*%SVS%*%Z-t(Z)%*%SVS%*%X%*%solve(t(X)%*%SVS%*%X)%*%t(X)%*%SVS%*%Z
  #P.[which(P.<10^(-6))] <- 10^(-6) Cambian considerablemente las estimaciones! X
  
  
  function.sigma <- function(sigma){
    
    nRand.previous <- 0
    d. <- NULL
    out <- NULL
    
    for (i in 1:nComp){      
      d. <- c(d.,rep(1/(sigma[i])^2,nRandComp[i]))
    }
    
    #Compute the trace
    D. <- diag(c(d.))
    P <- solve(P.+D.)
    trace <- diag(P)
    
    for(i in 1:nComp){
      
      # SIMPLIFIED
      out <- c(out,-nRandComp[i]*(sigma[i])^2+t(u[seq(nRand.previous+1,nRand.previous+nRandComp[i])])%*%u[seq(nRand.previous+1,nRand.previous+nRandComp[i])]+sum(trace[seq(nRand.previous+1,nRand.previous+nRandComp[i])]))
      
      # ORIGINAL
      #out <- c(out,-nRandComp[i]/(sigma[i])+t(u[seq(nRand.previous+1,nRand.previous+nRandComp[i])])%*%u[seq(nRand.previous+1,nRand.previous+nRandComp[i])]/(sigma[i])^3+sum(trace[seq(nRand.previous+1,nRand.previous+nRandComp[i])])/(sigma[i])^3)
      
      nRand.previous <- nRand.previous+nRandComp[i]
    }
    out
  }
  
  mle.sigma <- multiroot(function.sigma,start=OLDall.sigma,atol=0.01,ctol=0.01,rtol=0.01)
  all.sigma <- mle.sigma$root
  
  ############################################
  # VARIANCES
  ############################################
  #Variance of psi
  psi.var <- 1/diag(mle.psi$hessian)
  # Variance of sigma: in the uniroot equation the derivative of the score equation has been simplifyed, therefore we define the function without any simplification to calculate the derivative (2nd derivative of the score equation). We have simplifyed the derivative of the score equation because if not sometimes it gets errors as there are multipli solutions of sigma (0).
  function.sigma <- function(sigma){
    nRand.previous <- 0
    d. <- NULL
    out <- NULL
    for (i in 1:nComp){      
      d. <- c(d.,rep(1/(sigma[i])^2,nRandComp[i]))
    }
    #Compute the trace
    D. <- diag(c(d.))
    P <- solve(P.+D.)
    trace <- diag(P)
    for(i in 1:nComp){
      out <- c(out,-nRandComp[i]/(sigma[i])+t(u[seq(nRand.previous+1,nRand.previous+nRandComp[i])])%*%u[seq(nRand.previous+1,nRand.previous+nRandComp[i])]/(sigma[i])^3+sum(trace[seq(nRand.previous+1,nRand.previous+nRandComp[i])])/(sigma[i])^3)
      nRand.previous <- nRand.previous+nRandComp[i]
    }
    out
  }
  sigma.var <- -1/grad(function.sigma,all.sigma)
  
  #################################################
  # OUTPUT
  #################################################
  
  out <- list(phi=phi,all.sigma=all.sigma,psi=psi,psi.var=psi.var,all.sigma.var=sigma.var)
  return(out)
}
