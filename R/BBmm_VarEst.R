VarEst <- function(y,m,p,X,Z,u,nRand,nComp,nRandComp,OLDall.sigma,OLDphi,q,maxiter){
  
  #Number of observations
  nObs <- length(y)
  
  #The 'weight' matrix
  S <- diag(c(p*(1-p)))
  
  #Initial values
  phi <- OLDphi
  #psi <- log(phi)
  all.sigma <- OLDall.sigma
  theta <- c(phi,all.sigma)
  oldtheta <- rep(Inf,1+nComp)
  
  while(max(abs((theta-oldtheta)/theta))>0.001){
    
    oldtheta <- theta
    
    all.sigma2 <- all.sigma^2
    
    d <- d. <- NULL
    for (i in 1:nComp){
      d <- c(d,rep(all.sigma2[i],nRandComp[i]))
      d. <- c(d.,rep(1/(all.sigma2[i]),nRandComp[i]))
    }
    D <- diag(d)
    D. <- diag(d.)
    
    oldbetaphisigma <- rep(Inf,q+1+nComp)
    betaphisigma <- c(beta,phi,all.sigma) 
    
    ### ESTIMATION OF PSI
    function.psi <- function(psi){
      
      # Compute j and v
      j <- NULL
      v <- NULL
      hpsi <- NULL
      for (l in 1:nObs){
        j1 <- 0
        j2 <- 0
        v1 <- 0
        v2 <- 0
        hpsi1 <- 0
        hpsi2 <- 0
        hpsi3 <- 0
        
        if (y[l]==0){
          j1 <- 0
          v1 <- 0
          hpsi1 <- 0
        }else{
          for (k in 0:(y[l]-1)){
            j1 <- j1+k*exp(psi)/(p[l]+k*exp(psi))^3
            v1 <- v1+1/(p[l]+k*exp(phi))^2
            hpsi1 <- hpsi1+k*exp(psi)/(p[l]+k*exp(psi))
          }
        }
        if (y[l]==m[l]){
          j2 <- 0
          v2 <- 0
          hpsi2 <- 0
        }else{
          for (k in 0:(m[l]-y[l]-1)){
            j2 <- j2+k*exp(psi)/(1-p[l]+k*exp(psi))^3
            v2 <- v2+1/(1-p[l]+k*exp(psi))^2
            hpsi2 <- hpsi2+k*exp(psi)/(1-p[l]+k*exp(psi))
          }
        }
        for (k in 0:(m[l]-1)){
          hpsi3 <- hpsi3+k*exp(psi)/(1+k*exp(psi))
        }
        j <- c(j,-j1-j2)
        v <- c(v,v1+v2)
        hpsi <- c(hpsi,hpsi1+hpsi2-hpsi3)
      }
      
      V <- diag(c(v))
      J <- diag(c(j))
      
      W <- S%*%V%*%S
      WJ <- S%*%J%*%S
      #Compute H (Expected Hessian Matrix)
      H1 <- cbind(t(X)%*%W%*%X,t(X)%*%W%*%Z)
      H1psi <- cbind(t(X)%*%WJ%*%X,t(X)%*%WJ%*%Z)
      H2 <- cbind(t(Z)%*%W%*%X,t(Z)%*%W%*%Z+D.)
      H2psi <- cbind(t(Z)%*%WJ%*%X,t(Z)%*%WJ%*%Z)
      H <- rbind(H1,H2)
      Hpsi <- rbind(H1psi,H2psi)
      
      trace <- sum(diag(2*solve(H)%*%Hpsi))
      
      out <- sum(hpsi)-(1/2)*trace
      
      return(out)
      
    }
    
    mle.psi <- uniroot(function.psi,lower=-5,upper=3,tol=0.001)
    psi <- mle.psi$root
    phi <- exp(psi)
    ##########
    
    #ESTIMATION OF SIGMA
    
    # Compute v because it is the same for all sigma score equation
    v <- NULL
    for (l in 1:nObs){
      v1 <- 0
      v2 <- 0
      if (y[l]==0){
        v1 <- 0
      }else{
        for (k in 0:(y[l]-1)){
          v1 <- v1+1/(p[l]+k*phi)^2
        }
      }
      if (y[l]==m[l]){
        v2 <- 0
      }else{
        for (k in 0:(m[l]-y[l]-1)){
          v2 <- v2+1/(1-p[l]+k*phi)^2
        }
      }
      v <- c(v,v1+v2)
    }
    
    V <- diag(c(v))
    
    #We get the score equation terms where are the same for all sigma
    W <- S%*%V%*%S
    K. <- t(Z)%*%W%*%Z-t(Z)%*%W%*%X%*%solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%Z
    
    
    function.sigma <- function(sigma){
      
      nRand.previous <- 0
      d. <- NULL
      out <- NULL
      
      for (i in 1:nComp){      
        d. <- c(d.,rep(1/(sigma[i])^2,nRandComp[i]))
      }
      
      #Compute the trace
      D. <- diag(c(d.))
      K <- solve(K.+D.)
      trace <- diag(K)
      
      for(i in 1:nComp){
        out <- c(out,-nRandComp[i]*(sigma[i])^2+t(u[seq(nRand.previous+1,nRand.previous+nRandComp[i])])%*%u[seq(nRand.previous+1,nRand.previous+nRandComp[i])]+sum(trace[seq(nRand.previous+1,nRand.previous+nRandComp[i])]))
        
        nRand.previous <- nRand.previous+nRandComp[i]
      }
      out
    }
    
    mle.sigma <- multiroot(function.sigma,start=all.sigma)
    all.sigma <- mle.sigma$root
    ###################
    
    theta <- c(phi,all.sigma)
  }
  
  #The derivate in the point, the variance:
  psi.var <- -1/grad(function.psi,psi)
  sigma.var <- -1/grad(function.sigma,all.sigma)
  
  out <- list(phi=phi,all.sigma=all.sigma,psi=psi,psi.var=psi.var,all.sigma.var=sigma.var)
  return(out)
}