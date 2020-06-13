EffectsEst.BBDelta <- function(y,m,beta,u,p,phi,D.,X,Z,maxiter){
  
  #number of observations
  nObs <- length(y)
  
  #number of fixed effects
  q <- length(beta)
  
  #number of random effects
  nRand <- length(u)
  
  #Initial values
  iter <- 0
  eta <- X%*%beta+Z%*%u
  oldeta <- rep(Inf,nObs)
  
  while (max(abs((eta-oldeta)/eta))>0.001){
    
    oldeta <- eta
    
    oldbeta <- rep(Inf,q)
    beta.iter <- 1
    while (max(abs((beta-oldbeta)/beta))>0.001){
      
      oldbeta <- beta
      p <- 1/(1+exp(-(X%*%beta+Z%*%u)))
      
      # Compute S
      s <- p*(1-p)
      S <- diag(c(s))
      
      # Compute t and v
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
      
      V <- diag(v)
      
      sv <- s*v
      
      # Updating beta
      etab <- X%*%beta+(1/sv)*t
      Hb <- solve(t(X)%*%S%*%V%*%S%*%X)
      beta <- Hb%*%t(X)%*%as.vector(s*v*s*etab)
      
      
      if (beta.iter==1000){
        print("fixed effects estimation has failed")
        conv <- "no"
        out <- list(conv=conv)
        return(out)      
        }
      
      
      beta.iter <- beta.iter +1
    }
    
    oldu <- rep(Inf,nRand)
    u.iter <- 1
    while (max(abs((u-oldu)/u))>0.1){
      
      oldu <- u
      p <- 1/(1+exp(-(X%*%beta+Z%*%u)))
      
      # Compute S
      s <- p*(1-p)
      S <- diag(c(s))
      
      # Compute t and v
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
      
      V <- diag(v)
      
      sv <- s*v
      
      # Updating u
      etau <- Z%*%u+(1/sv)*t
      Hu <- solve(t(Z)%*%S%*%V%*%S%*%Z+D.)
      u <- Hu%*%t(Z)%*%as.vector(s*v*s*etau)
      
      if (u.iter==1000){
        print("random effects estimation has failed")
        conv <- "no"
        out <- list(conv=conv)
        return(out)
      }
      
      u.iter <- u.iter+1
      
    }
    
    iter <- iter+1
    
    if (iter==maxiter){
      print("The maximum number of iterations was reached without convergence")
      conv <- "no"
      out <- list(conv=conv)
      return(out)
    }
    
    eta <- X%*%beta+Z%*%u
  }
  conv <- "yes"
  
  #Random components names
  rownames(u) <- seq(1:length(u))
  
  #Variances
  vcov.fixed <- Hb
  var.random <- diag(Hu)
  names(var.random) <- seq(1:length(u))
  
  
  #return
  out <- list(fixed.est=beta,random.est=u,vcov.fixed=vcov.fixed,var.random=var.random,iter.fixrand=iter,conv.fixrand=conv)
  
}