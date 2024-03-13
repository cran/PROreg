EffectsEst.BBDelta <- function(y,m,beta,u,phi,D.,X,Z,maxiter=200,nDim){

  # Number of observations
  nObs <- length(y)/nDim

  # Number of fixed effects
  q <- length(beta)

  # Number of random effects
  nRand <- length(u)

  #Initial values
  iter <- 0
  eta <- X%*%beta+Z%*%u
  oldeta <- rep(Inf,nObs)

  tol <- 1
  while (tol>0.001){

    oldeta <- eta

    #----------------------------
    # Estimation of fixed effects
    #----------------------------
    oldbeta <- rep(Inf,q)
    beta.iter <- 1

    while (max(abs((beta-oldbeta)/beta))>0.01){

      oldbeta <- beta
      p <- 1/(1+exp(-(X%*%beta+Z%*%u)))
      p[which(p<10^(-6))] <- 10^(-6)
      p[which(p>(1-10^(-6)))] <- 1-10^(-6)

      # Compute s
      s <- c(p*(1-p))

      # Compute t and v
      t <- NULL
      v <- NULL
      for(i in 1:nDim){
        for (j in ((i-1)*nObs+1):(i*nObs)){
          t1 <- 0
          t2 <- 0
          v1 <- 0
          v2 <- 0
          if (y[j]==0){
            t1 <- 0
            v1 <- 0
          }else{
            for (k in 0:(y[j]-1)){
              t1 <- t1+1/(p[j]+k*phi[i])
              v1 <- v1+1/((p[j]+k*phi[i])^2)
            }
          }
          if (y[j]==m[j]){
            t2 <- 0
            v2 <- 0
          }else{
            for (k in 0:(m[j]-y[j]-1)){
              t2 <- t2+1/(1-p[j]+k*phi[i])
              v2 <- v2+1/((1-p[j]+k*phi[i])^2)
            }
          }
          t <- c(t,t1-t2)
          v <- c(v,v1+v2)
        }
      }

      #  Reducing to the working vector
      # etab <- X%*%beta+(1/s*v)*t
      # beta <- solve(t(X)%*%diag(s*v*s)%*%X)%*%t(X)%*%matrix(s*v*s*etau,ncol=1)

      # Without reducing to the working vector
      # beta <- Matrix::solve(t(X)%*%diag(s*v*s)%*%X)%*%t(X)%*%matrix(s*v*s*etab,ncol=1) # SLOW
      beta <- beta+crossprod(Matrix::solve(crossprod(X*s*sqrt(v))),crossprod(X,as.matrix(s*t,ncol=1))) # FAST

      if (beta.iter==500){
        cat("-Delta algorithm failed to converge fixed effects, NR is used instead-\n")
        conv <- "no"
        out <- list(conv=conv,beta=beta,u=u)
        return(out)
        }

      beta.iter <- beta.iter +1
    }

    oldu <- rep(Inf,nRand)
    u.iter <- 1
    while (max(abs(u-oldu))>0.01){
      # max(abs((u-oldu)/u)) sometimes there are problems because u can be 0 and then we cannot calculate the tolerance like that
      oldu <- u
      p <- 1/(1+exp(-(X%*%beta+Z%*%u)))
      p[which(p<10^(-6))] <- 10^(-6)
      p[which(p>(1-10^(-6)))] <- 1-10^(-6)

      # Compute S
      s <- c(p*(1-p))

      # Compute t and v
      t <- NULL
      v <- NULL
      for(i in 1:nDim){
        for (j in ((i-1)*nObs+1):(i*nObs)){
          t1 <- 0
          t2 <- 0
          v1 <- 0
          v2 <- 0
          if (y[j]==0){
            t1 <- 0
            v1 <- 0
          }else{
            for (k in 0:(y[j]-1)){
              t1 <- t1+1/(p[j]+k*phi[i])
              v1 <- v1+1/((p[j]+k*phi[i])^2)
            }
          }
          if (y[j]==m[j]){
            t2 <- 0
            v2 <- 0
          }else{
            for (k in 0:(m[j]-y[j]-1)){
              t2 <- t2+1/(1-p[j]+k*phi[i])
              v2 <- v2+1/((1-p[j]+k*phi[i])^2)
            }
          }
          t <- c(t,t1-t2)
          v <- c(v,v1+v2)
        }
      }


      #  Reducing to the working vector
      # etau <- Z%*%u+(1/s*v)*t
      # u <- solve(t(Z)%*%diag(s*v*s)%*%Z+D.)%*%t(Z)%*%matrix(s*v*s*etau,ncol=1)

      # Without reducing to the working vector
      # u <- u+solve(t(Z)%*%diag(s*v*s)%*%Z+D.)%*%(t(Z)%*%as.matrix(s*t,ncol=1)-D.%*%u) # SLOW
      u <- u+crossprod(Matrix::solve(crossprod(Z*s*sqrt(v))+D.),crossprod(Z,as.matrix(s*t,ncol=1))-crossprod(D.,u)) # FAST


      if (u.iter==500){
        cat("-Delta algorithm failed to converge random effects, NR is used instead-\n")
        conv <- "no"
        out <- list(conv=conv,beta=beta,u=u)
        return(out)
      }

      u.iter <- u.iter+1

    }

    # Calculate the tolerance
    eta <- X%*%beta+Z%*%u
    tol <- norm(eta-oldeta)/norm(eta)

    iter <- iter+1

    if (iter==maxiter){
      cat("-Delta algorithm failed to converge, NR is used instead-\n")
      conv <- "no"
      out <- list(conv=conv,beta=beta,u=u)
      return(out)
    }
  }
  conv <- "yes"

  #Random components names
  rownames(u) <- seq(1:length(u))

  #Variances
  Hb <- solve(crossprod(X*s*sqrt(v)))
  Hu <- solve(crossprod(Z*s*sqrt(v))+D.)
  vcov.fixed <- Hb
  var.random <- diag(Hu)
  names(var.random) <- seq(1:length(u))


  #return
  out <- list(fixed.est=beta,random.est=u,vcov.fixed=vcov.fixed,var.random=var.random,iter.fixrand=iter,conv.fixrand=conv)

}
