EffectsEst.NR <- function(y,m,beta,u,phi,D.,X,Z,nDim){

  # Number of observations
  nObs <- length(y)/nDim

  # Number of fixed effects
  q <- length(beta)

  # Number of random effects
  nRand <- length(u)

  # Initial values
  oldbeta <- beta
  oldu <- u

  # Defining the score equation functions
  EffectsEst.function <- function(x){
    beta <- x[1:q]
    u <- x[(q+1):length(x)]
    p <- 1/(1+exp(-(X%*%beta+Z%*%u)))

    p[which(p<10^(-6))] <- 10^(-6)
    p[which(p>(1-10^(-6)))] <- 1-10^(-6)

    # S <- diag(c(p*(1-p)))
    s <- c(p*(1-p))

    t <- NULL
    for(i in 1:nDim){
      for (j in ((i-1)*nObs+1):(i*nObs)){
        t1 <- 0
        t2 <- 0
        if (y[j]==0){
          t1 <- 0
        }else{
          for (k in 0:(y[j]-1)){
            t1 <- t1+1/(p[j]+k*phi[i])
          }
        }
        if (y[j]==m[j]){
          t2 <- 0
        }else{
          for (k in 0:(m[j]-y[j]-1)){
            t2 <- t2+1/(1-p[j]+k*phi[i])
          }
        }
        t <- c(t,t1-t2)
      }
    }

    # FASTEN CALCULATIONS
    # c(F1=t(X)%*%S%*%t,F2=t(Z)%*%S%*%t-D.%*%u)
    c(F1=crossprod(X,s*t),F2=crossprod(Z,s*t)-D.%*%u)
  }

  # Geting the mle
  start <- c(c(oldbeta),c(oldu))
  mle <- multiroot(f=EffectsEst.function,start=start,ctol=0.01,atol=0.01,rtol=0.01,useFortran = TRUE)

  # Calculing the values
  beta <- mle$root[1:q]
  names(beta) <- colnames(X)
  u <- mle$root[(q+1):(nRand+q)]

  # Output
  out <- list(fixed.est=beta,random.est=u)

}
