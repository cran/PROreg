rBI <- function(k,m,p,phi=1){
  
  if (k==as.integer(k)){
  } else {
    stop("k must be integer")
  }
  
  if (k>0){
  } else {
    stop("k must be positive")
  }
  
  if (sum(m==as.integer(m))==length(m)){
  } else {
    stop("m must be integer")
  }
  
  if ((length(m)>1) & (length(m)<k)){
    stop("m must a number or a vector of length k")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  if ((length(p)>1) & (length(p)<k)){
    stop("p must a number or a vector of length k")
  } 
  
  if (min(p)<0 | max(p) >1){
    stop("p must be bounded between 0 and 1")
  }
  
  if (phi < 0){
    stop("phi must be positive")
  }
  
  if (phi==1){
    out <- rbinom(k,m,p)
  } else{
    
    if (length(m)==1){
      m <- rep(m,k)
    }
    
    if (length(p)==1){
      p <- rep(p,k)
    }
    
    out <- NULL
    for (t in 1:k){
      #A little correction if p=0 or p=1.
      if (p[t]==0){
        p[t] <- 0.001
      }
      if(p[t]==1){
        p[t] <- 0.999
      }
      value. <- NULL
      
      #y=0
      y0 <- 0
      logout0 <- -(1/2)*log(2*pi*phi*m[t]*p[t]*(1-p[t]))-(2/(2*phi))*((m[t]-y0)*log((1-y0/m[t])/(1-p[t])))
      value0 <- exp(logout0)
      
      #y!=0 & y!=m
      value0m <- NULL
      if (m[t]>1){
        y0m <- 1:(m[t]-1)
        logout0m <- -(1/2)*log(2*pi*phi*m[t]*p[t]*(1-p[t]))-(2/(2*phi))*(y0m*log(y0m/(m[t]*p[t]))+(m[t]-y0m)*log((1-y0m/m[t])/(1-p[t])))
        value0m <- exp(logout0m)
      }
      value. <- c(value0,value0m)
      
      #y=m
      ym <- m[t]
      logoutm <- -(1/2)*log(2*pi*phi*m[t]*p[t]*(1-p[t]))-(2/(2*phi))*(ym*log(ym/(m[t]*p[t])))
      valuem <- exp(logoutm)
      
      #The probability of each y option.
      p. <-c(value.,valuem)/sum(value.,valuem)
      
      #We divide the [0,1] interval by those points.
      l <- NULL
      for (i in 1:(m[t]+1)){
        l[i] <- sum(p.[1:i])
      }
      
      u <- runif(1,0,1)
      for (i in (m[t]+1):1){
        if (u < l[i]){
          o <- i-1 
        }
      }
      out <- c(out,o)
      
    }
  }
  
  return(out)
}
