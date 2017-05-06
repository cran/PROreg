dBI <- function(m,p,phi=1){
  
  if ((length(m)>1) | (length(p)>1) | (length(phi)>1)){
    stop("parameters of the binomial distribution cannot change, i.e., m, p and phi must be numbers")
  }
  
  if (m==as.integer(m)){
  } else {
    stop("m must be integer")
  }
  
  if (m<=0){
    stop("m must be positive")
  }
  
  if (p<0 | p >1){
    stop("p must be bounded between 0 and 1")
  }
  
  if (phi < 0){
    stop("phi must be positive")
  }
  
  if (phi==1){
    s <- dbinom(0:m,m,p)
  } else{
    
    #y=0
    y0 <- 0
    logout0 <- -(1/2)*log(2*pi*phi*m*p*(1-p))-(2/(2*phi))*((m-y0)*log((1-y0/m)/(1-p)))
    out0 <- exp(logout0)
    
    out0m <- NULL
    #y!=0 and y!=m
    if (m>1){
      y0m <- 1:(m-1)
      logout0m <- -(1/2)*log(2*pi*phi*m*p*(1-p))-(2/(2*phi))*(y0m*log(y0m/(m*p))+(m-y0m)*log((1-y0m/m)/(1-p)))
      out0m <- exp(logout0m)
    }
    out0. <- c(out0,out0m)
    
    #y=m
    ym <- m
    logoutm <- -(1/2)*log(2*pi*phi*m*p*(1-p))-(2/(2*phi))*(ym*log(ym/(m*p)))
    outm <- exp(logoutm)
    
    s <-c(out0.,outm)/sum(out0.,outm)
  }
  
  s
}
