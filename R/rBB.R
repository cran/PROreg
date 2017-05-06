rBB <- function(k,m,p,phi){
  
  if (k==as.integer(k)){
  } else {
    stop("k must be integer")
  }
  
  if (k<=0){
    stop("k must be positive")
  }
  
  if (sum(as.integer(m)==m)==length(m)){
  } else{
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
  
  if (min(p)<0 | max(p)>1){
    stop("p must be bounded between 0 and 1")
  }
  
  if (phi<=0){
    stop("phi must be positive")
  }
  
  alpha <-p/phi
  beta <- (1-p)/phi
  
  p. <- rbeta(k,alpha,beta)
  
  out <- rbinom(k,m,p.)
  return(out)
}