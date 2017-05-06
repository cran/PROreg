dBB <- function(m,p,phi){
  
  if ((length(m)>1) | (length(p)>1) | (length(phi)>1)){
    stop("parameters of the beta-binomial distribution cannot change, i.e., m, p and phi must be numbers")
  }
  
  if (m==as.integer(m)){
  } else {
    stop("m must be integer")
  }
  
  if (m<=0){
    stop("m must be positive")
  }

  if (p<0 | p>1){
    stop("p must be bounded between 0 and 1")
  }
  
  if (phi<0){
    stop("phi must be positive")
  }
  
  if (phi==0){
    stop("phi cannot be zero, use the binomial distribution instead")
  }
  
  tt <- seq(0,m)
  tt1 <- gamma(m + 1)/(gamma(tt + 1) * gamma(m - tt + 1))
  tt2 <- gamma((p/phi) + tt)/gamma(p/phi)
  tt3 <- gamma(((1 - p)/phi) + m - tt)/gamma((1 - p)/phi)
  tt4 <- gamma(1/phi)/gamma((1/phi) + m)
  out <- tt1 * tt2 * tt3 * tt4
  
  out  
}