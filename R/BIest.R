BIest <- function(y,m,disp=FALSE,method="MM"){
  
  if (sum(as.integer(y)==y)==length(y)){
  } else {
    stop("y must be integer")
  }
  
  if (sum(m==as.integer(m))==length(m)){
  } else {
    stop("m must be integer")
  }
  
  if ((length(m)==1) | (length(m)==length(y))){
  } else{
    stop("m must be a number, or a vector of the length of y")
  }
  
  if (max(y-m)>0 | min(y) < 0){
    stop("y must be bounded between 0 and m")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  if ((disp==FALSE) | (disp==TRUE)){
  } else {
    stop("disp is logical")
  }
  
  if ((method=="MM") | (method=="MLE")){
  } else{
    stop("method option must be MM or MLE")
  }
  
  #Number of observations
  n <- length(y)
  #Balanced data
  if (length(m)==1){
    m <- rep(m,n)
    balanced <- "yes"
  } else{
    if(sum(m[1]==m)==length(m)){
      balanced <- "yes"
    } else {
      balanced <- "no"
    }
  }
  
  # The estimating equation for p is the same with both methodologies!
  p <- sum(y)/sum(m)
  # Ones we have computed the mle, we replace it in the Fisher information formula
  I <- sum(m)/(p*(1-p))
  
  if (disp){ 
    if (n==1){
      stop("Warning: Cannot calculate overdipersion with one observation")
    }
    a <- y==y[1]
    if (sum(a)==n){
      stop("Warning: Cannot calculate overdispersion if all the observations are equal")
    }
    
    if (method=="MLE"){
        #Estimate the dispersion parameter
        y0. <- which(y==0)
        if (length(y0.)==0){
          phi0 <- 0
        } else{
          y0 <- y[y0.]
          m0 <- m[y0.]
          phi0 <- sum(m0*log(1/(1-p)))
        }
        
        ym. <- which(y==m)
        if (length(ym.)==0){
          phim <- 0
        } else {
          ym <- y[ym.]
          mm <- m[ym.]
          phim <- sum(mm*log(1/p))
        }
        
        y0m. <- which((y==0) | (y==m))
        if (length(y0m.)==0){
          y0m <- y
          m0m <- m
          phi0m <- sum(y0m*log((y0m/m0m)/p)+(m0m-y0m)*log((1-y0m/m0m)/(1-p)))
        } else {
          y0m <- y[-y0m.]
          m0m <- m[-y0m.]
          phi0m <- sum(y0m*log((y0m/m0m)/p)+(m0m-y0m)*log((1-y0m/m0m)/(1-p)))
        }
        
        phi <- (2/(n-1))*(phi0+phi0m+phim)
    } else {
        phi <- sum((y-m*p)^2/(m*p*(1-p)))/(n-1)
      }
  } else {
    phi <- 1
  }
  
  #Stand. Errors of p
  pVar <- phi/I
  pSE <- sqrt(pVar)
  pIC <- c(p-1.96*pSE,p+1.96*pSE)

  
  out <- list(p=p, pVar=pVar, pIC.low=pIC[1],pIC.up=pIC[2],
              phi=phi,
              m=m,balanced=balanced,
              method=method)
  
  class(out) <- "BIest"
  
  out$call <- match.call()
  
  out
}

print.BIest <- function(x,...){
  
  cat("The probability parameter\n")
  cat("Estimation:",x$p,"\n")
  cat("Standard deviation:",sqrt(x$pVar),"\n")
  cat("95% confidence interval:", c(x$pIC.low,x$pIC.up),"\n")
  if (x$phi!=1){
    cat("\nThe dispersion parameter\n")
    cat("Estimation:",x$phi,"\n")
  }
  if (x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:",x$m[1])
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nEstimation approach:",x$method)
}
