BBest <- function(y,m,method="MM"){
  
  if (sum(as.integer(m)==m)==length(m)){
  } else{
    stop("m must be integer")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  if (is.numeric(y)){
  } else{
    stop("y must be numeric")
  }
  
  if (sum(y%%1==0)==length(y)){
  } else {
    stop("y must be integer")
  }
  
  if ((length(m)==1) | (length(m)==length(y))){
  } else{
    stop("m must be a number, or a vector of the length of y")
  }
  
  if (min(m-y)<0 | min(y) < 0){
    stop("y must be bounded between 0 and m")
  }
  
  if (is.character(method)){
  } else {
    stop("method must be MM or MLE.")
  }
  
  if ((method=="MM") | (method=="MLE")){
  } else {
    stop("method must be MM or MLE.")
  }
  
  #Number of observations
  nObs <- length(y)
  #Balanced
  if (length(m)==1){
    m. <- rep(m,nObs)
    balanced <- "yes"
  } else{
    m. <- m
    if (sum(m[1]==m)==length(m)){
      balanced <- "yes"
    } else{
      balanced <- "no"
    }
  }
  
  #Method of moments
  if (method=="MM"){
    
    p <- sum(y)/sum(m.)
    psi.est <- function(psi){
      out <- (1/(nObs-1))*sum(((y-m.*p)^2)/(m.*p*(1-p)*(1+(m.-1)*(exp(psi)/(1+exp(psi))))))-1
      out
    }
    psi <- uniroot(psi.est,lower = -5,upper =3)$root
    phi <- exp(psi)
  } else{
    
    p <- 0.5
    phi <- 0.5
    psi <- log(phi)
    param <- c(p,phi)
    old.param <- c(Inf,Inf)
    iter <- 0
    
    while(max(abs(param-old.param))>0.001){
      
      old.param <- param
      
      partial.p <- function(p){
        out <- 0
        
        for (j in 1:nObs){
          out1 <- 0
          out2 <- 0
          
          if (y[j]==0){
            out1 <- 0
          }else{
            for (k in 0:(y[j]-1)){
              out1 <- out1+1/(p+k*phi)
            }
          }
          if (y[j]==m.[j]){
            out2 <- 0
          }else{
            for (k in 0:(m.[j]-y[j]-1)){
              out2 <- out2+1/(1-p+k*phi)
            }
          }
          out <- out + (out1-out2)
        }
        return(out)
      }
      
      p <- uniroot(partial.p,lower = 0.001,upper = 0.999)$root
      
      
      partial.psi <- function(psi){
        out <- 0
        for (j in 1:nObs){
          out1 <- 0
          out2 <- 0
          out3 <- 0
          
          if (y[j]==0){
            out1 <- 0
          }else{
            for (k in 0:(y[j]-1)){
              out1 <- out1+(k*exp(psi))/(p+k*exp(psi))
            }
          }
          if (y[j]==m.[j]){
            out2 <- 0
          }else{
            for (k in 0:(m.[j]-y[j]-1)){
              out2 <- out2+(k*exp(psi))/(1-p+k*exp(psi))
            }
          }
          for (k in 0:(m.[j]-1)){
            out3 <- out3 + (k*exp(psi))/(1+k*exp(psi))
          }
          out <- out + (out1+out2-out3)
        }
        return(out)
      }
      
      psi <- uniroot(partial.psi,lower = -5,upper = 3)$root
      
      phi <- exp(psi)
      param <- c(p,phi)
      
      iter <- iter+1
    }
    p <- param[1]
    phi <- param[2]
  }
  
  
  # The variance is estimated through the Fisher Information Matrix:
  partial.p <- function(p){
    out <- 0
    
    for (j in 1:nObs){
      out1 <- 0
      out2 <- 0
      
      if (y[j]==0){
        out1 <- 0
      }else{
        for (k in 0:(y[j]-1)){
          out1 <- out1+1/(p+k*phi)
        }
      }
      if (y[j]==m.[j]){
        out2 <- 0
      }else{
        for (k in 0:(m.[j]-y[j]-1)){
          out2 <- out2+1/(1-p+k*phi)
        }
      }
      out <- out + (out1-out2)
    }
    return(out)
  }
  
  partial.psi <- function(psi){
    out <- 0
    for (j in 1:nObs){
      out1 <- 0
      out2 <- 0
      out3 <- 0
      
      if (y[j]==0){
        out1 <- 0
      }else{
        for (k in 0:(y[j]-1)){
          out1 <- out1+(k*exp(psi))/(p+k*exp(psi))
        }
      }
      if (y[j]==m.[j]){
        out2 <- 0
      }else{
        for (k in 0:(m.[j]-y[j]-1)){
          out2 <- out2+(k*exp(psi))/(1-p+k*exp(psi))
        }
      }
      for (k in 0:(m.[j]-1)){
        out3 <- out3 + (k*exp(psi))/(1+k*exp(psi))
      }
      out <- out + (out1+out2-out3)
    }
    return(out)
  }
  
  var.p <- -1/grad(partial.p,p)
  var.psi <- -1/grad(partial.psi,psi)
  
  out <- list(p=p, phi=phi,
              pVar=var.p,psi=psi,psiVar=var.psi,
              m=m,balanced=balanced,method=method)
  
  class(out) <- "BBest"
  
  out$call <- match.call()
  
  out
}


print.BBest <- function(x,...){
  
  cat("The probability parameter of the beta-binomial distribution:",x$p,"\n")
  cat("The dispersion parameter of the beta-binomial distribution:",x$phi,"\n")
  if (length(x$m)==1){
    cat("\nBalanced data, maximum score number:",x$m,"\n")
  } else {
    cat("\nNo balanced data.\n")
  }
}


summary.BBest <- function(object,...){
  p.se <- sqrt(object$pVar)
  p.coef <- cbind(object$p,p.se)
  colnames(p.coef) <- c("Estimate","StdErr")
  
  psi.se <- sqrt(object$psiVar)
  psi.coef <- cbind(object$psi,psi.se)
  colnames(psi.coef) <- c("Estimate","StdErr")
  
  
  res <- list(call=object$call,coefficients=cbind(p=object$p,phi=object$phi),
              p.coefficients=p.coef, psi.coefficients=psi.coef,
              m=object$m)
  
  class(res) <- "summary.BBest"
  res
}


print.summary.BBest <- function(x,...){
  cat("Call:\t")
  print(x$call)
  cat("\n")
  cat("Probability parameter estimation:\n")
  print(x$p.coefficients)
  cat("\n")
  cat("Log dispersion parameter estimation:\n")
  print(x$psi.coefficients)
  
  if (length(x$m)==1){
    cat("\nBalanced data, maximum score number:",x$m,"\n")
  } else {
    cat("\nNo balanced data.\n")
  }
}
