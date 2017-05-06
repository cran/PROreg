BBreg <- function(formula,m,data=list(),maxiter=100){
  
  if (sum(as.integer(m)==m)==length(m)){
  } else{
    stop("m must be integer")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  if (maxiter!=as.integer(maxiter)){
    stop("maxiter must be integer")
  }
  
  if(maxiter<=0){
    stop("maxiter must be positive")
  }
  
  #Get the response, covariates and model matrices
  mf <- model.frame(formula=formula, data=data)
  X <- model.matrix(attr(mf, "terms"), data=mf)
  y <- model.response(mf)
  #Number of observations
  n <- length(y)
  #Maximum number of scores and balanced data.
  if(length(m)==1){
    balanced <- "yes"
    m. <- rep(m,n)
  } else{
    m. <- m
    if (sum(m[1]==m)==length(m)){
      balanced <- "yes"
    } else{
      balanced <- "no"
    }
  }
  
  if (sum(as.integer(y)==y)==length(y)){
  } else {
    stop("y must be integer")
  }
  
  if ((length(m)==1) | (length(m)==length(y))){
  } else{
    stop("m must be a number, or a vector of the length of y")
  }
  
  if (max(y-m)>0 | min(y) < 0){
    stop("y must be bounded between 0 and m")
  }
  
  # We calculate which y are 0 or m.
  t0 <- which(y==0)
  if (length(t0!=0)){
    y0 <- y[-t0]
  }else{
    y0 <- y
  }
  tm <- which(y==m.)
  if (length(tm!=0)){
    ym <- y[-tm]
  }else{
    ym <- y
  }
  
  #BIG LOOP
  # Initial values
  iter <- 0
  phi <- 1
  psi <- log(phi)
  oldphi <- 10000
  beta <- c(1,rep(0,dim(X)[2]-1))
  
  while (abs(phi-oldphi)>0.001){
    
    oldphi <- phi
    oldbeta <- rep(Inf,dim(X)[2])
    
    #Beta LOOP
    beta.iter <- 0
    while (max(abs(beta-oldbeta))>0.001){
      oldbeta <- beta
      p <- 1/(1+exp(-(X%*%beta)))
      p1 <- which(p==1)
      p0 <- which(p==0)
      if (length(p1)!=0){
        p[p1] <- 0.99999
      }
      if(length(p0)!=0){
        p[p0] <- 0.00001
      }
      
      # Compute S
      s <- p*(1-p)
      
      # Compute u and v
      u <- NULL
      v <- NULL
      for (j in 1:n){
        u1 <- 0
        u2 <- 0
        v1 <- 0
        v2 <- 0
        if (y[j]==0){
          u1 <- 0
          v1 <- 0
        }else{
          for (k in 0:(y[j]-1)){
            u1 <- u1+1/(p[j]+k*phi)
            v1 <- v1+1/(p[j]+k*phi)^2
          }
        }
        if (y[j]==m.[j]){
          u2 <- 0
          v2 <- 0
        }else{
          for (k in 0:(m.[j]-y[j]-1)){
            u2 <- u2+1/(1-p[j]+k*phi)
            v2 <- v2+1/(1-p[j]+k*phi)^2
          }
        }
        u <- c(u,u1-u2)
        v <- c(v,v1+v2)
      }
      
      # Computing y
      sv <- s*v
      Y <- X%*%beta+(1/sv)*u
      
      # Updating beta
      XSV <- sweep(X,MARGIN=1,s*sqrt(v),'*')
      H <- solve(crossprod(XSV))
      beta <- H%*%t(X)%*%as.vector(s*v*s*Y)
      
      beta.iter <- beta.iter+1
      
      if (beta.iter>maxiter){
        print("The maximum number of iterations was has been reached, the method has not converged")
        conv <- "no"
        out <- list(conv=conv)
        return(out)
      }
    }
    
    conv <- "yes"
    
    # Once we reach convergence, we have to update phi using the profile likelihood
    if (length(t0)!=0){
      p0 <- p[-t0]
      m0 <- m.[-t0]
    }else{
      p0 <- p
      m0 <- m.
    }
    if (length(tm)!=0){
      pm <- p[-tm]
      mm <- m.[-tm]
    }else{
      pm <- p
      mm <- m.
    }
    
    # We define the profile likelihood of log(phi)
    Lpsi <- function(psi){ 
      Lp1 <- 0
      Lp2 <- 0
      Lp3 <- 0
      
      for (j in 1:length(y0)){
        for (k in 0:(y0[j]-1)){
          Lp1 <- Lp1+k*exp(psi)/(p0[j]+k*exp(psi))
        }
      }
      for (j in 1:length(ym)){
        for (k in 0:(mm[j]-ym[j]-1)){
          Lp2 <- Lp2+k*exp(psi)/(1-pm[j]+k*exp(psi))
        }
      }
      for (j in 1:n){ 
        for (k in 0:(m.[j]-1)){
          Lp3 <- Lp3+k*exp(psi)/(1+k*exp(psi))
        }
      }
      out <- Lp1+Lp2-Lp3
      return(out)
    }
    # We find out the mle, psi, based on the profile likelihood
    psi.mle <- uniroot(Lpsi,lower=-20,upper=4,maxiter=maxiter)
    psi <- psi.mle$root
    phi <- exp(psi)
    
    #Number of iterations
    iter <- iter +1
    
    if (iter>maxiter){
      print("The maximum number of iterations was has been reached, the method has not converged")
      conv <- "no"
      out <- list(conv=conv)
      return(out)
    }
  }
  
  #VARIANCES ESTIMATION
  #beta
  vcov.b <- H
  coef.b <- beta
  if (length(beta)==1){
    rownames(coef.b) <- c("Intercept")
  } else{
    rownames(coef.b) <- c("Intercept",rownames(beta)[2:length(beta)])
  }
  #log(phi)=psi
  psi.var <- -1/grad(Lpsi,psi)
  
  
  #Fitted value and residuals
  fitted.values <- p

  #DEVIANCE
    dev1 <- dev1.null <- 0
    dev2 <- dev2.null <- 0
    dev3 <- dev3.null <- 0
    
    e <- sum(y)/sum(m.)
    
    for (j in 1:length(y0)){
      for (k in 0:(y0[j]-1)){
        dev1 <- dev1+log(p0[j]+k*phi)
        dev1.null <- dev1.null+log(e+k*phi)
      }
    }
    for (j in 1:length(ym)){
      for (k in 0:(mm[j]-ym[j]-1)){
        dev2 <- dev2+log(1-pm[j]+k*phi)
        dev2.null <- dev2.null+log(1-e+k*phi)
      }
    }
    for (j in 1:n){ 
      for (k in 0:(m.[j]-1)){
        dev3 <- dev3+log(1+k*phi)
        dev3.null <- dev3.null+log(1+k*phi)
      }
    }
    dev <- dev1+dev2-dev3
    dev.null <- dev1.null+dev2.null-dev3.null
  #The deviance
  deviance <- as.numeric(2*(-sum(lgamma(m.+1)-(lgamma(y+1)+lgamma(m.-y+1)))-dev))
  deviance.null <- as.numeric(2*(-sum(lgamma(m.+1)-(lgamma(y+1)+lgamma(m.-y+1)))-dev.null))
  
  #Degrees of freedom
  df <- n-length(beta)-1
  df.null <- n-1
  
  #The output
  out <- list(coefficients=coef.b,vcov=vcov.b,
              phi=phi,psi=psi,psi.var=psi.var,
              conv=conv,fitted.values=fitted.values,
              deviance=deviance,df=df,null.deviance=deviance.null,null.df=df.null,
              iter=iter,X=X,y=y,m=m,balanced=balanced,nObs=n)
  
  class(out) <- "BBreg"
  
  out$call <- match.call()
  out$formula <- formula
  
  out
  
}


print.BBreg <- function(x,...){
  cat("Call:\t")
  print(x$call)
  cat("\nBeta coefficients:\n")
  print(t(x$coefficients))
  cat("\nDispersion parameter:",x$phi,"\n")
  cat("\nDeviance:",x$deviance, " on ", x$df, " degrees of freedom\n")
  cat("Null deviance:",x$null.deviance,"on",x$null.df," degrees of freedom\n")
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m,"\n")
  } else {
    cat("\nNo balanced data.\n")
  }
}


summary.BBreg <- function(object,...){
  beta.se <- sqrt(diag(object$vcov))
  beta.tval <- object$coefficients/beta.se
  beta.TAB <- cbind(object$coefficients,beta.se,beta.tval,2*pt(-abs(beta.tval),df=object$df))
  colnames(beta.TAB) <- c("Estimate","StdErr","t.value","p.value")
  
  coefficients.psi=cbind(object$psi,sqrt(object$psi.var))
  colnames(coefficients.psi) <- c("Estimate","StdErr")
  rownames(coefficients.psi) <- c("log(phi)")
  
  Chi <- object$null.deviance-object$deviance
  Chi.p.value <- 1-pchisq(Chi,object$null.df-object$df)
  
  out <- list(call=object$call,coefficients=beta.TAB,
              psi.table=coefficients.psi,
              deviance=object$deviance,df=object$df,null.deviance=object$null.deviance,null.df=object$null.df,
              Goodness.of.fit=Chi.p.value,
              iter=object$iter,
              X=object$X,y=object$y,Z=object$Z,nObs=object$nObs,
              m=object$m,balanced=object$balanced,
              conv=object$conv)
  
  class(out) <- "summary.BBreg"
  out
}


print.summary.BBreg <- function(x,...){
  cat("Call:\t")
  print(x$call)
  cat("\n")
  cat("Beta coefficients:\n")
  cat("\n")
  printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("Dispersion parameter coefficients:\n")
  cat("\n")
  print(x$psi.table)
  cat("\n---------------------------------------------------------------\n")
  cat("Deviance:",x$deviance," on ", x$df, " degrees of freedom\n")
  cat("Null deviance:",x$null.deviance," on ", x$null.df, " degrees of freedom\n")
  cat("Deviance test p-value:",x$Goodness.of.fit,"\n")
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m)
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nNumber of iterations:", x$iter,"\n")
}