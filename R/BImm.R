
BImm <- function(fixed.formula,random.formula=NULL,Z=NULL,nRandComp=NULL,m,data=list(),maxiter=100){
  
  if ((is.null(random.formula))&(is.null(Z))){
    stop("Random part of the model must be especified")
  }
  
  if ((is.null(random.formula)==FALSE)&(is.null(Z))==FALSE){
    stop("Random part of the model has been specified twice")
  }
  
  if ((is.null(Z)==FALSE)&(is.null(nRandComp))){
    stop("Number of random components must be specified")
  }
  
  if ((is.null(Z))&(is.null(nRandComp)==FALSE)){
    stop("Number of random components must be specified only when Z is defined")
  }
  
  if (is.null(Z)==FALSE){
    if (dim(Z)[2]!=sum(nRandComp)){
      stop("The number of random effects in each random component must match with the number of columns of the design matrix Z")
    }
  }
  
  if (maxiter!=as.integer(maxiter)){
    stop("maxiter must be integer")
  }
  
  if(maxiter<=0){
    stop("maxiter must be positive")
  }
  
  if (sum(as.integer(m)==m)==length(m)){
  } else{
    stop("m must be integer")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  #Fixed effects design matrix:
  fixed.mf <- model.frame(formula=fixed.formula, data=data)
  X <- model.matrix(attr(fixed.mf, "terms"), data=fixed.mf)
  q <- dim(X)[2]
  
  #Outcome variable:
  y <- model.response(fixed.mf)
  nObs <- length(y)
  #Balanced data
  if(length(m)==1){
    balanced <- "yes"
    m. <- rep(m,nObs)
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
  
  #Specifiying random effects structure
  if (is.null(random.formula)){
    nComp <- length(nRandComp) #Number of random components, number of u,v...
    nRand <- dim(Z)[2] #Number of random effects, number of u_i,v_i,..
    namesRand <- as.character(seq(1,nComp,1))
  } else {
    #Random effects design matrix:
    random.mf <- model.frame(formula=update(random.formula,~.-1), data=data)
    nComp <- dim(random.mf)[2] #Number of random components, number of u,v...
    nRandComp <- NULL #Number of random effects in each random components, number of u_i in u.
    Z <- NULL
    for(i in 1:nComp){
      z <- model.matrix(~random.mf[,i]-1)
      Z <- cbind(Z,z)
      nRandComp <- c(nRandComp,dim(z)[2])
    }
    nRand <- dim(Z)[2] #Number of random effects, number of u_i,v_i,..
    namesRand <- names(random.mf)
  }
  
  # Initial values
  iter=0
  while.iter <- 0
  beta <- BIreg(fixed.formula,m,data)$coefficients
  
  all.sigma <- rep(1,nRand)
  oldall.sigma <- rep(Inf,nRand)
  
  eta <- X%*%beta
  oldeta <- rep(Inf,nObs)
  
  while (max(abs(oldall.sigma-all.sigma))>0.001){
    
    oldall.sigma<- all.sigma
    oldeta <- rep(Inf,nObs)
    
    while (max(abs(eta-oldeta))>0.001){
      
      oldeta <- eta
      
      # Relevant to the working vector
      p <- 1/(1+exp(-(eta)))
      mu <- m.*p
      g. <- m./(mu*(m.-mu))
      yw <- eta + (g.*(y-mu))
      W <- diag(c((m.*p*(1-p))))
      
      
      #Variance-covariance matrix of random effects
      d <- d. <- NULL
      for (i in 1:nComp){
        d <- c(d,rep(all.sigma[i],nRandComp[i]))
        d. <- c(d.,rep(1/(all.sigma[i]),nRandComp[i]))
      }
      D <- diag(d)
      D. <- diag(d.)
      
      u <- solve(t(Z)%*%W%*%Z+D.)%*%t(Z)%*%W%*%(yw-X%*%beta)
      beta <- solve(t(X)%*%W%*%X)%*%t(X)%*%W%*%(yw-Z%*%u)
      
      eta <- X%*%beta+Z%*%u       
    }
    
    # Calculating the variance component sigmau
    t <- 0
    t1 <- NULL    
    t2 <- NULL
    sig. <- NULL
    for (l in 1:nComp){
      t1 <- t+1
      t2 <- t+nRandComp[l]
    
      D. <- (1/all.sigma[l])*diag(nRandComp[l])
      
      sig.. <- as.numeric((1/nRandComp[l])*(t(u[seq(t1,t2)])%*%u[seq(t1,t2)]+sum(diag(solve(t(Z[,seq(t1,t2)])%*%W%*%Z[,seq(t1,t2)]+D.)))))
      sig. <- cbind(sig.,sqrt(sig..))
      
      t <- t+nRandComp[l]
    }
    all.sigma <- sig.
    
    iter=iter+1
    cat("Iteration number:", iter,"\n")
    if (iter==maxiter){
      print("The maximum number of iterations was reached without convergence")
      conv <- "no"
      out <- list(conv=conv)
      return(out)
    }
  }
  
  conv <- "yes"
  # Returning the final values
  d <- NULL
  for (l in 1:nComp){
    d. <- all.sigma[l]*rep(1,nRandComp[l])
    d <- c(d,d.)
  }
  D <- diag(d)
  
  W. <- diag(1/diag(W))
  V <- W.+Z%*%D%*%t(Z)
  
  # Beta
  fixed.coef <- beta
  V. <- solve(V)
  fixed.vcov <- solve(t(X)%*%V.%*%X)
  
  # Random efffects u
  random.coef <- u
  sigma.coef <- all.sigma
  P <- V.- V.%*%X%*%fixed.vcov%*%t(X)%*%V.
  Fisher.randvar <- -(1/2)*4*(sigma.coef^2)*(sum(diag(P%*%Z%*%t(Z)%*%P%*%Z%*%t(Z))))
  sigma.var <- -1/Fisher.randvar
  
  # fitted values and residuals
  fitted <- 1/(1+exp(-eta))

  # Deviance
  # See Pawitan book
  #deviance <- 2*log(2*pi)+log(det(solve(W)))+t(yw-eta)%*%W%*%(yw-eta)+log(det(solve(D.)))+t(u)%*%D.%*%u+log(det(t(Z)%*%W%*%Z+D.))
  # Breslow's approach
  deviance <- sum(log(diag(V)))+log(det(t(X)%*%V.%*%X))+t(yw-X%*%beta)%*%V.%*%(yw-X%*%beta)
  df <- nObs-length(fixed.coef)-1
  
  # output 
  out <- list(fixed.coef=fixed.coef,fixed.vcov=fixed.vcov,
              random.coef=random.coef,sigma.coef=sigma.coef,sigma.var=sigma.var,
              fitted.values=fitted,conv=conv,
              deviance=deviance, df=df,conv=conv,
              nRand=nRand,nRandComp=nRandComp,namesRand=namesRand,
              iter=iter,nObs=nObs,
              y=y,X=X,Z=Z,
              balanced=balanced,m=m)
  
  
  class(out) <- "BImm"
  
  out$call <- match.call()
  out$formula <- formula
  
  out
}



print.BImm <- function(x,...){
  cat("Call: \t")
  print(x$call)
  cat("\nFixed effects estimation:\n")
  print(t(x$fixed.coef))
  cat("\n")
  cat("\nStandard deviation of normal random effects:\n")
  for (i in 1:length(x$namesRand)){
    cat(x$namesRand[i], x$sigma.coef[i])
    cat("\n")
  }
  cat("\nDeviance of the model:",x$deviance)
  cat("\nNumber of iterations:",x$iter)
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m,"\n")
  } else {
    cat("\nNo balanced data.\n")
  }
}


summary.BImm <- function(object,...){
  
  # Fixed
  se <- sqrt(diag(object$fixed.vcov))
  tval <- object$fixed.coef/se
  TAB <- cbind(object$fixed.coef,se,tval,2*pnorm(-abs(tval)))
  colnames(TAB) <- c("Estimate","StdErr","t.value","p.value")
  
  #Sigma
  sigma.table <- cbind(object$sigma.coef,sqrt(object$sigma.var))
  rownames(sigma.table) <- object$namesRand
  colnames(sigma.table) <- c("Estimate","StdErr")
  
  
  out <- list(call=object$call,fixed.coefficients=TAB,
              random.coef=object$random.coef,sigma.table=sigma.table,
              fitted.values=object$fitted.values,residuals=object$residuals,
              deviance=object$deviance, df=object$df,
              nRand=object$nRand,nComp=object$nComp,nRandComp=object$nRandComp,namesRand=object$namesRand,
              iter=object$iter,nObs=object$nObs,
              y=object$y,X=object$X,Z=object$Z,
              balanced=object$balanced,m=object$m,
              conv=object$conv)
  
  class(out) <- "summary.BImm"
  out
}

print.summary.BImm <- function(x,...){
  cat("Call:\t")
  print(x$call)
  
  cat("\nFixed effects coefficients:\n")
  cat("\n")
  printCoefmat(x$fixed.coefficients,P.values=TRUE,has.Pvalue=TRUE)
  cat("\n---------------------------------------------------------------\n")
  cat("Random effects dispersion parameter(s):\n")
  cat("\n")
  print(x$sigma.table)
  cat("\n---------------------------------------------------------------\n")
  
  cat("Deviance of the model:",x$deviance,"; with", x$df,"degrees of freedom.")
  cat("\nNumber of observations:",x$nObs)
  cat("\nNumber of iterations:", x$iter)
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m)
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nNumber of random effects in each random component:",x$nRandComp,"\n")
  cat("\n")
}
