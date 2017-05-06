BIreg <- function(formula,m,data=list(),disp=FALSE,maxiter=20){
  
  if (sum(m==as.integer(m))==length(m)){
  } else {
    stop("m must be integer")
  }
  
  if (min(m)<=0){
    stop("m must be positive")
  }
  
  if ((disp==FALSE) | (disp==TRUE)){
  } else {
    stop("disp is logical")
  }
  
  if (maxiter!=as.integer(maxiter)){
    stop("maxiter must be integer")
  }
  
  if(maxiter<=0){
    stop("maxiter must be positive")
  }
  
  #Decomposing the model matrices.
  mf <- model.frame(formula=formula,data=data)
  x <- model.matrix(attr(mf,"terms"),data=mf)
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
  
  # Calculing the coefficients
  test <- BIiwls(y,x,m.,maxiter)
  conv <- test$conv
  if (conv=="no"){
    stop("The iterative least squares method has not converged")
  }
  coef <- test$beta
  df <- nrow(x)-ncol(x) #degrees of freedom
  vcov <- test$vcov
  colnames(vcov) <- rownames(vcov) <- colnames(x)
  
  # Compute the overdispersion
  h <- (1/(1+exp(-as.vector(x%*%coef))))
  i0 <- which(y==0)
  im <- which(y==m.)
  i0m <- which((y==0) | (y==m.))
  if (disp==FALSE){
    phi=1
  } else{
    phi <- sum((y-m.*h)^2/(m.*h*(1-h)))/df
  }
  
  #Variance of the regression parameters
  vcov <- phi*vcov
  
  # Calculing the fitted values, residuals and degrees of freedom
  fitted.values <- h
  residuals <- y-m.*h
  
  # Calculating deviance / escaled deviance
  e <- mean(y/m.) #Observed probability parameter
  #y=0
  if (length(i0)==0){
    Deviance0 <- 0
    nullDeviance0 <- 0
  } else{
    y0 <- y[i0]
    h0 <- h[i0]
    m0 <- m.[i0]
    Deviance0 <- sum((m0-y0)*log((1-y0/m0)/(1-h0)))
    nullDeviance0 <- sum((m0-y0)*log((1-y0/m0)/(1-e)))
  }
  #y=m
  if (length(im)==0){
    Deviancem <- 0
    nullDeviancem <- 0
  } else {
    ym <- y[im]
    hm <- h[im]
    mm <- m.[im]
    Deviancem <- sum(ym*log((ym/mm)/hm))
    nullDeviancem <- sum(ym*log((ym/mm)/e))
  }
  #y!=0,y!=m
  if (length(i0m)==0){
    y0m <- y
    h0m <- h
    m0m <- m.
    Deviance0m <- sum(y0m*log((y0m/m0m)/h0m)+(m0m-y0m)*log((1-y0m/m0m)/(1-h0m)))
    nullDeviance0m <- sum(y0m*log((y0m/m0m)/e)+(m0m-y0m)*log((1-y0m/m0m)/(1-e)))
  } else {
    y0m <- y[-i0m]
    h0m <- h[-i0m]
    m0m <- m.[-i0m]
    Deviance0m <- sum(y0m*log((y0m/m0m)/h0m)+(m0m-y0m)*log((1-y0m/m0m)/(1-h0m)))
    nullDeviance0m <- sum(y0m*log((y0m/m0m)/e)+(m0m-y0m)*log((1-y0m/m0m)/(1-e)))
  }
  #The deviance
  deviance <- 2*(Deviance0+Deviancem+Deviance0m)
  null.df <- n-1
  null.deviance <- 2*(nullDeviance0+nullDeviancem+nullDeviance0m)
  
  # Number of iterations in IWLS
  iter <- test$iter
  
  #The output
  out <- list(coefficients=coef,vcov=vcov,phi=phi,
              fitted.values=fitted.values,residuas=residuals,
              deviance=deviance,df=df,null.deviance=null.deviance,null.df=null.df,
              iter=iter,conv=conv,X=x,y=y,balanced=balanced,m=m,nObs=n)
  
  class(out) <- "BIreg"
  
  out$call <- match.call()
  out$formula <- formula
  
  out
  
}

print.BIreg <- function(x,...){
  cat("Call:")
  print(x$call)
  cat("\nCoefficients:\n")
  print(t(x$coefficients))
  cat("\nDispersion parameter:",x$phi,"\n")
  cat("\nDegrees of freedom:", x$null.df, "Total ;", x$df, "Model")
  cat("\nNull deviance:",x$null.deviance)
  cat("\nModel Deviance:",x$deviance,"\n")
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m,"\n")
  } else {
    cat("\nNo balanced data.\n")
  }
}


summary.BIreg <- function(object,...){
  se <- sqrt(diag(object$vcov))
  tval <- object$coefficients/se
  TAB <- cbind(object$coefficients,se,tval,2*pt(-abs(tval),df=object$df))
  colnames(TAB) <- c("Estimate","StdErr","t.value","p.value")
  
  null.df <- object$null.df
  df <- object$df
  
  Chi <- (object$null.deviance-object$deviance)/object$phi
  Chi.p.value <- 1-pchisq(Chi,null.df-df)
  
  out <- list(call=object$call,coefficients=TAB,phi=object$phi,
              deviance=object$deviance,df=df,null.deviance=object$null.deviance,null.df=null.df,
              Goodness.of.fit=Chi.p.value,iter=object$iter,conv=object$conv,
              X=object$X,y=object$y,
              m=object$m, nObs=object$nObs,
              balanced=object$balanced)
  
  class(out) <- "summary.BIreg"
  out
}


print.summary.BIreg <- function(x,...){
  cat("Call:\n")
  print(x$call)
  cat("\n")
  
  printCoefmat(x$coefficients,P.values=TRUE,has.Pvalue=TRUE)  
  cat("\nDispersion parameter:",x$phi,"\n")
  cat("\nNull deviance:",x$null.deviance," on ",x$null.df, "degrees of freedom")
  cat("\nModel Deviance:",x$deviance," on ", x$df, "degrees of freedom") 
  cat("\nDeviance test:",x$Goodness.of.fit,"\n")
  
  if(x$balanced=="yes"){
    cat("\nBalanced data, maximum score in the trials:", x$m)
  } else {
    cat("\nNo balanced data.")
  }
  cat("\nNumber of iterations in IWLS:", x$iter,"\n")
  
}


# predict.BIreg <- function(object,newdata=NULL,...){
#   
#   if(is.null(newdata))
#     y <- fitted(object)
#   else{
#     if(!is.null(object$formula)){
#       x <- model.matrix(object$formula,newdata)
#     } else{
#       x <- newdata
#     }
#     y <- as.vector(x%*%coef(object))
#   }
#   y
# }