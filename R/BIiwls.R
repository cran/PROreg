BIiwls <- function(y,X,m,maxiter){
  
  # Initial values
  beta <- c(1,rep(0,dim(X)[2]-1))
  oldbeta <- rep(100,dim(X)[2])
  iter <- 0
  
  # The bucle
  while (max(abs(beta-oldbeta))>0.001){
    oldbeta <- beta
    p <- 1/(1+exp(-(X%*%beta)))
    
    #The method breaks when p=1 or p=0.
    p1 <- which(p==1)
    p0 <- which(p==0)
    if (length(p1)!=0){
      p[p1] <- 0.99999
    }
    if(length(p0)!=0){
      p[p0] <- 0.00001
    }
    
    # Working vector
    Y <- X%*%beta+(y-m*p)/(m*p*(1-p))
    # Weights
    W <- diag(c(m*p*(1-p)))
    
    U <- solve(t(X)%*%W%*%X)
    
    beta <-U%*%t(X)%*%as.vector(diag(W)*Y)
    iter <- iter+1
    if (iter > maxiter){
      conv <- "no"
      return(list(conv=conv))
    }
  }
  
  conv <- "yes"
  
  return(list(beta=beta,vcov=U,iter=iter,conv=conv))
}