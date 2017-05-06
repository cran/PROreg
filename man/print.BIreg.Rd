\name{print.BIreg}
\alias{print.BIreg}
\title{
Print a BIreg class model.
}
\description{
\code{print.BIreg} is the BIreg specific method fot the generic function print which prints objects returned by modelling functions.
}

\usage{
\method{print}{BIreg}(x, ...)
}

\arguments{
\item{x}{ a BIreg class model.}
\item{...}{ for extra arguments.}
}

\value{
Prints a BIreg object.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Pawitan Y. (2001): In All Likelihood: Statistical Modelling and Inference Using Likelihood, \emph{Oxford University Press}

Williams D. A. (1982): Extra-Binomial Variation in Logistic Linear Regression, \emph{Journal of the Royal Statistical Society. Series C}, \bold{31}, 144-148
}

\seealso{
\code{\link{BIreg}}
}

\examples{
  set.seed(1234)
  # We simulate a covariable and construct the outcome variable applying
  # an inverse logit link function on it.
  
  m <- 10                       
  k <-100                       
  covariate <- rnorm(k,2,0.5)  
  
  beta <- c(-6,4)
  p <- 1/(1+exp(-(beta[1]+beta[2]*covariate)))
  
  # without dispersion parameter
  outcome <- rBI(k,m,p)
  model <- BIreg(outcome~covariate,m,disp=FALSE)
  print(model) # or just model
}
