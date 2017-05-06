\name{summary.BIreg}
\alias{summary.BIreg}

\title{
Summarizes a BIreg class model.
}

\description{
\code{summary.BIreg} si the BIreg specific method for the generic function \code{summary} which summarizes objects returned by modelling functions.
}

\usage{
\method{summary}{BIreg}(object, ...)
}

\arguments{
\item{object}{ a BIreg class model.}
\item{...}{ for extra arguments.}
}

\details{
\code{summary.BIreg} summarizes all the relevant information about the estimation of the parameters in a BIreg class model.

The function performs statistical significance hypothesis about the estimated regression parameters based on the normal distribution of the estimates. It also performs a goodness of fit test based on the difference between the calculated deviance of the model and the null deviance or deviance of the null model, which it is suppose to follow a Chi-square distribution with degrees of freedom equal to the difference in degrees of freedom of the models.
}

\value{
  \code{summary.BIreg} returns an object of class "\code{summary.BIreg}".
  
  \item{coefficients}{ a table with all the relevant information about the significance of the regression coefficients of the model. It includes the estimations, the standard errors of the estimations, the test-statistics and the p-values.}
  \item{phi}{ the estimated value of the dispersion parameter. If \code{disp} FALSE, then the estimated value is 1.}
  \item{deviance}{ the deviance of the model.}
  \item{df}{ the degrees of freedom of the model.}
  \item{null.deviance}{ the deviance for the null model. The null model will include only an intercept.}
  \item{df.null}{ the degrees of freedom for the null model.}
  \item{Goodness.of.fit}{ p-value of the goodness of fit test.}
  \item{iter}{ number of iterations in the IWLS method.}
  \item{conv}{ convergence of the methodology. If the algorithm has converged it returns "yes", otherwise "no".}
  \item{X}{ the model matrix.}
  \item{y}{ the dependent variable in the model.}
  \item{balanced}{ if the response variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ the number of trials in each observation.}
  \item{nObs}{ number of observations.}
  \item{balanced}{ if the response binomial variable is balanced it returns "yes", otherwise "no".}
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
\code{\link{BIreg}}, \code{\link{BIiwls}}
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
  outcome <- rBI(k,m,p)
  
  model <- BIreg(outcome~covariate,m,disp=FALSE)
  sum.model <- summary(model)
}