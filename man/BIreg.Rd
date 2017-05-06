\name{BIreg}
\alias{BIreg}
\title{
Fit a binomial logistic regression model
}
\description{
\code{BIreg} function fits a binomial logistic regression model, i.e., it links the probability parameter of a binomial distribution with the given covariates by means of a logistic link function. There is the option to include a dispersion parameter in the binomial distribution, which will be estimated by the bias corrected method of moments.
}
\usage{
  BIreg(formula,m,data,disp=FALSE,maxiter=20)
}
\arguments{
\item{formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
\item{m}{ number of trials in each binomial observation.}
\item{data}{ an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
\item{disp}{ if TRUE a dispersion parameter will be estimated. Default FALSE.}
\item{maxiter}{ the maximum number of iterations in IWLS method. Default 20.}
}


\value{
  \code{BIreg} returns an object of class "\code{BIreg}".
  
  The function \code{summary} (i.e., \code{summary.BIreg}) can be used to obtain or print a summary of the results.
  
  \item{coefficients}{ the estimated value of the regression coefficients.}
  \item{vcov}{ the variance and covariance matrix of the estimated regression coefficients.}
  \item{phi}{ if \code{disp} TRUE, it returns the estimated value of the dispersion parameter. If \code{disp} FALSE, then the estimated value is 1. Default FALSE.}
  \item{fitted.values}{ the fitted mean values of the model.}
  \item{residuals}{ working residuals, i.e. the residuals in the final iteration of the IWLS method.}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{null.deviance}{ null-deviance, deviance for the null model. The null model will include only an intercept.}
  \item{df.null}{ degrees of freedom for the null model.}
  \item{iter}{ number of iterations in the IWLS method.}
  \item{conv}{ if the algorithm has converged it returns "yes", otherwise "no".}
  \item{X}{ model matrix.}
  \item{y}{ dependent variable in the model.}
  \item{balanced}{ if the response binomial variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ number of trials in each binomial observation.}
  \item{nObs}{ number of observations.}
  \item{call}{ the matched call.}
  \item{formula}{ the formula supplied.}
}
\details{
\code{BIreg} function performs a regression model linking by a logistic function the probability paramater of a binomial distribution with a linear predictor that consists of the given covariates. Following the exponential family theory, the binomial distribution with dispersion parameter has the following log-likelihood function:
\deqn{l=[y*log(p/(1-p))+m*log(1-p)]/phi+c(y,phi)}
where \eqn{c()} is a known function. If we any dispersion parameter is not considered the usual density function of the binomial distribution will be used,
\deqn{l=y*log(p)+(m-y)*log((1-p)).}

As explained before we link the probablity parameter with the given covariates by
\deqn{logit(p)=log(p/(1-p))=x_i'*beta}
where \eqn{beta} are the regression coefficients and \eqn{x_i} is the \eqn{i}th row of a full rank design matrix \eqn{X} composed by the given covariables.

The estimation of the regression parameters \eqn{beta} is done via maximum likelihood approach, where the iterative weighted least square (IWLS) method is applied.

If \code{disp} is TRUE, a dispersion parameter will be added in the binomial distribution and, consequently, the method will deal with the general definition of the log-likelihood formula, otherwise the usual and simpler one will be used. In case the dispersion parameter is included, the estimation will be done with a bias-corrected method of moments:
\deqn{phi=Var[y]/[(m-q)*p*(1-p)]}
where \eqn{q} is the number of estimated regression paramters, and \eqn{p} is the estimated probability parameter.

The deviance of the model is defined by the ratio between the log-likelihood of the estimated model and saturated or null model. If the dispersion paramter is included the scaled deviance is obtained dividing the deviance by the dispersion parameter. 

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

\keyword{ HRQoL }
\seealso{
Iterative weighted least squares method function \code{\link{BIiwls}} in R-package \code{HRQoL}.
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
  model
  
  # with dispersion parameter
  phi <- 2
  outcome.disp <- rBI(k,m,p,phi)
  model.disp <- BIreg(outcome.disp~covariate,m,disp=TRUE)
  model.disp
}