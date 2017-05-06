\name{BBreg}
\alias{BBreg}

\title{
Fit a beta-binomial logistic regression model
}

\description{
\code{BBreg} function fits a beta-binomial logistic regression model, i.e., it links the probability parameter of a beta-binomial distribution with the given covariates by means of a logistic link function. The estimation of the parameters in the model is done via maximum likelihood estimation.
}

\usage{
  BBreg(formula,m,data,maxiter=100)
}

\arguments{
  \item{formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}
  \item{m}{ maximum score number in each beta-binomial observation.}
  \item{data}{ an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
  \item{maxiter}{ the maximum number of iterations in the estimation process. Default 100.}
}

\value{
  \code{BBreg} returns an object of class "\code{BBreg}".
  
  The function \code{summary} (i.e., \code{summary.BBreg}) can be used to obtain or print a summary of the results.
  
  \item{coefficients}{ the estimated value of the regression coefficients.}
  \item{vcov}{ the variance-covariance matrix of the estimated coefficients of the regression.}
  \item{phi}{ the estimation of the dispersion parameter of the beta-binomial distribution.}
  \item{psi}{ the estimation of the logarithm of the dispersion parameter of the beta-binomial distribution. }
  \item{psi.var}{ the variance of the estimated logarithm of the dispersion parameter of the beta-binomial distribution.}
  \item{conv}{ convergence of the methodology. If the method has converged it returns "yes", otherwise "no".}
  \item{fitted.values}{ the fitted mean values of the model.}
  \item{deviance}{ the deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{null.deviance}{ null-deviance, the deviance for the null model. The null model will only include an intercept as the estimation of the probability parameter.}
  \item{null.df}{ the degrees of freedom for the null model.}
  \item{iter}{ number of iterations in the estimation process.}
  \item{X}{ the model matrix.}
  \item{y}{ the dependent response variable in the model.}
  \item{m}{ maximum score number in each beta-binomial observation.}
  \item{balanced}{ if the response beta-binomial variable is balanced it returns "yes", otherwise "no".}
  \item{nObs}{ number of observations.}
  \item{call}{ the matched call.}
  \item{formula}{ the formula supplied.}
}
\details{
There are two different ways of defining a regression model based on the beta-binomial distribution: (i) the marginal regression approach, (ii) hierarchical generalized linear model approach. \cite{Najera-Zuloaga et al. (2017)} proved that the first approach is more adequate when the interest consists of the interpretation of the regression coefficients. Consequenlty, this function is based on the first approach, i.e., the marginal regression approach.

Once the marginal density function of the beta-binomial distribution is explicity calculated, we connect the probability parameter with the given covariates by means of a logistic link function:
\deqn{logit(p)=log(p/(1-p))=X*beta}
where \eqn{X} a model matrix composed by the given covariates and \eqn{beta} are the regression coefficients of the model.

Remplacing the previous linear predictor in the marginal density function of the beta-binomial distribution, we can derive maximum likelihood estimations of both regression and dispersion parameters. Forcina and Franconi (1988) presented an estimation algorithm based on the Newton-Raphson approach. This function performs the estimation of the parameters following the presented methodology. 
}

\author{
  J. Najera-Zuloaga
  
  D.-J. Lee
  
  I. Arostegui
}
\references{
Forcina A. & Franconi L. (1988): Regression analysis with Beta-Binomial distribution, \emph{Revista di Statistica Applicata}, \bold{21}, 7-12

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2017): Comparison of beta-binomial regression model approaches to analyze health related quality of life data, \emph{Statistical Methods in Medical Research}, DOI: 10.1177/0962280217690413 
}

\examples{
# We simulate a covariate, fix the paramters of the beta-binomial 
# distribution and simulate a response variable.

# Then we apply the model, and try to get the same values.
set.seed(18)
k <- 1000
m <- 10
x <- rnorm(k,5,3)

beta <- c(-10,2)
p <- 1/(1+exp(-(beta[1]+beta[2]*x)))
phi <- 1.2

y <- rBB(k,m,p,phi)

model <- BBreg(y~x,m)
model
}
