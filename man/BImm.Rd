\name{BImm}
\alias{BImm}
\title{
Binomial Logistic Mixed-Effects Model Regression.}
\description{
\code{BImm} function performs binomial logistic mixed-effects models, i.e., it allows the inclusion of gaussian random effects in the linear predictor of a logistic binomial regression model.

The structure of the random part of the model can be expecified by two different ways: (i) determining the \code{random.formula} argument, or (ii) especifying the model matrix of the random effects, \code{Z}, and determining the number of random effects in each random component, \code{nRandComp}.
}

\usage{
BImm(fixed.formula,random.formula,Z=NULL,nRandComp=NULL,m,data,maxiter=100)
}

\arguments{
\item{fixed.formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the fixed part of the model to be fitted.}
\item{random.formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the random part of the model to be fitted. It must be specified in cases where the model matrix of the random effects \code{Z} is not determined.}
\item{Z}{ the design matrix of the random effects. If the \code{random.formula} argument is specified this argument should not be specified, as we will be specifying twice the random structure of the model.}
\item{nRandComp}{ the number of random components/levels in each random effect of the model. It must be especified as a vector, where the 'i'th value must correspond with the number of random components of the 'i'th random effect. It must be only determined when we specify the random structure of the model by the model matrix of the random effects, \code{Z}.}
\item{m}{ number of trials in each binomial observation.}
\item{data}{ an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
\item{maxiter}{ the maximum number of iterations in the parameters estimation algorithm. Default 100.}
}

\details{
The model that is performed by this function is a especial case of generalized linear mixed models (GLMMs), in which conditioned on some random components the response variable has a binomial distribution. As in the binomial (logistic) regression a logit link function is applied to the probability parameter of the conditioned binomial distribution, allowing the inclusion of random effects in the linear predictor,

\deqn{logit(p)=X*beta+Z*u,}
where \eqn{p} is the probability parameter, \eqn{X} a full rank matrix composed by the covariables, \eqn{beta} the fixed effects, \eqn{Z} the design matrix for the random effects ang \eqn{u} are the random effects. These random effects are independent and have a normal distribution with the same variance and mean 0.

The model estimates the fixed effects, predicts the random effects, and gets the estimation of the random effects variance parameters.

The estimation procedure is done by likelihood approximation, via iteartive weighted least squares method. The process is performed in two steps:  (i) fixed and random parameters are estimated for some given values of the random effects variance parameters, and (ii) random effects variance parameters are estimated for some given regression and random coefficients using a penalized profile likelihood. The estimation approach iterates between (i) and (ii) until convergence is obtained.
}

\value{
  \code{BImm} returns an object of class "\code{BImm}".
  
  The function \code{summary} (i.e., \code{summary.BImm}) can be used to obtain or print a summary of the results.

  \item{fixed.coef}{ estimated value of the fixed coefficients in the regression.}
  \item{fixed.vcov}{ variance and covariance matrix of the estimated fixed coefficients in the regression.}
  \item{random.coef}{ predicted random effects of the regression.}
  \item{sigma.coef}{ estimated value of the random effects variance parameters.}
  \item{sigma.var}{ variance of the estimated value of the random effects variance parameters.}
  \item{fitted.values}{ the fitted mean values of the probability parameter of the conditional beta-binomial distribution.}
  \item{conv}{ convergence of the methodology. If the method has converged it returns "yes", otherwise "no".}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{nRand}{ number of random effects.}
  \item{nComp}{ number of random components.}
  \item{nRandComp}{ number of random effects in each random component of the model.}
  \item{namesRand}{ names of the random components.}
  \item{iter}{ number of iterations in the estimation method.}
  \item{nObs}{ number of observations in the data.}
  \item{y}{ dependent response variable in the model.}
  \item{X}{ model matrix of the fixed effects.}
  \item{Z}{ model matrix of the random effects.}
  \item{balanced}{ if the conditional beta-binomial response variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ maximum score number in each binomial observation.}
  \item{call}{ the matched call.}
  \item{formula}{ the formula supplied.}
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25

McCulloch C. E. & Searle S. R. (2001): Generalized, Linear, and Mixed Models, \emph{Jhon Wiley & Sons}

Pawitan Y. (2001): In All Likelihood: Statistical Modelling and Inference Using Likelihood, \emph{Oxford University Press}

}

\keyword{matrixcalc, rootSolve}

\examples{
set.seed(5)
# Fixing parameters for the simulation:
nObs <- 1000
m <- 10
beta <- c(1.5,-1.1)
sigma <- 0.8

# Simulating the covariate:
x <- runif(nObs,-5,5)

# Simulating the random effects:
z <- as.factor(rBI(nObs,5,0.5,2))
u <- rnorm(6,0,sigma)

# Getting the linear predictor and probability parameter.
X <- model.matrix(~x)
Z <- model.matrix(~z-1)
eta <- beta[1]+beta[2]*x+crossprod(t(Z),u)
p <- 1/(1+exp(-eta))

# Simulating the response variable
y <- rBI(nObs,m,p)

# Apply the model
model <- BImm(fixed.formula = y~x,random.formula = ~z,m=m)
model
}