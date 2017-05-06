\name{summary.BImm}
\alias{summary.BImm}
\title{
Summarizes a BImm class model.
}

\description{
\code{summary.BImm} si the BImm specific method for the generic function \code{summary} which summarizes objects returned by modelling functions.
}

\usage{
\method{summary}{BImm}(object, ...)
}

\arguments{
\item{object}{ a BImm class model.}
\item{...}{ for extra arguments.}
}

\details{
\code{summary.BImm} summarizes all the relevant information about the estimation of the parameters in a BImm class model.

The function performs statistical significance hypothesis about the estimated fixed parameters based on the normal distribution of the estimates.
}

\value{
  \code{summary.BImm} returns an object of class "\code{summary.BImm}".
  
  \item{fixed.coefficients}{ a table with all the relevant information about the significance of the fixed effects of the model. It includes the estimations, the standard errors of the estimations, the test-statistics and the p-values.}
  \item{random.coef}{ predicted random effects of the regression.}
  \item{sigma.table}{ a table which inlcudes the estimation and standard errors of the parameters which the variance-covariance matrix of the random effects consists of.}
  \item{fitted.values}{ the fitted mean values of the probability parameter of the conditional beta-binomial distribution.}
  \item{residuals}{ residuals of the model.}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{nRand}{ number of random effects.}
  \item{nComp}{ number of random components.}
  \item{nRandComp}{ number of random components in each random effect of the model.}
  \item{namesRand}{ names of the random components.}
  \item{iter}{ number of iterations in the estimation method.}
  \item{nObs}{ number of observations in the data.}
  \item{y}{ dependent response variable in the model.}
  \item{X}{ model matrix of the fixed effects.}
  \item{Z}{ model matrix of the random effects.}
  \item{balanced}{ if the conditional binomial response variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ number of trials in each binomial observation.}
  \item{conv}{ convergence of the methodology. If the algorithm has converged it returns "yes", otherwise "no".}
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

\seealso{
\code{\link{BImm}}
}


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
sum.model <- summary(model)
}