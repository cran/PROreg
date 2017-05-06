\name{summary.BBmm}
\alias{summary.BBmm}

\title{
Summarizes a BBmm class model.}

\description{
\code{summary.BBmm} si the BBmm specific method for the generic function \code{summary} which summarizes objects returned by modelling functions.
}

\usage{
\method{summary}{BBmm}(object, ...)
}

\arguments{
\item{object}{ a BBmm class model.}
\item{...}{ for extra arguments.}
}

\details{
\code{summary.BBmm} summarizes all the relevant information about the estimation of the parameters in a BBmm class model.

The function performs statistical significance hypothesis about the estimated fixed parameters based on the normal distribution of the estimates. It also performs a goodness of fit test based on the difference between the calculated deviance of the model and the null deviance or deviance of the null model, which it is suppose to follow a Chi-square distribution with degrees of freedom equal to the difference in degrees of freedom of the models.
}

\value{
  \code{summary.BBmm} returns an object of class "\code{summary.BBmm}".
  
  \item{fixed.coefficients}{ a table with all the relevant information about the significance of the fixed effects estimates in the model. It includes the estimates, the standard errors of the estimates, the test-statistics and the p-values.}
  \item{sigma.table}{ a table which inlcudes the estimates and the standard errors of the estimates of the random effects variance parameters.}
  \item{psi.table}{ a table which includes the estimate and the standard errors of the estimate of the logarithm of the dispersion parameter of the conditional beta-binomial distribution.}
  \item{random.coef}{ predicted random effects of the regression.}
  \item{iter}{ number of iterations in the estimation method.}
  \item{nObs}{ number of observations in the data.}
  \item{nRand}{ number of random effects.}
  \item{nComp}{ number of random components.}
  \item{nRandComp}{ number of random effects in each random component of the model.}
  \item{namesRand}{ names of the random components.}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{null.deviance}{ null-deviance, deviance of the null model. The null model will only include an intercept as the estimation of the probability parameter of the conditinal beta-binomial distribution.}
  \item{null.df}{ degrees of freedom of the null model.}
  \item{Goodness.of.fit}{ p-value of the goodness of fit test.}
  \item{balanced}{ if the conditional beta-binomial response variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ maximum score number in each beta-binomial observation.}
  \item{conv}{ convergence of the methodology. If the algorithm has converged it returns "yes", otherwise "no".}
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25

Lee Y. & Nelder J. A. (1996): Hierarchical generalized linear models, \emph{Journal of the Royal Statistical
Society. Series B}, \bold{58}, 619-678

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2017): Comparison of beta-binomial regression model approaches to analyze health related quality of life data, \emph{Statistical Methods in Medical Research}, DOI: 10.1177/0962280217690413 
}

\keyword{rootSolve, numDeriv }
\seealso{
The \code{\link{multiroot}} and \code{\link{uniroot}} functions of the R-package \code{rootSolve} for the general Newton-Raphson algorithm.

\code{\link{BBmm}}.
}

\examples{
set.seed(14)

# Defining the parameters
k <- 100
m <- 10
phi <- 0.5
beta <- c(1.5,-1.1)
sigma <- 0.5

# Simulating the covariate and random effects
x <- runif(k,0,10)
X <- model.matrix(~x)
z <- as.factor(rBI(k,4,0.5,2))
Z <- model.matrix(~z-1)
u <- rnorm(5,0,sigma)


# The linear predictor and simulated response variable
eta <- beta[1]+beta[2]*x+crossprod(t(Z),u)
p <- 1/(1+exp(-eta))
y <- rBB(k,m,p,phi)
dat <- data.frame(cbind(y,x,z))
dat$z <- as.factor(dat$z)

# Apply the model
model <- BBmm(fixed.formula = y~x,random.formula = ~z,m=m,data=dat)
sum.model <- summary(model)

}