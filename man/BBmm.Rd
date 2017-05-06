\name{BBmm}
\alias{BBmm}

\title{
Beta-binomial mixed-effects model
}

\description{
\code{BBmm} function performs beta-binomial mixed-effects models, i.e., it allows the inclusion of gaussian random effects in the linear predictor of a logitic beta-binomial regression model.

The structure of the random part of the model can be expecified by two different ways: (i) determining the \code{random.formula} argument, or (ii) especifying the model matrix of the random effects, \code{Z}, and determining the number of random effects in each random component, \code{nRandComp}.

The estimation of the parameters can also be done by means of two approaches: (i) BB-NR, special Newton-Raphson algorithm developed for beta-binomial mixed-effect models, and (ii) using the \code{rootSolve} R-package.
}

\usage{
BBmm(fixed.formula,random.formula,Z=NULL,nRandComp=NULL,m,data,
      method="BB-NR",maxiter=100)
}

\arguments{
\item{fixed.formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the fixed part of the model to be fitted.}
\item{random.formula}{ an object of class \code{"formula"} (or one that can be coerced to that class): a symbolic description of the random part of the model to be fitted. It must be specified in cases where the model matrix of the random effects \code{Z} is not determined.}
\item{Z}{ the design matrix of the random effects. If the \code{random.formula} argument is specified this argument should not be specified, as we will be specifying twice the random structure of the model.}
\item{nRandComp}{ the number of random effects in each random component of the model. It must be especified as a vector, where the 'i'th value must correspond with the number of random effects of the 'i'th random component. It must be only determined when we specify the random structure of the model by the model matrix of the random effects, \code{Z}.}
\item{m}{ maximum score number in each beta-binomial observation.}
\item{data}{ an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula).}
\item{method}{ the methodology for performing the estimation of the parameters. Options "BB-NR" or "rootSolve". Default "BB-NR".}
\item{maxiter}{ the maximum number of iterations in the parameters estimation algorithm. Default 100.}
}

\details{
\code{BBmm} function performs beta-binomial mixed effects models. It extends the beta-binomial logistic regression to the inclusion of random effects in the linear predictor of the model. The model is defined as, conditional on some gaussian random effects \eqn{u} the response variable \eqn{y} follows a beta-binomial distribution of parameters \eqn{m}, \eqn{p} and \eqn{phi},
\deqn{y|u \sim BB(m,p,phi), u \sim N(0,D)}
where
\deqn{log(p/(1-p))=X*beta+Z*u}
and \eqn{D} is determined by some dispersion parameters icluded in the parameter vector \eqn{theta}.

The estimation of the regression paramters \eqn{beta} and the prediction of the random effects \eqn{u} is done via the extended likelihood, where the marginal likelihood is approximated to the h-likelihood by a first order Laplace approximation,
\deqn{h=f(y|beta,u,theta)+f(u|theta)}
The previous formula do not have a closed form and numerical methods are needed for the estimation precedure. Two approches are available in the function in order to perform the fixed and random effects estimation: (i) A special case of a Newton-Raphson algorithm developed for beta-binomial mixed-effects model estimations, and (ii) the general Newton-Raphson algorithm available in R-package \code{rootSolve}.


On the other hand, the estimation of dispersion parameters by the h-likelihood can be bias due to the estimation of the regression parameters. Consequenlty, a penalization of the h-likelihood must be performed to get an unbiased profile h-likelihood of the dispersion parameters. Lee and Nelder (1996) proposed the adjusted profile h-likelihood, where the following penalization is applied,
\deqn{h(theta)=h+(1/2)*log[det(2\pi H^{-1})]}
where \eqn{H} is the Hessian matrix of the model, and the terms \eqn{beta} and \eqn{u} involved in the previous formula are replaced by their estimates.

The method iterates between both estimation processes until convergence is reached.
}

\value{
  \code{BBmm} returns an object of class "\code{BBmm}".
  
  The function \code{summary} (i.e., \code{summary.BBmm}) can be used to obtain or print a summary of the results..

  \item{fixed.coef}{ estimated value of the fixed coefficients of the regression.}
  \item{fixed.vcov}{ variance and covariance matrix of the estimated fixed coefficients of the regression.}
  \item{random.coef}{ predicted random effects of the regression.}
  \item{sigma.coef}{ estimated value of the random effects variance parameters.}
  \item{sigma.var}{ variance of the estimated value of the random effects variance parameters.}
  \item{phi.coef}{ estimated value of the dispersion parameter of the conditional beta-binomial distribution.}
  \item{psi.coef}{ estimated value of the logrithm of the dispersion parameter of the conditional beta-binomial distribution.}
  \item{psi.var}{ variance of the estimation of the logarithm of the conditional beta-binomial distribution.}
  \item{fitted.values}{ the fitted mean values of the probability parameter of the conditional beta-binomial distribution.}
  \item{conv}{ convergence of the methodology. If the method has converged it returns "yes", otherwise "no".}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{null.deviance}{ null-deviance, deviance of the null model. The null model will only include an intercept as the estimation of the probability parameter.}
  \item{null.df}{ degrees of freedom of the null model.}
  \item{nRand}{ number of random effects.}
  \item{nComp}{ number of random components.}
  \item{nRandComp}{ number of random effects in each random component of the model.}
  \item{namesRand}{ names of the random components.}
  \item{iter}{ number of iterations in the estimation method.}
  \item{nObs}{ number of observations in the data.}
  \item{y}{ dependent response variable in the model.}
  \item{X}{ model matrix of the fixed effects.}
  \item{Z}{ model matrix of the random effects.}
  \item{D}{ variance and covariance matrix of the random effects.}
  \item{balanced}{ if the response dependent variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ maximum score number in each beta-binomial observation.}
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

Lee Y. & Nelder J. A. (1996): Hierarchical generalized linear models, \emph{Journal of the Royal Statistical
Society. Series B}, \bold{58}, 619-678

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2017): Comparison of beta-binomial regression model approaches to analyze health related quality of life data, \emph{Statistical Methods in Medical Research}. DOI: 10.1177/0962280217690413 
}

\keyword{rootSolve, numDeriv }
\seealso{
The \code{\link{multiroot}} and \code{\link{uniroot}} functions of the R-package \code{rootSolve} for the general Newton-Raphson algorithm.
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
model
}