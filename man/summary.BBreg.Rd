\name{summary.BBreg}
\alias{summary.BBreg}

\title{
Summarizes a BBreg class model.
}

\description{
\code{summary.BBreg} si the BBreg specific method for the generic function \code{summary} which summarizes objects returned by modelling functions.
}

\usage{
\method{summary}{BBreg}(object, ...)
}

\arguments{
\item{object}{ a BBreg class model.}
\item{...}{ for extra arguments.}
}

\details{
\code{summary.BBreg} summarizes all the relevant information about the estimation of the parameters in a BBreg class model.

The function performs statistical significance hypothesis about the estimated regression parameters based on the normal distribution of the estimates. It also performs a goodness of fit test based on the difference between the calculated deviance of the model and the null deviance or deviance of the null model, which it is suppose to follow a Chi-square distribution with degrees of freedom equal to the difference in degrees of freedom of the models.
}

\value{
  \code{summary.BBreg} returns an object of class "\code{summary.BBreg}".
  
  \item{coefficients}{ a table with all the relevant information about the significance of the regression coefficients of the model. It includes the estimations, the standard errors of the estimations, the test-statistics and the p-values.}
  \item{psi.table}{ a table which includes the estimation and standard errors of the logarithm of the dispersion parameter of the conditional beta-binomial distribution.}
  \item{deviance}{ deviance of the model.}
  \item{df}{ degrees of freedom of the model.}
  \item{null.deviance}{ null-deviance, deviance of the null model.}
  \item{null.df}{ degrees of freedom of the null model.}
  \item{Goodness.of.fit}{ p-value of the goodness of fit test.}
  \item{iter}{ number of iterations in the estimation method.}
  \item{X}{ the model matrix.}
  \item{y}{ the dependent variable in the model.}
  \item{balanced}{ if the response variable is balanced it returns "yes", otherwise "no".}
  \item{m}{ number of trials in each binomial observation.}
  \item{nObs}{ number of observations.}
  \item{m}{ number of trials in each observation.}
  \item{balanced}{ if the response beta-binomial variable is balanced it returns "yes", otherwise "no".}
  \item{conv}{ convergence of the methodology. If the algorithm has converged it returns "yes", otherwise "no".}
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
sum.model <- summary(model)
}
