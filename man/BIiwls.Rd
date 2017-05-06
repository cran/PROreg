\name{BIiwls}
\alias{BIiwls}

\title{
Iterative Weighted Least Squares (IWLS) method for binomial logistic regression
}

\description{
\code{BIiwls} performs estimation of the coefficients of binomial logistic regressions by iterative wieghted least squares method.
}

\usage{
BIiwls(y,X,m,maxiter)
}

\arguments{
\item{y}{ response dependent variable which follows a binomial distribution.}
\item{X}{ model matrix.}
\item{m}{ number of trials in each binomial observation.}
\item{maxiter}{ maximum number of iterations in the method.}
}

\value{
\item{beta}{ maximum likelihood estimations of the logistic regression coefficients.}
\item{vcov}{ variance-covariance matrix of the estimated regression coefficients.}
\item{iter}{ number of iterations of the algorithm.}
}

\details{
The iterative weighted least squares (IWLS) is a general algorithm to find the maximum likelihood estimations (mle) and standard deviations in generalized linear models. There are several ways to derive it, but the one that has been developed in this function is via the Newton-Raphson method. It consists of making a Taylor expansion in the score function, the first derivate of the log-likelihood, around the mle. This especific IWLS, \code{BIiwls}, has been developed to find out the mle and the standard errors in logistic regression by the introduction of a dependent variable, a matrix model of the regression covariates and the number of trials of the binomial dependent variable.
}
\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Pawitan Y. (2001): In All Likelihood: Statistical Modelling and Inference Using Likelihood, \emph{Oxford University Press}
}

\examples{
# We create a variable and a matrix model that consists of
# two covariates.
k=1000                    
m=10 
maxiter=100

y <- rbinom(k,m,0.87) 
x1 <- rnorm(k,1,50)
x2 <- rnorm(k,30,9) 
X <- cbind(1,x1,x2)

BIiwls(y,X,m,maxiter)
}
