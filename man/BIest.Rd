\name{BIest}
\alias{BIest}
\title{
Estimation of the parameters of a binomial distribution with optional dispersion parameter.
}
\description{
\code{BIest} function estimates the probability parameter of a binomial distribution for the given data and number of trials. It is possible to include a dispersion parameter in the binomial distribution and get the estimation by the method of moments or maximum quasi-likelihood approach. This function also returns the standard deviation of the estimated probability parameter and the upper and lower bounds of the 95\% confidence interval.
}
\usage{
BIest(y,m,disp=FALSE,method="MM")
}
\arguments{
\item{y}{ response variable wich follows a binomial distribution.}
\item{m}{ number of trials in each binomial observation.}
\item{disp}{ dispersion parameter of the binomial distribution. If \code{phi}=FALSE, then the binomial distribution without dispersion parameter will be considered for estimation. Default FALSE.}
\item{method}{ the method used for estimating the parameters, "MM" for the method of moments and "MLE" for maximum quasi-likelihood. Default "MM".}
}

\value{
  \code{BIest} returns an object of class "\code{BIest}".

\item{p}{ estimation of the probability parameter. Both estimating approaches, the method of moments and the maximum likelihood estimation, perform the same estimation.}
\item{pVar}{ the variance of the estimated probability parameter.}
\item{pIC.low}{ the lower bound of the 95\% confidence interval of the estimated probability parameter.}
\item{pIC.up}{ the upper bound of the 95\% confidence interval of the estimated probability parameter.}
\item{phi}{ if the \code{disp} option is TRUE, it returns the estimated value of the dispersion parameter. Default FALSE.}
\item{m}{ number of trials in each binomial observation.}
\item{balanced}{ if the data is balanced it returns "yes", otherwise "no".}
\item{method}{ the used methodology for performing the estimation of the parameters.}
}

\details{This function performs the estimation of the parameters involved in a binomial distribution for a given data.

The estimation of the probability parameter is done by either maximum likelihood approach or method of moments due to the fact that both approaches give the same estimation,
\deqn{p=sum(y)/(m*n),}
where \eqn{m} is the number of trials and \eqn{n} is the number of observations. 

If the dispersion parameter is included in the model, \code{BIest} function performs its estimation by the method of moments or maximum quasi-likelihood methodology. The method of moments is based on the variance equation of a binomial distribution with dispersion parameter
\deqn{Var[y]=phi*mp(1-p).}
The maximum quasi-likelihood approach is based on the quadratic approximation of the log-likelihood function of a binomial distribution with dispersion parameter, where the unknown term involving \eqn{phi} is estimated with the deviance of the model.

The standard deviation of the estimated probability parameter is calculated by the Fisher information, i.e., the negative of the second derivative of the log-likelihood (log-quasi-likelihood) function.}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Pawitan Y. (2001): \emph{In All Likelihood: Statistical Modelling and Inference Using Likelihood}, Oxford University Press}

\seealso{
The \code{\link{rBI}} and \code{\link{dBI}} functions of package \code{HRQoL}. The functions perform simulations and estimate the density of a binomial distribution with optional dispersion parameter for a given set of parameters.
}

\examples{
set.seed(9)
# We simulate the binomial data with some fixed parameters and
# then try to reach the same estimations.
m <- 10   
k <- 100
p <- 0.654
y <- rBI(k,m,p) #Simulations

# without dispersion parameter
BIest(y,m)

# with dispersion parameter
# estimation by method of moments.
BIest(y,m,disp=TRUE,method="MM")
# estimation by maximum quasi-likelihood.
BIest(y,m,disp=TRUE,method="MLE")
}
