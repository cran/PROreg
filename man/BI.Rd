\name{BI}
\alias{dBI}
\alias{rBI}
\title{
The Binomial distribution with optional Dispersion Parameter
}
\description{
Density and random generation for the binomial distribution with optional dispersion parameter.
}
\usage{
dBI(m,p,phi)
rBI(k,m,p,phi)}

\arguments{
\item{k}{ number of simulations.}
\item{m}{ number of trials in each binomial observation.}
\item{p}{ probability parameter of the binomial distribution.}
\item{phi}{ dispersion parameter of the binomial distribution. If \code{phi}=1, the binomial distribution without dispersion parameter will be considered. Default 1.}
}

\value{
\code{dBI} gives the density of the binomial distribution for those \code{m}, \code{p} and \code{phi} parameters.

\code{rBI} generates \code{k} random observations based on a binomial distribution for those \code{m}, \code{p} and \code{phi} parameters.}

\details{The binomial distribution belongs to the exponential family of distributions. Consequenlty, although the usual binomial distribution only consists of two paramters, an additional dispersion parameter can be included. The inclusion of a dispersion parameter softens the relationship between the expectation and variance that the binomial distribution keeps, i.e. the model allows overdispersion to be included,
\deqn{E[y]=mp, Var[y]=phi*mp(1-p).}
The density function of the binomial distribution with dispersion parameter is based on the exponential family approach and it is defined as
\deqn{f(y)=exp\{[y*log(p/(1-p))+m*log(1-p)]/phi+c(y,phi)\},}
where \eqn{c()} is a function that it is approximated with the deviance of the model by quadratic approximations of the log-likelihood function.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Pawitan Y. (2001): \emph{In All Likelihood: Statistical Modelling and Inference Using Likelihood}. Oxford University Press
}

\keyword{ stats }
\seealso{
The \code{\link{rbinom}} functions of package \code{stats}. This function performs simulations based on a binomial distribution without dispersion parameter.
}
\examples{
k <- 1000
m <- 10
p <- 0.765
phi <- 4.35

#simulating
y <- rBI(k,m,p,phi)
y

#density function
d <- dBI(m,p,phi)
d

#plot the simulated variable and fit the density
hist(y,col="lightgrey")
lines(seq(0,m),k*d,col="blue",lwd=2)
}
