\name{BB}
\alias{dBB}
\alias{rBB}
\title{
The Beta-Binomial Distribution
}
\description{
Density and random generation for the beta-binomial distribution.
}
\usage{
dBB(m,p,phi)
rBB(k,m,p,phi)}

\arguments{
\item{k}{ number of simulations.}
\item{m}{ maximum socre number in each beta-binomial observation..}
\item{p}{ probability parameter of the beta-binomial distribution.}
\item{phi}{ dispersion parameter of the beta-binomial distribution.}
}

\value{
\code{dBB} gives the density of a beta-binomial distribution with the defined \code{m}, \code{p} and \code{phi} parameters.

\code{rBB} generates \code{k} random observations based on a beta-binomial distribution with the defined \code{m}, \code{p} and \code{phi} parameters.
}

\details{The beta-binomial distribution consists of a finite sum of Bernoulli dependent variables whose probability parameter is random and follows a beta distribution. Assume that we have \eqn{y_j} a set of variables, \eqn{j=1,...,m}, with m integer, that conditioned on a random variable \eqn{u}, are independent and follow a Bernoulli distribution with probability parameter \eqn{u}. On the other hand, the random variable \eqn{u} follows a beta distribution with parameter \eqn{p/phi} and \eqn{(1-p)/phi}. Namely,
\deqn{y_j \sim Ber(u), u \sim Beta(p/phi,(1-p)/phi),}
where \eqn{0<p<1} and \eqn{phi>0}. The first and second order marginal moments of this distribution are defined as
\deqn{E[y_j]=p, Var[y_j]=p(1-p),}
and correlation between observations is defined as
\deqn{Corr[y_j,y_k]=phi/(1+phi),}
where \eqn{j,k=1,...,m} are different. Consequently, \eqn{phi} can be considered as a dispersion parameter.

If we sum up all the variables we will define a new variable which follows a new distribution that is called beta-binomial distribution, and it is defined as follows. The variable \eqn{y} follows a beta-binomial distribution with parameters \eqn{m}, \eqn{p} and \eqn{phi} if
\deqn{y|u \sim Bin(m,u), u\sim Beta(p/phi,(1-p)/phi).}
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}

\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2006): Analysis of short-form-36 (SF-36): The
beta-binomial distribution approach, \emph{Statistics in Medicine}, \bold{26}, 1318-1342
}

\keyword{ stats }
\seealso{
The \code{\link{rbeta}} and \code{\link{rbinom}} functions of package \code{stats}.
}

\examples{
set.seed(12)
# We define
m <- 10     
p <- 0.4    
phi <- 1.8  

# We perform k beta-binomial simulations for those parameters.
k <- 100
bb <- rBB(k,m,p,phi)
bb
dd <- dBB(m,p,phi)

# We are going to plot the histogram of the created variable,
# and using dBB() function we are going to fit the distribution:
hist(bb,col="grey",breaks=seq(-0.5,m+0.5,1),probability=TRUE,
  main="Histogram",xlab="Beta-binomial random variable")
lines(seq(0,m),dd,col="red",lwd=4)

}
