\name{BBest}
\alias{BBest}
\title{Estimation of the parameters of a beta-binomial distribution}

\description{
This function performs the estimation of the parameters of a beta-binomial distribution for the given data and maximum score number in each observation. 

There are two different approaches available for performing the estimation of the parameters: (i) Method of moments, and, (ii) maximum likelihood approach.
}

\usage{
BBest(y,m,method="MM")
}

\arguments{
\item{y}{ response variable which folloes a beta-binomial distribution.}
\item{m}{ maximum score number in each beta-binomial observation.}
\item{method}{ the method used for performing the estimation of the probability and dispersion parameters of a beta-binomial distribution, "MM" for the method of moments and "MLE" for maximum likelihood estimation. Default "MM".}
}

\value{
\code{BBest} returns an object of class "\code{BBest}".

The function \code{summary} (i.e., \code{summary.BBest}) can be used to obtain or print a summary of the results.

\item{p}{ estimated probability parameter of the beta-binomial distribution.}
\item{phi}{ estimated dispersion parameter of the beta-binomial distribution.}
\item{pVar}{ variance of the estimation of the probability parameter \code{p}.}
\item{psi}{ estimation of the logarithm of the dispersion parameter \code{phi}.}
\item{psiVar}{ variance of the estimation of the logarithm of the dispersion parameter \code{psi}.}
\item{m}{ maximum score number in each beta-binomial observation.}
\item{balanced}{ if the response variable is balanced it returns "yes", otherwise "no".}
\item{method}{ the used approach for performing the estimations.}
}

\details{
\code{BBest} function performs estimations in the parameters of a beta-binomial distribution for the given data. The estimations can be performed using two different approaches, the methods of moments and the maximum likelihood estimation approach.

The density function of a given observation \eqn{y}{y} that follows a beta-binomial distribution with paramters \eqn{m}, \eqn{p} and \eqn{phi} is defined as
\deqn{f(y)=[\Gamma(1/phi)*\Gamma(p/phi+m)*\Gamma((1-p)/phi+m)]/[\Gamma(1/phi+m)*\Gamma(p/phi)*\Gamma((1-p)/phi)].}

The first and second order moments are defined as
\deqn{E[y]=mp}
\deqn{Var[y]=mp(1-p)[1+(n-1)phi/(1+phi)].}

Hence, if \eqn{y=(y_1,...,y_n)} is the given data, we can conclude the method of moments from the previous as
\deqn{p=E/m,}
\deqn{phi=[V-mp(1-p)]/[mp(1-p)m-V],}
where \eqn{E} is the sample mean and \eqn{V} is the sample variance.

On the other hand, the maximum likelihood estimation of both parameters consits of solving the derivative of the log-likelihood defined by the density function with respect to each parameter and equaling them to zero. An iterative algorithm is needed for both parameter estimation as the score equations the parameters depend each other.

The variance of the estimation of the probability parameter of the beta-binomial distribution for the given data set is computed by the inverse of the Fisher information, i.e., the inverse of the negative second derivate of the log-likelihodd remplacing \eqn{p} by its estimation.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2006): Analysis of short-form-36 (SF-36): The beta-binomial distribution approach, \emph{Statistics in Medicine}, \bold{26}, 1318-1342
}

\examples{
# We simulate 1000 observations of a beta-binomial distribution
# for the fixed paramters.
m <- 10
k <- 1000
p <- 0.7
phi <- 1.6

set.seed(5)
y <- rBB(k,m,p,phi)

# Performing the estimation of the parameters

# Method of moments:
MM <- BBest(y,m)
MM

# Maximum likelihood approach
MLE <- BBest(y,m,method="MLE")
MLE

}
