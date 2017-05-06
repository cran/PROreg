\name{summary.BBest}
\alias{summary.BBest}
\title{Summarizes a BBest class model.}

\description{
\code{summary.BBest} si the BBest specific method for the generic function \code{summary} which summarizes objects returned by modelling functions.
}

\usage{
\method{summary}{BBest}(object, ...)
}

\arguments{
\item{object}{ a BBest class model.}
\item{...}{ for extra arguments.}
}

\value{
\code{summary.BBest} returns an object of class "\code{summary.BBest}".
  
\item{coefficients}{ a table with the estimated parameters is in the BBest class model.}
\item{p.coefficients}{ a summarized table of the estimation of the probability parameter of the beta-binomial distribution. The table contents the estimation of the probability parameter and the standard errors of the estimations.}
\item{psi.coefficients}{ a summarized table of the estimation of the logarithm of the dispersion parameter of the beta-binomial distribution. The table contents the estimation of the logarithm of the dispersion parameter and the standard errors of the estimations.}
\item{m}{ the maximum score number in each beta-binomial observation.}
}

\details{
\code{summary.BBest} summarizes all the relevant information about the estimation of the parameters in a BBest class model. 
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2006): Analysis of short-form-36 (SF-36): The beta-binomial distribution approach, \emph{Statistics in Medicine}, \bold{26}, 1318-1342
}

\seealso{
\code{\link{BBest}}
}


\examples{
set.seed(9)
# Simulate a binomial distribution
y <- rBB(100,10,0.5,2)

# Apply the model
model <- BBest(y,10)
sum.model <- summary(model)
}
