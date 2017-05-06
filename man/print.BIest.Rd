\name{print.BIest}
\alias{print.BIest}
\title{
Print a BIest class model.
}
\description{
\code{print.BIest} is the BIest specific method fot the generic function print which prints objects returned by modelling functions.
}

\usage{
\method{print}{BIest}(x, ...)
}

\arguments{
\item{x}{ a BIest class model.}
\item{...}{ for extra arguments.}
}

\value{
Prints a BIest object.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Pawitan Y. (2001): \emph{In All Likelihood: Statistical Modelling and Inference Using Likelihood}, Oxford University Press}

\seealso{
\code{\link{BIest}}
}

\examples{
set.seed(9)
# Simulate a binomial distribution
y <- rBI(100,10,0.5)

# Apply the model
model <- BIest(y,10)
print(model)   # or just model

}
