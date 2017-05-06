\name{print.summary.BBreg}
\alias{print.summary.BBreg}
\title{
Print a summary.BBreg class model.
}
\description{
\code{print.summary.BBreg} is the summary.BBreg specific method fot the generic function print which prints objects returned by modelling functions.
}

\usage{
\method{print}{summary.BBreg}(x, ...)
}

\arguments{
\item{x}{ a summary.BBreg class model.}
\item{...}{ for extra arguments.}
}

\value{
Prints a summary.BBreg object.
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

\seealso{
\code{\link{BBreg}},\code{\link{summary.BBreg}}
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
print(sum.model) # or just sum.model
}
