\name{print.summary.BImm}
\alias{print.summary.BImm}
\title{
Print a summary.BImm class model.
}
\description{
\code{print.summary.BImm} is the summary.BImm specific method fot the generic function print which prints objects returned by modelling functions.
}

\usage{
\method{print}{summary.BImm}(x, ...)
}

\arguments{
\item{x}{ a summary.BImm class model.}
\item{...}{ for extra arguments.}
}

\value{
Prints a summary.BImm object.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25

McCulloch C. E. & Searle S. R. (2001): Generalized, Linear, and Mixed Models, \emph{Jhon Wiley & Sons}

Pawitan Y. (2001): In All Likelihood: Statistical Modelling and Inference Using Likelihood, \emph{Oxford University Press}
}

\seealso{
\code{\link{BImm}}, \code{\link{summary.BImm}}
}

\examples{
set.seed(5)
# Fixing parameters for the simulation:
nObs <- 1000
m <- 10
beta <- c(1.5,-1.1)
sigma <- 0.8

# Simulating the covariate:
x <- runif(nObs,-5,5)

# Simulating the random effects:
z <- as.factor(rBI(nObs,5,0.5,2))
u <- rnorm(6,0,sigma)

# Getting the linear predictor and probability parameter.
X <- model.matrix(~x)
Z <- model.matrix(~z-1)
eta <- beta[1]+beta[2]*x+crossprod(t(Z),u)
p <- 1/(1+exp(-eta))

# Simulating the response variable
y <- rBI(nObs,m,p)

# Apply the model
model <- BImm(fixed.formula = y~x,random.formula = ~z,m=m)
sum.model <- summary(model)
print(sum.model) # or just sum.model
}
