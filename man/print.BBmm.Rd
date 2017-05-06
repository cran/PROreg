\name{print.BBmm}
\alias{print.BBmm}
\title{
Print a BBmm class model.
}
\description{
\code{print.BBmm} is the BBmm specific method fot the generic function print which prints objects returned by modelling functions.
}

\usage{
\method{print}{BBmm}(x, ...)
}

\arguments{
\item{x}{ a BBmm class model.}
\item{...}{ for extra arguments.}
}

\value{
Prints a BBmm object.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Breslow N. E. & Calyton D. G. (1993): Approximate Inference in Generalized Linear Mixed Models, \emph{Journal of the American Statistical Association}, \bold{88}, 9-25

Lee Y. & Nelder J. A. (1996): Hierarchical generalized linear models, \emph{Journal of the Royal Statistical
Society. Series B}, \bold{58}, 619-678

Najera-Zuloaga J., Lee D.-J. & Arostegui I. (2017): Comparison of beta-binomial regression model approaches to analyze health related quality of life data, \emph{Statistical Methods in Medical Research}, DOI: 10.1177/0962280217690413 
}

\seealso{
\code{\link{BBmm}}
}

\examples{
set.seed(14)

# Defining the parameters
k <- 100
m <- 10
phi <- 0.5
beta <- c(1.5,-1.1)
sigma <- 0.5

# Simulating the covariate and random effects
x <- runif(k,0,10)
X <- model.matrix(~x)
z <- as.factor(rBI(k,4,0.5,2))
Z <- model.matrix(~z-1)
u <- rnorm(5,0,sigma)


# The linear predictor and simulated response variable
eta <- beta[1]+beta[2]*x+crossprod(t(Z),u)
p <- 1/(1+exp(-eta))
y <- rBB(k,m,p,phi)
dat <- data.frame(cbind(y,x,z))
dat$z <- as.factor(dat$z)

# Apply the model
model <- BBmm(fixed.formula = y~x,random.formula = ~z,m=m,data=dat)
print(model) # or just model
}
