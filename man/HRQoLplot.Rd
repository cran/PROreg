\name{HRQoLplot}
\alias{HRQoLplot}
\title{
Spider plot of the dimensions of the Short Form-36 Health Survey
}
\description{
This function creates a spider plot with the 8 health realted quality of life dimensions provided by the Short Form-36 Health Survey.
}
\usage{
HRQoLplot(data,legend=FALSE,title="Short Form-36 Health Survey",
  dimlabel.cex=NULL,legend.cex=1,linewidth=3,title.cex=1,lty=1)
}
\arguments{
\item{data}{ a data frame with each column relative to the observations of each SF-36 dimension. The columns of the data frame must be introduced in the following order:
\enumerate{
\item column -> Physical Functioning
\item column -> Role Physical
\item column -> Body Pain
\item column -> General Health
\item column -> Vitality
\item column -> Social Functioning
\item column -> Role Emotional
\item column -> Mental Health
}
}
\item{legend}{ logical parameter, if TRUE the legend with the name of the rows of the data will appear. Default FALSE.}
\item{title}{ the title of the plot. Default "Short Form-36 Health Survey".}
\item{dimlabel.cex}{ font size magnification for the labels of the dimension in the plot. If NULL, the font size is fixed at text()'s default. Default NULL.}
\item{legend.cex}{ font size of legend text(). Default 1.}
\item{linewidth}{ the width of the lines of the plot. Default 3.}
\item{title.cex}{ the font size of the title. Default 1.}
\item{lty}{ the line type of the plot and the legend. Default 1.}
}

\details{
The Short Form-36 Health Survey is a commonly used technique to measure the Health Related Quality of Life (HRQoL) in chronich diseases. It was developed within the Medical Outcomes Study \cite{(Ware et al. (1993))}. It measures generic HRQoL concepts and provides an objective way to measure HRQoL from the patients point of view by scoring standardized responses to standardized questions. The validity and reability of this instrument has been broadly tested \cite{(Stansfeld et al. (1997))}. The SF-36 has 36 items, with different answer options. It was constructed to respresent eight health dimensions, which are \emph{physical functioning} (PF), \emph{role physical} (RP), \emph{body pain} (BP), \emph{general health} (GH), \emph{vitality} (VT), \emph{social functioning} (SF), \emph{role emotional} (RE) and \emph{ mental health} (MH). Each item is assigned to a unique helath dimension. Each of the multi-item dimensions contains two to ten items. The first four dimensions are mainly physical, whereas the last four measure mental aspects of HRQoL. The resulting raw scores are tipically transformed to standardized scale scores from 0 to 100, where a higher score indicates a better health status.

\cite{Arostegui et al. (2013)} proposed a recoding methodology for the Short Form-36 Health Survey (SF-36) dimensions in order to apply a beta-binomial distribution. The \code{HRQoLplot} function plots the SF-36 dimensions scores in a spider plot. Each axis of the plot refers to an expecific SF-36 dimension. Hence, the order of the dimensions in the data frame object of the function has been stablished as it has been explained in Arguments section. Each observation of the data frame, the value of each observation in all the dimensions, is drawn with a line of a different color in the plot. The plot shows the name of each dimension and the maximum number of scores each dimension can obtain in each axis of the plot.
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui

This function depends on the function \code{\link[fmsb]{radarchart}} of the package \code{fmsb} created by Minato Nakazawa.

}

\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2013): On the recoding of continuous and bounded indexes to a binomial form: an application to quality-of-life scores, \emph{Journal of Applied Statistics}, \bold{40}, 563-583
}

\keyword{ fmsb }
\seealso{
As it is said in the author section, the function depends on the function \code{\link[fmsb]{radarchart}} of the package \code{fmsb}
}
\examples{
set.seed(5)
# We insert the columns in the order that has been determined:
n <- c(20,4,9,20,20,8,3,13)
k=3
p=runif(8,0,1)
phi <- runif(8,1,3)
dat <- data.frame(
  PF=rBB(k,n[1],p[1],phi[1]),
  RP=rBB(k,n[2],p[2],phi[2]),
  BP=rBB(k,n[3],p[3],phi[3]),
  GH=rBB(k,n[4],p[4],phi[4]),
  VT=rBB(k,n[5],p[5],phi[5]),
  SF=rBB(k,n[6],p[6],phi[6]),
  RE=rBB(k,n[7],p[7],phi[7]),
  MH=rBB(k,n[8],p[8],phi[8]))

rownames(dat) <- c("ID1", "ID2", "ID3")
HRQoLplot(dat,TRUE)
}
