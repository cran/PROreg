\name{SF36rec}
\alias{SF36rec}
\title{
Short Form-36 Health Survey recode
}

\description{
The \code{SF36rec} function recodes to a binomial form the 0-100 original standardized scores of the dimensions provided by the Short Form-36 Health Survey (SF-36) based on \cite{Arostegui et al. (2013)}.
}

\usage{
SF36rec(x,k)
}

\arguments{
\item{x}{ the 0-100 scale standardized dimension of the SF-36. It must be numeric and bounded between 0 and 100.}
\item{k}{ an integer from 1 to 8 that defines which SF-36 dimension belongs \code{x}. These are the dimensions depending on the \code{k} value:

\code{k=1} -> Physical functioning

\code{k=2} -> Role physical

\code{k=3} -> Body pain

\code{k=4} -> General health

\code{k=5} -> Vitality

\code{k=6} -> Social functioning

\code{k=7} -> Role emotional

\code{k=8} -> Mental health
}
}

\value{
The score values of the recoded dimension.
}

\details{
The Short Form-36 Health Survey is a commonly used technique to measure the Health Related Quality of Life (HRQoL) in chronich diseases. It was developed within the Medical Outcomes Study \cite{(Ware et al. (1993))}. It measures generic HRQoL concepts and provides an objective way to measure HRQoL from the patients point of view by scoring standardized responses to standardized questions. The validity and reability of this instrument has been broadly tested \cite{(Stansfeld et al. (1997))}. The SF-36 has 36 items, with different answer options. It was constructed to respresent eight health dimensions, which are \emph{physical functioning} (PF), \emph{role physical} (RP), \emph{body pain} (BP), \emph{general health} (GH), \emph{vitality} (VT), \emph{social functioning} (SF), \emph{role emotional} (RE) and \emph{ mental health} (MH). Each item is assigned to a unique helath dimension. Each of the multi-item dimensions contains two to ten items. The first four dimensions are mainly physical, whereas the last four measure mental aspects of HRQoL. The resulting raw scores are tipically transformed to standardized scale scores from 0 to 100, where a higher score indicates a better health status.

\cite{Arostegui et al. (2013)} proposed a recoding methodology of the SF-36 standardized scores to a binomial form in order to apply the beta-binomial distribution. The method was mainly based on the possible number of values each dimension can obtain, which comes from the number of items related with the construction of each dimension.

The \code{SF36rec} function performes the cited recoding methodology to the specified SF-36 dimension. It has two inputs. The first one is the dimension that will be recoded, and the second one identifies which SF-36 dimension is. 
}

\author{
J. Najera-Zuloaga

D.-J. Lee

I. Arostegui
}
\references{
Arostegui I., Nunez-Anton V. & Quintana J. M. (2013): On the recoding of continuous and bounded indexes to a binomial form: an application to quality-of-life scores, \emph{Journal of Applied Statistics}, \bold{40}, 563-583

Ware J. E., Snow K. K., Kosinski M. A. & Gandek B. (1993): \emph{SF-36 Health Survey, Manual and Interpretation Guides}. The Health Institute, New England Medical Center.

Stansfeld S. A., Roberts R. & Foot S. P. (1997): Assessing the validity of the SF-36 general health survey.
Quality of Life Research, \bold{6}, 217-224.
}


\examples{
set.seed(2)
# We simulate a variable bounded between 0 and 100.
BodyPain <- runif(1000,0,100)

# We specify that the simulated dimension corresponds
# with body pain dimension.
k <- 3

# We perform the recoding.
BodyPain.rec <- SF36rec(BodyPain,k)
}
