\name{accept_rate}
\alias{accept_rate}

\title{
Compute Acceptance Rates for Metropolis-Hastings and Reversible Jump Algorithms
}
\description{
This function computes the acceptance rates of the Metropolis-Hastings and reversible 
jump algorithms from the MCMC output of \code{bcct} and \code{bict} objects.
}
\usage{
accept_rate(object)
}

\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"}. 
}
}
\details{
Acceptance rates can be used to assess the performance of MCMC methods (in particular the 
peformance of the reversible jump method, Brooks et al, 2003).
}
\value{
This function will return an object of class \code{"acceptrate"} which is a list with the 
following components.
\item{rj_ar}{Acceptance rate (as a \%) of the reversible jump algorithm.}
\item{mh_ar}{Acceptance rate (as a \%) of the Metropolis-Hastings algorithm.}
}
\references{
Brooks, S.P., Giudici, P., & Roberts, G.O. (2003) Efficient construction of reversible jump 
Markov chain Monte Carlo proposal distributions. \emph{Journal of the Royal Statistical 
Society, Series B}, \bold{65 (1)}, 3--55.
}
\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk}.
}

\seealso{
\code{\link{print.acceptrate}},
\code{\link{bcct}},
\code{\link{bict}}.

}
\examples{
set.seed(1) ## set a seed for reproducibility
data(AOH)
test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=500,prior="UIP")
## Create a bcct object for the AOH dataset for a very small number of
## iterations  (500).
accept_rate(test1)
## Calculate accept rates. Will get:

#Acceptance rate of reversible jump proposals =  32.5581 % 
#Acceptance rate of Metropolis-Hastings proposals =  76.8595 %

}
