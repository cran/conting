\name{inter_probs}
\alias{inter_probs}

\title{
Calculate Posterior Probability of Each Term
}
\description{
This function computes the posterior probability of each 
term using the MCMC output of \code{"bcct"} and \code{"bict"} objects. 
}
\usage{
inter_probs(object, cutoff = 0.75, n.burnin = 0, thin = 1)
}

\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"}.
}
  \item{cutoff}{
An optional argument giving the cutoff posterior probability for displaying posterior 
summary statistics of the log-linear parameters. Only those log-linear parameters with 
a posterior probability greater than \code{cutoff} will be returned as part of the output. 
The default value is 0.75.
}
  \item{n.burnin}{
An optional argument giving the number of iterations to use as burn-in. 
The default value is 0.
}
  \item{thin}{
An optional argument giving the amount of thinning to use, i.e. the computations are 
based on every \code{thin}-th value in the MCMC sample. The default value is 1, i.e. no
thinning.
}
}
\details{
This function provides a scaled back version of what \code{\link{inter_stats}} provides.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
computer memory storage.
}
\value{
This function returns an object of class \code{"interprob"} which is a list with the 
following components. 

\item{term}{A vector of term labels.}

\item{prob}{A vector of posterior probabilities.}

\item{thin}{The value of the argument \code{thin}.}

The function will only return elements in the above list if \code{prob} > \code{cutoff}.

}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}

\seealso{
\code{\link{bcct}},
\code{\link{bict}},
\code{\link{print.interprob}},
\code{\link{inter_stats}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility
data(AOH)
## Load AOH data

test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=100,prior="UIP")
## Starting from maximal model of saturated model do 100 iterations of MCMC
## algorithm.

inter_probs(test1,n.burnin=10,cutoff=0)
## Calculate posterior probabilities having used a burn-in phase of 
## 10 iterations and a cutoff of 0 (i.e. display all terms with 
## non-zero posterior probability). Will get the following:

#Posterior probabilities of log-linear parameters:
#            post_prob
#(Intercept)    1.0000
#alc            1.0000
#hyp            1.0000
#obe            1.0000
#alc:hyp        0.1778
#alc:obe        0.0000
#hyp:obe        0.4444
#alc:hyp:obe    0.0000

## Note that the MCMC chain (after burn-in) does not visit any models 
## with the alc:obe or alc:hyp:obe interactions.
}