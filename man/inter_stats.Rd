\name{inter_stats}
\alias{inter_stats}

\title{
Compute Posterior Summary Statistics of the Log-Linear Parameters.
}
\description{
This function computes the posterior summary statistics of the log-linear parameters
using the MCMC output of \code{"bcct"} and \code{"bict"} objects. The posterior summary 
statistics are posterior probability, posterior mean, posterior variance and lower and 
upper limits highest posterior density intervals (HPDIs).
}
\usage{
inter_stats(object, cutoff = 0.75, n.burnin = 0, thin = 1, prob.level = 0.95)
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
  \item{prob.level}{
An optional argument giving the probability content of the HPDIs. The default value is 0.95.
}
}
\details{
This function provides an expanded version of what \code{\link{inter_probs}} provides.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
comuter memory storage.
}
\value{

This function will return an object of class \code{"interstat"} which is a list with the following 
components:

\item{term}{A vector of term labels for each parameter.}
\item{prob}{A vector of posterior probabilities for each parameter.}
\item{post_mean}{A vector of posterior means for each parameter.}
\item{post_var}{A vector of posterior variances for each parameter.}
\item{lower}{A vector of lower limits for the 100*\code{prob.level}\% HPDI for each parameter.}
\item{upper}{A vector of upper limits for the 100*\code{prob.level}\% HPDI for each parameter.}
\item{prob.level}{The argument \code{prob.level}.}

The function will only return elements in the above list if \code{prob} > \code{cutoff}.

}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\seealso{
\code{\link{bcct}},
\code{\link{bict}},
\code{\link{print.interstat}}
\code{\link{inter_probs}}
}

\examples{
set.seed(1)
## Set seed for reproducibility
data(AOH)
## Load AOH data

test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=100,prior="UIP")
## Starting from maximal model of saturated model do 100 iterations of MCMC
## algorithm.

inter_stats(test1,n.burnin=10,cutoff=0.5)
## Calculate posterior summary statistics having used a burn-in phase of 
## 10 iterations and a cutoff of 0 (i.e. display all terms with 
## non-zero posterior probability. Will get the following:

#Posterior summary statistics of log-linear parameters:
#            post_prob post_mean post_var lower_lim upper_lim
#(Intercept)         1   2.88291 0.002565   2.78778   2.97185
#alc1                1  -0.05246 0.008762  -0.27772   0.06655
#alc2                1  -0.05644 0.006407  -0.20596   0.11786
#alc3                1   0.06822 0.005950  -0.09635   0.18596
#hyp1                1  -0.53895 0.003452  -0.63301  -0.39888
#obe1                1  -0.04686 0.007661  -0.20929   0.12031
#obe2                1   0.01395 0.004024  -0.11024   0.11783
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively

}
