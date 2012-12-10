\name{total_pop}
\alias{total_pop}

\title{
Evaluate Posterior Distribution of Total Population Size
}
\description{
This function uses the MCMC output of a "\code{bict}" object to derive an MCMC sample from the 
posterior distribution of the total population size.
}
\usage{
total_pop(object, n.burnin = 0, thin = 1, prob.level = 0.95)
}

\arguments{
  \item{object}{
An object of class "\code{bict}".
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
An optional argument giving the target probability content of the highest posterior density 
intervals for the total population size. The default value is 0.95.
}
}

\value{
This function will return an object of class "\code{totpop}" which is a list with the following components.

\item{TOT}{A vector of length (\code{n.sample}-\code{n.burnin}) giving the MCMC sample from the posterior distribution of the 
total population size.}
\item{int}{The 100*\code{prob.level}\% highest posterior density interval (HPDI) for the total population size.}
\item{meanTOT}{The posterior mean of the total population size.}
\item{prob.level}{The argument \code{prob.level}.}

}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\seealso{
\code{\link{bict}},
\code{\link{print.totpop}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility

data(rabbit)
## Load rabbit data

test1<-bict(formula=y~(first+second+hunter)^2,data=rabbit,n.sample=100)

## For the rabbit dataset. We do 100 iterations under the unit information 
## prior. The maximal model is the saturated model and we start from this
## model at the posterior model

tp<-total_pop(test1,n.burnin=10)
## Use a burn-in phase of 10 iterations
tp
## Print out results. Will get:
#Posterior mean of total population size = 263.1333 
#95 % highest posterior density interval for total population size = ( 207 318 )

## Could do a plot
\dontrun{plot(tp)}

## Do a summary of MCMC sample from total population size
summary(tp$TOT)
## Will get
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  207.0   238.5   256.5   263.1   289.0   327.0 

}
