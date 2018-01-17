\name{total_pop}
\alias{total_pop}

\title{
Evaluate Posterior Distribution of Total Population Size
}
\description{
This function uses the MCMC output of a \code{"bict"} object to derive an MCMC sample from the 
posterior distribution of the total population size.
}
\usage{
total_pop(object, n.burnin = 0, thin = 1, prob.level = 0.95)
}

\arguments{
  \item{object}{
An object of class \code{"bict"}.
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

\details{
The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
comuter memory storage.
}

\value{
This function will return an object of class \code{"totpop"} which is a list with the following components.

\item{TOT}{A vector of length (\code{n.sample}-\code{n.burnin}) giving the MCMC sample from the posterior distribution of the 
total population size.}
\item{int}{The 100*\code{prob.level}\% highest posterior density interval (HPDI) for the total population size.}
\item{meanTOT}{The posterior mean of the total population size.}
\item{prob.level}{The argument \code{prob.level}.}

}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\seealso{
\code{\link{bict}},
\code{\link{print.totpop}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility

data(spina)
## Load spina data

test1<-bict(formula=y~(S1+S2+S3+eth)^2,data=spina,n.sample=100,prior="UIP")

## For the spina dataset. We do 100 iterations under the unit information 
## prior. The maximal model is the model with two-way interactions and we 
## start from this model at the posterior model

tp<-total_pop(test1,n.burnin=10)
## Use a burn-in phase of 10 iterations
tp
## Print out results. Will get:

#Posterior mean of total population size = 727.0667 
#95 % highest posterior density interval for total population size = ( 706 757 )

## Could do a plot
\dontrun{plot(tp)}

## Do a summary of MCMC sample from total population size
summary(tp$TOT)
## Will get

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  697.0   716.2   727.0   727.1   735.8   763.0

}
