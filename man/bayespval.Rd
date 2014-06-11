\name{bayespval}
\alias{bayespval}

\title{
Compute Bayesian p-value
}
\description{
This function will compute the Bayesian (or posterior predictive) p-value. This can be used 
as a diagnostic tool to check model adequacy. Additionally this function outputs 
predictions from the model which can also be used in other assessments of model adequacy.
}
\usage{
bayespval(object, n.burnin = 0, thin = 1, statistic = "X2")
}

\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"} object.
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
  \item{statistic}{
An optional argument giving the discrepancy statistic to use for calculating the Bayesian p-value. It can be one of 
\code{c("X2","FreemanTukey","deviance")} which correspond to the different statistics: 
\code{"X2"} = Chi-squared statistic, \code{"FreemanTukey"} = Freeman-Tukey statistic, 
\code{"deviance"} = deviance statistic. See Overstall & King (2014), and references 
therein, for descriptions of these statistics.
}
}

\details{See Gelman et al (2004, Chapter 6) for more details on Bayesian p-values and see 
Overstall & King (2014), and references therein, for details of their application to 
contingency tables.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
comuter memory storage.
}

\value{
The function will produce an object of class \code{"pval"} which is a list 
with the following components.
\item{PRED}{An (\code{n.sample}-\code{n.burnin}) by n* matrix where (n* is the number of observed cell counts) containing the
predictions of the observed cell counts.}
\item{Tpred}{A vector of length (\code{n.sample}-\code{n.burnin}) containing the discrepancies between the predicted cell counts 
and their means.}
\item{Tobs}{A vector of length (\code{n.sample}-\code{n.burnin}) containing the discrepancies between the observed cell counts 
and their means.}
\item{pval}{A scalar giving the Bayesian p-value, i.e. the proportion of \code{Tpred}>\code{Tobs}.}
\item{statnum}{A numeric scalar identifying which statistic is used.}
\item{statistic}{A character string identifying which statistic is used.}
\item{thin}{The value of the argument \code{thin}.}
}
\references{
Gelman, A., Carlin, J.B., Stern, H.S. & Rubin, D.B. (2004) \emph{Bayesian Data Analysis}, 2nd edition, Chapman & Hall.

Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Journal of Statistical Software}, \bold{58 (7)}, 
1--27. \url{http://www.jstatsoft.org/v58/i07/}
}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}

\seealso{
\code{\link{bict}},
\code{\link{bcct}},
\code{\link{print.pval}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility
data(spina)
## Load spina data

test1<-bict(formula=y~(S1+S2+S3+eth)^2,data=spina,n.sample=50,prior="UIP")
## Do 50 iterations starting at maximal model containing all two-way interactions.

test1p<-bayespval(object=test1,statistic="FreemanTukey",n.burnin=5)
## Use the Freeman-Tukey statistic and a burn-in phase of 5 iterations.
test1p
## Will get following output

#Under the Freeman-Tukey statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  2.812   4.695   5.190   5.777   6.405  14.490 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  4.566   4.861   5.197   5.430   6.108   6.460 
#
#Bayesian p-value =  0.4667

## Can do a plot

\dontrun{plot(test1p)}


}
