\name{sub_model}
\alias{sub_model}

\title{
Compute Posterior Summary Statistics for (Sub-) Models
}
\encoding{UTF-8}
\description{
This function computes posterior summary statistics for (sub-) models 
using the MCMC output of \code{"bcct"} and \code{"bict"} objects. 
}
\usage{
sub_model(object, formula = NULL, order = 1, n.burnin = 0, thin = 1, 
prob.level = 0.95, statistic = "X2")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"}.
}
  \item{formula}{
An optional argument of class \code{"formula"}: a symbolic description of the 
model of interest. The default value is \code{NULL}. If not \code{NULL} then 
this argument takes precedent over \code{order}.
}
  \item{order}{
A scalar argument identifying the model for which to compute summary 
statistics. The function will compute statistics for the model with the 
\code{order}-th largest posterior model probability. The default value is 1, 
meaning that, by default, the function will compute summary statistics for 
the posterior modal model.  
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
An optional argument giving the probability content of the highest posterior density intervals (HPDIs). 
The default value is 0.95.
}
  \item{statistic}{
An optional argument giving the discrepancy statistic to use for calculating the Bayesian p-value. It can be one of 
\code{c("deviance","X2","FreemanTukey")} which correspond to the different statistics: 
\code{"deviance"} = deviance statistic, \code{"X2"} = Chi-squared statistic, 
\code{"FreemanTukey"} = Freeman-Tukey statistic. See Overstall & King (2013), and references 
therein, for descriptions of these statistics.
}
}
\details{
If the MCMC algorithm does not visit the model of interest in the thinned MCMC sample, 
after burn-in, then an error message will be returned.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
computer memory storage.
}
\value{
This function will return an object of class \code{"submod"} which is a list with the following 
components. Note that, unless otherwise stated, all components are conditional on the model of 
interest. 

\item{term}{A vector of term labels for each log-linear parameter.}


\item{post_prob}{A scalar giving the posterior model probability for the model of interest.}

\item{post_mean}{A vector of posterior means for each of the log-linear parameters.}

\item{post_var}{A vector of posterior variances for each of the log-linear parameters.}

\item{lower}{A vector of lower limits for the 100*\code{prob.level}\% HPDI for each log-linear parameter.}

\item{upper}{A vector of upper limits for the 100*\code{prob.level}\% HPDI for each log-linear parameter.}

\item{prob.level}{The argument \code{prob.level}.}

\item{order}{The ranking of the model of interest in terms of posterior model probabilities.}

\item{formula}{The formula of the model of interest.}

\item{BETA}{A matrix containing the sampled values of the log-linear parameters, where the 
number of columns is the number of log-linear parameters in the model of interest.} 

\item{SIG}{A vector containing the sampled values of sigma^2 under the Sabanes-Bové & Held 
prior. If the unit information prior is used then the components of this vector will be one.} 

If \code{object} is of class \code{"bict"}, then \code{sub_model} will also return the following 
component.

\item{Y0}{A matrix (with k columns) containing the sampled values of the missing and censored 
cell counts, where k is the total number of missing and censored cell counts.}
}
\references{
Overstall, A.M. & King, R. (2013) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}

\seealso{
\code{\link{bcct}},
\code{\link{bict}},
}
\examples{
set.seed(1)
## Set seed for reproducibility.

data(AOH)
## Load the AOH data

test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=100,prior="UIP")
## Let the maximal model be the saturated model. Starting from the 
## posterior mode of the maximal model do 100 iterations under the unit 
## information prior.

test1sm<-sub_model(object=test1,order=1,n.burnin=10)
## Obtain posterior summary statistics for posterior modal model using a 
## burnin of 10.

test1sm

#Posterior model probability =  0.5
#
#Posterior summary statistics of log-linear parameters:
#            post_mean post_var lower_lim upper_lim
#(Intercept)  2.907059 0.002311   2.81725   2.97185
#alc1        -0.023605 0.004009  -0.20058   0.06655
#alc2        -0.073832 0.005949  -0.22995   0.10845
#alc3         0.062491 0.006252  -0.09635   0.18596
#hyp1        -0.529329 0.002452  -0.63301  -0.43178
#obe1         0.005441 0.004742  -0.12638   0.12031
#obe2        -0.002783 0.004098  -0.17082   0.07727
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively
#
#Under the X2 statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  11.07   19.76   23.34   24.47   29.04   50.37 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  30.82   34.78   35.74   36.28   37.45   42.49 
#
#Bayesian p-value =  0.0444

}

