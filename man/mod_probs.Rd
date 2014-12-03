\name{mod_probs}
\alias{mod_probs}

\title{
Compute Posterior Model Probabilities
}
\description{
This function computes the posterior model probabilities using the MCMC output of 
\code{"bcct"} and \code{"bict"} objects. 
}
\usage{
mod_probs(object, n.burnin = 0, scale = 0.1, best = NULL, thin = 1)
}

\arguments{
  \item{object}{
An object of class \code{"bcct"} or \code{"bict"}.
}
  \item{n.burnin}{
An optional argument giving the number of iterations to use as burn-in. 
The default value is 0.
}
  \item{scale}{
An optional argument for controlling how the posterior model probabilities are returned 
as output. The function will return details on the models with the posterior model probability 
larger than \code{scale} times the probability of the posterior modal model. The default value is
0.1.
}
  \item{best}{
An optional argument for controlling how the posterior model probabilities are returned 
as output. The function will return details on the \code{best} models with the highest posterior 
model probabilities. For example, if \code{best=4}, then details on the four models with the highest 
posterior model probabilities will be returned. The default value is \code{NULL}. If not \code{NULL}, then this argument 
takes precedent over \code{scale}.
}
  \item{thin}{
An optional argument giving the amount of thinning to use, i.e. the computations are 
based on every \code{thin}-th value in the MCMC sample. The default value is 1, i.e. no
thinning.
}
}
\details{
It will output 
only the probabilities of the "best" models, as defined by the user specifying either the \code{best} or 
\code{scale} arguments.

The use of thinning is recommended when the number of MCMC iterations and/or the number of 
log-linear parameters in the maximal model are/is large, which may cause problems with 
comuter memory storage.
}
\value{
The function will return an object of class \code{"modprobs"} which is a list containing the following 
components.

\item{table}{An object of class \code{"data.frame"} with number of rows defined by \code{scale} or 
\code{best} and columns: \code{model_formula}; giving the model (in terms of a printed formula), 
and \code{prob}; giving the posterior model probability.}

\item{totmodsvisit}{A numeric scalar giving the total number of models visited after the burn-in iterations.}

}

\author{
Antony M. Overstall \email{Antony.Overstall@glasgow.ac.uk}.
}

\seealso{
\code{\link{bcct}},
\code{\link{bict}},
\code{\link{print.modprobs}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility
data(AOH)
## Load AOH data

test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=100,prior="UIP")
## Starting from maximal model of saturated model do 100 iterations of MCMC
## algorithm.

mod_probs(object=test1,n.burnin=10,best=6)
## Using a burn-in of 10 iterations find the posterior model probabilities 
## of the 6 models with the highest posterior model probability. Will get:

#Posterior model probabilities:
#  prob    model_formula                       
#1 0.50000 ~alc + hyp + obe                    
#2 0.32222 ~alc + hyp + obe + hyp:obe          
#3 0.12222 ~alc + hyp + obe + alc:hyp + hyp:obe
#4 0.05556 ~alc + hyp + obe + alc:hyp          
#
#Total number of models visited =  4

## Note that since the chain only visited 4 models we only get probabilities
## for 4 models not 6.
}
