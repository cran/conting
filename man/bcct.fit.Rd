\name{bcct.fit}
\alias{bcct.fit}

\title{
Bayesian Analysis of Complete Contingency Tables
}
\description{
This function is the workhorse behind \code{\link{bcct}} and \code{\link{bcctu}}.
}
\usage{
bcct.fit(priornum, maximal.mod, IP, eta.hat, ini.index, ini.beta, ini.sig, 
iters, save, name, null.move.prob, a, b, progress)
}

\arguments{
  \item{priornum}{
A numeric scalar indicating which prior is to be used: 1 = \code{"UIP"}, 2 = \code{"SBH"}.
}
  \item{maximal.mod}{
An object of class \code{"glm"} giving the fit of the maximal model.
}
  \item{IP}{
A p by p matrix giving the inverse of the prior scale matrix for the maximal model.
}
  \item{eta.hat}{
A vector of length n (number of cells) giving the posterior mode of the linear predictor under 
the maximal model.
}
  \item{ini.index}{
A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the initial model.
}
  \item{ini.beta}{
A numeric vector giving the starting values of the log-linear parameters for the MCMC algorithm.
}
  \item{ini.sig}{
A numeric scalar giving the starting value of sigma^2 for the MCMC algorithm.
}
  \item{iters}{
The number of iterations of the MCMC algorithm to peform.
}
  \item{save}{
If positive, the function will save the MCMC output to external text 
files every \code{save} iterations. If zero , the function will not save 
the MCMC output to external files. 
}
  \item{name}{
A prefix to the external files saved if the argument \code{save} is 
positive. If \code{NULL}, then the external files will have no prefix. 
}
  \item{null.move.prob}{
A scalar argument giving the probability of performing a null move, i.e. proposing a move to the 
current model.
}
  \item{a}{
The shape hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014).
}
  \item{b}{
The scale hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014).
}
  \item{progress}{
Logical argument. If \code{TRUE}, then a progress bar will be displayed.
}

}

\value{
The function will return a list with the following components:

\item{BETA}{An \code{iters} by p matrix containing the sampled values of the log-linear parameters, where p 
is the number of log-linear parameters in the maximal model. For elements of this matrix which correspond 
to a log-linear parameter which is not present for the current model a zero is returned.}

\item{MODEL}{A vector of length \code{iters} giving the sampled model indicators in hexadecimal form.}

\item{SIG}{A vector of length \code{iters} giving the sampled values for sigma^2 under the Sabanes-Bove & Held prior. If the 
unit information prior is used then the components of this vector will be one.}

\item{rj_acc}{A binary vector of the same length as the number of reversible jump moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{mh_acc}{A binary vector of the same length as the number of Metropolis-Hastings moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

}

\references{
Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Journal of Statistical Software}, \bold{58 (7)}, 
1--27. \url{http://www.jstatsoft.org/v58/i07/}
}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
This function will not typically be called by the user.
}

\seealso{
\code{\link{bcct}},
\code{\link{bcctu}}.
}

\examples{
data(AOH)
## Load the AOH data.

maximal.mod<-glm(formula=y~(alc+hyp+obe)^3,data=AOH,x=TRUE,y=TRUE,
contrasts=list(alc="contr.sum",hyp="contr.sum",obe="contr.sum"))
## Set up the maximal model which in this case is the saturated 
## model.

curr.index<-formula2index(big.X=maximal.mod$x,formula=y~alc+hyp+obe+hyp:obe,data=AOH)
## Set up the binary vector for the model containing all main effects and the 
## hyp:obe interaction.

IP<-t(maximal.mod$x)\%*\%maximal.mod$x/length(maximal.mod$y)
IP[,1]<-0
IP[1,]<-0
## Set up the inverse scale matrix for the prior distribution under
## the maximal model. 

bmod<-beta_mode(X=maximal.mod$x,prior="UIP",y=maximal.mod$y,IP=IP)
## Find the posterior mode under the maximal model

eta.hat<-as.vector(maximal.mod$x\%*\%bmod)
## Find the posterior mode of the linear predictor
## under the maximal model.

set.seed(1)
## Set seed for reproducibility

test1<-bcct.fit(priornum=1, maximal.mod=maximal.mod, IP=IP, eta.hat=eta.hat, 
ini.index=curr.index, ini.beta=bmod[curr.index==1], ini.sig=1, iters=5, save=0, 
name=NULL,null.move.prob=0.5, a=0.001, b=0.001, progress=TRUE) 
## Run for 5 iterations starting at model defined by curr.index.

test1$MODEL
## Look at sampled model indicators. Should be:
## [1] "fe00c0" "fe0000" "fe0000" "fe0000" "fe0000"

model2index(test1$MODEL,dig=24)
## Convert these to binary indicators of the log-linear parameters.
## Will get: 

#       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11] [,12] [,13]
#fe00c0    1    1    1    1    1    1    1    0    0     0     0     0     0
#fe0000    1    1    1    1    1    1    1    0    0     0     0     0     0
#fe0000    1    1    1    1    1    1    1    0    0     0     0     0     0
#fe0000    1    1    1    1    1    1    1    0    0     0     0     0     0
#fe0000    1    1    1    1    1    1    1    0    0     0     0     0     0
#       [,14] [,15] [,16] [,17] [,18] [,19] [,20] [,21] [,22] [,23] [,24]
#fe00c0     0     0     0     1     1     0     0     0     0     0     0
#fe0000     0     0     0     0     0     0     0     0     0     0     0
#fe0000     0     0     0     0     0     0     0     0     0     0     0
#fe0000     0     0     0     0     0     0     0     0     0     0     0
#fe0000     0     0     0     0     0     0     0     0     0     0     0

## See how the hyp:obe interactions in columns 17 and 18 gets dropped after 
## the 1st iteration. 

}


