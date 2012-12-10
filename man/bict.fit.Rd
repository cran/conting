\name{bict.fit}
\alias{bict.fit}

\title{
Bayesian Analysis of Incomplete Contingency Tables
}
\description{
This function is the workhorse behind \code{\link{bict}} and \code{\link{bictu}}.
}
\usage{
bict.fit(priornum, missing1, missing2, maximal.mod, IP, eta.hat, ini.index, 
ini.beta, ini.sig, ini.y0, iters, save, name, null.move.prob)
}

\arguments{
  \item{priornum}{
A numeric scalar indicating which prior is to be used: 1=\code{"UIP"}, 2=\code{"SBH"}.
}
  \item{missing1}{
A vector of the same length as the number of missing cell counts giving the row numbers of 
the data.frame in \code{data} which correspond to the missing cell counts.
}
  \item{missing2}{
A vector of the same length as the number of censored cell counts giving the row numbers of 
the data.frame in \code{data} which correspond to the censored cell counts.
}
  \item{maximal.mod}{
An object of class "\code{glm}" giving the fit of the maximal model.
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
  \item{ini.y0}{
A numeric vector giving the starting values of the missing and censored cell entries for the MCMC algorithm.
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
positive. If NULL, then the external files will have no prefix. 
}
  \item{null.move.prob}{
A scalar argument giving the probability of performing a null move, i.e. proposing a move to the 
current model.
}
}

\value{
The function will return a list with the following components.

\item{BETA}{An \code{iters} by p matrix containing the sampled values of the log-linear parameters, where p 
is the number of log-linear parameters in the maximal model. For elements of this matrix which correspond 
to a log-linear parameter which is not present for the current model a zero is returned.}

\item{MODEL}{A vector of length \code{iters} giving the sampled model indicators in hexadecimal format.}

\item{SIG}{A vector of length \code{iters} giving the sampled values for sigma^2 under the Sabanes-Bove & Held prior. 
If the unit information prior is used then the components of this vector will be one.}

\item{Y0}{An \code{iters} by k matrix giving the sampled values of the missing and censored cell counts, where k is 
the total number of missing and censored cell counts.}

\item{rj_acc}{A binary vector of the same length as the number of reversible jump moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{mh_acc}{A binary vector of the same length as the number of Metropolis-Hastings moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

}

\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
This function will not typically be called by the user.
}

\seealso{
\code{\link{bict}},
\code{\link{bictu}}.
}
\examples{
data(rabbit)
## Load rabbit data.

rabbit$z<-rabbit$y
rabbit$z[is.na(rabbit$y)]<-0
## Define a new variable in rabbit data.frame which is equal to y except where
## y is NA, in which case z=0. This is just so we can fit maximal model to the
## complete contingency table.

maximal.mod<-glm(formula=z~(first+second+hunter)^2,data=rabbit,x=TRUE,y=TRUE,
contrasts=list(first="contr.sum",second="contr.sum",hunter="contr.sum"))
## Fit maximal model to complete contingency table.

curr.index<-formula2index(big.X=maximal.mod$x,formula=z~first+second+hunter,data=rabbit)
## Set up binary vector for independence model.

IP<-t(maximal.mod$x)\%*\%maximal.mod$x/length(maximal.mod$y)
IP[,1]<-0
IP[1,0]<-0
## Set up the inverse scale matrix for the prior distribution under
## the maximal model. 

bmod<-beta_mode(X=maximal.mod$x[!is.na(rabbit$y),],prior="UIP",
y=maximal.mod$y[!is.na(rabbit$y)],IP=IP)
## Find the posterior mode under the maximal model fitted to observed cell counts.
eta.hat<-as.vector(maximal.mod$x\%*\%bmod)
## Find the posterior mode of the linear predictor
## under the maximal model.

set.seed(1)
## Set seed for reproducibility
test1<-bict.fit(priornum=1, missing1=(1:length(maximal.mod$y))[is.na(rabbit$y)],
missing2=NULL,maximal.mod=maximal.mod, IP=IP, eta.hat=eta.hat, ini.index=curr.index, 
ini.beta=bmod[curr.index==1], ini.sig=1, ini.y0=222,iters=5, save=0, name=NULL,
null.move.prob=0.5) 
## Run for 5 iterations starting at model defined by curr.index.
test1$MODEL
## Look at sampled model indicators. Should be:
## [1] "78" "78" "78" "78" "78"

}