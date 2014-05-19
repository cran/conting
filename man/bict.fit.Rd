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
ini.beta, ini.sig, ini.y0, iters, save, name, null.move.prob, a, b, progress)
}

\arguments{
  \item{priornum}{
A numeric scalar indicating which prior is to be used: 1 = \code{"UIP"}, 2 = \code{"SBH"}.
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

\references{
Overstall, A.M. & King, R. (2014) conting: An R package for the Bayesian analysis of 
complete and incomplete contingency tables. \emph{Technical Report}, University of St Andrews.
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
data(spina)
## Load spina data.

spina$z<-spina$y
spina$z[is.na(spina$y)]<-0
## Define a new variable in spina data.frame which is equal to y except where
## y is NA, in which case z=0. This is just so we can fit maximal model to the
## complete contingency table.

maximal.mod<-glm(formula=z~(S1+S2+S3+eth)^2,data=spina,x=TRUE,y=TRUE,
contrasts=list(S1="contr.sum",S2="contr.sum",S3="contr.sum",
eth="contr.sum"))
## Fit maximal model to complete contingency table.

curr.index<-formula2index(big.X=maximal.mod$x,formula=z~S1+S2+S3+eth,data=spina)
## Set up binary vector for independence model.

IP<-t(maximal.mod$x)\%*\%maximal.mod$x/length(maximal.mod$y)
IP[,1]<-0
IP[1,]<-0
## Set up the inverse scale matrix for the prior distribution under
## the maximal model. 

bmod<-beta_mode(X=maximal.mod$x[!is.na(spina$y),],prior="UIP",
y=maximal.mod$y[!is.na(spina$y)],IP=IP)
## Find the posterior mode under the maximal model fitted to observed cell counts.
eta.hat<-as.vector(maximal.mod$x\%*\%bmod)
## Find the posterior mode of the linear predictor
## under the maximal model.

set.seed(1)
## Set seed for reproducibility
test1<-bict.fit(priornum=1, missing1=(1:length(maximal.mod$y))[is.na(spina$y)],
missing2=NULL,maximal.mod=maximal.mod, IP=IP, eta.hat=eta.hat, ini.index=curr.index, 
ini.beta=bmod[curr.index==1], ini.sig=1, ini.y0=c(500,200,20),iters=10, save=0, name=NULL,
null.move.prob=0.5, a=0.001, b=0.001, progress = FALSE) 
## Run for 10 iterations starting at model defined by curr.index.
test1$MODEL
## Look at sampled model indicators. Should be:
#  [1] "7e00" "7e00" "7e00" "7e00" "7e00" "7e00" "7e00" "7e00" "7f00" "7f00"

model2index(test1$MODEL,dig=15)
## Convert these to binary indicators of the log-linear parameters.
## Will get: 

#     [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10] [,11]
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7e00    1    1    1    1    1    1    0    0    0     0     0
#7f00    1    1    1    1    1    1    1    0    0     0     0
#7f00    1    1    1    1    1    1    1    0    0     0     0
#     [,12] [,13] [,14] [,15]
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7e00     0     0     0     0
#7f00     0     0     0     0
#7f00     0     0     0     0

}