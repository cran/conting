\name{iwls_mh}
\alias{iwls_mh}

\title{
Iterated Weighted Least Square Metropolis Hastings Algorithm
}
\description{
This function implements one iteration of the Iterated Weight Least Square Metropolis 
Hastings Algorithm as proposed by Gamerman (1997) for generalised linear models as applied 
to log-linear models.
}
\usage{
iwls_mh(curr.y, curr.X, curr.beta, iprior.var)
}

\arguments{
  \item{curr.y}{
A vector of length n giving the cell counts.
}
  \item{curr.X}{
An n by p design matrix for the current model, where p is the number of 
log-linear parameters.
}
  \item{curr.beta}{
A vector of length p giving the current log-linear parameters.
}
  \item{iprior.var}{
A p by p matrix giving the inverse of the prior variance matrix.
}
}
\details{
For details of the original algorithm see Gamerman (1997). For its application to 
log-linear models see Overstall & King (2012), and the references therein.
}
\value{
The function will output a vector of length p giving the new values of the log-linear 
parameters.
}
\references{
Gamerman, D. (1997) Sampling from the posterior distribution in generalised linear mixed
models. \emph{Statistics and Computing}, \bold{7}, 57--68.

Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
This function will not typically be called by the user.
}
\examples{
set.seed(1)
## Set seed for reproducibility
data(AOH)
## Load AOH data

maximal.mod<-glm(y~alc+hyp+obe,family=poisson,x=TRUE,contrasts=list(alc="contr.sum",
hyp="contr.sum",obe="contr.sum"),data=AOH)
## Fit independence model to get a design matrix

IP<-t(maximal.mod$x)\%*\%maximal.mod$x/length(AOH$y)
IP[,1]<-0
IP[1,]<-0
## Set up inverse prior variance matrix under the UIP

## Let the current parameters be the MLE under the independence model
as.vector(coef(maximal.mod))
#[1]  2.89365105 -0.04594959 -0.07192507  0.08971628 -0.50545335  0.00818037
#[7] -0.01636074

## Update parameters using MH algorithm
iwls_mh(curr.y=AOH$y,curr.X=maximal.mod$x,curr.beta=coef(maximal.mod),iprior.var=IP)

## Will get:
#[1]  2.86468919 -0.04218623 -0.16376055  0.21656167 -0.49528676 -0.05026597
#[7]  0.02726671
}
