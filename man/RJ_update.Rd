\name{RJ_update}
\alias{RJ_update}

\title{
Reversible Jump Algorithm
}
\description{
This function implements one iteration of the orthogonal projection reversible jump algorithm for generalised 
linear models proposed by Forster et al (2012) applied to log-linear models.
}
\usage{
RJ_update(prop.index, curr.index, curr.beta, eta.hat, curr.y, big.X, 
proposal.probs, i.prop.prior.var, i.curr.prior.var)
}

\arguments{
  \item{prop.index}{
A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the proposed model.
}
  \item{curr.index}{
A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the current model.
}
  \item{curr.beta}{
A vector of length \code{sum(curr.index)} giving the log-linear parameters under the current model.
}
  \item{eta.hat}{
A vector of length n (number of cells) giving the posterior mode of the linear predictor under 
the maximal model.
}
  \item{curr.y}{
A vector of length n giving the cell counts.
}
  \item{big.X}{
The design matrix under the maximal model.
}
  \item{proposal.probs}{
A numeric vector of length 2. The first element gives the probability of proposing a move from the proposed model to the 
current model. The second element gives the probability of proposing a move from the current model to the 
proposed model.
}
  \item{i.prop.prior.var}{
A matrix giving the inverse of the prior variance matrix for the log-linear parameters under the proposed model.
}
  \item{i.curr.prior.var}{
A matrix giving the inverse of the prior variance matrix for the log-linear parameters under the current model.
}
}
\details{
For the original algorithm see Forster et al (2012). For details on its application to log-linear models 
see Overstall & King (2013), and the references therein.
}
\value{
The function will return a list with the following components:
\item{new.beta}{A vector giving the new log-linear parameters.}
\item{new.index}{A binary vector indicating which log-linear parameters are present in the new model.}
}
\references{
Forster, J.J., Gill, R.C. & Overstall, A.M. (2012) Reversible jump methods for generalised linear models and generalised 
linear mixed models. \emph{Statistics and Computing}, \bold{22}, 107--120.

Overstall, A.M. & King, R. (2013) conting: An R package for Bayesian analysis of
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
set.seed(4)
## Set seed for reproducibility
data(AOH)
## Load data
maximal.mod<-glm(y~(alc+hyp+obe)^3,family=poisson,x=TRUE,contrasts=list(alc="contr.sum",
hyp="contr.sum",obe="contr.sum"),data=AOH)
## Fit maximal model to get a design matrix

IP<-t(maximal.mod$x)\%*\%maximal.mod$x/length(AOH$y)
IP[,1]<-0
IP[1,]<-0
## Calculate inverse prior scale matrix under maximal model. Under the UIP this
## is the inverse prior variance matrix. Under the SBH prior, we need to divide
## this matrix by the current value of SIG.

bmod<-beta_mode(X=maximal.mod$x,y=AOH$y,IP=IP)
## Find posterior mode under maximal model with UIP
eta.hat<-as.vector(maximal.mod$x\%*\%bmod)
## Find posterior mode of linear predictor.

curr.index<-formula2index(big.X=maximal.mod$x,formula=y~alc+hyp+obe+alc:hyp,data=AOH)
## Calculate index for current model including alc:hyp interaction
curr.index
## Print out current index
#[1] 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0

pm<-prop_mod(curr.index=curr.index,data=AOH,maximal.mod=maximal.mod)
## Propose a model
p2<-(1-pm$null.move.prob)/pm$total.choices
p2
## Calculate probability of proposing proposed model from current model
#[1] 0.1666667

prop.index<-pm$new.index
prop.index
## Assign and print out proposal index
# [1] 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

dm<-prop_mod(curr.index=prop.index,data=AOH,maximal.mod=maximal.mod,null.move.prob=0)
p1<-(1-pm$null.move.prob)/dm$total.choices
p1
## Calculate probability of proposing current model from proposed model
#[1] 0.1666667

RJ_update(prop.index=prop.index,curr.index=curr.index,
curr.beta=coef(maximal.mod)[curr.index==1],eta.hat=eta.hat,curr.y=AOH$y,big.X=maximal.mod$x,
proposal.probs=c(p1,p2),
i.prop.prior.var=IP[prop.index==1,prop.index==1],
i.curr.prior.var=IP[curr.index==1,curr.index==1])

## Do one iteration of reversible jump algorithm. Will get:

#$new.beta
#(Intercept)        alc1        alc2        alc3        hyp1        obe1 
# 2.87128918 -0.07098006 -0.07221330  0.08748803 -0.51899802 -0.07855115 
#       obe2 
#-0.02474727 
#
#$new.index
# [1] 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0

}
