\name{beta_mode}
\alias{beta_mode}

\title{
Posterior Mode
}

\description{
This function finds the posterior mode of the log-linear parameters of a 
log-linear model with a given design matrix and prior distribution.
}
\usage{
beta_mode(X, prior = "SBH", y, IP , a = 0.001 , b = 0.001)
}

\arguments{
  \item{X}{
The n by p design matrix where n is the number of cells and p is the number 
of log-linear parameters.
}
  \item{prior}{
The prior distribution. It can be one of \code{c("UIP","SBH")}, where \code{"UIP"} = 
unit information prior; and \code{"SBH"} = Sabanes-Bove & Held prior. The default 
value is \code{"SBH"}.
}
  \item{y}{
The n by 1 vector of cell counts.
}
  \item{IP}{
A p by p matrix giving the inverse of the prior scale matrix.
}
  \item{a}{
The shape hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014).
}
  \item{b}{
The scale hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014).
}

}
\details{
The posterior mode is found by maximising the log unnormalised posterior pdf given by the sum of the 
log-likelihood and the log of the prior pdf. This optimisation is achieved using a quasi Newton-Raphson
method. 

For details on the unit information and Sabanes-Bove & Held priors for generalised linear models see Ntzoufras 
et al (2003) and Sabanes-Bove & Held (2011), respectively. See Overstall & King (2014), and the references therein, 
for their application to log-linear models and contingency tables.

The posterior mode is required for the reversible jump algorithm 
implemented from Forster et al (2012).
}
\value{
\code{beta_mode} will return a p by 1 vector containing the posterior mode of the log
linear parameters.}

\references{
Sabanes-Bove, D. & Held, L. (2011) Hyper-g priors for generalized linear models. \emph{Bayesian
Analysis}, \bold{6 (3)}, 387--410.

Forster, J.J., Gill, R.C. & Overstall, A.M. (2012) Reversible jump methods for 
generalised linear models and generalised linear mixed models. \emph{Statistics and Computing}, 
\bold{22 (1)}, 107--120.

Nztoufras, I., Dellaportas, P. & Forster, J.J. (2003) Bayesian variable and link 
determination for generalised linear models. \emph{Journal of Statistical Planning and 
Inference}, \bold{111 (1)}, 165--180.

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

\examples{
data(AOH) ## loads the AOH data

X<-model.matrix(~alc+hyp+obe,data=AOH,
contrasts=list(alc="contr.sum",hyp="contr.sum",obe="contr.sum"))
## Sets up the design matrix for the independence model

IP<-(t(X)\%*\%X)/dim(X)[1]
## Set up inverse of prior scale matrix

beta_mode(X=X,prior="UIP",y=AOH$y,IP=IP)
## Finds the posterior mode of the log-linear parameters under the
## independence model with the unit information prior. Will get:
#X(Intercept)        Xalc1        Xalc2        Xalc3        Xhyp1        Xobe1 
# 2.894270420 -0.045859743 -0.071775824  0.089541068 -0.504141954  0.008163604 
#       Xobe2 
#-0.016327209

beta_mode(X=X,prior="SBH",y=AOH$y,IP=IP)
## Finds the posterior mode of the log-linear parameters under the
## independence model with the Sabanes-Bove & Held prior. Will get:
#X(Intercept)        Xalc1        Xalc2        Xalc3        Xhyp1        Xobe1 
# 2.908298763 -0.043704371 -0.068212247  0.085338704 -0.473628107  0.007762839 
#       Xobe2 
#-0.015525678


}