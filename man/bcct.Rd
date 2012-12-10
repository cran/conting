\name{bcct}
\alias{bcct}
\alias{bcctu}
\title{
Bayesian Analysis of Complete Contingency Tables
}
\description{
These functions implement a Bayesian analysis of complete contingency tables. This is 
accomplished using an MCMC algorithm where the null moves are performed using 
a Metropolis-Hastings algorithm and the between models moves are performed using a 
reversible jump algorithm.

\code{bcct} should be used initially, and \code{bcctu} should be used to do additional 
MCMC iterations, if required.
}
\usage{
bcct(formula, data, n.sample, prior = "UIP", start.formula = NULL, 
start.beta = NULL, start.sig = NULL, save = 0, name = NULL, null.move.prob=0.5)

bcctu(object, n.sample, save = NULL, name = NULL)
}

\arguments{
  
  \item{formula}{An object of class "\code{formula}": a symbolic description of the 
maximal model.}

  \item{object}{An object of class "\code{bcct}" produced as a previous call to \code{bcct} or \code{bcctu}.}

  \item{data}{An object of class "\code{data.frame}" containing the variables in the model.}

  \item{n.sample}{
A numeric scalar giving the number of (additional, in the case of \code{bcctu}) MCMC iterations to perform.
}
  \item{prior}{An optional argument giving the prior to be used in the analysis. It can be one of 
\code{c("UIP","SBH")}, where \code{"UIP"}=unit information prior; and \code{"SBH"}=Sabanes-Bove & Held prior. The 
default value is \code{"UIP"}.
}

  \item{start.formula}{
An optional argument giving an object of class "\code{formula}": a symbolic description of the 
starting model in the MCMC algorithm. If NULL (the default) the starting model will be the maximal 
model. 
}
  \item{start.beta}{
An optional argument giving the starting values of the log-linear parameters for the MCMC algorithm. 
It should be a vector of the same length as the number of log-linear parameters in the staring model 
implied by the argument \code{start.formula}. If NULL (the default) the starting value will be the 
posterior mode under the maximal model.
}
  \item{start.sig}{
An optional argument giving the starting value of sigma^2 for the MCMC algorithm when the argument 
prior is "SBH". If NULL (the default) the starting value will be one.
}
  \item{save}{
An optional argument for saving the MCMC output mid-algorithm. 

For \code{bcct}, if positive, the function will save the MCMC output to external text 
files every \code{save} iterations. If zero (the default), the function will not save 
the MCMC output to external files. 

For \code{bcctu}, if non-NULL, the function will save the MCMC output to external text 
files every \code{save} iterations. If NULL (the default), it will inherit the value of 
\code{save} from the previous call to \code{bcct} or \code{bcctu}.
}
  \item{name}{
An optional argument giving a prefix to the file name of the external files saved if the argument \code{save} is 
positive. For \code{bcct}, a value of NULL means the external files will not have a prefix. For 
\code{bcctu}, a value of NULL, means the prefix will be inherited from the previous call to 
\code{bcct} or \code{bcctu}. 
}
  \item{null.move.prob}{
An optional scalar argument giving the probability of performing a null move in the reversible jump 
algorithm, i.e. proposing a move to the current model. The default value is 0.5.
}
}
\details{
For identifiability, the parameters are constrained. The \code{\link{conting-package}} uses sum-to-zero constraints. 
See Overstall & King (2012), and the references therein, for more details. 

The Metropolis-Hastings algorithm employed is the iterated weighted least squares method for 
generalised linear models (GLMs) proposed by Gamerman (1997). The reversible jump algorithm employed is 
that orthogonal projections method for GLMs proposed by Forster et al (2012). For details on these 
methods applied to log-linear models see Overstall & King (2012), and the references therein. 

For details on the unit information and Sabanes-Bove & Held priors for generalised linear models see Ntzoufras 
et al (2003) and Sabanes-Bove & Held (2011), respectively. See Overstall & King (2012), and the references therein, 
for their application to log-linear models and contingency tables.

}
\value{
The functions will return an object of class "\code{bcct}" which is a list with the following components:

\item{BETA}{An \code{n.sample} by p matrix containing the sampled values of the log-linear parameters, where p 
is the number of log-linear parameters in the maximal model. For elements of this matrix which correspond 
to a log-linear parameter which is not present for the current model a zero is returned.}

\item{MODEL}{A vector of length \code{n.sample} giving the sampled model indicators in hexadecimal format.}

\item{SIG}{A vector of length \code{n.sample} giving the sampled values for sigma^2 under the Sabanes-Bove & Held prior. If the 
unit information prior is used then the components of this vector will be one.}

\item{rj_acc}{A binary vector of the same length as the number of reversible jump moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{mh_acc}{A binary vector of the same length as the number of Metropolis-Hastings moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.}

\item{priornum}{A numeric scalar indicating which prior was used: 1=\code{"UIP"}, 2=\code{"SBH"}.}

\item{maximal.mod}{An object of class "\code{glm}" giving the fit of the maximal model.}

\item{IP}{A p by p matrix giving the inverse of the prior scale matrix for the maximal model.}

\item{eta.hat}{A vector of length n (number of cells) giving the posterior mode of the linear predictor under 
the maximal model.}

\item{save}{The argument \code{save}.}

\item{name}{The argument \code{name}.}

\item{null.move.prob}{
The argument \code{null.move.prob}.}

\item{time}{
The total computer time (in seconds) used for the MCMC computations.}
}

\references{
Sabanes-Bove, D. & Held, L. (2011) Hyper-g priors for generalized linear models. \emph{Bayesian
Analysis}, \bold{6}, 387--410.

Forster, J.J., Gill, R.C. & Overstall, A.M. (2012) Reversible jump methods for generalised linear models and generalised 
linear mixed models. \emph{Statistics and Computing}, \bold{22}, 107--120.

Gamerman, D. (1997) Sampling from the posterior distribution in generalised linear mixed
models. \emph{Statistics and Computing}, \bold{7}, 57--68.

Nztoufras, I., Dellaportas, P. & Forster, J.J. (2003) Bayesian variable and link 
determination for generalised linear models. \emph{Journal of Statistical Planning and 
Inference}, \bold{111}, 165--180.

Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
These functions are wrappers for \code{\link{bcct.fit}}.
}

\seealso{
\code{\link{bcct.fit}},
\code{\link{AOH}},
\code{\link{heart}}.

}
\examples{
set.seed(1)
## Set seed for reproducibility.

data(AOH)
## Load the AOH data

test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=50)
## Let the maximal model be the saturated model. Starting from the 
## posterior mode of the maximal model do 500 iterations under the unit 
## information prior.

test1<-bcctu(object=test1,n.sample=50)
## Do another 500 iterations

test1
## Printing out a bcct object produces this simple summary

#Number of cells in table = 24 
#
#Maximal model =
#y ~ (alc + hyp + obe)^3
#
#Number of log-linear parameters in maximal model = 24 
#
#Number of MCMC iterations = 100 
#
#Computer time for MCMC = 00:00:01 
#
#Prior distribution for log-linear parameters = UIP

summary(test1)
## Printing out a summary produces a bit more:

#Posterior summary statistics of log-linear parameters:
#            post_prob post_mean post_var lower_lim upper_lim
#(Intercept)         1  2.877924 0.002574   2.78778   2.97185
#alc1                1 -0.060274 0.008845  -0.27772   0.06655
#alc2                1 -0.049450 0.006940  -0.20157   0.11786
#alc3                1  0.073111 0.005673  -0.05929   0.20185
#hyp1                1 -0.544988 0.003485  -0.65004  -0.42620
#obe1                1 -0.054672 0.007812  -0.19623   0.12031
#obe2                1  0.007809 0.004127  -0.11024   0.11783
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively
#
#Posterior model probabilities:
#  prob model_formula                                 
#1 0.45 ~alc + hyp + obe                              
#2 0.30 ~alc + hyp + obe + hyp:obe                    
#3 0.11 ~alc + hyp + obe + alc:hyp + hyp:obe          
#4 0.06 ~alc + hyp + obe + alc:hyp + alc:obe + hyp:obe
#5 0.05 ~alc + hyp + obe + alc:hyp                    
#
#Total number of models visited =  7 
#
#Under the X2 statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  11.79   20.16   23.98   24.70   28.77   52.40 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   8.18   24.22   31.51   30.12   35.63   42.49 
#
#Bayesian p-value =  0.28

## For more examples see Overstall & King (2012).


}
