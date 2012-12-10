\name{bict}
\alias{bict}
\alias{bictu}

\title{
Bayesian Analysis of Incomplete Contingency Tables
}
\description{
These functions implement a Bayesian analysis of incomplete contingency tables. This is 
accomplished using a data augmentation MCMC algorithm where the null moves 
are performed using the Metropolis-Hastings algorithm and the between models moves are performed using the 
reversible jump algorithm. This function can also accomodate cases where one of the sources observes 
a mixture of individuals from target and non-target populations. This results in the some of the cell 
counts being censored.

\code{bict} should be used initially, and \code{bictu} should be used to do additional 
MCMC iterations, if needed.
}
\usage{
bict(formula, data, n.sample, prior = "UIP", cens = NULL, start.formula = NULL, 
start.beta = NULL, start.sig = NULL, start.y0 = NULL, save = 0, name = NULL, 
null.move.prob=0.5)

bictu(object, n.sample, save = NULL, name = NULL)
}

\arguments{
  \item{formula}{
An object of class "\code{formula}": a symbolic description of the 
maximal model.}

  \item{object}{An object of class \code{bict} produced as a previous call to \code{bict} or \code{bictu}.}

  \item{data}{
An object of class "\code{data.frame}" containing the variables in the model.
}
  \item{n.sample}{
A numeric scalar giving the number of MCMC iterations to peform.
}
  \item{prior}{
An optional argument giving the prior to be used in the analysis. It can be one of 
\code{c("UIP","SBH")}, where \code{"UIP"}=unit information prior; and \code{"SBH"}=Sabanes-Bove & Held 
prior. The default value is \code{"UIP"}.
}
  \item{cens}{
A numeric vector indicating the row numbers of the data.frame in \code{data} which correspond to the 
censored cells. This can be found using the function \code{\link{find_cens}}.
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
prior is \code{"SBH"}. If NULL (the default) the starting value will be one.
}
  \item{start.y0}{
An optional argument giving the starting values of the missing and censored cell counts. This 
should have the same length as the number of missing and censored cell counts. 
}
  \item{save}{
An optional argument for saving the MCMC output mid-algorithm. 

For \code{bict}, if positive, the function will save the MCMC output to external text 
files every \code{save} iterations. If zero (the default), the function will not save 
the MCMC output to external files. 

For \code{bictu}, if non-NULL, the function will save the MCMC output to external text 
files every \code{save} iterations. If NULL (the default), it will inherit the value of 
\code{save} from the previous call to \code{bict} or \code{bictu}.
}
  \item{name}{
An optional argument giving a prefix to the external files saved if the argument \code{save} is 
positive. For \code{bict}, a value of NULL means the external files will not have a prefix. For 
\code{bictu}, a value of NULL, means the prefix will be inherited from the previous call to 
\code{bict} or \code{bictu}. 
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
the orthogonal projections method for GLMs proposed by Forster et al (2012). For details on these 
methods applied to log-linear models through the data-augmentation algorithm see Overstall & King (2012), 
and the references therein. For details on the censored approach see Overstall et al (2012).

For details on the unit information and Sabanes-Bove & Held priors for generalised linear models see Ntzoufras 
et al (2003) and Sabanes-Bove & Held (2011), respectively. See Overstall & King (2012), and the references therein, 
for their application to log-linear models and contingency tables.
}
\value{
The functions will return an object of class "\code{bict}" which is a list with the following components.

\item{BETA}{
An \code{n.sample} by p matrix containing the sampled values of the log-linear parameters, where p 
is the number of log-linear parameters in the maximal model. For elements of this matrix which correspond 
to a log-linear parameter which is not present for the current model a zero is returned.
}

\item{MODEL}{
A vector of length \code{n.sample} giving the sampled model indicators in hexadecimal format.
}

\item{SIG}{
A vector of length \code{n.sample} giving the sampled values for sigma^2 under the Sabanes-Bove & Held prior. 
If the unit information prior is used then the components of this vector will be one.
}

\item{Y0}{
An \code{n.sample} by k matrix giving the sampled values of the missing and censored cell counts, where k is 
the total number of missing and censored cell counts.
}

\item{missing1}{
A vector of the same length as the number of missing cell counts giving the row numbers of 
the data.frame in \code{data} which correspond to the missing cell counts.
}

\item{missing2}{
A vector of the same length as the number of censored cell counts giving the row numbers of 
the data.frame in \code{data} which correspond to the censored cell counts.
}

\item{missing_details}{
The rows of the data.frame in \code{data} corresponding to the missing cell counts.
}

\item{censored_details}{
The rows of the data.frame in \code{data} corresponding to the censored cell counts.
}

\item{rj_acc}{
A binary vector of the same length as the number of reversible jump moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.
}

\item{mh_acc}{
A binary vector of the same length as the number of Metropolis-Hastings moves attempted. A 0 indicates 
that the proposal was rejected, and a 1 that the proposal was accepted.
}

\item{priornum}{
A numeric scalar indicating which prior was used: 1=\code{"UIP"}, 2=\code{"SBH"}.
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

\item{save}{
The argument \code{save}.
}

\item{name}{
The argument \code{name}.
}
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

Overstall, A.M., King, R., Bird, S.M., Hutchinson, S.J. & Hay, G. (2012) Estimating the 
Number of People Who Inject Drugs in Scotland using Multi-List Data with Left Censoring. 
\emph{Technical Report}, University of St Andrews.

}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
These functions are wrappers for \code{\link{bict.fit}}.
}

\seealso{
\code{\link{bict.fit}},
\code{\link{rabbit}},
\code{\link{ScotPWID}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility.

data(rabbit)
## Load the rabbit data

test1<-bict(formula=y~(first+second+hunter)^2,data=rabbit,n.sample=50,prior="SBH")
## Let the maximal model be the model with two-way interactions. Starting from the 
## posterior mode of the maximal model do 50 iterations under the Sabanes-Bove & Held 
## prior.

test1<-bictu(object=test1,n.sample=50)
## Do another 50 iterations

test1
## Printing out a bict object produces this simple summary

#Number of cells in table = 8 
#
#Maximal model =
#y ~ (first + second + hunter)^2
#
#Number of log-linear parameters in maximal model = 7 
#
#Number of MCMC iterations = 100 
#
#Computer time for MCMC = 00:00:00 
#
#Prior distribution for log-linear parameters = SBH 
#
#Number of missing cells = 1 
#
#Number of censored cells = 0

summary(test1)
## Printing out a summary produces a bit more:

#Posterior summary statistics of log-linear parameters:
#                post_prob post_mean post_var lower_lim upper_lim
#(Intercept)          1.00   2.65380 0.021194    2.4632   2.90036
#first1               1.00  -0.12840 0.011542   -0.2899   0.09038
#second1              1.00  -0.07974 0.009306   -0.2552   0.04317
#hunter1              1.00  -0.36272 0.009595   -0.4925  -0.18710
#first1:second1       0.79   0.23423 0.008526    0.1137   0.39316
#second1:hunter1      0.86  -0.29970 0.011300   -0.5783  -0.13356
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively
#
#Posterior model probabilities:
#  prob model_formula                                                         
#1 0.44 ~first + second + hunter + first:second + first:hunter + second:hunter
#2 0.21 ~first + second + hunter + first:hunter + second:hunter               
#3 0.21 ~first + second + hunter + first:second + second:hunter               
#4 0.10 ~first + second + hunter + first:second                               
#
#Total number of models visited =  5 
#
#Posterior mean of total population size = 137.96 
#95 % highest posterior density interval for total population size = ( 122 156 ) 
#
#Under the X2 statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.435   4.108   6.774   7.320   8.955  26.000 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  5.836   9.616  13.580  14.150  19.430  25.950 
#
#Bayesian p-value =  0.12

}
