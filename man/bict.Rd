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
bict(formula, data, n.sample, prior = "SBH", cens = NULL, start.formula = NULL, 
start.beta = NULL, start.sig = NULL, start.y0 = NULL, save = 0, name = NULL, 
null.move.prob=0.5, a = 0.001, b = 0.001, progress = FALSE)

bictu(object, n.sample, save = NULL, name = NULL, progress = FALSE)
}

\arguments{
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the 
maximal model.}

  \item{object}{An object of class \code{"bict"} produced as a previous call to \code{bict} or \code{bictu}.}

  \item{data}{An object of class \code{"data.frame"} (or \code{"table"}) containing the variables in the 
model. If the model variables are not found in \code{data}, the variables are taken from 
\code{environment(formula)}, typically the environment from which \code{bict} is called.}

  \item{n.sample}{
A numeric scalar giving the number of MCMC iterations to peform.
}
  \item{prior}{
An optional argument giving the prior to be used in the analysis. It can be one of 
\code{c("UIP","SBH")}, where \code{"UIP"} = unit information prior; and \code{"SBH"} = Sabanes-Bove & Held 
prior. The default value is \code{"SBH"}.
}
  \item{cens}{
A numeric vector indicating the row numbers of the data.frame in \code{data} which correspond to the 
censored cells. This can be found using the function \code{\link{find_cens}}.
}
  \item{start.formula}{
An optional argument giving an object of class \code{"formula"}: a symbolic description of the 
starting model in the MCMC algorithm. If \code{NULL} (the default) the starting model will be the maximal 
model. 
}
  \item{start.beta}{
An optional argument giving the starting values of the log-linear parameters for the MCMC algorithm. 
It should be a vector of the same length as the number of log-linear parameters in the starting model 
implied by the argument \code{start.formula}. If \code{NULL} (the default) the starting value will be the 
posterior mode under the maximal model.
}
  \item{start.sig}{
An optional argument giving the starting value of sigma^2 (under the Sabanes-Bove & Held prior) for 
the MCMC algorithm when the argument of prior is \code{"SBH"}. If \code{NULL} (the default) the 
starting value will be one.
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

For \code{bictu}, if non-\code{NULL}, the function will save the MCMC output to external text 
files every \code{save} iterations. If \code{NULL} (the default), it will inherit the value of 
\code{save} from the previous call to \code{bict} or \code{bictu}.
}
  \item{name}{
An optional argument giving a prefix to the external files saved if the argument \code{save} is 
positive. For \code{bict}, a value of \code{NULL} means the external files will not have a prefix. For 
\code{bictu}, a value of \code{NULL}, means the prefix will be inherited from the previous call to 
\code{bict} or \code{bictu}. 
}
  \item{null.move.prob}{
An optional scalar argument giving the probability of performing a null move in the reversible jump 
algorithm, i.e. proposing a move to the current model. The default value is 0.5.
}
  \item{a}{
The shape hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014). The default 
value is 0.001.
}
  \item{b}{
The scale hyperparameter of the Sabanes-Bove & Held prior, see Overstall & King (2014). The default 
value is 0.001.
}
  \item{progress}{
Logical argument. If \code{TRUE}, then a progress bar will be displayed. The default value is 
\code{FALSE}.
}
}
\details{
For identifiability, the parameters are constrained. The \code{\link{conting-package}} uses sum-to-zero constraints. 
See Overstall & King (2014), and the references therein, for more details. 

The Metropolis-Hastings algorithm employed is the iterated weighted least squares method for 
generalised linear models (GLMs) proposed by Gamerman (1997). The reversible jump algorithm employed is 
the orthogonal projections method for GLMs proposed by Forster et al (2012). For details on these 
methods applied to log-linear models through the data-augmentation algorithm see Overstall & King (2014), 
and the references therein. For details on the censored approach see Overstall et al (2014).

For details on the unit information and Sabanes-Bove & Held priors for generalised linear models see Ntzoufras 
et al (2003) and Sabanes-Bove & Held (2011), respectively. See Overstall & King (2014), and the references therein, 
for their application to log-linear models and contingency tables.
}
\value{
The functions will return an object of class \code{"bict"} which is a list with the following components.

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
the \code{data.frame} in \code{data} (or the elements of the variables) which correspond 
to the missing cell counts.
}

\item{missing2}{
A vector of the same length as the number of censored cell counts giving the row numbers of 
the \code{data.frame} in \code{data} (or the elements of the variables) which correspond to 
the censored cell counts.
}

\item{missing_details}{
The rows of the \code{data.frame} in \code{data} (or the elements of the variables) which 
correspond to the missing cell counts.
}

\item{censored_details}{
The rows of the \code{data.frame} in \code{data} (or the elements of the variables) which 
correspond to the censored cell counts.
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
A numeric scalar indicating which prior was used: 1 = \code{"UIP"}, 2 = \code{"SBH"}.
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
  \item{a}{
The argument \code{a}.
}
  \item{b}{
The argument \code{b}.
}
}
\references{
Sabanes-Bove, D. & Held, L. (2011) Hyper-g priors for generalized linear models. \emph{Bayesian
Analysis}, \bold{6 (3)}, 387--410.

Forster, J.J., Gill, R.C. & Overstall, A.M. (2012) Reversible jump methods for generalised linear models and generalised 
linear mixed models. \emph{Statistics and Computing}, \bold{22 (1)}, 107--120.

Gamerman, D. (1997) Sampling from the posterior distribution in generalised linear mixed
models. \emph{Statistics and Computing}, \bold{7 (1)}, 57--68.

Nztoufras, I., Dellaportas, P. & Forster, J.J. (2003) Bayesian variable and link 
determination for generalised linear models. \emph{Journal of Statistical Planning and 
Inference}, \bold{111 (1)}, 165--180.

Overstall, A.M. & King, R. (2014) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Journal of Statistical Software}, \bold{58 (7)}, 
1--27. \url{http://www.jstatsoft.org/v58/i07/}

Overstall, A.M., King, R., Bird, S.M., Hutchinson, S.J. & Hay, G. (2014) Incomplete contingency 
tables with censored cells with application to estimating the number of people who inject drugs 
in Scotland. \emph{Statistics in Medicine}, \bold{33 (9)}, 1564--1579.
}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}.
}
\note{
These functions are wrappers for \code{\link{bict.fit}}.

In Version 1.0 of \code{\link{conting-package}}, note that the default value for \code{prior} was \code{"UIP"}. From 
Version 1.1 onwards, the default value is \code{"SBH"}. 
}

\seealso{
\code{\link{bict.fit}},
\code{\link{spina}},
\code{\link{ScotPWID}}.
}
\examples{
set.seed(1)
## Set seed for reproducibility.

data(spina)
## Load the spina data

test1<-bict(formula=y~(S1 + S2 + S3 + eth)^2,data=spina,n.sample=50, prior="UIP")
## Let the maximal model be the model with two-way interactions. Starting from the 
## posterior mode of the model with two-way interactions, do 50 iterations under the 
## unit information prior.

test1<-bictu(object=test1,n.sample=50)
## Do another 50 iterations

test1

#Number of cells in table = 24
#
#Maximal model =
#y ~ (S1 + S2 + S3 + eth)^2
#
#Number of log-linear parameters in maximal model = 15 
#
#Number of MCMC iterations = 100 
#
#Computer time for MCMC = 00:00:01 
#
#Prior distribution for log-linear parameters = UIP 
#
#Number of missing cells = 3 
#
#Number of censored cells = 0

summary(test1)
## Summarise the result. Will get:

#Posterior summary statistics of log-linear parameters:
#            post_prob post_mean post_var lower_lim upper_lim
#(Intercept)         1    1.0427 0.033967    0.6498    1.4213
#S11                 1   -0.3159 0.015785   -0.4477   -0.1203
#S21                 1    0.8030 0.018797    0.6127    1.1865
#S31                 1    0.7951 0.003890    0.6703    0.8818
#eth1                1    2.8502 0.033455    2.4075    3.1764
#eth2                1    0.1435 0.072437   -0.4084    0.5048
#S21:S31             1   -0.4725 0.002416   -0.5555   -0.3928
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively
#
#Posterior model probabilities:
#  prob model_formula                                                         
#1 0.36 ~S1 + S2 + S3 + eth + S2:S3                                           
#2 0.19 ~S1 + S2 + S3 + eth + S2:S3 + S2:eth                                  
#3 0.12 ~S1 + S2 + S3 + eth + S1:eth + S2:S3                                  
#4 0.12 ~S1 + S2 + S3 + eth + S1:S2 + S1:S3 + S1:eth + S2:S3 + S2:eth + S3:eth
#5 0.10 ~S1 + S2 + S3 + eth + S1:S3 + S1:eth + S2:S3                          
#6 0.06 ~S1 + S2 + S3 + eth + S1:S3 + S1:eth + S2:S3 + S2:eth                 
#
#Total number of models visited =  8 
#
#Posterior mean of total population size = 726.75 
#95 % highest posterior density interval for total population size = ( 706 758 ) 
#
#Under the X2 statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  8.329  15.190  20.040  22.550  24.180 105.200 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  5.329  18.270  22.580  21.290  24.110  37.940 
#
#Bayesian p-value =  0.45


}
