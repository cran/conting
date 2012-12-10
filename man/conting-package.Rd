\name{conting-package}
\alias{conting-package}
\alias{conting}
\docType{package}
\title{
Bayesian Analysis of Complete and Incomplete Contingency Tables
}
\description{
Performs Bayesian analysis of complete and incomplete contingency tables incorporating 
model uncertainty using log-linear models. These analyses can be used to identify 
associations/interactions between categorical factors and to estimate unknown closed 
populations.
}
\details{
\tabular{ll}{
Package: \tab conting\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-11-17\cr
License: \tab GPL-2\cr
}
For the Bayesian analysis of complete contingency tables the key function is \code{\link{bcct}} 
which uses MCMC methods to generate a sample from the joint posterior distribution of the model parameters and model 
indicator. Further MCMC iterations can be performed by using \code{\link{bcctu}}.

For the Bayesian analysis of incomplete contingency tables the key function is \code{\link{bict}} 
which uses MCMC methods to generate a sample from the joint posterior distribution of the model parameters, model 
indicator and the missing and/or censored cell entries. Further MCMC iterations can be performed 
by using \code{\link{bictu}}.

In both cases see Overstall & King (2012), and the references therein, for details on the statistical 
and computational methods, as well as detailed examples.

}
\author{
Antony M. Overstall \email{antony@mcs.st-and.ac.uk}

Maintainer: Antony M. Overstall \email{antony@mcs.st-and.ac.uk}
}
\references{
Overstall, A.M. & King, R. (2012) conting: An R package for Bayesian analysis of
complete and incomplete contingency tables. \emph{Technical Report}, University 
of St Andrews, UK.
}

\keyword{ package }

\examples{
set.seed(1)
## Set seed for reproducibility
data(AOH)
## Load AOH data
test1<-bcct(formula=y~(alc+hyp+obe)^3,data=AOH,n.sample=100)
## Bayesian analysis of complete contingency table. Let the saturated model
## be the maximal model and do 100 iterations.

summary(test1)
## Summarise the result. Will get:
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

set.seed(1)
## Set seed for reproducibility
data(rabbit)
## Load rabbit data
test2<-bict(formula=y~(first+second+hunter)^2,data=rabbit,n.sample=100)
## Bayesian analysis of incomplete contingency table. Let the model with two-way 
## interactions be the maximal model and do 100 iterations.

summary(test2)
## Summarise the result. Will get:
#Posterior summary statistics of log-linear parameters:
#               post_prob post_mean post_var lower_lim upper_lim
#(Intercept)            1    2.8491 0.009506    2.6741    2.9971
#first1                 1   -0.4172 0.009777   -0.5679   -0.2332
#second1                1   -0.3432 0.008762   -0.5568   -0.2133
#hunter1                1   -0.7493 0.007090   -0.9160   -0.6181
#first1:second1         1    0.5711 0.006318    0.4333    0.7056
#NB: lower_lim and upper_lim refer to the lower and upper values of the
#95 % highest posterior density intervals, respectively
#
#Posterior model probabilities:
#  prob model_formula                                          
#1 0.71 ~first + second + hunter + first:second                
#2 0.14 ~first + second + hunter + first:second + second:hunter
#3 0.10 ~first + second + hunter + first:second + first:hunter 
#
#Total number of models visited =  4 
#
#Posterior mean of total population size = 261.19 
#95 % highest posterior density interval for total population size = ( 207 317 ) 
#
#Under the X2 statistic 
#
#Summary statistics for T_pred 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#  1.284   4.240   6.584   7.238   9.653  16.970 
#
#Summary statistics for T_obs 
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.7725  3.3950  5.0520  5.6300  6.1670 14.7200 
#
#Bayesian p-value =  0.68 

}
