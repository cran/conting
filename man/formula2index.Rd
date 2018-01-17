\name{formula2index}
\alias{formula2index}
\alias{index2formula}

\title{
Convert Between Formula and Index
}
\description{
These functions will convert a formula object to a binary index and vice versa.
}
\usage{
formula2index(big.X, formula, data)

index2formula(index, maximal.mod)
}

\arguments{
  \item{big.X}{
The design matrix under the maximal model.
}
  \item{formula}{
An object of class \code{"formula"}: a symbolic description of the 
model to convert to a binary index.
}
  \item{data}{
An object of class \code{"data.frame"} containing the variables in the model.
}

  \item{maximal.mod}{
An object of class \code{"glm"} giving the fit of the maximal model.
}

  \item{index}{
A binary vector, of the same length as the number of log-linear parameters in the maximal model, indicating 
which parameters are present in the model to be converted to a formula.
}
}

\value{
The function \code{formula2index} will produce a binary vector of the same length as the number of log-linear 
parameters in the maximal model, indicating which parameters are present in the model represented by the 
argument \code{formula}.

The function \code{index2formula} will produce an object of class \code{"formula"}: a symbolic description of the model given by 
the argument \code{index}.
}

\author{
Antony M. Overstall \email{A.M.Overstall@soton.ac.uk}.
}
\note{
This function will not typically be called by the user.
}

\examples{
data(ScotPWID)
## Load the ScotPWID data

maximal.mod<-glm(y~(S1+S2+S3+S4+Region+Gender+Age)^2,family=poisson,contrasts=list(
S1="contr.sum",S2="contr.sum",S3="contr.sum",S4="contr.sum",
Region="contr.sum",Gender="contr.sum",Age="contr.sum"),data=ScotPWID,x=TRUE)
## Fit the maximal model containing all two-way interactions.

big.X<-maximal.mod$x
## Set the design matrix under the maximal model

index<-formula2index(big.X=big.X,
formula=~S1+S2+S3+S4+Region+Gender+Age+S1:S2+S1:Age+S2:Gender+S3:S4+S4:Age,
data=ScotPWID)
## Find the index under the model with the following interactions:
## S1:S2
## S1:Age
## S2:Gender
## S3:S4
## S4:Age

index
## Print index
# [1] 1 1 1 1 1 1 1 1 1 0 0 0 0 1 0 0 0 1 0 1 0 0 0 0 0 1 0 0 0

index2formula(index=index,maximal.mod=maximal.mod)
## Go back to formula
#y ~ S1 + S2 + S3 + S4 + Region + Gender + Age + S1:S2 + S1:Age + 
#    S2:Gender + S3:S4 + S4:Age
}