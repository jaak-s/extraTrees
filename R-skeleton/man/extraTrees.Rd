\name{extraTrees}
\alias{extraTrees}
\title{Function for training ExtraTree classifier.}
\description{
  This function executes ExtraTree building method (implemented in Java).
}
\usage{
  extraTrees(x, y, ntree=500,
             mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1
             )
}
\arguments{
  \item{x}{ a numberic input data matrix, each row is an input. }
  \item{y}{ a vector of output values: if vector of numbers then regression, if vector of factors then classification. }
  \item{ntree}{ the number of trees (default 500). }
  \item{mtry}{ the number of features tried at each node (default is ncol(x)/3 for regression and sqrt(ncol(x)) for classification). }
  \item{nodesize}{ the size of leaves of the tree (default is 5 for regression and 1 for classification) }
  \item{numRandomCuts}{ the number of random cuts for each (randomly chosen) feature (default 1, which corresponds to the official ExtraTrees method). The higher the number of cuts the higher the chance of a good cut. }
}
%\details{
%  Details are given in the package vignette.
%}
\value{
  The trained model  from input x and output values y, stored in ExtraTree object.
}
\author{Jaak Simm}
\examples{
  ## Regression with ExtraTrees:
  n = 1000  ## number of samples
  p = 5     ## number of dimensions
  x = matrix(runif(n*p), n, p)
  y = (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
  et = extraTrees(x, y, nodesize=3, mtry=p, numRandomCuts=2)
}
\keyword{regression,classification,trees}

