\name{predict.extraTrees}
\alias{predict.extraTrees}
\title{Function for making predictions from trained ExtraTree object.}
\description{
  This function makes predictions for regression/classification using the given trained ExtraTree object and provided input matrix (newdata).
}
\usage{
  \method{predict}{extraTrees}(object, newdata, quantile=NULL, allValues=F, probability=F, newtasks=NULL, ...)
}
\arguments{
  \item{object}{ extraTree (S3) object, created by extraTrees(). }
  \item{newdata}{ a new numberic input data matrix, for each row a prediction is made. }
  \item{quantile}{ the quantile value between 0.0 and 1.0 for quantile regression, or NULL (default) for standard predictions. }
  \item{allValues}{ whether or not to return outputs of all trees (default FALSE). }
  \item{probability}{ whether to return a matrix of class (factor) probabilities, default FALSE. Can only be used in the case of classification. Calculated as the proportion of trees voting for particular class. }
  \item{newtasks}{ list of tasks, for each input in newdata (default NULL). Must be NULL if no multi-task learning was used at training. }
  \item{...}{ not used currently. }
}
\value{
  The vector of predictions from the ExtraTree et.
  The length of the vector is equal to the the number of rows in newdata.
}
\author{Jaak Simm}
\examples{
  ## Regression with ExtraTrees:
  n <- 1000  ## number of samples
  p <- 5     ## number of dimensions
  x <- matrix(runif(n*p), n, p)
  y <- (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4) + 0.1*runif(nrow(x))
  et <- extraTrees(x, y, nodesize=3, mtry=p, numRandomCuts=2)
  yhat <- predict(et, x)
}
\keyword{regression,classification,trees}

