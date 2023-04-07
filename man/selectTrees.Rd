\name{selectTrees}
\alias{selectTrees}
\title{Makes a sub-ExtraTrees object by keeping only selected trees.}
\description{
  This function creates a sub-ExtraTrees object by keeping only selected trees specified by selection.
}
\usage{
  selectTrees(object, selection)
}
\arguments{
  \item{object}{ extraTrees (S3) object, created by extraTrees(). }
  \item{selection}{ a list of logicals (T/F) of length object$ntree. }
}
\value{
  A new ExtraTrees (S3) object based on the existing object by keeping only the trees
  present in the selection.
}
\author{Jaak Simm}
\examples{
  ## Regression with ExtraTrees:
  n <- 1000  ## number of samples
  p <- 5     ## number of dimensions
  x <- matrix(runif(n*p), n, p)
  y <- (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4) + 0.1*runif(nrow(x))
  et <- extraTrees(x, y, nodesize=3, mtry=p, numRandomCuts=2, ntree=500)
  ## random selection of trees:
  trees <- sample(c(FALSE, TRUE), replace=TRUE, et$ntree)
  et2   <- selectTrees(et, selection=trees)
}
\keyword{regression,classification,trees}
