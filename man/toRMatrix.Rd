\name{toRMatrix}
\alias{toRMatrix}
\title{Utility function for converting Java matrix to R matrix (matrix of doubles).}
\description{
  Internal function used for converting a Matrix object from Java to an R matrix. Matrix class is a custom Java class used for storing matrices by the implementation of ExtraTrees in Java.
}
\usage{
  toRMatrix( javam )
}
\arguments{
  \item{javam}{ Java matrix (Matrix class). }
}
\value{
  R (double) matrix with the same contents as the input.
}
\author{Jaak Simm}
\keyword{java,matrix,conversion}

