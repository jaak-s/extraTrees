\name{toJavaMatrix}
\alias{toJavaMatrix}
\title{Utility function for converting an R matrix (numeric matrix) to Java matrix.}
\description{
  Internal function used for converting an R matrix to a Matrix object in Java. Matrix class is a custom Java class used for storing matrices by the implementation of ExtraTrees in Java.
}
\usage{
  toJavaMatrix( m )
}
\arguments{
  \item{m}{ matrix of numeric values. }
}
\value{
  reference to Java matrix with the same contents as the input R matrix.
}
\author{Jaak Simm}
\keyword{java,matrix,conversion}

