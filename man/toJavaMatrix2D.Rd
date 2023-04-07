\name{toJavaMatrix2D}
\alias{toJavaMatrix2D}
\title{Utility function for converting an R matrix (standard matrix or SparseMatrix) to appropriate Java matrix object.}
\description{
  Internal function used for converting an R matrix to an appropriate object in Java.
  It uses toJavaMatrix() and toJavaCSMatrix() underneath and returns a reference to 
  general matrix representation in Java of type Array2D (interface).
}
\usage{
  toJavaMatrix2D( m )
}
\arguments{
  \item{m}{ matrix of numeric values. }
}
\value{
  reference to Java matrix (dense or sparse) with the same contents as the input R matrix.
}
\author{Jaak Simm}
\keyword{java,matrix,conversion}

