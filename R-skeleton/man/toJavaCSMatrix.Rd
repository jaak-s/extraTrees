\name{toJavaCSMatrix}
\alias{toJavaCSMatrix}
\title{Utility function for converting an R SparseMatrix (package Matrix) to Java (column) sparse matrix.}
\description{
  Internal function used for converting an R SparseMatrix (package Matrix) to a CSparseMatrix object in Java. CSparseMatrix class is a custom Java class used for storing sparse matrices by the implementation of ExtraTrees in Java.
}
\usage{
  toJavaCSMatrix( m )
}
\arguments{
  \item{m}{ matrix of numeric values. }
}
\value{
  reference to Java matrix with the same contents as the input R matrix.
}
\author{Jaak Simm}
\keyword{java,matrix,conversion}

