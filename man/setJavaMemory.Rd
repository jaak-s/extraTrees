\name{setJavaMemory}
\alias{setJavaMemory}
\title{Utility function for setting Java memory.}
\description{
  Function for setting JVM memory, specified in MB. If you get java.lang.OutOfMemoryError you can use this function to increase the memory available to ExtraTrees.
}
\usage{
  setJavaMemory( memoryInMB )
}
\arguments{
  \item{memoryInMB}{Integer specifying the amount of memory (MB)}
}
%\value{
%}
\author{Jaak Simm}
\examples{
  ## use 2G memory
  setJavaMemory(2000)
}
\keyword{java.lang.OutOfMemoryError,memory,JVM}

