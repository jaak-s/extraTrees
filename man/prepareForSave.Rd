\name{prepareForSave}
\alias{prepareForSave}
\title{Prepares ExtraTrees object for save() function}
\description{
  This function prepares ExtraTrees for saving by serializing the trees in Java VM.
  It is equivalent to calling \code{.jcache(et$jobject)}.
  Afterwards the object can be saved by \code{save} (or automatic R session saving) and will be fully recovered after \code{load}.

  Note: the object can still be used as usual after \code{prepareForSave}.
}
\usage{
  prepareForSave(object)
}
\arguments{
  \item{object}{ extraTrees (S3) object, created by extraTrees(). }
}
\value{
  Nothing is returned.
}
\author{Jaak Simm}
\examples{
  et <- extraTrees(iris[,1:4], iris$Species)
  prepareForSave(et)
  ## saving to a file
  save(et, file="temp.Rdata")

  ## testing: remove et and load it back from file
  rm(list = "et")
  load("temp.Rdata")
  predict(et, iris[,1:4])
}
\keyword{save,load,extraTrees}
