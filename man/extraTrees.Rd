\name{extraTrees}
\alias{extraTrees}
\alias{extraTrees.default}
\title{Function for training ExtraTree classifier or regression.}
\description{
  This function executes ExtraTree building method (implemented in Java).
}
\usage{
  \method{extraTrees}{default}(x, y, 
             ntree=500,
             mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             numRandomCuts = 1,
             evenCuts = FALSE,
             numThreads = 1,
             quantile = F,
             weights = NULL,
             subsetSizes = NULL,
             subsetGroups = NULL,
             tasks = NULL,
             probOfTaskCuts = mtry / ncol(x),
             numRandomTaskCuts = 1,
             na.action = "stop",
             ...)
}
\arguments{
  \item{x}{ a numberic input data matrix, each row is an input. }
  \item{y}{ a vector of output values: if vector of numbers then regression, if vector of factors then classification. }
  \item{ntree}{ the number of trees (default 500). }
  \item{mtry}{ the number of features tried at each node (default is ncol(x)/3 for regression and sqrt(ncol(x)) for classification). }
  \item{nodesize}{ the size of leaves of the tree (default is 5 for regression and 1 for classification) }
  \item{numRandomCuts}{ the number of random cuts for each (randomly chosen) feature (default 1, which corresponds to the official ExtraTrees method). The higher the number of cuts the higher the chance of a good cut. }
  \item{evenCuts}{ if FALSE then cutting thresholds are uniformly sampled (default). If TRUE then the range is split into even intervals (the number of intervals is numRandomCuts) and a cut is uniformly sampled from each interval. }
  \item{numThreads}{ the number of CPU threads to use (default is 1). }
  \item{quantile}{ if TRUE then quantile regression is performed (default is FALSE), only for regression data.
                   Then use predict(et, newdata, quantile=k) to make predictions for k quantile. }
  \item{weights}{ a vector of sample weights, one positive real value for each sample. NULL means standard learning, i.e. equal weights. }
  \item{subsetSizes}{ subset size (one integer) or subset sizes (vector of integers, requires subsetGroups), if supplied every tree is built from a random subset of size subsetSizes. NULL means no subsetting, i.e. all samples are used. }
  \item{subsetGroups}{ list specifying subset group for each sample: from samples in group g, each tree will randomly select subsetSizes[g] samples. }
  \item{tasks}{ vector of tasks, integers from 1 and up. NULL if no multi-task learning }
  \item{probOfTaskCuts}{ probability of performing task cut at a node (default mtry / ncol(x)). Used only if tasks is specified. } 
  \item{numRandomTaskCuts}{ number of times task cut is performed at a node (default 1). Used only if tasks is specified. }
  \item{na.action}{ specifies how to handle NA in x: "stop" (default) will give error is any NA present, "zero" will set all NA to zero and "fuse" will build trees by skipping samples when the chosen feature is NA for them.} 
  \item{...}{ not used currently. }
}
\details{
  For classification ExtraTrees at each node chooses the cut based on minimizing the Gini impurity index and for regression the variance.

  For more details see the package vignette, i.e. vignette("extraTrees").
  
  If Java runs out of memory: \code{java.lang.OutOfMemoryError: Java heap space},
  then  (assuming you have free memory) you can increase the heap size by:
  \code{options( java.parameters = "-Xmx2g" )} before calling
  \code{library( "extraTrees" )}, where 2g defines 2GB of heap size. 
  Change it as necessary.
}
\value{
  The trained model from input x and output values y, stored in ExtraTree object.
}
\seealso{
   \code{\link{predict.extraTrees}} for predicting and \code{\link{prepareForSave}} for saving ExtraTrees models to disk.
}
\author{Jaak Simm}
\examples{
  ## Regression with ExtraTrees:
  n <- 1000  ## number of samples
  p <- 5     ## number of dimensions
  x <- matrix(runif(n*p), n, p)
  y <- (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4) +
       0.1*runif(nrow(x))
  et <- extraTrees(x, y, nodesize=3, mtry=p, numRandomCuts=2)
  yhat <- predict(et, x)
  
  #######################################
  ## Multi-task regression with ExtraTrees:
  n <- 1000  ## number of samples
  p <- 5     ## number of dimensions
  x <- matrix(runif(n*p), n, p)
  task <- sample(1:10, size=n, replace=TRUE)
  ## y depends on the task: 
  y <- 0.5*(x[,1]>0.5) + 0.6*(x[,2]>0.6) + 0.8*(x[cbind(1:n,(task \%\% 2) + 3)]>0.4)
  et <- extraTrees(x, y, nodesize=3, mtry=p-1, numRandomCuts=2, tasks=task)
  yhat <- predict(et, x, newtasks=task)
  
  #######################################
  ## Classification with ExtraTrees (with test data)
  make.data <- function(n) {
    p <- 4
    f <- function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
    x <- matrix(runif(n*p), n, p)
    y <- as.factor(f(x))
    return(list(x=x, y=y))
  }
  train <- make.data(800)
  test  <- make.data(500)
  et    <- extraTrees(train$x, train$y)
  yhat  <- predict(et, test$x)
  ## accuracy
  mean(test$y == yhat)
  ## class probabilities
  yprob = predict(et, test$x, probability=TRUE)
  head(yprob)
  
  #######################################
  ## Quantile regression with ExtraTrees (with test data)
  make.qdata <- function(n) {
    p <- 4
    f <- function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
    x <- matrix(runif(n*p), n, p)
    y <- as.numeric(f(x))
    return(list(x=x, y=y))
  }
  train <- make.qdata(400)
  test  <- make.qdata(200)
  
  ## learning extra trees:
  et <- extraTrees(train$x, train$y, quantile=TRUE)
  ## estimate median (0.5 quantile)
  yhat0.5 <- predict(et, test$x, quantile = 0.5)
  ## estimate 0.8 quantile (80\%)
  yhat0.8 <- predict(et, test$x, quantile = 0.8)


  #######################################
  ## Weighted regression with ExtraTrees 
  make.wdata <- function(n) {
    p <- 4
    f <- function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
    x <- matrix(runif(n*p), n, p)
    y <- as.numeric(f(x))
    return(list(x=x, y=y))
  }
  train <- make.wdata(400)
  test  <- make.wdata(200)
  
  ## first half of the samples have weight 1, rest 0.3
  weights <- rep(c(1, 0.3), each = nrow(train$x) / 2)
  et <- extraTrees(train$x, train$y, weights = weights, numRandomCuts = 2)
  ## estimates of the weighted model
  yhat <- predict(et, test$x)
}
\keyword{regression,classification,trees}

