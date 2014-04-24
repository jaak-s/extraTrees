library("extraTrees")

get.data <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
  y = factor(f(x), levels=0:3)
  return(list(x=x, y=y))
}

test.classification <- function() {
  train = get.data(100)
  test  = get.data(100)
  
  et = extraTrees(train$x, train$y, nodesize=1, mTry=2, numRandomCuts=2, ntree=50)
  yhat = predict(et, test$x)
  checkEquals( length(yhat), length(test$y) )
  checkEquals( levels(yhat), levels(test$y) )
  checkEquals( 50, et$ntree )
  checkEquals( TRUE, et$factor )
  
  yall = predict(et, test$x, allValues=T)
  checkEquals( nrow(yall), nrow(test$x) )
  checkEquals( ncol(yall), 50 )
  checkEquals( FALSE, is.numeric(yall) )
}

test.regression <- function() {
  train = get.data(100)
  test  = get.data(100)
  
  et = extraTrees(train$x, as.numeric(train$y), nodesize=1, mTry=2, numRandomCuts=2, ntree=50)
  yhat = predict(et, test$x)
  checkEquals( length(yhat), length(test$y) )
  checkEquals( 50,    et$ntree )
  checkEquals( TRUE,  is.numeric(yhat) )
  checkEquals( TRUE,  is.double(yhat) )
  checkEquals( FALSE, et$factor )
  
  yall = predict(et, test$x, allValues=T)
  checkEquals( nrow(yall), nrow(test$x) )
  checkEquals( ncol(yall), 50 )
  checkEquals( TRUE, is.double(yall) )
}

test.regressionint <- function() {
  train = get.data(100)
  test  = get.data(100)
  
  ## make inputs to int
  et = extraTrees(train$x, as.integer(train$y), ntree=50)
  yhat = predict(et, test$x)
  
  checkEquals( length(yhat), length(test$y) )
  checkEquals( TRUE,  is.double(yhat) )
}
