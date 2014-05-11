
context("Regression basics")

get.data.regr <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
  return(list(x=x, y=f(x)))
}

train <- get.data.regr(100)
test  <- get.data.regr(100)

test_that("basic regression and prediction", {  
  et    <- extraTrees(train$x, train$y, nodesize=1, mTry=2, numRandomCuts=2, ntree=50)
  yhat  <- predict(et, test$x)

  expect_is   ( et, "extraTrees")
  expect_equal( length(yhat), length(test$y) )
  expect_equal( 50,    et$ntree )
  expect_true(  is.numeric(yhat) )
  expect_true(  is.double(yhat) )
  expect_false( et$factor )
  
  yall <- predict(et, test$x, allValues=T)
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_true(  is.double(yall) )
})

test_that("integer y is used as double for regression", {  
  et   <- extraTrees(train$x, as.integer(train$y), ntree=50)
  yhat <- predict(et, train$x)
  
  expect_equal( length(yhat), length(train$y) )
  expect_true(  is.double(yhat) )
})

test_that("multi-threaded regression works", {
  et   <- extraTrees(train$x, train$y, numRandomCuts=2, ntree=50, numThreads=2)
  yhat <- predict(et, test$x)
  expect_equal( length(yhat), length(test$y) )
  expect_equal( 50, et$ntree )
  expect_equal( 2,  et$numThreads )
  expect_true(  is.numeric(yhat) )
  expect_true(  is.double(yhat) )
  expect_false( et$factor )
  
  yall = predict(et, test$x, allValues=T)
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_true(  is.double(yall) )
  
  expect_error( extraTrees(train$x, train$y, numThreads=0 ) )
})

test_that("training creates a usable call object", {
  et   <- extraTrees(train$x, train$y, ntree=50)
  expect_true( is.call(et$call) )
  et2  <- eval(et$call)
  expect_is(et2, "extraTrees")
})

