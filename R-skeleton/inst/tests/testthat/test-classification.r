
context("Classification basics")

get.data.class <- function(n=800) {
  p = 4
  x = matrix(runif(n*p), n, p)
  f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
  y = factor(f(x), levels=0:3)
  return(list(x=x, y=y))
}

train <- get.data.class(100)
test  <- get.data.class(100)

test_that("basic classification and prediction", {  
  et   <- extraTrees(train$x, train$y, nodesize=1, mTry=2, numRandomCuts=2, ntree=50)
  yhat <- predict(et, test$x)
  
  expect_equal( length(yhat), length(test$y) )
  expect_equal( levels(yhat), levels(test$y) )
  expect_equal( 50, et$ntree )
  expect_true(  et$factor )  

  ## allValues prediction
  yall = predict(et, test$x, allValues=T)
  expect_equal( nrow(yall), nrow(test$x) )
  expect_equal( ncol(yall), 50 )
  expect_false( is.numeric(yall) )
})
