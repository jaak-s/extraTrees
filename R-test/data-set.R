library(extraTrees)
library(randomForest)

data(iris)
x = iris[,1:4]
y = iris[,5]
et   <- extraTrees(x, y)
yhat <- predict(et, x)
print( sprintf( "accuracy(extraTrees): %f", mean(y==yhat) ) )

