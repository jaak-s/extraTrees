library(extraTrees)
library(randomForest)

testing <- function(f, x, y, xtest, ytest) {
    time.f = system.time( {m = f(x, y)} )
    yhat   = predict(m, xtest)
    return( mean( ytest==yhat ) )
    ##return( c(sqError=yerr.f, time=time.f[3]) )
}

methods = c(
    extraTrees.1c = function(x,y) extraTrees(x,y),
    extraTrees.1c4try = function(x,y) extraTrees(x,y, mtry=4),
    extraTrees.2c = function(x,y) extraTrees(x,y, numRandomCuts=2),
    extraTrees.2e = function(x,y) extraTrees(x,y, numRandomCuts=2, evenCuts=TRUE),
    extraTrees.5c = function(x,y) extraTrees(x,y, numRandomCuts=5),
    extraTrees.5e = function(x,y) extraTrees(x,y, numRandomCuts=5, evenCuts=TRUE),
    randomForest  = function(x,y) randomForest(x,y)
)

data(iris)
X = sapply( 1:40, function(i) {
    train = rep(TRUE, nrow(iris))
    ntrain = 35 # nrow(iris)/3
    train[ sample.int( nrow(iris), nrow(iris)-ntrain ) ] = FALSE
    x = iris[train, 1:4]
    y = iris[train, 5]
    xtest = iris[!train, 1:4]
    ytest = iris[!train, 5]

    results = sapply( methods, testing, x=x, y=y, xtest=xtest, ytest=ytest )

    return ( results )
})
print( rowMeans(X) )

