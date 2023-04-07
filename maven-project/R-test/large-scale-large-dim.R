library(extraTrees)
library(randomForest)

## making data:
n = 1000
p = 100
f = function(x) (x[,1]<0.4) + 0.8*rowSums(x[,2:20]>0.6) + 0.5*(x[,50]>0.4) + 0.1*runif(nrow(x))
x = matrix(runif(n*p), n, p)
y = as.numeric(f(x))
xtest = matrix(runif(n*p), n, p)
ytest = f(xtest)

methods = c(
    extraTrees.1c = function(x,y) extraTrees(x,y, nodesize=5),
    extraTrees.1c2t = function(x,y) extraTrees(x,y, nodesize=5, numThreads=2),
    extraTrees.2c = function(x,y) extraTrees(x,y, nodesize=5, numRandomCuts=2),
    extraTrees.5c = function(x,y) extraTrees(x,y, nodesize=5, numRandomCuts=5)
    #randomForest  = function(x,y) randomForest(x,y)
)

## tests a learning function f: (returns sqError and time)
testing <- function(f, x, y, xtest, ytest) {
    time.f = system.time( {m = f(x, y)} )
    yhat   = predict(m, xtest)
    yerr.f = mean( (ytest-yhat)^2 )
    return( c(sqError=yerr.f, time=time.f[3]) )
}

results = t(sapply( methods, testing, x=x, y=y, xtest=xtest, ytest=ytest ))

print( results )

