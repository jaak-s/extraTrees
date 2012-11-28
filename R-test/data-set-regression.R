library(extraTrees)
library(randomForest)


testing <- function(f, x, y, xtest, ytest) {
    time.f = system.time( {m = f(x, y)} )
    yhat   = predict(m, xtest)
    return( mean( (ytest-yhat)^2 ) )
    ##return( c(sqError=yerr.f, time=time.f[3]) )
}

methods = c(
    extraTrees.1c = function(x,y) extraTrees(x,y),
    extraTrees.2c = function(x,y) extraTrees(x,y, numRandomCuts=2),
    #extraTrees.2e = function(x,y) extraTrees(x,y, numRandomCuts=2, evenCuts=TRUE),
    extraTrees.5c = function(x,y) extraTrees(x,y, numRandomCuts=5),
    #extraTrees.5e = function(x,y) extraTrees(x,y, numRandomCuts=5, evenCuts=TRUE),
    randomForest  = function(x,y) randomForest(x,y)
)

datafilename = "http://personality-project.org/R/datasets/psychometrics.prob2.txt"
dataset = read.table(datafilename,header=TRUE)
datax = dataset[,2:8]
datay = dataset[,9]

X = sapply( 1:2, function(i) {
    train = rep(TRUE, nrow(dataset))
    ntrain = 500 # nrow(iris)/3
    train[ sample.int( nrow(dataset), nrow(dataset)-ntrain ) ] = FALSE
    x = datax[train,]
    y = datay[train]
    xtest = datax[!train, ]
    ytest = datay[!train]

    results = sapply( methods, testing, x=x, y=y, xtest=xtest, ytest=ytest )

    return ( results )
})
print( rowMeans(X) )

