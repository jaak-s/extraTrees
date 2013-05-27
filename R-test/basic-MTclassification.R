library("extraTrees")

## making data:
makeData <- function(n = 1000, p = 10, ntasks = 50) {
    out = list()
    x = matrix(runif(n*p), n, p)
    out$x    = x
    out$tasks = rep( 1:ntasks, length.out=nrow(x) )
    y = rep( 1, n)
    for (i in 1:n) {
        if (x[i,1] < 0.5) {
            if (out$tasks[i] %% 2 == 0) {
                y[i] = runif(1)<0.05
            } else {
                y[i] = runif(1)<0.95
            }
        } else {
            y[i] = x[i,2]<0.5
        }
    }
    out$y = as.factor(y)
    return(out)
}
Dtrain = makeData()
Dtest  = makeData()

## learning extra trees:
et = extraTrees(Dtrain$x, Dtrain$y, nodesize=1, mTry=p, numRandomCuts=10, tasks=Dtrain$tasks)
yhat = predict(et, Dtest$x, newtasks=Dtest$tasks)
print( sprintf( "accuracy(extraTrees): %f", mean(ytest==yhat) ) )
