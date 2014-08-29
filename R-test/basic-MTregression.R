library("extraTrees")

## making data:
makeDataRegr <- function(n = 1000, p = 10, ntasks = 50) {
    out = list()
    x = matrix(runif(n*p), n, p)
    out$x    = x
    out$tasks = rep( 1:ntasks, length.out=nrow(x) )
    y = rep( 1, n)
    for (i in 1:n) {
        if (x[i,1] < 0.5) {
            if (out$tasks[i] %% 2 == 0) {
                ## for even numbered tasks
                y[i] = 0 + 0.1*rnorm(1) + x[i,3]
            } else {
                ## for odd numbered tasks
                y[i] = 1 + 0.1*rnorm(1) + x[i,4]
            }
        } else {
            ## no task effect if x1 < 0.5
            y[i] = 0.5 + 0.1*rnorm(1)
        }
    }
    out$y = y
    return(out)
}
Dtrain = makeDataRegr()
Dtest  = makeDataRegr()

## learning extra trees with multi-task learning:
et = extraTrees(Dtrain$x, Dtrain$y, nodesize=1, mtry=p, numRandomCuts=1, tasks=Dtrain$tasks)
et0 = extraTrees(Dtrain$x, Dtrain$y, nodesize=1, mtry=p, numRandomCuts=1)
yhat = predict(et, Dtest$x, newtasks=Dtest$tasks)
## testing allValues
yhatAll = predict(et, Dtest$x[1:10,], newtasks=Dtest$tasks[1:10], allValues=T)

yhat0 = predict(et0, Dtest$x)
print( sprintf( "meanSqError(extraTreesMT): %f",     mean((Dtest$y-yhat)^2) ) )
print( sprintf( "meanSqError(plain extraTrees): %f", mean((Dtest$y-yhat0)^2) ) )

