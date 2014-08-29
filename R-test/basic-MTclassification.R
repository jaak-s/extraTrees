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

## single-task learning, all separated
train.sep <- function(X, y, tasks) {
    x = tapply( 1:nrow(X), tasks, function(i) {
        ## training separate extraTrees
        et = extraTrees( X[i,,drop=FALSE], y[i], nodesize=1, mtry=p, numRandomCuts=1 )
        return(et)
    })
    return(x)
}

## prediction with separately learned trees
predict.sep <- function(etsep, newdata, newtasks) {
    res = sapply( 1:length(newtasks), function(i) {
        task  = newtasks[i]
        input = newdata[i,,drop=FALSE]
        et    = etsep[[ task ]]
        yhat  = predict(et, input)
        return( yhat )
    })
    if (etsep[[1]]$factor) {
        return( as.factor(res) )
    }
    return( res )
}


## learning extra trees with multi-task learning:
et = extraTrees(Dtrain$x, Dtrain$y, nodesize=1, mtry=p, numRandomCuts=1, tasks=Dtrain$tasks)
et0 = extraTrees(Dtrain$x, Dtrain$y, nodesize=1, mtry=p, numRandomCuts=1)
ets = train.sep(Dtrain$x, Dtrain$y, Dtrain$tasks)

yhat = predict(et, Dtest$x, newtasks=Dtest$tasks)
## testing allValues
yhatAll = predict(et, Dtest$x[1:10,], newtasks=Dtest$tasks[1:10], allValues=T)

## pooled
yhat0 = predict(et0, Dtest$x)
## separated
yhatSep = predict.sep(ets, Dtest$x, newtasks=Dtest$tasks)

print( sprintf( "accuracy(extraTreesMT): %f", mean(Dtest$y==yhat) ) )
print( sprintf( "accuracy(pooled extraTrees): %f", mean(Dtest$y==yhat0) ) )
print( sprintf( "accuracy(separ. extraTrees): %f", mean(Dtest$y==yhatSep) ) )


