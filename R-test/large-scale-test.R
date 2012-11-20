library("extraTrees")

## making data:
n = 10000
p = 4
x = matrix(runif(n*p), n, p)
f = function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
y = as.numeric(f(x))

## extraTrees
time.et = system.time( {et = extraTrees(x, y, nodesize=3, mtry=p)} )
yhat    = predict(et, x)
yerr.et = mean( (y-yhat)^2 )

## randomForest
library(randomForest)
time.rf = system.time( {rf = randomForest(x, y)} )
yhat    = predict(rf, newdata=x)
yerr.rf = mean( (y-yhat)^2 )

## results:
results = data.frame( 
    errorSq=c(extraTrees=yerr.et, randomForest=yerr.rf),
    time   =c(extraTrees=time.et[3], randomForest=time.rf[3])
)
print( results )

