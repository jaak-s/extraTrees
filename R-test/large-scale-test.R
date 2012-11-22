library("extraTrees")

## making data:
n = 5000
p = 4
x = matrix(runif(n*p), n, p)
f = function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
y = as.numeric(f(x))

xtest = matrix(runif(n*p), n, p)
ytest = f(xtest)

## extraTrees
time.et = system.time( {et = extraTrees(x, y, nodesize=5, mtry=p)} )
yhat    = predict(et, xtest)
yerr.et = mean( (ytest-yhat)^2 )

## randomForest
library(randomForest)
time.rf = system.time( {rf = randomForest(x, y)} )
yhat    = predict(rf, newdata=xtest)
yerr.rf = mean( (ytest-yhat)^2 )

## results:
results = data.frame( 
    errorSq=c(extraTrees=yerr.et, randomForest=yerr.rf),
    time   =c(extraTrees=time.et[3], randomForest=time.rf[3])
)
print( results )

