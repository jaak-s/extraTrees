
library("extraTrees")
library("randomForest")

## making data:
n = 800
p = 4
x = matrix(runif(n*p), n, p)
f = function(x) (x[,1]>0.5) + (x[,2]>0.6) + (x[,3]>0.4)
y = as.factor(f(x))

xtest = matrix(runif(n*p), n, p)
ytest = as.factor( f(xtest) )

## learning extra trees:
et = extraTrees(x, y, nodesize=1, mtry=p, numRandomCuts=4)
yhat = predict(et, xtest)
print( sprintf( "accuracy(extraTrees): %f", mean(ytest==yhat) ) )

rf = randomForest(x, y)
yhat.rf = predict(rf, xtest)
print( sprintf( "accuracy(randomForest): %f", mean(ytest==yhat.rf) ) )

y2 = predict(et, x)
print( mean(y2==y) )
#yerr = mean( (y-yhat)^2 )


