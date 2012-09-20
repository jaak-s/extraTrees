
library("extraTrees")

## making data:
n = 200
p = 10
x = matrix(runif(n*p), n, p)
f = function(x) (x[,1]>0.5) + x[,2] - x[,3] + 0.5*x[,4]
y = f(x)

xtest = matrix(runif(n*p), n, p)
ytest = f(xtest)

## learning extra trees:
et = extraTrees(x, y, nodesize=3)
yhat = predict(et, x)
yerr = mean( (y-yhat)^2 )

ytest_hat = predict(et, xtest)
ytest_err = mean( (ytest-ytest_hat)^2 )

print( sprintf("ytrain err: %f",yerr) )
print( sprintf("ytest err:  %f",ytest_err) )

