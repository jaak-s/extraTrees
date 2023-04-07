
library("extraTrees")

## making data:
n = 400
p = 4
x = matrix(runif(n*p), n, p)
f = function(x) (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4)
#f = function(x) (x[,1]>0.5) + (x[,2]>0.7) - 0.5*(x[,3]<0.3) + (0.1<x[,4] & x[,4]<0.7)
#f = function(x) (x[,1]>0.5) + x[,2] - x[,3] + 0.5*x[,4]
y = as.numeric(f(x))

xtest = matrix(runif(n*p), n, p)
ytest = f(xtest)

## learning extra trees:
et = extraTrees(x, y, mtry=p, quantile=T)
yhat = predict(et, x)
yqhat = predict(et, x, quantile=0.5)
yerr = mean( (y-yhat)^2 )
yqerr = mean( (y-yqhat)^2 )

ytest_et = predict(et, xtest)
ytest_qet0.4 = predict(et, xtest, quantile=0.4)
ytest_qet0.5 = predict(et, xtest, quantile=0.5)
ytest_qet0.6 = predict(et, xtest, quantile=0.6)

## lm:
m = lm(y~., data=data.frame(x,y=y))
ytest_lm = predict(m, newdata=data.frame(xtest, y=ytest))

## randomForest:
library(randomForest)
rf = randomForest(x, y)
ytest_rf = predict(rf, newdata=xtest)

print( sprintf("train error(et): %f", yerr) )
print( sprintf("test error(et):  %f", mean((ytest-ytest_et)^2) ) )
print( sprintf("test er(qet0.4): %f", mean((ytest-ytest_qet0.4)^2) ) )
print( sprintf("test er(qet0.5): %f", mean((ytest-ytest_qet0.5)^2) ) )
print( sprintf("test er(qet0.6): %f", mean((ytest-ytest_qet0.6)^2) ) )
print( sprintf("test error(lm):  %f", mean((ytest-ytest_lm)^2) ) )
print( sprintf("test error(rf):  %f", mean((ytest-ytest_rf)^2) ) )


