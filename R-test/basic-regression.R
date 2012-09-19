
library("extraTrees")

## making data:
n = 200
p = 10
x = matrix(runif(n*p), n, p)
y = (x[,1]>0.5) + x[,2] - x[,3] + 0.5*x[,4]

## learning extra trees:
et = extraTrees(x, y, nodesize=5)
yhat = predict(et, x)
yerr = mean( (y-yhat)^2 )

