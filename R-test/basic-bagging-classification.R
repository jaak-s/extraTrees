
library("extraTrees")

## test setup: positive points + unlabelled points.
## we treat the unlabelled points as negative and weak weight.
nPos = 50
nUnlabelled = 7*nPos
n = nPos + nUnlabelled
p = 4
fx.pos = function(n) matrix(runif(n*p, min=c(0.7, 0.75, rep(0,p-2)), max=1), n, p, byrow=TRUE)
fx.all = function(n) matrix(runif(n*p), n, p)
x = rbind( fx.pos(nPos), fx.all(nUnlabelled) )
f = function(x) (x[,1]>0.7 & x[,2]>0.75) + 0
y = as.factor( c(rep(1, nPos), rep(0, nUnlabelled) ) )
w = c(rep(1, nPos), rep(0.1, nUnlabelled) )
bagLabels = c(rep(1, nPos), rep(2, nUnlabelled))
bagSizes  = c(nPos, nPos)

xtest = fx.all( 10000 )
ytest = as.factor( f(xtest) )

## learning extra trees:
methods = list()
methods$et = extraTrees(x, y, nodesize=2, mTry=p-1, numRandomCuts=2)
methods$et2 = extraTrees(x[1:(2*nPos),], y[1:(2*nPos)], nodesize=2, mTry=p-1, numRandomCuts=2)
methods$etw = extraTrees(x, y, nodesize=8, mTry=p-1, numRandomCuts=2, weights=w )
methods$etb = extraTrees(x, y, nodesize=2, mTry=p-1, numRandomCuts=2, bagLabels=bagLabels, bagSizes=bagSizes )
methods$etwb = extraTrees(x, y, nodesize=8, mTry=p-1, numRandomCuts=2, bagLabels=bagLabels, bagSizes=c(nPos, 4*nPos), 
                          weights=c(rep(1, nPos), rep(0.25, nUnlabelled)) )

yhat = lapply( methods, predict, xtest)

pind = which(ytest==1)

results = data.frame(
  accuracy = sapply(yhat, function(yh) mean(ytest==yh) ),
  recall   = sapply(yhat, function(yh) mean(ytest[pind]==yh[pind]) )
)
print(results)


