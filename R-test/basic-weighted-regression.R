
library("extraTrees")

## test setup: positive points (1.0) + unlabelled points.
## we treat the unlabelled points as 0.0 and use weak weight.
nPos = 100
nUnlabelled = 700
n = nPos + nUnlabelled
p = 4
fx.pos = function(n) matrix(runif(n*p, min=c(0.7, 0.75, rep(0,p-2)), max=1), n, p, byrow=TRUE)
fx.all = function(n) matrix(runif(n*p), n, p)
x = rbind( fx.pos(nPos), fx.all(nUnlabelled) )
f = function(x) (x[,1]>0.7 & x[,2]>0.75) + 0
y = c(rep(1, nPos), rep(0, nUnlabelled) )
w = c(rep(1, nPos), rep(0.1, nUnlabelled) )

xtest = fx.all( 10000 )
ytest = f(xtest)

## learning extra trees:
et = extraTrees(x, y, nodesize=5, mtry=p-1, numRandomCuts=2)
et2 = extraTrees(x[1:200,], y[1:200], nodesize=5, mtry=p-1, numRandomCuts=2)
etw = extraTrees(x, y, nodesize=15, mtry=p-1, numRandomCuts=2, weights=w )
etw2 = extraTrees(x[1:200,], y[1:200], nodesize=5, mtry=p-1, numRandomCuts=2, weights=w[1:200])
yhat.et = predict(et, xtest)
yhat.et2 = predict(et2, xtest)
yhat.etw = predict(etw, xtest)
yhat.etw2 = predict(etw2, xtest)

pind = which(ytest==1)
results = data.frame(
	mse = c( 
		et  = mean( (ytest-yhat.et)^2 ),
		et2 = mean( (ytest-yhat.et2)^2 ),
		etw = mean( (ytest-yhat.etw)^2 ),
		etw2= mean( (ytest-yhat.etw2)^2 )
	),
	mse1 = c(
		et  = mean( (ytest[pind]-yhat.et[pind])^2 ),
		et2 = mean( (ytest[pind]-yhat.et2[pind])^2 ),
		etw = mean( (ytest[pind]-yhat.etw[pind])^2 ),
		etw2= mean( (ytest[pind]-yhat.etw2[pind])^2 )
	)
)
print( results )
