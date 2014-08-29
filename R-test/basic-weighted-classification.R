
library("extraTrees")

## test setup: positive points + unlabelled points.
## we treat the unlabelled points as negative and weak weight.
nPos = 100
nUnlabelled = 700
n = nPos + nUnlabelled
p = 4
fx.pos = function(n) matrix(runif(n*p, min=c(0.7, 0.75, rep(0,p-2)), max=1), n, p, byrow=TRUE)
fx.all = function(n) matrix(runif(n*p), n, p)
x = rbind( fx.pos(nPos), fx.all(nUnlabelled) )
f = function(x) (x[,1]>0.7 & x[,2]>0.75) + 0
y = as.factor( c(rep(1, nPos), rep(0, nUnlabelled) ) )
w = c(rep(1, nPos), rep(0.1, nUnlabelled) )

xtest = fx.all( 10000 )
ytest = as.factor( f(xtest) )

## learning extra trees:
et = extraTrees(x, y, nodesize=2, mtry=p-1, numRandomCuts=2)
et2 = extraTrees(x[1:200,], y[1:200], nodesize=2, mtry=p-1, numRandomCuts=2)
etw = extraTrees(x, y, nodesize=8, mtry=p-1, numRandomCuts=2, weights=w )
yhat.et = predict(et, xtest)
yhat.et2 = predict(et2, xtest)
yhat.etw = predict(etw, xtest)

pind = which(ytest==1)
results = data.frame(
	accuracy = c( 
		et  = mean(ytest==yhat.et),
		et  = mean(ytest==yhat.et2),
		etw = mean(ytest==yhat.etw)
	),
	recall = c(
		et  = mean(ytest[pind]==yhat.et[pind]),
		et2 = mean(ytest[pind]==yhat.et2[pind]),
		etw = mean(ytest[pind]==yhat.etw[pind])
	)
)
print( results )
