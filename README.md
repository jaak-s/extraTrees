extraTrees
==========

ExtraTrees method for Java and R. ExtraTrees trains an ensemble of binary decision trees for *classification* and *regression*. ExtraTrees is very closely related to RandomForest. 

The software is available in R (2.15.2 and up).

```R
## Installing:
install.package("extraTrees")

## Loading ExtraTrees:
library("extraTrees")

## Generating artificial data:
n <- 1000  ## number of samples
p <- 5     ## number of dimensions
## x - inputs
## y - outputs to predict
x <- matrix(runif(n*p), n, p)
y <- (x[,1]>0.5) + 0.8*(x[,2]>0.6) + 0.5*(x[,3]>0.4) +
     0.1*runif(nrow(x))

## Training and predicting with ExtraTrees:
et <- extraTrees(x, y)
yhat <- predict(et, x)
```

## Development setup
For Java development checkout the git repository and follow [development.md](development.md).

