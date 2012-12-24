

## converts Java matrix to R matrix (doubles):
toRMatrix <- function( javam ) {
    return( matrix( 
        .jfield(javam, "[D", "v" ),
        .jfield(javam, "I", "nrows" ),
        .jfield(javam, "I", "ncols" )
    ))
}

## converts R matrix (or data.frame) into Java matrix (doubles)
toJavaMatrix <- function( m ) {
    if ( !is.matrix(m) ) { m = as.matrix(m) }
    if ( !is.double(m) ) { m = m + 0.0 }
    return(.jnew(
        "org.extratrees.Matrix", .jarray(m), nrow(m), ncol(m)
    ))
}

## creates a new extraTrees object based on selection
selectTrees <- function( object, selection ) {
    ## checking if right object:
    if (!inherits(object, "extraTrees")) {
        stop("object not of class extraTrees")
    }
    ## checking selection is logical:
    if (!is.logical(selection)) {
        stop("selection should be list of logical (T/F) values.")
    }
    ## checking selection has correct length:
    if (length(selection)!=object$ntree) {
        stop(sprintf("Length of selection (%d) should be equal to the number of trees (%d)", length(selection), object$ntree))
    }
    ## choosing the correct class:
    if (object$factor) {
        etClass = "Lorg/extratrees/FactorExtraTrees;"
    } else {
        etClass = "Lorg/extratrees/ExtraTrees;"
    }
    ## copying S3 object and creating new java object:
    etNew = object
    etNew$jobject = .jcall(object$jobject, etClass, "selectTrees", .jarray(selection) )
    etNew$ntree = .jcall(etNew$jobject, "I", "getNumTrees" )
    return(etNew)
}

extraTrees.default <- function(x, y, 
             #xtest=NULL, ytest=NULL, 
             ntree=500,
             mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1,
             numRandomCuts = 1,
             evenCuts = FALSE,
             numThreads = 1,
             quantile = F,
             ...) {
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## making sure no NAs:
    if ( any(is.na(y)) ) {
        stop("Output vector y contains NAs.")
    }
    if ( any(is.na(x)) ) {
        stop("Input matrix x contains NAs.")
    }

    ## uncomment when xtest/ytest are used:
    #testdat <- !is.null(xtest)
    #if (testdat) {
    #    if (ncol(x) != ncol(xtest))
    #        stop("x and xtest must have same number of columns")
    #    ntest <- nrow(xtest)
    #    xts.row.names <- rownames(xtest)
    #}
    
    et <- list()
    et$ntree    = ntree
    et$nodesize = nodesize
    et$ndim     = p
    et$mtry     = mtry
    et$factor   = is.factor(y)
    et$numRandomCuts = numRandomCuts
    et$evenCuts = evenCuts
    et$numThreads = numThreads
    et$quantile = quantile
    class(et) = "extraTrees"

    if (et$factor && length(unique(y)) < 2) {
        stop("Need at least two classes to do classification.")
    }
    if (et$factor) {
        if (quantile) {
            stop("Option quantile cannot be used for classification.")
        }
        ## classification:
        et$levels = levels(y)
        #stop("classification with extraTrees is not yet implemented.")
        ## creating FactorExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.FactorExtraTrees",
            toJavaMatrix(x),
            .jarray( as.integer( as.integer(y)-1 ) )
        )
        .jcall( et$jobject, "V", "setnFactors", as.integer(length(et$levels)) )
    } else if (et$quantile) {
        ## quantile regression:
        et$jobject = .jnew(
            "org.extratrees.QuantileExtraTrees",
            toJavaMatrix(x),
            .jarray(y)
        )
    } else {
        ## regression:
        ## creating ExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.ExtraTrees",
            toJavaMatrix(x),
            .jarray(y)
        )
    }
    ## setting variables:
    .jcall( et$jobject, "V", "setNumRandomCuts", as.integer(et$numRandomCuts) )
    .jcall( et$jobject, "V", "setEvenCuts", et$evenCuts )
    .jcall( et$jobject, "V", "setNumThreads", as.integer(et$numThreads) )
    
    ## learning the trees (stored at the et$jobject)
    .jcall( et$jobject, "V", "learnTrees", as.integer(et$nodesize), as.integer(et$mtry), as.integer(et$ntree) )
    
    return( et )
}

