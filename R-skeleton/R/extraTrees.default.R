

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
        "org.extratrees.data.Matrix", .jarray(m), nrow(m), ncol(m)
    ))
}

toJavaMatrix2D <- function( m ) {
  .jcast(toJavaMatrix(m), new.class="org/extratrees/data/Array2D")
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

## main extraTree training function
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
             weights = NULL,
             bagSizes = NULL,
             bagLabels = NULL,
             tasks = NULL,
             probOfTaskCuts = mtry / ncol(x),
             numRandomTaskCuts = 1,
             na.action = "stop",
             ...) {
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)
    
    ## making sure no NAs:
    if ( any(is.na(y)) ) stop("Output vector y contains NAs.")
    if ( !is.null(tasks) && any(is.na(tasks)) ) stop("Task vector contains NAs.")
    
    if ( numThreads < 1 ) stop("numThreads has to be 1 or bigger.")

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
    et$quantile   = quantile
    et$useWeights = ! is.null(weights)
    et$useBagging = ! is.null(bagSizes) && sum(bagSizes) != nrow(x)
    et$multitask  = ! is.null(tasks)
    et$probOfTaskCuts = probOfTaskCuts
    et$numRandomTaskCuts = numRandomTaskCuts
    et$call <- match.call()
    et$call[[1]] <- as.name("extraTrees")
    
    class(et) = "extraTrees"

    if (nrow(x) != length(y)) {
        stop(sprintf("Length of y (%d) is not equal to the number of samples in x (%d).", length(y), nrow(x) ) )
    }
    
    if ( ! et$factor && ! is.numeric(y) ) {
      stop("y values have to be either factor (classification) or numeric (regression).")
    }

    et$xHasNA = FALSE
    if (na.action == "stop") {
      if ( any(is.na(x)) ) stop("Input matrix x contains NAs. Change na.action to 'zero' or 'fuse' to allow NA to be used or remove samples with NA.")
    } else if (na.action == "zero") {
      x[ is.na(x) ] = 0
    } else if (na.action == "fuse") {
      et$xHasNA = any(is.na(x))
    } else {
      stop("na.action should be either 'stop', 'zero' or 'fuse'. See manual for details.")
    }
    
    if ( ! is.null(weights) ) {
      if (nrow(x) != length(weights)) {
        stop(sprintf("Length of weights (%d) is not equal to the number of samples in x (%d).", length(weights), nrow(x) ) )
      }
      if (any(weights <= 0)) {
        stop("Weights have to be positive.")
      }
    }
    
    if ( et$useBagging ) {
      if (sum(bagSizes) > nrow(x)) {
        stop(sprintf("Total of bagSizes (%d) should not be bigger than the number of samples in x (%d).", sum(bagSizes), nrow(x) ))
      }
      ## if only one bag size then bagLabels are not used
      if (length(bagSizes) >= 2) {
        if (nrow(x) != length(bagLabels)) {
          stop(sprintf("Length of bagLabels (%d) is not equal to the number of samples in x (%d).", length(weights), nrow(x) ) )
        }
        if ( ! is.factor(bagLabels) ) {
          bagLabels = as.factor( bagLabels )
        }
        
        numUnique = length(levels(bagLabels))
        if (numUnique != length(bagSizes)) {
          stop(sprintf("Number of unique bagLabels (%d) has to be the same as length of number of bagSizes (%d).", numUnique, length(bagSizes)  ))
        }
      }
    }
    
    ## making sure if tasks is present there are only two factors
    if ( ! is.null(tasks) ) {
        if (nrow(x)!=length(tasks)) {
            stop(sprintf("Length of tasks (%d) is not equal to the number of inputs in x (%d).", length(tasks), nrow(x) ) )
        }
        if ( et$factor && length(unique(y)) != 2 ) {
            stop("Multi-task learning only works with 2 factors (binary classification). 3 or more classes is not supported.")
        }
        if (min(tasks) < 1) {
            stop("Tasks should be positive integers.")
        }
        if ( et$quantile ) {
            stop("Quantile regression is not (yet) supported with multi-task learning.")
        }
        if ( et$useWeights ) {
            stop("Weights are not (yet) supported with multi-task learning.")
        }
    }

    if (et$factor && length(unique(y)) < 2) {
        stop("Need at least two classes to do classification.")
    }
    if (et$factor) {
        if (quantile) {
            stop("Option quantile cannot be used for classification.")
        }
        ## classification:
        et$levels = levels(y)
        ## creating FactorExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.FactorExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.integer( as.integer(y)-1 ) )
        )
        .jcall( et$jobject, "V", "setnFactors", as.integer(length(et$levels)) )
    } else if (et$quantile) {
        ## quantile regression:
        et$jobject = .jnew(
            "org.extratrees.QuantileExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.double(y) )
        )
    } else {
        ## regression:
        ## creating ExtraTree object with the data
        et$jobject = .jnew(
            "org.extratrees.ExtraTrees",
            toJavaMatrix2D(x),
            .jarray( as.double(y) )
        )
    }
    
    #if (et$xHasNA) {
    #  print("Note: Input matrix has NA. See ?extraTrees how ET builds trees with NAs.")
    #}
    
    ## setting variables:
    .jcall( et$jobject, "V", "setNumRandomCuts", as.integer(et$numRandomCuts) )
    .jcall( et$jobject, "V", "setEvenCuts", et$evenCuts )
    .jcall( et$jobject, "V", "setNumThreads", as.integer(et$numThreads) )
    .jcall( et$jobject, "V", "setHasNaN", et$xHasNA)
    
    ## if present set weights:
    if (et$useWeights) {
      .jcall( et$jobject, "V", "setWeights", .jarray(as.double(weights)) )
    }
    
    ## if given set bagging:
    if (et$useBagging) {
      if (length(bagSizes) == 1) {
        .jcall( et$jobject, "V", "setBagging", as.integer(bagSizes[1]) )
      } else {
        .jcall( et$jobject, "V", "setBagging", 
                .jarray(as.integer( bagSizes )), 
                .jarray(as.integer( as.integer(bagLabels)-1 )) 
              )
      }
    }
    
    ## multitask variables:
    if (et$multitask) {
        .jcall( et$jobject, "V", "setTasks", .jarray(as.integer(tasks-1)) )
        .jcall( et$jobject, "V", "setProbOfTaskCuts", et$probOfTaskCuts )
        .jcall( et$jobject, "V", "setNumRandomTaskCuts", as.integer(et$numRandomTaskCuts) )
    }
    
    ## learning the trees (stored at the et$jobject)
    .jcall( et$jobject, "V", "learnTrees", as.integer(et$nodesize), as.integer(et$mtry), as.integer(et$ntree) )
    
    return( et )
}

