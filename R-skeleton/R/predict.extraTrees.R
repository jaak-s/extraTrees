

predict.extraTrees <- function( object, newdata, quantile=NULL, allValues=F, newtasks=NULL, ... )
{
    if (!inherits(object, "extraTrees")) {
        stop("Object not of class extraTrees")
    }
    et = object
    
    if (et$multitask && is.null(newtasks)) {
        stop("Prediction requires list of newtasks.")
    }
    
    if (!et$multitask && ! is.null(newtasks)) {
        stop("to predict with newtasks extraTrees must be trained with tasks")
    }
        
    if (ncol(newdata)!=et$ndim) {
        stop( sprintf("newdata(ncol=%d) does not have the same dimensions as the original x (ncol=%d)", ncol(newdata), et$ndim) )
    }
    if (!is.null(quantile)) {
        if (!et$quantile) {
            stop("to predict quantiles, extraTrees must be trained with quantile=T")
        }
        if (quantile[1] < 0.0 || quantile[1] > 1.0) {
            stop("quantile has to be between 0.0 and 1.0.")
        }
        if (allValues) {
            stop("Can't use allValues=T with quantile.")
        }
        ## quantile regression:
        return( .jcall( et$jobject, "[D", "getQuantiles", toJavaMatrix(newdata), quantile[1] ) )
    }
    if (allValues) {
        ## returning allValues prediction:
        if (et$multitask) {
            ## multi-task version:
            m = toRMatrix( .jcall( 
                et$jobject,
                "Lorg/extratrees/Matrix;", 
                "getAllValuesMT", 
                toJavaMatrix(newdata),
                .jarray(as.integer(newtasks-1))
            ))
        } else {
            m = toRMatrix( .jcall( et$jobject, "Lorg/extratrees/Matrix;", "getAllValues", toJavaMatrix(newdata) ) )
        }
        if (!et$factor) {
            ## regression model:
            return(m)
        }
        ## factor model: convert double matrix into data.frame of factors
        lvls = 0:(length(et$levels)-1)
        m = round(m)
        mlist = lapply( 1:ncol(m), function(j) factor(round(m[,j]), levels=lvls, labels=et$levels) )
        ## changing list to data.frame:
        attributes(mlist) <- list(
            row.names=c(NA_integer_,nrow(m)), 
            class="data.frame",
            names=make.names(names(mlist)),
            unique=TRUE
        )
        return(mlist)
    }
    if (!et$factor) {
        ## multi-task regression:
        if (et$multitask) {
            return( .jcall( et$jobject, "[D", "getValuesMT", toJavaMatrix(newdata), .jarray(as.integer(newtasks-1)) ) )
        }
        ## single-task regression:
        return( .jcall( et$jobject, "[D", "getValues", toJavaMatrix(newdata) ) )
    }
    if (et$multitask) {
        ## multi-task classification:
        yhat = .jcall( et$jobject, "[I", "getValuesMT", toJavaMatrix(newdata), .jarray(as.integer(newtasks-1)) )
    } else {
        ## single-task classification:
        yhat = .jcall( et$jobject, "[I", "getValues", toJavaMatrix(newdata) )
    }
    return( factor(yhat, levels=0:(length(et$levels)-1), labels=et$levels ) )
}
