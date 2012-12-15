

predict.extraTrees <- function( object, newdata, quantile=NULL, ... )
{
    if (!inherits(object, "extraTrees")) {
        stop("object not of class extraTrees")
    }
    et = object
        
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
        ## quantile regression:
        return( .jcall( et$jobject, "[D", "getQuantiles", toJavaMatrix(newdata), quantile[1] ) )
    }
    if (!et$factor) {
        ## regression:
        return( .jcall( et$jobject, "[D", "getValues", toJavaMatrix(newdata) ) )
    }
    ## classification:
    yhat = .jcall( et$jobject, "[I", "getValues", toJavaMatrix(newdata) )
    return( factor(yhat, levels=0:(length(et$levels)-1), labels=et$levels ) )
}
