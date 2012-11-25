

predict.extraTrees <- function( object, newdata, ... )
{
    if (!inherits(object, "extraTrees")) {
        stop("object not of class extraTrees")
    }
    et = object
        
    #getValues(Matrix input)
    if (ncol(newdata)!=et$ndim) {
        stop( sprintf("newdata(ncol=%d) does not have the same dimensions as the original x (ncol=%d)", ncol(newdata), et$ndim) )
    }
    if (!et$factor) {
        ## regression:
        return( .jcall( et$jobject, "[D", "getValues", toJavaMatrix(newdata) )
 )
    }
    ## classification:
    yhat = .jcall( et$jobject, "[I", "getValues", toJavaMatrix(newdata) )
    #return(yhat)
    return( factor(yhat, levels=0:(length(et$levels)-1), labels=et$levels ) )
}
