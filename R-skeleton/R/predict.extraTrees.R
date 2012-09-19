

predict.extraTrees <- function( et, newdata )
{
    #getValues(Matrix input)
    if (ncol(newdata)!=et$ndim) {
        stop( sprintf("newdata(ncol=%d) does not have the same dimensions as the original x (ncol=%d)", ncol(newdata), et$ndim) )
    }
    
    return( .jcall( et$jobject, "[D", "getValues", toJavaMatrix(newdata) ) )
}
