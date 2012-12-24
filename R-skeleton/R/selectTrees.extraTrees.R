## creates a new extraTrees object based on selection
selectTrees.extraTrees <- function( object, selection ) {
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
