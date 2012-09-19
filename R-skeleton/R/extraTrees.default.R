

## converts Java matrix to R matrix (doubles):
toRMatrix <- function( javam ) {
    return( matrix( 
        .jfield(javam, "[D", "v" ),
        .jfield(javam, "I", "nrows" ),
        .jfield(javam, "I", "ncols" )
    ))
}

## converts R matrix into Java matrix (doubles)
toJavaMatrix <- function( m ) {
    return(.jnew(
        "Matrix", .jarray(m), nrow(m), ncol(m)
    ))
}


extraTrees.default <- function(x, y, xtest=NULL, ytest=NULL, ntree=500,
             mtry = if (!is.null(y) && !is.factor(y))
                    max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))),
             nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1
){
   #hjw <- .jnew("ExtraTrees") # create instance of HelloJavaWorld class
   #out <- .jcall(hjw, "S", "sayHello") # invoke sayHello method
   #x <- as.matrix(x)
   #y <- as.numeric(y)
    n <- nrow(x)
    p <- ncol(x)
    if (n == 0) stop("data (x) has 0 rows")
    x.row.names <- rownames(x)
    x.col.names <- if (is.null(colnames(x))) 1:ncol(x) else colnames(x)

    testdat <- !is.null(xtest)
    if (testdat) {
        if (ncol(x) != ncol(xtest))
            stop("x and xtest must have same number of columns")
        ntest <- nrow(xtest)
        xts.row.names <- rownames(xtest)
    }
    
    ## remove when factor classification is implemented:
    if (is.factor(y)) {
        stop("classification with extraTrees is not yet implemented.")
    }
    
    et <- list()
    et$ntree    = ntree
    et$nodesize = nodesize
    et$ndim     = p
    et$mtry     = mtry
    class(et) = "extraTrees"
    
    ## creating ExtraTree object with the data
    et$jobject  = .jnew(
        "ExtraTrees",
        toJavaMatrix(x),
        .jarray(y)
    )
    ## learning the trees (stored at the et$jobject)
    .jcall( et$jobject, "V", "learnTrees", as.integer(et$nodesize), as.integer(et$mtry), as.integer(et$ntree) )
    
    return( et )
}

