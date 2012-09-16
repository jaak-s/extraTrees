extraTrees.default <- function(x, y, ...){
   #hjw <- .jnew("ExtraTrees") # create instance of HelloJavaWorld class
   #out <- .jcall(hjw, "S", "sayHello") # invoke sayHello method
   #x <- as.matrix(x)
   #y <- as.numeric(y)
   
   x <- list()
   x$nmin = 2
   class(x) = "extraTrees"
   return( x )
}

