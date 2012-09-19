print.extraTrees <- function(x, ...) {
    cat( "ExtraTrees:\n" )
    cat( sprintf(" - # of trees: %d\n", x$ntree) )
    cat( sprintf(" - node size:  %d\n", x$nodesize) )
    cat( sprintf(" - # of dim:   %d\n", x$ndim) )
}

