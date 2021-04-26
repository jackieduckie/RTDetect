#' Combining matching transcripts 
#' @details
#' This is an internal function used to merge all overlapping transcripts of a breakpoint into one vector.
#' @param gr A GRanges object
#' @param names A vector of granges names.
#' @return A list of vectors. Each vector is named with the name of the corresponding granges.
.combineMatchingTranscripts <- function(gr, names){
    names <- unique(names)
    txs.list <- vector(mode="list", length=length(names))
    names(txs.list) <- names
    for (name in names) {
        #txs.list[[name]] <- name
        txs.list[[name]] <- Reduce(union, gr[names(gr) == name]$txs)
    }
    return(txs.list)
}