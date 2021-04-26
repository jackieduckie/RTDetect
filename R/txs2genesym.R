#' Adding gene symbol annotations
#' @details 
#' This is an internal function which takes a list of txs in UCSC id format as input 
#' and convert the txs to gene symbol.
#' @param txs A list of transcript ids in UCSC format.
#' @param unique.genesyms TRUE or FALSE. If TRUE, the converted gene symbols will remove duplicates.
#' @return A list of names in gene symbols
.txs2genesym <- function(txs, unique.genesyms=TRUE){
    assertthat::assert_that(class(txs)=="list", msg = "txs should be a list object")
    #ucsc id to gene symbol look up table
    gene_symbol <- utils::read.delim(system.file("extdata", "gene_symbol.txt", package = "RTDetect"), 
                                     header=TRUE, comment.char="#")
    gene_symbol <- dplyr::bind_rows(gene_symbol, data.frame(kgID=NA, geneSymbol=NA))
    
    #txs <- RT$insSite$txs
    gene_syms <- lapply(txs, function(x) gene_symbol$geneSymbol[gene_symbol$kgID %in% x])
    if (unique.genesyms) {
        gene_syms <- lapply(gene_syms, unique)
    }
    return(gene_syms)
}
