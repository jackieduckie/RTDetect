#' Ranking matching transcripts
#' @details
#' This is an internal function which returns overlapping transcript names with ranking scores. 
#' The ranking score is the proportion of exon-exon fusions (intronic deletion events) detected for a given transcript.
#' @param genes TxDb object of genes. hg19 and hg38 are supported in the current version.
#' @param transcripts.col A vector of transcript names.
#' @return A dataframe with two columns, tx_name and score. 
.scoreByTranscripts <- function(genes, transcripts.col){
    overlapIntron.df <- as.data.frame(table(transcripts.col)/2)
    colnames(overlapIntron.df) <- c("tx_name", "count")
    overlapIntron.df <- merge(overlapIntron.df, 
                              mcols(GenomicFeatures::transcripts(genes, columns=c("tx_name","exon_rank"), 
                                                                 filter=list(tx_name=overlapIntron.df[,1]))))
    return(data.frame(tx_name=overlapIntron.df$tx_name, 
                      score= overlapIntron.df$count / (sapply(overlapIntron.df$exon_rank, length)-1)))
}

