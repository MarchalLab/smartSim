#' @title getTranscriptExpression()
#'
#' @description get transcript expression from cluster-specific PSI values and fractions of unspliced transcripts, and cell-specific gene expression values
#' 
#' @param clusterPSI transcript PSI values per cluster
#' @param unspliced fractions of unspliced transcripts per cluster
#' @param geneExpression gene expression for each cell
#' @param outFile output file to write transcript expression values to
#'
#' @return matrix containing expression values for each transcript (incl. unspliced transcripts) for each gene and cell 
#'
#' @export
getTranscriptExpression <- function(clusterPSI, unspliced, geneExpression, outFile){
  
  #rescale spliced transcripts
  rescaled <- merge(clusterPSI, unspliced , by = c("cluster", "gene_id"))
  rescaled$PSI <- rescaled$PSI * (1-rescaled$nUnspliced)
  rescaled$nUnspliced <- NULL
  
  #add unspliced transcripts
  rowsToAdd <- unspliced[unspliced$nUnspliced > 0,]
  if (nrow(rowsToAdd) > 0){
    rowsToAdd$transcript_id <- paste0(rowsToAdd$gene_id, "_unspliced")
    names(rowsToAdd)[names(rowsToAdd) == 'nUnspliced'] <- 'PSI'
  }
  transcriptPSI <- rbind(rescaled, rowsToAdd)
  
  #combine cluster specific transcript expression with gene expression
  transcriptExpression <- merge(transcriptPSI, geneExpression, by = c("gene_id", "cluster"))
  
  #round transcript expression to nearest integer
  transcriptExpression$transcriptExpression <- round(transcriptExpression$geneExpression * transcriptExpression$PSI)
  
  #write to output
  write.table(transcriptExpression, outFile, 
              sep = "\t")
  
  return(transcriptExpression)
  
}
