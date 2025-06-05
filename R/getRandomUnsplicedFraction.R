#' @title getRandomUnspliced()
#'
#' @description generate random fractions of unspliced transcripts based on the distribution of unspliced transcripts from data.
#' 
#' @param nClusters number of clusters
#' @param fixedGenes genes that have the unspliced fraction across clusters
#' @param diffGenes genes that have different unspliced fraction across clusters
#' @param dataChar data characteristics object
#' @param seed for reproducability [default = 1234]
#'
#' @return matrix containing fractions of unspliced transcripts per cluster
#'
#' @export
getRandomUnspliced <- function(nClusters, fixedGenes, diffGenes, dataChar, seed = 1234){
  
  print("getting random fractions of unspliced transcript")
  set.seed(seed)
  
  if (nClusters < 1 | nClusters%%1 != 0)
    stop("nClusters = ", nClusters, " is not a valid number")
  
  unspliced_fixed <- data.frame(cluster = rep(1:nClusters, each = length(fixedGenes)), 
                                gene_id = rep(fixedGenes, times = nClusters), 
                                nUnspliced = rep(sampleUnspliced(dataChar, n = length(fixedGenes)), times = nClusters))
  
  unspliced_diff <- data.frame(cluster = rep(1:nClusters, each = length(diffGenes)),
                               gene_id = rep(diffGenes, times = nClusters), 
                               nUnspliced = sampleUnspliced(dataChar, n = length(diffGenes) * nClusters))
  
  unspliced = rbind(unspliced_fixed, unspliced_diff)
  
  return(unspliced)
}
