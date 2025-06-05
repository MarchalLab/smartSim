#' @title DataCharacteristics
#'
#' @description data characteristics based on existing Smart-seq3 data
#'
#' @param inDir directory containing output from characterizeData
#' @return object of class "DataCharacteristics"
#'
#' @importFrom dplyr %>% group_by n_distinct
#' @importFrom utils read.table 
#' @importFrom stats quantile 
#' @importFrom graphics hist
#' @export
DataCharacteristics <- function(inDir){
  cell <- gene <- UMI <- numReads <- NULL #binding variable locally to function
  
  #read data
  UMIcounts <- read.table(paste0(inDir, "/", "UMIcounts.tsv"), sep = "\t", header = TRUE)
  fragLen <- read.table(paste0(inDir, "/", "fragmentLength.tsv"), sep = "\t", header= TRUE)
  fragLen_UMI <- read.table(paste0(inDir, "/", "fragmentLength_UMI.tsv"), sep = "\t", header = TRUE)
  fragLen_nonUMI <- read.table(paste0(inDir, "/", "fragmentLength_nonUMI.tsv"), sep = "\t", header = TRUE)
  readQual <- read.table(paste0(inDir, "/", "readQual.tsv"), sep = "\t", header = TRUE, comment.char="", quote = NULL)
  bcQual <- read.table(paste0(inDir, "/", "readQual.tsv"), sep = "\t", header = TRUE, comment.char="", quote = NULL)
  intronicReads <- read.table(paste0(inDir, "/", "intronicReads.tsv"), sep = "\t", header = TRUE, comment.char = "", quote = NULL)

  # number of unique UMIs per gene
  UMIsPerGene <- UMIcounts[UMIcounts$UMI != "noUMI",] %>%
    group_by(cell, gene) %>%
    summarize(numUMI = n_distinct(UMI))
  q01 <- quantile(UMIsPerGene$numUMI, 0.01)
  q99 <- quantile(UMIsPerGene$numUMI, 0.99)
  UMIsPerGene_hist <- hist(UMIsPerGene$numUMI[UMIsPerGene$numUMI >= q01 & UMIsPerGene$numUMI <= q99],
                           breaks = seq(1, max(UMIsPerGene$numUMI)),
                           plot = FALSE)

  #number of reads per UMI
  readsPerUMI <- UMIcounts[UMIcounts$UMI != "noUMI",]
  q01 <- quantile(readsPerUMI$numReads, 0.01)
  q99 <- quantile(readsPerUMI$numReads, 0.99)
  readsPerUMI_hist <- hist(readsPerUMI$numReads[readsPerUMI$numReads >= q01 & readsPerUMI$numReads <= q99],
                           breaks = seq(1, max(readsPerUMI$numReads)),
                           plot = FALSE)

  #fragment length for UMI-containing reads
  fragLen_UMI <- fragLen_UMI[order(fragLen_UMI$val), ]
  fragLen_UMI$cumSum <- cumsum(fragLen_UMI$count)
  tot <- fragLen_UMI$cumSum[nrow(fragLen_UMI)]
  fragLen_UMI <- fragLen_UMI[fragLen_UMI$cumSum >= 0.01*tot, ] #q01
  fragLen_UMI <- fragLen_UMI[fragLen_UMI$cumSum <= 0.99*tot, ] #q99
  fragLen_UMI <- smoothData(fragLen_UMI)
  fragLen_UMI_hist <- list("mids" = fragLen_UMI$val, "counts" = fragLen_UMI$smooth)

  #fragment length for non UMI-containing reads
  fragLen_nonUMI <- fragLen_nonUMI[order(fragLen_nonUMI$val), ]
  fragLen_nonUMI$cumSum <- cumsum(fragLen_nonUMI$count)
  tot <- fragLen_nonUMI$cumSum[nrow(fragLen_nonUMI)]
  fragLen_nonUMI <- fragLen_nonUMI[fragLen_nonUMI$cumSum >= 0.01*tot, ] #q01
  fragLen_nonUMI <- fragLen_nonUMI[fragLen_nonUMI$cumSum <= 0.99*tot, ] #q99
  fragLen_nonUMI <- smoothData(fragLen_nonUMI)
  fragLen_nonUMI_hist <- list("mids" = fragLen_nonUMI$val, "counts" = fragLen_nonUMI$smooth)


  #distribution of number of nonUMI reads vs. UMI reads
  numNonUMIvsUMI <- UMIcounts %>%
    group_by(cell, gene) %>%
    summarize(numNonUMIvsUMI = sum(numReads[UMI == "noUMI"]) / sum(numReads[UMI != "noUMI"]))
  numNonUMIvsUMI <- numNonUMIvsUMI$numNonUMIvsUMI[! is.infinite(numNonUMIvsUMI$numNonUMIvsUMI)]
  q01 <- quantile(numNonUMIvsUMI, 0.01)
  q99 <- quantile(numNonUMIvsUMI, 0.99)
  nonUMIvsUMI_hist <- hist(numNonUMIvsUMI[numNonUMIvsUMI >= q01 & numNonUMIvsUMI <= q99],
                           breaks = seq(0, max(numNonUMIvsUMI), by=0.1),
                           plot = FALSE)

  #distribution of number of intronic vs exonic reads
  intronicFrac <- intronicReads[,"intronicReads"]/(intronicReads[,"intronicReads"] + intronicReads[,"exonicReads"])
  intronicFrac_hist <- hist(intronicFrac,
                            breaks = seq(0, 1, by=0.01),
                            plot = FALSE)
  
  #quality
  readQual$pos <- readQual$pos + 1 #R uses 1-based indexing
  bcQual$pos <- bcQual$pos + 1 #R uses 1-based indexing

  #create object
  obj <- list("fragLen_UMI" = fragLen_UMI_hist,
              "fragLen_nonUMI" = fragLen_nonUMI_hist,
              "UMIsPerGene" = UMIsPerGene_hist,
              "readsPerUMI" = readsPerUMI_hist,
              "nonUMIvsUMI" = nonUMIvsUMI_hist,
              "intronicFrac" = intronicFrac_hist,
              "readQual" = readQual,
              "bcQual" = bcQual)
  class(obj) <- "DataCharacteristics"
  return(obj)

}

#helper functions
ma <- function(x, n = 11){
  smoothVal <- as.vector(filter(x, rep(1 / n, n), sides = 2))
  
  smoothVal[1:as.integer(n/2)] = mean(x[1:as.integer(n/2)]) #set first values
  smoothVal[(length(x)-(as.integer(n/2)+1)):length(x)] = mean(x[(length(x)-(as.integer(n/2)+1)):length(x)]) #set last values
  return(smoothVal)
}

smoothData <- function(countTable){
  #create table with all possible values
  fullTable<- data.frame(count = rep(0, max(countTable$val)-min(countTable$val)+1), 
                         row.names = min(countTable$val):max(countTable$val))
  fullTable[as.character(countTable$val), "count"] = countTable$count

  fullTable$smooth <- ma(fullTable$count, 11)
  fullTable$val <- as.numeric(row.names(fullTable))
  return(fullTable[, c("val", "smooth")])
}


#' @title sampleUMIcounts
#'
#' @description sample n elements from the distribution of UMIcounts
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minCount minimum number of counts
#'
#' @export
sampleUMIcounts <- function(dataChar, n, minCount = 0){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$UMIsPerGene$mids[dataChar$UMIsPerGene$mids > minCount] - 0.5
  counts <- dataChar$UMIsPerGene$counts[dataChar$UMIsPerGene$mids > minCount]
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
}

#' @title sampleNumReadsPerUMI
#'
#' @description sample n elements from the distribution of number of reads per UMI
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleNumReadsPerUMI <- function(dataChar, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$readsPerUMI$mids - 0.5
  counts <- dataChar$readsPerUMI$counts

  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
}

#' @title sampleUMIfragLen
#'
#' @description sample n elements from the distribution of fragment length for reads containing UMIs
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minLen minumum fragment length
#'
#' @export
sampleUMIfragLen <- function(dataChar, n, minLen){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$fragLen_UMI$mids
  counts <- dataChar$fragLen_UMI$counts
  
  return(sample(midpoints[midpoints >= minLen], size = n, replace = TRUE, prob = counts[midpoints >= minLen]))

}

#' @title sampleNonUMIfragLen
#'
#' @description sample n elements from the distribution of fragment length for reads not containing UMIs
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#' @param minLen minumum fragment length
#'
#' @export
sampleNonUMIfragLen <- function(dataChar, n, minLen){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$fragLen_nonUMI$mids
  counts <- dataChar$fragLen_nonUMI$counts

  return(sample(midpoints[midpoints >= minLen], size = n, replace = TRUE, prob = counts[midpoints >= minLen]))

}

#' @title sampleNonUMIvsUMI
#'
#' @description sample n elements from the distribution of number of nonUMI reads vs. UMI reads
#'
#' @param dataChar DataCharacteristics object
#' @param n number of elements to sample
#'
#' @export
sampleNonUMIvsUMI <- function(dataChar, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")

  midpoints <- dataChar$nonUMIvsUMI$mids
  counts <- dataChar$nonUMIvsUMI$counts

  return(sample(midpoints, size = n, replace = TRUE, prob = counts))

}

#' @title sampleReadQual
#'
#' @description sample read quality for read of length rLen
#'
#' @param dataChar DataCharacteristics object
#' @param rLen read length
#' @param n number of read qualities to sample
#'
#' @export
sampleReadQual <- function(dataChar, rLen, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(rLen) | ! rLen%%1 == 0 | ! rLen>0) stop("value pos = ", rLen, " is not valid")

  maxPos <- max(dataChar$readQual$pos)

  if (rLen > maxPos) warning("original data did not contain reads with length ", rLen)

  qualStr <- matrix("", nrow = rLen, ncol = n)
  for (pos in 1:rLen){

    if (pos <= maxPos){
      midpoints <- dataChar$readQual$qual[dataChar$readQual$pos == pos]
      counts <- dataChar$readQual$count[dataChar$readQual$pos == pos]
    }else{
      midpoints <- dataChar$readQual$qual
      counts <- dataChar$readQual$count
    }
    qualStr[pos,] = sample(midpoints, size = n, prob = counts, replace = TRUE)
  }

  return(apply(qualStr, 2, paste0, collapse = ""))

}


#' @title sampleBCQual
#'
#' @description sample barcode quality for barcode of length rLen
#'
#' @param dataChar DataCharacteristics object
#' @param rLen length of barcode
#' @param n number of barcode qualities to sample
#'
#' @export
sampleBCQual <- function(dataChar, rLen, n){

  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(rLen) | ! rLen%%1 == 0 | ! rLen>0) stop("value pos = ", rLen, " is not valid")

  maxPos <- max(dataChar$bcQual$pos)

  if (rLen > maxPos) warning("original data did not contain barcodes with length ", rLen)

  qualStr <- matrix("", nrow = rLen, ncol = n)
  for (pos in 1:rLen){

    if (pos <= maxPos){
      midpoints <- dataChar$bcQual[dataChar$bcQual$pos == pos, "qual"]
      counts <- dataChar$bcQual[dataChar$bcQual$pos == pos, "count"]
    }else{
      midpoints <- dataChar$bcQual$qual
      counts <- dataChar$bcQual$count
    }
    qualStr[pos, ] = sample(midpoints, size = n, prob = counts, replace = TRUE)
  }
  return(apply(qualStr, 2, paste0, collapse = ""))

}

#' @title sampleUnspliced
#'
#' @description sample fraction of unspliced transcripts
#'
#' @param dataChar DataCharacteristics object
#' @param n number of sample fractions to sample
#'
#' @export
sampleUnspliced <- function(dataChar, n){
  
  #param checks
  if (!inherits(dataChar, "DataCharacteristics")) stop("Not a DataCharacteristics object")
  if (! is.numeric(n) | ! n%%1 == 0 | ! n>0) stop("value n = ", n, " is not valid")
  
  midpoints <- dataChar$intronicFrac$mids - 0.005
  counts <- dataChar$intronicFrac$counts
  
  return(sample(midpoints, size = n, replace = TRUE, prob = counts))
}
