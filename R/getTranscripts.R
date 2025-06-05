
### helper functions

filtSmallTranscripts <- function(gtf_attributes, minTxLen){
  transcript_id <- NULL #binding variable locally to function
  
  length <- gtf_attributes[gtf_attributes$type == "exon",] %>%
    group_by(transcript_id) %>%
    summarize(length = sum(width))
  transcripts_toRemove <- length$transcript_id[length$length < minTxLen]
  return (gtf_attributes[! gtf_attributes$transcript_id %in% transcripts_toRemove, ])
}

filtEnoughTx <- function(gtf_attributes, nTx){
  gene_id <- transcript_id <- NULL #binding variable locally to function
  
  nTxPerGene <- gtf_attributes[gtf_attributes$type == "transcript",] %>%
    group_by(gene_id) %>%
    summarize(nTx = n_distinct(transcript_id))
  genes_toKeep <- nTxPerGene$gene_id[nTxPerGene$nTx >= nTx]
  return (gtf_attributes[gtf_attributes$gene_id %in% genes_toKeep, ])
}

selectGenes <- function(gtf_attributes, genes, nRandomGenes){
  if (! is.null(nRandomGenes)){
    allGenes <- unique(gtf_attributes$gene_id)
    if (nRandomGenes > length(allGenes)){
      warning("could not select ", nRandomGenes, " genes from gtf, not enough genes. All genes (=",
              length(allGenes), ") genes selected instead.")
      nRandomGenes <- length(allGenes)
    }
    genes <- sample(unique(gtf_attributes$gene_id), nRandomGenes, replace = F)
  }
  if (any(! genes %in% gtf_attributes$gene_id)){
    warning("The following genes are ignored as they are not present in the reference genome, not long enough, or do not have enough annotated transcripts: ",
            paste(genes[! genes %in% gtf_attributes$gene_id], collapse = " "))
  }

  return(gtf_attributes[gtf_attributes$gene_id %in% genes, ])
}

writeGeneInfo <- function(gtf_attributes, out){
  geneInfo <- gtf_attributes[gtf_attributes$type == "gene",
                             c("seqnames", "start", "end", "gene_id")]
  write.table(geneInfo, file = paste0(out, "/geneInfo.tsv"), sep = "\t", row.names = F, col.names = F, quote = F)
}

selectTranscripts <- function(gtf_attributes, nTx){
  gene_id <- NULL #binding variable locally to function
  
  #select n transcripts for each gene
  transcripts <- gtf_attributes[gtf_attributes$type == "transcript",c("gene_id", "transcript_id")]
  transcripts <- transcripts %>%
    group_by(gene_id) %>%
    slice_sample(n = nTx)
  rowsToKeep = ((gtf_attributes$transcript_id %in% transcripts$transcript_id) & (gtf_attributes$type %in% c("transcript", "exon"))) | (gtf_attributes$type == "gene")
  return(gtf_attributes[rowsToKeep,])
}

###

#' @title getTranscripts()
#'
#' @description get transcripts to simulate
#'
#' @param gtf gtfFile containing the reference annotation
#' @param genes list of genes
#' @param transcripts list of transcripts
#' @param nRandomGenes number of random genes for which transcripts are selected (not used if list of genes or transcripts given)
#' @param nTx number of transcripts per gene (not used if list of transcripts given
#' @param out outdir to write geneInfo file
#' @param minTxLen minimum transcript length [default = 1000]
#' @param seed seed for reproducability [default = 1234]
#'
#' @return gtf table containing all selected genes and transcripts
#'
#' @importFrom rtracklayer import
#' @importFrom dplyr %>% group_by summarize n_distinct slice_sample
#' @importFrom utils write.table
#'
#' @export
getTranscripts <- function(gtf, nTx = NULL, out, genes = NULL, transcripts = NULL, nRandomGenes = NULL, minTxLen = 1000, seed = 1234){
  
  #check arguments
  if (is.null(nTx) & is.null(transcripts)){
    stop("provide list of transcripts or number of transcripts per gene")
  } 
  if (is.null(transcripts) & is.null(genes) & is.null(nRandomGenes)){
    stop("provide list of genes or number of random genes")
  }
  if (! is.null(transcripts)){
    if (! is.null(nTx)){
      warning("argument nTx is ignored, using predefined list of transcripts provided in transcripts")
    }
    if (! is.null(genes)){
      warning("argument genes is ignored, using predefined list of transcripts provided in transcripts")
    }
    if (! is.null(nRandomGenes)){
      warning("argument nRandomGenes is ignored, using predefined list of transcripts provided in transcripts")
    }
  }
  if (! is.null(genes) & ! is.null(nRandomGenes)){
    warning("argument nRandomGenes is ignored, using predefined list of genes provided in genes")
  }

  print("Selecting transcripts for simulation")
  set.seed(seed)

  #create output directory
  dir.create(out, showWarnings = FALSE)

  #import reference data
  gtf_attributes <- as.data.frame(import(gtf))

  #select transcripts that are large enough (> 1kb)
  gtf_attributes <- filtSmallTranscripts(gtf_attributes, minTxLen)

  if (! is.null(transcripts)){
    
    genesToKeep = unique(gtf_attributes$gene_id[gtf_attributes$transcript_id %in% transcripts])
    rowsToKeep = ((gtf_attributes$transcript_id %in% transcripts) & (gtf_attributes$type %in% c("transcript", "exon"))) | 
      ((gtf_attributes$gene_id %in% genesToKeep) & (gtf_attributes$type == "gene"))
    
    gtf_attributes <- gtf_attributes[rowsToKeep,]

  }else{
    #select genes with enough transcripts
    gtf_attributes <- filtEnoughTx(gtf_attributes, nTx)
    
    #select random genes, or use predefined list
    gtf_attributes <- selectGenes(gtf_attributes, genes, nRandomGenes)
    
    #select random transcripts
    gtf_attributes <- selectTranscripts(gtf_attributes, nTx)
  }
 
  writeGeneInfo(gtf_attributes, out)
  return(gtf_attributes)

}




